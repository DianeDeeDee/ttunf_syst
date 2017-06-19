#include "WjetsCorrections/HFsys.h"

#include<iostream>
#include<cmath>
#include<cstdlib>

using namespace std;

void SetWflavors_elec(int idsys, TString sysname, int ijet, double Wjets_in[5],double Wjets_out[5], double& canorm, int _mode) {
	//-
	// This function corrects the Wjet ligth and heavy flavor components and return the normalization based from Charge Assymmetry.
	//
	// NOTE: this function shoud always be called with PRETAGGED event counts. NEVER WITH TAGGED counts.
	//
	// INPUTS:
	// for idsys there are two input modes: the individul systematics (long list), or the combined systematics (4 checks).
	// Long list:
	// input idsys>0..999 or sysname are passed to the function GetFFactors_elec (see below) that returns the relevant HF factors.
	// This allows to use the full list of systematics.
	// Short list:
	// most analysers will use the total sum of all systematics. In that case only four checks up and two down are needed.
	// Complete set: idsys=2000 (up), -2000 (down) ---> checks Fbb + Fcc + Fc versus Fl
	//                     2001 (up), -2001 (down) ---> checks Fbb+Fcc versus Fc
	//                     2002 (up), -2002 (down) ---> checks Fll  OBSOLETE, not needed anymore, return nominal
	//                     2003 (up), -2003 (down) ---> check Fbb+25% in 1,3,4,5 jetbin. --> see remark *
	//                     2004 (up), -2004 (down) ---> check Fc+25% in 1,3,4,5 jetbin.  --> see remark *
	//                     2005 (up), -2005 (down) ---> check canorm. Please note that canorm is returned and NOT yet applied to Wjets_out.
	// * please note that the 25% per jetbin change should be applied UNCORRELATED to all jetbins. So, for example, when you use simultaneoulsy 
	//   jetbin 4 and 5 in your analysis, then first use the changed 4jetbin HF factors and 5 nominal and then change 5jetbin and keep 4jet nominal.
	//  
	//
	// input ijet = jetbin you want to set. Available ijet values: 1=1ex, 2=2ex, 3=3ex, 4=4ex, -5=5=5in, -3=3in, -4=4in.
	// input Wjets_in are the PRETAGGED event numbers (for ijet) for 0:Wbb, 1:Wcc, 2:Wc, 3:Wlight
	// The flavor fractions defined as bb, cc, c and light flavor jets correspond to the HFOR flag equals to 0, 1, 2 and 3 respectively
	//--
	// OUTPUTS:
	// output Wjets_out are the PRETAGGED event counts (for ijet) for 0:Wbb, 1:Wcc, 2:Wc, 3:Wlight 
	// output canorm, which is  NOT yet applied to Wjets_out. This normlization factor can (or should be) be used for pretagged and tagged events.
	//
	// NOTE: if you want to use the HF corrections for your tagged sample OR after any other cut, then you still have to input
	//       the PRETAGGED number of events and use the pretagged fraction: Wjets_out[iflavor]/Wjets_in[iflavor] to reweigh your events.
	//
	//
	// TYPICAL USE
	// LONG:
	// Typical use for the full systematics list:
	//      .....
	//      int idsys=0; // for nominal, but you have to loop over many idsys, see in explanation on www.nikhef.nl/~h73/factors.html
	//      int ijet=3; // in case you input the 3jetbin
	//      double Wjets_in[4]={ 6888, 14281, 22937, 102347} ; 
	//      double Wjets_out[4]; 
	//      SetWflavors(idsys, "" , ijet, Wjets_in, Wjets_out, canorm);
	//      Wjets_out[i]=canorm*Wjets_out[i] .....
	//      please note that you have to do the 25% change per jetbin per flavor, which is NOT available in the long list. For these check you can use idsys=+/-2003,2004.
	//      Another note: you can only use this method when you have all the systematic changes to MC samples available. For example, when you apply the JESup systematic to your Wjets samples, you use:
	//      SetWflavors(-1, "jes_up" , ijet, Wjets_in_jesup, Wjets_out_jesup, canorm_jesup);
	//      Wjets_out_jesup[i]=canorm_jesup*Wjets_out_jesup[i];
	//      and thus (symbolically), the result is:
	//      Wjets_out_jesup= HFjesup x CAjesup x Wjets_in_jesup
	//      These individual systematics should not be applied to Wjets_nominal --> then you can use the short version as approximation.
	//      
	// SHORT:
	// Typical use when the combined uncertainties (short list)  are used:
	//      .....
	//      int ijet=3; // in case you input the 3jetbin
	//      int idsys=2000; // 
	//      double Wjets_in[4]={ 6888, 14281, 22937, 102347} ;
	//      double Wjets_out[4];
	//      SetWflavors(idsys, "" , ijet, Wjets_in, Wjets_out, canorm);
	//      Wjets_out[i]=canorm*Wjets_out[i] .....
	//      do this also for idsys=2001, 2002, 2003, 2004 and all negative: -2000, -2001, -2002, -2003, -2004. 
	//      if you also need to check the normalization (when you use canorm), then use 2005, -2005 as check.
	//      add all positive effects on you analysis in quadrature, and all negative effect also, as usual. Please note that negative effects on your analsysis
	//      do not neccessarily correspond with the negative idsys.
	//      See also remark * for checks 2003 and 2004.
	//      .....
	//
	// FAQ: I don't have the PRETAGGED numbers available at all times, but I want to apply the fraction to my TAGGED results.
	//      DON'T use this function 'SetWflavors' with your number of TAGGED events.
	//      If possible, once call the function SetWflavor with your PRETAGGED numbers and store the ratio's Wjets_in/Wjets_out -->
	//      these ratio's should be applied to your TAGGED samples.
	//
	// get HF Kfactors
	_mode=0; // should be ZERO
	if (_mode!=0) cout << " SetWFlavors: mode should be zero " << _mode << endl;
	double Hbb[9], Hcc[9], Hc[9],  Hll[9];
	double Fbb[9], Fcc[9], Fc[9],  Fll[9];
	double cafac[9];
	double fbb[9], fcc[9], fc[9],  fll[9];
	double winpretag[9], wintag[9], woutpretag[9],  wouttag[9];
        if (ijet==-3) { ijet=6; }
        if (ijet==-4) { ijet=7; }
        if (ijet==-5) { ijet=5;  } // default inc.
        if (ijet<=0 || ijet>7) cout << "FATAL  HF KFACTORS dont know which jetbin you are going to use " << endl;

	if (idsys*idsys<999*999) { // go for the full list of systematics
		GetFFactors_elec( idsys, sysname,  Hbb,Hcc, Hc,  Hll, cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
	} else if (idsys*idsys<1999*1999) {
		cout << " FATAL  HF KFACTORS probabbly called with old defintion of systematics ids " << endl;
	} else  {                  // go for the reduced errors.
		GetFFactors_elec( 0, "",  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
	}

	// errors on fraction, not on K factors
	double sigma_hf[2]={ 0.21 ,  - 0.20}; // error on ttoal HF
	double sigma_as[2] ={ 0.21 ,  - 0.21}; // error on fbb (and fcc) wrt to fc
	double sigma_ll[2] ={0.,  - 0.}; // -100% correlated with error on total HF
	double sigma_canorm_up[9]=  {0,  0.066, 0.060, 0.09, 0.11, 0.17, 0.08,0.099,0};
	double sigma_canorm_down[9]={0, -0.061, -0.061,-0.085,-0.1,-0.17, -0.076,-0.091,0};


	// use 2jetbin as a base for the default user.
	Fbb[ijet]=Hbb[2];
	Fcc[ijet]=Hcc[2];
	Fc[ijet]=Hc[2];
	Fll[ijet]=Hll[2];
	canorm=cafac[ijet];
	Wjets_out[0]=Fbb[ijet]*Wjets_in[0];
	Wjets_out[1]=Fcc[ijet]*Wjets_in[1];
	Wjets_out[2]=Fc[ijet] *Wjets_in[2];
	Wjets_out[3]=Fll[ijet]*Wjets_in[3];
	double nw_in=0;
	double nw_raw=0;
	for (int iw=0;iw<4;iw++) {
		nw_in+=Wjets_in[iw];
		nw_raw+=Wjets_out[iw];
	}
	if (nw_raw==0) { nw_in=1; nw_raw=1;} // some protection
	for (int iw=0;iw<4;iw++) {
		Wjets_out[iw]*=nw_in/nw_raw;
	}
	// done setting HFs
	// now systematics if combined mode is used.
	double fbb_in=Wjets_in[0]/nw_in;
	double fcc_in=Wjets_in[1]/nw_in;
	double fc_in= Wjets_in[2]/nw_in;
	double fll_in=Wjets_in[3]/nw_in;
	double fhf_in=fbb_in+fcc_in+fc_in;
	//   cout << " fraction in " << fbb_in <<  " " << fcc_in << " " << fc_in << " " << fll_in << endl;
	double fbb_out=Wjets_out[0]/nw_in;
	double fcc_out=Wjets_out[1]/nw_in;
	double fc_out= Wjets_out[2]/nw_in;
	double fll_out=Wjets_out[3]/nw_in;
	double fhf_out=fbb_out+fcc_out+fc_out;

	// total HF change
	if (idsys==2000) {
		fbb_out*=1+sigma_hf[0];
		fcc_out*=1+sigma_hf[0];
		fc_out*=1+sigma_hf[0];
		fhf_out=fbb_out+fcc_out+fc_out;
		fll_out=1-fhf_out;
	} else if (idsys==-2000) {
		fbb_out*=1+sigma_hf[1];
		fcc_out*=1+sigma_hf[1];
		fc_out*=1+sigma_hf[1];
		fhf_out=fbb_out+fcc_out+fc_out;
		fll_out=1-fhf_out;
	}

	// total HF constant, but bb+cc change wrt c
	if (idsys==2001) {
		fbb_out*=1+sigma_as[0];
		fcc_out*=1+sigma_as[0];
		fc_out=fhf_out-fbb_out-fcc_out; 
		// fll_out unchanged
	} else if (idsys==-2001) {
		fbb_out*=1+sigma_as[1];
		fcc_out*=1+sigma_as[1];
		fc_out=fhf_out-fbb_out-fcc_out;
		// fll_out unchanged
	}

	if (idsys==2002) {
		// dummy
	} else if (idsys==-2002) {
		// dummy
	} 
	TString cname;
	if (ijet!=2 && idsys==2003) { // change bb,cc at the cost of fc+fl
		fbb_out*=1+0.25;
		fcc_out*=1+0.25;
		double rest=fc_out+fll_out;
		double rest_new=1-fbb_out-fcc_out;
		fc_out*=rest_new/rest;
		fll_out*=rest_new/rest;
		if (ijet==1) cname="WbbWccjet1_up";
		if (ijet==3 || ijet==6) cname="WbbWccjet3_up";
		if (ijet==4 || ijet==7) cname="WbbWccjet4_up";
		if (ijet==5) cname="WbbWccjet5_up"; 
		GetFFactors_elec( -1, cname,  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		canorm=cafac[ijet];
	} else if (ijet!=2 && idsys==-2003) {
		fbb_out*=1-0.25;
		fcc_out*=1-0.25;
		double rest=fc_out+fll_out;
		double rest_new=1-fbb_out-fcc_out;
		fc_out*=rest_new/rest;
		fll_out*=rest_new/rest;
		if (ijet==1) cname="WbbWccjet1_down";
		if (ijet==3 || ijet==6) cname="WbbWccjet3_down";
		if (ijet==4 || ijet==7) cname="WbbWccjet4_down";
		if (ijet==5) cname="WbbWccjet5_down"; 
		GetFFactors_elec( -1, cname,  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		canorm=cafac[ijet];
	}
	if (ijet!=2 && idsys==2004) { // change c. at the cost of fbb+fcc+fl
		fc_out*=1+0.25;
		double rest=fbb_out+fcc_out+fll_out;
		double rest_new=1-fc_out; 
		fbb_out*=rest_new/rest;
		fcc_out*=rest_new/rest;
		fll_out*=rest_new/rest;
		if (ijet==1) cname="Wcjet1_up";
		if (ijet==3 || ijet==6) cname="Wcjet3_up";
		if (ijet==4 || ijet==7) cname="Wcjet4_up";
		if (ijet==5) cname="Wcjet5_up"; // 
		GetFFactors_elec( -1, cname,  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		canorm=cafac[ijet];
	} else if (ijet!=2 && idsys==-2004) {
		fc_out*=1-0.25;
		double rest=fbb_out+fcc_out+fll_out;
		double rest_new=1-fc_out; 
		fbb_out*=rest_new/rest;
		fcc_out*=rest_new/rest;
		fll_out*=rest_new/rest;
		if (ijet==1) cname="Wcjet1_down";
		if (ijet==3 || ijet==6) cname="Wcjet3_down";
		if (ijet==4 || ijet==7) cname="Wcjet4_down";
		if (ijet==5) cname="Wcjet5_down"; 
		GetFFactors_elec( -1, cname,  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		canorm=cafac[ijet];
	}

	fhf_out=fbb_out+fcc_out+fc_out;
	if (fabs(fll_out+fhf_out -1) > 0.001) cout << " FATAL this shoud be one " << fll_out+fhf_out << endl;
	if (fbb_out<0 || fcc_out<0 || fc_out<0 || fll_out<0) cout << " FATAL negative fraction " << fbb_out << " " << fcc_out << " " << fc_out << " " << fll_out << endl;

	if (idsys==2005) { // change c. at the cost of fbb+fcc+fl
		canorm*=1+sigma_canorm_up[ijet]; }
	else if (idsys==-2005) {
		canorm*=1+sigma_canorm_down[ijet];
	}
	if (abs(idsys)>2005) cout << " FATAL: unknown systematic check for HF factors " ;

	// calculate new Wjets_out (keeps normalization constant)
	Wjets_out[0]=fbb_out*nw_in;
	Wjets_out[1]=fcc_out*nw_in;
	Wjets_out[2]=fc_out*nw_in;
	Wjets_out[3]=fll_out*nw_in;

	if (_mode==1) {
		// mode for experts.
		// check PRETAGGED systematics long list
		// special mode to test long list systematics by defining: Wsys= HF x CA x Wnominal (default: HF x CA x Wsys)
		double nwsys=winpretag[ijet];
		double fbbsys=fbb[ijet];
		double fccsys=fcc[ijet];
		double fcsys=fc[ijet];
		double fllsys=fll[ijet];
		GetFFactors_elec( 0, "",  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		double nwnom=winpretag[ijet];	
		double fbbnom=fbb[ijet];
		double fccnom=fcc[ijet];
		double fcnom=fc[ijet];
		double fllnom=fll[ijet];
		canorm*=nwsys/nwnom;
		if (idsys<-1999 || idsys>1999) {
			fbbsys=fbb_out;
			fccsys=fcc_out;
			fcsys=fc_out;
			fllsys=fll_out;
		}
		Wjets_out[0]=fbbsys*nw_in;
		Wjets_out[1]=fccsys*nw_in;
		Wjets_out[2]=fcsys*nw_in;
		Wjets_out[3]=fllsys*nw_in;
		//  renormalize  to ensure that the overal normalization remains constant
		 nw_in=0;
		 nw_raw=0;
		for (int iw=0;iw<4;iw++) {
			nw_in+=Wjets_in[iw];
			nw_raw+=Wjets_out[iw];
		}
		if (nw_raw==0) { nw_in=0; nw_raw=1;} // some protection
		for (int iw=0;iw<4;iw++) {
			Wjets_out[iw]*=nw_in/nw_raw;
		}
	} // mode=1

	if (_mode==2) { 
		// mode for experts
		// check TAGGED systematics
		// special mode to test long list systematics by defining: Wsys= HF x CA x Wnominal (default: HF x CA x Wsys)
		double nwsys=wintag[ijet];
		double nwpresys=winpretag[ijet];
		double fbbsys=fbb[ijet];
		double fccsys=fcc[ijet];
		double fcsys=fc[ijet];
		double fllsys=fll[ijet];
		GetFFactors_elec( 0, "",  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		double nwnom=wintag[ijet];
		double nwprenom=winpretag[ijet];

		double fbbnom=fbb[ijet];
		double fccnom=fcc[ijet];
		double fcnom=fc[ijet];
		double fllnom=fll[ijet];
		canorm*=nwpresys/nwprenom/nwpresys * nwsys;
		//		canorm x winpretag_sys  --> canorm x winpretag_nom      x winpretag_sys/winpretag_nom      /nwpresys*nwsys
		if (idsys<-1999 || idsys>1999) {
			fbbsys=fbb_out;
			fccsys=fcc_out;
			fcsys=fc_out;
			fllsys=fll_out;
		}
		Wjets_out[0]=fbbsys*nw_in;
		Wjets_out[1]=fccsys*nw_in;
		Wjets_out[2]=fcsys*nw_in;
		Wjets_out[3]=fllsys*nw_in;
		//  renormalize  to ensure that the overal normalization remains constant
		 nw_in=0;
		 nw_raw=0;
		for (int iw=0;iw<4;iw++) {
			nw_in+=Wjets_in[iw];
			nw_raw+=Wjets_out[iw];
		}
		if (nw_raw==0) { nw_in=0; nw_raw=1;} // some protection
		for (int iw=0;iw<4;iw++) {
			Wjets_out[iw]*=nw_in/nw_raw;
		}
	} //mode=2
}





void GetFFactors_elec(int idsys, TString sysname, double Fbb[9],double Fcc[9],double Fc[9], double Fll[9], double cafac[9],
	double fbb[9],double fcc[9],double fc[9], double fll[9],
	double winpretag[9],double wintag[9],double woutpretag[9], double wouttag[9]) {
		//--
		// This function should NOT be called directly by the user. It is called via SetWflavors
		//
		// Ffactors for r17. ttbar cuts april 2012.  MUON
		// the Ffactors Fbb, Fcc, Fc, Fll contain weight factors for Wbb events, Wcc events, Wc event and Wlight events. 
		// These factors do NOT contains the normalization of MC Wjets events to data.
		//--
		// input: int idsys --> the id of the systematic variation, then call with sysname="" 
		//        TString sysname --> the name of the systematic variation, then call with idsys=-1 (NOT: idsys=0, because that is nominal).
		// output: Fbb[9],Fcc[9],Fc[9],Fll[9]
		//         position Fxx[0] is not used, position 1,2,3,4 are for the jetbins 1,2,3,4exclusive.
		//         position Fxx[5,6,7] for 5,3,4inclusive respectively. 
		//         HOWEVER, ONLY the 2jetbin should be used in all jetbins while keeping the pretagged normalization constant.
		//         cafac[9] --> the normalization obtained from Charge Assymmetry to be used in each jetbin (same for tagged and pretagged)
		// typical usage case:
		// You call this function once for nominal and then again for each systematic variation under study.
		// The input is either the integer idsys or the TString sysname. For the first mode with idsys, you typically may use:
		//      .....
		//      int idsys; 
		//      double Fbb[9], Fcc[9], Fc[9],  Fll[9];
		//      idsys=0; // for nominal
		//      GetFFactors_elec(idsys, "" , Fbb,Fcc,Fc,Fll); 
		//      .....
		//--
		// If you want to read this function by 'eye', i have two hints:
		// The HF factors are Fbb[2], Fcc[2] , Fc[2], Fll[2] for EACH jetbin, defined on the PRETAGGED sample.
		// The normalization factors are cafac[1], cafac[2], cafac[3], cafac[4], cafac[5] for jetbins 1,2,3,4,5 resp.
		// Ignore the other quantities which are used for additional checks and studies.
		//--
		double fmcbb[9];double fmccc[9];double fmcc[9];double fmcll[9];	



      Fbb[1]= 1.13522;   Fcc[1]= 1.13522;   Fc[1]= 0.886292;   Fll[1]= 1.01156; 
         cafac[1]= 0.972142;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08049e+06;   wouttag[1]= 38391.7;  
       fbb[1]= 0.0135665;   fcc[1]= 0.0396548;   fc[1]= 0.122801;   fll[1]= 0.823978;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.880444;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18364.2;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12441;   Fcc[3]= 1.12441;   Fc[3]= 0.87785;   Fll[3]= 1.00192; 
         cafac[3]= 0.795526;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48498.9;   wouttag[3]= 5791.92;  
       fbb[3]= 0.0553975;   fcc[3]= 0.10794;   fc[3]= 0.139503;   fll[3]= 0.697159;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.1176;   Fcc[4]= 1.1176;   Fc[4]= 0.872539;   Fll[4]= 0.995863; 
         cafac[4]= 0.787955;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.1;   wouttag[4]= 1840.07;  
       fbb[4]= 0.0767186;   fcc[4]= 0.127934;   fc[4]= 0.128458;   fll[4]= 0.666889;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.1103;   Fcc[5]= 1.1103;   Fc[5]= 0.866836;   Fll[5]= 0.989354; 
         cafac[5]= 0.722298;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.68;   wouttag[5]= 562.777;  
       fbb[5]= 0.0942562;   fcc[5]= 0.149158;   fc[5]= 0.112274;   fll[5]= 0.644312;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12254;   Fcc[6]= 1.12254;   Fc[6]= 0.876395;   Fll[6]= 1.00026; 
         cafac[6]= 0.789205;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61736;   wouttag[6]= 8209.2;  
       fbb[6]= 0.0609791;   fcc[6]= 0.113406;   fc[6]= 0.136268;   fll[6]= 0.689347;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.11601;   Fcc[7]= 1.11601;   Fc[7]= 0.871295;   Fll[7]= 0.994443; 
         cafac[7]= 0.771294;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13313.2;   wouttag[7]= 2402.64;  
       fbb[7]= 0.0805456;   fcc[7]= 0.132566;   fc[7]= 0.124926;   fll[7]= 0.661962;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 if ( idsys==0   || sysname=="norminal" ) { 
// OK
 }  
 else if ( idsys==1   || sysname=="tt_up" ) { 
      Fbb[1]= 1.13298;   Fcc[1]= 1.13298;   Fc[1]= 0.830933;   Fll[1]= 1.0211; 
         cafac[1]= 0.962193;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.06943e+06;   wouttag[1]= 36839.6;  
       fbb[1]= 0.0135397;   fcc[1]= 0.0395767;   fc[1]= 0.115131;   fll[1]= 0.831753;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.13053;   Fcc[2]= 1.13053;   Fc[2]= 0.829135;   Fll[2]= 1.0189; 
         cafac[2]= 0.870039;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 234107;   wouttag[2]= 17792.2;  
       fbb[2]= 0.0331329;   fcc[2]= 0.0799364;   fc[2]= 0.131347;   fll[2]= 0.755584;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12495;   Fcc[3]= 1.12495;   Fc[3]= 0.825046;   Fll[3]= 1.01387; 
         cafac[3]= 0.785743;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 47902.5;   wouttag[3]= 5644.31;  
       fbb[3]= 0.0554245;   fcc[3]= 0.107993;   fc[3]= 0.131112;   fll[3]= 0.705471;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.11779;   Fcc[4]= 1.11779;   Fc[4]= 0.819789;   Fll[4]= 1.00741; 
         cafac[4]= 0.779135;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10528.9;   wouttag[4]= 1804.47;  
       fbb[4]= 0.0767311;   fcc[4]= 0.127955;   fc[4]= 0.120692;   fll[4]= 0.674622;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.10968;   Fcc[5]= 1.10968;   Fc[5]= 0.813842;   Fll[5]= 1.0001; 
         cafac[5]= 0.714488;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2677.42;   wouttag[5]= 553.04;  
       fbb[5]= 0.0942035;   fcc[5]= 0.149074;   fc[5]= 0.10541;   fll[5]= 0.651312;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12297;   Fcc[6]= 1.12297;   Fc[6]= 0.823591;   Fll[6]= 1.01208; 
         cafac[6]= 0.779724;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 60994.3;   wouttag[6]= 8016.16;  
       fbb[6]= 0.0610023;   fcc[6]= 0.113449;   fc[6]= 0.128058;   fll[6]= 0.697491;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.11602;   Fcc[7]= 1.11602;   Fc[7]= 0.818491;   Fll[7]= 1.00581; 
         cafac[7]= 0.762716;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13165.1;   wouttag[7]= 2357.39;  
       fbb[7]= 0.0805461;   fcc[7]= 0.132566;   fc[7]= 0.117355;   fll[7]= 0.669532;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==2   || sysname=="tt_down" ) { 
      Fbb[1]= 1.13757;   Fcc[1]= 1.13757;   Fc[1]= 0.944359;   Fll[1]= 1.00155; 
         cafac[1]= 0.982802;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.09234e+06;   wouttag[1]= 40054.6;  
       fbb[1]= 0.0135945;   fcc[1]= 0.0397368;   fc[1]= 0.130847;   fll[1]= 0.815822;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.13068;   Fcc[2]= 1.13068;   Fc[2]= 0.938641;   Fll[2]= 0.995483; 
         cafac[2]= 0.891585;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 239904;   wouttag[2]= 18976.6;  
       fbb[2]= 0.0331372;   fcc[2]= 0.0799469;   fc[2]= 0.148694;   fll[2]= 0.738222;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12384;   Fcc[3]= 1.12384;   Fc[3]= 0.93296;   Fll[3]= 0.989458; 
         cafac[3]= 0.805998;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 49137.4;   wouttag[3]= 5949.95;  
       fbb[3]= 0.0553694;   fcc[3]= 0.107886;   fc[3]= 0.148261;   fll[3]= 0.688484;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.11741;   Fcc[4]= 1.11741;   Fc[4]= 0.927629;   Fll[4]= 0.983804; 
         cafac[4]= 0.797382;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10775.4;   wouttag[4]= 1878.12;  
       fbb[4]= 0.0767056;   fcc[4]= 0.127913;   fc[4]= 0.136568;   fll[4]= 0.658813;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.11095;   Fcc[5]= 1.11095;   Fc[5]= 0.922262;   Fll[5]= 0.978112; 
         cafac[5]= 0.730651;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2737.99;   wouttag[5]= 573.191;  
       fbb[5]= 0.0943114;   fcc[5]= 0.149245;   fc[5]= 0.119453;   fll[5]= 0.636991;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.1221;   Fcc[6]= 1.1221;   Fc[6]= 0.931518;   Fll[6]= 0.987928; 
         cafac[6]= 0.799351;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 62529.7;   wouttag[6]= 8415.79;  
       fbb[6]= 0.0609549;   fcc[6]= 0.113361;   fc[6]= 0.144839;   fll[6]= 0.680845;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.116;   Fcc[7]= 1.116;   Fc[7]= 0.926459;   Fll[7]= 0.982563; 
         cafac[7]= 0.780463;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13471.5;   wouttag[7]= 2451.01;  
       fbb[7]= 0.0805452;   fcc[7]= 0.132565;   fc[7]= 0.132836;   fll[7]= 0.654054;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==19   || sysname=="Kbb_up" ) { 
      Fbb[1]= 1.27993;   Fcc[1]= 1.27993;   Fc[1]= 0.879941;   Fll[1]= 1.00431; 
         cafac[1]= 0.965176;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08049e+06;   wouttag[1]= 38985;  
       fbb[1]= 0.0152958;   fcc[1]= 0.0447099;   fc[1]= 0.121921;   fll[1]= 0.818073;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.26454;   Fcc[2]= 1.26454;   Fc[2]= 0.869357;   Fll[2]= 0.99223; 
         cafac[2]= 0.867147;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18940.7;  
       fbb[2]= 0.0370602;   fcc[2]= 0.0894116;   fc[2]= 0.137718;   fll[2]= 0.73581;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.24922;   Fcc[3]= 1.24922;   Fc[3]= 0.858827;   Fll[3]= 0.980212; 
         cafac[3]= 0.778286;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48498.9;   wouttag[3]= 6016.54;  
       fbb[3]= 0.0615469;   fcc[3]= 0.119922;   fc[3]= 0.13648;   fll[3]= 0.682051;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.23489;   Fcc[4]= 1.23489;   Fc[4]= 0.848977;   Fll[4]= 0.968971; 
         cafac[4]= 0.766677;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.1;   wouttag[4]= 1909.72;  
       fbb[4]= 0.08477;   fcc[4]= 0.141361;   fc[4]= 0.124989;   fll[4]= 0.64888;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.22058;   Fcc[5]= 1.22058;   Fc[5]= 0.839136;   Fll[5]= 0.957739; 
         cafac[5]= 0.699217;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.68;   wouttag[5]= 585.917;  
       fbb[5]= 0.103618;   fcc[5]= 0.163973;   fc[5]= 0.108686;   fll[5]= 0.623723;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.24532;   Fcc[6]= 1.24532;   Fc[6]= 0.856149;   Fll[6]= 0.977156; 
         cafac[6]= 0.770972;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61736;   wouttag[6]= 8532.2;  
       fbb[6]= 0.0676488;   fcc[6]= 0.12581;   fc[6]= 0.13312;   fll[6]= 0.673421;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.23176;   Fcc[7]= 1.23176;   Fc[7]= 0.846821;   Fll[7]= 0.96651; 
         cafac[7]= 0.749629;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13313.2;   wouttag[7]= 2495.98;  
       fbb[7]= 0.0888994;   fcc[7]= 0.146315;   fc[7]= 0.121417;   fll[7]= 0.643369;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==20   || sysname=="Kbb_down" ) { 
      Fbb[1]= 0.988402;   Fcc[1]= 0.988402;   Fc[1]= 0.892735;   Fll[1]= 1.01891; 
         cafac[1]= 0.97921;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08049e+06;   wouttag[1]= 37789.8;  
       fbb[1]= 0.0118119;   fcc[1]= 0.0345263;   fc[1]= 0.123694;   fll[1]= 0.829968;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 0.992498;   Fcc[2]= 0.992498;   Fc[2]= 0.896434;   Fll[2]= 1.02314; 
         cafac[2]= 0.894156;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 17769.8;  
       fbb[2]= 0.0290875;   fcc[2]= 0.0701766;   fc[2]= 0.142008;   fll[2]= 0.758728;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 0.993939;   Fcc[3]= 0.993939;   Fc[3]= 0.897736;   Fll[3]= 1.02462; 
         cafac[3]= 0.813546;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48498.9;   wouttag[3]= 5557.13;  
       fbb[3]= 0.0489696;   fcc[3]= 0.0954158;   fc[3]= 0.142663;   fll[3]= 0.712951;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 0.993619;   Fcc[4]= 0.993619;   Fc[4]= 0.897447;   Fll[4]= 1.02429; 
         cafac[4]= 0.810448;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.1;   wouttag[4]= 1766.45;  
       fbb[4]= 0.0682076;   fcc[4]= 0.113742;   fc[4]= 0.132125;   fll[4]= 0.685926;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 0.99249;   Fcc[5]= 0.99249;   Fc[5]= 0.896427;   Fll[5]= 1.02313; 
         cafac[5]= 0.746955;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.68;   wouttag[5]= 538.058;  
       fbb[5]= 0.0842551;   fcc[5]= 0.133331;   fc[5]= 0.116107;   fll[5]= 0.666307;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 0.993814;   Fcc[6]= 0.993814;   Fc[6]= 0.897623;   Fll[6]= 1.02449; 
         cafac[6]= 0.808321;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61736;   wouttag[6]= 7870.57;  
       fbb[6]= 0.0539862;   fcc[6]= 0.100401;   fc[6]= 0.139569;   fll[6]= 0.706044;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 0.993373;   Fcc[7]= 0.993373;   Fc[7]= 0.897225;   Fll[7]= 1.02404; 
         cafac[7]= 0.794248;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13313.2;   wouttag[7]= 2303.74;  
       fbb[7]= 0.0716946;   fcc[7]= 0.117998;   fc[7]= 0.128644;   fll[7]= 0.681663;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==21   || sysname=="Kc_up" ) { 
      Fbb[1]= 1.12343;   Fcc[1]= 1.12343;   Fc[1]= 0.952019;   Fll[1]= 1.00106; 
         cafac[1]= 0.96205;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08049e+06;   wouttag[1]= 39709.7;  
       fbb[1]= 0.0134256;   fcc[1]= 0.0392432;   fc[1]= 0.131908;   fll[1]= 0.815423;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.11726;   Fcc[2]= 1.11726;   Fc[2]= 0.946784;   Fll[2]= 0.995554; 
         cafac[2]= 0.870051;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18727.3;  
       fbb[2]= 0.0327439;   fcc[2]= 0.0789979;   fc[2]= 0.149984;   fll[2]= 0.738274;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.11116;   Fcc[3]= 1.11116;   Fc[3]= 0.941622;   Fll[3]= 0.990125; 
         cafac[3]= 0.786156;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48498.9;   wouttag[3]= 5859.83;  
       fbb[3]= 0.0547451;   fcc[3]= 0.106669;   fc[3]= 0.149637;   fll[3]= 0.688948;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.10547;   Fcc[4]= 1.10547;   Fc[4]= 0.936798;   Fll[4]= 0.985053; 
         cafac[4]= 0.779402;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.1;   wouttag[4]= 1850.79;  
       fbb[4]= 0.0758859;   fcc[4]= 0.126546;   fc[4]= 0.137918;   fll[4]= 0.65965;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.09975;   Fcc[5]= 1.09975;   Fc[5]= 0.931949;   Fll[5]= 0.979955; 
         cafac[5]= 0.715436;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.68;   wouttag[5]= 564.75;  
       fbb[5]= 0.0933608;   fcc[5]= 0.147741;   fc[5]= 0.120708;   fll[5]= 0.638191;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.10963;   Fcc[6]= 1.10963;   Fc[6]= 0.940318;   Fll[6]= 0.988754; 
         cafac[6]= 0.780124;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61736;   wouttag[6]= 8288.52;  
       fbb[6]= 0.0602774;   fcc[6]= 0.112101;   fc[6]= 0.146207;   fll[6]= 0.681415;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.10423;   Fcc[7]= 1.10423;   Fc[7]= 0.935741;   Fll[7]= 0.983942; 
         cafac[7]= 0.763149;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13313.2;   wouttag[7]= 2415.12;  
       fbb[7]= 0.0796951;   fcc[7]= 0.131166;   fc[7]= 0.134167;   fll[7]= 0.654972;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==22   || sysname=="Kc_down" ) { 
      Fbb[1]= 1.14725;   Fcc[1]= 1.14725;   Fc[1]= 0.819171;   Fll[1]= 1.02228; 
         cafac[1]= 0.982449;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08049e+06;   wouttag[1]= 37045.7;  
       fbb[1]= 0.0137103;   fcc[1]= 0.0400753;   fc[1]= 0.113501;   fll[1]= 0.832713;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.14427;   Fcc[2]= 1.14427;   Fc[2]= 0.817041;   Fll[2]= 1.01963; 
         cafac[2]= 0.891089;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 17992.3;  
       fbb[2]= 0.0335356;   fcc[2]= 0.080908;   fc[2]= 0.129431;   fll[2]= 0.756126;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.13797;   Fcc[3]= 1.13797;   Fc[3]= 0.812541;   Fll[3]= 1.01401; 
         cafac[3]= 0.805121;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48498.9;   wouttag[3]= 5722.37;  
       fbb[3]= 0.0560657;   fcc[3]= 0.109242;   fc[3]= 0.129125;   fll[3]= 0.705567;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.13;   Fcc[4]= 1.13;   Fc[4]= 0.806855;   Fll[4]= 1.00691; 
         cafac[4]= 0.796698;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.1;   wouttag[4]= 1829.12;  
       fbb[4]= 0.0775699;   fcc[4]= 0.129354;   fc[4]= 0.118788;   fll[4]= 0.674289;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.12105;   Fcc[5]= 1.12105;   Fc[5]= 0.800461;   Fll[5]= 0.998935; 
         cafac[5]= 0.729293;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.68;   wouttag[5]= 560.766;  
       fbb[5]= 0.095169;   fcc[5]= 0.150602;   fc[5]= 0.103677;   fll[5]= 0.650552;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.13576;   Fcc[6]= 1.13576;   Fc[6]= 0.810967;   Fll[6]= 1.01205; 
         cafac[6]= 0.7985;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61736;   wouttag[6]= 8128.02;  
       fbb[6]= 0.0616973;   fcc[6]= 0.114742;   fc[6]= 0.126095;   fll[6]= 0.697466;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.12805;   Fcc[7]= 1.12805;   Fc[7]= 0.805458;   Fll[7]= 1.00517; 
         cafac[7]= 0.779614;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13313.2;   wouttag[7]= 2389.88;  
       fbb[7]= 0.0814145;   fcc[7]= 0.133996;   fc[7]= 0.115487;   fll[7]= 0.669103;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==23   || sysname=="Kll_up" ) { 
      Fbb[1]= 1.12779;   Fcc[1]= 1.12779;   Fc[1]= 0.880494;   Fll[1]= 1.01297; 
         cafac[1]= 0.973501;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08049e+06;   wouttag[1]= 38232.9;  
       fbb[1]= 0.0134777;   fcc[1]= 0.0393954;   fc[1]= 0.121998;   fll[1]= 0.825129;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.12389;   Fcc[2]= 1.12389;   Fc[2]= 0.877449;   Fll[2]= 1.00947; 
         cafac[2]= 0.882213;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18295.7;  
       fbb[2]= 0.0329383;   fcc[2]= 0.0794671;   fc[2]= 0.139;   fll[2]= 0.748594;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.11818;   Fcc[3]= 1.11818;   Fc[3]= 0.872987;   Fll[3]= 1.00434; 
         cafac[3]= 0.79744;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48498.9;   wouttag[3]= 5772.07;  
       fbb[3]= 0.0550906;   fcc[3]= 0.107342;   fc[3]= 0.13873;   fll[3]= 0.698837;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.11168;   Fcc[4]= 1.11168;   Fc[4]= 0.867914;   Fll[4]= 0.9985; 
         cafac[4]= 0.790042;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.1;   wouttag[4]= 1834.87;  
       fbb[4]= 0.0763119;   fcc[4]= 0.127256;   fc[4]= 0.127777;   fll[4]= 0.668655;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.10461;   Fcc[5]= 1.10461;   Fc[5]= 0.862395;   Fll[5]= 0.992152; 
         cafac[5]= 0.724341;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.68;   wouttag[5]= 561.19;  
       fbb[5]= 0.0937734;   fcc[5]= 0.148394;   fc[5]= 0.111699;   fll[5]= 0.646134;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.11639;   Fcc[6]= 1.11639;   Fc[6]= 0.871594;   Fll[6]= 1.00273; 
         cafac[6]= 0.791154;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61736;   wouttag[6]= 8182.36;  
       fbb[6]= 0.060645;   fcc[6]= 0.112785;   fc[6]= 0.135521;   fll[6]= 0.691049;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.11014;   Fcc[7]= 1.11014;   Fc[7]= 0.86671;   Fll[7]= 0.997115; 
         cafac[7]= 0.773366;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13313.2;   wouttag[7]= 2395.83;  
       fbb[7]= 0.0801217;   fcc[7]= 0.131868;   fc[7]= 0.124269;   fll[7]= 0.663741;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==24   || sysname=="Kll_down" ) { 
      Fbb[1]= 1.14274;   Fcc[1]= 1.14274;   Fc[1]= 0.892167;   Fll[1]= 1.01013; 
         cafac[1]= 0.970766;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08049e+06;   wouttag[1]= 38552.5;  
       fbb[1]= 0.0136564;   fcc[1]= 0.0399177;   fc[1]= 0.123615;   fll[1]= 0.822811;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.13739;   Fcc[2]= 1.13739;   Fc[2]= 0.88799;   Fll[2]= 1.0054; 
         cafac[2]= 0.878654;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18433.5;  
       fbb[2]= 0.033334;   fcc[2]= 0.0804217;   fc[2]= 0.14067;   fll[2]= 0.745574;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.13071;   Fcc[3]= 1.13071;   Fc[3]= 0.882769;   Fll[3]= 0.999486; 
         cafac[3]= 0.793589;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48498.9;   wouttag[3]= 5811.99;  
       fbb[3]= 0.0557079;   fcc[3]= 0.108545;   fc[3]= 0.140285;   fll[3]= 0.695462;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.12359;   Fcc[4]= 1.12359;   Fc[4]= 0.877215;   Fll[4]= 0.993198; 
         cafac[4]= 0.785846;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.1;   wouttag[4]= 1845.33;  
       fbb[4]= 0.0771297;   fcc[4]= 0.12862;   fc[4]= 0.129146;   fll[4]= 0.665104;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.11605;   Fcc[5]= 1.11605;   Fc[5]= 0.871323;   Fll[5]= 0.986527; 
         cafac[5]= 0.720234;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.68;   wouttag[5]= 564.38;  
       fbb[5]= 0.0947441;   fcc[5]= 0.14993;   fc[5]= 0.112855;   fll[5]= 0.642471;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12876;   Fcc[6]= 1.12876;   Fc[6]= 0.88125;   Fll[6]= 0.997767; 
         cafac[6]= 0.787235;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61736;   wouttag[6]= 8236.35;  
       fbb[6]= 0.0613169;   fcc[6]= 0.114034;   fc[6]= 0.137023;   fll[6]= 0.687626;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.12194;   Fcc[7]= 1.12194;   Fc[7]= 0.875929;   Fll[7]= 0.991742; 
         cafac[7]= 0.769199;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13313.2;   wouttag[7]= 2409.51;  
       fbb[7]= 0.080974;   fcc[7]= 0.133271;   fc[7]= 0.125591;   fll[7]= 0.660165;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==33   || sysname=="WbbWccjet1_up" ) { 
      Fbb[1]= 1.13333;   Fcc[1]= 1.13333;   Fc[1]= 0.88482;   Fll[1]= 1.00988; 
         cafac[1]= 0.968404;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 42270.2 ;   woutpretag[1]= 1.07634e+06;   wouttag[1]= 39393.5;  
       fbb[1]= 0.0169299;   fcc[1]= 0.0494863;   fc[1]= 0.12109;   fll[1]= 0.812494;  
       fmcbb[1]= 0.0149382;   fmccc[1]= 0.0436643;   fmcc[1]= 0.136852;   fmcll[1]= 0.804545;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.880444;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18364.2;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12441;   Fcc[3]= 1.12441;   Fc[3]= 0.87785;   Fll[3]= 1.00192; 
         cafac[3]= 0.795526;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48498.9;   wouttag[3]= 5791.92;  
       fbb[3]= 0.0553975;   fcc[3]= 0.10794;   fc[3]= 0.139503;   fll[3]= 0.697159;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.1176;   Fcc[4]= 1.1176;   Fc[4]= 0.872539;   Fll[4]= 0.995863; 
         cafac[4]= 0.787955;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.1;   wouttag[4]= 1840.07;  
       fbb[4]= 0.0767186;   fcc[4]= 0.127934;   fc[4]= 0.128458;   fll[4]= 0.666889;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.1103;   Fcc[5]= 1.1103;   Fc[5]= 0.866836;   Fll[5]= 0.989354; 
         cafac[5]= 0.722298;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.68;   wouttag[5]= 562.777;  
       fbb[5]= 0.0942562;   fcc[5]= 0.149158;   fc[5]= 0.112274;   fll[5]= 0.644312;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12254;   Fcc[6]= 1.12254;   Fc[6]= 0.876395;   Fll[6]= 1.00026; 
         cafac[6]= 0.789205;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61736;   wouttag[6]= 8209.2;  
       fbb[6]= 0.0609791;   fcc[6]= 0.113406;   fc[6]= 0.136268;   fll[6]= 0.689347;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.11601;   Fcc[7]= 1.11601;   Fc[7]= 0.871295;   Fll[7]= 0.994443; 
         cafac[7]= 0.771294;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13313.2;   wouttag[7]= 2402.64;  
       fbb[7]= 0.0805456;   fcc[7]= 0.132566;   fc[7]= 0.124926;   fll[7]= 0.661962;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==34   || sysname=="WbbWccjet1_down" ) { 
      Fbb[1]= 1.13711;   Fcc[1]= 1.13711;   Fc[1]= 0.887768;   Fll[1]= 1.01324; 
         cafac[1]= 0.975922;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 40218.7 ;   woutpretag[1]= 1.08469e+06;   wouttag[1]= 37378.7;  
       fbb[1]= 0.0101918;   fcc[1]= 0.0297907;   fc[1]= 0.124518;   fll[1]= 0.835499;  
       fmcbb[1]= 0.00896289;   fmccc[1]= 0.0261986;   fmcc[1]= 0.14026;   fmcll[1]= 0.824579;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.880444;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18364.2;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12441;   Fcc[3]= 1.12441;   Fc[3]= 0.87785;   Fll[3]= 1.00192; 
         cafac[3]= 0.795526;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48498.9;   wouttag[3]= 5791.92;  
       fbb[3]= 0.0553975;   fcc[3]= 0.10794;   fc[3]= 0.139503;   fll[3]= 0.697159;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.1176;   Fcc[4]= 1.1176;   Fc[4]= 0.872539;   Fll[4]= 0.995863; 
         cafac[4]= 0.787955;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.1;   wouttag[4]= 1840.07;  
       fbb[4]= 0.0767186;   fcc[4]= 0.127934;   fc[4]= 0.128458;   fll[4]= 0.666889;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.1103;   Fcc[5]= 1.1103;   Fc[5]= 0.866836;   Fll[5]= 0.989354; 
         cafac[5]= 0.722298;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.68;   wouttag[5]= 562.777;  
       fbb[5]= 0.0942562;   fcc[5]= 0.149158;   fc[5]= 0.112274;   fll[5]= 0.644312;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12254;   Fcc[6]= 1.12254;   Fc[6]= 0.876395;   Fll[6]= 1.00026; 
         cafac[6]= 0.789205;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61736;   wouttag[6]= 8209.2;  
       fbb[6]= 0.0609791;   fcc[6]= 0.113406;   fc[6]= 0.136268;   fll[6]= 0.689347;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.11601;   Fcc[7]= 1.11601;   Fc[7]= 0.871295;   Fll[7]= 0.994443; 
         cafac[7]= 0.771294;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13313.2;   wouttag[7]= 2402.64;  
       fbb[7]= 0.0805456;   fcc[7]= 0.132566;   fc[7]= 0.124926;   fll[7]= 0.661962;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==35   || sysname=="WbbWccjet3_up" ) { 
      Fbb[1]= 1.13522;   Fcc[1]= 1.13522;   Fc[1]= 0.886292;   Fll[1]= 1.01156; 
         cafac[1]= 0.972142;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08049e+06;   wouttag[1]= 38391.7;  
       fbb[1]= 0.0135665;   fcc[1]= 0.0396548;   fc[1]= 0.122801;   fll[1]= 0.823978;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.880444;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18364.2;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.11849;   Fcc[3]= 1.11849;   Fc[3]= 0.873235;   Fll[3]= 0.996657; 
         cafac[3]= 0.789672;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7749.26 ;   woutpretag[3]= 48142.1;   wouttag[3]= 6238.26;  
       fbb[3]= 0.0688828;   fcc[3]= 0.134216;   fc[3]= 0.132874;   fll[3]= 0.664027;  
       fmcbb[3]= 0.0615853;   fmccc[3]= 0.119997;   fmcc[3]= 0.152163;   fmcll[3]= 0.666255;  
      Fbb[4]= 1.1176;   Fcc[4]= 1.1176;   Fc[4]= 0.872539;   Fll[4]= 0.995863; 
         cafac[4]= 0.787955;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.1;   wouttag[4]= 1840.07;  
       fbb[4]= 0.0767186;   fcc[4]= 0.127934;   fc[4]= 0.128458;   fll[4]= 0.666889;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.1103;   Fcc[5]= 1.1103;   Fc[5]= 0.866836;   Fll[5]= 0.989354; 
         cafac[5]= 0.722298;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.68;   wouttag[5]= 562.777;  
       fbb[5]= 0.0942562;   fcc[5]= 0.149158;   fc[5]= 0.112274;   fll[5]= 0.644312;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.11794;   Fcc[6]= 1.11794;   Fc[6]= 0.872806;   Fll[6]= 0.996167; 
         cafac[6]= 0.784995;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10791 ;   woutpretag[6]= 61406.6;   wouttag[6]= 8647.82;  
       fbb[6]= 0.0714607;   fcc[6]= 0.133851;   fc[6]= 0.131117;   fll[6]= 0.663571;  
       fmcbb[6]= 0.0639215;   fmccc[6]= 0.11973;   fmcc[6]= 0.150225;   fmcll[6]= 0.666124;  
      Fbb[7]= 1.11601;   Fcc[7]= 1.11601;   Fc[7]= 0.871295;   Fll[7]= 0.994443; 
         cafac[7]= 0.771294;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13313.2;   wouttag[7]= 2402.64;  
       fbb[7]= 0.0805456;   fcc[7]= 0.132566;   fc[7]= 0.124926;   fll[7]= 0.661962;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==36   || sysname=="WbbWccjet3_down" ) { 
      Fbb[1]= 1.13522;   Fcc[1]= 1.13522;   Fc[1]= 0.886292;   Fll[1]= 1.01156; 
         cafac[1]= 0.972142;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08049e+06;   wouttag[1]= 38391.7;  
       fbb[1]= 0.0135665;   fcc[1]= 0.0396548;   fc[1]= 0.122801;   fll[1]= 0.823978;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.880444;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18364.2;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.13038;   Fcc[3]= 1.13038;   Fc[3]= 0.882515;   Fll[3]= 1.00725; 
         cafac[3]= 0.80153;  
         winpretag[3]= 60964.7 ;   wintag[3]= 6634.54 ;   woutpretag[3]= 48865;   wouttag[3]= 5334.07;  
       fbb[3]= 0.0417689;   fcc[3]= 0.0813855;   fc[3]= 0.146203;   fll[3]= 0.730642;  
       fmcbb[3]= 0.0369512;   fmccc[3]= 0.0719983;   fmcc[3]= 0.165667;   fmcll[3]= 0.725384;  
      Fbb[4]= 1.1176;   Fcc[4]= 1.1176;   Fc[4]= 0.872539;   Fll[4]= 0.995863; 
         cafac[4]= 0.787955;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.1;   wouttag[4]= 1840.07;  
       fbb[4]= 0.0767186;   fcc[4]= 0.127934;   fc[4]= 0.128458;   fll[4]= 0.666889;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.1103;   Fcc[5]= 1.1103;   Fc[5]= 0.866836;   Fll[5]= 0.989354; 
         cafac[5]= 0.722298;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.68;   wouttag[5]= 562.777;  
       fbb[5]= 0.0942562;   fcc[5]= 0.149158;   fc[5]= 0.112274;   fll[5]= 0.644312;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12718;   Fcc[6]= 1.12718;   Fc[6]= 0.880015;   Fll[6]= 1.00439; 
         cafac[6]= 0.793496;  
         winpretag[6]= 78225.5 ;   wintag[6]= 9676.26 ;   woutpretag[6]= 62071.6;   wouttag[6]= 7762.17;  
       fbb[6]= 0.0504108;   fcc[6]= 0.0927917;   fc[6]= 0.141462;   fll[6]= 0.715336;  
       fmcbb[6]= 0.044723;   fmccc[6]= 0.0823221;   fmcc[6]= 0.160749;   fmcll[6]= 0.712206;  
      Fbb[7]= 1.11601;   Fcc[7]= 1.11601;   Fc[7]= 0.871295;   Fll[7]= 0.994443; 
         cafac[7]= 0.771294;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13313.2;   wouttag[7]= 2402.64;  
       fbb[7]= 0.0805456;   fcc[7]= 0.132566;   fc[7]= 0.124926;   fll[7]= 0.661962;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==37   || sysname=="WbbWccjet4_up" ) { 
      Fbb[1]= 1.13522;   Fcc[1]= 1.13522;   Fc[1]= 0.886292;   Fll[1]= 1.01156; 
         cafac[1]= 0.972142;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08049e+06;   wouttag[1]= 38391.7;  
       fbb[1]= 0.0135665;   fcc[1]= 0.0396548;   fc[1]= 0.122801;   fll[1]= 0.823978;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.880444;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18364.2;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12441;   Fcc[3]= 1.12441;   Fc[3]= 0.87785;   Fll[3]= 1.00192; 
         cafac[3]= 0.795526;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48498.9;   wouttag[3]= 5791.92;  
       fbb[3]= 0.0553975;   fcc[3]= 0.10794;   fc[3]= 0.139503;   fll[3]= 0.697159;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.11029;   Fcc[4]= 1.11029;   Fc[4]= 0.866826;   Fll[4]= 0.989343; 
         cafac[4]= 0.78338;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2471.09 ;   woutpretag[4]= 10586.2;   wouttag[4]= 1988.93;  
       fbb[4]= 0.0952704;   fcc[4]= 0.158871;   fc[4]= 0.120465;   fll[4]= 0.625394;  
       fmcbb[4]= 0.085807;   fmccc[4]= 0.14309;   fmcc[4]= 0.138972;   fmcll[4]= 0.63213;  
      Fbb[5]= 1.1103;   Fcc[5]= 1.1103;   Fc[5]= 0.866836;   Fll[5]= 0.989354; 
         cafac[5]= 0.722298;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.68;   wouttag[5]= 562.777;  
       fbb[5]= 0.0942562;   fcc[5]= 0.149158;   fc[5]= 0.112274;   fll[5]= 0.644312;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12126;   Fcc[6]= 1.12126;   Fc[6]= 0.875394;   Fll[6]= 0.999122; 
         cafac[6]= 0.788239;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10419.9 ;   woutpretag[6]= 61660.4;   wouttag[6]= 8364.08;  
       fbb[6]= 0.0642336;   fcc[6]= 0.11882;   fc[6]= 0.134865;   fll[6]= 0.682082;  
       fmcbb[6]= 0.0572869;   fmccc[6]= 0.10597;   fmcc[6]= 0.154062;   fmcll[6]= 0.682682;  
      Fbb[7]= 1.11029;   Fcc[7]= 1.11029;   Fc[7]= 0.866829;   Fll[7]= 0.989345; 
         cafac[7]= 0.768107;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3227.98 ;   woutpretag[7]= 13258.2;   wouttag[7]= 2548.63;  
       fbb[7]= 0.0950502;   fcc[7]= 0.156762;   fc[7]= 0.118687;   fll[7]= 0.629501;  
       fmcbb[7]= 0.0856085;   fmccc[7]= 0.141191;   fmcc[7]= 0.136921;   fmcll[7]= 0.63628;  
 }  
 else if ( idsys==38   || sysname=="WbbWccjet4_down" ) { 
      Fbb[1]= 1.13522;   Fcc[1]= 1.13522;   Fc[1]= 0.886292;   Fll[1]= 1.01156; 
         cafac[1]= 0.972142;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08049e+06;   wouttag[1]= 38391.7;  
       fbb[1]= 0.0135665;   fcc[1]= 0.0396548;   fc[1]= 0.122801;   fll[1]= 0.823978;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.880444;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18364.2;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12441;   Fcc[3]= 1.12441;   Fc[3]= 0.87785;   Fll[3]= 1.00192; 
         cafac[3]= 0.795526;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48498.9;   wouttag[3]= 5791.92;  
       fbb[3]= 0.0553975;   fcc[3]= 0.10794;   fc[3]= 0.139503;   fll[3]= 0.697159;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.12502;   Fcc[4]= 1.12502;   Fc[4]= 0.878328;   Fll[4]= 1.00247; 
         cafac[4]= 0.792645;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2098.57 ;   woutpretag[4]= 10711.4;   wouttag[4]= 1687.45;  
       fbb[4]= 0.0579207;   fcc[4]= 0.0965875;   fc[4]= 0.136557;   fll[4]= 0.708935;  
       fmcbb[4]= 0.0514842;   fmccc[4]= 0.0858541;   fmcc[4]= 0.155474;   fmcll[4]= 0.707188;  
      Fbb[5]= 1.1103;   Fcc[5]= 1.1103;   Fc[5]= 0.866836;   Fll[5]= 0.989354; 
         cafac[5]= 0.722298;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.68;   wouttag[5]= 562.777;  
       fbb[5]= 0.0942562;   fcc[5]= 0.149158;   fc[5]= 0.112274;   fll[5]= 0.644312;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12383;   Fcc[6]= 1.12383;   Fc[6]= 0.877399;   Fll[6]= 1.00141; 
         cafac[6]= 0.790176;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10047.4 ;   woutpretag[6]= 61811.9;   wouttag[6]= 8053.59;  
       fbb[6]= 0.0577171;   fcc[6]= 0.10798;   fc[6]= 0.137675;   fll[6]= 0.696629;  
       fmcbb[6]= 0.0513576;   fmccc[6]= 0.0960821;   fmcc[6]= 0.156912;   fmcll[6]= 0.695648;  
      Fbb[7]= 1.12179;   Fcc[7]= 1.12179;   Fc[7]= 0.875807;   Fll[7]= 0.999593; 
         cafac[7]= 0.77454;  
         winpretag[7]= 17260.9 ;   wintag[7]= 2855.46 ;   woutpretag[7]= 13369.2;   wouttag[7]= 2253.9;  
       fbb[7]= 0.0658908;   fcc[7]= 0.108119;   fc[7]= 0.131231;   fll[7]= 0.69476;  
       fmcbb[7]= 0.0587372;   fmccc[7]= 0.0963804;   fmcc[7]= 0.149839;   fmcll[7]= 0.695043;  
 }  
 else if ( idsys==73   || sysname=="WbbWccjet5_up" ) { 
      Fbb[1]= 1.13522;   Fcc[1]= 1.13522;   Fc[1]= 0.886292;   Fll[1]= 1.01156; 
         cafac[1]= 0.972142;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08049e+06;   wouttag[1]= 38391.7;  
       fbb[1]= 0.0135665;   fcc[1]= 0.0396548;   fc[1]= 0.122801;   fll[1]= 0.823978;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.880444;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18364.2;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12441;   Fcc[3]= 1.12441;   Fc[3]= 0.87785;   Fll[3]= 1.00192; 
         cafac[3]= 0.795526;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48498.9;   wouttag[3]= 5791.92;  
       fbb[3]= 0.0553975;   fcc[3]= 0.10794;   fc[3]= 0.139503;   fll[3]= 0.697159;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.1176;   Fcc[4]= 1.1176;   Fc[4]= 0.872539;   Fll[4]= 0.995863; 
         cafac[4]= 0.787955;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.1;   wouttag[4]= 1840.07;  
       fbb[4]= 0.0767186;   fcc[4]= 0.127934;   fc[4]= 0.128458;   fll[4]= 0.666889;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.10177;   Fcc[5]= 1.10177;   Fc[5]= 0.860176;   Fll[5]= 0.981752; 
         cafac[5]= 0.715266;  
         winpretag[5]= 3747.32 ;   wintag[5]= 828.828 ;   woutpretag[5]= 2680.33;   wouttag[5]= 612.757;  
       fbb[5]= 0.116915;   fcc[5]= 0.185015;   fc[5]= 0.103591;   fll[5]= 0.59448;  
       fmcbb[5]= 0.106116;   fmccc[5]= 0.167925;   fmcc[5]= 0.12043;   fmcll[5]= 0.605529;  
      Fbb[6]= 1.12212;   Fcc[6]= 1.12212;   Fc[6]= 0.876067;   Fll[6]= 0.999889; 
         cafac[6]= 0.788607;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10305.6 ;   woutpretag[6]= 61689.2;   wouttag[6]= 8267;  
       fbb[6]= 0.0620971;   fcc[6]= 0.115169;   fc[6]= 0.135835;   fll[6]= 0.686899;  
       fmcbb[6]= 0.055339;   fmccc[6]= 0.102635;   fmcc[6]= 0.155051;   fmcll[6]= 0.686975;  
      Fbb[7]= 1.11413;   Fcc[7]= 1.11413;   Fc[7]= 0.869825;   Fll[7]= 0.992765; 
         cafac[7]= 0.769171;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3113.66 ;   woutpretag[7]= 13276.5;   wouttag[7]= 2456.95;  
       fbb[7]= 0.0855431;   fcc[7]= 0.140466;   fc[7]= 0.122999;   fll[7]= 0.650993;  
       fmcbb[7]= 0.0767804;   fmccc[7]= 0.126077;   fmcc[7]= 0.141406;   fmcll[7]= 0.655737;  
 }  
 else if ( idsys==74   || sysname=="WbbWccjet5_down" ) { 
      Fbb[1]= 1.13522;   Fcc[1]= 1.13522;   Fc[1]= 0.886292;   Fll[1]= 1.01156; 
         cafac[1]= 0.972142;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08049e+06;   wouttag[1]= 38391.7;  
       fbb[1]= 0.0135665;   fcc[1]= 0.0396548;   fc[1]= 0.122801;   fll[1]= 0.823978;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.880444;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18364.2;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12441;   Fcc[3]= 1.12441;   Fc[3]= 0.87785;   Fll[3]= 1.00192; 
         cafac[3]= 0.795526;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48498.9;   wouttag[3]= 5791.92;  
       fbb[3]= 0.0553975;   fcc[3]= 0.10794;   fc[3]= 0.139503;   fll[3]= 0.697159;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.1176;   Fcc[4]= 1.1176;   Fc[4]= 0.872539;   Fll[4]= 0.995863; 
         cafac[4]= 0.787955;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.1;   wouttag[4]= 1840.07;  
       fbb[4]= 0.0767186;   fcc[4]= 0.127934;   fc[4]= 0.128458;   fll[4]= 0.666889;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.11896;   Fcc[5]= 1.11896;   Fc[5]= 0.8736;   Fll[5]= 0.997074; 
         cafac[5]= 0.729583;  
         winpretag[5]= 3747.32 ;   wintag[5]= 684.954 ;   woutpretag[5]= 2733.98;   wouttag[5]= 511.001;  
       fbb[5]= 0.0712438;   fcc[5]= 0.112741;   fc[5]= 0.121093;   fll[5]= 0.694922;  
       fmcbb[5]= 0.0636695;   fmccc[5]= 0.100755;   fmcc[5]= 0.138614;   fmcll[5]= 0.696961;  
      Fbb[6]= 1.12296;   Fcc[6]= 1.12296;   Fc[6]= 0.876724;   Fll[6]= 1.00064; 
         cafac[6]= 0.789804;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10161.7 ;   woutpretag[6]= 61782.8;   wouttag[6]= 8151.28;  
       fbb[6]= 0.0598603;   fcc[6]= 0.111642;   fc[6]= 0.136701;   fll[6]= 0.691797;  
       fmcbb[6]= 0.0533056;   fmccc[6]= 0.0994171;   fmcc[6]= 0.155922;   fmcll[6]= 0.691355;  
      Fbb[7]= 1.1179;   Fcc[7]= 1.1179;   Fc[7]= 0.87277;   Fll[7]= 0.996126; 
         cafac[7]= 0.773436;  
         winpretag[7]= 17260.9 ;   wintag[7]= 2969.78 ;   woutpretag[7]= 13350.2;   wouttag[7]= 2347.84;  
       fbb[7]= 0.0755312;   fcc[7]= 0.124639;   fc[7]= 0.126861;   fll[7]= 0.672969;  
       fmcbb[7]= 0.0675653;   fmccc[7]= 0.111494;   fmcc[7]= 0.145354;   fmcll[7]= 0.675587;  
 }  
 else if ( idsys==41   || sysname=="Wcjet1_up" ) { 
      Fbb[1]= 1.14043;   Fcc[1]= 1.14043;   Fc[1]= 0.890363;   Fll[1]= 1.01621; 
         cafac[1]= 1.01541;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 46425.2 ;   woutpretag[1]= 1.12858e+06;   wouttag[1]= 44848;  
       fbb[1]= 0.0130808;   fcc[1]= 0.0382351;   fc[1]= 0.154206;   fll[1]= 0.794478;  
       fmcbb[1]= 0.01147;   fmccc[1]= 0.0335269;   fmcc[1]= 0.173195;   fmcll[1]= 0.781808;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.880444;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18364.2;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12441;   Fcc[3]= 1.12441;   Fc[3]= 0.87785;   Fll[3]= 1.00192; 
         cafac[3]= 0.795526;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48498.9;   wouttag[3]= 5791.92;  
       fbb[3]= 0.0553975;   fcc[3]= 0.10794;   fc[3]= 0.139503;   fll[3]= 0.697159;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.1176;   Fcc[4]= 1.1176;   Fc[4]= 0.872539;   Fll[4]= 0.995863; 
         cafac[4]= 0.787955;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.1;   wouttag[4]= 1840.07;  
       fbb[4]= 0.0767186;   fcc[4]= 0.127934;   fc[4]= 0.128458;   fll[4]= 0.666889;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.1103;   Fcc[5]= 1.1103;   Fc[5]= 0.866836;   Fll[5]= 0.989354; 
         cafac[5]= 0.722298;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.68;   wouttag[5]= 562.777;  
       fbb[5]= 0.0942562;   fcc[5]= 0.149158;   fc[5]= 0.112274;   fll[5]= 0.644312;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12254;   Fcc[6]= 1.12254;   Fc[6]= 0.876395;   Fll[6]= 1.00026; 
         cafac[6]= 0.789205;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61736;   wouttag[6]= 8209.2;  
       fbb[6]= 0.0609791;   fcc[6]= 0.113406;   fc[6]= 0.136268;   fll[6]= 0.689347;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.11601;   Fcc[7]= 1.11601;   Fc[7]= 0.871295;   Fll[7]= 0.994443; 
         cafac[7]= 0.771294;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13313.2;   wouttag[7]= 2402.64;  
       fbb[7]= 0.0805456;   fcc[7]= 0.132566;   fc[7]= 0.124926;   fll[7]= 0.661962;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==42   || sysname=="Wcjet1_down" ) { 
      Fbb[1]= 1.13005;   Fcc[1]= 1.13005;   Fc[1]= 0.882258;   Fll[1]= 1.00696; 
         cafac[1]= 0.932757;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 36063.7 ;   woutpretag[1]= 1.03672e+06;   wouttag[1]= 32514.9;  
       fbb[1]= 0.0140477;   fcc[1]= 0.0410616;   fc[1]= 0.0916816;   fll[1]= 0.853209;  
       fmcbb[1]= 0.0124311;   fmccc[1]= 0.0363361;   fmcc[1]= 0.103917;   fmcll[1]= 0.847316;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.880444;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18364.2;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12441;   Fcc[3]= 1.12441;   Fc[3]= 0.87785;   Fll[3]= 1.00192; 
         cafac[3]= 0.795526;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48498.9;   wouttag[3]= 5791.92;  
       fbb[3]= 0.0553975;   fcc[3]= 0.10794;   fc[3]= 0.139503;   fll[3]= 0.697159;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.1176;   Fcc[4]= 1.1176;   Fc[4]= 0.872539;   Fll[4]= 0.995863; 
         cafac[4]= 0.787955;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.1;   wouttag[4]= 1840.07;  
       fbb[4]= 0.0767186;   fcc[4]= 0.127934;   fc[4]= 0.128458;   fll[4]= 0.666889;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.1103;   Fcc[5]= 1.1103;   Fc[5]= 0.866836;   Fll[5]= 0.989354; 
         cafac[5]= 0.722298;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.68;   wouttag[5]= 562.777;  
       fbb[5]= 0.0942562;   fcc[5]= 0.149158;   fc[5]= 0.112274;   fll[5]= 0.644312;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12254;   Fcc[6]= 1.12254;   Fc[6]= 0.876395;   Fll[6]= 1.00026; 
         cafac[6]= 0.789205;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61736;   wouttag[6]= 8209.2;  
       fbb[6]= 0.0609791;   fcc[6]= 0.113406;   fc[6]= 0.136268;   fll[6]= 0.689347;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.11601;   Fcc[7]= 1.11601;   Fc[7]= 0.871295;   Fll[7]= 0.994443; 
         cafac[7]= 0.771294;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13313.2;   wouttag[7]= 2402.64;  
       fbb[7]= 0.0805456;   fcc[7]= 0.132566;   fc[7]= 0.124926;   fll[7]= 0.661962;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==43   || sysname=="Wcjet3_up" ) { 
      Fbb[1]= 1.13522;   Fcc[1]= 1.13522;   Fc[1]= 0.886292;   Fll[1]= 1.01156; 
         cafac[1]= 0.972142;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08049e+06;   wouttag[1]= 38391.7;  
       fbb[1]= 0.0135665;   fcc[1]= 0.0396548;   fc[1]= 0.122801;   fll[1]= 0.823978;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.880444;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18364.2;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.13093;   Fcc[3]= 1.13093;   Fc[3]= 0.882945;   Fll[3]= 1.00774; 
         cafac[3]= 0.839916;  
         winpretag[3]= 60964.6 ;   wintag[3]= 7538.49 ;   woutpretag[3]= 51205.2;   wouttag[3]= 6369.03;  
       fbb[3]= 0.0530871;   fcc[3]= 0.103439;   fc[3]= 0.175391;   fll[3]= 0.668083;  
       fmcbb[3]= 0.0469411;   fmccc[3]= 0.0914632;   fmcc[3]= 0.198643;   fmcll[3]= 0.662952;  
      Fbb[4]= 1.1176;   Fcc[4]= 1.1176;   Fc[4]= 0.872539;   Fll[4]= 0.995863; 
         cafac[4]= 0.787955;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.1;   wouttag[4]= 1840.07;  
       fbb[4]= 0.0767186;   fcc[4]= 0.127934;   fc[4]= 0.128458;   fll[4]= 0.666889;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.1103;   Fcc[5]= 1.1103;   Fc[5]= 0.866836;   Fll[5]= 0.989354; 
         cafac[5]= 0.722298;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.68;   wouttag[5]= 562.777;  
       fbb[5]= 0.0942562;   fcc[5]= 0.149158;   fc[5]= 0.112274;   fll[5]= 0.644312;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.1276;   Fcc[6]= 1.1276;   Fc[6]= 0.880347;   Fll[6]= 1.00477; 
         cafac[6]= 0.821378;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10580.2 ;   woutpretag[6]= 64252.7;   wouttag[6]= 8795.38;  
       fbb[6]= 0.0592089;   fcc[6]= 0.109932;   fc[6]= 0.16414;   fll[6]= 0.666718;  
       fmcbb[6]= 0.0525086;   fmccc[6]= 0.097492;   fmcc[6]= 0.186449;   fmcll[6]= 0.66355;  
      Fbb[7]= 1.11601;   Fcc[7]= 1.11601;   Fc[7]= 0.871295;   Fll[7]= 0.994443; 
         cafac[7]= 0.771294;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13313.2;   wouttag[7]= 2402.64;  
       fbb[7]= 0.0805456;   fcc[7]= 0.132566;   fc[7]= 0.124926;   fll[7]= 0.661962;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==44   || sysname=="Wcjet3_down" ) { 
      Fbb[1]= 1.13522;   Fcc[1]= 1.13522;   Fc[1]= 0.886292;   Fll[1]= 1.01156; 
         cafac[1]= 0.972142;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08049e+06;   wouttag[1]= 38391.7;  
       fbb[1]= 0.0135665;   fcc[1]= 0.0396548;   fc[1]= 0.122801;   fll[1]= 0.823978;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.880444;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18364.2;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.11796;   Fcc[3]= 1.11796;   Fc[3]= 0.872815;   Fll[3]= 0.996177; 
         cafac[3]= 0.756027;  
         winpretag[3]= 60964.7 ;   wintag[3]= 6845.31 ;   woutpretag[3]= 46090.9;   wouttag[3]= 5278.41;  
       fbb[3]= 0.0576814;   fcc[3]= 0.112391;   fc[3]= 0.104027;   fll[3]= 0.725901;  
       fmcbb[3]= 0.0515954;   fmccc[3]= 0.100532;   fmcc[3]= 0.119186;   fmcll[3]= 0.728686;  
      Fbb[4]= 1.1176;   Fcc[4]= 1.1176;   Fc[4]= 0.872539;   Fll[4]= 0.995863; 
         cafac[4]= 0.787955;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.1;   wouttag[4]= 1840.07;  
       fbb[4]= 0.0767186;   fcc[4]= 0.127934;   fc[4]= 0.128458;   fll[4]= 0.666889;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.1103;   Fcc[5]= 1.1103;   Fc[5]= 0.866836;   Fll[5]= 0.989354; 
         cafac[5]= 0.722298;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.68;   wouttag[5]= 562.777;  
       fbb[5]= 0.0942562;   fcc[5]= 0.149158;   fc[5]= 0.112274;   fll[5]= 0.644312;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.11753;   Fcc[6]= 1.11753;   Fc[6]= 0.872479;   Fll[6]= 0.995794; 
         cafac[6]= 0.759713;  
         winpretag[6]= 78225.5 ;   wintag[6]= 9887.03 ;   woutpretag[6]= 59429;   wouttag[6]= 7671.89;  
       fbb[6]= 0.0627334;   fcc[6]= 0.116848;   fc[6]= 0.108645;   fll[6]= 0.711773;  
       fmcbb[6]= 0.056136;   fmccc[6]= 0.10456;   fmcc[6]= 0.124525;   fmcll[6]= 0.71478;  
      Fbb[7]= 1.11601;   Fcc[7]= 1.11601;   Fc[7]= 0.871295;   Fll[7]= 0.994443; 
         cafac[7]= 0.771294;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13313.2;   wouttag[7]= 2402.64;  
       fbb[7]= 0.0805456;   fcc[7]= 0.132566;   fc[7]= 0.124926;   fll[7]= 0.661962;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==45   || sysname=="Wcjet4_up" ) { 
      Fbb[1]= 1.13522;   Fcc[1]= 1.13522;   Fc[1]= 0.886292;   Fll[1]= 1.01156; 
         cafac[1]= 0.972142;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08049e+06;   wouttag[1]= 38391.7;  
       fbb[1]= 0.0135665;   fcc[1]= 0.0396548;   fc[1]= 0.122801;   fll[1]= 0.823978;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.880444;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18364.2;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12441;   Fcc[3]= 1.12441;   Fc[3]= 0.87785;   Fll[3]= 1.00192; 
         cafac[3]= 0.795526;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48498.9;   wouttag[3]= 5791.92;  
       fbb[3]= 0.0553975;   fcc[3]= 0.10794;   fc[3]= 0.139503;   fll[3]= 0.697159;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.12379;   Fcc[4]= 1.12379;   Fc[4]= 0.877366;   Fll[4]= 1.00137; 
         cafac[4]= 0.826998;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2341.09 ;   woutpretag[4]= 11175.7;   wouttag[4]= 1970.49;  
       fbb[4]= 0.0738135;   fcc[4]= 0.12309;   fc[4]= 0.161461;   fll[4]= 0.641636;  
       fmcbb[4]= 0.0656829;   fmccc[4]= 0.109532;   fmcc[4]= 0.184029;   fmcll[4]= 0.640757;  
      Fbb[5]= 1.1103;   Fcc[5]= 1.1103;   Fc[5]= 0.866836;   Fll[5]= 0.989354; 
         cafac[5]= 0.722298;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.68;   wouttag[5]= 562.777;  
       fbb[5]= 0.0942562;   fcc[5]= 0.149158;   fc[5]= 0.112274;   fll[5]= 0.644312;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12362;   Fcc[6]= 1.12362;   Fc[6]= 0.877233;   Fll[6]= 1.00122; 
         cafac[6]= 0.796395;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10289.9 ;   woutpretag[6]= 62298.4;   wouttag[6]= 8319.41;  
       fbb[6]= 0.0604623;   fcc[6]= 0.112555;   fc[6]= 0.141976;   fll[6]= 0.685007;  
       fmcbb[6]= 0.0538105;   fmccc[6]= 0.100172;   fmcc[6]= 0.161845;   fmcll[6]= 0.684172;  
      Fbb[7]= 1.12083;   Fcc[7]= 1.12083;   Fc[7]= 0.875058;   Fll[7]= 0.998738; 
         cafac[7]= 0.799322;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3097.98 ;   woutpretag[7]= 13797;   wouttag[7]= 2528.24;  
       fbb[7]= 0.0782937;   fcc[7]= 0.128803;   fc[7]= 0.150681;   fll[7]= 0.642222;  
       fmcbb[7]= 0.0698533;   fmccc[7]= 0.114917;   fmcc[7]= 0.172195;   fmcll[7]= 0.643034;  
 }  
 else if ( idsys==46   || sysname=="Wcjet4_down" ) { 
      Fbb[1]= 1.13522;   Fcc[1]= 1.13522;   Fc[1]= 0.886292;   Fll[1]= 1.01156; 
         cafac[1]= 0.972142;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08049e+06;   wouttag[1]= 38391.7;  
       fbb[1]= 0.0135665;   fcc[1]= 0.0396548;   fc[1]= 0.122801;   fll[1]= 0.823978;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.880444;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18364.2;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12441;   Fcc[3]= 1.12441;   Fc[3]= 0.87785;   Fll[3]= 1.00192; 
         cafac[3]= 0.795526;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48498.9;   wouttag[3]= 5791.92;  
       fbb[3]= 0.0553975;   fcc[3]= 0.10794;   fc[3]= 0.139503;   fll[3]= 0.697159;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.11149;   Fcc[4]= 1.11149;   Fc[4]= 0.867766;   Fll[4]= 0.990415; 
         cafac[4]= 0.752804;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2228.56 ;   woutpretag[4]= 10173;   wouttag[4]= 1722.65;  
       fbb[4]= 0.0795919;   fcc[4]= 0.132726;   fc[4]= 0.0958163;   fll[4]= 0.691866;  
       fmcbb[4]= 0.0716084;   fmccc[4]= 0.119413;   fmcc[4]= 0.110417;   fmcll[4]= 0.698562;  
      Fbb[5]= 1.1103;   Fcc[5]= 1.1103;   Fc[5]= 0.866836;   Fll[5]= 0.989354; 
         cafac[5]= 0.722298;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.68;   wouttag[5]= 562.777;  
       fbb[5]= 0.0942562;   fcc[5]= 0.149158;   fc[5]= 0.112274;   fll[5]= 0.644312;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12147;   Fcc[6]= 1.12147;   Fc[6]= 0.87556;   Fll[6]= 0.99931; 
         cafac[6]= 0.782157;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10177.4 ;   woutpretag[6]= 61184.6;   wouttag[6]= 8101.17;  
       fbb[6]= 0.0614949;   fcc[6]= 0.114255;   fc[6]= 0.130571;   fll[6]= 0.693679;  
       fmcbb[6]= 0.0548341;   fmccc[6]= 0.101879;   fmcc[6]= 0.149129;   fmcll[6]= 0.694158;  
      Fbb[7]= 1.11123;   Fcc[7]= 1.11123;   Fc[7]= 0.867564;   Fll[7]= 0.990184; 
         cafac[7]= 0.745381;  
         winpretag[7]= 17260.9 ;   wintag[7]= 2985.45 ;   woutpretag[7]= 12865.9;   wouttag[7]= 2286.52;  
       fbb[7]= 0.0827782;   fcc[7]= 0.136296;   fc[7]= 0.0993923;   fll[7]= 0.681533;  
       fmcbb[7]= 0.0744924;   fmccc[7]= 0.122653;   fmcc[7]= 0.114565;   fmcll[7]= 0.688289;  
 }  
 else if ( idsys==75   || sysname=="Wcjet5_up" ) { 
      Fbb[1]= 1.13522;   Fcc[1]= 1.13522;   Fc[1]= 0.886292;   Fll[1]= 1.01156; 
         cafac[1]= 0.972142;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08049e+06;   wouttag[1]= 38391.7;  
       fbb[1]= 0.0135665;   fcc[1]= 0.0396548;   fc[1]= 0.122801;   fll[1]= 0.823978;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.880444;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18364.2;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12441;   Fcc[3]= 1.12441;   Fc[3]= 0.87785;   Fll[3]= 1.00192; 
         cafac[3]= 0.795526;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48498.9;   wouttag[3]= 5791.92;  
       fbb[3]= 0.0553975;   fcc[3]= 0.10794;   fc[3]= 0.139503;   fll[3]= 0.697159;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.1176;   Fcc[4]= 1.1176;   Fc[4]= 0.872539;   Fll[4]= 0.995863; 
         cafac[4]= 0.787955;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.1;   wouttag[4]= 1840.07;  
       fbb[4]= 0.0767186;   fcc[4]= 0.127934;   fc[4]= 0.128458;   fll[4]= 0.666889;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.11583;   Fcc[5]= 1.11583;   Fc[5]= 0.871151;   Fll[5]= 0.994279; 
         cafac[5]= 0.756519;  
         winpretag[5]= 3747.32 ;   wintag[5]= 768.412 ;   woutpretag[5]= 2834.92;   wouttag[5]= 596.487;  
       fbb[5]= 0.0912018;   fcc[5]= 0.144324;   fc[5]= 0.141041;   fll[5]= 0.623433;  
       fmcbb[5]= 0.0817348;   fmccc[5]= 0.129343;   fmcc[5]= 0.161902;   fmcll[5]= 0.62702;  
      Fbb[6]= 1.12281;   Fcc[6]= 1.12281;   Fc[6]= 0.876606;   Fll[6]= 1.0005; 
         cafac[6]= 0.791554;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10245.1 ;   woutpretag[6]= 61919.8;   wouttag[6]= 8239.95;  
       fbb[6]= 0.0608239;   fcc[6]= 0.113164;   fc[6]= 0.13766;   fll[6]= 0.688351;  
       fmcbb[6]= 0.054171;   fmccc[6]= 0.100787;   fmcc[6]= 0.157038;   fmcll[6]= 0.688004;  
      Fbb[7]= 1.11722;   Fcc[7]= 1.11722;   Fc[7]= 0.872238;   Fll[7]= 0.995519; 
         cafac[7]= 0.780299;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3053.24 ;   woutpretag[7]= 13468.6;   wouttag[7]= 2437.57;  
       fbb[7]= 0.0798668;   fcc[7]= 0.131497;   fc[7]= 0.131193;   fll[7]= 0.657443;  
       fmcbb[7]= 0.0714873;   fmccc[7]= 0.117701;   fmcc[7]= 0.15041;   fmcll[7]= 0.660402;  
 }  
 else if ( idsys==76   || sysname=="Wcjet5_down" ) { 
      Fbb[1]= 1.13522;   Fcc[1]= 1.13522;   Fc[1]= 0.886292;   Fll[1]= 1.01156; 
         cafac[1]= 0.972142;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08049e+06;   wouttag[1]= 38391.7;  
       fbb[1]= 0.0135665;   fcc[1]= 0.0396548;   fc[1]= 0.122801;   fll[1]= 0.823978;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.880444;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236907;   wouttag[2]= 18364.2;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12441;   Fcc[3]= 1.12441;   Fc[3]= 0.87785;   Fll[3]= 1.00192; 
         cafac[3]= 0.795526;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48498.9;   wouttag[3]= 5791.92;  
       fbb[3]= 0.0553975;   fcc[3]= 0.10794;   fc[3]= 0.139503;   fll[3]= 0.697159;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.1176;   Fcc[4]= 1.1176;   Fc[4]= 0.872539;   Fll[4]= 0.995863; 
         cafac[4]= 0.787955;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.1;   wouttag[4]= 1840.07;  
       fbb[4]= 0.0767186;   fcc[4]= 0.127934;   fc[4]= 0.128458;   fll[4]= 0.666889;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.10483;   Fcc[5]= 1.10483;   Fc[5]= 0.862563;   Fll[5]= 0.984477; 
         cafac[5]= 0.691334;  
         winpretag[5]= 3747.32 ;   wintag[5]= 745.37 ;   woutpretag[5]= 2590.65;   wouttag[5]= 532.275;  
       fbb[5]= 0.0972805;   fcc[5]= 0.153944;   fc[5]= 0.0837905;   fll[5]= 0.664985;  
       fmcbb[5]= 0.0880506;   fmccc[5]= 0.139337;   fmcc[5]= 0.0971413;   fmcll[5]= 0.675471;  
      Fbb[6]= 1.12227;   Fcc[6]= 1.12227;   Fc[6]= 0.876185;   Fll[6]= 1.00002; 
         cafac[6]= 0.786871;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10222.1 ;   woutpretag[6]= 61553.4;   wouttag[6]= 8178.66;  
       fbb[6]= 0.0611342;   fcc[6]= 0.113647;   fc[6]= 0.134876;   fll[6]= 0.690342;  
       fmcbb[6]= 0.0544736;   fmccc[6]= 0.101265;   fmcc[6]= 0.153936;   fmcll[6]= 0.690325;  
      Fbb[7]= 1.1148;   Fcc[7]= 1.1148;   Fc[7]= 0.870354;   Fll[7]= 0.993369; 
         cafac[7]= 0.762513;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3030.2 ;   woutpretag[7]= 13161.6;   wouttag[7]= 2368.58;  
       fbb[7]= 0.0812229;   fcc[7]= 0.133632;   fc[7]= 0.118673;   fll[7]= 0.666472;  
       fmcbb[7]= 0.0728584;   fmccc[7]= 0.11987;   fmcc[7]= 0.13635;   fmcll[7]= 0.670921;  
 }  
 else if ( idsys==3   || sysname=="wt_up" ) { 
      Fbb[1]= 1.13473;   Fcc[1]= 1.13473;   Fc[1]= 0.880523;   Fll[1]= 1.01257; 
         cafac[1]= 0.971099;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.07933e+06;   wouttag[1]= 38227.3;  
       fbb[1]= 0.0135607;   fcc[1]= 0.0396379;   fc[1]= 0.122002;   fll[1]= 0.8248;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.13036;   Fcc[2]= 1.13036;   Fc[2]= 0.87713;   Fll[2]= 1.00867; 
         cafac[2]= 0.879354;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236613;   wouttag[2]= 18303;  
       fbb[2]= 0.0331279;   fcc[2]= 0.0799244;   fc[2]= 0.13895;   fll[2]= 0.747998;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12424;   Fcc[3]= 1.12424;   Fc[3]= 0.872383;   Fll[3]= 1.00321; 
         cafac[3]= 0.794477;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48435;   wouttag[3]= 5775.85;  
       fbb[3]= 0.0553895;   fcc[3]= 0.107925;   fc[3]= 0.138634;   fll[3]= 0.698051;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.11741;   Fcc[4]= 1.11741;   Fc[4]= 0.867085;   Fll[4]= 0.997114; 
         cafac[4]= 0.786917;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10634;   wouttag[4]= 1835.95;  
       fbb[4]= 0.0767056;   fcc[4]= 0.127913;   fc[4]= 0.127655;   fll[4]= 0.667727;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.11004;   Fcc[5]= 1.11004;   Fc[5]= 0.86136;   Fll[5]= 0.990531; 
         cafac[5]= 0.721384;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2703.26;   wouttag[5]= 561.639;  
       fbb[5]= 0.094234;   fcc[5]= 0.149123;   fc[5]= 0.111565;   fll[5]= 0.645079;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12237;   Fcc[6]= 1.12237;   Fc[6]= 0.87093;   Fll[6]= 1.00154; 
         cafac[6]= 0.788167;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61654.8;   wouttag[6]= 8187.89;  
       fbb[6]= 0.0609697;   fcc[6]= 0.113388;   fc[6]= 0.135418;   fll[6]= 0.690224;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.1158;   Fcc[7]= 1.1158;   Fc[7]= 0.865835;   Fll[7]= 0.995678; 
         cafac[7]= 0.770286;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13295.8;   wouttag[7]= 2397.38;  
       fbb[7]= 0.0805308;   fcc[7]= 0.132541;   fc[7]= 0.124144;   fll[7]= 0.662784;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==4   || sysname=="wt_down" ) { 
      Fbb[1]= 1.13571;   Fcc[1]= 1.13571;   Fc[1]= 0.892048;   Fll[1]= 1.01055; 
         cafac[1]= 0.973186;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08165e+06;   wouttag[1]= 38556.1;  
       fbb[1]= 0.0135723;   fcc[1]= 0.0396718;   fc[1]= 0.123599;   fll[1]= 0.823157;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.13085;   Fcc[2]= 1.13085;   Fc[2]= 0.888232;   Fll[2]= 1.00623; 
         cafac[2]= 0.881534;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 237200;   wouttag[2]= 18425.3;  
       fbb[2]= 0.0331422;   fcc[2]= 0.0799588;   fc[2]= 0.140708;   fll[2]= 0.746191;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12457;   Fcc[3]= 1.12457;   Fc[3]= 0.883303;   Fll[3]= 1.00065; 
         cafac[3]= 0.796574;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48562.8;   wouttag[3]= 5808;  
       fbb[3]= 0.0554057;   fcc[3]= 0.107956;   fc[3]= 0.14037;   fll[3]= 0.696268;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.11779;   Fcc[4]= 1.11779;   Fc[4]= 0.87798;   Fll[4]= 0.994615; 
         cafac[4]= 0.788993;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10662.1;   wouttag[4]= 1844.2;  
       fbb[4]= 0.0767317;   fcc[4]= 0.127956;   fc[4]= 0.129259;   fll[4]= 0.666053;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.11056;   Fcc[5]= 1.11056;   Fc[5]= 0.872298;   Fll[5]= 0.988179; 
         cafac[5]= 0.723212;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2710.11;   wouttag[5]= 563.916;  
       fbb[5]= 0.0942785;   fcc[5]= 0.149193;   fc[5]= 0.112982;   fll[5]= 0.643547;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12272;   Fcc[6]= 1.12272;   Fc[6]= 0.881846;   Fll[6]= 0.998995; 
         cafac[6]= 0.790243;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61817.2;   wouttag[6]= 8230.52;  
       fbb[6]= 0.0609885;   fcc[6]= 0.113424;   fc[6]= 0.137116;   fll[6]= 0.688472;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.11622;   Fcc[7]= 1.11622;   Fc[7]= 0.87674;   Fll[7]= 0.993211; 
         cafac[7]= 0.772302;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13330.6;   wouttag[7]= 2407.9;  
       fbb[7]= 0.0805605;   fcc[7]= 0.13259;   fc[7]= 0.125707;   fll[7]= 0.661142;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==5   || sysname=="st_up" ) { 
      Fbb[1]= 1.11116;   Fcc[1]= 1.11116;   Fc[1]= 0.887422;   Fll[1]= 1.01275; 
         cafac[1]= 0.972275;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08064e+06;   wouttag[1]= 38299.8;  
       fbb[1]= 0.0132789;   fcc[1]= 0.0388144;   fc[1]= 0.122958;   fll[1]= 0.824949;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.10812;   Fcc[2]= 1.10812;   Fc[2]= 0.884998;   Fll[2]= 1.00998; 
         cafac[2]= 0.880234;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236850;   wouttag[2]= 18263.5;  
       fbb[2]= 0.0324761;   fcc[2]= 0.078352;   fc[2]= 0.140196;   fll[2]= 0.748976;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.10329;   Fcc[3]= 1.10329;   Fc[3]= 0.88114;   Fll[3]= 1.00558; 
         cafac[3]= 0.794401;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48430.4;   wouttag[3]= 5745.89;  
       fbb[3]= 0.0543572;   fcc[3]= 0.105913;   fc[3]= 0.140026;   fll[3]= 0.699703;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.09763;   Fcc[4]= 1.09763;   Fc[4]= 0.876622;   Fll[4]= 1.00043; 
         cafac[4]= 0.786171;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10624;   wouttag[4]= 1824.1;  
       fbb[4]= 0.0753478;   fcc[4]= 0.125649;   fc[4]= 0.129059;   fll[4]= 0.669945;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.09141;   Fcc[5]= 1.09141;   Fc[5]= 0.871652;   Fll[5]= 0.994754; 
         cafac[5]= 0.720396;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2699.56;   wouttag[5]= 557.347;  
       fbb[5]= 0.0926528;   fcc[5]= 0.14662;   fc[5]= 0.112898;   fll[5]= 0.647829;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.10174;   Fcc[6]= 1.10174;   Fc[6]= 0.879898;   Fll[6]= 1.00416; 
         cafac[6]= 0.787942;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61637.1;   wouttag[6]= 8141.55;  
       fbb[6]= 0.0598488;   fcc[6]= 0.111304;   fc[6]= 0.136813;   fll[6]= 0.692035;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.09628;   Fcc[7]= 1.09628;   Fc[7]= 0.875538;   Fll[7]= 0.999189; 
         cafac[7]= 0.769507;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13282.3;   wouttag[7]= 2381.22;  
       fbb[7]= 0.0791215;   fcc[7]= 0.130222;   fc[7]= 0.125535;   fll[7]= 0.665122;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==6   || sysname=="st_down" ) { 
      Fbb[1]= 1.15709;   Fcc[1]= 1.15709;   Fc[1]= 0.885264;   Fll[1]= 1.01048; 
         cafac[1]= 0.972021;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08036e+06;   wouttag[1]= 38475.1;  
       fbb[1]= 0.0138278;   fcc[1]= 0.0404188;   fc[1]= 0.122659;   fll[1]= 0.823095;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.15098;   Fcc[2]= 1.15098;   Fc[2]= 0.880594;   Fll[2]= 1.00514; 
         cafac[2]= 0.880635;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236958;   wouttag[2]= 18455.5;  
       fbb[2]= 0.0337324;   fcc[2]= 0.0813827;   fc[2]= 0.139498;   fll[2]= 0.745387;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.14351;   Fcc[3]= 1.14351;   Fc[3]= 0.874875;   Fll[3]= 0.998616; 
         cafac[3]= 0.796545;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48561.1;   wouttag[3]= 5833.67;  
       fbb[3]= 0.0563387;   fcc[3]= 0.109774;   fc[3]= 0.13903;   fll[3]= 0.694857;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.13564;   Fcc[4]= 1.13564;   Fc[4]= 0.868853;   Fll[4]= 0.991742; 
         cafac[4]= 0.789573;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10669.9;   wouttag[4]= 1854.56;  
       fbb[4]= 0.0779566;   fcc[4]= 0.129999;   fc[4]= 0.127915;   fll[4]= 0.664129;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.12733;   Fcc[5]= 1.12733;   Fc[5]= 0.862494;   Fll[5]= 0.984485; 
         cafac[5]= 0.724024;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2713.15;   wouttag[5]= 567.699;  
       fbb[5]= 0.0957019;   fcc[5]= 0.151445;   fc[5]= 0.111712;   fll[5]= 0.641141;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.14136;   Fcc[6]= 1.14136;   Fc[6]= 0.873229;   Fll[6]= 0.996737; 
         cafac[6]= 0.790351;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61825.6;   wouttag[6]= 8270.57;  
       fbb[6]= 0.0620012;   fcc[6]= 0.115307;   fc[6]= 0.135776;   fll[6]= 0.686916;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.13382;   Fcc[7]= 1.13382;   Fc[7]= 0.867464;   Fll[7]= 0.990158; 
         cafac[7]= 0.772914;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13341.2;   wouttag[7]= 2422.06;  
       fbb[7]= 0.0818313;   fcc[7]= 0.134682;   fc[7]= 0.124377;   fll[7]= 0.65911;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==7   || sysname=="z_up" ) { 
      Fbb[1]= 1.128;   Fcc[1]= 1.128;   Fc[1]= 0.827806;   Fll[1]= 1.02192; 
         cafac[1]= 0.959815;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.06679e+06;   wouttag[1]= 36659.1;  
       fbb[1]= 0.0134802;   fcc[1]= 0.0394026;   fc[1]= 0.114698;   fll[1]= 0.83242;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.12599;   Fcc[2]= 1.12599;   Fc[2]= 0.826334;   Fll[2]= 1.02011; 
         cafac[2]= 0.86693;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 233270;   wouttag[2]= 17687.9;  
       fbb[2]= 0.0329999;   fcc[2]= 0.0796155;   fc[2]= 0.130903;   fll[2]= 0.756482;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12073;   Fcc[3]= 1.12073;   Fc[3]= 0.822473;   Fll[3]= 1.01534; 
         cafac[3]= 0.785982;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 47917.1;   wouttag[3]= 5633.79;  
       fbb[3]= 0.0552164;   fcc[3]= 0.107587;   fc[3]= 0.130703;   fll[3]= 0.706493;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.11378;   Fcc[4]= 1.11378;   Fc[4]= 0.817369;   Fll[4]= 1.00904; 
         cafac[4]= 0.774339;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10464.1;   wouttag[4]= 1790.1;  
       fbb[4]= 0.0764558;   fcc[4]= 0.127496;   fc[4]= 0.120336;   fll[4]= 0.675712;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.10584;   Fcc[5]= 1.10584;   Fc[5]= 0.811546;   Fll[5]= 1.00185; 
         cafac[5]= 0.72385;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2712.5;   wouttag[5]= 559.261;  
       fbb[5]= 0.0938778;   fcc[5]= 0.148559;   fc[5]= 0.105113;   fll[5]= 0.652451;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.1188;   Fcc[6]= 1.1188;   Fc[6]= 0.821058;   Fll[6]= 1.01359; 
         cafac[6]= 0.779595;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 60984.2;   wouttag[6]= 7998.17;  
       fbb[6]= 0.0607758;   fcc[6]= 0.113028;   fc[6]= 0.127664;   fll[6]= 0.698533;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.11204;   Fcc[7]= 1.11204;   Fc[7]= 0.816098;   Fll[7]= 1.00747; 
         cafac[7]= 0.761478;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13143.8;   wouttag[7]= 2349.27;  
       fbb[7]= 0.0802593;   fcc[7]= 0.132095;   fc[7]= 0.117012;   fll[7]= 0.670634;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==8   || sysname=="z_down" ) { 
      Fbb[1]= 1.14226;   Fcc[1]= 1.14226;   Fc[1]= 0.943266;   Fll[1]= 1.00146; 
         cafac[1]= 0.984461;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.09418e+06;   wouttag[1]= 40122.9;  
       fbb[1]= 0.0136506;   fcc[1]= 0.0399007;   fc[1]= 0.130695;   fll[1]= 0.815754;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.13508;   Fcc[2]= 1.13508;   Fc[2]= 0.937337;   Fll[2]= 0.995168; 
         cafac[2]= 0.893959;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 240543;   wouttag[2]= 19040.4;  
       fbb[2]= 0.0332661;   fcc[2]= 0.0802579;   fc[2]= 0.148487;   fll[2]= 0.737989;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12797;   Fcc[3]= 1.12797;   Fc[3]= 0.931467;   Fll[3]= 0.988936; 
         cafac[3]= 0.804983;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 49075.5;   wouttag[3]= 5948.68;  
       fbb[3]= 0.055573;   fcc[3]= 0.108282;   fc[3]= 0.148024;   fll[3]= 0.688121;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.12131;   Fcc[4]= 1.12131;   Fc[4]= 0.925972;   Fll[4]= 0.983102; 
         cafac[4]= 0.801596;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10832.4;   wouttag[4]= 1890.15;  
       fbb[4]= 0.0769733;   fcc[4]= 0.128359;   fc[4]= 0.136324;   fll[4]= 0.658343;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.11462;   Fcc[5]= 1.11462;   Fc[5]= 0.920447;   Fll[5]= 0.977236; 
         cafac[5]= 0.720457;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2699.78;   wouttag[5]= 565.904;  
       fbb[5]= 0.0946233;   fcc[5]= 0.149739;   fc[5]= 0.119218;   fll[5]= 0.63642;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12617;   Fcc[6]= 1.12617;   Fc[6]= 0.92998;   Fll[6]= 0.987357; 
         cafac[6]= 0.798737;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 62481.6;   wouttag[6]= 8418.55;  
       fbb[6]= 0.061176;   fcc[6]= 0.113772;   fc[6]= 0.1446;   fll[6]= 0.680452;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.11985;   Fcc[7]= 1.11985;   Fc[7]= 0.924767;   Fll[7]= 0.981822; 
         cafac[7]= 0.781052;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13481.6;   wouttag[7]= 2455.69;  
       fbb[7]= 0.0808231;   fcc[7]= 0.133022;   fc[7]= 0.132593;   fll[7]= 0.653561;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==9   || sysname=="di_up" ) { 
      Fbb[1]= 1.13385;   Fcc[1]= 1.13385;   Fc[1]= 0.885137;   Fll[1]= 1.01183; 
         cafac[1]= 0.971901;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08022e+06;   wouttag[1]= 38351.1;  
       fbb[1]= 0.0135501;   fcc[1]= 0.0396071;   fc[1]= 0.122641;   fll[1]= 0.824202;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.12937;   Fcc[2]= 1.12937;   Fc[2]= 0.88164;   Fll[2]= 1.00784; 
         cafac[2]= 0.88005;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236801;   wouttag[2]= 18342.8;  
       fbb[2]= 0.0330989;   fcc[2]= 0.0798544;   fc[2]= 0.139664;   fll[2]= 0.747383;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12326;   Fcc[3]= 1.12326;   Fc[3]= 0.876873;   Fll[3]= 1.00239; 
         cafac[3]= 0.795069;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48471.1;   wouttag[3]= 5784.83;  
       fbb[3]= 0.0553413;   fcc[3]= 0.107831;   fc[3]= 0.139348;   fll[3]= 0.69748;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.11652;   Fcc[4]= 1.11652;   Fc[4]= 0.871606;   Fll[4]= 0.996366; 
         cafac[4]= 0.787536;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10642.4;   wouttag[4]= 1838.12;  
       fbb[4]= 0.076644;   fcc[4]= 0.12781;   fc[4]= 0.12832;   fll[4]= 0.667226;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.10925;   Fcc[5]= 1.10925;   Fc[5]= 0.865936;   Fll[5]= 0.989884; 
         cafac[5]= 0.722106;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2705.97;   wouttag[5]= 562.331;  
       fbb[5]= 0.0941676;   fcc[5]= 0.149017;   fc[5]= 0.112158;   fll[5]= 0.644658;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12141;   Fcc[6]= 1.12141;   Fc[6]= 0.875429;   Fll[6]= 1.00074; 
         cafac[6]= 0.788774;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61702.3;   wouttag[6]= 8199.65;  
       fbb[6]= 0.0609178;   fcc[6]= 0.113292;   fc[6]= 0.136118;   fll[6]= 0.689673;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.11493;   Fcc[7]= 1.11493;   Fc[7]= 0.870368;   Fll[7]= 0.994951; 
         cafac[7]= 0.770933;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13307;   wouttag[7]= 2400.24;  
       fbb[7]= 0.0804678;   fcc[7]= 0.132438;   fc[7]= 0.124793;   fll[7]= 0.662301;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==10   || sysname=="di_down" ) { 
      Fbb[1]= 1.13659;   Fcc[1]= 1.13659;   Fc[1]= 0.887445;   Fll[1]= 1.01128; 
         cafac[1]= 0.972383;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08076e+06;   wouttag[1]= 38432.2;  
       fbb[1]= 0.0135828;   fcc[1]= 0.0397027;   fc[1]= 0.122961;   fll[1]= 0.823754;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.13184;   Fcc[2]= 1.13184;   Fc[2]= 0.883735;   Fll[2]= 1.00706; 
         cafac[2]= 0.880839;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 237013;   wouttag[2]= 18385.6;  
       fbb[2]= 0.0331711;   fcc[2]= 0.0800287;   fc[2]= 0.139996;   fll[2]= 0.746804;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12555;   Fcc[3]= 1.12555;   Fc[3]= 0.878827;   Fll[3]= 1.00146; 
         cafac[3]= 0.795982;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48526.8;   wouttag[3]= 5799.01;  
       fbb[3]= 0.0554539;   fcc[3]= 0.10805;   fc[3]= 0.139658;   fll[3]= 0.696838;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.11869;   Fcc[4]= 1.11869;   Fc[4]= 0.873472;   Fll[4]= 0.995361; 
         cafac[4]= 0.788374;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10653.7;   wouttag[4]= 1842.03;  
       fbb[4]= 0.0767933;   fcc[4]= 0.128059;   fc[4]= 0.128595;   fll[4]= 0.666553;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.11134;   Fcc[5]= 1.11134;   Fc[5]= 0.867735;   Fll[5]= 0.988823; 
         cafac[5]= 0.72249;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2707.4;   wouttag[5]= 563.224;  
       fbb[5]= 0.094345;   fcc[5]= 0.149298;   fc[5]= 0.112391;   fll[5]= 0.643966;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12367;   Fcc[6]= 1.12367;   Fc[6]= 0.87736;   Fll[6]= 0.999792; 
         cafac[6]= 0.789636;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61769.7;   wouttag[6]= 8218.77;  
       fbb[6]= 0.0610404;   fcc[6]= 0.11352;   fc[6]= 0.136418;   fll[6]= 0.689022;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.11709;   Fcc[7]= 1.11709;   Fc[7]= 0.87222;   Fll[7]= 0.993934; 
         cafac[7]= 0.771654;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13319.4;   wouttag[7]= 2405.04;  
       fbb[7]= 0.0806234;   fcc[7]= 0.132694;   fc[7]= 0.125059;   fll[7]= 0.661624;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==11   || sysname=="qcd_up" ) { 
      Fbb[1]= 1.122;   Fcc[1]= 1.122;   Fc[1]= 0.72366;   Fll[1]= 1.03998; 
         cafac[1]= 0.943515;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.04867e+06;   wouttag[1]= 33889.5;  
       fbb[1]= 0.0134085;   fcc[1]= 0.0391932;   fc[1]= 0.100268;   fll[1]= 0.847131;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.12417;   Fcc[2]= 1.12417;   Fc[2]= 0.725055;   Fll[2]= 1.04199; 
         cafac[2]= 0.850544;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 228861;   wouttag[2]= 16689.3;  
       fbb[2]= 0.0329463;   fcc[2]= 0.0794864;   fc[2]= 0.114859;   fll[2]= 0.772708;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12018;   Fcc[3]= 1.12018;   Fc[3]= 0.722483;   Fll[3]= 1.03829; 
         cafac[3]= 0.767375;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 46782.8;   wouttag[3]= 5356.04;  
       fbb[3]= 0.0551892;   fcc[3]= 0.107534;   fc[3]= 0.114813;   fll[3]= 0.722463;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.11262;   Fcc[4]= 1.11262;   Fc[4]= 0.717609;   Fll[4]= 1.03129; 
         cafac[4]= 0.762513;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10304.3;   wouttag[4]= 1734.05;  
       fbb[4]= 0.0763765;   fcc[4]= 0.127364;   fc[4]= 0.105649;   fll[4]= 0.690611;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.10325;   Fcc[5]= 1.10325;   Fc[5]= 0.711568;   Fll[5]= 1.02261; 
         cafac[5]= 0.699831;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2622.49;   wouttag[5]= 533.644;  
       fbb[5]= 0.0936582;   fcc[5]= 0.148211;   fc[5]= 0.0921636;   fll[5]= 0.665967;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.11804;   Fcc[6]= 1.11804;   Fc[6]= 0.721107;   Fll[6]= 1.03631; 
         cafac[6]= 0.761919;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 59601.5;   wouttag[6]= 7637.75;  
       fbb[6]= 0.0607347;   fcc[6]= 0.112951;   fc[6]= 0.112123;   fll[6]= 0.714191;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.11057;   Fcc[7]= 1.11057;   Fc[7]= 0.716289;   Fll[7]= 1.02939; 
         cafac[7]= 0.746572;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 12886.5;   wouttag[7]= 2267.74;  
       fbb[7]= 0.0801533;   fcc[7]= 0.13192;   fc[7]= 0.102702;   fll[7]= 0.685225;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==12   || sysname=="qcd_down" ) { 
      Fbb[1]= 1.14771;   Fcc[1]= 1.14771;   Fc[1]= 1.04002;   Fll[1]= 0.984691; 
         cafac[1]= 1.00085;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.1124e+06;   wouttag[1]= 42906;  
       fbb[1]= 0.0137158;   fcc[1]= 0.0400913;   fc[1]= 0.144101;   fll[1]= 0.802092;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.13662;   Fcc[2]= 1.13662;   Fc[2]= 1.02997;   Fll[2]= 0.975173; 
         cafac[2]= 0.910345;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 244952;   wouttag[2]= 20039.1;  
       fbb[2]= 0.0333113;   fcc[2]= 0.0803669;   fc[2]= 0.163161;   fll[2]= 0.723161;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12834;   Fcc[3]= 1.12834;   Fc[3]= 1.02247;   Fll[3]= 0.968074; 
         cafac[3]= 0.82365;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 50213.5;   wouttag[3]= 6227.39;  
       fbb[3]= 0.0555915;   fcc[3]= 0.108318;   fc[3]= 0.162485;   fll[3]= 0.673605;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.12225;   Fcc[4]= 1.12225;   Fc[4]= 1.01695;   Fll[4]= 0.962846; 
         cafac[4]= 0.813246;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10989.8;   wouttag[4]= 1945.47;  
       fbb[4]= 0.0770375;   fcc[4]= 0.128466;   fc[4]= 0.149718;   fll[4]= 0.644778;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.11689;   Fcc[5]= 1.11689;   Fc[5]= 1.01209;   Fll[5]= 0.958247; 
         cafac[5]= 0.744662;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2790.49;   wouttag[5]= 591.777;  
       fbb[5]= 0.0948158;   fcc[5]= 0.150043;   fc[5]= 0.131087;   fll[5]= 0.624054;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12673;   Fcc[6]= 1.12673;   Fc[6]= 1.02101;   Fll[6]= 0.966692; 
         cafac[6]= 0.816434;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 63865.9;   wouttag[6]= 8779.45;  
       fbb[6]= 0.0612067;   fcc[6]= 0.113829;   fc[6]= 0.158753;   fll[6]= 0.666211;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.12108;   Fcc[7]= 1.12108;   Fc[7]= 1.01589;   Fll[7]= 0.961843; 
         cafac[7]= 0.795877;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13737.5;   wouttag[7]= 2536.79;  
       fbb[7]= 0.0809116;   fcc[7]= 0.133168;   fc[7]= 0.145658;   fll[7]= 0.640262;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==13   || sysname=="jes_up" ) { 
      Fbb[1]= 1.13151;   Fcc[1]= 1.13151;   Fc[1]= 0.941565;   Fll[1]= 1.00257; 
         cafac[1]= 0.913043;  
         winpretag[1]= 1.19744e+06 ;   wintag[1]= 43166.2 ;   woutpretag[1]= 1.09332e+06;   wouttag[1]= 38851;  
       fbb[1]= 0.0126207;   fcc[1]= 0.037121;   fc[1]= 0.127092;   fll[1]= 0.823166;  
       fmcbb[1]= 0.0111539;   fmccc[1]= 0.0328066;   fmcc[1]= 0.13498;   fmcll[1]= 0.82106;  
      Fbb[2]= 1.12563;   Fcc[2]= 1.12563;   Fc[2]= 0.936669;   Fll[2]= 0.997352; 
         cafac[2]= 0.814676;  
         winpretag[2]= 298038 ;   wintag[2]= 22921.3 ;   woutpretag[2]= 242804;   wouttag[2]= 18845.6;  
       fbb[2]= 0.0309907;   fcc[2]= 0.0765733;   fc[2]= 0.148337;   fll[2]= 0.744099;  
       fmcbb[2]= 0.027532;   fmccc[2]= 0.0680272;   fmcc[2]= 0.158366;   fmcll[2]= 0.746075;  
      Fbb[3]= 1.11951;   Fcc[3]= 1.11951;   Fc[3]= 0.931583;   Fll[3]= 0.991936; 
         cafac[3]= 0.704281;  
         winpretag[3]= 70845.6 ;   wintag[3]= 8229.4 ;   woutpretag[3]= 49895.2;   wouttag[3]= 5928.11;  
       fbb[3]= 0.0522235;   fcc[3]= 0.103755;   fc[3]= 0.149897;   fll[3]= 0.694125;  
       fmcbb[3]= 0.0466483;   fmccc[3]= 0.0926782;   fmcc[3]= 0.160906;   fmcll[3]= 0.699768;  
      Fbb[4]= 1.11328;   Fcc[4]= 1.11328;   Fc[4]= 0.926391;   Fll[4]= 0.986409; 
         cafac[4]= 0.670308;  
         winpretag[4]= 16524.1 ;   wintag[4]= 2679.86 ;   woutpretag[4]= 11076.2;   wouttag[4]= 1846.91;  
       fbb[4]= 0.0741152;   fcc[4]= 0.123469;   fc[4]= 0.137763;   fll[4]= 0.664652;  
       fmcbb[4]= 0.066574;   fmccc[4]= 0.110906;   fmcc[4]= 0.14871;   fmcll[4]= 0.67381;  
      Fbb[5]= 1.1069;   Fcc[5]= 1.1069;   Fc[5]= 0.921089;   Fll[5]= 0.980763; 
         cafac[5]= 0.562124;  
         winpretag[5]= 4798.16 ;   wintag[5]= 958.955 ;   woutpretag[5]= 2697.16;   wouttag[5]= 557.274;  
       fbb[5]= 0.0912081;   fcc[5]= 0.147675;   fc[5]= 0.12326;   fll[5]= 0.637857;  
       fmcbb[5]= 0.0823992;   fmccc[5]= 0.133412;   fmcc[5]= 0.13382;   fmcll[5]= 0.650369;  
      Fbb[6]= 1.11773;   Fcc[6]= 1.11773;   Fc[6]= 0.930097;   Fll[6]= 0.990354; 
         cafac[6]= 0.687642;  
         winpretag[6]= 92167.9 ;   wintag[6]= 11868.2 ;   woutpretag[6]= 63378.5;   wouttag[6]= 8369.44;  
       fbb[6]= 0.0582133;   fcc[6]= 0.109612;   fc[6]= 0.146313;   fll[6]= 0.685862;  
       fmcbb[6]= 0.0520818;   fmccc[6]= 0.0980668;   fmcc[6]= 0.157309;   fmcll[6]= 0.692542;  
      Fbb[7]= 1.11184;   Fcc[7]= 1.11184;   Fc[7]= 0.925193;   Fll[7]= 0.985132; 
         cafac[7]= 0.641897;  
         winpretag[7]= 21322.2 ;   wintag[7]= 3638.82 ;   woutpretag[7]= 13686.7;   wouttag[7]= 2405.53;  
       fbb[7]= 0.0779788;   fcc[7]= 0.128941;   fc[7]= 0.134485;   fll[7]= 0.658596;  
       fmcbb[7]= 0.0701352;   fmccc[7]= 0.115971;   fmcc[7]= 0.145359;   fmcll[7]= 0.668535;  
 }  
 else if ( idsys==14   || sysname=="jes_down" ) { 
      Fbb[1]= 1.17417;   Fcc[1]= 1.17417;   Fc[1]= 0.820453;   Fll[1]= 1.02073; 
         cafac[1]= 1.02916;  
         winpretag[1]= 1.032e+06 ;   wintag[1]= 39303 ;   woutpretag[1]= 1.06209e+06;   wouttag[1]= 37490.1;  
       fbb[1]= 0.0149757;   fcc[1]= 0.0432764;   fc[1]= 0.116104;   fll[1]= 0.825644;  
       fmcbb[1]= 0.0127542;   fmccc[1]= 0.0368569;   fmcc[1]= 0.141512;   fmcll[1]= 0.808877;  
      Fbb[2]= 1.16829;   Fcc[2]= 1.16829;   Fc[2]= 0.816343;   Fll[2]= 1.01562; 
         cafac[2]= 0.942554;  
         winpretag[2]= 242652 ;   wintag[2]= 19398.9 ;   woutpretag[2]= 228713;   wouttag[2]= 17923.7;  
       fbb[2]= 0.0359217;   fcc[2]= 0.0857721;   fc[2]= 0.129119;   fll[2]= 0.749188;  
       fmcbb[2]= 0.0307472;   fmccc[2]= 0.0734167;   fmcc[2]= 0.158167;   fmcll[2]= 0.737669;  
      Fbb[3]= 1.15989;   Fcc[3]= 1.15989;   Fc[3]= 0.810474;   Fll[3]= 1.00831; 
         cafac[3]= 0.902474;  
         winpretag[3]= 52071.7 ;   wintag[3]= 6266.02 ;   woutpretag[3]= 46993.4;   wouttag[3]= 5713.03;  
       fbb[3]= 0.0602083;   fcc[3]= 0.114679;   fc[3]= 0.127687;   fll[3]= 0.697426;  
       fmcbb[3]= 0.0519085;   fmccc[3]= 0.0988703;   fmcc[3]= 0.157546;   fmcll[3]= 0.691675;  
      Fbb[4]= 1.14982;   Fcc[4]= 1.14982;   Fc[4]= 0.803437;   Fll[4]= 0.999559; 
         cafac[4]= 0.91493;  
         winpretag[4]= 11286.4 ;   wintag[4]= 1988.11 ;   woutpretag[4]= 10326.3;   wouttag[4]= 1866.92;  
       fbb[4]= 0.0831197;   fcc[4]= 0.135031;   fc[4]= 0.114984;   fll[4]= 0.666865;  
       fmcbb[4]= 0.0722892;   fmccc[4]= 0.117437;   fmcc[4]= 0.143115;   fmcll[4]= 0.667159;  
      Fbb[5]= 1.13996;   Fcc[5]= 1.13996;   Fc[5]= 0.796543;   Fll[5]= 0.990983; 
         cafac[5]= 0.966187;  
         winpretag[5]= 2922.38 ;   wintag[5]= 594.552 ;   woutpretag[5]= 2823.57;   wouttag[5]= 593.782;  
       fbb[5]= 0.0979943;   fcc[5]= 0.157921;   fc[5]= 0.100067;   fll[5]= 0.644018;  
       fmcbb[5]= 0.0859633;   fmccc[5]= 0.138533;   fmcc[5]= 0.125626;   fmcll[5]= 0.649878;  
      Fbb[6]= 1.15727;   Fcc[6]= 1.15727;   Fc[6]= 0.808645;   Fll[6]= 1.00604; 
         cafac[6]= 0.908022;  
         winpretag[6]= 66280.5 ;   wintag[6]= 8848.67 ;   woutpretag[6]= 60184.2;   wouttag[6]= 8166.52;  
       fbb[6]= 0.0658264;   fcc[6]= 0.120103;   fc[6]= 0.124274;   fll[6]= 0.689797;  
       fmcbb[6]= 0.0568805;   fmccc[6]= 0.103781;   fmcc[6]= 0.153682;   fmcll[6]= 0.685657;  
      Fbb[7]= 1.14778;   Fcc[7]= 1.14778;   Fc[7]= 0.80201;   Fll[7]= 0.997783; 
         cafac[7]= 0.926361;  
         winpretag[7]= 14208.8 ;   wintag[7]= 2582.66 ;   woutpretag[7]= 13162.5;   wouttag[7]= 2460.1;  
       fbb[7]= 0.0862001;   fcc[7]= 0.139771;   fc[7]= 0.111895;   fll[7]= 0.662133;  
       fmcbb[7]= 0.0751016;   fmccc[7]= 0.121776;   fmcc[7]= 0.139518;   fmcll[7]= 0.663604;  
 }  
 else if ( idsys==79   || sysname=="leff_up" ) { 
      Fbb[1]= 1.11951;   Fcc[1]= 1.11951;   Fc[1]= 0.868518;   Fll[1]= 1.01548; 
         cafac[1]= 0.945029;  
         winpretag[1]= 1.13918e+06 ;   wintag[1]= 42270 ;   woutpretag[1]= 1.07656e+06;   wouttag[1]= 37799.8;  
       fbb[1]= 0.0133802;   fcc[1]= 0.0391106;   fc[1]= 0.120312;   fll[1]= 0.827198;  
       fmcbb[1]= 0.0119518;   fmccc[1]= 0.0349353;   fmcc[1]= 0.138525;   fmcll[1]= 0.814588;  
      Fbb[2]= 1.1166;   Fcc[2]= 1.1166;   Fc[2]= 0.866256;   Fll[2]= 1.01284; 
         cafac[2]= 0.855203;  
         winpretag[2]= 275816 ;   wintag[2]= 21538.9 ;   woutpretag[2]= 235879;   wouttag[2]= 18105.3;  
       fbb[2]= 0.0327272;   fcc[2]= 0.0789598;   fc[2]= 0.137193;   fll[2]= 0.75112;  
       fmcbb[2]= 0.0293098;   fmccc[2]= 0.0707147;   fmcc[2]= 0.158374;   fmcll[2]= 0.741601;  
      Fbb[3]= 1.11146;   Fcc[3]= 1.11146;   Fc[3]= 0.862269;   Fll[3]= 1.00817; 
         cafac[3]= 0.772198;  
         winpretag[3]= 62497.8 ;   wintag[3]= 7371.86 ;   woutpretag[3]= 48260.7;   wouttag[3]= 5713.71;  
       fbb[3]= 0.0547599;   fcc[3]= 0.106702;   fc[3]= 0.136982;   fll[3]= 0.701557;  
       fmcbb[3]= 0.0492685;   fmccc[3]= 0.0960014;   fmcc[3]= 0.158862;   fmcll[3]= 0.695868;  
      Fbb[4]= 1.10525;   Fcc[4]= 1.10525;   Fc[4]= 0.857455;   Fll[4]= 1.00255; 
         cafac[4]= 0.764385;  
         winpretag[4]= 13855.2 ;   wintag[4]= 2342.6 ;   woutpretag[4]= 10590.7;   wouttag[4]= 1817.82;  
       fbb[4]= 0.0758744;   fcc[4]= 0.126524;   fc[4]= 0.126195;   fll[4]= 0.671407;  
       fmcbb[4]= 0.0686489;   fmccc[4]= 0.114475;   fmcc[4]= 0.147174;   fmcll[4]= 0.669702;  
      Fbb[5]= 1.09836;   Fcc[5]= 1.09836;   Fc[5]= 0.852106;   Fll[5]= 0.996291; 
         cafac[5]= 0.701581;  
         winpretag[5]= 3842.42 ;   wintag[5]= 776.279 ;   woutpretag[5]= 2695.77;   wouttag[5]= 556.965;  
       fbb[5]= 0.0932686;   fcc[5]= 0.147542;   fc[5]= 0.110329;   fll[5]= 0.64886;  
       fmcbb[5]= 0.0849164;   fmccc[5]= 0.13433;   fmcc[5]= 0.129478;   fmcll[5]= 0.651276;  
      Fbb[6]= 1.10975;   Fcc[6]= 1.10975;   Fc[6]= 0.860942;   Fll[6]= 1.00662; 
         cafac[6]= 0.766024;  
         winpretag[6]= 80195.5 ;   wintag[6]= 10490.7 ;   woutpretag[6]= 61431.6;   wouttag[6]= 8102.85;  
       fbb[6]= 0.0602869;   fcc[6]= 0.112117;   fc[6]= 0.13382;   fll[6]= 0.693775;  
       fmcbb[6]= 0.0543248;   fmccc[6]= 0.10103;   fmcc[6]= 0.155435;   fmcll[6]= 0.689211;  
      Fbb[7]= 1.10375;   Fcc[7]= 1.10375;   Fc[7]= 0.856288;   Fll[7]= 1.00118; 
         cafac[7]= 0.748461;  
         winpretag[7]= 17697.6 ;   wintag[7]= 3118.88 ;   woutpretag[7]= 13246;   wouttag[7]= 2374.62;  
       fbb[7]= 0.0796695;   fcc[7]= 0.13111;   fc[7]= 0.122733;   fll[7]= 0.666487;  
       fmcbb[7]= 0.0721808;   fmccc[7]= 0.118786;   fmcc[7]= 0.143332;   fmcll[7]= 0.665701;  
 }  
 else if ( idsys==80   || sysname=="leff_down" ) { 
      Fbb[1]= 1.15083;   Fcc[1]= 1.15083;   Fc[1]= 0.903957;   Fll[1]= 1.00766; 
         cafac[1]= 1.00066;  
         winpretag[1]= 1.08373e+06 ;   wintag[1]= 40218.9 ;   woutpretag[1]= 1.08444e+06;   wouttag[1]= 38984.9;  
       fbb[1]= 0.0137515;   fcc[1]= 0.0401955;   fc[1]= 0.125278;   fll[1]= 0.820775;  
       fmcbb[1]= 0.0119492;   fmccc[1]= 0.0349274;   fmcc[1]= 0.138589;   fmcll[1]= 0.814535;  
      Fbb[2]= 1.14448;   Fcc[2]= 1.14448;   Fc[2]= 0.898971;   Fll[2]= 1.0021; 
         cafac[2]= 0.906983;  
         winpretag[2]= 262336 ;   wintag[2]= 20487.5 ;   woutpretag[2]= 237934;   wouttag[2]= 18623.1;  
       fbb[2]= 0.033539;   fcc[2]= 0.0809139;   fc[2]= 0.142447;   fll[2]= 0.7431;  
       fmcbb[2]= 0.0293049;   fmccc[2]= 0.070699;   fmcc[2]= 0.158456;   fmcll[2]= 0.74154;  
      Fbb[3]= 1.13722;   Fcc[3]= 1.13722;   Fc[3]= 0.893262;   Fll[3]= 0.99574; 
         cafac[3]= 0.820051;  
         winpretag[3]= 59431.5 ;   wintag[3]= 7011.94 ;   woutpretag[3]= 48736.8;   wouttag[3]= 5870.11;  
       fbb[3]= 0.0560284;   fcc[3]= 0.109166;   fc[3]= 0.142002;   fll[3]= 0.692804;  
       fmcbb[3]= 0.049268;   fmccc[3]= 0.0959937;   fmcc[3]= 0.15897;   fmcll[3]= 0.695768;  
      Fbb[4]= 1.12981;   Fcc[4]= 1.12981;   Fc[4]= 0.887447;   Fll[4]= 0.989257; 
         cafac[4]= 0.81274;  
         winpretag[4]= 13171.9 ;   wintag[4]= 2227.05 ;   woutpretag[4]= 10705.3;   wouttag[4]= 1862.3;  
       fbb[4]= 0.0775528;   fcc[4]= 0.129328;   fc[4]= 0.130699;   fll[4]= 0.662421;  
       fmcbb[4]= 0.0686422;   fmccc[4]= 0.114469;   fmcc[4]= 0.147275;   fmcll[4]= 0.669614;  
      Fbb[5]= 1.12209;   Fcc[5]= 1.12209;   Fc[5]= 0.881385;   Fll[5]= 0.982499; 
         cafac[5]= 0.744074;  
         winpretag[5]= 3652.22 ;   wintag[5]= 737.503 ;   woutpretag[5]= 2717.52;   wouttag[5]= 568.552;  
       fbb[5]= 0.0952297;   fcc[5]= 0.150754;   fc[5]= 0.114199;   fll[5]= 0.639817;  
       fmcbb[5]= 0.0848678;   fmccc[5]= 0.134351;   fmcc[5]= 0.129568;   fmcll[5]= 0.651213;  
      Fbb[6]= 1.1352;   Fcc[6]= 1.1352;   Fc[6]= 0.891677;   Fll[6]= 0.993973; 
         cafac[6]= 0.813578;  
         winpretag[6]= 76255.5 ;   wintag[6]= 9976.49 ;   woutpretag[6]= 62039.8;   wouttag[6]= 8315.47;  
       fbb[6]= 0.0616635;   fcc[6]= 0.11468;   fc[6]= 0.138693;   fll[6]= 0.684963;  
       fmcbb[6]= 0.0543196;   fmccc[6]= 0.101022;   fmcc[6]= 0.155542;   fmcll[6]= 0.689117;  
      Fbb[7]= 1.12813;   Fcc[7]= 1.12813;   Fc[7]= 0.886124;   Fll[7]= 0.987782; 
         cafac[7]= 0.795301;  
         winpretag[7]= 16824.1 ;   wintag[7]= 2964.56 ;   woutpretag[7]= 13380.2;   wouttag[7]= 2430.59;  
       fbb[7]= 0.0814108;   fcc[7]= 0.134004;   fc[7]= 0.127098;   fll[7]= 0.657487;  
       fmcbb[7]= 0.0721645;   fmccc[7]= 0.118785;   fmcc[7]= 0.143431;   fmcll[7]= 0.66562;  
 }  
 else if ( idsys==81   || sysname=="ca_up" ) { 
      Fbb[1]= 1.13522;   Fcc[1]= 1.13522;   Fc[1]= 0.886292;   Fll[1]= 1.01156; 
         cafac[1]= 0.977977;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08698e+06;   wouttag[1]= 38622.1;  
       fbb[1]= 0.0135665;   fcc[1]= 0.0396548;   fc[1]= 0.122801;   fll[1]= 0.823978;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.893841;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 240512;   wouttag[2]= 18643.6;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12441;   Fcc[3]= 1.12441;   Fc[3]= 0.87785;   Fll[3]= 1.00192; 
         cafac[3]= 0.818815;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 49918.8;   wouttag[3]= 5961.48;  
       fbb[3]= 0.0553975;   fcc[3]= 0.10794;   fc[3]= 0.139503;   fll[3]= 0.697159;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.1176;   Fcc[4]= 1.1176;   Fc[4]= 0.872539;   Fll[4]= 0.995863; 
         cafac[4]= 0.840949;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 11364.2;   wouttag[4]= 1963.83;  
       fbb[4]= 0.0767186;   fcc[4]= 0.127934;   fc[4]= 0.128458;   fll[4]= 0.666889;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.1103;   Fcc[5]= 1.1103;   Fc[5]= 0.866836;   Fll[5]= 0.989354; 
         cafac[5]= 0.830743;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 3113.06;   wouttag[5]= 647.272;  
       fbb[5]= 0.0942562;   fcc[5]= 0.149158;   fc[5]= 0.112274;   fll[5]= 0.644312;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12254;   Fcc[6]= 1.12254;   Fc[6]= 0.876395;   Fll[6]= 1.00026; 
         cafac[6]= 0.810433;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 63396.6;   wouttag[6]= 8430.02;  
       fbb[6]= 0.0609791;   fcc[6]= 0.113406;   fc[6]= 0.136268;   fll[6]= 0.689347;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.11601;   Fcc[7]= 1.11601;   Fc[7]= 0.871295;   Fll[7]= 0.994443; 
         cafac[7]= 0.819414;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 14143.8;   wouttag[7]= 2552.54;  
       fbb[7]= 0.0805456;   fcc[7]= 0.132566;   fc[7]= 0.124926;   fll[7]= 0.661962;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==82   || sysname=="ca_down" ) { 
      Fbb[1]= 1.13522;   Fcc[1]= 1.13522;   Fc[1]= 0.886292;   Fll[1]= 1.01156; 
         cafac[1]= 0.966308;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.07401e+06;   wouttag[1]= 38161.3;  
       fbb[1]= 0.0135665;   fcc[1]= 0.0396548;   fc[1]= 0.122801;   fll[1]= 0.823978;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.867047;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 233302;   wouttag[2]= 18084.8;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12441;   Fcc[3]= 1.12441;   Fc[3]= 0.87785;   Fll[3]= 1.00192; 
         cafac[3]= 0.772236;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 47079.1;   wouttag[3]= 5622.36;  
       fbb[3]= 0.0553975;   fcc[3]= 0.10794;   fc[3]= 0.139503;   fll[3]= 0.697159;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.1176;   Fcc[4]= 1.1176;   Fc[4]= 0.872539;   Fll[4]= 0.995863; 
         cafac[4]= 0.734961;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 9931.92;   wouttag[4]= 1716.32;  
       fbb[4]= 0.0767186;   fcc[4]= 0.127934;   fc[4]= 0.128458;   fll[4]= 0.666889;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.1103;   Fcc[5]= 1.1103;   Fc[5]= 0.866836;   Fll[5]= 0.989354; 
         cafac[5]= 0.613853;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2300.3;   wouttag[5]= 478.282;  
       fbb[5]= 0.0942562;   fcc[5]= 0.149158;   fc[5]= 0.112274;   fll[5]= 0.644312;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12254;   Fcc[6]= 1.12254;   Fc[6]= 0.876395;   Fll[6]= 1.00026; 
         cafac[6]= 0.767977;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 60075.4;   wouttag[6]= 7988.39;  
       fbb[6]= 0.0609791;   fcc[6]= 0.113406;   fc[6]= 0.136268;   fll[6]= 0.689347;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.11601;   Fcc[7]= 1.11601;   Fc[7]= 0.871295;   Fll[7]= 0.994443; 
         cafac[7]= 0.723173;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 12482.6;   wouttag[7]= 2252.74;  
       fbb[7]= 0.0805456;   fcc[7]= 0.132566;   fc[7]= 0.124926;   fll[7]= 0.661962;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==83   || sysname=="chmis_up" ) { 
      Fbb[1]= 1.13522;   Fcc[1]= 1.13522;   Fc[1]= 0.886292;   Fll[1]= 1.01156; 
         cafac[1]= 0.977975;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08698e+06;   wouttag[1]= 38622;  
       fbb[1]= 0.0135665;   fcc[1]= 0.0396548;   fc[1]= 0.122801;   fll[1]= 0.823978;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.884847;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 238091;   wouttag[2]= 18456;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12441;   Fcc[3]= 1.12441;   Fc[3]= 0.87785;   Fll[3]= 1.00192; 
         cafac[3]= 0.801094;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48838.4;   wouttag[3]= 5832.46;  
       fbb[3]= 0.0553975;   fcc[3]= 0.10794;   fc[3]= 0.139503;   fll[3]= 0.697159;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.1176;   Fcc[4]= 1.1176;   Fc[4]= 0.872539;   Fll[4]= 0.995863; 
         cafac[4]= 0.793471;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10722.6;   wouttag[4]= 1852.95;  
       fbb[4]= 0.0767186;   fcc[4]= 0.127934;   fc[4]= 0.128458;   fll[4]= 0.666889;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.1103;   Fcc[5]= 1.1103;   Fc[5]= 0.866836;   Fll[5]= 0.989354; 
         cafac[5]= 0.728076;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2728.34;   wouttag[5]= 567.279;  
       fbb[5]= 0.0942562;   fcc[5]= 0.149158;   fc[5]= 0.112274;   fll[5]= 0.644312;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12254;   Fcc[6]= 1.12254;   Fc[6]= 0.876395;   Fll[6]= 1.00026; 
         cafac[6]= 0.795519;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 62229.9;   wouttag[6]= 8274.88;  
       fbb[6]= 0.0609791;   fcc[6]= 0.113406;   fc[6]= 0.136268;   fll[6]= 0.689347;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.11601;   Fcc[7]= 1.11601;   Fc[7]= 0.871295;   Fll[7]= 0.994443; 
         cafac[7]= 0.777464;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13419.7;   wouttag[7]= 2421.86;  
       fbb[7]= 0.0805456;   fcc[7]= 0.132566;   fc[7]= 0.124926;   fll[7]= 0.661962;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==84   || sysname=="chmis_down" ) { 
      Fbb[1]= 1.13522;   Fcc[1]= 1.13522;   Fc[1]= 0.886292;   Fll[1]= 1.01156; 
         cafac[1]= 0.966309;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.07401e+06;   wouttag[1]= 38161.3;  
       fbb[1]= 0.0135665;   fcc[1]= 0.0396548;   fc[1]= 0.122801;   fll[1]= 0.823978;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.1306;   Fcc[2]= 1.1306;   Fc[2]= 0.882688;   Fll[2]= 1.00745; 
         cafac[2]= 0.876042;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 235722;   wouttag[2]= 18272.4;  
       fbb[2]= 0.033135;   fcc[2]= 0.0799415;   fc[2]= 0.13983;   fll[2]= 0.747093;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.12441;   Fcc[3]= 1.12441;   Fc[3]= 0.87785;   Fll[3]= 1.00192; 
         cafac[3]= 0.789957;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48159.4;   wouttag[3]= 5751.38;  
       fbb[3]= 0.0553975;   fcc[3]= 0.10794;   fc[3]= 0.139503;   fll[3]= 0.697159;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.1176;   Fcc[4]= 1.1176;   Fc[4]= 0.872539;   Fll[4]= 0.995863; 
         cafac[4]= 0.782439;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10573.5;   wouttag[4]= 1827.19;  
       fbb[4]= 0.0767186;   fcc[4]= 0.127934;   fc[4]= 0.128458;   fll[4]= 0.666889;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.1103;   Fcc[5]= 1.1103;   Fc[5]= 0.866836;   Fll[5]= 0.989354; 
         cafac[5]= 0.71652;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2685.03;   wouttag[5]= 558.275;  
       fbb[5]= 0.0942562;   fcc[5]= 0.149158;   fc[5]= 0.112274;   fll[5]= 0.644312;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12254;   Fcc[6]= 1.12254;   Fc[6]= 0.876395;   Fll[6]= 1.00026; 
         cafac[6]= 0.782891;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61242.1;   wouttag[6]= 8143.53;  
       fbb[6]= 0.0609791;   fcc[6]= 0.113406;   fc[6]= 0.136268;   fll[6]= 0.689347;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.11601;   Fcc[7]= 1.11601;   Fc[7]= 0.871295;   Fll[7]= 0.994443; 
         cafac[7]= 0.765123;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13206.7;   wouttag[7]= 2383.42;  
       fbb[7]= 0.0805456;   fcc[7]= 0.132566;   fc[7]= 0.124926;   fll[7]= 0.661962;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==67   || sysname=="eer_up" ) { 
      Fbb[1]= 1.13514;   Fcc[1]= 1.13514;   Fc[1]= 0.886769;   Fll[1]= 1.01148; 
         cafac[1]= 0.971972;  
         winpretag[1]= 1.1112e+06 ;   wintag[1]= 41227.1 ;   woutpretag[1]= 1.08006e+06;   wouttag[1]= 38378.5;  
       fbb[1]= 0.0135702;   fcc[1]= 0.0396389;   fc[1]= 0.122868;   fll[1]= 0.823923;  
       fmcbb[1]= 0.0119546;   fmccc[1]= 0.0349196;   fmcc[1]= 0.138557;   fmcll[1]= 0.814569;  
      Fbb[2]= 1.13052;   Fcc[2]= 1.13052;   Fc[2]= 0.883157;   Fll[2]= 1.00736; 
         cafac[2]= 0.880775;  
         winpretag[2]= 269009 ;   wintag[2]= 21003.5 ;   woutpretag[2]= 236937;   wouttag[2]= 18364.6;  
       fbb[2]= 0.0331357;   fcc[2]= 0.0799301;   fc[2]= 0.139934;   fll[2]= 0.747;  
       fmcbb[2]= 0.0293101;   fmccc[2]= 0.070702;   fmcc[2]= 0.158448;   fmcll[2]= 0.74154;  
      Fbb[3]= 1.12434;   Fcc[3]= 1.12434;   Fc[3]= 0.878325;   Fll[3]= 1.00185; 
         cafac[3]= 0.795859;  
         winpretag[3]= 60949.5 ;   wintag[3]= 7185.13 ;   woutpretag[3]= 48507.3;   wouttag[3]= 5789.29;  
       fbb[3]= 0.055363;   fcc[3]= 0.10787;   fc[3]= 0.139606;   fll[3]= 0.697161;  
       fmcbb[3]= 0.0492406;   fmccc[3]= 0.0959414;   fmcc[3]= 0.158946;   fmcll[3]= 0.695872;  
      Fbb[4]= 1.11749;   Fcc[4]= 1.11749;   Fc[4]= 0.872975;   Fll[4]= 0.995749; 
         cafac[4]= 0.787298;  
         winpretag[4]= 13512.2 ;   wintag[4]= 2283.48 ;   woutpretag[4]= 10638.1;   wouttag[4]= 1837.74;  
       fbb[4]= 0.0767679;   fcc[4]= 0.12793;   fc[4]= 0.128334;   fll[4]= 0.666968;  
       fmcbb[4]= 0.0686969;   fmccc[4]= 0.11448;   fmcc[4]= 0.147008;   fmcll[4]= 0.669815;  
      Fbb[5]= 1.11025;   Fcc[5]= 1.11025;   Fc[5]= 0.867324;   Fll[5]= 0.989304; 
         cafac[5]= 0.723888;  
         winpretag[5]= 3749.67 ;   wintag[5]= 758.329 ;   woutpretag[5]= 2714.34;   wouttag[5]= 565.053;  
       fbb[5]= 0.0941931;   fcc[5]= 0.149079;   fc[5]= 0.112388;   fll[5]= 0.64434;  
       fmcbb[5]= 0.0848392;   fmccc[5]= 0.134274;   fmcc[5]= 0.12958;   fmcll[5]= 0.651306;  
      Fbb[6]= 1.12246;   Fcc[6]= 1.12246;   Fc[6]= 0.876863;   Fll[6]= 1.00018; 
         cafac[6]= 0.789435;  
         winpretag[6]= 78211.4 ;   wintag[6]= 10226.9 ;   woutpretag[6]= 61742.8;   wouttag[6]= 8206.93;  
       fbb[6]= 0.0609596;   fcc[6]= 0.113349;   fc[6]= 0.136331;   fll[6]= 0.689361;  
       fmcbb[6]= 0.0543087;   fmccc[6]= 0.100982;   fmcc[6]= 0.155476;   fmcll[6]= 0.689234;  
      Fbb[7]= 1.11591;   Fcc[7]= 1.11591;   Fc[7]= 0.871741;   Fll[7]= 0.994342; 
         cafac[7]= 0.771237;  
         winpretag[7]= 17261.8 ;   wintag[7]= 3041.81 ;   woutpretag[7]= 13313;   wouttag[7]= 2402.78;  
       fbb[7]= 0.0805723;   fcc[7]= 0.132547;   fc[7]= 0.124853;   fll[7]= 0.662027;  
       fmcbb[7]= 0.0722034;   fmccc[7]= 0.11878;   fmcc[7]= 0.143222;   fmcll[7]= 0.665794;  
 }  
 else if ( idsys==68   || sysname=="eer_down" ) { 
      Fbb[1]= 1.14719;   Fcc[1]= 1.14719;   Fc[1]= 0.873593;   Fll[1]= 1.01302; 
         cafac[1]= 0.967929;  
         winpretag[1]= 1.11174e+06 ;   wintag[1]= 41283.8 ;   woutpretag[1]= 1.07609e+06;   wouttag[1]= 38052.4;  
       fbb[1]= 0.0137317;   fcc[1]= 0.0400448;   fc[1]= 0.120993;   fll[1]= 0.825231;  
       fmcbb[1]= 0.0119699;   fmccc[1]= 0.0349069;   fmcc[1]= 0.1385;   fmcll[1]= 0.814623;  
      Fbb[2]= 1.14219;   Fcc[2]= 1.14219;   Fc[2]= 0.869786;   Fll[2]= 1.00861; 
         cafac[2]= 0.878848;  
         winpretag[2]= 269149 ;   wintag[2]= 21093.5 ;   woutpretag[2]= 236541;   wouttag[2]= 18370.4;  
       fbb[2]= 0.033592;   fcc[2]= 0.0807829;   fc[2]= 0.137736;   fll[2]= 0.747889;  
       fmcbb[2]= 0.0294102;   fmccc[2]= 0.0707264;   fmcc[2]= 0.158356;   fmcll[2]= 0.741507;  
      Fbb[3]= 1.13547;   Fcc[3]= 1.13547;   Fc[3]= 0.864668;   Fll[3]= 1.00267; 
         cafac[3]= 0.789923;  
         winpretag[3]= 60961.4 ;   wintag[3]= 7190.1 ;   woutpretag[3]= 48154.8;   wouttag[3]= 5752.3;  
       fbb[3]= 0.0558665;   fcc[3]= 0.108717;   fc[3]= 0.137344;   fll[3]= 0.698072;  
       fmcbb[3]= 0.0492013;   fmccc[3]= 0.0957467;   fmcc[3]= 0.15884;   fmcll[3]= 0.696212;  
      Fbb[4]= 1.12799;   Fcc[4]= 1.12799;   Fc[4]= 0.858971;   Fll[4]= 0.996066; 
         cafac[4]= 0.783069;  
         winpretag[4]= 13490 ;   wintag[4]= 2275.95 ;   woutpretag[4]= 10563.6;   wouttag[4]= 1824.17;  
       fbb[4]= 0.077358;   fcc[4]= 0.128741;   fc[4]= 0.126372;   fll[4]= 0.667529;  
       fmcbb[4]= 0.0685806;   fmccc[4]= 0.114134;   fmcc[4]= 0.14712;   fmcll[4]= 0.670165;  
      Fbb[5]= 1.11966;   Fcc[5]= 1.11966;   Fc[5]= 0.852634;   Fll[5]= 0.988717; 
         cafac[5]= 0.719237;  
         winpretag[5]= 3742.57 ;   wintag[5]= 759.788 ;   woutpretag[5]= 2691.8;   wouttag[5]= 563.669;  
       fbb[5]= 0.0949357;   fcc[5]= 0.151582;   fc[5]= 0.109948;   fll[5]= 0.643534;  
       fmcbb[5]= 0.0847894;   fmccc[5]= 0.135382;   fmcc[5]= 0.128951;   fmcll[5]= 0.650878;  
      Fbb[6]= 1.1334;   Fcc[6]= 1.1334;   Fc[6]= 0.863098;   Fll[6]= 1.00085; 
         cafac[6]= 0.783885;  
         winpretag[6]= 78194 ;   wintag[6]= 10225.8 ;   woutpretag[6]= 61295.2;   wouttag[6]= 8154.69;  
       fbb[6]= 0.0614849;   fcc[6]= 0.114265;   fc[6]= 0.134115;   fll[6]= 0.690135;  
       fmcbb[6]= 0.0542479;   fmccc[6]= 0.100816;   fmcc[6]= 0.155388;   fmcll[6]= 0.689549;  
      Fbb[7]= 1.12617;   Fcc[7]= 1.12617;   Fc[7]= 0.857587;   Fll[7]= 0.99446; 
         cafac[7]= 0.766877;  
         winpretag[7]= 17232.6 ;   wintag[7]= 3035.74 ;   woutpretag[7]= 13215.3;   wouttag[7]= 2388.07;  
       fbb[7]= 0.0811977;   fcc[7]= 0.133731;   fc[7]= 0.122785;   fll[7]= 0.662287;  
       fmcbb[7]= 0.0721008;   fmccc[7]= 0.118748;   fmcc[7]= 0.143175;   fmcll[7]= 0.665976;  
 }  
 else if ( idsys==65   || sysname=="ees_up" ) { 
      Fbb[1]= 1.13497;   Fcc[1]= 1.13497;   Fc[1]= 0.886405;   Fll[1]= 1.01155; 
         cafac[1]= 0.972097;  
         winpretag[1]= 1.11147e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08045e+06;   wouttag[1]= 38391;  
       fbb[1]= 0.0135633;   fcc[1]= 0.0396456;   fc[1]= 0.122816;   fll[1]= 0.823976;  
       fmcbb[1]= 0.0119504;   fmccc[1]= 0.0349311;   fmcc[1]= 0.138555;   fmcll[1]= 0.814564;  
      Fbb[2]= 1.13036;   Fcc[2]= 1.13036;   Fc[2]= 0.88281;   Fll[2]= 1.00745; 
         cafac[2]= 0.880465;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236912;   wouttag[2]= 18364.3;  
       fbb[2]= 0.033128;   fcc[2]= 0.0799247;   fc[2]= 0.13985;   fll[2]= 0.747098;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.0707071;   fmcc[2]= 0.158414;   fmcll[2]= 0.741571;  
      Fbb[3]= 1.12418;   Fcc[3]= 1.12418;   Fc[3]= 0.877982;   Fll[3]= 1.00194; 
         cafac[3]= 0.795558;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48500.9;   wouttag[3]= 5791.89;  
       fbb[3]= 0.0553864;   fcc[3]= 0.107919;   fc[3]= 0.139524;   fll[3]= 0.697171;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.11739;   Fcc[4]= 1.11739;   Fc[4]= 0.872679;   Fll[4]= 0.995891; 
         cafac[4]= 0.78806;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10649.5;   wouttag[4]= 1840.22;  
       fbb[4]= 0.076704;   fcc[4]= 0.12791;   fc[4]= 0.128478;   fll[4]= 0.666907;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.1101;   Fcc[5]= 1.1101;   Fc[5]= 0.866984;   Fll[5]= 0.989392; 
         cafac[5]= 0.722318;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.76;   wouttag[5]= 562.757;  
       fbb[5]= 0.0942393;   fcc[5]= 0.149131;   fc[5]= 0.112293;   fll[5]= 0.644337;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12232;   Fcc[6]= 1.12232;   Fc[6]= 0.876529;   Fll[6]= 1.00028; 
         cafac[6]= 0.789251;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61739.5;   wouttag[6]= 8209.27;  
       fbb[6]= 0.060967;   fcc[6]= 0.113383;   fc[6]= 0.136289;   fll[6]= 0.689361;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.1158;   Fcc[7]= 1.1158;   Fc[7]= 0.871437;   Fll[7]= 0.994473; 
         cafac[7]= 0.771377;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13314.6;   wouttag[7]= 2402.76;  
       fbb[7]= 0.0805304;   fcc[7]= 0.132541;   fc[7]= 0.124947;   fll[7]= 0.661982;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==66   || sysname=="ees_down" ) { 
      Fbb[1]= 1.15989;   Fcc[1]= 1.15989;   Fc[1]= 0.873837;   Fll[1]= 1.01229; 
         cafac[1]= 0.974705;  
         winpretag[1]= 1.1067e+06 ;   wintag[1]= 41072.2 ;   woutpretag[1]= 1.07871e+06;   wouttag[1]= 38187.7;  
       fbb[1]= 0.0138486;   fcc[1]= 0.0405201;   fc[1]= 0.121211;   fll[1]= 0.82442;  
       fmcbb[1]= 0.0119395;   fmccc[1]= 0.0349344;   fmcc[1]= 0.138712;   fmcll[1]= 0.814414;  
      Fbb[2]= 1.15401;   Fcc[2]= 1.15401;   Fc[2]= 0.869406;   Fll[2]= 1.00715; 
         cafac[2]= 0.882115;  
         winpretag[2]= 268197 ;   wintag[2]= 20985.9 ;   woutpretag[2]= 236580;   wouttag[2]= 18399.7;  
       fbb[2]= 0.0338983;   fcc[2]= 0.0815866;   fc[2]= 0.137906;   fll[2]= 0.746609;  
       fmcbb[2]= 0.0293744;   fmccc[2]= 0.0706984;   fmcc[2]= 0.158621;   fmcll[2]= 0.741306;  
      Fbb[3]= 1.14646;   Fcc[3]= 1.14646;   Fc[3]= 0.863715;   Fll[3]= 1.00056; 
         cafac[3]= 0.795041;  
         winpretag[3]= 60727.5 ;   wintag[3]= 7176.47 ;   woutpretag[3]= 48280.9;   wouttag[3]= 5799.76;  
       fbb[3]= 0.0564944;   fcc[3]= 0.109919;   fc[3]= 0.137201;   fll[3]= 0.696385;  
       fmcbb[3]= 0.0492774;   fmccc[3]= 0.0958773;   fmcc[3]= 0.15885;   fmcll[3]= 0.695995;  
      Fbb[4]= 1.13842;   Fcc[4]= 1.13842;   Fc[4]= 0.857661;   Fll[4]= 0.993547; 
         cafac[4]= 0.786089;  
         winpretag[4]= 13488.7 ;   wintag[4]= 2279.22 ;   woutpretag[4]= 10603.3;   wouttag[4]= 1840.12;  
       fbb[4]= 0.078283;   fcc[4]= 0.129851;   fc[4]= 0.126442;   fll[4]= 0.665424;  
       fmcbb[4]= 0.0687646;   fmccc[4]= 0.114062;   fmcc[4]= 0.147427;   fmcll[4]= 0.669746;  
      Fbb[5]= 1.1296;   Fcc[5]= 1.1296;   Fc[5]= 0.851017;   Fll[5]= 0.98585; 
         cafac[5]= 0.733294;  
         winpretag[5]= 3730.21 ;   wintag[5]= 757.239 ;   woutpretag[5]= 2735.34;   wouttag[5]= 575.225;  
       fbb[5]= 0.0967154;   fcc[5]= 0.151429;   fc[5]= 0.110001;   fll[5]= 0.641855;  
       fmcbb[5]= 0.0856191;   fmccc[5]= 0.134055;   fmcc[5]= 0.129258;   fmcll[5]= 0.651067;  
      Fbb[6]= 1.14424;   Fcc[6]= 1.14424;   Fc[6]= 0.862047;   Fll[6]= 0.998627; 
         cafac[6]= 0.78919;  
         winpretag[6]= 77946.4 ;   wintag[6]= 10212.9 ;   woutpretag[6]= 61514.5;   wouttag[6]= 8229.87;  
       fbb[6]= 0.062234;   fcc[6]= 0.115398;   fc[6]= 0.134011;   fll[6]= 0.688357;  
       fmcbb[6]= 0.0543889;   fmccc[6]= 0.100851;   fmcc[6]= 0.155457;   fmcll[6]= 0.689303;  
      Fbb[7]= 1.1365;   Fcc[7]= 1.1365;   Fc[7]= 0.856213;   Fll[7]= 0.991869; 
         cafac[7]= 0.772783;  
         winpretag[7]= 17218.9 ;   wintag[7]= 3036.45 ;   woutpretag[7]= 13306.5;   wouttag[7]= 2415.83;  
       fbb[7]= 0.0823004;   fcc[7]= 0.134554;   fc[7]= 0.122859;   fll[7]= 0.660287;  
       fmcbb[7]= 0.0724159;   fmccc[7]= 0.118394;   fmcc[7]= 0.143491;   fmcll[7]= 0.6657;  
 }  
 else if ( idsys==63   || sysname=="jer_one" ) { 
      Fbb[1]= 1.20212;   Fcc[1]= 1.20212;   Fc[1]= 0.859135;   Fll[1]= 1.01198; 
         cafac[1]= 0.90827;  
         winpretag[1]= 1.17816e+06 ;   wintag[1]= 42400.7 ;   woutpretag[1]= 1.07009e+06;   wouttag[1]= 36687.8;  
       fbb[1]= 0.0137586;   fcc[1]= 0.0399002;   fc[1]= 0.115038;   fll[1]= 0.831303;  
       fmcbb[1]= 0.0114453;   fmccc[1]= 0.0331917;   fmcc[1]= 0.1339;   fmcll[1]= 0.821463;  
      Fbb[2]= 1.19468;   Fcc[2]= 1.19468;   Fc[2]= 0.853817;   Fll[2]= 1.00571; 
         cafac[2]= 0.851333;  
         winpretag[2]= 281373 ;   wintag[2]= 21415.6 ;   woutpretag[2]= 239542;   wouttag[2]= 18198.6;  
       fbb[2]= 0.0336685;   fcc[2]= 0.082451;   fc[2]= 0.135362;   fll[2]= 0.748518;  
       fmcbb[2]= 0.0281821;   fmccc[2]= 0.0690154;   fmcc[2]= 0.158537;   fmcll[2]= 0.744265;  
      Fbb[3]= 1.18505;   Fcc[3]= 1.18505;   Fc[3]= 0.846939;   Fll[3]= 0.997613; 
         cafac[3]= 0.77626;  
         winpretag[3]= 64542.9 ;   wintag[3]= 7558.08 ;   woutpretag[3]= 50102;   wouttag[3]= 6002.69;  
       fbb[3]= 0.0571256;   fcc[3]= 0.111581;   fc[3]= 0.136573;   fll[3]= 0.69472;  
       fmcbb[3]= 0.0482052;   fmccc[3]= 0.0941572;   fmcc[3]= 0.161255;   fmcll[3]= 0.696382;  
      Fbb[4]= 1.17509;   Fcc[4]= 1.17509;   Fc[4]= 0.839822;   Fll[4]= 0.98923; 
         cafac[4]= 0.742838;  
         winpretag[4]= 14893.7 ;   wintag[4]= 2546.98 ;   woutpretag[4]= 11063.6;   wouttag[4]= 1957.23;  
       fbb[4]= 0.0789482;   fcc[4]= 0.129863;   fc[4]= 0.125111;   fll[4]= 0.666078;  
       fmcbb[4]= 0.0671846;   fmccc[4]= 0.110513;   fmcc[4]= 0.148974;   fmcll[4]= 0.673329;  
      Fbb[5]= 1.16357;   Fcc[5]= 1.16357;   Fc[5]= 0.831588;   Fll[5]= 0.979531; 
         cafac[5]= 0.679695;  
         winpretag[5]= 4111.66 ;   wintag[5]= 819.876 ;   woutpretag[5]= 2794.67;   wouttag[5]= 581.229;  
       fbb[5]= 0.100207;   fcc[5]= 0.156678;   fc[5]= 0.113334;   fll[5]= 0.629781;  
       fmcbb[5]= 0.0861199;   fmccc[5]= 0.134653;   fmcc[5]= 0.136286;   fmcll[5]= 0.642942;  
      Fbb[6]= 1.18219;   Fcc[6]= 1.18219;   Fc[6]= 0.844895;   Fll[6]= 0.995205; 
         cafac[6]= 0.763359;  
         winpretag[6]= 83548.2 ;   wintag[6]= 10924.9 ;   woutpretag[6]= 63777.3;   wouttag[6]= 8575.35;  
       fbb[6]= 0.0631934;   fcc[6]= 0.117115;   fc[6]= 0.133356;   fll[6]= 0.686336;  
       fmcbb[6]= 0.0534544;   fmccc[6]= 0.0990657;   fmcc[6]= 0.157837;   fmcll[6]= 0.689643;  
      Fbb[7]= 1.17258;   Fcc[7]= 1.17258;   Fc[7]= 0.838027;   Fll[7]= 0.987116; 
         cafac[7]= 0.727216;  
         winpretag[7]= 19005.3 ;   wintag[7]= 3366.86 ;   woutpretag[7]= 13821;   wouttag[7]= 2538.66;  
       fbb[7]= 0.0835829;   fcc[7]= 0.135709;   fc[7]= 0.122544;   fll[7]= 0.658164;  
       fmcbb[7]= 0.0712811;   fmccc[7]= 0.115735;   fmcc[7]= 0.146229;   fmcll[7]= 0.666755;  
 }  
 else if ( idsys==64   || sysname=="jef_one" ) { 
      Fbb[1]= 1.10799;   Fcc[1]= 1.10799;   Fc[1]= 0.889769;   Fll[1]= 1.01255; 
         cafac[1]= 0.971864;  
         winpretag[1]= 1.11279e+06 ;   wintag[1]= 41279.1 ;   woutpretag[1]= 1.08148e+06;   wouttag[1]= 38343.1;  
       fbb[1]= 0.0132693;   fcc[1]= 0.0387244;   fc[1]= 0.123389;   fll[1]= 0.824617;  
       fmcbb[1]= 0.011976;   fmccc[1]= 0.0349502;   fmcc[1]= 0.138675;   fmcll[1]= 0.814399;  
      Fbb[2]= 1.10507;   Fcc[2]= 1.10507;   Fc[2]= 0.887421;   Fll[2]= 1.00988; 
         cafac[2]= 0.880045;  
         winpretag[2]= 269282 ;   wintag[2]= 21066.4 ;   woutpretag[2]= 236981;   wouttag[2]= 18306.8;  
       fbb[2]= 0.0324953;   fcc[2]= 0.0782296;   fc[2]= 0.140692;   fll[2]= 0.748583;  
       fmcbb[2]= 0.0294057;   fmccc[2]= 0.0707918;   fmcc[2]= 0.158541;   fmcll[2]= 0.741262;  
      Fbb[3]= 1.10041;   Fcc[3]= 1.10041;   Fc[3]= 0.883686;   Fll[3]= 1.00563; 
         cafac[3]= 0.800789;  
         winpretag[3]= 60925.5 ;   wintag[3]= 7219.89 ;   woutpretag[3]= 48788.5;   wouttag[3]= 5811.75;  
       fbb[3]= 0.0540866;   fcc[3]= 0.105639;   fc[3]= 0.140474;   fll[3]= 0.6998;  
       fmcbb[3]= 0.0491511;   fmccc[3]= 0.0959996;   fmcc[3]= 0.158964;   fmcll[3]= 0.695885;  
      Fbb[4]= 1.09483;   Fcc[4]= 1.09483;   Fc[4]= 0.879201;   Fll[4]= 1.00052; 
         cafac[4]= 0.798089;  
         winpretag[4]= 13538.8 ;   wintag[4]= 2294.47 ;   woutpretag[4]= 10805.2;   wouttag[4]= 1858.63;  
       fbb[4]= 0.0755976;   fcc[4]= 0.125148;   fc[4]= 0.129097;   fll[4]= 0.670157;  
       fmcbb[4]= 0.0690496;   fmccc[4]= 0.114308;   fmcc[4]= 0.146835;   fmcll[4]= 0.669808;  
      Fbb[5]= 1.08882;   Fcc[5]= 1.08882;   Fc[5]= 0.874379;   Fll[5]= 0.995034; 
         cafac[5]= 0.704839;  
         winpretag[5]= 3765.78 ;   wintag[5]= 764.142 ;   woutpretag[5]= 2654.27;   wouttag[5]= 550.342;  
       fbb[5]= 0.093091;   fcc[5]= 0.145635;   fc[5]= 0.113036;   fll[5]= 0.648238;  
       fmcbb[5]= 0.0854968;   fmccc[5]= 0.133754;   fmcc[5]= 0.129276;   fmcll[5]= 0.651473;  
      Fbb[6]= 1.09888;   Fcc[6]= 1.09888;   Fc[6]= 0.882455;   Fll[6]= 1.00422; 
         cafac[6]= 0.793853;  
         winpretag[6]= 78230.1 ;   wintag[6]= 10278.5 ;   woutpretag[6]= 62103.2;   wouttag[6]= 8234.55;  
       fbb[6]= 0.059718;   fcc[6]= 0.110971;   fc[6]= 0.137165;   fll[6]= 0.692146;  
       fmcbb[6]= 0.0543444;   fmccc[6]= 0.100986;   fmcc[6]= 0.155436;   fmcll[6]= 0.689234;  
      Fbb[7]= 1.09352;   Fcc[7]= 1.09352;   Fc[7]= 0.878148;   Fll[7]= 0.999323; 
         cafac[7]= 0.774047;  
         winpretag[7]= 17304.6 ;   wintag[7]= 3058.62 ;   woutpretag[7]= 13394.6;   wouttag[7]= 2407.46;  
       fbb[7]= 0.0794208;   fcc[7]= 0.129625;   fc[7]= 0.125587;   fll[7]= 0.665367;  
       fmcbb[7]= 0.0726288;   fmccc[7]= 0.11854;   fmcc[7]= 0.143014;   fmcll[7]= 0.665818;  
 }  
 else if ( idsys==61   || sysname=="musms_up" ) { 
      Fbb[1]= 1.13604;   Fcc[1]= 1.13604;   Fc[1]= 0.88604;   Fll[1]= 1.01155; 
         cafac[1]= 0.972022;  
         winpretag[1]= 1.11147e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08037e+06;   wouttag[1]= 38385.8;  
       fbb[1]= 0.0135761;   fcc[1]= 0.0396831;   fc[1]= 0.122764;   fll[1]= 0.823976;  
       fmcbb[1]= 0.0119504;   fmccc[1]= 0.0349311;   fmcc[1]= 0.138554;   fmcll[1]= 0.814564;  
      Fbb[2]= 1.13137;   Fcc[2]= 1.13137;   Fc[2]= 0.882403;   Fll[2]= 1.0074; 
         cafac[2]= 0.880423;  
         winpretag[2]= 269073 ;   wintag[2]= 21012.1 ;   woutpretag[2]= 236898;   wouttag[2]= 18364.6;  
       fbb[2]= 0.0331503;   fcc[2]= 0.0799972;   fc[2]= 0.139782;   fll[2]= 0.747071;  
       fmcbb[2]= 0.0293009;   fmccc[2]= 0.070708;   fmcc[2]= 0.158411;   fmcll[2]= 0.741581;  
      Fbb[3]= 1.12513;   Fcc[3]= 1.12513;   Fc[3]= 0.877536;   Fll[3]= 1.00185; 
         cafac[3]= 0.795432;  
         winpretag[3]= 60964.1 ;   wintag[3]= 7191.31 ;   woutpretag[3]= 48492.8;   wouttag[3]= 5791.72;  
       fbb[3]= 0.0554339;   fcc[3]= 0.108002;   fc[3]= 0.139454;   fll[3]= 0.69711;  
       fmcbb[3]= 0.0492687;   fmccc[3]= 0.09599;   fmcc[3]= 0.158916;   fmcll[3]= 0.695825;  
      Fbb[4]= 1.11829;   Fcc[4]= 1.11829;   Fc[4]= 0.872197;   Fll[4]= 0.995751; 
         cafac[4]= 0.787994;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.6;   wouttag[4]= 1840.51;  
       fbb[4]= 0.0767657;   fcc[4]= 0.128013;   fc[4]= 0.128407;   fll[4]= 0.666814;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.11094;   Fcc[5]= 1.11094;   Fc[5]= 0.866467;   Fll[5]= 0.98921; 
         cafac[5]= 0.72225;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.51;   wouttag[5]= 562.861;  
       fbb[5]= 0.094311;   fcc[5]= 0.149244;   fc[5]= 0.112226;   fll[5]= 0.644218;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12326;   Fcc[6]= 1.12326;   Fc[6]= 0.876073;   Fll[6]= 1.00018; 
         cafac[6]= 0.789139;  
         winpretag[6]= 78225 ;   wintag[6]= 10233 ;   woutpretag[6]= 61730.4;   wouttag[6]= 8209.5;  
       fbb[6]= 0.0610184;   fcc[6]= 0.113472;   fc[6]= 0.136219;   fll[6]= 0.689291;  
       fmcbb[6]= 0.0543226;   fmccc[6]= 0.10102;   fmcc[6]= 0.155488;   fmcll[6]= 0.689169;  
      Fbb[7]= 1.11669;   Fcc[7]= 1.11669;   Fc[7]= 0.870947;   Fll[7]= 0.994324; 
         cafac[7]= 0.77131;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13313.5;   wouttag[7]= 2403.16;  
       fbb[7]= 0.0805944;   fcc[7]= 0.132646;   fc[7]= 0.124876;   fll[7]= 0.661883;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==62   || sysname=="musms_down" ) { 
      Fbb[1]= 1.13588;   Fcc[1]= 1.13588;   Fc[1]= 0.886052;   Fll[1]= 1.01156; 
         cafac[1]= 0.972027;  
         winpretag[1]= 1.11147e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08038e+06;   wouttag[1]= 38385.4;  
       fbb[1]= 0.0135742;   fcc[1]= 0.0396776;   fc[1]= 0.122766;   fll[1]= 0.823982;  
       fmcbb[1]= 0.0119504;   fmccc[1]= 0.0349311;   fmcc[1]= 0.138555;   fmcll[1]= 0.814564;  
      Fbb[2]= 1.13123;   Fcc[2]= 1.13123;   Fc[2]= 0.882422;   Fll[2]= 1.00742; 
         cafac[2]= 0.880384;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236890;   wouttag[2]= 18364.3;  
       fbb[2]= 0.0331534;   fcc[2]= 0.0799859;   fc[2]= 0.139788;   fll[2]= 0.747072;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.0707071;   fmcc[2]= 0.158414;   fmcll[2]= 0.741571;  
      Fbb[3]= 1.125;   Fcc[3]= 1.125;   Fc[3]= 0.877561;   Fll[3]= 1.00187; 
         cafac[3]= 0.795465;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48495.3;   wouttag[3]= 5792.25;  
       fbb[3]= 0.0554266;   fcc[3]= 0.107997;   fc[3]= 0.139457;   fll[3]= 0.697119;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.11816;   Fcc[4]= 1.11816;   Fc[4]= 0.872228;   Fll[4]= 0.99578; 
         cafac[4]= 0.788003;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.7;   wouttag[4]= 1840.46;  
       fbb[4]= 0.0767568;   fcc[4]= 0.127998;   fc[4]= 0.128412;   fll[4]= 0.666833;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.11082;   Fcc[5]= 1.11082;   Fc[5]= 0.866502;   Fll[5]= 0.989242; 
         cafac[5]= 0.722651;  
         winpretag[5]= 3746.74 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2707.58;   wouttag[5]= 563.147;  
       fbb[5]= 0.0943152;   fcc[5]= 0.149251;   fc[5]= 0.112248;   fll[5]= 0.644186;  
       fmcbb[5]= 0.084906;   fmccc[5]= 0.134361;   fmcc[5]= 0.129542;   fmcll[5]= 0.651191;  
      Fbb[6]= 1.12312;   Fcc[6]= 1.12312;   Fc[6]= 0.8761;   Fll[6]= 1.0002; 
         cafac[6]= 0.789193;  
         winpretag[6]= 78224.9 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61734.6;   wouttag[6]= 8210.24;  
       fbb[6]= 0.0610111;   fcc[6]= 0.113465;   fc[6]= 0.136223;   fll[6]= 0.6893;  
       fmcbb[6]= 0.0543227;   fmccc[6]= 0.101027;   fmcc[6]= 0.155488;   fmcll[6]= 0.689163;  
      Fbb[7]= 1.11656;   Fcc[7]= 1.11656;   Fc[7]= 0.870979;   Fll[7]= 0.994353; 
         cafac[7]= 0.771425;  
         winpretag[7]= 17260.3 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13315;   wouttag[7]= 2403.42;  
       fbb[7]= 0.0805879;   fcc[7]= 0.132635;   fc[7]= 0.124885;   fll[7]= 0.661891;  
       fmcbb[7]= 0.0721753;   fmccc[7]= 0.11879;   fmcc[7]= 0.143385;   fmcll[7]= 0.66565;  
 }  
 else if ( idsys==59   || sysname=="musid_up" ) { 
      Fbb[1]= 1.1351;   Fcc[1]= 1.1351;   Fc[1]= 0.886398;   Fll[1]= 1.01155; 
         cafac[1]= 0.972097;  
         winpretag[1]= 1.11147e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08045e+06;   wouttag[1]= 38391.5;  
       fbb[1]= 0.0135649;   fcc[1]= 0.0396502;   fc[1]= 0.122815;   fll[1]= 0.82397;  
       fmcbb[1]= 0.0119504;   fmccc[1]= 0.0349311;   fmcc[1]= 0.138555;   fmcll[1]= 0.814564;  
      Fbb[2]= 1.13049;   Fcc[2]= 1.13049;   Fc[2]= 0.882797;   Fll[2]= 1.00744; 
         cafac[2]= 0.880462;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236911;   wouttag[2]= 18364.7;  
       fbb[2]= 0.0331316;   fcc[2]= 0.0799334;   fc[2]= 0.139847;   fll[2]= 0.747088;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.0707071;   fmcc[2]= 0.158414;   fmcll[2]= 0.741571;  
      Fbb[3]= 1.1243;   Fcc[3]= 1.1243;   Fc[3]= 0.877964;   Fll[3]= 1.00192; 
         cafac[3]= 0.795507;  
         winpretag[3]= 60964.1 ;   wintag[3]= 7191.31 ;   woutpretag[3]= 48497.4;   wouttag[3]= 5791.2;  
       fbb[3]= 0.0553926;   fcc[3]= 0.107921;   fc[3]= 0.139523;   fll[3]= 0.697164;  
       fmcbb[3]= 0.0492687;   fmccc[3]= 0.09599;   fmcc[3]= 0.158916;   fmcll[3]= 0.695825;  
      Fbb[4]= 1.1175;   Fcc[4]= 1.1175;   Fc[4]= 0.872656;   Fll[4]= 0.995866; 
         cafac[4]= 0.788058;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10649.4;   wouttag[4]= 1840.28;  
       fbb[4]= 0.0767115;   fcc[4]= 0.127923;   fc[4]= 0.128475;   fll[4]= 0.666891;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.1102;   Fcc[5]= 1.1102;   Fc[5]= 0.866955;   Fll[5]= 0.989361; 
         cafac[5]= 0.722712;  
         winpretag[5]= 3746.74 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2707.81;   wouttag[5]= 563.085;  
       fbb[5]= 0.0942626;   fcc[5]= 0.149168;   fc[5]= 0.112307;   fll[5]= 0.644263;  
       fmcbb[5]= 0.084906;   fmccc[5]= 0.134361;   fmcc[5]= 0.129542;   fmcll[5]= 0.651191;  
      Fbb[6]= 1.12243;   Fcc[6]= 1.12243;   Fc[6]= 0.87651;   Fll[6]= 1.00026; 
         cafac[6]= 0.78924;  
         winpretag[6]= 78224.4 ;   wintag[6]= 10233 ;   woutpretag[6]= 61737.8;   wouttag[6]= 8208.94;  
       fbb[6]= 0.0609741;   fcc[6]= 0.113389;   fc[6]= 0.136288;   fll[6]= 0.689349;  
       fmcbb[6]= 0.054323;   fmccc[6]= 0.101021;   fmcc[6]= 0.155489;   fmcll[6]= 0.689167;  
      Fbb[7]= 1.11591;   Fcc[7]= 1.11591;   Fc[7]= 0.871412;   Fll[7]= 0.994447; 
         cafac[7]= 0.771482;  
         winpretag[7]= 17260.3 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13316;   wouttag[7]= 2403.17;  
       fbb[7]= 0.0805409;   fcc[7]= 0.132558;   fc[7]= 0.124947;   fll[7]= 0.661954;  
       fmcbb[7]= 0.0721753;   fmccc[7]= 0.11879;   fmcc[7]= 0.143385;   fmcll[7]= 0.66565;  
 }  
 else if ( idsys==60   || sysname=="musid_down" ) { 
      Fbb[1]= 1.13535;   Fcc[1]= 1.13535;   Fc[1]= 0.886175;   Fll[1]= 1.01157; 
         cafac[1]= 0.972053;  
         winpretag[1]= 1.11147e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.0804e+06;   wouttag[1]= 38386.4;  
       fbb[1]= 0.0135679;   fcc[1]= 0.039659;   fc[1]= 0.122784;   fll[1]= 0.823989;  
       fmcbb[1]= 0.0119504;   fmccc[1]= 0.0349311;   fmcc[1]= 0.138555;   fmcll[1]= 0.814564;  
      Fbb[2]= 1.13073;   Fcc[2]= 1.13073;   Fc[2]= 0.882569;   Fll[2]= 1.00745; 
         cafac[2]= 0.880441;  
         winpretag[2]= 269075 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236904;   wouttag[2]= 18363.9;  
       fbb[2]= 0.0331362;   fcc[2]= 0.079951;   fc[2]= 0.13981;   fll[2]= 0.747103;  
       fmcbb[2]= 0.0293052;   fmccc[2]= 0.0707075;   fmcc[2]= 0.158412;   fmcll[2]= 0.741575;  
      Fbb[3]= 1.12453;   Fcc[3]= 1.12453;   Fc[3]= 0.877727;   Fll[3]= 1.00193; 
         cafac[3]= 0.795513;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48498.1;   wouttag[3]= 5791.89;  
       fbb[3]= 0.0554034;   fcc[3]= 0.107952;   fc[3]= 0.139484;   fll[3]= 0.697161;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.11772;   Fcc[4]= 1.11772;   Fc[4]= 0.872412;   Fll[4]= 0.995861; 
         cafac[4]= 0.78802;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10648.9;   wouttag[4]= 1840.26;  
       fbb[4]= 0.0767263;   fcc[4]= 0.127947;   fc[4]= 0.128439;   fll[4]= 0.666887;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.1104;   Fcc[5]= 1.1104;   Fc[5]= 0.866704;   Fll[5]= 0.989345; 
         cafac[5]= 0.72228;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.62;   wouttag[5]= 562.778;  
       fbb[5]= 0.0942651;   fcc[5]= 0.149172;   fc[5]= 0.112257;   fll[5]= 0.644306;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12266;   Fcc[6]= 1.12266;   Fc[6]= 0.876271;   Fll[6]= 1.00027; 
         cafac[6]= 0.789206;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61736.1;   wouttag[6]= 8209.33;  
       fbb[6]= 0.0609855;   fcc[6]= 0.113418;   fc[6]= 0.136249;   fll[6]= 0.689348;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.11612;   Fcc[7]= 1.11612;   Fc[7]= 0.871166;   Fll[7]= 0.994439; 
         cafac[7]= 0.771338;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13314;   wouttag[7]= 2402.83;  
       fbb[7]= 0.0805536;   fcc[7]= 0.132579;   fc[7]= 0.124908;   fll[7]= 0.66196;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==57   || sysname=="mscale_one" ) { 
      Fbb[1]= 1.13722;   Fcc[1]= 1.13722;   Fc[1]= 0.885632;   Fll[1]= 1.01156; 
         cafac[1]= 0.971928;  
         winpretag[1]= 1.11147e+06 ;   wintag[1]= 41244 ;   woutpretag[1]= 1.08027e+06;   wouttag[1]= 38378.9;  
       fbb[1]= 0.0135897;   fcc[1]= 0.0397244;   fc[1]= 0.122708;   fll[1]= 0.823978;  
       fmcbb[1]= 0.0119499;   fmccc[1]= 0.0349311;   fmcc[1]= 0.138554;   fmcll[1]= 0.814565;  
      Fbb[2]= 1.13249;   Fcc[2]= 1.13249;   Fc[2]= 0.881949;   Fll[2]= 1.00735; 
         cafac[2]= 0.88036;  
         winpretag[2]= 269074 ;   wintag[2]= 21012.1 ;   woutpretag[2]= 236882;   wouttag[2]= 18365.8;  
       fbb[2]= 0.033183;   fcc[2]= 0.0800759;   fc[2]= 0.139712;   fll[2]= 0.747029;  
       fmcbb[2]= 0.0293008;   fmccc[2]= 0.0707076;   fmcc[2]= 0.158413;   fmcll[2]= 0.741579;  
      Fbb[3]= 1.12619;   Fcc[3]= 1.12619;   Fc[3]= 0.877038;   Fll[3]= 1.00174; 
         cafac[3]= 0.795402;  
         winpretag[3]= 60964.1 ;   wintag[3]= 7191.31 ;   woutpretag[3]= 48491;   wouttag[3]= 5792.9;  
       fbb[3]= 0.0554858;   fcc[3]= 0.108103;   fc[3]= 0.139375;   fll[3]= 0.697036;  
       fmcbb[3]= 0.0492687;   fmccc[3]= 0.09599;   fmcc[3]= 0.158916;   fmcll[3]= 0.695825;  
      Fbb[4]= 1.11928;   Fcc[4]= 1.11928;   Fc[4]= 0.87166;   Fll[4]= 0.995598; 
         cafac[4]= 0.787926;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10647.7;   wouttag[4]= 1840.84;  
       fbb[4]= 0.0768338;   fcc[4]= 0.128127;   fc[4]= 0.128328;   fll[4]= 0.666711;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.11188;   Fcc[5]= 1.11188;   Fc[5]= 0.865893;   Fll[5]= 0.98901; 
         cafac[5]= 0.722175;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.22;   wouttag[5]= 562.975;  
       fbb[5]= 0.0943902;   fcc[5]= 0.14937;   fc[5]= 0.112152;   fll[5]= 0.644088;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.1243;   Fcc[6]= 1.1243;   Fc[6]= 0.875565;   Fll[6]= 1.00006; 
         cafac[6]= 0.789097;  
         winpretag[6]= 78225 ;   wintag[6]= 10233 ;   woutpretag[6]= 61727.1;   wouttag[6]= 8211.17;  
       fbb[6]= 0.0610747;   fcc[6]= 0.113576;   fc[6]= 0.13614;   fll[6]= 0.689209;  
       fmcbb[6]= 0.0543226;   fmccc[6]= 0.10102;   fmcc[6]= 0.155488;   fmcll[6]= 0.689169;  
      Fbb[7]= 1.11767;   Fcc[7]= 1.11767;   Fc[7]= 0.870402;   Fll[7]= 0.99416; 
         cafac[7]= 0.771239;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13312.2;   wouttag[7]= 2403.61;  
       fbb[7]= 0.0806651;   fcc[7]= 0.132762;   fc[7]= 0.124798;   fll[7]= 0.661774;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==55   || sysname=="metpileup_up" ) { 
      Fbb[1]= 1.11471;   Fcc[1]= 1.11471;   Fc[1]= 0.887253;   Fll[1]= 1.01258; 
         cafac[1]= 0.971991;  
         winpretag[1]= 1.11479e+06 ;   wintag[1]= 41419.1 ;   woutpretag[1]= 1.08356e+06;   wouttag[1]= 38462.6;  
       fbb[1]= 0.0132995;   fcc[1]= 0.0389621;   fc[1]= 0.122968;   fll[1]= 0.824771;  
       fmcbb[1]= 0.0119309;   fmccc[1]= 0.0349525;   fmcc[1]= 0.138594;   fmcll[1]= 0.814523;  
      Fbb[2]= 1.11143;   Fcc[2]= 1.11143;   Fc[2]= 0.88464;   Fll[2]= 1.0096; 
         cafac[2]= 0.881878;  
         winpretag[2]= 269542 ;   wintag[2]= 21067.2 ;   woutpretag[2]= 237703;   wouttag[2]= 18358.1;  
       fbb[2]= 0.0326084;   fcc[2]= 0.0785704;   fc[2]= 0.140075;   fll[2]= 0.748746;  
       fmcbb[2]= 0.029339;   fmccc[2]= 0.070693;   fmcc[2]= 0.158341;   fmcll[2]= 0.741627;  
      Fbb[3]= 1.1064;   Fcc[3]= 1.1064;   Fc[3]= 0.880633;   Fll[3]= 1.00503; 
         cafac[3]= 0.796436;  
         winpretag[3]= 60980.1 ;   wintag[3]= 7204.3 ;   woutpretag[3]= 48566.7;   wouttag[3]= 5776.24;  
       fbb[3]= 0.0543841;   fcc[3]= 0.106146;   fc[3]= 0.139713;   fll[3]= 0.699757;  
       fmcbb[3]= 0.0491542;   fmccc[3]= 0.0959384;   fmcc[3]= 0.15865;   fmcll[3]= 0.696257;  
      Fbb[4]= 1.10066;   Fcc[4]= 1.10066;   Fc[4]= 0.876069;   Fll[4]= 0.999818; 
         cafac[4]= 0.793574;  
         winpretag[4]= 13550 ;   wintag[4]= 2282.57 ;   woutpretag[4]= 10753;   wouttag[4]= 1840.71;  
       fbb[4]= 0.0750593;   fcc[4]= 0.125637;   fc[4]= 0.128887;   fll[4]= 0.670416;  
       fmcbb[4]= 0.0681946;   fmccc[4]= 0.114146;   fmcc[4]= 0.14712;   fmcll[4]= 0.670539;  
      Fbb[5]= 1.09408;   Fcc[5]= 1.09408;   Fc[5]= 0.870832;   Fll[5]= 0.993842; 
         cafac[5]= 0.728599;  
         winpretag[5]= 3748.19 ;   wintag[5]= 759.096 ;   woutpretag[5]= 2730.92;   wouttag[5]= 565.954;  
       fbb[5]= 0.0934343;   fcc[5]= 0.147074;   fc[5]= 0.112406;   fll[5]= 0.647086;  
       fmcbb[5]= 0.0853995;   fmccc[5]= 0.134426;   fmcc[5]= 0.129078;   fmcll[5]= 0.651096;  
      Fbb[6]= 1.10481;   Fcc[6]= 1.10481;   Fc[6]= 0.879366;   Fll[6]= 1.00358; 
         cafac[6]= 0.791402;  
         winpretag[6]= 78278.3 ;   wintag[6]= 10246 ;   woutpretag[6]= 61949.6;   wouttag[6]= 8194.81;  
       fbb[6]= 0.0598646;   fcc[6]= 0.111512;   fc[6]= 0.136511;   fll[6]= 0.692112;  
       fmcbb[6]= 0.0541856;   fmccc[6]= 0.100933;   fmcc[6]= 0.155239;   fmcll[6]= 0.689643;  
      Fbb[7]= 1.09923;   Fcc[7]= 1.09923;   Fc[7]= 0.874929;   Fll[7]= 0.998517; 
         cafac[7]= 0.777104;  
         winpretag[7]= 17298.2 ;   wintag[7]= 3041.67 ;   woutpretag[7]= 13442.5;   wouttag[7]= 2406.63;  
       fbb[7]= 0.0790595;   fcc[7]= 0.130304;   fc[7]= 0.125299;   fll[7]= 0.665337;  
       fmcbb[7]= 0.0719225;   fmccc[7]= 0.118541;   fmcc[7]= 0.143211;   fmcll[7]= 0.666326;  
 }  
 else if ( idsys==56   || sysname=="metpileup_down" ) { 
      Fbb[1]= 1.11578;   Fcc[1]= 1.11578;   Fc[1]= 0.895551;   Fll[1]= 1.01109; 
         cafac[1]= 0.972584;  
         winpretag[1]= 1.10852e+06 ;   wintag[1]= 41037.2 ;   woutpretag[1]= 1.07813e+06;   wouttag[1]= 38311.5;  
       fbb[1]= 0.0133816;   fcc[1]= 0.0390129;   fc[1]= 0.124046;   fll[1]= 0.82356;  
       fmcbb[1]= 0.0119931;   fmccc[1]= 0.0349648;   fmcc[1]= 0.138513;   fmcll[1]= 0.814529;  
      Fbb[2]= 1.11214;   Fcc[2]= 1.11214;   Fc[2]= 0.892634;   Fll[2]= 1.00779; 
         cafac[2]= 0.882857;  
         winpretag[2]= 268951 ;   wintag[2]= 20989.1 ;   woutpretag[2]= 237445;   wouttag[2]= 18369.8;  
       fbb[2]= 0.0326023;   fcc[2]= 0.0786958;   fc[2]= 0.141355;   fll[2]= 0.747347;  
       fmcbb[2]= 0.0293149;   fmccc[2]= 0.0707605;   fmcc[2]= 0.158358;   fmcll[2]= 0.741567;  
      Fbb[3]= 1.10697;   Fcc[3]= 1.10697;   Fc[3]= 0.888486;   Fll[3]= 1.00311; 
         cafac[3]= 0.795189;  
         winpretag[3]= 60964.6 ;   wintag[3]= 7181.01 ;   woutpretag[3]= 48478.4;   wouttag[3]= 5761.61;  
       fbb[3]= 0.0544633;   fcc[3]= 0.106294;   fc[3]= 0.141023;   fll[3]= 0.698219;  
       fmcbb[3]= 0.0492002;   fmccc[3]= 0.0960224;   fmcc[3]= 0.158723;   fmcll[3]= 0.696054;  
      Fbb[4]= 1.1012;   Fcc[4]= 1.1012;   Fc[4]= 0.883852;   Fll[4]= 0.997879; 
         cafac[4]= 0.795039;  
         winpretag[4]= 13536.6 ;   wintag[4]= 2286.54 ;   woutpretag[4]= 10762.1;   wouttag[4]= 1850.5;  
       fbb[4]= 0.0754519;   fcc[4]= 0.125873;   fc[4]= 0.129975;   fll[4]= 0.6687;  
       fmcbb[4]= 0.0685179;   fmccc[4]= 0.114305;   fmcc[4]= 0.147055;   fmcll[4]= 0.670122;  
      Fbb[5]= 1.09492;   Fcc[5]= 1.09492;   Fc[5]= 0.878813;   Fll[5]= 0.99219; 
         cafac[5]= 0.727919;  
         winpretag[5]= 3745.28 ;   wintag[5]= 760.461 ;   woutpretag[5]= 2726.26;   wouttag[5]= 567.114;  
       fbb[5]= 0.0926579;   fcc[5]= 0.147231;   fc[5]= 0.113927;   fll[5]= 0.646184;  
       fmcbb[5]= 0.0846251;   fmccc[5]= 0.134467;   fmcc[5]= 0.129638;   fmcll[5]= 0.65127;  
      Fbb[6]= 1.10539;   Fcc[6]= 1.10539;   Fc[6]= 0.887214;   Fll[6]= 1.00167; 
         cafac[6]= 0.790726;  
         winpretag[6]= 78246.5 ;   wintag[6]= 10228 ;   woutpretag[6]= 61871.5;   wouttag[6]= 8190.47;  
       fbb[6]= 0.0599538;   fcc[6]= 0.111672;   fc[6]= 0.137795;   fll[6]= 0.690579;  
       fmcbb[6]= 0.0542377;   fmccc[6]= 0.101025;   fmcc[6]= 0.155312;   fmcll[6]= 0.689425;  
      Fbb[7]= 1.09983;   Fcc[7]= 1.09983;   Fc[7]= 0.882755;   Fll[7]= 0.99664; 
         cafac[7]= 0.778047;  
         winpretag[7]= 17281.9 ;   wintag[7]= 3047 ;   woutpretag[7]= 13446.1;   wouttag[7]= 2417.59;  
       fbb[7]= 0.0791975;   fcc[7]= 0.130522;   fc[7]= 0.126481;   fll[7]= 0.663799;  
       fmcbb[7]= 0.0720086;   fmccc[7]= 0.118675;   fmcc[7]= 0.14328;   fmcll[7]= 0.666036;  
 }  
 else if ( idsys==53   || sysname=="met_up" ) { 
      Fbb[1]= 1.13496;   Fcc[1]= 1.13496;   Fc[1]= 0.886405;   Fll[1]= 1.01155; 
         cafac[1]= 0.972097;  
         winpretag[1]= 1.11147e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08045e+06;   wouttag[1]= 38391;  
       fbb[1]= 0.0135633;   fcc[1]= 0.0396456;   fc[1]= 0.122816;   fll[1]= 0.823976;  
       fmcbb[1]= 0.0119504;   fmccc[1]= 0.0349311;   fmcc[1]= 0.138555;   fmcll[1]= 0.814564;  
      Fbb[2]= 1.13036;   Fcc[2]= 1.13036;   Fc[2]= 0.882811;   Fll[2]= 1.00745; 
         cafac[2]= 0.880465;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 236912;   wouttag[2]= 18364.3;  
       fbb[2]= 0.033128;   fcc[2]= 0.0799246;   fc[2]= 0.13985;   fll[2]= 0.747098;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.0707071;   fmcc[2]= 0.158414;   fmcll[2]= 0.741571;  
      Fbb[3]= 1.12418;   Fcc[3]= 1.12418;   Fc[3]= 0.877982;   Fll[3]= 1.00194; 
         cafac[3]= 0.795558;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48500.9;   wouttag[3]= 5791.89;  
       fbb[3]= 0.0553864;   fcc[3]= 0.107919;   fc[3]= 0.139524;   fll[3]= 0.697171;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.11739;   Fcc[4]= 1.11739;   Fc[4]= 0.87268;   Fll[4]= 0.995891; 
         cafac[4]= 0.78806;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10649.5;   wouttag[4]= 1840.22;  
       fbb[4]= 0.076704;   fcc[4]= 0.12791;   fc[4]= 0.128479;   fll[4]= 0.666907;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.1101;   Fcc[5]= 1.1101;   Fc[5]= 0.866984;   Fll[5]= 0.989392; 
         cafac[5]= 0.722318;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2706.76;   wouttag[5]= 562.757;  
       fbb[5]= 0.0942392;   fcc[5]= 0.149131;   fc[5]= 0.112293;   fll[5]= 0.644337;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.12232;   Fcc[6]= 1.12232;   Fc[6]= 0.876529;   Fll[6]= 1.00028; 
         cafac[6]= 0.789251;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61739.5;   wouttag[6]= 8209.27;  
       fbb[6]= 0.060967;   fcc[6]= 0.113383;   fc[6]= 0.136289;   fll[6]= 0.689361;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.1158;   Fcc[7]= 1.1158;   Fc[7]= 0.871437;   Fll[7]= 0.994473; 
         cafac[7]= 0.771378;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13314.6;   wouttag[7]= 2402.76;  
       fbb[7]= 0.0805304;   fcc[7]= 0.132541;   fc[7]= 0.124947;   fll[7]= 0.661982;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==54   || sysname=="met_down" ) { 
      Fbb[1]= 1.1362;   Fcc[1]= 1.1362;   Fc[1]= 0.889372;   Fll[1]= 1.01095; 
         cafac[1]= 0.974225;  
         winpretag[1]= 1.10546e+06 ;   wintag[1]= 40878 ;   woutpretag[1]= 1.07696e+06;   wouttag[1]= 38201.2;  
       fbb[1]= 0.0136401;   fcc[1]= 0.0397633;   fc[1]= 0.12318;   fll[1]= 0.823416;  
       fmcbb[1]= 0.012005;   fmccc[1]= 0.0349966;   fmcc[1]= 0.138503;   fmcll[1]= 0.814496;  
      Fbb[2]= 1.13141;   Fcc[2]= 1.13141;   Fc[2]= 0.885617;   Fll[2]= 1.00668; 
         cafac[2]= 0.88211;  
         winpretag[2]= 268550 ;   wintag[2]= 20963.4 ;   woutpretag[2]= 236891;   wouttag[2]= 18381.5;  
       fbb[2]= 0.0332107;   fcc[2]= 0.080067;   fc[2]= 0.140238;   fll[2]= 0.746485;  
       fmcbb[2]= 0.0293535;   fmccc[2]= 0.0707678;   fmcc[2]= 0.15835;   fmcll[2]= 0.741529;  
      Fbb[3]= 1.12515;   Fcc[3]= 1.12515;   Fc[3]= 0.880721;   Fll[3]= 1.00112; 
         cafac[3]= 0.792317;  
         winpretag[3]= 60881 ;   wintag[3]= 7167.56 ;   woutpretag[3]= 48237.1;   wouttag[3]= 5753.91;  
       fbb[3]= 0.0553894;   fcc[3]= 0.107981;   fc[3]= 0.139927;   fll[3]= 0.696702;  
       fmcbb[3]= 0.0492284;   fmccc[3]= 0.0959705;   fmcc[3]= 0.158878;   fmcll[3]= 0.695923;  
      Fbb[4]= 1.11834;   Fcc[4]= 1.11834;   Fc[4]= 0.875388;   Fll[4]= 0.995057; 
         cafac[4]= 0.792593;  
         winpretag[4]= 13534.1 ;   wintag[4]= 2282.62 ;   woutpretag[4]= 10727;   wouttag[4]= 1850.4;  
       fbb[4]= 0.0766467;   fcc[4]= 0.127828;   fc[4]= 0.128726;   fll[4]= 0.666799;  
       fmcbb[4]= 0.0685363;   fmccc[4]= 0.114301;   fmcc[4]= 0.14705;   fmcll[4]= 0.670112;  
      Fbb[5]= 1.11123;   Fcc[5]= 1.11123;   Fc[5]= 0.869825;   Fll[5]= 0.988733; 
         cafac[5]= 0.731103;  
         winpretag[5]= 3733 ;   wintag[5]= 760.61 ;   woutpretag[5]= 2729.21;   wouttag[5]= 572.615;  
       fbb[5]= 0.0940465;   fcc[5]= 0.149789;   fc[5]= 0.114209;   fll[5]= 0.641955;  
       fmcbb[5]= 0.0846328;   fmccc[5]= 0.134796;   fmcc[5]= 0.131301;   fmcll[5]= 0.64927;  
      Fbb[6]= 1.12329;   Fcc[6]= 1.12329;   Fc[6]= 0.879268;   Fll[6]= 0.999467; 
         cafac[6]= 0.788309;  
         winpretag[6]= 78148.2 ;   wintag[6]= 10210.8 ;   woutpretag[6]= 61604.9;   wouttag[6]= 8188.02;  
       fbb[6]= 0.0609538;   fcc[6]= 0.113452;   fc[6]= 0.136737;   fll[6]= 0.688857;  
       fmcbb[6]= 0.0542634;   fmccc[6]= 0.101;   fmcc[6]= 0.155512;   fmcll[6]= 0.689224;  
      Fbb[7]= 1.11679;   Fcc[7]= 1.11679;   Fc[7]= 0.87418;   Fll[7]= 0.993683; 
         cafac[7]= 0.777057;  
         winpretag[7]= 17267.1 ;   wintag[7]= 3043.23 ;   woutpretag[7]= 13417.5;   wouttag[7]= 2423.28;  
       fbb[7]= 0.0804272;   fcc[7]= 0.132599;   fc[7]= 0.125572;   fll[7]= 0.661401;  
       fmcbb[7]= 0.0720162;   fmccc[7]= 0.118732;   fmcc[7]= 0.143646;   fmcll[7]= 0.665606;  
 }  
 else if ( idsys==98   || sysname=="ifsr_nom" ) { 
      Fbb[1]= 1.04839;   Fcc[1]= 1.04839;   Fc[1]= 0.898615;   Fll[1]= 1.01446; 
         cafac[1]= 0.974366;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08296e+06;   wouttag[1]= 38301.5;  
       fbb[1]= 0.0125288;   fcc[1]= 0.0366218;   fc[1]= 0.124509;   fll[1]= 0.826341;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.04891;   Fcc[2]= 1.04891;   Fc[2]= 0.899063;   Fll[2]= 1.01497; 
         cafac[2]= 0.882995;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 237593;   wouttag[2]= 18120.2;  
       fbb[2]= 0.0307409;   fcc[2]= 0.0741655;   fc[2]= 0.142424;   fll[2]= 0.752669;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.04736;   Fcc[3]= 1.04736;   Fc[3]= 0.897736;   Fll[3]= 1.01347; 
         cafac[3]= 0.799767;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48757.5;   wouttag[3]= 5695.57;  
       fbb[3]= 0.0516018;   fcc[3]= 0.100545;   fc[3]= 0.142663;   fll[3]= 0.70519;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.04461;   Fcc[4]= 1.04461;   Fc[4]= 0.895375;   Fll[4]= 1.0108; 
         cafac[4]= 0.79647;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10763.1;   wouttag[4]= 1818.54;  
       fbb[4]= 0.071708;   fcc[4]= 0.119579;   fc[4]= 0.13182;   fll[4]= 0.676893;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.04121;   Fcc[5]= 1.04121;   Fc[5]= 0.892462;   Fll[5]= 1.00751; 
         cafac[5]= 0.735006;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2754.3;   wouttag[5]= 558.491;  
       fbb[5]= 0.0883913;   fcc[5]= 0.139877;   fc[5]= 0.115593;   fll[5]= 0.656139;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.04659;   Fcc[6]= 1.04659;   Fc[6]= 0.897073;   Fll[6]= 1.01272; 
         cafac[6]= 0.794922;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 62183.2;   wouttag[6]= 8082.49;  
       fbb[6]= 0.0568532;   fcc[6]= 0.105733;   fc[6]= 0.139483;   fll[6]= 0.697931;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.04387;   Fcc[7]= 1.04387;   Fc[7]= 0.894741;   Fll[7]= 1.01009; 
         cafac[7]= 0.780959;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13480;   wouttag[7]= 2376.79;  
       fbb[7]= 0.0753392;   fcc[7]= 0.123997;   fc[7]= 0.128288;   fll[7]= 0.672376;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==47   || sysname=="ifsr_up" ) { 
      Fbb[1]= 1.09961;   Fcc[1]= 1.09961;   Fc[1]= 0.857655;   Fll[1]= 1.01848; 
         cafac[1]= 0.966891;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.07466e+06;   wouttag[1]= 37409.7;  
       fbb[1]= 0.0131409;   fcc[1]= 0.0384109;   fc[1]= 0.118833;   fll[1]= 0.829615;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.09838;   Fcc[2]= 1.09838;   Fc[2]= 0.856698;   Fll[2]= 1.01734; 
         cafac[2]= 0.875087;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 235465;   wouttag[2]= 17920.2;  
       fbb[2]= 0.0321906;   fcc[2]= 0.0776632;   fc[2]= 0.135713;   fll[2]= 0.754433;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.09445;   Fcc[3]= 1.09445;   Fc[3]= 0.853637;   Fll[3]= 1.01371; 
         cafac[3]= 0.790821;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48212.1;   wouttag[3]= 5661.58;  
       fbb[3]= 0.0539219;   fcc[3]= 0.105065;   fc[3]= 0.135655;   fll[3]= 0.705358;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.08909;   Fcc[4]= 1.08909;   Fc[4]= 0.849451;   Fll[4]= 1.00874; 
         cafac[4]= 0.790135;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10677.5;   wouttag[4]= 1819.81;  
       fbb[4]= 0.0747611;   fcc[4]= 0.12467;   fc[4]= 0.125059;   fll[4]= 0.67551;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.08289;   Fcc[5]= 1.08289;   Fc[5]= 0.844618;   Fll[5]= 1.003; 
         cafac[5]= 0.722515;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2707.5;   wouttag[5]= 555.243;  
       fbb[5]= 0.0919296;   fcc[5]= 0.145476;   fc[5]= 0.109396;   fll[5]= 0.653198;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.09296;   Fcc[6]= 1.09296;   Fc[6]= 0.852475;   Fll[6]= 1.01233; 
         cafac[6]= 0.786162;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61497.9;   wouttag[6]= 8047.44;  
       fbb[6]= 0.0593723;   fcc[6]= 0.110418;   fc[6]= 0.132549;   fll[6]= 0.697661;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.08774;   Fcc[7]= 1.08774;   Fc[7]= 0.848397;   Fll[7]= 1.00749; 
         cafac[7]= 0.773011;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13342.8;   wouttag[7]= 2374.86;  
       fbb[7]= 0.078505;   fcc[7]= 0.129207;   fc[7]= 0.121643;   fll[7]= 0.670644;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==98   || sysname=="ifsr_nom" ) { 
      Fbb[1]= 1.11896;   Fcc[1]= 1.11896;   Fc[1]= 0.954858;   Fll[1]= 1.00083; 
         cafac[1]= 0.984897;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.09467e+06;   wouttag[1]= 40268.2;  
       fbb[1]= 0.0133722;   fcc[1]= 0.0390871;   fc[1]= 0.132301;   fll[1]= 0.815239;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.11299;   Fcc[2]= 1.11299;   Fc[2]= 0.949764;   Fll[2]= 0.995492; 
         cafac[2]= 0.89371;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 240476;   wouttag[2]= 19008.2;  
       fbb[2]= 0.032619;   fcc[2]= 0.0786965;   fc[2]= 0.150456;   fll[2]= 0.738229;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.10713;   Fcc[3]= 1.10713;   Fc[3]= 0.944762;   Fll[3]= 0.990249; 
         cafac[3]= 0.809165;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 49330.5;   wouttag[3]= 5956.74;  
       fbb[3]= 0.0545465;   fcc[3]= 0.106282;   fc[3]= 0.150137;   fll[3]= 0.689035;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.10167;   Fcc[4]= 1.10167;   Fc[4]= 0.940103;   Fll[4]= 0.985366; 
         cafac[4]= 0.791868;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10700.9;   wouttag[4]= 1858.46;  
       fbb[4]= 0.0756251;   fcc[4]= 0.126111;   fc[4]= 0.138405;   fll[4]= 0.659859;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.09619;   Fcc[5]= 1.09619;   Fc[5]= 0.935424;   Fll[5]= 0.980462; 
         cafac[5]= 0.734467;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2752.28;   wouttag[5]= 573.681;  
       fbb[5]= 0.0930586;   fcc[5]= 0.147262;   fc[5]= 0.121158;   fll[5]= 0.638521;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.10566;   Fcc[6]= 1.10566;   Fc[6]= 0.943503;   Fll[6]= 0.98893; 
         cafac[6]= 0.800943;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 62654.2;   wouttag[6]= 8405.92;  
       fbb[6]= 0.0600619;   fcc[6]= 0.1117;   fc[6]= 0.146702;   fll[6]= 0.681536;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.10048;   Fcc[7]= 1.10048;   Fc[7]= 0.939083;   Fll[7]= 0.984297; 
         cafac[7]= 0.777318;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13417.2;   wouttag[7]= 2431.85;  
       fbb[7]= 0.0794247;   fcc[7]= 0.130721;   fc[7]= 0.134646;   fll[7]= 0.655209;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==48   || sysname=="ifsr_down" ) { 
      Fbb[1]= 1.17061;   Fcc[1]= 1.17061;   Fc[1]= 0.91475;   Fll[1]= 1.00468; 
         cafac[1]= 0.977416;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08635e+06;   wouttag[1]= 39378.1;  
       fbb[1]= 0.0139894;   fcc[1]= 0.0408911;   fc[1]= 0.126744;   fll[1]= 0.818375;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.16244;   Fcc[2]= 1.16244;   Fc[2]= 0.908364;   Fll[2]= 0.997668; 
         cafac[2]= 0.885802;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 238348;   wouttag[2]= 18808.2;  
       fbb[2]= 0.034068;   fcc[2]= 0.0821924;   fc[2]= 0.143898;   fll[2]= 0.739842;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.15388;   Fcc[3]= 1.15388;   Fc[3]= 0.901681;   Fll[3]= 0.990328; 
         cafac[3]= 0.800211;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48784.6;   wouttag[3]= 5921.72;  
       fbb[3]= 0.0568499;   fcc[3]= 0.11077;   fc[3]= 0.14329;   fll[3]= 0.689089;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.14561;   Fcc[4]= 1.14561;   Fc[4]= 0.895211;   Fll[4]= 0.983222; 
         cafac[4]= 0.785671;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10617.2;   wouttag[4]= 1859.5;  
       fbb[4]= 0.0786408;   fcc[4]= 0.13114;   fc[4]= 0.131796;   fll[4]= 0.658424;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.13716;   Fcc[5]= 1.13716;   Fc[5]= 0.888616;   Fll[5]= 0.975978; 
         cafac[5]= 0.722007;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2705.6;   wouttag[5]= 570.097;  
       fbb[5]= 0.096537;   fcc[5]= 0.152767;   fc[5]= 0.115095;   fll[5]= 0.635601;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.15164;   Fcc[6]= 1.15164;   Fc[6]= 0.899924;   Fll[6]= 0.988398; 
         cafac[6]= 0.792208;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61970.9;   wouttag[6]= 8369.4;  
       fbb[6]= 0.0625595;   fcc[6]= 0.116345;   fc[6]= 0.139926;   fll[6]= 0.681169;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.14376;   Fcc[7]= 1.14376;   Fc[7]= 0.893771;   Fll[7]= 0.98164; 
         cafac[7]= 0.769481;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13281.9;   wouttag[7]= 2429.36;  
       fbb[7]= 0.0825486;   fcc[7]= 0.135862;   fc[7]= 0.128149;   fll[7]= 0.65344;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==72   || sysname=="iqopt3_one" ) { 
      Fbb[1]= 1.13134;   Fcc[1]= 1.13134;   Fc[1]= 0.888055;   Fll[1]= 1.0115; 
         cafac[1]= 0.978296;  
         winpretag[1]= 1.10548e+06 ;   wintag[1]= 40954.6 ;   woutpretag[1]= 1.08149e+06;   wouttag[1]= 38382.9;  
       fbb[1]= 0.0134852;   fcc[1]= 0.0394281;   fc[1]= 0.123032;   fll[1]= 0.824054;  
       fmcbb[1]= 0.0119197;   fmccc[1]= 0.0348507;   fmcc[1]= 0.138541;   fmcll[1]= 0.814688;  
      Fbb[2]= 1.12688;   Fcc[2]= 1.12688;   Fc[2]= 0.884555;   Fll[2]= 1.00751; 
         cafac[2]= 0.878731;  
         winpretag[2]= 269484 ;   wintag[2]= 21060.1 ;   woutpretag[2]= 236804;   wouttag[2]= 18364.2;  
       fbb[2]= 0.0330868;   fcc[2]= 0.0797779;   fc[2]= 0.140037;   fll[2]= 0.747098;  
       fmcbb[2]= 0.0293613;   fmccc[2]= 0.0707952;   fmcc[2]= 0.158314;   fmcll[2]= 0.74153;  
      Fbb[3]= 1.12078;   Fcc[3]= 1.12078;   Fc[3]= 0.879765;   Fll[3]= 1.00205; 
         cafac[3]= 0.788296;  
         winpretag[3]= 61379.2 ;   wintag[3]= 7248.94 ;   woutpretag[3]= 48385;   wouttag[3]= 5781.42;  
       fbb[3]= 0.0554512;   fcc[3]= 0.107956;   fc[3]= 0.139303;   fll[3]= 0.697289;  
       fmcbb[3]= 0.0494756;   fmccc[3]= 0.0963225;   fmcc[3]= 0.158341;   fmcll[3]= 0.695861;  
      Fbb[4]= 1.1141;   Fcc[4]= 1.1141;   Fc[4]= 0.874523;   Fll[4]= 0.996083; 
         cafac[4]= 0.786892;  
         winpretag[4]= 13523.6 ;   wintag[4]= 2298.36 ;   woutpretag[4]= 10641.6;   wouttag[4]= 1847.06;  
       fbb[4]= 0.0768694;   fcc[4]= 0.127833;   fc[4]= 0.127822;   fll[4]= 0.667475;  
       fmcbb[4]= 0.0689967;   fmccc[4]= 0.114741;   fmcc[4]= 0.146162;   fmcll[4]= 0.6701;  
      Fbb[5]= 1.10661;   Fcc[5]= 1.10661;   Fc[5]= 0.86864;   Fll[5]= 0.989382; 
         cafac[5]= 0.661976;  
         winpretag[5]= 4110.69 ;   wintag[5]= 825.336 ;   woutpretag[5]= 2721.18;   wouttag[5]= 561.905;  
       fbb[5]= 0.0946611;   fcc[5]= 0.151354;   fc[5]= 0.111102;   fll[5]= 0.642883;  
       fmcbb[5]= 0.0855417;   fmccc[5]= 0.136773;   fmcc[5]= 0.127903;   fmcll[5]= 0.649783;  
      Fbb[6]= 1.11889;   Fcc[6]= 1.11889;   Fc[6]= 0.878278;   Fll[6]= 1.00036; 
         cafac[6]= 0.779153;  
         winpretag[6]= 79013.5 ;   wintag[6]= 10372.6 ;   woutpretag[6]= 61563.6;   wouttag[6]= 8210.17;  
       fbb[6]= 0.0611953;   fcc[6]= 0.113656;   fc[6]= 0.135846;   fll[6]= 0.689303;  
       fmcbb[6]= 0.0546931;   fmccc[6]= 0.101579;   fmcc[6]= 0.154673;   fmcll[6]= 0.689054;  
      Fbb[7]= 1.11235;   Fcc[7]= 1.11235;   Fc[7]= 0.873144;   Fll[7]= 0.994513; 
         cafac[7]= 0.753285;  
         winpretag[7]= 17634.3 ;   wintag[7]= 3123.7 ;   woutpretag[7]= 13283.7;   wouttag[7]= 2408.12;  
       fbb[7]= 0.0810383;   fcc[7]= 0.133345;   fc[7]= 0.123904;   fll[7]= 0.661713;  
       fmcbb[7]= 0.0728535;   fmccc[7]= 0.119877;   fmcc[7]= 0.141906;   fmcll[7]= 0.665364;  
 }  
 else if ( idsys==71   || sysname=="ptjmin_one" ) { 
      Fbb[1]= 1.12811;   Fcc[1]= 1.12811;   Fc[1]= 0.887965;   Fll[1]= 1.01168; 
         cafac[1]= 0.972485;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08087e+06;   wouttag[1]= 38404.6;  
       fbb[1]= 0.0134815;   fcc[1]= 0.0394065;   fc[1]= 0.123033;   fll[1]= 0.824079;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.12385;   Fcc[2]= 1.12385;   Fc[2]= 0.884611;   Fll[2]= 1.00786; 
         cafac[2]= 0.879066;  
         winpretag[2]= 269187 ;   wintag[2]= 21068.6 ;   woutpretag[2]= 236633;   wouttag[2]= 18364.2;  
       fbb[2]= 0.033094;   fcc[2]= 0.0797515;   fc[2]= 0.14002;   fll[2]= 0.747134;  
       fmcbb[2]= 0.029447;   fmccc[2]= 0.0709628;   fmcc[2]= 0.158285;   fmcll[2]= 0.741306;  
      Fbb[3]= 1.11805;   Fcc[3]= 1.11805;   Fc[3]= 0.880047;   Fll[3]= 1.00266; 
         cafac[3]= 0.79625;  
         winpretag[3]= 60869.8 ;   wintag[3]= 7186.34 ;   woutpretag[3]= 48467.5;   wouttag[3]= 5783.62;  
       fbb[3]= 0.0552161;   fcc[3]= 0.107491;   fc[3]= 0.13963;   fll[3]= 0.697663;  
       fmcbb[3]= 0.0493861;   fmccc[3]= 0.0961412;   fmcc[3]= 0.158662;   fmcll[3]= 0.695811;  
      Fbb[4]= 1.11152;   Fcc[4]= 1.11152;   Fc[4]= 0.874908;   Fll[4]= 0.996808; 
         cafac[4]= 0.788595;  
         winpretag[4]= 13493.5 ;   wintag[4]= 2289 ;   woutpretag[4]= 10640.9;   wouttag[4]= 1841.7;  
       fbb[4]= 0.0765592;   fcc[4]= 0.127348;   fc[4]= 0.128131;   fll[4]= 0.667961;  
       fmcbb[4]= 0.0688777;   fmccc[4]= 0.114571;   fmcc[4]= 0.146451;   fmcll[4]= 0.6701;  
      Fbb[5]= 1.10469;   Fcc[5]= 1.10469;   Fc[5]= 0.869533;   Fll[5]= 0.990684; 
         cafac[5]= 0.733758;  
         winpretag[5]= 3690.46 ;   wintag[5]= 747.061 ;   woutpretag[5]= 2707.91;   wouttag[5]= 563.008;  
       fbb[5]= 0.0941475;   fcc[5]= 0.148478;   fc[5]= 0.112853;   fll[5]= 0.644522;  
       fmcbb[5]= 0.085225;   fmccc[5]= 0.134406;   fmcc[5]= 0.129786;   fmcll[5]= 0.650583;  
      Fbb[6]= 1.11628;   Fcc[6]= 1.11628;   Fc[6]= 0.878652;   Fll[6]= 1.00107; 
         cafac[6]= 0.790661;  
         winpretag[6]= 78053.7 ;   wintag[6]= 10222.4 ;   woutpretag[6]= 61714;   wouttag[6]= 8201.38;  
       fbb[6]= 0.0607816;   fcc[6]= 0.112896;   fc[6]= 0.136354;   fll[6]= 0.689968;  
       fmcbb[6]= 0.0544502;   fmccc[6]= 0.101136;   fmcc[6]= 0.155186;   fmcll[6]= 0.689228;  
      Fbb[7]= 1.11005;   Fcc[7]= 1.11005;   Fc[7]= 0.873748;   Fll[7]= 0.995486; 
         cafac[7]= 0.774816;  
         winpretag[7]= 17183.9 ;   wintag[7]= 3036.06 ;   woutpretag[7]= 13314.4;   wouttag[7]= 2404.51;  
       fbb[7]= 0.0803548;   fcc[7]= 0.131908;   fc[7]= 0.124834;   fll[7]= 0.662903;  
       fmcbb[7]= 0.0723885;   fmccc[7]= 0.118831;   fmcc[7]= 0.142872;   fmcll[7]= 0.665909;  
 }  
 else if ( idsys==69   || sysname=="powhe_one" ) { 
      Fbb[1]= 1.11657;   Fcc[1]= 1.11657;   Fc[1]= 0.902958;   Fll[1]= 1.0098; 
         cafac[1]= 0.975282;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.08398e+06;   wouttag[1]= 38772.1;  
       fbb[1]= 0.0133436;   fcc[1]= 0.0390034;   fc[1]= 0.12511;   fll[1]= 0.822543;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.11262;   Fcc[2]= 1.11262;   Fc[2]= 0.899762;   Fll[2]= 1.00622; 
         cafac[2]= 0.883939;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 237847;   wouttag[2]= 18462.2;  
       fbb[2]= 0.0326079;   fcc[2]= 0.0786698;   fc[2]= 0.142535;   fll[2]= 0.746187;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.10734;   Fcc[3]= 1.10734;   Fc[3]= 0.895499;   Fll[3]= 1.00146; 
         cafac[3]= 0.79873;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 48694.3;   wouttag[3]= 5806.87;  
       fbb[3]= 0.054557;   fcc[3]= 0.106303;   fc[3]= 0.142308;   fll[3]= 0.696833;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.10156;   Fcc[4]= 1.10156;   Fc[4]= 0.890825;   Fll[4]= 0.996229; 
         cafac[4]= 0.790786;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 10686.3;   wouttag[4]= 1841.51;  
       fbb[4]= 0.0756176;   fcc[4]= 0.126098;   fc[4]= 0.13115;   fll[4]= 0.667134;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.09535;   Fcc[5]= 1.09535;   Fc[5]= 0.885802;   Fll[5]= 0.990612; 
         cafac[5]= 0.72497;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2716.7;   wouttag[5]= 562.735;  
       fbb[5]= 0.0929876;   fcc[5]= 0.14715;   fc[5]= 0.114731;   fll[5]= 0.645132;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.10576;   Fcc[6]= 1.10576;   Fc[6]= 0.894219;   Fll[6]= 1.00003; 
         cafac[6]= 0.792335;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 61980.8;   wouttag[6]= 8225.16;  
       fbb[6]= 0.0600675;   fcc[6]= 0.111711;   fc[6]= 0.139039;   fll[6]= 0.689182;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.10021;   Fcc[7]= 1.10021;   Fc[7]= 0.889729;   Fll[7]= 0.995004; 
         cafac[7]= 0.77411;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13361.8;   wouttag[7]= 2404;  
       fbb[7]= 0.0794054;   fcc[7]= 0.130689;   fc[7]= 0.127569;   fll[7]= 0.662336;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==70   || sysname=="powpy_one" ) { 
      Fbb[1]= 1.03824;   Fcc[1]= 1.03824;   Fc[1]= 0.942134;   Fll[1]= 1.00764; 
         cafac[1]= 0.981764;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41244.5 ;   woutpretag[1]= 1.09119e+06;   wouttag[1]= 39461.7;  
       fbb[1]= 0.0124076;   fcc[1]= 0.0362674;   fc[1]= 0.130538;   fll[1]= 0.820787;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.03791;   Fcc[2]= 1.03791;   Fc[2]= 0.941828;   Fll[2]= 1.00731; 
         cafac[2]= 0.891584;  
         winpretag[2]= 269076 ;   wintag[2]= 21013.2 ;   woutpretag[2]= 239904;   wouttag[2]= 18531.2;  
       fbb[2]= 0.0304183;   fcc[2]= 0.0733872;   fc[2]= 0.149199;   fll[2]= 0.746996;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.0365;   Fcc[3]= 1.0365;   Fc[3]= 0.940556;   Fll[3]= 1.00595; 
         cafac[3]= 0.817017;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7191.9 ;   woutpretag[3]= 49809.1;   wouttag[3]= 5861.16;  
       fbb[3]= 0.0510668;   fcc[3]= 0.0995021;   fc[3]= 0.149468;   fll[3]= 0.699963;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.03452;   Fcc[4]= 1.03452;   Fc[4]= 0.938753;   Fll[4]= 1.00403; 
         cafac[4]= 0.82272;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2284.83 ;   woutpretag[4]= 11117.8;   wouttag[4]= 1884.73;  
       fbb[4]= 0.0710151;   fcc[4]= 0.118423;   fc[4]= 0.138206;   fll[4]= 0.672355;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.03219;   Fcc[5]= 1.03219;   Fc[5]= 0.936639;   Fll[5]= 1.00177; 
         cafac[5]= 0.771692;  
         winpretag[5]= 3747.32 ;   wintag[5]= 756.891 ;   woutpretag[5]= 2891.78;   wouttag[5]= 587.345;  
       fbb[5]= 0.0876253;   fcc[5]= 0.138664;   fc[5]= 0.121315;   fll[5]= 0.652395;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.03595;   Fcc[6]= 1.03595;   Fc[6]= 0.940056;   Fll[6]= 1.00542; 
         cafac[6]= 0.815138;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10233.6 ;   woutpretag[6]= 63764.6;   wouttag[6]= 8337.21;  
       fbb[6]= 0.0562754;   fcc[6]= 0.104658;   fc[6]= 0.146166;   fll[6]= 0.6929;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.03401;   Fcc[7]= 1.03401;   Fc[7]= 0.938293;   Fll[7]= 1.00353; 
         cafac[7]= 0.809851;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3041.72 ;   woutpretag[7]= 13978.7;   wouttag[7]= 2471.82;  
       fbb[7]= 0.0746275;   fcc[7]= 0.122825;   fc[7]= 0.134533;   fll[7]= 0.668014;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==15   || sysname=="btagb_up" ) { 
      Fbb[1]= 1.05097;   Fcc[1]= 1.05097;   Fc[1]= 0.848546;   Fll[1]= 1.02283; 
         cafac[1]= 0.965773;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41486.3 ;   woutpretag[1]= 1.07341e+06;   wouttag[1]= 37177.5;  
       fbb[1]= 0.0125596;   fcc[1]= 0.0367119;   fc[1]= 0.117571;   fll[1]= 0.833157;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.05304;   Fcc[2]= 1.05304;   Fc[2]= 0.850217;   Fll[2]= 1.02484; 
         cafac[2]= 0.874631;  
         winpretag[2]= 269076 ;   wintag[2]= 21237 ;   woutpretag[2]= 235343;   wouttag[2]= 17849;  
       fbb[2]= 0.0308618;   fcc[2]= 0.0744573;   fc[2]= 0.134686;   fll[2]= 0.759995;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.05179;   Fcc[3]= 1.05179;   Fc[3]= 0.849208;   Fll[3]= 1.02363; 
         cafac[3]= 0.789682;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7293.26 ;   woutpretag[3]= 48142.7;   wouttag[3]= 5645.11;  
       fbb[3]= 0.0518198;   fcc[3]= 0.100969;   fc[3]= 0.134952;   fll[3]= 0.712259;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.04853;   Fcc[4]= 1.04853;   Fc[4]= 0.846579;   Fll[4]= 1.02046; 
         cafac[4]= 0.782432;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2316.12 ;   woutpretag[4]= 10573.4;   wouttag[4]= 1800.62;  
       fbb[4]= 0.0719773;   fcc[4]= 0.120028;   fc[4]= 0.124636;   fll[4]= 0.683359;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.04426;   Fcc[5]= 1.04426;   Fc[5]= 0.843129;   Fll[5]= 1.0163; 
         cafac[5]= 0.718082;  
         winpretag[5]= 3747.32 ;   wintag[5]= 770.705 ;   woutpretag[5]= 2690.89;   wouttag[5]= 553.362;  
       fbb[5]= 0.0886501;   fcc[5]= 0.140286;   fc[5]= 0.109204;   fll[5]= 0.66186;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.05086;   Fcc[6]= 1.05086;   Fc[6]= 0.84846;   Fll[6]= 1.02272; 
         cafac[6]= 0.78364;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10380.1 ;   woutpretag[6]= 61300.6;   wouttag[6]= 8012.09;  
       fbb[6]= 0.0570853;   fcc[6]= 0.106164;   fc[6]= 0.131924;   fll[6]= 0.704826;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.0476;   Fcc[7]= 1.0476;   Fc[7]= 0.845828;   Fll[7]= 1.01955; 
         cafac[7]= 0.76618;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3086.83 ;   woutpretag[7]= 13224.9;   wouttag[7]= 2353.97;  
       fbb[7]= 0.0756085;   fcc[7]= 0.12444;   fc[7]= 0.121275;   fll[7]= 0.678677;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==16   || sysname=="btagb_down" ) { 
      Fbb[1]= 1.22867;   Fcc[1]= 1.22867;   Fc[1]= 0.925554;   Fll[1]= 0.999502; 
         cafac[1]= 0.978826;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 41002.6 ;   woutpretag[1]= 1.08792e+06;   wouttag[1]= 39680.8;  
       fbb[1]= 0.0146832;   fcc[1]= 0.0429191;   fc[1]= 0.128241;   fll[1]= 0.814157;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.21565;   Fcc[2]= 1.21565;   Fc[2]= 0.915749;   Fll[2]= 0.988913; 
         cafac[2]= 0.886401;  
         winpretag[2]= 269076 ;   wintag[2]= 20787.7 ;   woutpretag[2]= 238509;   wouttag[2]= 18898.5;  
       fbb[2]= 0.0356275;   fcc[2]= 0.0859549;   fc[2]= 0.145067;   fll[2]= 0.73335;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.20335;   Fcc[3]= 1.20335;   Fc[3]= 0.906481;   Fll[3]= 0.978906; 
         cafac[3]= 0.801499;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7088.6 ;   woutpretag[3]= 48863.1;   wouttag[3]= 5940.64;  
       fbb[3]= 0.0592868;   fcc[3]= 0.115519;   fc[3]= 0.144053;   fll[3]= 0.681142;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.19221;   Fcc[4]= 1.19221;   Fc[4]= 0.898091;   Fll[4]= 0.969845; 
         cafac[4]= 0.793582;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2252.17 ;   woutpretag[4]= 10724.1;   wouttag[4]= 1878.94;  
       fbb[4]= 0.0818399;   fcc[4]= 0.136475;   fc[4]= 0.13222;   fll[4]= 0.649466;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.18122;   Fcc[5]= 1.18122;   Fc[5]= 0.889815;   Fll[5]= 0.960908; 
         cafac[5]= 0.726508;  
         winpretag[5]= 3747.32 ;   wintag[5]= 742.287 ;   woutpretag[5]= 2722.46;   wouttag[5]= 571.342;  
       fbb[5]= 0.100277;   fcc[5]= 0.158686;   fc[5]= 0.11525;   fll[5]= 0.625787;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.20033;   Fcc[6]= 1.20033;   Fc[6]= 0.904211;   Fll[6]= 0.976454; 
         cafac[6]= 0.794873;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10083.1 ;   woutpretag[6]= 62179.3;   wouttag[6]= 8406.81;  
       fbb[6]= 0.0652048;   fcc[6]= 0.121265;   fc[6]= 0.140593;   fll[6]= 0.672938;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.18981;   Fcc[7]= 1.18981;   Fc[7]= 0.896281;   Fll[7]= 0.967891; 
         cafac[7]= 0.776473;  
         winpretag[7]= 17260.9 ;   wintag[7]= 2994.46 ;   woutpretag[7]= 13402.6;   wouttag[7]= 2449.8;  
       fbb[7]= 0.0858717;   fcc[7]= 0.141332;   fc[7]= 0.128509;   fll[7]= 0.644288;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==85   || sysname=="btagc_up" ) { 
      Fbb[1]= 1.06043;   Fcc[1]= 1.06043;   Fc[1]= 0.772034;   Fll[1]= 1.0353; 
         cafac[1]= 0.952188;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 45514.2 ;   woutpretag[1]= 1.05831e+06;   wouttag[1]= 38182.7;  
       fbb[1]= 0.0126727;   fcc[1]= 0.0370425;   fc[1]= 0.10697;   fll[1]= 0.843315;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.06458;   Fcc[2]= 1.06458;   Fc[2]= 0.77505;   Fll[2]= 1.03934; 
         cafac[2]= 0.860209;  
         winpretag[2]= 269076 ;   wintag[2]= 22811.5 ;   woutpretag[2]= 231462;   wouttag[2]= 18234;  
       fbb[2]= 0.0312;   fcc[2]= 0.0752731;   fc[2]= 0.122779;   fll[2]= 0.770748;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.0635;   Fcc[3]= 1.0635;   Fc[3]= 0.774269;   Fll[3]= 1.0383; 
         cafac[3]= 0.776144;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7728.65 ;   woutpretag[3]= 47317.3;   wouttag[3]= 5747.88;  
       fbb[3]= 0.0523969;   fcc[3]= 0.102094;   fc[3]= 0.123043;   fll[3]= 0.722467;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.05922;   Fcc[4]= 1.05922;   Fc[4]= 0.771152;   Fll[4]= 1.03412; 
         cafac[4]= 0.770224;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2420.25 ;   woutpretag[4]= 10408.5;   wouttag[4]= 1826.71;  
       fbb[4]= 0.072711;   fcc[4]= 0.121251;   fc[4]= 0.113531;   fll[4]= 0.692506;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.05336;   Fcc[5]= 1.05336;   Fc[5]= 0.766887;   Fll[5]= 1.0284; 
         cafac[5]= 0.707132;  
         winpretag[5]= 3747.32 ;   wintag[5]= 801.103 ;   woutpretag[5]= 2649.85;   wouttag[5]= 560.453;  
       fbb[5]= 0.089423;   fcc[5]= 0.141509;   fc[5]= 0.0993286;   fll[5]= 0.669739;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.06227;   Fcc[6]= 1.06227;   Fc[6]= 0.773372;   Fll[6]= 1.03709; 
         cafac[6]= 0.770488;  
         winpretag[6]= 78225.5 ;   wintag[6]= 10950 ;   woutpretag[6]= 60271.9;   wouttag[6]= 8147.83;  
       fbb[6]= 0.057705;   fcc[6]= 0.107317;   fc[6]= 0.120249;   fll[6]= 0.714729;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.05795;   Fcc[7]= 1.05795;   Fc[7]= 0.770223;   Fll[7]= 1.03287; 
         cafac[7]= 0.754256;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3221.36 ;   woutpretag[7]= 13019.1;   wouttag[7]= 2387.09;  
       fbb[7]= 0.0763549;   fcc[7]= 0.125669;   fc[7]= 0.110435;   fll[7]= 0.687542;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==86   || sysname=="btagc_down" ) { 
      Fbb[1]= 1.21331;   Fcc[1]= 1.21331;   Fc[1]= 1.0343;   Fll[1]= 0.981888; 
         cafac[1]= 0.999375;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 36974.7 ;   woutpretag[1]= 1.11076e+06;   wouttag[1]= 38575.2;  
       fbb[1]= 0.0144997;   fcc[1]= 0.0423827;   fc[1]= 0.143309;   fll[1]= 0.799809;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.19734;   Fcc[2]= 1.19734;   Fc[2]= 1.02069;   Fll[2]= 0.968965; 
         cafac[2]= 0.907955;  
         winpretag[2]= 269076 ;   wintag[2]= 19209 ;   woutpretag[2]= 244309;   wouttag[2]= 18494.8;  
       fbb[2]= 0.035091;   fcc[2]= 0.0846605;   fc[2]= 0.161691;   fll[2]= 0.718557;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.18506;   Fcc[3]= 1.18506;   Fc[3]= 1.01022;   Fll[3]= 0.959029; 
         cafac[3]= 0.821621;  
         winpretag[3]= 60964.7 ;   wintag[3]= 6652.25 ;   woutpretag[3]= 50089.8;   wouttag[3]= 5836.68;  
       fbb[3]= 0.058386;   fcc[3]= 0.113763;   fc[3]= 0.160539;   fll[3]= 0.667311;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.17571;   Fcc[4]= 1.17571;   Fc[4]= 1.00225;   Fll[4]= 0.951458; 
         cafac[4]= 0.811567;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2148.1 ;   woutpretag[4]= 10967.1;   wouttag[4]= 1853.23;  
       fbb[4]= 0.0807072;   fcc[4]= 0.134586;   fc[4]= 0.147554;   fll[4]= 0.637153;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.1673;   Fcc[5]= 1.1673;   Fc[5]= 0.995084;   Fll[5]= 0.944657; 
         cafac[5]= 0.742555;  
         winpretag[5]= 3747.32 ;   wintag[5]= 711.837 ;   woutpretag[5]= 2782.6;   wouttag[5]= 564.42;  
       fbb[5]= 0.0990956;   fcc[5]= 0.156816;   fc[5]= 0.128885;   fll[5]= 0.615204;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.18258;   Fcc[6]= 1.18258;   Fc[6]= 1.0081;   Fll[6]= 0.957016; 
         cafac[6]= 0.814367;  
         winpretag[6]= 78225.5 ;   wintag[6]= 9512.19 ;   woutpretag[6]= 63704.2;   wouttag[6]= 8270.6;  
       fbb[6]= 0.0642402;   fcc[6]= 0.119471;   fc[6]= 0.156747;   fll[6]= 0.659542;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.17387;   Fcc[7]= 1.17387;   Fc[7]= 1.00068;   Fll[7]= 0.949974; 
         cafac[7]= 0.794015;  
         winpretag[7]= 17260.9 ;   wintag[7]= 2859.94 ;   woutpretag[7]= 13705.4;   wouttag[7]= 2417.26;  
       fbb[7]= 0.0847218;   fcc[7]= 0.139439;   fc[7]= 0.143478;   fll[7]= 0.632361;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==17   || sysname=="bmtag_up" ) { 
      Fbb[1]= 0.876056;   Fcc[1]= 0.876056;   Fc[1]= 0.852045;   Fll[1]= 1.0323; 
         cafac[1]= 0.967328;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 43910.6 ;   woutpretag[1]= 1.07514e+06;   wouttag[1]= 38818.9;  
       fbb[1]= 0.0104693;   fcc[1]= 0.0306019;   fc[1]= 0.118056;   fll[1]= 0.840873;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 0.88659;   Fcc[2]= 0.88659;   Fc[2]= 0.86229;   Fll[2]= 1.04471; 
         cafac[2]= 0.878117;  
         winpretag[2]= 269076 ;   wintag[2]= 22455.7 ;   woutpretag[2]= 236280;   wouttag[2]= 18260.7;  
       fbb[2]= 0.0259836;   fcc[2]= 0.0626881;   fc[2]= 0.136599;   fll[2]= 0.774729;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 0.893061;   Fcc[3]= 0.893061;   Fc[3]= 0.868585;   Fll[3]= 1.05234; 
         cafac[3]= 0.792169;  
         winpretag[3]= 60964.7 ;   wintag[3]= 7751.76 ;   woutpretag[3]= 48294.3;   wouttag[3]= 5741.32;  
       fbb[3]= 0.0439996;   fcc[3]= 0.0857318;   fc[3]= 0.138031;   fll[3]= 0.732238;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 0.89654;   Fcc[4]= 0.89654;   Fc[4]= 0.871968;   Fll[4]= 1.05644; 
         cafac[4]= 0.784129;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2504.59 ;   woutpretag[4]= 10596.4;   wouttag[4]= 1864.97;  
       fbb[4]= 0.0615436;   fcc[4]= 0.102629;   fc[4]= 0.128374;   fll[4]= 0.707454;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 0.898796;   Fcc[5]= 0.898796;   Fc[5]= 0.874162;   Fll[5]= 1.0591; 
         cafac[5]= 0.72111;  
         winpretag[5]= 3747.32 ;   wintag[5]= 846.339 ;   woutpretag[5]= 2702.23;   wouttag[5]= 581.498;  
       fbb[5]= 0.0763012;   fcc[5]= 0.120744;   fc[5]= 0.113223;   fll[5]= 0.689731;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 0.893934;   Fcc[6]= 0.893934;   Fc[6]= 0.869433;   Fll[6]= 1.05337; 
         cafac[6]= 0.786273;  
         winpretag[6]= 78225.5 ;   wintag[6]= 11102.7 ;   woutpretag[6]= 61506.6;   wouttag[6]= 8199.41;  
       fbb[6]= 0.0485605;   fcc[6]= 0.0903105;   fc[6]= 0.135185;   fll[6]= 0.725944;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 0.897029;   Fcc[7]= 0.897029;   Fc[7]= 0.872444;   Fll[7]= 1.05701; 
         cafac[7]= 0.768405;  
         winpretag[7]= 17260.9 ;   wintag[7]= 3350.93 ;   woutpretag[7]= 13263.3;   wouttag[7]= 2446.99;  
       fbb[7]= 0.0647412;   fcc[7]= 0.106554;   fc[7]= 0.125091;   fll[7]= 0.703614;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
 else if ( idsys==18   || sysname=="bmtag_down" ) { 
      Fbb[1]= 1.39858;   Fcc[1]= 1.39858;   Fc[1]= 0.919144;   Fll[1]= 0.990813; 
         cafac[1]= 0.976727;  
         winpretag[1]= 1.11146e+06 ;   wintag[1]= 38578.3 ;   woutpretag[1]= 1.08559e+06;   wouttag[1]= 37982;  
       fbb[1]= 0.0167137;   fcc[1]= 0.0488544;   fc[1]= 0.127353;   fll[1]= 0.807079;  
       fmcbb[1]= 0.0119505;   fmccc[1]= 0.0349315;   fmcc[1]= 0.138556;   fmcll[1]= 0.814562;  
      Fbb[2]= 1.37083;   Fcc[2]= 1.37083;   Fc[2]= 0.900908;   Fll[2]= 0.971155; 
         cafac[2]= 0.88238;  
         winpretag[2]= 269076 ;   wintag[2]= 19579 ;   woutpretag[2]= 237428;   wouttag[2]= 18467.8;  
       fbb[2]= 0.0401754;   fcc[2]= 0.0969272;   fc[2]= 0.142716;   fll[2]= 0.720181;  
       fmcbb[2]= 0.0293074;   fmccc[2]= 0.070707;   fmcc[2]= 0.158414;   fmcll[2]= 0.741572;  
      Fbb[3]= 1.34652;   Fcc[3]= 1.34652;   Fc[3]= 0.884934;   Fll[3]= 0.953936; 
         cafac[3]= 0.798432;  
         winpretag[3]= 60964.7 ;   wintag[3]= 6648.62 ;   woutpretag[3]= 48676.1;   wouttag[3]= 5851;  
       fbb[3]= 0.0663408;   fcc[3]= 0.129263;   fc[3]= 0.140629;   fll[3]= 0.663767;  
       fmcbb[3]= 0.0492683;   fmccc[3]= 0.0959977;   fmcc[3]= 0.158915;   fmcll[3]= 0.695819;  
      Fbb[4]= 1.32575;   Fcc[4]= 1.32575;   Fc[4]= 0.871284;   Fll[4]= 0.939221; 
         cafac[4]= 0.791286;  
         winpretag[4]= 13513.5 ;   wintag[4]= 2072.29 ;   woutpretag[4]= 10693.1;   wouttag[4]= 1821.56;  
       fbb[4]= 0.0910071;   fcc[4]= 0.151762;   fc[4]= 0.128273;   fll[4]= 0.628958;  
       fmcbb[4]= 0.0686456;   fmccc[4]= 0.114472;   fmcc[4]= 0.147223;   fmcll[4]= 0.669659;  
      Fbb[5]= 1.30595;   Fcc[5]= 1.30595;   Fc[5]= 0.858271;   Fll[5]= 0.925193; 
         cafac[5]= 0.723134;  
         winpretag[5]= 3747.32 ;   wintag[5]= 676.942 ;   woutpretag[5]= 2709.82;   wouttag[5]= 552.169;  
       fbb[5]= 0.110866;   fcc[5]= 0.175442;   fc[5]= 0.111165;   fll[5]= 0.602528;  
       fmcbb[5]= 0.0848927;   fmccc[5]= 0.13434;   fmcc[5]= 0.129522;   fmcll[5]= 0.651245;  
      Fbb[6]= 1.3409;   Fcc[6]= 1.3409;   Fc[6]= 0.881238;   Fll[6]= 0.949951; 
         cafac[6]= 0.791696;  
         winpretag[6]= 78225.5 ;   wintag[6]= 9397.86 ;   woutpretag[6]= 61930.8;   wouttag[6]= 8241.43;  
       fbb[6]= 0.0728407;   fcc[6]= 0.135465;   fc[6]= 0.137021;   fll[6]= 0.654673;  
       fmcbb[6]= 0.0543223;   fmccc[6]= 0.101026;   fmcc[6]= 0.155487;   fmcll[6]= 0.689165;  
      Fbb[7]= 1.3214;   Fcc[7]= 1.3214;   Fc[7]= 0.868425;   Fll[7]= 0.93614; 
         cafac[7]= 0.773726;  
         winpretag[7]= 17260.9 ;   wintag[7]= 2749.24 ;   woutpretag[7]= 13355.2;   wouttag[7]= 2373.09;  
       fbb[7]= 0.0953694;   fcc[7]= 0.156963;   fc[7]= 0.124515;   fll[7]= 0.623152;  
       fmcbb[7]= 0.0721729;   fmccc[7]= 0.118785;   fmcc[7]= 0.14338;   fmcll[7]= 0.665662;  
 }  
       else { cout << " WARNING: Ffactors request for unknown variation. Set to Nominal! " << idsys << "  " << sysname <<  endl; }

     return;


}

// This is the function to use:
void SetWflavors_muon(int idsys, TString sysname, int ijet, double Wjets_in[4],double Wjets_out[4], double& canorm, int _mode) {
	//-
	// This function corrects the Wjet ligth and heavy flavor components and return the normalization based from Charge Assymmetry.
	//
	// NOTE: this function shoud always be called with PRETAGGED event counts. NEVER WITH TAGGED counts.
	//
	// INPUTS:
	// for idsys there are two input modes: the individul systematics (long list), or the combined systematics (4 checks).
	// Long list:
	// input idsys>0..999 or sysname are passed to the function GetFFactors_muon (see below) that returns the relevant HF factors.
	// This allows to use the full list of systematics.
	// Short list:
	// most analysers will use the total sum of all systematics. In that case only four checks up and two down are needed.
	// Complete set: idsys=2000 (up), -2000 (down) ---> checks Fbb + Fcc + Fc versus Fl
	//                     2001 (up), -2001 (down) ---> checks Fbb+Fcc versus Fc
	//                     2002 (up), -2002 (down) ---> checks Fll  OBSOLETE, not needed anymore, return nominal
	//                     2003 (up), -2003 (down) ---> check Fbb+25% in 1,3,4,5 jetbin. --> see remark *
	//                     2004 (up), -2004 (down) ---> check Fc+25% in 1,3,4,5 jetbin.  --> see remark *
	//                     2005 (up), -2005 (down) ---> check canorm. Please note that canorm is returned and NOT yet applied to Wjets_out.
	// * please note that the 25% per jetbin change should be applied UNCORRELATED to all jetbins. So, for example, when you use simultaneoulsy 
	//   jetbin 4 and 5 in your analysis, then first use the changed 4jetbin HF factors and 5 nominal and then change 5jetbin and keep 4jet nominal.
	//  
	//
	// input ijet = jetbin you want to set. Available ijet values: 1=1ex, 2=2ex, 3=3ex, 4=4ex, -5=5=5in, -3=3in, -4=4in.
	// input Wjets_in are the PRETAGGED event numbers (for ijet) for 0:Wbb, 1:Wcc, 2:Wc, 3:Wlight
	// The flavor fractions defined as bb, cc, c and light flavor jets correspond to the HFOR flag equals to 0, 1, 2 and 3 respectively
	//--
	// OUTPUTS:
	// output Wjets_out are the PRETAGGED event counts (for ijet) for 0:Wbb, 1:Wcc, 2:Wc, 3:Wlight 
	// output canorm, which is  NOT yet applied to Wjets_out. This normlization factor can (or should be) be used for pretagged and tagged events.
	//
	// NOTE: if you want to use the HF corrections for your tagged sample OR after any other cut, then you still have to input
	//       the PRETAGGED number of events and use the pretagged fraction: Wjets_out[iflavor]/Wjets_in[iflavor] to reweigh your events.
	//
	//
	// TYPICAL USE
	// LONG:
	// Typical use for the full systematics list:
	//      .....
	//      int idsys=0; // for nominal, but you have to loop over many idsys, see in explanation on www.nikhef.nl/~h73/factors.html
	//      int ijet=3; // in case you input the 3jetbin
	//      double Wjets_in[4]={ 6888, 14281, 22937, 102347} ; 
	//      double Wjets_out[4]; 
	//      SetWflavors(idsys, "" , ijet, Wjets_in, Wjets_out, canorm);
	//      Wjets_out[i]=canorm*Wjets_out[i] .....
	//      please note that you have to do the 25% change per jetbin per flavor, which is NOT available in the long list. For these check you can use idsys=+/-2003,2004.
	//      Another note: you can only use this method when you have all the systematic changes to MC samples available. For example, when you apply the JESup systematic to your Wjets samples, you use:
	//      SetWflavors(-1, "jes_up" , ijet, Wjets_in_jesup, Wjets_out_jesup, canorm_jesup);
	//      Wjets_out_jesup[i]=canorm_jesup*Wjets_out_jesup[i];
	//      and thus (symbolically), the result is:
	//      Wjets_out_jesup= HFjesup x CAjesup x Wjets_in_jesup
	//      These individual systematics should not be applied to Wjets_nominal --> then you can use the short version as approximation.
	//      
	// SHORT:
	// Typical use when the combined uncertainties (short list)  are used:
	//      .....
	//      int ijet=3; // in case you input the 3jetbin
	//      int idsys=2000; // 
	//      double Wjets_in[4]={ 6888, 14281, 22937, 102347} ;
	//      double Wjets_out[4];
	//      SetWflavors(idsys, "" , ijet, Wjets_in, Wjets_out, canorm);
	//      Wjets_out[i]=canorm*Wjets_out[i] .....
	//      do this also for idsys=2001, 2002, 2003, 2004 and all negative: -2000, -2001, -2002, -2003, -2004. 
	//      if you also need to check the normalization (when you use canorm), then use 2005, -2005 as check.
	//      add all positive effects on you analysis in quadrature, and all negative effect also, as usual. Please note that negative effects on your analsysis
	//      do not neccessarily correspond with the negative idsys.
	//      See also remark * for checks 2003 and 2004.
	//      .....
	//
	// FAQ: I don't have the PRETAGGED numbers available at all times, but I want to apply the fraction to my TAGGED results.
	//      DON'T use this function 'SetWflavors' with your number of TAGGED events.
	//      If possible, once call the function SetWflavor with your PRETAGGED numbers and store the ratio's Wjets_in/Wjets_out -->
	//      these ratio's should be applied to your TAGGED samples.
	//
	// get HF Kfactors
	_mode=0; // should be ZERO
	if (_mode!=0) cout << " SetWFlavors: mode should be zero " << _mode << endl;
	double Hbb[9], Hcc[9], Hc[9],  Hll[9];
	double Fbb[9], Fcc[9], Fc[9],  Fll[9];
	double cafac[9];
	double fbb[9], fcc[9], fc[9],  fll[9];
	double winpretag[9], wintag[9], woutpretag[9],  wouttag[9];
        if (ijet==-3) { ijet=6; }
        if (ijet==-4) { ijet=7; }
        if (ijet==-5) { ijet=5;  } // default inc.
        if (ijet<=0 || ijet>7) cout << "FATAL  HF KFACTORS dont know which jetbin you are going to use " << endl;

	if (idsys*idsys<999*999) { // go for the full list of systematics
		GetFFactors_muon( idsys, sysname,  Hbb,Hcc, Hc,  Hll, cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
	} else if (idsys*idsys<1999*1999) {
		cout << " FATAL  HF KFACTORS probabbly called with old defintion of systematics ids " << endl;
	} else  {                  // go for the reduced errors.
		GetFFactors_muon( 0, "",  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
	}

	// errors on fraction, not on K factors
	double sigma_hf[2]={ 0.19 ,  - 0.18}; // error on ttoal HF
	double sigma_as[2] ={ 0.17 ,  - 0.18}; // error on fbb (and fcc) wrt to fc
	double sigma_ll[2] ={0.,  - 0.}; // -100% correlated with error on total HF
	double sigma_canorm_up[9]=  {0,  0.064, 0.053, 0.083, 0.085, 0.13, 0.071,0.077,0};
	double sigma_canorm_down[9]={0, -0.058, -0.049,-0.072,-0.079,-0.15, -0.065,-0.075,0};


	// use 2jetbin as a base for the default user.
	Fbb[ijet]=Hbb[2];
	Fcc[ijet]=Hcc[2];
	Fc[ijet]=Hc[2];
	Fll[ijet]=Hll[2];
	canorm=cafac[ijet];
	Wjets_out[0]=Fbb[ijet]*Wjets_in[0];
	Wjets_out[1]=Fcc[ijet]*Wjets_in[1];
	Wjets_out[2]=Fc[ijet] *Wjets_in[2];
	Wjets_out[3]=Fll[ijet]*Wjets_in[3];
	double nw_in=0;
	double nw_raw=0;
	for (int iw=0;iw<4;iw++) {
		nw_in+=Wjets_in[iw];
		nw_raw+=Wjets_out[iw];
	}
	if (nw_raw==0) { nw_in=1; nw_raw=1;} // some protection
	for (int iw=0;iw<4;iw++) {
		Wjets_out[iw]*=nw_in/nw_raw;
	}
	// done setting HFs
	// now systematics if combined mode is used.
	double fbb_in=Wjets_in[0]/nw_in;
	double fcc_in=Wjets_in[1]/nw_in;
	double fc_in= Wjets_in[2]/nw_in;
	double fll_in=Wjets_in[3]/nw_in;
	double fhf_in=fbb_in+fcc_in+fc_in;
	//   cout << " fraction in " << fbb_in <<  " " << fcc_in << " " << fc_in << " " << fll_in << endl;
	double fbb_out=Wjets_out[0]/nw_in;
	double fcc_out=Wjets_out[1]/nw_in;
	double fc_out= Wjets_out[2]/nw_in;
	double fll_out=Wjets_out[3]/nw_in;
	double fhf_out=fbb_out+fcc_out+fc_out;

	// total HF change
	if (idsys==2000) {
		fbb_out*=1+sigma_hf[0];
		fcc_out*=1+sigma_hf[0];
		fc_out*=1+sigma_hf[0];
		fhf_out=fbb_out+fcc_out+fc_out;
		fll_out=1-fhf_out;
	} else if (idsys==-2000) {
		fbb_out*=1+sigma_hf[1];
		fcc_out*=1+sigma_hf[1];
		fc_out*=1+sigma_hf[1];
		fhf_out=fbb_out+fcc_out+fc_out;
		fll_out=1-fhf_out;
	}

	// total HF constant, but bb+cc change wrt c
	if (idsys==2001) {
		fbb_out*=1+sigma_as[0];
		fcc_out*=1+sigma_as[0];
		fc_out=fhf_out-fbb_out-fcc_out; 
		// fll_out unchanged
	} else if (idsys==-2001) {
		fbb_out*=1+sigma_as[1];
		fcc_out*=1+sigma_as[1];
		fc_out=fhf_out-fbb_out-fcc_out;
		// fll_out unchanged
	}

	if (idsys==2002) {
		// dummy
	} else if (idsys==-2002) {
		// dummy
	} 
	TString cname;
	if (ijet!=2 && idsys==2003) { // change bb,cc at the cost of fc+fl
		fbb_out*=1+0.25;
		fcc_out*=1+0.25;
		double rest=fc_out+fll_out;
		double rest_new=1-fbb_out-fcc_out;
		fc_out*=rest_new/rest;
		fll_out*=rest_new/rest;
		if (ijet==1) cname="WbbWccjet1_up";
		if (ijet==3 || ijet==6) cname="WbbWccjet3_up";
		if (ijet==4 || ijet==7) cname="WbbWccjet4_up";
		if (ijet==5) cname="WbbWccjet5_up"; 
		GetFFactors_muon( -1, cname,  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		canorm=cafac[ijet];
	} else if (ijet!=2 && idsys==-2003) {
		fbb_out*=1-0.25;
		fcc_out*=1-0.25;
		double rest=fc_out+fll_out;
		double rest_new=1-fbb_out-fcc_out;
		fc_out*=rest_new/rest;
		fll_out*=rest_new/rest;
		if (ijet==1) cname="WbbWccjet1_down";
		if (ijet==3 || ijet==6) cname="WbbWccjet3_down";
		if (ijet==4 || ijet==7) cname="WbbWccjet4_down";
		if (ijet==5) cname="WbbWccjet5_down"; 
		GetFFactors_muon( -1, cname,  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		canorm=cafac[ijet];
	}
	if (ijet!=2 && idsys==2004) { // change c. at the cost of fbb+fcc+fl
		fc_out*=1+0.25;
		double rest=fbb_out+fcc_out+fll_out;
		double rest_new=1-fc_out; 
		fbb_out*=rest_new/rest;
		fcc_out*=rest_new/rest;
		fll_out*=rest_new/rest;
		if (ijet==1) cname="Wcjet1_up";
		if (ijet==3 || ijet==6) cname="Wcjet3_up";
		if (ijet==4 || ijet==7) cname="Wcjet4_up";
		if (ijet==5) cname="Wcjet5_up"; // 
		GetFFactors_muon( -1, cname,  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		canorm=cafac[ijet];
	} else if (ijet!=2 && idsys==-2004) {
		fc_out*=1-0.25;
		double rest=fbb_out+fcc_out+fll_out;
		double rest_new=1-fc_out; 
		fbb_out*=rest_new/rest;
		fcc_out*=rest_new/rest;
		fll_out*=rest_new/rest;
		if (ijet==1) cname="Wcjet1_down";
		if (ijet==3 || ijet==6) cname="Wcjet3_down";
		if (ijet==4 || ijet==7) cname="Wcjet4_down";
		if (ijet==5) cname="Wcjet5_down"; 
		GetFFactors_muon( -1, cname,  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		canorm=cafac[ijet];
	}

	fhf_out=fbb_out+fcc_out+fc_out;
	if (fabs(fll_out+fhf_out -1) > 0.001) cout << " FATAL this shoud be one " << fll_out+fhf_out << endl;
	if (fbb_out<0 || fcc_out<0 || fc_out<0 || fll_out<0) cout << " FATAL negative fraction " << fbb_out << " " << fcc_out << " " << fc_out << " " << fll_out << endl;

	if (idsys==2005) { // change c. at the cost of fbb+fcc+fl
		canorm*=1+sigma_canorm_up[ijet]; }
	else if (idsys==-2005) {
		canorm*=1+sigma_canorm_down[ijet];
	}
	if (abs(idsys)>2005) cout << " FATAL: unknown systematic check for HF factors " ;

	// calculate new Wjets_out (keeps normalization constant)
	Wjets_out[0]=fbb_out*nw_in;
	Wjets_out[1]=fcc_out*nw_in;
	Wjets_out[2]=fc_out*nw_in;
	Wjets_out[3]=fll_out*nw_in;

	if (_mode==1) {
		// mode for experts.
		// check PRETAGGED systematics long list
		// special mode to test long list systematics by defining: Wsys= HF x CA x Wnominal (default: HF x CA x Wsys)
		double nwsys=winpretag[ijet];
		double fbbsys=fbb[ijet];
		double fccsys=fcc[ijet];
		double fcsys=fc[ijet];
		double fllsys=fll[ijet];
		GetFFactors_muon( 0, "",  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		double nwnom=winpretag[ijet];	
		double fbbnom=fbb[ijet];
		double fccnom=fcc[ijet];
		double fcnom=fc[ijet];
		double fllnom=fll[ijet];
		canorm*=nwsys/nwnom;
		if (idsys<-1999 || idsys>1999) {
			fbbsys=fbb_out;
			fccsys=fcc_out;
			fcsys=fc_out;
			fllsys=fll_out;
		}
		Wjets_out[0]=fbbsys*nw_in;
		Wjets_out[1]=fccsys*nw_in;
		Wjets_out[2]=fcsys*nw_in;
		Wjets_out[3]=fllsys*nw_in;
		//  renormalize  to ensure that the overal normalization remains constant
		 nw_in=0;
		 nw_raw=0;
		for (int iw=0;iw<4;iw++) {
			nw_in+=Wjets_in[iw];
			nw_raw+=Wjets_out[iw];
		}
		if (nw_raw==0) { nw_in=0; nw_raw=1;} // some protection
		for (int iw=0;iw<4;iw++) {
			Wjets_out[iw]*=nw_in/nw_raw;
		}
	} // mode=1

	if (_mode==2) { 
		// mode for experts
		// check TAGGED systematics
		// special mode to test long list systematics by defining: Wsys= HF x CA x Wnominal (default: HF x CA x Wsys)
		double nwsys=wintag[ijet];
		double nwpresys=winpretag[ijet];
		double fbbsys=fbb[ijet];
		double fccsys=fcc[ijet];
		double fcsys=fc[ijet];
		double fllsys=fll[ijet];
		GetFFactors_muon( 0, "",  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		double nwnom=wintag[ijet];
		double nwprenom=winpretag[ijet];

		double fbbnom=fbb[ijet];
		double fccnom=fcc[ijet];
		double fcnom=fc[ijet];
		double fllnom=fll[ijet];
		canorm*=nwpresys/nwprenom/nwpresys * nwsys;
		//		canorm x winpretag_sys  --> canorm x winpretag_nom      x winpretag_sys/winpretag_nom      /nwpresys*nwsys
		if (idsys<-1999 || idsys>1999) {
			fbbsys=fbb_out;
			fccsys=fcc_out;
			fcsys=fc_out;
			fllsys=fll_out;
		}
		Wjets_out[0]=fbbsys*nw_in;
		Wjets_out[1]=fccsys*nw_in;
		Wjets_out[2]=fcsys*nw_in;
		Wjets_out[3]=fllsys*nw_in;
		//  renormalize  to ensure that the overal normalization remains constant
		 nw_in=0;
		 nw_raw=0;
		for (int iw=0;iw<4;iw++) {
			nw_in+=Wjets_in[iw];
			nw_raw+=Wjets_out[iw];
		}
		if (nw_raw==0) { nw_in=0; nw_raw=1;} // some protection
		for (int iw=0;iw<4;iw++) {
			Wjets_out[iw]*=nw_in/nw_raw;
		}
	} //mode=2
}





void GetFFactors_muon(int idsys, TString sysname, double Fbb[9],double Fcc[9],double Fc[9], double Fll[9], double cafac[9],
	double fbb[9],double fcc[9],double fc[9], double fll[9],
	double winpretag[9],double wintag[9],double woutpretag[9], double wouttag[9]) {
		//--
		// This function should NOT be called directly by the user. It is called via SetWflavors
		//
		// Ffactors for r17. ttbar cuts april 2012.  MUON
		// the Ffactors Fbb, Fcc, Fc, Fll contain weight factors for Wbb events, Wcc events, Wc event and Wlight events. 
		// These factors do NOT contains the normalization of MC Wjets events to data.
		//--
		// input: int idsys --> the id of the systematic variation, then call with sysname="" 
		//        TString sysname --> the name of the systematic variation, then call with idsys=-1 (NOT: idsys=0, because that is nominal).
		// output: Fbb[9],Fcc[9],Fc[9],Fll[9]
		//         position Fxx[0] is not used, position 1,2,3,4 are for the jetbins 1,2,3,4exclusive.
		//         position Fxx[5,6,7] for 5,3,4inclusive respectively. 
		//         HOWEVER, ONLY the 2jetbin should be used in all jetbins while keeping the pretagged normalization constant.
		//         cafac[9] --> the normalization obtained from Charge Assymmetry to be used in each jetbin (same for tagged and pretagged)
		// typical usage case:
		// You call this function once for nominal and then again for each systematic variation under study.
		// The input is either the integer idsys or the TString sysname. For the first mode with idsys, you typically may use:
		//      .....
		//      int idsys; 
		//      double Fbb[9], Fcc[9], Fc[9],  Fll[9];
		//      idsys=0; // for nominal
		//      GetFFactors_muon(idsys, "" , Fbb,Fcc,Fc,Fll); 
		//      .....
		//--
		// If you want to read this function by 'eye', i have two hints:
		// The HF factors are Fbb[2], Fcc[2] , Fc[2], Fll[2] for EACH jetbin, defined on the PRETAGGED sample.
		// The normalization factors are cafac[1], cafac[2], cafac[3], cafac[4], cafac[5] for jetbins 1,2,3,4,5 resp.
		// Ignore the other quantities which are used for additional checks and studies.
		//--
		double fmcbb[9];double fmccc[9];double fmcc[9];double fmcll[9];	



      Fbb[1]= 1.23358;   Fcc[1]= 1.23358;   Fc[1]= 0.961029;   Fll[1]= 0.993059; 
         cafac[1]= 0.9976;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.2387e+06;   wouttag[1]= 80628.3;  
       fbb[1]= 0.01471;   fcc[1]= 0.0426761;   fc[1]= 0.127449;   fll[1]= 0.815165;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.91873;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37108.9;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20558;   Fcc[3]= 1.20558;   Fc[3]= 0.93922;   Fll[3]= 0.970523; 
         cafac[3]= 0.832384;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95780.3;   wouttag[3]= 11676.5;  
       fbb[3]= 0.05821;   fcc[3]= 0.117277;   fc[3]= 0.142185;   fll[3]= 0.682328;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19366;   Fcc[4]= 1.19366;   Fc[4]= 0.929929;   Fll[4]= 0.960923; 
         cafac[4]= 0.888562;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22242.5;   wouttag[4]= 3569.19;  
       fbb[4]= 0.0827553;   fcc[4]= 0.140318;   fc[4]= 0.132517;   fll[4]= 0.644409;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18386;   Fcc[5]= 1.18386;   Fc[5]= 0.922294;   Fll[5]= 0.953033; 
         cafac[5]= 0.824143;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.68;   wouttag[5]= 1216.94;  
       fbb[5]= 0.104276;   fcc[5]= 0.156627;   fc[5]= 0.117107;   fll[5]= 0.621991;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.2025;   Fcc[6]= 1.2025;   Fc[6]= 0.93682;   Fll[6]= 0.968043; 
         cafac[6]= 0.842063;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123775;   wouttag[6]= 16452.6;  
       fbb[6]= 0.0646144;   fcc[6]= 0.123104;   fc[6]= 0.139332;   fll[6]= 0.67295;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.19153;   Fcc[7]= 1.19153;   Fc[7]= 0.928271;   Fll[7]= 0.959209; 
         cafac[7]= 0.874032;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27901.2;   wouttag[7]= 4803.54;  
       fbb[7]= 0.0874306;   fcc[7]= 0.143861;   fc[7]= 0.129169;   fll[7]= 0.639539;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
  if ( idsys==0   || sysname=="norminal" ) { 
// OK
 }  
 else if ( idsys==1   || sysname=="tt_up" ) { 
      Fbb[1]= 1.23137;   Fcc[1]= 1.23137;   Fc[1]= 0.915944;   Fll[1]= 1.00047; 
         cafac[1]= 0.989717;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.22101e+06;   wouttag[1]= 78090;  
       fbb[1]= 0.0146837;   fcc[1]= 0.0425999;   fc[1]= 0.121469;   fll[1]= 0.821247;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.2182;   Fcc[2]= 1.2182;   Fc[2]= 0.906142;   Fll[2]= 0.989762; 
         cafac[2]= 0.91017;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 468774;   wouttag[2]= 36215.4;  
       fbb[2]= 0.0360601;   fcc[2]= 0.0866484;   fc[2]= 0.138403;   fll[2]= 0.738888;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20571;   Fcc[3]= 1.20571;   Fc[3]= 0.896856;   Fll[3]= 0.979618; 
         cafac[3]= 0.824755;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 94902.5;   wouttag[3]= 11448.6;  
       fbb[3]= 0.0582162;   fcc[3]= 0.11729;   fc[3]= 0.135771;   fll[3]= 0.688723;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19368;   Fcc[4]= 1.19368;   Fc[4]= 0.887908;   Fll[4]= 0.969845; 
         cafac[4]= 0.880673;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22045;   wouttag[4]= 3508.34;  
       fbb[4]= 0.082757;   fcc[4]= 0.140321;   fc[4]= 0.126529;   fll[4]= 0.650393;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.1833;   Fcc[5]= 1.1833;   Fc[5]= 0.880187;   Fll[5]= 0.961412; 
         cafac[5]= 0.817208;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5630.89;   wouttag[5]= 1199.73;  
       fbb[5]= 0.104227;   fcc[5]= 0.156554;   fc[5]= 0.11176;   fll[5]= 0.627459;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.20258;   Fcc[6]= 1.20258;   Fc[6]= 0.894527;   Fll[6]= 0.977074; 
         cafac[6]= 0.834421;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 122651;   wouttag[6]= 16146.5;  
       fbb[6]= 0.0646186;   fcc[6]= 0.123112;   fc[6]= 0.133042;   fll[6]= 0.679228;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.19143;   Fcc[7]= 1.19143;   Fc[7]= 0.88623;   Fll[7]= 0.968012; 
         cafac[7]= 0.86636;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27656.3;   wouttag[7]= 4725.42;  
       fbb[7]= 0.0874231;   fcc[7]= 0.143849;   fc[7]= 0.123319;   fll[7]= 0.645409;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==2   || sysname=="tt_down" ) { 
      Fbb[1]= 1.2359;   Fcc[1]= 1.2359;   Fc[1]= 1.00853;   Fll[1]= 0.985253; 
         cafac[1]= 1.00604;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.25764e+06;   wouttag[1]= 83346.8;  
       fbb[1]= 0.0147376;   fcc[1]= 0.0427564;   fc[1]= 0.133748;   fll[1]= 0.808758;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21877;   Fcc[2]= 1.21877;   Fc[2]= 0.994553;   Fll[2]= 0.971596; 
         cafac[2]= 0.927895;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 477903;   wouttag[2]= 38065.6;  
       fbb[2]= 0.036077;   fcc[2]= 0.086689;   fc[2]= 0.151907;   fll[2]= 0.725327;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20545;   Fcc[3]= 1.20545;   Fc[3]= 0.983683;   Fll[3]= 0.960977; 
         cafac[3]= 0.840544;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 96719.3;   wouttag[3]= 11920.4;  
       fbb[3]= 0.0582034;   fcc[3]= 0.117264;   fc[3]= 0.148916;   fll[3]= 0.675617;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19363;   Fcc[4]= 1.19363;   Fc[4]= 0.97404;   Fll[4]= 0.951557; 
         cafac[4]= 0.896997;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22453.6;   wouttag[4]= 3634.25;  
       fbb[4]= 0.0827533;   fcc[4]= 0.140315;   fc[4]= 0.138803;   fll[4]= 0.638129;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18444;   Fcc[5]= 1.18444;   Fc[5]= 0.966539;   Fll[5]= 0.944229; 
         cafac[5]= 0.831557;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5729.76;   wouttag[5]= 1235.34;  
       fbb[5]= 0.104327;   fcc[5]= 0.156704;   fc[5]= 0.122724;   fll[5]= 0.616245;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.20242;   Fcc[6]= 1.20242;   Fc[6]= 0.981213;   Fll[6]= 0.958564; 
         cafac[6]= 0.850238;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 124976;   wouttag[6]= 16780;  
       fbb[6]= 0.0646099;   fcc[6]= 0.123095;   fc[6]= 0.145935;   fll[6]= 0.66636;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.19163;   Fcc[7]= 1.19163;   Fc[7]= 0.972411;   Fll[7]= 0.949966; 
         cafac[7]= 0.882235;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 28163;   wouttag[7]= 4887.06;  
       fbb[7]= 0.0874383;   fcc[7]= 0.143874;   fc[7]= 0.135312;   fll[7]= 0.633376;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==19   || sysname=="Kbb_up" ) { 
      Fbb[1]= 1.33357;   Fcc[1]= 1.33357;   Fc[1]= 0.956287;   Fll[1]= 0.988158; 
         cafac[1]= 0.992677;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.2387e+06;   wouttag[1]= 81467.3;  
       fbb[1]= 0.0159023;   fcc[1]= 0.0461353;   fc[1]= 0.12682;   fll[1]= 0.811143;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.30988;   Fcc[2]= 1.30988;   Fc[2]= 0.939301;   Fll[2]= 0.970607; 
         cafac[2]= 0.909087;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37889.8;  
       fbb[2]= 0.0387741;   fcc[2]= 0.0931698;   fc[2]= 0.143468;   fll[2]= 0.724588;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.2902;   Fcc[3]= 1.2902;   Fc[3]= 0.925189;   Fll[3]= 0.956024; 
         cafac[3]= 0.819949;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95780.3;   wouttag[3]= 11971.9;  
       fbb[3]= 0.0622956;   fcc[3]= 0.125509;   fc[3]= 0.140061;   fll[3]= 0.672135;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.27229;   Fcc[4]= 1.27229;   Fc[4]= 0.912342;   Fll[4]= 0.942748; 
         cafac[4]= 0.871757;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22242.5;   wouttag[4]= 3671.98;  
       fbb[4]= 0.0882064;   fcc[4]= 0.149561;   fc[4]= 0.130011;   fll[4]= 0.632222;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.25781;   Fcc[5]= 1.25781;   Fc[5]= 0.901958;   Fll[5]= 0.932019; 
         cafac[5]= 0.805971;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.68;   wouttag[5]= 1248.32;  
       fbb[5]= 0.110789;   fcc[5]= 0.16641;   fc[5]= 0.114524;   fll[5]= 0.608276;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.28557;   Fcc[6]= 1.28557;   Fc[6]= 0.921866;   Fll[6]= 0.95259; 
         cafac[6]= 0.828621;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123775;   wouttag[6]= 16887.3;  
       fbb[6]= 0.0690776;   fcc[6]= 0.131607;   fc[6]= 0.137108;   fll[6]= 0.662207;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.26913;   Fcc[7]= 1.26913;   Fc[7]= 0.91008;   Fll[7]= 0.940412; 
         cafac[7]= 0.856905;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27901.2;   wouttag[7]= 4938.79;  
       fbb[7]= 0.0931248;   fcc[7]= 0.153231;   fc[7]= 0.126638;   fll[7]= 0.627006;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==20   || sysname=="Kbb_down" ) { 
      Fbb[1]= 1.13259;   Fcc[1]= 1.13259;   Fc[1]= 0.965819;   Fll[1]= 0.998008; 
         cafac[1]= 1.00257;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.2387e+06;   wouttag[1]= 79780.9;  
       fbb[1]= 0.0135057;   fcc[1]= 0.0391824;   fc[1]= 0.128084;   fll[1]= 0.819228;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.12511;   Fcc[2]= 1.12511;   Fc[2]= 0.959441;   Fll[2]= 0.991417; 
         cafac[2]= 0.928579;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 36311.3;  
       fbb[2]= 0.0333046;   fcc[2]= 0.0800273;   fc[2]= 0.146544;   fll[2]= 0.740124;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.11836;   Fcc[3]= 1.11836;   Fc[3]= 0.953683;   Fll[3]= 0.985468; 
         cafac[3]= 0.845202;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95780.3;   wouttag[3]= 11372.1;  
       fbb[3]= 0.0539984;   fcc[3]= 0.108792;   fc[3]= 0.144374;   fll[3]= 0.692835;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.11194;   Fcc[4]= 1.11194;   Fc[4]= 0.948209;   Fll[4]= 0.979811; 
         cafac[4]= 0.906028;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22242.5;   wouttag[4]= 3462.37;  
       fbb[4]= 0.0770898;   fcc[4]= 0.130712;   fc[4]= 0.135122;   fll[4]= 0.657076;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.1065;   Fcc[5]= 1.1065;   Fc[5]= 0.943569;   Fll[5]= 0.975017; 
         cafac[5]= 0.843153;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.68;   wouttag[5]= 1184.11;  
       fbb[5]= 0.0974619;   fcc[5]= 0.146392;   fc[5]= 0.119808;   fll[5]= 0.636338;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.1167;   Fcc[6]= 1.1167;   Fc[6]= 0.952268;   Fll[6]= 0.984006; 
         cafac[6]= 0.855949;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123775;   wouttag[6]= 16003.6;  
       fbb[6]= 0.0600039;   fcc[6]= 0.11432;   fc[6]= 0.14163;   fll[6]= 0.684047;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.11076;   Fcc[7]= 1.11076;   Fc[7]= 0.947203;   Fll[7]= 0.978772; 
         cafac[7]= 0.891859;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27901.2;   wouttag[7]= 4662.77;  
       fbb[7]= 0.081504;   fcc[7]= 0.134109;   fc[7]= 0.131804;   fll[7]= 0.652583;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==21   || sysname=="Kc_up" ) { 
      Fbb[1]= 1.22464;   Fcc[1]= 1.22464;   Fc[1]= 1.0087;   Fll[1]= 0.985864; 
         cafac[1]= 0.990372;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.2387e+06;   wouttag[1]= 82537.5;  
       fbb[1]= 0.0146034;   fcc[1]= 0.0423669;   fc[1]= 0.13377;   fll[1]= 0.809259;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.20844;   Fcc[2]= 1.20844;   Fc[2]= 0.995358;   Fll[2]= 0.972825; 
         cafac[2]= 0.911165;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37597.3;  
       fbb[2]= 0.0357714;   fcc[2]= 0.0859546;   fc[2]= 0.15203;   fll[2]= 0.726244;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.19585;   Fcc[3]= 1.19585;   Fc[3]= 0.984983;   Fll[3]= 0.962685; 
         cafac[3]= 0.825661;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95780.3;   wouttag[3]= 11770.3;  
       fbb[3]= 0.0577398;   fcc[3]= 0.11633;   fc[3]= 0.149113;   fll[3]= 0.676818;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.18467;   Fcc[4]= 1.18467;   Fc[4]= 0.975776;   Fll[4]= 0.953686; 
         cafac[4]= 0.88187;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22242.5;   wouttag[4]= 3588.17;  
       fbb[4]= 0.082132;   fcc[4]= 0.139261;   fc[4]= 0.13905;   fll[4]= 0.639556;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.17597;   Fcc[5]= 1.17597;   Fc[5]= 0.968613;   Fll[5]= 0.946685; 
         cafac[5]= 0.818653;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.68;   wouttag[5]= 1220.69;  
       fbb[5]= 0.103581;   fcc[5]= 0.155584;   fc[5]= 0.122988;   fll[5]= 0.617848;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.19298;   Fcc[6]= 1.19298;   Fc[6]= 0.982626;   Fll[6]= 0.960381; 
         cafac[6]= 0.835398;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123775;   wouttag[6]= 16568.3;  
       fbb[6]= 0.0641029;   fcc[6]= 0.122129;   fc[6]= 0.146145;   fll[6]= 0.667623;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.18278;   Fcc[7]= 1.18278;   Fc[7]= 0.974221;   Fll[7]= 0.952166; 
         cafac[7]= 0.867615;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27901.2;   wouttag[7]= 4825.96;  
       fbb[7]= 0.0867886;   fcc[7]= 0.142805;   fc[7]= 0.135563;   fll[7]= 0.634843;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==22   || sysname=="Kc_down" ) { 
      Fbb[1]= 1.24265;   Fcc[1]= 1.24265;   Fc[1]= 0.912658;   Fll[1]= 1.00036; 
         cafac[1]= 1.00493;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.2387e+06;   wouttag[1]= 78691.1;  
       fbb[1]= 0.0148181;   fcc[1]= 0.0429898;   fc[1]= 0.121034;   fll[1]= 0.821158;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.22868;   Fcc[2]= 1.22868;   Fc[2]= 0.902398;   Fll[2]= 0.989114; 
         cafac[2]= 0.926422;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 36612.4;  
       fbb[2]= 0.0363703;   fcc[2]= 0.0873938;   fc[2]= 0.137832;   fll[2]= 0.738404;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.21548;   Fcc[3]= 1.21548;   Fc[3]= 0.892706;   Fll[3]= 0.97849; 
         cafac[3]= 0.839217;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95780.3;   wouttag[3]= 11581.2;  
       fbb[3]= 0.0586878;   fcc[3]= 0.11824;   fc[3]= 0.135143;   fll[3]= 0.687929;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.20279;   Fcc[4]= 1.20279;   Fc[4]= 0.883382;   Fll[4]= 0.96827; 
         cafac[4]= 0.895356;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22242.5;   wouttag[4]= 3549.92;  
       fbb[4]= 0.083388;   fcc[4]= 0.141391;   fc[4]= 0.125884;   fll[4]= 0.649337;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.19185;   Fcc[5]= 1.19185;   Fc[5]= 0.875351;   Fll[5]= 0.959467; 
         cafac[5]= 0.829707;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.68;   wouttag[5]= 1213.14;  
       fbb[5]= 0.10498;   fcc[5]= 0.157684;   fc[5]= 0.111146;   fll[5]= 0.62619;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.21217;   Fcc[6]= 1.21217;   Fc[6]= 0.890278;   Fll[6]= 0.975829; 
         cafac[6]= 0.848836;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123775;   wouttag[6]= 16335.1;  
       fbb[6]= 0.065134;   fcc[6]= 0.124094;   fc[6]= 0.13241;   fll[6]= 0.678362;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.20041;   Fcc[7]= 1.20041;   Fc[7]= 0.881636;   Fll[7]= 0.966356; 
         cafac[7]= 0.880545;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27901.2;   wouttag[7]= 4780.78;  
       fbb[7]= 0.0880821;   fcc[7]= 0.144933;   fc[7]= 0.12268;   fll[7]= 0.644305;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==23   || sysname=="Kll_up" ) { 
      Fbb[1]= 1.2281;   Fcc[1]= 1.2281;   Fc[1]= 0.956759;   Fll[1]= 0.994059; 
         cafac[1]= 0.998605;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.2387e+06;   wouttag[1]= 80392.1;  
       fbb[1]= 0.0146446;   fcc[1]= 0.0424865;   fc[1]= 0.126882;   fll[1]= 0.815987;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21361;   Fcc[2]= 1.21361;   Fc[2]= 0.945473;   Fll[2]= 0.982334; 
         cafac[2]= 0.920071;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37012.2;  
       fbb[2]= 0.0359243;   fcc[2]= 0.0863221;   fc[2]= 0.144411;   fll[2]= 0.733343;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.2011;   Fcc[3]= 1.2011;   Fc[3]= 0.935724;   Fll[3]= 0.972205; 
         cafac[3]= 0.833826;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95780.3;   wouttag[3]= 11648.6;  
       fbb[3]= 0.0579933;   fcc[3]= 0.116841;   fc[3]= 0.141655;   fll[3]= 0.683511;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.18946;   Fcc[4]= 1.18946;   Fc[4]= 0.92666;   Fll[4]= 0.962787; 
         cafac[4]= 0.890286;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22242.5;   wouttag[4]= 3560.76;  
       fbb[4]= 0.0824643;   fcc[4]= 0.139825;   fc[4]= 0.132051;   fll[4]= 0.64566;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.17984;   Fcc[5]= 1.17984;   Fc[5]= 0.919164;   Fll[5]= 0.954999; 
         cafac[5]= 0.825843;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.68;   wouttag[5]= 1214.55;  
       fbb[5]= 0.103922;   fcc[5]= 0.156095;   fc[5]= 0.116709;   fll[5]= 0.623274;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.19809;   Fcc[6]= 1.19809;   Fc[6]= 0.933381;   Fll[6]= 0.96977; 
         cafac[6]= 0.843566;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123775;   wouttag[6]= 16413.6;  
       fbb[6]= 0.0643772;   fcc[6]= 0.122652;   fc[6]= 0.138821;   fll[6]= 0.67415;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.18737;   Fcc[7]= 1.18737;   Fc[7]= 0.925032;   Fll[7]= 0.961095; 
         cafac[7]= 0.875751;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27901.2;   wouttag[7]= 4792.68;  
       fbb[7]= 0.0871255;   fcc[7]= 0.143359;   fc[7]= 0.128719;   fll[7]= 0.640797;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==24   || sysname=="Kll_down" ) { 
      Fbb[1]= 1.23911;   Fcc[1]= 1.23911;   Fc[1]= 0.965338;   Fll[1]= 0.992049; 
         cafac[1]= 0.996586;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.2387e+06;   wouttag[1]= 80866.7;  
       fbb[1]= 0.0147759;   fcc[1]= 0.0428674;   fc[1]= 0.12802;   fll[1]= 0.814337;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.22338;   Fcc[2]= 1.22338;   Fc[2]= 0.953085;   Fll[2]= 0.979458; 
         cafac[2]= 0.917378;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37206.4;  
       fbb[2]= 0.0362135;   fcc[2]= 0.0870171;   fc[2]= 0.145574;   fll[2]= 0.731196;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.2101;   Fcc[3]= 1.2101;   Fc[3]= 0.942742;   Fll[3]= 0.968829; 
         cafac[3]= 0.830931;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95780.3;   wouttag[3]= 11704.7;  
       fbb[3]= 0.0584282;   fcc[3]= 0.117717;   fc[3]= 0.142718;   fll[3]= 0.681137;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19788;   Fcc[4]= 1.19788;   Fc[4]= 0.933222;   Fll[4]= 0.959045; 
         cafac[4]= 0.886826;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22242.5;   wouttag[4]= 3577.69;  
       fbb[4]= 0.0830483;   fcc[4]= 0.140815;   fc[4]= 0.132986;   fll[4]= 0.64315;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.1879;   Fcc[5]= 1.1879;   Fc[5]= 0.925446;   Fll[5]= 0.951054; 
         cafac[5]= 0.822431;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.68;   wouttag[5]= 1219.34;  
       fbb[5]= 0.104632;   fcc[5]= 0.157162;   fc[5]= 0.117507;   fll[5]= 0.620699;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.20695;   Fcc[6]= 1.20695;   Fc[6]= 0.940285;   Fll[6]= 0.966303; 
         cafac[6]= 0.84055;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123775;   wouttag[6]= 16491.9;  
       fbb[6]= 0.0648533;   fcc[6]= 0.123559;   fc[6]= 0.139847;   fll[6]= 0.67174;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.19572;   Fcc[7]= 1.19572;   Fc[7]= 0.931533;   Fll[7]= 0.957309; 
         cafac[7]= 0.872301;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27901.2;   wouttag[7]= 4814.48;  
       fbb[7]= 0.0877378;   fcc[7]= 0.144367;   fc[7]= 0.129623;   fll[7]= 0.638272;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==33   || sysname=="WbbWccjet1_up" ) { 
      Fbb[1]= 1.23007;   Fcc[1]= 1.23007;   Fc[1]= 0.958299;   Fll[1]= 0.990237; 
         cafac[1]= 0.993943;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 82120.3 ;   woutpretag[1]= 2.23049e+06;   wouttag[1]= 82874.4;  
       fbb[1]= 0.0183352;   fcc[1]= 0.0531936;   fc[1]= 0.125536;   fll[1]= 0.802935;  
       fmcbb[1]= 0.0149058;   fmccc[1]= 0.0432442;   fmcc[1]= 0.130999;   fmcll[1]= 0.810851;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.91873;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37108.9;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20558;   Fcc[3]= 1.20558;   Fc[3]= 0.93922;   Fll[3]= 0.970523; 
         cafac[3]= 0.832384;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95780.3;   wouttag[3]= 11676.5;  
       fbb[3]= 0.05821;   fcc[3]= 0.117277;   fc[3]= 0.142185;   fll[3]= 0.682328;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19366;   Fcc[4]= 1.19366;   Fc[4]= 0.929929;   Fll[4]= 0.960923; 
         cafac[4]= 0.888562;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22242.5;   wouttag[4]= 3569.19;  
       fbb[4]= 0.0827553;   fcc[4]= 0.140318;   fc[4]= 0.132517;   fll[4]= 0.644409;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18386;   Fcc[5]= 1.18386;   Fc[5]= 0.922294;   Fll[5]= 0.953033; 
         cafac[5]= 0.824143;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.68;   wouttag[5]= 1216.94;  
       fbb[5]= 0.104276;   fcc[5]= 0.156627;   fc[5]= 0.117107;   fll[5]= 0.621991;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.2025;   Fcc[6]= 1.2025;   Fc[6]= 0.93682;   Fll[6]= 0.968043; 
         cafac[6]= 0.842063;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123775;   wouttag[6]= 16452.6;  
       fbb[6]= 0.0646144;   fcc[6]= 0.123104;   fc[6]= 0.139332;   fll[6]= 0.67295;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.19153;   Fcc[7]= 1.19153;   Fc[7]= 0.928271;   Fll[7]= 0.959209; 
         cafac[7]= 0.874032;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27901.2;   wouttag[7]= 4803.54;  
       fbb[7]= 0.0874306;   fcc[7]= 0.143861;   fc[7]= 0.129169;   fll[7]= 0.639539;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==34   || sysname=="WbbWccjet1_down" ) { 
      Fbb[1]= 1.2371;   Fcc[1]= 1.2371;   Fc[1]= 0.963775;   Fll[1]= 0.995896; 
         cafac[1]= 1.0013;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 77943.4 ;   woutpretag[1]= 2.24701e+06;   wouttag[1]= 78352.7;  
       fbb[1]= 0.011064;   fcc[1]= 0.0320985;   fc[1]= 0.129372;   fll[1]= 0.827466;  
       fmcbb[1]= 0.00894348;   fmccc[1]= 0.0259465;   fmcc[1]= 0.134234;   fmcll[1]= 0.830876;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.91873;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37108.9;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20558;   Fcc[3]= 1.20558;   Fc[3]= 0.93922;   Fll[3]= 0.970523; 
         cafac[3]= 0.832384;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95780.3;   wouttag[3]= 11676.5;  
       fbb[3]= 0.05821;   fcc[3]= 0.117277;   fc[3]= 0.142185;   fll[3]= 0.682328;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19366;   Fcc[4]= 1.19366;   Fc[4]= 0.929929;   Fll[4]= 0.960923; 
         cafac[4]= 0.888562;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22242.5;   wouttag[4]= 3569.19;  
       fbb[4]= 0.0827553;   fcc[4]= 0.140318;   fc[4]= 0.132517;   fll[4]= 0.644409;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18386;   Fcc[5]= 1.18386;   Fc[5]= 0.922294;   Fll[5]= 0.953033; 
         cafac[5]= 0.824143;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.68;   wouttag[5]= 1216.94;  
       fbb[5]= 0.104276;   fcc[5]= 0.156627;   fc[5]= 0.117107;   fll[5]= 0.621991;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.2025;   Fcc[6]= 1.2025;   Fc[6]= 0.93682;   Fll[6]= 0.968043; 
         cafac[6]= 0.842063;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123775;   wouttag[6]= 16452.6;  
       fbb[6]= 0.0646144;   fcc[6]= 0.123104;   fc[6]= 0.139332;   fll[6]= 0.67295;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.19153;   Fcc[7]= 1.19153;   Fc[7]= 0.928271;   Fll[7]= 0.959209; 
         cafac[7]= 0.874032;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27901.2;   wouttag[7]= 4803.54;  
       fbb[7]= 0.0874306;   fcc[7]= 0.143861;   fc[7]= 0.129169;   fll[7]= 0.639539;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==35   || sysname=="WbbWccjet3_up" ) { 
      Fbb[1]= 1.23358;   Fcc[1]= 1.23358;   Fc[1]= 0.961029;   Fll[1]= 0.993059; 
         cafac[1]= 0.9976;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.2387e+06;   wouttag[1]= 80628.3;  
       fbb[1]= 0.01471;   fcc[1]= 0.0426761;   fc[1]= 0.127449;   fll[1]= 0.815165;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.91873;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37108.9;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.19512;   Fcc[3]= 1.19512;   Fc[3]= 0.931068;   Fll[3]= 0.962099; 
         cafac[3]= 0.824675;  
         winpretag[3]= 115067 ;   wintag[3]= 14302.9 ;   woutpretag[3]= 94893.3;   wouttag[3]= 12565.3;  
       fbb[3]= 0.0721309;   fcc[3]= 0.145324;   fc[3]= 0.134947;   fll[3]= 0.647598;  
       fmcbb[3]= 0.0603545;   fmccc[3]= 0.121598;   fmcc[3]= 0.144938;   fmcll[3]= 0.673109;  
      Fbb[4]= 1.19366;   Fcc[4]= 1.19366;   Fc[4]= 0.929929;   Fll[4]= 0.960923; 
         cafac[4]= 0.888562;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22242.5;   wouttag[4]= 3569.19;  
       fbb[4]= 0.0827553;   fcc[4]= 0.140318;   fc[4]= 0.132517;   fll[4]= 0.644409;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18386;   Fcc[5]= 1.18386;   Fc[5]= 0.922294;   Fll[5]= 0.953033; 
         cafac[5]= 0.824143;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.68;   wouttag[5]= 1216.94;  
       fbb[5]= 0.104276;   fcc[5]= 0.156627;   fc[5]= 0.117107;   fll[5]= 0.621991;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.19434;   Fcc[6]= 1.19434;   Fc[6]= 0.930459;   Fll[6]= 0.96147; 
         cafac[6]= 0.836252;  
         winpretag[6]= 146990 ;   wintag[6]= 19443.2 ;   woutpretag[6]= 122921;   wouttag[6]= 17340.1;  
       fbb[6]= 0.0754614;   fcc[6]= 0.145006;   fc[6]= 0.13369;   fll[6]= 0.645843;  
       fmcbb[6]= 0.0631826;   fmccc[6]= 0.121411;   fmcc[6]= 0.143681;   fmcll[6]= 0.671725;  
      Fbb[7]= 1.19153;   Fcc[7]= 1.19153;   Fc[7]= 0.928271;   Fll[7]= 0.959209; 
         cafac[7]= 0.874032;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27901.2;   wouttag[7]= 4803.54;  
       fbb[7]= 0.0874306;   fcc[7]= 0.143861;   fc[7]= 0.129169;   fll[7]= 0.639539;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==36   || sysname=="WbbWccjet3_down" ) { 
      Fbb[1]= 1.23358;   Fcc[1]= 1.23358;   Fc[1]= 0.961029;   Fll[1]= 0.993059; 
         cafac[1]= 0.9976;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.2387e+06;   wouttag[1]= 80628.3;  
       fbb[1]= 0.01471;   fcc[1]= 0.0426761;   fc[1]= 0.127449;   fll[1]= 0.815165;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.91873;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37108.9;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.21623;   Fcc[3]= 1.21623;   Fc[3]= 0.947516;   Fll[3]= 0.979096; 
         cafac[3]= 0.840378;  
         winpretag[3]= 115067 ;   wintag[3]= 12214.5 ;   woutpretag[3]= 96700.2;   wouttag[3]= 10754.9;  
       fbb[3]= 0.0440431;   fcc[3]= 0.0887349;   fc[3]= 0.14955;   fll[3]= 0.717672;  
       fmcbb[3]= 0.0362127;   fmccc[3]= 0.0729588;   fmcc[3]= 0.157833;   fmcll[3]= 0.732995;  
      Fbb[4]= 1.19366;   Fcc[4]= 1.19366;   Fc[4]= 0.929929;   Fll[4]= 0.960923; 
         cafac[4]= 0.888562;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22242.5;   wouttag[4]= 3569.19;  
       fbb[4]= 0.0827553;   fcc[4]= 0.140318;   fc[4]= 0.132517;   fll[4]= 0.644409;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18386;   Fcc[5]= 1.18386;   Fc[5]= 0.922294;   Fll[5]= 0.953033; 
         cafac[5]= 0.824143;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.68;   wouttag[5]= 1216.94;  
       fbb[5]= 0.104276;   fcc[5]= 0.156627;   fc[5]= 0.117107;   fll[5]= 0.621991;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.21078;   Fcc[6]= 1.21078;   Fc[6]= 0.943269;   Fll[6]= 0.974707; 
         cafac[6]= 0.848038;  
         winpretag[6]= 146990 ;   wintag[6]= 17354.7 ;   woutpretag[6]= 124653;   wouttag[6]= 15540.2;  
       fbb[6]= 0.053618;   fcc[6]= 0.1009;   fc[6]= 0.145052;   fll[6]= 0.70043;  
       fmcbb[6]= 0.0442838;   fmccc[6]= 0.0833349;   fmcc[6]= 0.153776;   fmcll[6]= 0.718605;  
      Fbb[7]= 1.19153;   Fcc[7]= 1.19153;   Fc[7]= 0.928271;   Fll[7]= 0.959209; 
         cafac[7]= 0.874032;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27901.2;   wouttag[7]= 4803.54;  
       fbb[7]= 0.0874306;   fcc[7]= 0.143861;   fc[7]= 0.129169;   fll[7]= 0.639539;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==37   || sysname=="WbbWccjet4_up" ) { 
      Fbb[1]= 1.23358;   Fcc[1]= 1.23358;   Fc[1]= 0.961029;   Fll[1]= 0.993059; 
         cafac[1]= 0.9976;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.2387e+06;   wouttag[1]= 80628.3;  
       fbb[1]= 0.01471;   fcc[1]= 0.0426761;   fc[1]= 0.127449;   fll[1]= 0.815165;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.91873;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37108.9;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20558;   Fcc[3]= 1.20558;   Fc[3]= 0.93922;   Fll[3]= 0.970523; 
         cafac[3]= 0.832384;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95780.3;   wouttag[3]= 11676.5;  
       fbb[3]= 0.05821;   fcc[3]= 0.117277;   fc[3]= 0.142185;   fll[3]= 0.682328;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.18052;   Fcc[4]= 1.18052;   Fc[4]= 0.919696;   Fll[4]= 0.950348; 
         cafac[4]= 0.879134;  
         winpretag[4]= 25032 ;   wintag[4]= 4119.6 ;   woutpretag[4]= 22006.5;   wouttag[4]= 3896.05;  
       fbb[4]= 0.102306;   fcc[4]= 0.173467;   fc[4]= 0.123529;   fll[4]= 0.600699;  
       fmcbb[4]= 0.0866614;   fmccc[4]= 0.146941;   fmcc[4]= 0.134315;   fmcll[4]= 0.632083;  
      Fbb[5]= 1.18386;   Fcc[5]= 1.18386;   Fc[5]= 0.922294;   Fll[5]= 0.953033; 
         cafac[5]= 0.824143;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.68;   wouttag[5]= 1216.94;  
       fbb[5]= 0.104276;   fcc[5]= 0.156627;   fc[5]= 0.117107;   fll[5]= 0.621991;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.20021;   Fcc[6]= 1.20021;   Fc[6]= 0.935035;   Fll[6]= 0.966199; 
         cafac[6]= 0.840252;  
         winpretag[6]= 146990 ;   wintag[6]= 18765.5 ;   woutpretag[6]= 123509;   wouttag[6]= 16778.1;  
       fbb[6]= 0.0680338;   fcc[6]= 0.128876;   fc[6]= 0.137763;   fll[6]= 0.665327;  
       fmcbb[6]= 0.0566849;   fmccc[6]= 0.107378;   fmcc[6]= 0.147334;   fmcll[6]= 0.688603;  
      Fbb[7]= 1.18124;   Fcc[7]= 1.18124;   Fc[7]= 0.920255;   Fll[7]= 0.950926; 
         cafac[7]= 0.866904;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5506.85 ;   woutpretag[7]= 27673.6;   wouttag[7]= 5121.44;  
       fbb[7]= 0.10273;   fcc[7]= 0.16984;   fc[7]= 0.122145;   fll[7]= 0.605284;  
       fmcbb[7]= 0.0869679;   fmccc[7]= 0.143781;   fmcc[7]= 0.13273;   fmcll[7]= 0.636521;  
 }  
 else if ( idsys==38   || sysname=="WbbWccjet4_down" ) { 
      Fbb[1]= 1.23358;   Fcc[1]= 1.23358;   Fc[1]= 0.961029;   Fll[1]= 0.993059; 
         cafac[1]= 0.9976;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.2387e+06;   wouttag[1]= 80628.3;  
       fbb[1]= 0.01471;   fcc[1]= 0.0426761;   fc[1]= 0.127449;   fll[1]= 0.815165;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.91873;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37108.9;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20558;   Fcc[3]= 1.20558;   Fc[3]= 0.93922;   Fll[3]= 0.970523; 
         cafac[3]= 0.832384;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95780.3;   wouttag[3]= 11676.5;  
       fbb[3]= 0.05821;   fcc[3]= 0.117277;   fc[3]= 0.142185;   fll[3]= 0.682328;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.20709;   Fcc[4]= 1.20709;   Fc[4]= 0.940393;   Fll[4]= 0.971735; 
         cafac[4]= 0.898413;  
         winpretag[4]= 25032 ;   wintag[4]= 3386.43 ;   woutpretag[4]= 22489;   wouttag[4]= 3227.65;  
       fbb[4]= 0.0627648;   fcc[4]= 0.106423;   fc[4]= 0.141708;   fll[4]= 0.689104;  
       fmcbb[4]= 0.0519968;   fmccc[4]= 0.0881647;   fmcc[4]= 0.15069;   fmcll[4]= 0.709148;  
      Fbb[5]= 1.18386;   Fcc[5]= 1.18386;   Fc[5]= 0.922294;   Fll[5]= 0.953033; 
         cafac[5]= 0.824143;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.68;   wouttag[5]= 1216.94;  
       fbb[5]= 0.104276;   fcc[5]= 0.156627;   fc[5]= 0.117107;   fll[5]= 0.621991;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.2048;   Fcc[6]= 1.2048;   Fc[6]= 0.938612;   Fll[6]= 0.969895; 
         cafac[6]= 0.84389;  
         winpretag[6]= 146990 ;   wintag[6]= 18032.4 ;   woutpretag[6]= 124043;   wouttag[6]= 16124.5;  
       fbb[6]= 0.0611818;   fcc[6]= 0.117309;   fc[6]= 0.140907;   fll[6]= 0.680601;  
       fmcbb[6]= 0.0507816;   fmccc[6]= 0.0973681;   fmcc[6]= 0.150123;   fmcll[6]= 0.701727;  
      Fbb[7]= 1.202;   Fcc[7]= 1.202;   Fc[7]= 0.936427;   Fll[7]= 0.967637; 
         cafac[7]= 0.881408;  
         winpretag[7]= 31922.4 ;   wintag[7]= 4773.68 ;   woutpretag[7]= 28136.6;   wouttag[7]= 4474.64;  
       fbb[7]= 0.0718622;   fcc[7]= 0.117425;   fc[7]= 0.136317;   fll[7]= 0.674396;  
       fmcbb[7]= 0.0597857;   fmccc[7]= 0.0976917;   fmcc[7]= 0.145571;   fmcll[7]= 0.696952;  
 }  
 else if ( idsys==73   || sysname=="WbbWccjet5_up" ) { 
      Fbb[1]= 1.23358;   Fcc[1]= 1.23358;   Fc[1]= 0.961029;   Fll[1]= 0.993059; 
         cafac[1]= 0.9976;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.2387e+06;   wouttag[1]= 80628.3;  
       fbb[1]= 0.01471;   fcc[1]= 0.0426761;   fc[1]= 0.127449;   fll[1]= 0.815165;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.91873;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37108.9;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20558;   Fcc[3]= 1.20558;   Fc[3]= 0.93922;   Fll[3]= 0.970523; 
         cafac[3]= 0.832384;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95780.3;   wouttag[3]= 11676.5;  
       fbb[3]= 0.05821;   fcc[3]= 0.117277;   fc[3]= 0.142185;   fll[3]= 0.682328;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19366;   Fcc[4]= 1.19366;   Fc[4]= 0.929929;   Fll[4]= 0.960923; 
         cafac[4]= 0.888562;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22242.5;   wouttag[4]= 3569.19;  
       fbb[4]= 0.0827553;   fcc[4]= 0.140318;   fc[4]= 0.132517;   fll[4]= 0.644409;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.16867;   Fcc[5]= 1.16867;   Fc[5]= 0.910464;   Fll[5]= 0.940809; 
         cafac[5]= 0.81101;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1515.59 ;   woutpretag[5]= 5588.19;   wouttag[5]= 1313.21;  
       fbb[5]= 0.128673;   fcc[5]= 0.193272;   fc[5]= 0.107435;   fll[5]= 0.57062;  
       fmcbb[5]= 0.110102;   fmccc[5]= 0.165378;   fmcc[5]= 0.118;   fmcll[5]= 0.606521;  
      Fbb[6]= 1.20176;   Fcc[6]= 1.20176;   Fc[6]= 0.936241;   Fll[6]= 0.967445; 
         cafac[6]= 0.841247;  
         winpretag[6]= 146990 ;   wintag[6]= 18527.3 ;   woutpretag[6]= 123655;   wouttag[6]= 16566.3;  
       fbb[6]= 0.0658149;   fcc[6]= 0.124891;   fc[6]= 0.138852;   fll[6]= 0.670442;  
       fmcbb[6]= 0.0547655;   fmccc[6]= 0.103923;   fmcc[6]= 0.148308;   fmcll[6]= 0.693003;  
      Fbb[7]= 1.18817;   Fcc[7]= 1.18817;   Fc[7]= 0.925658;   Fll[7]= 0.956509; 
         cafac[7]= 0.870714;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5268.61 ;   woutpretag[7]= 27795.3;   wouttag[7]= 4914.85;  
       fbb[7]= 0.0928319;   fcc[7]= 0.151939;   fc[7]= 0.127013;   fll[7]= 0.628216;  
       fmcbb[7]= 0.0781299;   fmccc[7]= 0.127876;   fmcc[7]= 0.137214;   fmcll[7]= 0.656781;  
 }  
 else if ( idsys==74   || sysname=="WbbWccjet5_down" ) { 
      Fbb[1]= 1.23358;   Fcc[1]= 1.23358;   Fc[1]= 0.961029;   Fll[1]= 0.993059; 
         cafac[1]= 0.9976;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.2387e+06;   wouttag[1]= 80628.3;  
       fbb[1]= 0.01471;   fcc[1]= 0.0426761;   fc[1]= 0.127449;   fll[1]= 0.815165;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.91873;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37108.9;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20558;   Fcc[3]= 1.20558;   Fc[3]= 0.93922;   Fll[3]= 0.970523; 
         cafac[3]= 0.832384;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95780.3;   wouttag[3]= 11676.5;  
       fbb[3]= 0.05821;   fcc[3]= 0.117277;   fc[3]= 0.142185;   fll[3]= 0.682328;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19366;   Fcc[4]= 1.19366;   Fc[4]= 0.929929;   Fll[4]= 0.960923; 
         cafac[4]= 0.888562;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22242.5;   wouttag[4]= 3569.19;  
       fbb[4]= 0.0827553;   fcc[4]= 0.140318;   fc[4]= 0.132517;   fll[4]= 0.644409;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.19944;   Fcc[5]= 1.19944;   Fc[5]= 0.934436;   Fll[5]= 0.965579; 
         cafac[5]= 0.838071;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1258.9 ;   woutpretag[5]= 5774.65;   wouttag[5]= 1114.84;  
       fbb[5]= 0.0792364;   fcc[5]= 0.119017;   fc[5]= 0.127033;   fll[5]= 0.674714;  
       fmcbb[5]= 0.0660611;   fmccc[5]= 0.0992266;   fmcc[5]= 0.135946;   fmcll[5]= 0.698766;  
      Fbb[6]= 1.20325;   Fcc[6]= 1.20325;   Fc[6]= 0.9374;   Fll[6]= 0.968642; 
         cafac[6]= 0.842883;  
         winpretag[6]= 146990 ;   wintag[6]= 18270.6 ;   woutpretag[6]= 123895;   wouttag[6]= 16338.6;  
       fbb[6]= 0.0634123;   fcc[6]= 0.121314;   fc[6]= 0.139813;   fll[6]= 0.675461;  
       fmcbb[6]= 0.052701;   fmccc[6]= 0.100822;   fmcc[6]= 0.149149;   fmcll[6]= 0.697327;  
      Fbb[7]= 1.1949;   Fcc[7]= 1.1949;   Fc[7]= 0.930898;   Fll[7]= 0.961924; 
         cafac[7]= 0.877395;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5011.92 ;   woutpretag[7]= 28008.5;   wouttag[7]= 4690.74;  
       fbb[7]= 0.0819986;   fcc[7]= 0.135738;   fc[7]= 0.131338;   fll[7]= 0.650926;  
       fmcbb[7]= 0.0686237;   fmccc[7]= 0.113597;   fmcc[7]= 0.141087;   fmcll[7]= 0.676692;  
 }  
 else if ( idsys==41   || sysname=="Wcjet1_up" ) { 
      Fbb[1]= 1.23542;   Fcc[1]= 1.23542;   Fc[1]= 0.962463;   Fll[1]= 0.99454; 
         cafac[1]= 1.0424;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 90158.2 ;   woutpretag[1]= 2.33922e+06;   wouttag[1]= 94377.9;  
       fbb[1]= 0.0141688;   fcc[1]= 0.0411061;   fc[1]= 0.159548;   fll[1]= 0.785177;  
       fmcbb[1]= 0.0114688;   fmccc[1]= 0.033273;   fmcc[1]= 0.165771;   fmcll[1]= 0.789487;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.91873;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37108.9;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20558;   Fcc[3]= 1.20558;   Fc[3]= 0.93922;   Fll[3]= 0.970523; 
         cafac[3]= 0.832384;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95780.3;   wouttag[3]= 11676.5;  
       fbb[3]= 0.05821;   fcc[3]= 0.117277;   fc[3]= 0.142185;   fll[3]= 0.682328;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19366;   Fcc[4]= 1.19366;   Fc[4]= 0.929929;   Fll[4]= 0.960923; 
         cafac[4]= 0.888562;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22242.5;   wouttag[4]= 3569.19;  
       fbb[4]= 0.0827553;   fcc[4]= 0.140318;   fc[4]= 0.132517;   fll[4]= 0.644409;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18386;   Fcc[5]= 1.18386;   Fc[5]= 0.922294;   Fll[5]= 0.953033; 
         cafac[5]= 0.824143;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.68;   wouttag[5]= 1216.94;  
       fbb[5]= 0.104276;   fcc[5]= 0.156627;   fc[5]= 0.117107;   fll[5]= 0.621991;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.2025;   Fcc[6]= 1.2025;   Fc[6]= 0.93682;   Fll[6]= 0.968043; 
         cafac[6]= 0.842063;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123775;   wouttag[6]= 16452.6;  
       fbb[6]= 0.0646144;   fcc[6]= 0.123104;   fc[6]= 0.139332;   fll[6]= 0.67295;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.19153;   Fcc[7]= 1.19153;   Fc[7]= 0.928271;   Fll[7]= 0.959209; 
         cafac[7]= 0.874032;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27901.2;   wouttag[7]= 4803.54;  
       fbb[7]= 0.0874306;   fcc[7]= 0.143861;   fc[7]= 0.129169;   fll[7]= 0.639539;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==42   || sysname=="Wcjet1_down" ) { 
      Fbb[1]= 1.23174;   Fcc[1]= 1.23174;   Fc[1]= 0.9596;   Fll[1]= 0.991582; 
         cafac[1]= 0.956613;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 69905.5 ;   woutpretag[1]= 2.14672e+06;   wouttag[1]= 68047.7;  
       fbb[1]= 0.0152495;   fcc[1]= 0.0442414;   fc[1]= 0.0954442;   fll[1]= 0.845065;  
       fmcbb[1]= 0.0123804;   fmccc[1]= 0.0359177;   fmcc[1]= 0.0994626;   fmcll[1]= 0.852239;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.91873;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37108.9;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20558;   Fcc[3]= 1.20558;   Fc[3]= 0.93922;   Fll[3]= 0.970523; 
         cafac[3]= 0.832384;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95780.3;   wouttag[3]= 11676.5;  
       fbb[3]= 0.05821;   fcc[3]= 0.117277;   fc[3]= 0.142185;   fll[3]= 0.682328;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19366;   Fcc[4]= 1.19366;   Fc[4]= 0.929929;   Fll[4]= 0.960923; 
         cafac[4]= 0.888562;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22242.5;   wouttag[4]= 3569.19;  
       fbb[4]= 0.0827553;   fcc[4]= 0.140318;   fc[4]= 0.132517;   fll[4]= 0.644409;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18386;   Fcc[5]= 1.18386;   Fc[5]= 0.922294;   Fll[5]= 0.953033; 
         cafac[5]= 0.824143;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.68;   wouttag[5]= 1216.94;  
       fbb[5]= 0.104276;   fcc[5]= 0.156627;   fc[5]= 0.117107;   fll[5]= 0.621991;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.2025;   Fcc[6]= 1.2025;   Fc[6]= 0.93682;   Fll[6]= 0.968043; 
         cafac[6]= 0.842063;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123775;   wouttag[6]= 16452.6;  
       fbb[6]= 0.0646144;   fcc[6]= 0.123104;   fc[6]= 0.139332;   fll[6]= 0.67295;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.19153;   Fcc[7]= 1.19153;   Fc[7]= 0.928271;   Fll[7]= 0.959209; 
         cafac[7]= 0.874032;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27901.2;   wouttag[7]= 4803.54;  
       fbb[7]= 0.0874306;   fcc[7]= 0.143861;   fc[7]= 0.129169;   fll[7]= 0.639539;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==43   || sysname=="Wcjet3_up" ) { 
      Fbb[1]= 1.23358;   Fcc[1]= 1.23358;   Fc[1]= 0.961029;   Fll[1]= 0.993059; 
         cafac[1]= 0.9976;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.2387e+06;   wouttag[1]= 80628.3;  
       fbb[1]= 0.01471;   fcc[1]= 0.0426761;   fc[1]= 0.127449;   fll[1]= 0.815165;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.91873;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37108.9;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20886;   Fcc[3]= 1.20886;   Fc[3]= 0.941773;   Fll[3]= 0.973161; 
         cafac[3]= 0.877835;  
         winpretag[3]= 115067 ;   wintag[3]= 13915.1 ;   woutpretag[3]= 101010;   wouttag[3]= 12828.5;  
       fbb[3]= 0.0557651;   fcc[3]= 0.112351;   fc[3]= 0.178214;   fll[3]= 0.65367;  
       fmcbb[3]= 0.0461303;   fmccc[3]= 0.09294;   fmcc[3]= 0.189232;   fmcll[3]= 0.671697;  
      Fbb[4]= 1.19366;   Fcc[4]= 1.19366;   Fc[4]= 0.929929;   Fll[4]= 0.960923; 
         cafac[4]= 0.888562;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22242.5;   wouttag[4]= 3569.19;  
       fbb[4]= 0.0827553;   fcc[4]= 0.140318;   fc[4]= 0.132517;   fll[4]= 0.644409;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18386;   Fcc[5]= 1.18386;   Fc[5]= 0.922294;   Fll[5]= 0.953033; 
         cafac[5]= 0.824143;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.68;   wouttag[5]= 1216.94;  
       fbb[5]= 0.104276;   fcc[5]= 0.156627;   fc[5]= 0.117107;   fll[5]= 0.621991;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.20505;   Fcc[6]= 1.20505;   Fc[6]= 0.938807;   Fll[6]= 0.970096; 
         cafac[6]= 0.876518;  
         winpretag[6]= 146990 ;   wintag[6]= 19055.4 ;   woutpretag[6]= 128839;   wouttag[6]= 17640.8;  
       fbb[6]= 0.06272;   fcc[6]= 0.119272;   fc[6]= 0.167442;   fll[6]= 0.650566;  
       fmcbb[6]= 0.0520475;   fmccc[6]= 0.0989767;   fmcc[6]= 0.178356;   fmcll[6]= 0.67062;  
      Fbb[7]= 1.19153;   Fcc[7]= 1.19153;   Fc[7]= 0.928271;   Fll[7]= 0.959209; 
         cafac[7]= 0.874032;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27901.2;   wouttag[7]= 4803.54;  
       fbb[7]= 0.0874306;   fcc[7]= 0.143861;   fc[7]= 0.129169;   fll[7]= 0.639539;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==44   || sysname=="Wcjet3_down" ) { 
      Fbb[1]= 1.23358;   Fcc[1]= 1.23358;   Fc[1]= 0.961029;   Fll[1]= 0.993059; 
         cafac[1]= 0.9976;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.2387e+06;   wouttag[1]= 80628.3;  
       fbb[1]= 0.01471;   fcc[1]= 0.0426761;   fc[1]= 0.127449;   fll[1]= 0.815165;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.91873;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37108.9;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20232;   Fcc[3]= 1.20232;   Fc[3]= 0.936681;   Fll[3]= 0.967899; 
         cafac[3]= 0.791619;  
         winpretag[3]= 115067 ;   wintag[3]= 12602.2 ;   woutpretag[3]= 91089.6;   wouttag[3]= 10643.4;  
       fbb[3]= 0.0606416;   fcc[3]= 0.122176;   fc[3]= 0.10635;   fll[3]= 0.710832;  
       fmcbb[3]= 0.050437;   fmccc[3]= 0.101617;   fmcc[3]= 0.113539;   fmcll[3]= 0.734407;  
      Fbb[4]= 1.19366;   Fcc[4]= 1.19366;   Fc[4]= 0.929929;   Fll[4]= 0.960923; 
         cafac[4]= 0.888562;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22242.5;   wouttag[4]= 3569.19;  
       fbb[4]= 0.0827553;   fcc[4]= 0.140318;   fc[4]= 0.132517;   fll[4]= 0.644409;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18386;   Fcc[5]= 1.18386;   Fc[5]= 0.922294;   Fll[5]= 0.953033; 
         cafac[5]= 0.824143;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.68;   wouttag[5]= 1216.94;  
       fbb[5]= 0.104276;   fcc[5]= 0.156627;   fc[5]= 0.117107;   fll[5]= 0.621991;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.19996;   Fcc[6]= 1.19996;   Fc[6]= 0.934842;   Fll[6]= 0.965999; 
         cafac[6]= 0.810345;  
         winpretag[6]= 146990 ;   wintag[6]= 17742.5 ;   woutpretag[6]= 119112;   wouttag[6]= 15358.8;  
       fbb[6]= 0.0665007;   fcc[6]= 0.126919;   fc[6]= 0.111341;   fll[6]= 0.695239;  
       fmcbb[6]= 0.0554189;   fmccc[6]= 0.105769;   fmcc[6]= 0.119101;   fmcll[6]= 0.71971;  
      Fbb[7]= 1.19153;   Fcc[7]= 1.19153;   Fc[7]= 0.928271;   Fll[7]= 0.959209; 
         cafac[7]= 0.874032;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27901.2;   wouttag[7]= 4803.54;  
       fbb[7]= 0.0874306;   fcc[7]= 0.143861;   fc[7]= 0.129169;   fll[7]= 0.639539;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==45   || sysname=="Wcjet4_up" ) { 
      Fbb[1]= 1.23358;   Fcc[1]= 1.23358;   Fc[1]= 0.961029;   Fll[1]= 0.993059; 
         cafac[1]= 0.9976;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.2387e+06;   wouttag[1]= 80628.3;  
       fbb[1]= 0.01471;   fcc[1]= 0.0426761;   fc[1]= 0.127449;   fll[1]= 0.815165;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.91873;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37108.9;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20558;   Fcc[3]= 1.20558;   Fc[3]= 0.93922;   Fll[3]= 0.970523; 
         cafac[3]= 0.832384;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95780.3;   wouttag[3]= 11676.5;  
       fbb[3]= 0.05821;   fcc[3]= 0.117277;   fc[3]= 0.142185;   fll[3]= 0.682328;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19714;   Fcc[4]= 1.19714;   Fc[4]= 0.932645;   Fll[4]= 0.963728; 
         cafac[4]= 0.935168;  
         winpretag[4]= 25032 ;   wintag[4]= 3881.82 ;   woutpretag[4]= 23409.1;   wouttag[4]= 3859.18;  
       fbb[4]= 0.0795487;   fcc[4]= 0.134881;   fc[4]= 0.16613;   fll[4]= 0.61944;  
       fmcbb[4]= 0.0664488;   fmccc[4]= 0.112669;   fmcc[4]= 0.178128;   fmcll[4]= 0.642754;  
      Fbb[5]= 1.18386;   Fcc[5]= 1.18386;   Fc[5]= 0.922294;   Fll[5]= 0.953033; 
         cafac[5]= 0.824143;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.68;   wouttag[5]= 1216.94;  
       fbb[5]= 0.104276;   fcc[5]= 0.156627;   fc[5]= 0.117107;   fll[5]= 0.621991;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.2031;   Fcc[6]= 1.2031;   Fc[6]= 0.937288;   Fll[6]= 0.968527; 
         cafac[6]= 0.849923;  
         winpretag[6]= 146990 ;   wintag[6]= 18527.8 ;   woutpretag[6]= 124930;   wouttag[6]= 16698.3;  
       fbb[6]= 0.0640565;   fcc[6]= 0.122165;   fc[6]= 0.145088;   fll[6]= 0.668691;  
       fmcbb[6]= 0.0532427;   fmccc[6]= 0.101541;   fmcc[6]= 0.154796;   fmcll[6]= 0.69042;  
      Fbb[7]= 1.19425;   Fcc[7]= 1.19425;   Fc[7]= 0.930391;   Fll[7]= 0.961399; 
         cafac[7]= 0.909052;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5269.07 ;   woutpretag[7]= 29019.1;   wouttag[7]= 5096.45;  
       fbb[7]= 0.0849329;   fcc[7]= 0.139616;   fc[7]= 0.155456;   fll[7]= 0.619995;  
       fmcbb[7]= 0.0711182;   fmccc[7]= 0.116907;   fmcc[7]= 0.167086;   fmcll[7]= 0.644889;  
 }  
 else if ( idsys==46   || sysname=="Wcjet4_down" ) { 
      Fbb[1]= 1.23358;   Fcc[1]= 1.23358;   Fc[1]= 0.961029;   Fll[1]= 0.993059; 
         cafac[1]= 0.9976;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.2387e+06;   wouttag[1]= 80628.3;  
       fbb[1]= 0.01471;   fcc[1]= 0.0426761;   fc[1]= 0.127449;   fll[1]= 0.815165;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.91873;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37108.9;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20558;   Fcc[3]= 1.20558;   Fc[3]= 0.93922;   Fll[3]= 0.970523; 
         cafac[3]= 0.832384;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95780.3;   wouttag[3]= 11676.5;  
       fbb[3]= 0.05821;   fcc[3]= 0.117277;   fc[3]= 0.142185;   fll[3]= 0.682328;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19019;   Fcc[4]= 1.19019;   Fc[4]= 0.92723;   Fll[4]= 0.958133; 
         cafac[4]= 0.846615;  
         winpretag[4]= 25032 ;   wintag[4]= 3624.21 ;   woutpretag[4]= 21192.4;   wouttag[4]= 3308.19;  
       fbb[4]= 0.0859432;   fcc[4]= 0.145723;   fc[4]= 0.0990995;   fll[4]= 0.669234;  
       fmcbb[4]= 0.0722095;   fmccc[4]= 0.122437;   fmcc[4]= 0.106877;   fmcll[4]= 0.698477;  
      Fbb[5]= 1.18386;   Fcc[5]= 1.18386;   Fc[5]= 0.922294;   Fll[5]= 0.953033; 
         cafac[5]= 0.824143;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.68;   wouttag[5]= 1216.94;  
       fbb[5]= 0.104276;   fcc[5]= 0.156627;   fc[5]= 0.117107;   fll[5]= 0.621991;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.2019;   Fcc[6]= 1.2019;   Fc[6]= 0.936353;   Fll[6]= 0.96756; 
         cafac[6]= 0.834355;  
         winpretag[6]= 146990 ;   wintag[6]= 18270.1 ;   woutpretag[6]= 122642;   wouttag[6]= 16211.6;  
       fbb[6]= 0.0651717;   fcc[6]= 0.124042;   fc[6]= 0.133582;   fll[6]= 0.677205;  
       fmcbb[6]= 0.0542237;   fmccc[6]= 0.103205;   fmcc[6]= 0.142662;   fmcll[6]= 0.69991;  
      Fbb[7]= 1.18882;   Fcc[7]= 1.18882;   Fc[7]= 0.92616;   Fll[7]= 0.957028; 
         cafac[7]= 0.841753;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5011.46 ;   woutpretag[7]= 26870.7;   wouttag[7]= 4533.55;  
       fbb[7]= 0.0899169;   fcc[7]= 0.148087;   fc[7]= 0.103003;   fll[7]= 0.658994;  
       fmcbb[7]= 0.0756354;   fmccc[7]= 0.124566;   fmcc[7]= 0.111215;   fmcll[7]= 0.688584;  
 }  
 else if ( idsys==75   || sysname=="Wcjet5_up" ) { 
      Fbb[1]= 1.23358;   Fcc[1]= 1.23358;   Fc[1]= 0.961029;   Fll[1]= 0.993059; 
         cafac[1]= 0.9976;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.2387e+06;   wouttag[1]= 80628.3;  
       fbb[1]= 0.01471;   fcc[1]= 0.0426761;   fc[1]= 0.127449;   fll[1]= 0.815165;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.91873;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37108.9;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20558;   Fcc[3]= 1.20558;   Fc[3]= 0.93922;   Fll[3]= 0.970523; 
         cafac[3]= 0.832384;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95780.3;   wouttag[3]= 11676.5;  
       fbb[3]= 0.05821;   fcc[3]= 0.117277;   fc[3]= 0.142185;   fll[3]= 0.682328;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19366;   Fcc[4]= 1.19366;   Fc[4]= 0.929929;   Fll[4]= 0.960923; 
         cafac[4]= 0.888562;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22242.5;   wouttag[4]= 3569.19;  
       fbb[4]= 0.0827553;   fcc[4]= 0.140318;   fc[4]= 0.132517;   fll[4]= 0.644409;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18721;   Fcc[5]= 1.18721;   Fc[5]= 0.924908;   Fll[5]= 0.955733; 
         cafac[5]= 0.865019;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1415.36 ;   woutpretag[5]= 5960.33;   wouttag[5]= 1297.19;  
       fbb[5]= 0.100769;   fcc[5]= 0.15136;   fc[5]= 0.146798;   fll[5]= 0.601073;  
       fmcbb[5]= 0.0848788;   fmccc[5]= 0.127492;   fmcc[5]= 0.158716;   fmcll[5]= 0.628913;  
      Fbb[6]= 1.20266;   Fcc[6]= 1.20266;   Fc[6]= 0.936946;   Fll[6]= 0.968173; 
         cafac[6]= 0.844231;  
         winpretag[6]= 146990 ;   wintag[6]= 18427.1 ;   woutpretag[6]= 124093;   wouttag[6]= 16513.3;  
       fbb[6]= 0.0644425;   fcc[6]= 0.122849;   fc[6]= 0.140745;   fll[6]= 0.671963;  
       fmcbb[6]= 0.0535831;   fmccc[6]= 0.102147;   fmcc[6]= 0.150217;   fmcll[6]= 0.694053;  
      Fbb[7]= 1.19226;   Fcc[7]= 1.19226;   Fc[7]= 0.928841;   Fll[7]= 0.959798; 
         cafac[7]= 0.883469;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5168.38 ;   woutpretag[7]= 28202.4;   wouttag[7]= 4875.07;  
       fbb[7]= 0.0866601;   fcc[7]= 0.142711;   fc[7]= 0.135613;   fll[7]= 0.635016;  
       fmcbb[7]= 0.0726855;   fmccc[7]= 0.119698;   fmcc[7]= 0.146002;   fmcll[7]= 0.661614;  
 }  
 else if ( idsys==76   || sysname=="Wcjet5_down" ) { 
      Fbb[1]= 1.23358;   Fcc[1]= 1.23358;   Fc[1]= 0.961029;   Fll[1]= 0.993059; 
         cafac[1]= 0.9976;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.2387e+06;   wouttag[1]= 80628.3;  
       fbb[1]= 0.01471;   fcc[1]= 0.0426761;   fc[1]= 0.127449;   fll[1]= 0.815165;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.91873;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473183;   wouttag[2]= 37108.9;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20558;   Fcc[3]= 1.20558;   Fc[3]= 0.93922;   Fll[3]= 0.970523; 
         cafac[3]= 0.832384;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95780.3;   wouttag[3]= 11676.5;  
       fbb[3]= 0.05821;   fcc[3]= 0.117277;   fc[3]= 0.142185;   fll[3]= 0.682328;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19366;   Fcc[4]= 1.19366;   Fc[4]= 0.929929;   Fll[4]= 0.960923; 
         cafac[4]= 0.888562;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22242.5;   wouttag[4]= 3569.19;  
       fbb[4]= 0.0827553;   fcc[4]= 0.140318;   fc[4]= 0.132517;   fll[4]= 0.644409;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18052;   Fcc[5]= 1.18052;   Fc[5]= 0.919696;   Fll[5]= 0.950348; 
         cafac[5]= 0.787156;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1359.13 ;   woutpretag[5]= 5423.82;   wouttag[5]= 1144.33;  
       fbb[5]= 0.107763;   fcc[5]= 0.161865;   fc[5]= 0.0875824;   fll[5]= 0.64279;  
       fmcbb[5]= 0.0912841;   fmccc[5]= 0.137113;   fmcc[5]= 0.0952298;   fmcll[5]= 0.676374;  
      Fbb[6]= 1.20234;   Fcc[6]= 1.20234;   Fc[6]= 0.936694;   Fll[6]= 0.967913; 
         cafac[6]= 0.839908;  
         winpretag[6]= 146990 ;   wintag[6]= 18370.8 ;   woutpretag[6]= 123458;   wouttag[6]= 16392.3;  
       fbb[6]= 0.0647862;   fcc[6]= 0.123358;   fc[6]= 0.13792;   fll[6]= 0.673936;  
       fmcbb[6]= 0.0538833;   fmccc[6]= 0.102598;   fmcc[6]= 0.147241;   fmcll[6]= 0.696278;  
      Fbb[7]= 1.1908;   Fcc[7]= 1.1908;   Fc[7]= 0.927701;   Fll[7]= 0.95862; 
         cafac[7]= 0.864807;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5112.15 ;   woutpretag[7]= 27606.7;   wouttag[7]= 4733.6;  
       fbb[7]= 0.0882001;   fcc[7]= 0.145009;   fc[7]= 0.122734;   fll[7]= 0.644057;  
       fmcbb[7]= 0.0740681;   fmccc[7]= 0.121775;   fmcc[7]= 0.132299;   fmcll[7]= 0.671858;  
 }  
 else if ( idsys==3   || sysname=="wt_up" ) { 
      Fbb[1]= 1.23342;   Fcc[1]= 1.23342;   Fc[1]= 0.956097;   Fll[1]= 0.993864; 
         cafac[1]= 0.99673;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.23674e+06;   wouttag[1]= 80349.5;  
       fbb[1]= 0.0147081;   fcc[1]= 0.0426708;   fc[1]= 0.126795;   fll[1]= 0.815827;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21852;   Fcc[2]= 1.21852;   Fc[2]= 0.944549;   Fll[2]= 0.98186; 
         cafac[2]= 0.917792;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 472700;   wouttag[2]= 37011.4;  
       fbb[2]= 0.0360698;   fcc[2]= 0.0866717;   fc[2]= 0.14427;   fll[2]= 0.732989;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20567;   Fcc[3]= 1.20567;   Fc[3]= 0.934585;   Fll[3]= 0.971503; 
         cafac[3]= 0.831533;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95682.4;   wouttag[3]= 11651.6;  
       fbb[3]= 0.0582142;   fcc[3]= 0.117286;   fc[3]= 0.141483;   fll[3]= 0.683017;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19373;   Fcc[4]= 1.19373;   Fc[4]= 0.925329;   Fll[4]= 0.96188; 
         cafac[4]= 0.887825;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22224;   wouttag[4]= 3563.12;  
       fbb[4]= 0.0827602;   fcc[4]= 0.140326;   fc[4]= 0.131862;   fll[4]= 0.645052;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18386;   Fcc[5]= 1.18386;   Fc[5]= 0.91768;   Fll[5]= 0.953929; 
         cafac[5]= 0.823333;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5673.09;   wouttag[5]= 1215.01;  
       fbb[5]= 0.104276;   fcc[5]= 0.156628;   fc[5]= 0.116521;   fll[5]= 0.622576;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.20258;   Fcc[6]= 1.20258;   Fc[6]= 0.932192;   Fll[6]= 0.969015; 
         cafac[6]= 0.841235;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123653;   wouttag[6]= 16419.6;  
       fbb[6]= 0.0646187;   fcc[6]= 0.123112;   fc[6]= 0.138644;   fll[6]= 0.673626;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.19159;   Fcc[7]= 1.19159;   Fc[7]= 0.923667;   Fll[7]= 0.960153; 
         cafac[7]= 0.873279;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27877.1;   wouttag[7]= 4795.56;  
       fbb[7]= 0.0874347;   fcc[7]= 0.143868;   fc[7]= 0.128529;   fll[7]= 0.640169;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==4   || sysname=="wt_down" ) { 
      Fbb[1]= 1.23373;   Fcc[1]= 1.23373;   Fc[1]= 0.965953;   Fll[1]= 0.992254; 
         cafac[1]= 0.99847;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.24065e+06;   wouttag[1]= 80907.2;  
       fbb[1]= 0.0147118;   fcc[1]= 0.0426814;   fc[1]= 0.128102;   fll[1]= 0.814505;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21843;   Fcc[2]= 1.21843;   Fc[2]= 0.95397;   Fll[2]= 0.979945; 
         cafac[2]= 0.919668;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473666;   wouttag[2]= 37206.4;  
       fbb[2]= 0.0360669;   fcc[2]= 0.0866648;   fc[2]= 0.145709;   fll[2]= 0.73156;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.2055;   Fcc[3]= 1.2055;   Fc[3]= 0.943846;   Fll[3]= 0.969545; 
         cafac[3]= 0.833235;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95878.3;   wouttag[3]= 11701.5;  
       fbb[3]= 0.0582057;   fcc[3]= 0.117269;   fc[3]= 0.142885;   fll[3]= 0.681641;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19359;   Fcc[4]= 1.19359;   Fc[4]= 0.934521;   Fll[4]= 0.959967; 
         cafac[4]= 0.889298;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22260.9;   wouttag[4]= 3575.27;  
       fbb[4]= 0.0827503;   fcc[4]= 0.14031;   fc[4]= 0.133172;   fll[4]= 0.643768;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18385;   Fcc[5]= 1.18385;   Fc[5]= 0.926901;   Fll[5]= 0.952139; 
         cafac[5]= 0.824953;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5684.26;   wouttag[5]= 1218.88;  
       fbb[5]= 0.104275;   fcc[5]= 0.156626;   fc[5]= 0.117691;   fll[5]= 0.621407;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.20242;   Fcc[6]= 1.20242;   Fc[6]= 0.941439;   Fll[6]= 0.967073; 
         cafac[6]= 0.842892;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123897;   wouttag[6]= 16485.7;  
       fbb[6]= 0.06461;   fcc[6]= 0.123095;   fc[6]= 0.140019;   fll[6]= 0.672276;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.19147;   Fcc[7]= 1.19147;   Fc[7]= 0.932866;   Fll[7]= 0.958266; 
         cafac[7]= 0.874785;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27925.2;   wouttag[7]= 4811.51;  
       fbb[7]= 0.0874264;   fcc[7]= 0.143854;   fc[7]= 0.129809;   fll[7]= 0.638911;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==5   || sysname=="st_up" ) { 
      Fbb[1]= 1.21309;   Fcc[1]= 1.21309;   Fc[1]= 0.96189;   Fll[1]= 0.994081; 
         cafac[1]= 0.997664;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.23884e+06;   wouttag[1]= 80456.9;  
       fbb[1]= 0.0144656;   fcc[1]= 0.0419672;   fc[1]= 0.127563;   fll[1]= 0.816004;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.19962;   Fcc[2]= 1.19962;   Fc[2]= 0.951213;   Fll[2]= 0.983047; 
         cafac[2]= 0.918499;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473064;   wouttag[2]= 36937.2;  
       fbb[2]= 0.0355102;   fcc[2]= 0.0853271;   fc[2]= 0.145288;   fll[2]= 0.733875;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.18803;   Fcc[3]= 1.18803;   Fc[3]= 0.942026;   Fll[3]= 0.973552; 
         cafac[3]= 0.831452;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95673.1;   wouttag[3]= 11602;  
       fbb[3]= 0.0573626;   fcc[3]= 0.11557;   fc[3]= 0.142609;   fll[3]= 0.684458;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.17727;   Fcc[4]= 1.17727;   Fc[4]= 0.933491;   Fll[4]= 0.964732; 
         cafac[4]= 0.887043;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22204.4;   wouttag[4]= 3541.63;  
       fbb[4]= 0.0816192;   fcc[4]= 0.138392;   fc[4]= 0.133025;   fll[4]= 0.646964;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.16839;   Fcc[5]= 1.16839;   Fc[5]= 0.926445;   Fll[5]= 0.95745; 
         cafac[5]= 0.822317;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5666.1;   wouttag[5]= 1207.68;  
       fbb[5]= 0.102913;   fcc[5]= 0.15458;   fc[5]= 0.117634;   fll[5]= 0.624873;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.18526;   Fcc[6]= 1.18526;   Fc[6]= 0.939822;   Fll[6]= 0.971274; 
         cafac[6]= 0.840994;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123618;   wouttag[6]= 16341.2;  
       fbb[6]= 0.0636876;   fcc[6]= 0.121338;   fc[6]= 0.139778;   fll[6]= 0.675196;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.17534;   Fcc[7]= 1.17534;   Fc[7]= 0.931961;   Fll[7]= 0.963151; 
         cafac[7]= 0.872453;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27850.8;   wouttag[7]= 4766.61;  
       fbb[7]= 0.0862428;   fcc[7]= 0.141907;   fc[7]= 0.129683;   fll[7]= 0.642168;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==6   || sysname=="st_down" ) { 
      Fbb[1]= 1.25219;   Fcc[1]= 1.25219;   Fc[1]= 0.960246;   Fll[1]= 0.99213; 
         cafac[1]= 0.997541;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.23856e+06;   wouttag[1]= 80784;  
       fbb[1]= 0.014932;   fcc[1]= 0.0433202;   fc[1]= 0.127345;   fll[1]= 0.814403;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.23557;   Fcc[2]= 1.23557;   Fc[2]= 0.947496;   Fll[2]= 0.978957; 
         cafac[2]= 0.918939;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473291;   wouttag[2]= 37264.6;  
       fbb[2]= 0.0365743;   fcc[2]= 0.087884;   fc[2]= 0.14472;   fll[2]= 0.730822;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.22146;   Fcc[3]= 1.22146;   Fc[3]= 0.93668;   Fll[3]= 0.967782; 
         cafac[3]= 0.833229;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95877.6;   wouttag[3]= 11744.2;  
       fbb[3]= 0.0589767;   fcc[3]= 0.118822;   fc[3]= 0.1418;   fll[3]= 0.680401;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.20846;   Fcc[4]= 1.20846;   Fc[4]= 0.92671;   Fll[4]= 0.957481; 
         cafac[4]= 0.88994;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22277;   wouttag[4]= 3594.18;  
       fbb[4]= 0.0837817;   fcc[4]= 0.142058;   fc[4]= 0.132059;   fll[4]= 0.642101;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.19782;   Fcc[5]= 1.19782;   Fc[5]= 0.918548;   Fll[5]= 0.949047; 
         cafac[5]= 0.825799;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5690.09;   wouttag[5]= 1225.34;  
       fbb[5]= 0.105506;   fcc[5]= 0.158474;   fc[5]= 0.116631;   fll[5]= 0.61939;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.21811;   Fcc[6]= 1.21811;   Fc[6]= 0.934105;   Fll[6]= 0.96512; 
         cafac[6]= 0.843033;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123917;   wouttag[6]= 16553.7;  
       fbb[6]= 0.0654527;   fcc[6]= 0.124701;   fc[6]= 0.138928;   fll[6]= 0.670918;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.20615;   Fcc[7]= 1.20615;   Fc[7]= 0.924936;   Fll[7]= 0.955648; 
         cafac[7]= 0.875465;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27946.9;   wouttag[7]= 4837.02;  
       fbb[7]= 0.0885034;   fcc[7]= 0.145626;   fc[7]= 0.128705;   fll[7]= 0.637165;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==7   || sysname=="z_up" ) { 
      Fbb[1]= 1.23278;   Fcc[1]= 1.23278;   Fc[1]= 0.904485;   Fll[1]= 1.00224; 
         cafac[1]= 0.981298;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.20211e+06;   wouttag[1]= 76966.8;  
       fbb[1]= 0.0147004;   fcc[1]= 0.0426484;   fc[1]= 0.11995;   fll[1]= 0.822701;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21993;   Fcc[2]= 1.21993;   Fc[2]= 0.89506;   Fll[2]= 0.991795; 
         cafac[2]= 0.903243;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 465207;   wouttag[2]= 35817.8;  
       fbb[2]= 0.0361114;   fcc[2]= 0.0867718;   fc[2]= 0.136711;   fll[2]= 0.740406;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20742;   Fcc[3]= 1.20742;   Fc[3]= 0.885883;   Fll[3]= 0.981627; 
         cafac[3]= 0.81847;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 94179.3;   wouttag[3]= 11336.8;  
       fbb[3]= 0.0582988;   fcc[3]= 0.117456;   fc[3]= 0.13411;   fll[3]= 0.690135;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19526;   Fcc[4]= 1.19526;   Fc[4]= 0.876955;   Fll[4]= 0.971734; 
         cafac[4]= 0.877026;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 21953.7;   wouttag[4]= 3488.5;  
       fbb[4]= 0.082866;   fcc[4]= 0.140506;   fc[4]= 0.124968;   fll[4]= 0.65166;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18464;   Fcc[5]= 1.18464;   Fc[5]= 0.869167;   Fll[5]= 0.963104; 
         cafac[5]= 0.818331;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5638.63;   wouttag[5]= 1200.24;  
       fbb[5]= 0.104345;   fcc[5]= 0.156731;   fc[5]= 0.110361;   fll[5]= 0.628564;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.20425;   Fcc[6]= 1.20425;   Fc[6]= 0.883554;   Fll[6]= 0.979047; 
         cafac[6]= 0.829013;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 121857;   wouttag[6]= 16011;  
       fbb[6]= 0.0647083;   fcc[6]= 0.123283;   fc[6]= 0.13141;   fll[6]= 0.680599;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.19295;   Fcc[7]= 1.19295;   Fc[7]= 0.875262;   Fll[7]= 0.969858; 
         cafac[7]= 0.863784;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27574;   wouttag[7]= 4704.99;  
       fbb[7]= 0.0875347;   fcc[7]= 0.144032;   fc[7]= 0.121793;   fll[7]= 0.64664;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==8   || sysname=="z_down" ) { 
      Fbb[1]= 1.23435;   Fcc[1]= 1.23435;   Fc[1]= 1.0159;   Fll[1]= 0.98415; 
         cafac[1]= 1.01393;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.27534e+06;   wouttag[1]= 84298.1;  
       fbb[1]= 0.0147192;   fcc[1]= 0.0427029;   fc[1]= 0.134725;   fll[1]= 0.807853;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21707;   Fcc[2]= 1.21707;   Fc[2]= 1.00167;   Fll[2]= 0.970369; 
         cafac[2]= 0.934216;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 481159;   wouttag[2]= 38400.1;  
       fbb[2]= 0.0360267;   fcc[2]= 0.0865681;   fc[2]= 0.152995;   fll[2]= 0.724411;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.2038;   Fcc[3]= 1.2038;   Fc[3]= 0.990757;   Fll[3]= 0.959794; 
         cafac[3]= 0.846289;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 97380.4;   wouttag[3]= 12016.1;  
       fbb[3]= 0.058124;   fcc[3]= 0.117104;   fc[3]= 0.149987;   fll[3]= 0.674785;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19211;   Fcc[4]= 1.19211;   Fc[4]= 0.981134;   Fll[4]= 0.950472; 
         cafac[4]= 0.900018;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22529.2;   wouttag[4]= 3649.28;  
       fbb[4]= 0.0826481;   fcc[4]= 0.140136;   fc[4]= 0.139814;   fll[4]= 0.637401;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.1831;   Fcc[5]= 1.1831;   Fc[5]= 0.973715;   Fll[5]= 0.943285; 
         cafac[5]= 0.829792;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5717.6;   wouttag[5]= 1233.25;  
       fbb[5]= 0.104209;   fcc[5]= 0.156526;   fc[5]= 0.123636;   fll[5]= 0.615629;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.20081;   Fcc[6]= 1.20081;   Fc[6]= 0.988295;   Fll[6]= 0.95741; 
         cafac[6]= 0.855083;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 125688;   wouttag[6]= 16893;  
       fbb[6]= 0.0645236;   fcc[6]= 0.122931;   fc[6]= 0.146988;   fll[6]= 0.665558;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.19016;   Fcc[7]= 1.19016;   Fc[7]= 0.979523;   Fll[7]= 0.948912; 
         cafac[7]= 0.884182;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 28225.2;   wouttag[7]= 4901.11;  
       fbb[7]= 0.0873298;   fcc[7]= 0.143695;   fc[7]= 0.136301;   fll[7]= 0.632674;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==9   || sysname=="di_up" ) { 
      Fbb[1]= 1.23207;   Fcc[1]= 1.23207;   Fc[1]= 0.960104;   Fll[1]= 0.993294; 
         cafac[1]= 0.997399;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.23825e+06;   wouttag[1]= 80557.6;  
       fbb[1]= 0.014692;   fcc[1]= 0.0426239;   fc[1]= 0.127326;   fll[1]= 0.815358;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21713;   Fcc[2]= 1.21713;   Fc[2]= 0.948461;   Fll[2]= 0.981248; 
         cafac[2]= 0.918359;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 472992;   wouttag[2]= 37070.3;  
       fbb[2]= 0.0360284;   fcc[2]= 0.0865723;   fc[2]= 0.144867;   fll[2]= 0.732532;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20434;   Fcc[3]= 1.20434;   Fc[3]= 0.938495;   Fll[3]= 0.970937; 
         cafac[3]= 0.831979;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95733.8;   wouttag[3]= 11663.8;  
       fbb[3]= 0.0581498;   fcc[3]= 0.117156;   fc[3]= 0.142075;   fll[3]= 0.682619;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19249;   Fcc[4]= 1.19249;   Fc[4]= 0.929265;   Fll[4]= 0.961388; 
         cafac[4]= 0.888185;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22233;   wouttag[4]= 3565.51;  
       fbb[4]= 0.0826745;   fcc[4]= 0.140181;   fc[4]= 0.132423;   fll[4]= 0.644722;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18275;   Fcc[5]= 1.18275;   Fc[5]= 0.921669;   Fll[5]= 0.95353; 
         cafac[5]= 0.823834;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5676.55;   wouttag[5]= 1215.86;  
       fbb[5]= 0.104178;   fcc[5]= 0.15648;   fc[5]= 0.117027;   fll[5]= 0.622315;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.20128;   Fcc[6]= 1.20128;   Fc[6]= 0.93611;   Fll[6]= 0.96847; 
         cafac[6]= 0.84167;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123717;   wouttag[6]= 16435;  
       fbb[6]= 0.0645485;   fcc[6]= 0.122978;   fc[6]= 0.139226;   fll[6]= 0.673247;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.19038;   Fcc[7]= 1.19038;   Fc[7]= 0.927615;   Fll[7]= 0.959681; 
         cafac[7]= 0.873671;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27889.6;   wouttag[7]= 4798.75;  
       fbb[7]= 0.0873459;   fcc[7]= 0.143722;   fc[7]= 0.129078;   fll[7]= 0.639854;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==10   || sysname=="di_down" ) { 
      Fbb[1]= 1.23509;   Fcc[1]= 1.23509;   Fc[1]= 0.961953;   Fll[1]= 0.992824; 
         cafac[1]= 0.9978;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.23915e+06;   wouttag[1]= 80699;  
       fbb[1]= 0.0147279;   fcc[1]= 0.0427283;   fc[1]= 0.127571;   fll[1]= 0.814973;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21982;   Fcc[2]= 1.21982;   Fc[2]= 0.950066;   Fll[2]= 0.980555; 
         cafac[2]= 0.919101;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473374;   wouttag[2]= 37147.5;  
       fbb[2]= 0.0361082;   fcc[2]= 0.0867641;   fc[2]= 0.145112;   fll[2]= 0.732015;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20683;   Fcc[3]= 1.20683;   Fc[3]= 0.939945;   Fll[3]= 0.970109; 
         cafac[3]= 0.832789;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95826.9;   wouttag[3]= 11689.3;  
       fbb[3]= 0.05827;   fcc[3]= 0.117398;   fc[3]= 0.142294;   fll[3]= 0.682037;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19482;   Fcc[4]= 1.19482;   Fc[4]= 0.930593;   Fll[4]= 0.960457; 
         cafac[4]= 0.888939;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22251.9;   wouttag[4]= 3572.88;  
       fbb[4]= 0.0828359;   fcc[4]= 0.140455;   fc[4]= 0.132612;   fll[4]= 0.644097;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18497;   Fcc[5]= 1.18497;   Fc[5]= 0.922919;   Fll[5]= 0.952537; 
         cafac[5]= 0.824451;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5680.8;   wouttag[5]= 1218.02;  
       fbb[5]= 0.104374;   fcc[5]= 0.156774;   fc[5]= 0.117186;   fll[5]= 0.621667;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.20373;   Fcc[6]= 1.20373;   Fc[6]= 0.93753;   Fll[6]= 0.967616; 
         cafac[6]= 0.842457;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123833;   wouttag[6]= 16470.2;  
       fbb[6]= 0.0646801;   fcc[6]= 0.123229;   fc[6]= 0.139438;   fll[6]= 0.672653;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.19268;   Fcc[7]= 1.19268;   Fc[7]= 0.928926;   Fll[7]= 0.958737; 
         cafac[7]= 0.874394;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27912.7;   wouttag[7]= 4808.33;  
       fbb[7]= 0.0875151;   fcc[7]= 0.144;   fc[7]= 0.12926;   fll[7]= 0.639224;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==11   || sysname=="qcd_up" ) { 
      Fbb[1]= 1.22751;   Fcc[1]= 1.22751;   Fc[1]= 0.802447;   Fll[1]= 1.01902; 
         cafac[1]= 0.970407;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.17767e+06;   wouttag[1]= 71891;  
       fbb[1]= 0.0146376;   fcc[1]= 0.0424662;   fc[1]= 0.106418;   fll[1]= 0.836478;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21904;   Fcc[2]= 1.21904;   Fc[2]= 0.796913;   Fll[2]= 1.012; 
         cafac[2]= 0.88918;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 457964;   wouttag[2]= 34039.8;  
       fbb[2]= 0.0360852;   fcc[2]= 0.0867087;   fc[2]= 0.12172;   fll[2]= 0.755486;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20748;   Fcc[3]= 1.20748;   Fc[3]= 0.789356;   Fll[3]= 1.0024; 
         cafac[3]= 0.806014;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 92746;   wouttag[3]= 10894;  
       fbb[3]= 0.0583017;   fcc[3]= 0.117462;   fc[3]= 0.119497;   fll[3]= 0.704739;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.1951;   Fcc[4]= 1.1951;   Fc[4]= 0.781257;   Fll[4]= 0.992114; 
         cafac[4]= 0.861277;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 21559.5;   wouttag[4]= 3360.6;  
       fbb[4]= 0.0828549;   fcc[4]= 0.140487;   fc[4]= 0.111331;   fll[4]= 0.665327;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18318;   Fcc[5]= 1.18318;   Fc[5]= 0.773465;   Fll[5]= 0.982219; 
         cafac[5]= 0.800141;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5513.3;   wouttag[5]= 1157.95;  
       fbb[5]= 0.104216;   fcc[5]= 0.156537;   fc[5]= 0.0982092;   fll[5]= 0.641038;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.2042;   Fcc[6]= 1.2042;   Fc[6]= 0.787208;   Fll[6]= 0.999671; 
         cafac[6]= 0.81564;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 119891;   wouttag[6]= 15402.2;  
       fbb[6]= 0.0647055;   fcc[6]= 0.123277;   fc[6]= 0.11708;   fll[6]= 0.694937;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.1925;   Fcc[7]= 1.1925;   Fc[7]= 0.779562;   Fll[7]= 0.989961; 
         cafac[7]= 0.847493;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27054;   wouttag[7]= 4535.78;  
       fbb[7]= 0.087502;   fcc[7]= 0.143979;   fc[7]= 0.108476;   fll[7]= 0.660043;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==12   || sysname=="qcd_down" ) { 
      Fbb[1]= 1.23933;   Fcc[1]= 1.23933;   Fc[1]= 1.1113;   Fll[1]= 0.968456; 
         cafac[1]= 1.02481;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.29976e+06;   wouttag[1]= 89371.6;  
       fbb[1]= 0.0147785;   fcc[1]= 0.042875;   fc[1]= 0.147376;   fll[1]= 0.79497;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21794;   Fcc[2]= 1.21794;   Fc[2]= 1.09212;   Fll[2]= 0.951745; 
         cafac[2]= 0.948279;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 488402;   wouttag[2]= 40178.1;  
       fbb[2]= 0.0360525;   fcc[2]= 0.0866302;   fc[2]= 0.16681;   fll[2]= 0.710508;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.2038;   Fcc[3]= 1.2038;   Fc[3]= 1.07944;   Fll[3]= 0.940698; 
         cafac[3]= 0.858669;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 98804.9;   wouttag[3]= 12456.6;  
       fbb[3]= 0.058124;   fcc[3]= 0.117104;   fc[3]= 0.163412;   fll[3]= 0.661359;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19231;   Fcc[4]= 1.19231;   Fc[4]= 1.06914;   Fll[4]= 0.931717; 
         cafac[4]= 0.915725;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22922.4;   wouttag[4]= 3776.85;  
       fbb[4]= 0.0826619;   fcc[4]= 0.14016;   fc[4]= 0.152355;   fll[4]= 0.624824;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.1845;   Fcc[5]= 1.1845;   Fc[5]= 1.06213;   Fll[5]= 0.925611; 
         cafac[5]= 0.848044;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5843.37;   wouttag[5]= 1275.68;  
       fbb[5]= 0.104332;   fcc[5]= 0.156712;   fc[5]= 0.134862;   fll[5]= 0.604094;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.20091;   Fcc[6]= 1.20091;   Fc[6]= 1.07685;   Fll[6]= 0.93844; 
         cafac[6]= 0.868395;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 127645;   wouttag[6]= 17499.3;  
       fbb[6]= 0.064529;   fcc[6]= 0.122941;   fc[6]= 0.160159;   fll[6]= 0.652371;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.19062;   Fcc[7]= 1.19062;   Fc[7]= 1.06762;   Fll[7]= 0.930392; 
         cafac[7]= 0.900455;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 28744.6;   wouttag[7]= 5070.11;  
       fbb[7]= 0.0873636;   fcc[7]= 0.143751;   fc[7]= 0.14856;   fll[7]= 0.620326;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==13   || sysname=="jes_up" ) { 
      Fbb[1]= 1.36059;   Fcc[1]= 1.36059;   Fc[1]= 1.02218;   Fll[1]= 0.977825; 
         cafac[1]= 0.911356;  
         winpretag[1]= 2.45765e+06 ;   wintag[1]= 85407 ;   woutpretag[1]= 2.23979e+06;   wouttag[1]= 82146;  
       fbb[1]= 0.0149394;   fcc[1]= 0.0436222;   fc[1]= 0.13138;   fll[1]= 0.810058;  
       fmcbb[1]= 0.0109802;   fmccc[1]= 0.0320613;   fmcc[1]= 0.128529;   fmcll[1]= 0.828429;  
      Fbb[2]= 1.33268;   Fcc[2]= 1.33268;   Fc[2]= 1.00121;   Fll[2]= 0.957766; 
         cafac[2]= 0.814736;  
         winpretag[2]= 585660 ;   wintag[2]= 42791.9 ;   woutpretag[2]= 477158;   wouttag[2]= 37945.8;  
       fbb[2]= 0.0366605;   fcc[2]= 0.0900177;   fc[2]= 0.152023;   fll[2]= 0.721299;  
       fmcbb[2]= 0.0275089;   fmccc[2]= 0.0675466;   fmcc[2]= 0.151839;   fmcll[2]= 0.753106;  
      Fbb[3]= 1.31105;   Fcc[3]= 1.31105;   Fc[3]= 0.984964;   Fll[3]= 0.942225; 
         cafac[3]= 0.717681;  
         winpretag[3]= 137477 ;   wintag[3]= 15373.1 ;   woutpretag[3]= 98664.6;   wouttag[3]= 12167;  
       fbb[3]= 0.0602631;   fcc[3]= 0.121744;   fc[3]= 0.151476;   fll[3]= 0.666516;  
       fmcbb[3]= 0.0459655;   fmccc[3]= 0.0928601;   fmcc[3]= 0.153789;   fmcll[3]= 0.707386;  
      Fbb[4]= 1.29257;   Fcc[4]= 1.29257;   Fc[4]= 0.971082;   Fll[4]= 0.928946; 
         cafac[4]= 0.710619;  
         winpretag[4]= 31736.6 ;   wintag[4]= 4704.78 ;   woutpretag[4]= 22552.6;   wouttag[4]= 3724.19;  
       fbb[4]= 0.0825076;   fcc[4]= 0.148341;   fc[4]= 0.140854;   fll[4]= 0.628297;  
       fmcbb[4]= 0.063832;   fmccc[4]= 0.114764;   fmcc[4]= 0.145049;   fmcll[4]= 0.676355;  
      Fbb[5]= 1.27624;   Fcc[5]= 1.27624;   Fc[5]= 0.958813;   Fll[5]= 0.917209; 
         cafac[5]= 0.582724;  
         winpretag[5]= 9184.03 ;   wintag[5]= 1782.93 ;   woutpretag[5]= 5351.75;   wouttag[5]= 1150.34;  
       fbb[5]= 0.110002;   fcc[5]= 0.16518;   fc[5]= 0.123908;   fll[5]= 0.600911;  
       fmcbb[5]= 0.0861919;   fmccc[5]= 0.129427;   fmcc[5]= 0.12923;   fmcll[5]= 0.655151;  
      Fbb[6]= 1.3059;   Fcc[6]= 1.3059;   Fc[6]= 0.981091;   Fll[6]= 0.93852; 
         cafac[6]= 0.707544;  
         winpretag[6]= 178398 ;   wintag[6]= 21860.8 ;   woutpretag[6]= 126224;   wouttag[6]= 17123.5;  
       fbb[6]= 0.0668813;   fcc[6]= 0.128813;   fc[6]= 0.148115;   fll[6]= 0.656191;  
       fmcbb[6]= 0.0512148;   fmccc[6]= 0.0986393;   fmcc[6]= 0.15097;   fmcll[6]= 0.699176;  
      Fbb[7]= 1.28887;   Fcc[7]= 1.28887;   Fc[7]= 0.968301;   Fll[7]= 0.926285; 
         cafac[7]= 0.678999;  
         winpretag[7]= 40920.7 ;   wintag[7]= 6487.7 ;   woutpretag[7]= 27785.1;   wouttag[7]= 4901.94;  
       fbb[7]= 0.0887393;   fcc[7]= 0.152158;   fc[7]= 0.137013;   fll[7]= 0.62209;  
       fmcbb[7]= 0.0688504;   fmccc[7]= 0.118055;   fmcc[7]= 0.141498;   fmcll[7]= 0.671596;  
 }  
 else if ( idsys==14   || sysname=="jes_down" ) { 
      Fbb[1]= 1.14963;   Fcc[1]= 1.14963;   Fc[1]= 0.903332;   Fll[1]= 1.00698; 
         cafac[1]= 1.09301;  
         winpretag[1]= 2.03912e+06 ;   wintag[1]= 75210.1 ;   woutpretag[1]= 2.22878e+06;   wouttag[1]= 79602.3;  
       fbb[1]= 0.0148029;   fcc[1]= 0.0425173;   fc[1]= 0.122808;   fll[1]= 0.819872;  
       fmcbb[1]= 0.0128762;   fmccc[1]= 0.0369835;   fmcc[1]= 0.13595;   fmcll[1]= 0.81419;  
      Fbb[2]= 1.14248;   Fcc[2]= 1.14248;   Fc[2]= 0.897716;   Fll[2]= 1.00072; 
         cafac[2]= 1.0329;  
         winpretag[2]= 451981 ;   wintag[2]= 34926 ;   woutpretag[2]= 466852;   wouttag[2]= 36251;  
       fbb[2]= 0.0361699;   fcc[2]= 0.084833;   fc[2]= 0.137116;   fll[2]= 0.741881;  
       fmcbb[2]= 0.031659;   fmccc[2]= 0.0742532;   fmcc[2]= 0.152738;   fmcll[2]= 0.741349;  
      Fbb[3]= 1.13481;   Fcc[3]= 1.13481;   Fc[3]= 0.891688;   Fll[3]= 0.993998; 
         cafac[3]= 0.988754;  
         winpretag[3]= 95481.3 ;   wintag[3]= 11239.5 ;   woutpretag[3]= 94407.5;   wouttag[3]= 11355.5;  
       fbb[3]= 0.0571529;   fcc[3]= 0.11534;   fc[3]= 0.134233;   fll[3]= 0.693275;  
       fmcbb[3]= 0.0503634;   fmccc[3]= 0.101638;   fmcc[3]= 0.150538;   fmcll[3]= 0.697461;  
      Fbb[4]= 1.12722;   Fcc[4]= 1.12722;   Fc[4]= 0.885725;   Fll[4]= 0.987351; 
         cafac[4]= 1.11852;  
         winpretag[4]= 19686.2 ;   wintag[4]= 3014.29 ;   woutpretag[4]= 22019.5;   wouttag[4]= 3489.56;  
       fbb[4]= 0.0817505;   fcc[4]= 0.134911;   fc[4]= 0.12407;   fll[4]= 0.659269;  
       fmcbb[4]= 0.0725238;   fmccc[4]= 0.119684;   fmcc[4]= 0.140078;   fmcll[4]= 0.667714;  
      Fbb[5]= 1.11888;   Fcc[5]= 1.11888;   Fc[5]= 0.879172;   Fll[5]= 0.980046; 
         cafac[5]= 1.05125;  
         winpretag[5]= 5174.01 ;   wintag[5]= 1086.93 ;   woutpretag[5]= 5439.19;   wouttag[5]= 1181.02;  
       fbb[5]= 0.103712;   fcc[5]= 0.154449;   fc[5]= 0.105281;   fll[5]= 0.636558;  
       fmcbb[5]= 0.0926924;   fmccc[5]= 0.138039;   fmcc[5]= 0.119751;   fmcll[5]= 0.649518;  
      Fbb[6]= 1.13287;   Fcc[6]= 1.13287;   Fc[6]= 0.890163;   Fll[6]= 0.992298; 
         cafac[6]= 1.01469;  
         winpretag[6]= 120342 ;   wintag[6]= 15340.7 ;   woutpretag[6]= 122110;   wouttag[6]= 15969.2;  
       fbb[6]= 0.0632236;   fcc[6]= 0.12026;   fc[6]= 0.131302;   fll[6]= 0.685215;  
       fmcbb[6]= 0.0558084;   fmccc[6]= 0.106155;   fmcc[6]= 0.147503;   fmcll[6]= 0.690534;  
      Fbb[7]= 1.12548;   Fcc[7]= 1.12548;   Fc[7]= 0.884353;   Fll[7]= 0.985822; 
         cafac[7]= 1.1033;  
         winpretag[7]= 24860.2 ;   wintag[7]= 4101.22 ;   woutpretag[7]= 27428.3;   wouttag[7]= 4683.52;  
       fbb[7]= 0.0863481;   fcc[7]= 0.139001;   fc[7]= 0.120137;   fll[7]= 0.654514;  
       fmcbb[7]= 0.0767214;   fmccc[7]= 0.123504;   fmcc[7]= 0.135847;   fmcll[7]= 0.663927;  
 }  
 else if ( idsys==79   || sysname=="leff_up" ) { 
      Fbb[1]= 1.22641;   Fcc[1]= 1.22641;   Fc[1]= 0.949515;   Fll[1]= 0.99533; 
         cafac[1]= 0.984158;  
         winpretag[1]= 2.27306e+06 ;   wintag[1]= 81110.5 ;   woutpretag[1]= 2.23706e+06;   wouttag[1]= 80054.7;  
       fbb[1]= 0.0146262;   fcc[1]= 0.0424266;   fc[1]= 0.126008;   fll[1]= 0.816939;  
       fmcbb[1]= 0.011926;   fmccc[1]= 0.0345942;   fmcc[1]= 0.132708;   fmcll[1]= 0.820772;  
      Fbb[2]= 1.21234;   Fcc[2]= 1.21234;   Fc[2]= 0.938623;   Fll[2]= 0.983913; 
         cafac[2]= 0.905807;  
         winpretag[2]= 521776 ;   wintag[2]= 39249.7 ;   woutpretag[2]= 472629;   wouttag[2]= 36877;  
       fbb[2]= 0.0358835;   fcc[2]= 0.0862324;   fc[2]= 0.143459;   fll[2]= 0.734425;  
       fmcbb[2]= 0.0295985;   fmccc[2]= 0.0711288;   fmcc[2]= 0.15284;   fmcll[2]= 0.746433;  
      Fbb[3]= 1.19998;   Fcc[3]= 1.19998;   Fc[3]= 0.929053;   Fll[3]= 0.973881; 
         cafac[3]= 0.820347;  
         winpretag[3]= 116582 ;   wintag[3]= 13435.2 ;   woutpretag[3]= 95637.3;   wouttag[3]= 11609.1;  
       fbb[3]= 0.057936;   fcc[3]= 0.116727;   fc[3]= 0.140734;   fll[3]= 0.684603;  
       fmcbb[3]= 0.0482808;   fmccc[3]= 0.0972744;   fmcc[3]= 0.151481;   fmcll[3]= 0.702964;  
      Fbb[4]= 1.1884;   Fcc[4]= 1.1884;   Fc[4]= 0.920085;   Fll[4]= 0.96448; 
         cafac[4]= 0.875544;  
         winpretag[4]= 25363.7 ;   wintag[4]= 3803.91 ;   woutpretag[4]= 22207.1;   wouttag[4]= 3549.96;  
       fbb[4]= 0.0823988;   fcc[4]= 0.139712;   fc[4]= 0.131195;   fll[4]= 0.646695;  
       fmcbb[4]= 0.0693361;   fmccc[4]= 0.117563;   fmcc[4]= 0.14259;   fmcll[4]= 0.670511;  
      Fbb[5]= 1.17873;   Fcc[5]= 1.17873;   Fc[5]= 0.912598;   Fll[5]= 0.956632; 
         cafac[5]= 0.811768;  
         winpretag[5]= 6981.5 ;   wintag[5]= 1405.52 ;   woutpretag[5]= 5667.36;   wouttag[5]= 1210.54;  
       fbb[5]= 0.103834;   fcc[5]= 0.156021;   fc[5]= 0.115932;   fll[5]= 0.624213;  
       fmcbb[5]= 0.0880897;   fmccc[5]= 0.132364;   fmcc[5]= 0.127035;   fmcll[5]= 0.652511;  
      Fbb[6]= 1.19698;   Fcc[6]= 1.19698;   Fc[6]= 0.926731;   Fll[6]= 0.971447; 
         cafac[6]= 0.829836;  
         winpretag[6]= 148927 ;   wintag[6]= 18644.7 ;   woutpretag[6]= 123585;   wouttag[6]= 16359.7;  
       fbb[6]= 0.0643173;   fcc[6]= 0.122541;   fc[6]= 0.137917;   fll[6]= 0.675225;  
       fmcbb[6]= 0.0537329;   fmccc[6]= 0.102375;   fmcc[6]= 0.148821;   fmcll[6]= 0.695072;  
      Fbb[7]= 1.1863;   Fcc[7]= 1.1863;   Fc[7]= 0.918459;   Fll[7]= 0.962776; 
         cafac[7]= 0.86116;  
         winpretag[7]= 32345.2 ;   wintag[7]= 5209.43 ;   woutpretag[7]= 27854.4;   wouttag[7]= 4777.91;  
       fbb[7]= 0.0870551;   fcc[7]= 0.143255;   fc[7]= 0.127879;   fll[7]= 0.641811;  
       fmcbb[7]= 0.0733839;   fmccc[7]= 0.120758;   fmcc[7]= 0.139232;   fmcll[7]= 0.666626;  
 }  
 else if ( idsys==80   || sysname=="leff_down" ) { 
      Fbb[1]= 1.24073;   Fcc[1]= 1.24073;   Fc[1]= 0.972636;   Fll[1]= 0.990776; 
         cafac[1]= 1.01137;  
         winpretag[1]= 2.2151e+06 ;   wintag[1]= 78953.2 ;   woutpretag[1]= 2.2403e+06;   wouttag[1]= 81202.9;  
       fbb[1]= 0.0147936;   fcc[1]= 0.0429252;   fc[1]= 0.128897;   fll[1]= 0.813385;  
       fmcbb[1]= 0.0119232;   fmccc[1]= 0.0345966;   fmcc[1]= 0.132523;   fmcll[1]= 0.820957;  
      Fbb[2]= 1.22459;   Fcc[2]= 1.22459;   Fc[2]= 0.959978;   Fll[2]= 0.977882; 
         cafac[2]= 0.931977;  
         winpretag[2]= 508305 ;   wintag[2]= 38218.2 ;   woutpretag[2]= 473728;   wouttag[2]= 37340.8;  
       fbb[2]= 0.0362526;   fcc[2]= 0.0871024;   fc[2]= 0.146527;   fll[2]= 0.730118;  
       fmcbb[2]= 0.029604;   fmccc[2]= 0.071128;   fmcc[2]= 0.152636;   fmcll[2]= 0.746632;  
      Fbb[3]= 1.21116;   Fcc[3]= 1.21116;   Fc[3]= 0.949453;   Fll[3]= 0.96716; 
         cafac[3]= 0.844726;  
         winpretag[3]= 113553 ;   wintag[3]= 13082.1 ;   woutpretag[3]= 95921.6;   wouttag[3]= 11743.9;  
       fbb[3]= 0.0584828;   fcc[3]= 0.117825;   fc[3]= 0.143641;   fll[3]= 0.680051;  
       fmcbb[3]= 0.0482866;   fmccc[3]= 0.0972825;   fmcc[3]= 0.151288;   fmcll[3]= 0.703143;  
      Fbb[4]= 1.19889;   Fcc[4]= 1.19889;   Fc[4]= 0.939835;   Fll[4]= 0.957363; 
         cafac[4]= 0.901913;  
         winpretag[4]= 24700.2 ;   wintag[4]= 3702.13 ;   woutpretag[4]= 22277.4;   wouttag[4]= 3588.36;  
       fbb[4]= 0.0831096;   fcc[4]= 0.140921;   fc[4]= 0.133845;   fll[4]= 0.642125;  
       fmcbb[4]= 0.069322;   fmccc[4]= 0.117543;   fmcc[4]= 0.142413;   fmcll[4]= 0.670722;  
      Fbb[5]= 1.18896;   Fcc[5]= 1.18896;   Fc[5]= 0.932052;   Fll[5]= 0.949435; 
         cafac[5]= 0.836831;  
         winpretag[5]= 6799.31 ;   wintag[5]= 1368.97 ;   woutpretag[5]= 5689.87;   wouttag[5]= 1223.33;  
       fbb[5]= 0.104716;   fcc[5]= 0.157227;   fc[5]= 0.118286;   fll[5]= 0.619771;  
       fmcbb[5]= 0.0880729;   fmccc[5]= 0.132239;   fmcc[5]= 0.12691;   fmcll[5]= 0.652779;  
      Fbb[6]= 1.208;   Fcc[6]= 1.208;   Fc[6]= 0.946974;   Fll[6]= 0.964635; 
         cafac[6]= 0.854601;  
         winpretag[6]= 145053 ;   wintag[6]= 18153.2 ;   woutpretag[6]= 123962;   wouttag[6]= 16545.4;  
       fbb[6]= 0.06491;   fcc[6]= 0.123664;   fc[6]= 0.140753;   fll[6]= 0.670673;  
       fmcbb[6]= 0.0537335;   fmccc[6]= 0.102371;   fmcc[6]= 0.148634;   fmcll[6]= 0.695261;  
      Fbb[7]= 1.19674;   Fcc[7]= 1.19674;   Fc[7]= 0.938144;   Fll[7]= 0.95564; 
         cafac[7]= 0.887234;  
         winpretag[7]= 31499.5 ;   wintag[7]= 5071.1 ;   woutpretag[7]= 27947.4;   wouttag[7]= 4829.09;  
       fbb[7]= 0.0878038;   fcc[7]= 0.144464;   fc[7]= 0.130465;   fll[7]= 0.637268;  
       fmcbb[7]= 0.0733694;   fmccc[7]= 0.120715;   fmcc[7]= 0.139067;   fmcll[7]= 0.666849;  
 }  
 else if ( idsys==81   || sysname=="ca_up" ) { 
      Fbb[1]= 1.23358;   Fcc[1]= 1.23358;   Fc[1]= 0.961029;   Fll[1]= 0.993059; 
         cafac[1]= 1.00135;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.24712e+06;   wouttag[1]= 80931.5;  
       fbb[1]= 0.01471;   fcc[1]= 0.0426761;   fc[1]= 0.127449;   fll[1]= 0.815165;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.928425;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 478176;   wouttag[2]= 37500.5;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20558;   Fcc[3]= 1.20558;   Fc[3]= 0.93922;   Fll[3]= 0.970523; 
         cafac[3]= 0.847804;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 97554.7;   wouttag[3]= 11892.9;  
       fbb[3]= 0.05821;   fcc[3]= 0.117277;   fc[3]= 0.142185;   fll[3]= 0.682328;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19366;   Fcc[4]= 1.19366;   Fc[4]= 0.929929;   Fll[4]= 0.960923; 
         cafac[4]= 0.924449;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 23140.8;   wouttag[4]= 3713.34;  
       fbb[4]= 0.0827553;   fcc[4]= 0.140318;   fc[4]= 0.132517;   fll[4]= 0.644409;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18386;   Fcc[5]= 1.18386;   Fc[5]= 0.922294;   Fll[5]= 0.953033; 
         cafac[5]= 0.907317;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 6251.78;   wouttag[5]= 1339.76;  
       fbb[5]= 0.104276;   fcc[5]= 0.156627;   fc[5]= 0.117107;   fll[5]= 0.621991;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.2025;   Fcc[6]= 1.2025;   Fc[6]= 0.93682;   Fll[6]= 0.968043; 
         cafac[6]= 0.856257;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 125861;   wouttag[6]= 16729.9;  
       fbb[6]= 0.0646144;   fcc[6]= 0.123104;   fc[6]= 0.139332;   fll[6]= 0.67295;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.19153;   Fcc[7]= 1.19153;   Fc[7]= 0.928271;   Fll[7]= 0.959209; 
         cafac[7]= 0.907545;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 28971;   wouttag[7]= 4987.72;  
       fbb[7]= 0.0874306;   fcc[7]= 0.143861;   fc[7]= 0.129169;   fll[7]= 0.639539;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==82   || sysname=="ca_down" ) { 
      Fbb[1]= 1.23358;   Fcc[1]= 1.23358;   Fc[1]= 0.961029;   Fll[1]= 0.993059; 
         cafac[1]= 0.993848;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.23028e+06;   wouttag[1]= 80325.1;  
       fbb[1]= 0.01471;   fcc[1]= 0.0426761;   fc[1]= 0.127449;   fll[1]= 0.815165;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.909035;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 468190;   wouttag[2]= 36717.3;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20558;   Fcc[3]= 1.20558;   Fc[3]= 0.93922;   Fll[3]= 0.970523; 
         cafac[3]= 0.816964;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 94006;   wouttag[3]= 11460.2;  
       fbb[3]= 0.05821;   fcc[3]= 0.117277;   fc[3]= 0.142185;   fll[3]= 0.682328;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19366;   Fcc[4]= 1.19366;   Fc[4]= 0.929929;   Fll[4]= 0.960923; 
         cafac[4]= 0.852675;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 21344.1;   wouttag[4]= 3425.04;  
       fbb[4]= 0.0827553;   fcc[4]= 0.140318;   fc[4]= 0.132517;   fll[4]= 0.644409;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18386;   Fcc[5]= 1.18386;   Fc[5]= 0.922294;   Fll[5]= 0.953033; 
         cafac[5]= 0.740969;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5105.57;   wouttag[5]= 1094.12;  
       fbb[5]= 0.104276;   fcc[5]= 0.156627;   fc[5]= 0.117107;   fll[5]= 0.621991;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.2025;   Fcc[6]= 1.2025;   Fc[6]= 0.93682;   Fll[6]= 0.968043; 
         cafac[6]= 0.82787;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 121689;   wouttag[6]= 16175.3;  
       fbb[6]= 0.0646144;   fcc[6]= 0.123104;   fc[6]= 0.139332;   fll[6]= 0.67295;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.19153;   Fcc[7]= 1.19153;   Fc[7]= 0.928271;   Fll[7]= 0.959209; 
         cafac[7]= 0.84052;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 26831.4;   wouttag[7]= 4619.36;  
       fbb[7]= 0.0874306;   fcc[7]= 0.143861;   fc[7]= 0.129169;   fll[7]= 0.639539;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==83   || sysname=="chmis_up" ) { 
      Fbb[1]= 1.23358;   Fcc[1]= 1.23358;   Fc[1]= 0.961029;   Fll[1]= 0.993059; 
         cafac[1]= 1.00359;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.25213e+06;   wouttag[1]= 81112.1;  
       fbb[1]= 0.01471;   fcc[1]= 0.0426761;   fc[1]= 0.127449;   fll[1]= 0.815165;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.923323;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 475549;   wouttag[2]= 37294.5;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20558;   Fcc[3]= 1.20558;   Fc[3]= 0.93922;   Fll[3]= 0.970523; 
         cafac[3]= 0.838211;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 96450.8;   wouttag[3]= 11758.3;  
       fbb[3]= 0.05821;   fcc[3]= 0.117277;   fc[3]= 0.142185;   fll[3]= 0.682328;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19366;   Fcc[4]= 1.19366;   Fc[4]= 0.929929;   Fll[4]= 0.960923; 
         cafac[4]= 0.894782;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22398.1;   wouttag[4]= 3594.18;  
       fbb[4]= 0.0827553;   fcc[4]= 0.140318;   fc[4]= 0.132517;   fll[4]= 0.644409;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18386;   Fcc[5]= 1.18386;   Fc[5]= 0.922294;   Fll[5]= 0.953033; 
         cafac[5]= 0.830736;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5724.11;   wouttag[5]= 1226.68;  
       fbb[5]= 0.104276;   fcc[5]= 0.156627;   fc[5]= 0.117107;   fll[5]= 0.621991;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.2025;   Fcc[6]= 1.2025;   Fc[6]= 0.93682;   Fll[6]= 0.968043; 
         cafac[6]= 0.8488;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 124765;   wouttag[6]= 16584.2;  
       fbb[6]= 0.0646144;   fcc[6]= 0.123104;   fc[6]= 0.139332;   fll[6]= 0.67295;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.19153;   Fcc[7]= 1.19153;   Fc[7]= 0.928271;   Fll[7]= 0.959209; 
         cafac[7]= 0.881025;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 28124.4;   wouttag[7]= 4841.97;  
       fbb[7]= 0.0874306;   fcc[7]= 0.143861;   fc[7]= 0.129169;   fll[7]= 0.639539;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==84   || sysname=="chmis_down" ) { 
      Fbb[1]= 1.23358;   Fcc[1]= 1.23358;   Fc[1]= 0.961029;   Fll[1]= 0.993059; 
         cafac[1]= 0.991614;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.22526e+06;   wouttag[1]= 80144.5;  
       fbb[1]= 0.01471;   fcc[1]= 0.0426761;   fc[1]= 0.127449;   fll[1]= 0.815165;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21848;   Fcc[2]= 1.21848;   Fc[2]= 0.949264;   Fll[2]= 0.980902; 
         cafac[2]= 0.914136;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 470817;   wouttag[2]= 36923.4;  
       fbb[2]= 0.0360684;   fcc[2]= 0.0866682;   fc[2]= 0.14499;   fll[2]= 0.732274;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20558;   Fcc[3]= 1.20558;   Fc[3]= 0.93922;   Fll[3]= 0.970523; 
         cafac[3]= 0.826557;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95109.9;   wouttag[3]= 11594.8;  
       fbb[3]= 0.05821;   fcc[3]= 0.117277;   fc[3]= 0.142185;   fll[3]= 0.682328;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19366;   Fcc[4]= 1.19366;   Fc[4]= 0.929929;   Fll[4]= 0.960923; 
         cafac[4]= 0.882342;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22086.8;   wouttag[4]= 3544.21;  
       fbb[4]= 0.0827553;   fcc[4]= 0.140318;   fc[4]= 0.132517;   fll[4]= 0.644409;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.18386;   Fcc[5]= 1.18386;   Fc[5]= 0.922294;   Fll[5]= 0.953033; 
         cafac[5]= 0.81755;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5633.25;   wouttag[5]= 1207.21;  
       fbb[5]= 0.104276;   fcc[5]= 0.156627;   fc[5]= 0.117107;   fll[5]= 0.621991;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.2025;   Fcc[6]= 1.2025;   Fc[6]= 0.93682;   Fll[6]= 0.968043; 
         cafac[6]= 0.835327;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 122785;   wouttag[6]= 16321;  
       fbb[6]= 0.0646144;   fcc[6]= 0.123104;   fc[6]= 0.139332;   fll[6]= 0.67295;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.19153;   Fcc[7]= 1.19153;   Fc[7]= 0.928271;   Fll[7]= 0.959209; 
         cafac[7]= 0.86704;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27678;   wouttag[7]= 4765.11;  
       fbb[7]= 0.0874306;   fcc[7]= 0.143861;   fc[7]= 0.129169;   fll[7]= 0.639539;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==67   || sysname=="eer_up" ) { 
      Fbb[1]= 1.25496;   Fcc[1]= 1.25496;   Fc[1]= 0.958625;   Fll[1]= 0.992645; 
         cafac[1]= 0.998577;  
         winpretag[1]= 2.24112e+06 ;   wintag[1]= 79886.4 ;   woutpretag[1]= 2.23793e+06;   wouttag[1]= 80655.5;  
       fbb[1]= 0.0149823;   fcc[1]= 0.0418166;   fc[1]= 0.127297;   fll[1]= 0.815904;  
       fmcbb[1]= 0.0119385;   fmccc[1]= 0.0333211;   fmcc[1]= 0.132792;   fmcll[1]= 0.821949;  
      Fbb[2]= 1.23819;   Fcc[2]= 1.23819;   Fc[2]= 0.945817;   Fll[2]= 0.979383; 
         cafac[2]= 0.91934;  
         winpretag[2]= 514333 ;   wintag[2]= 38575.6 ;   woutpretag[2]= 472847;   wouttag[2]= 37109.5;  
       fbb[2]= 0.0367;   fcc[2]= 0.0864975;   fc[2]= 0.144663;   fll[2]= 0.732139;  
       fmcbb[2]= 0.02964;   fmccc[2]= 0.069858;   fmcc[2]= 0.152951;   fmcll[2]= 0.747551;  
      Fbb[3]= 1.22378;   Fcc[3]= 1.22378;   Fc[3]= 0.934812;   Fll[3]= 0.967988; 
         cafac[3]= 0.832689;  
         winpretag[3]= 114964 ;   wintag[3]= 13218.9 ;   woutpretag[3]= 95729.5;   wouttag[3]= 11700.7;  
       fbb[3]= 0.0591419;   fcc[3]= 0.118061;   fc[3]= 0.141644;   fll[3]= 0.681152;  
       fmcbb[3]= 0.0483271;   fmccc[3]= 0.0964723;   fmcc[3]= 0.151522;   fmcll[3]= 0.703679;  
      Fbb[4]= 1.21049;   Fcc[4]= 1.21049;   Fc[4]= 0.924654;   Fll[4]= 0.957469; 
         cafac[4]= 0.887703;  
         winpretag[4]= 25023.7 ;   wintag[4]= 3746.92 ;   woutpretag[4]= 22213.6;   wouttag[4]= 3580.2;  
       fbb[4]= 0.0839497;   fcc[4]= 0.141909;   fc[4]= 0.131809;   fll[4]= 0.642333;  
       fmcbb[4]= 0.0693521;   fmccc[4]= 0.117233;   fmcc[4]= 0.14255;   fmcll[4]= 0.670865;  
      Fbb[5]= 1.19969;   Fcc[5]= 1.19969;   Fc[5]= 0.91641;   Fll[5]= 0.948932; 
         cafac[5]= 0.822718;  
         winpretag[5]= 6888.13 ;   wintag[5]= 1386.46 ;   woutpretag[5]= 5666.99;   wouttag[5]= 1220.55;  
       fbb[5]= 0.105705;   fcc[5]= 0.158378;   fc[5]= 0.116398;   fll[5]= 0.619519;  
       fmcbb[5]= 0.0881105;   fmccc[5]= 0.132016;   fmcc[5]= 0.127015;   fmcll[5]= 0.652859;  
      Fbb[6]= 1.22035;   Fcc[6]= 1.22035;   Fc[6]= 0.93219;   Fll[6]= 0.965272; 
         cafac[6]= 0.842049;  
         winpretag[6]= 146876 ;   wintag[6]= 18352.3 ;   woutpretag[6]= 123677;   wouttag[6]= 16493.5;  
       fbb[6]= 0.0656242;   fcc[6]= 0.124081;   fc[6]= 0.138751;   fll[6]= 0.671544;  
       fmcbb[6]= 0.0537749;   fmccc[6]= 0.101676;   fmcc[6]= 0.148844;   fmcll[6]= 0.695705;  
      Fbb[7]= 1.20814;   Fcc[7]= 1.20814;   Fc[7]= 0.922862;   Fll[7]= 0.955613; 
         cafac[7]= 0.873029;  
         winpretag[7]= 31911.8 ;   wintag[7]= 5133.38 ;   woutpretag[7]= 27859.9;   wouttag[7]= 4818.51;  
       fbb[7]= 0.0886787;   fcc[7]= 0.145489;   fc[7]= 0.128459;   fll[7]= 0.637374;  
       fmcbb[7]= 0.0734011;   fmccc[7]= 0.120424;   fmcc[7]= 0.139197;   fmcll[7]= 0.666979;  
 }  
 else if ( idsys==68   || sysname=="eer_down" ) { 
      Fbb[1]= 1.23406;   Fcc[1]= 1.23406;   Fc[1]= 0.960708;   Fll[1]= 0.993083; 
         cafac[1]= 0.997551;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.23859e+06;   wouttag[1]= 80615.9;  
       fbb[1]= 0.0147157;   fcc[1]= 0.042693;   fc[1]= 0.127406;   fll[1]= 0.815185;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345957;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21893;   Fcc[2]= 1.21893;   Fc[2]= 0.94893;   Fll[2]= 0.980908; 
         cafac[2]= 0.918744;  
         winpretag[2]= 515037 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473187;   wouttag[2]= 37109.7;  
       fbb[2]= 0.0360819;   fcc[2]= 0.0867021;   fc[2]= 0.144941;   fll[2]= 0.732275;  
       fmcbb[2]= 0.0296014;   fmccc[2]= 0.0711299;   fmcc[2]= 0.152742;   fmcll[2]= 0.746527;  
      Fbb[3]= 1.20601;   Fcc[3]= 1.20601;   Fc[3]= 0.938871;   Fll[3]= 0.97051; 
         cafac[3]= 0.832328;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95773.5;   wouttag[3]= 11676.4;  
       fbb[3]= 0.0582306;   fcc[3]= 0.117319;   fc[3]= 0.142132;   fll[3]= 0.682318;  
       fmcbb[3]= 0.0482839;   fmccc[3]= 0.0972788;   fmcc[3]= 0.151387;   fmcll[3]= 0.703051;  
      Fbb[4]= 1.19405;   Fcc[4]= 1.19405;   Fc[4]= 0.929565;   Fll[4]= 0.960891; 
         cafac[4]= 0.888355;  
         winpretag[4]= 25032.7 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22237.9;   wouttag[4]= 3568.69;  
       fbb[4]= 0.0827803;   fcc[4]= 0.140361;   fc[4]= 0.132462;   fll[4]= 0.644397;  
       fmcbb[4]= 0.0693272;   fmccc[4]= 0.11755;   fmcc[4]= 0.142499;   fmcll[4]= 0.670625;  
      Fbb[5]= 1.18423;   Fcc[5]= 1.18423;   Fc[5]= 0.921915;   Fll[5]= 0.952983; 
         cafac[5]= 0.824327;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5679.95;   wouttag[5]= 1217.32;  
       fbb[5]= 0.104308;   fcc[5]= 0.156676;   fc[5]= 0.117058;   fll[5]= 0.621958;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.20292;   Fcc[6]= 1.20292;   Fc[6]= 0.936467;   Fll[6]= 0.968026; 
         cafac[6]= 0.841993;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123765;   wouttag[6]= 16452.4;  
       fbb[6]= 0.0646366;   fcc[6]= 0.123146;   fc[6]= 0.139279;   fll[6]= 0.672938;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695166;  
      Fbb[7]= 1.19192;   Fcc[7]= 1.19192;   Fc[7]= 0.927903;   Fll[7]= 0.959173; 
         cafac[7]= 0.873915;  
         winpretag[7]= 31923.1 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27898;   wouttag[7]= 4803.33;  
       fbb[7]= 0.0874572;   fcc[7]= 0.143905;   fc[7]= 0.129115;   fll[7]= 0.639523;  
       fmcbb[7]= 0.0733752;   fmccc[7]= 0.120734;   fmcc[7]= 0.139147;   fmcll[7]= 0.666743;  
 }  
 else if ( idsys==65   || sysname=="ees_up" ) { 
      Fbb[1]= 1.23416;   Fcc[1]= 1.23416;   Fc[1]= 0.960736;   Fll[1]= 0.993073; 
         cafac[1]= 0.997541;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.23856e+06;   wouttag[1]= 80617.3;  
       fbb[1]= 0.0147169;   fcc[1]= 0.0426963;   fc[1]= 0.127409;   fll[1]= 0.815177;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21902;   Fcc[2]= 1.21902;   Fc[2]= 0.948951;   Fll[2]= 0.980892; 
         cafac[2]= 0.918724;  
         winpretag[2]= 515036 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473176;   wouttag[2]= 37110.1;  
       fbb[2]= 0.0360849;   fcc[2]= 0.086708;   fc[2]= 0.144942;   fll[2]= 0.732265;  
       fmcbb[2]= 0.0296015;   fmccc[2]= 0.0711291;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.20609;   Fcc[3]= 1.20609;   Fc[3]= 0.938887;   Fll[3]= 0.970489; 
         cafac[3]= 0.832336;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95774.8;   wouttag[3]= 11676.9;  
       fbb[3]= 0.0582346;   fcc[3]= 0.117327;   fc[3]= 0.142134;   fll[3]= 0.682304;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19414;   Fcc[4]= 1.19414;   Fc[4]= 0.929577;   Fll[4]= 0.960866; 
         cafac[4]= 0.888457;  
         winpretag[4]= 25032.7 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22240.4;   wouttag[4]= 3569.22;  
       fbb[4]= 0.0827861;   fcc[4]= 0.14037;   fc[4]= 0.132463;   fll[4]= 0.64438;  
       fmcbb[4]= 0.0693272;   fmccc[4]= 0.11755;   fmcc[4]= 0.142499;   fmcll[4]= 0.670625;  
      Fbb[5]= 1.1843;   Fcc[5]= 1.1843;   Fc[5]= 0.921924;   Fll[5]= 0.952955; 
         cafac[5]= 0.824223;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5679.23;   wouttag[5]= 1217.21;  
       fbb[5]= 0.104315;   fcc[5]= 0.156686;   fc[5]= 0.117059;   fll[5]= 0.62194;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.20301;   Fcc[6]= 1.20301;   Fc[6]= 0.936482;   Fll[6]= 0.968003; 
         cafac[6]= 0.842012;  
         winpretag[6]= 146991 ;   wintag[6]= 18399 ;   woutpretag[6]= 123768;   wouttag[6]= 16453.3;  
       fbb[6]= 0.064641;   fcc[6]= 0.123155;   fc[6]= 0.139281;   fll[6]= 0.672923;  
       fmcbb[6]= 0.053733;   fmccc[6]= 0.102372;   fmcc[6]= 0.148728;   fmcll[6]= 0.695167;  
      Fbb[7]= 1.192;   Fcc[7]= 1.192;   Fc[7]= 0.927914;   Fll[7]= 0.959147; 
         cafac[7]= 0.87397;  
         winpretag[7]= 31923.1 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27899.8;   wouttag[7]= 4803.81;  
       fbb[7]= 0.0874632;   fcc[7]= 0.143915;   fc[7]= 0.129117;   fll[7]= 0.639505;  
       fmcbb[7]= 0.0733752;   fmccc[7]= 0.120734;   fmcc[7]= 0.139147;   fmcll[7]= 0.666743;  
 }  
 else if ( idsys==66   || sysname=="ees_down" ) { 
      Fbb[1]= 1.23372;   Fcc[1]= 1.23372;   Fc[1]= 0.960777;   Fll[1]= 0.993092; 
         cafac[1]= 0.997546;  
         winpretag[1]= 2.24409e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.23858e+06;   wouttag[1]= 80614.9;  
       fbb[1]= 0.0147116;   fcc[1]= 0.0426809;   fc[1]= 0.127415;   fll[1]= 0.815193;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345953;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21861;   Fcc[2]= 1.21861;   Fc[2]= 0.949015;   Fll[2]= 0.980934; 
         cafac[2]= 0.91871;  
         winpretag[2]= 515042 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473174;   wouttag[2]= 37106.3;  
       fbb[2]= 0.0360724;   fcc[2]= 0.0866779;   fc[2]= 0.144953;   fll[2]= 0.732297;  
       fmcbb[2]= 0.0296011;   fmccc[2]= 0.0711282;   fmcc[2]= 0.15274;   fmcll[2]= 0.74653;  
      Fbb[3]= 1.20571;   Fcc[3]= 1.20571;   Fc[3]= 0.938969;   Fll[3]= 0.97055; 
         cafac[3]= 0.832247;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95764.1;   wouttag[3]= 11674.4;  
       fbb[3]= 0.0582165;   fcc[3]= 0.11729;   fc[3]= 0.142147;   fll[3]= 0.682346;  
       fmcbb[3]= 0.0482839;   fmccc[3]= 0.0972788;   fmcc[3]= 0.151387;   fmcll[3]= 0.703051;  
      Fbb[4]= 1.19378;   Fcc[4]= 1.19378;   Fc[4]= 0.929676;   Fll[4]= 0.960945; 
         cafac[4]= 0.888604;  
         winpretag[4]= 25033.2 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22244.6;   wouttag[4]= 3569.37;  
       fbb[4]= 0.0827596;   fcc[4]= 0.140325;   fc[4]= 0.132474;   fll[4]= 0.644441;  
       fmcbb[4]= 0.0693256;   fmccc[4]= 0.117547;   fmcc[4]= 0.142495;   fmcll[4]= 0.670632;  
      Fbb[5]= 1.18397;   Fcc[5]= 1.18397;   Fc[5]= 0.922034;   Fll[5]= 0.953046; 
         cafac[5]= 0.824071;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.18;   wouttag[5]= 1216.85;  
       fbb[5]= 0.104286;   fcc[5]= 0.156642;   fc[5]= 0.117074;   fll[5]= 0.621999;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.20263;   Fcc[6]= 1.20263;   Fc[6]= 0.936568;   Fll[6]= 0.968068; 
         cafac[6]= 0.841963;  
         winpretag[6]= 146991 ;   wintag[6]= 18399 ;   woutpretag[6]= 123761;   wouttag[6]= 16450.5;  
       fbb[6]= 0.0646209;   fcc[6]= 0.123116;   fc[6]= 0.139294;   fll[6]= 0.672969;  
       fmcbb[6]= 0.0537329;   fmccc[6]= 0.102372;   fmcc[6]= 0.148728;   fmcll[6]= 0.695167;  
      Fbb[7]= 1.19165;   Fcc[7]= 1.19165;   Fc[7]= 0.928016;   Fll[7]= 0.959229; 
         cafac[7]= 0.874049;  
         winpretag[7]= 31923.6 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27902.8;   wouttag[7]= 4803.65;  
       fbb[7]= 0.0874359;   fcc[7]= 0.14387;   fc[7]= 0.129129;   fll[7]= 0.639565;  
       fmcbb[7]= 0.0733739;   fmccc[7]= 0.120732;   fmcc[7]= 0.139145;   fmcll[7]= 0.666749;  
 }  
 else if ( idsys==63   || sysname=="jer_one" ) { 
      Fbb[1]= 1.31793;   Fcc[1]= 1.31793;   Fc[1]= 0.948624;   Fll[1]= 0.99092; 
         cafac[1]= 0.942121;  
         winpretag[1]= 2.3449e+06 ;   wintag[1]= 81512.9 ;   woutpretag[1]= 2.20918e+06;   wouttag[1]= 77859.2;  
       fbb[1]= 0.0149783;   fcc[1]= 0.0434908;   fc[1]= 0.121741;   fll[1]= 0.81979;  
       fmcbb[1]= 0.011365;   fmccc[1]= 0.0329993;   fmcc[1]= 0.128334;   fmcll[1]= 0.827301;  
      Fbb[2]= 1.29645;   Fcc[2]= 1.29645;   Fc[2]= 0.933159;   Fll[2]= 0.974766; 
         cafac[2]= 0.883663;  
         winpretag[2]= 533773 ;   wintag[2]= 39757.3 ;   woutpretag[2]= 471675;   wouttag[2]= 37108.9;  
       fbb[2]= 0.0370962;   fcc[2]= 0.0900331;   fc[2]= 0.141515;   fll[2]= 0.731356;  
       fmcbb[2]= 0.0286138;   fmccc[2]= 0.0694461;   fmcc[2]= 0.151651;   fmcll[2]= 0.750289;  
      Fbb[3]= 1.27845;   Fcc[3]= 1.27845;   Fc[3]= 0.920209;   Fll[3]= 0.961238; 
         cafac[3]= 0.791356;  
         winpretag[3]= 120740 ;   wintag[3]= 13673.1 ;   woutpretag[3]= 95548.1;   wouttag[3]= 11646.2;  
       fbb[3]= 0.0598349;   fcc[3]= 0.12162;   fc[3]= 0.140432;   fll[3]= 0.678113;  
       fmcbb[3]= 0.0468025;   fmccc[3]= 0.0951303;   fmcc[3]= 0.152609;   fmcll[3]= 0.705458;  
      Fbb[4]= 1.26168;   Fcc[4]= 1.26168;   Fc[4]= 0.908132;   Fll[4]= 0.948623; 
         cafac[4]= 0.810573;  
         winpretag[4]= 26888.1 ;   wintag[4]= 3949.75 ;   woutpretag[4]= 21794.8;   wouttag[4]= 3486.61;  
       fbb[4]= 0.0832409;   fcc[4]= 0.147321;   fc[4]= 0.130778;   fll[4]= 0.63866;  
       fmcbb[4]= 0.0659764;   fmccc[4]= 0.116766;   fmcc[4]= 0.144008;   fmcll[4]= 0.673249;  
      Fbb[5]= 1.24565;   Fcc[5]= 1.24565;   Fc[5]= 0.896597;   Fll[5]= 0.936573; 
         cafac[5]= 0.728406;  
         winpretag[5]= 7504.66 ;   wintag[5]= 1495.33 ;   woutpretag[5]= 5466.44;   wouttag[5]= 1184.24;  
       fbb[5]= 0.11013;   fcc[5]= 0.165845;   fc[5]= 0.113248;   fll[5]= 0.610778;  
       fmcbb[5]= 0.0884115;   fmccc[5]= 0.133139;   fmcc[5]= 0.126308;   fmcll[5]= 0.652141;  
      Fbb[6]= 1.2739;   Fcc[6]= 1.2739;   Fc[6]= 0.916927;   Fll[6]= 0.95781; 
         cafac[6]= 0.791083;  
         winpretag[6]= 155132 ;   wintag[6]= 19118.1 ;   woutpretag[6]= 122723;   wouttag[6]= 16351.7;  
       fbb[6]= 0.0664192;   fcc[6]= 0.128306;   fc[6]= 0.137398;   fll[6]= 0.667878;  
       fmcbb[6]= 0.0521387;   fmccc[6]= 0.100719;   fmcc[6]= 0.149846;   fmcll[6]= 0.697297;  
      Fbb[7]= 1.25814;   Fcc[7]= 1.25814;   Fc[7]= 0.90559;   Fll[7]= 0.945967; 
         cafac[7]= 0.791407;  
         winpretag[7]= 34392.8 ;   wintag[7]= 5445.08 ;   woutpretag[7]= 27218.7;   wouttag[7]= 4694.21;  
       fbb[7]= 0.089167;   fcc[7]= 0.151404;   fc[7]= 0.126915;   fll[7]= 0.632515;  
       fmcbb[7]= 0.0708719;   fmccc[7]= 0.120339;   fmcc[7]= 0.140146;   fmcll[7]= 0.668643;  
 }  
 else if ( idsys==64   || sysname=="jef_one" ) { 
      Fbb[1]= 1.23559;   Fcc[1]= 1.23559;   Fc[1]= 0.95853;   Fll[1]= 0.993346; 
         cafac[1]= 0.996095;  
         winpretag[1]= 2.24217e+06 ;   wintag[1]= 79963.5 ;   woutpretag[1]= 2.23341e+06;   wouttag[1]= 80356.7;  
       fbb[1]= 0.0147469;   fcc[1]= 0.0427497;   fc[1]= 0.127145;   fll[1]= 0.815358;  
       fmcbb[1]= 0.0119351;   fmccc[1]= 0.0345986;   fmcc[1]= 0.132646;   fmcll[1]= 0.82082;  
      Fbb[2]= 1.22042;   Fcc[2]= 1.22042;   Fc[2]= 0.946758;   Fll[2]= 0.981147; 
         cafac[2]= 0.91939;  
         winpretag[2]= 514670 ;   wintag[2]= 38680.7 ;   woutpretag[2]= 473182;   wouttag[2]= 37076.5;  
       fbb[2]= 0.0361274;   fcc[2]= 0.086817;   fc[2]= 0.144547;   fll[2]= 0.732509;  
       fmcbb[2]= 0.0296026;   fmccc[2]= 0.0711373;   fmcc[2]= 0.152676;   fmcll[2]= 0.746584;  
      Fbb[3]= 1.20734;   Fcc[3]= 1.20734;   Fc[3]= 0.936619;   Fll[3]= 0.970639; 
         cafac[3]= 0.831849;  
         winpretag[3]= 114920 ;   wintag[3]= 13253.3 ;   woutpretag[3]= 95595.7;   wouttag[3]= 11663.4;  
       fbb[3]= 0.0583676;   fcc[3]= 0.117664;   fc[3]= 0.141807;   fll[3]= 0.682161;  
       fmcbb[3]= 0.0483438;   fmccc[3]= 0.097457;   fmcc[3]= 0.151403;   fmcll[3]= 0.702796;  
      Fbb[4]= 1.19547;   Fcc[4]= 1.19547;   Fc[4]= 0.927408;   Fll[4]= 0.961093; 
         cafac[4]= 0.891417;  
         winpretag[4]= 25002.2 ;   wintag[4]= 3745.08 ;   woutpretag[4]= 22287.4;   wouttag[4]= 3573.71;  
       fbb[4]= 0.082886;   fcc[4]= 0.140105;   fc[4]= 0.132478;   fll[4]= 0.644531;  
       fmcbb[4]= 0.0693333;   fmccc[4]= 0.117196;   fmcc[4]= 0.142848;   fmcll[4]= 0.670622;  
      Fbb[5]= 1.18536;   Fcc[5]= 1.18536;   Fc[5]= 0.919561;   Fll[5]= 0.952961; 
         cafac[5]= 0.827478;  
         winpretag[5]= 6877.66 ;   wintag[5]= 1382.7 ;   woutpretag[5]= 5691.11;   wouttag[5]= 1218.03;  
       fbb[5]= 0.104147;   fcc[5]= 0.157413;   fc[5]= 0.116768;   fll[5]= 0.621672;  
       fmcbb[5]= 0.0878613;   fmccc[5]= 0.132798;   fmcc[5]= 0.126983;   fmcll[5]= 0.652358;  
      Fbb[6]= 1.20426;   Fcc[6]= 1.20426;   Fc[6]= 0.934227;   Fll[6]= 0.96816; 
         cafac[6]= 0.842338;  
         winpretag[6]= 146799 ;   wintag[6]= 18381 ;   woutpretag[6]= 123655;   wouttag[6]= 16441.7;  
       fbb[6]= 0.0647532;   fcc[6]= 0.123406;   fc[6]= 0.139015;   fll[6]= 0.672826;  
       fmcbb[6]= 0.0537701;   fmccc[6]= 0.102475;   fmcc[6]= 0.148802;   fmcll[6]= 0.694953;  
      Fbb[7]= 1.19327;   Fcc[7]= 1.19327;   Fc[7]= 0.925704;   Fll[7]= 0.959327; 
         cafac[7]= 0.877005;  
         winpretag[7]= 31879.9 ;   wintag[7]= 5127.78 ;   woutpretag[7]= 27958.8;   wouttag[7]= 4809.02;  
       fbb[7]= 0.0875034;   fcc[7]= 0.143864;   fc[7]= 0.129066;   fll[7]= 0.639566;  
       fmcbb[7]= 0.0733305;   fmccc[7]= 0.120562;   fmcc[7]= 0.139425;   fmcll[7]= 0.666682;  
 }  
 else if ( idsys==61   || sysname=="musms_up" ) { 
      Fbb[1]= 1.22708;   Fcc[1]= 1.22708;   Fc[1]= 0.976432;   Fll[1]= 0.990939; 
         cafac[1]= 0.999519;  
         winpretag[1]= 2.24378e+06 ;   wintag[1]= 80032.7 ;   woutpretag[1]= 2.2427e+06;   wouttag[1]= 81366.1;  
       fbb[1]= 0.0146113;   fcc[1]= 0.0424672;   fc[1]= 0.129456;   fll[1]= 0.813466;  
       fmcbb[1]= 0.0119074;   fmccc[1]= 0.0346083;   fmcc[1]= 0.13258;   fmcll[1]= 0.820904;  
      Fbb[2]= 1.2119;   Fcc[2]= 1.2119;   Fc[2]= 0.96435;   Fll[2]= 0.978678; 
         cafac[2]= 0.918538;  
         winpretag[2]= 514692 ;   wintag[2]= 38645.1 ;   woutpretag[2]= 472764;   wouttag[2]= 37148.6;  
       fbb[2]= 0.0359113;   fcc[2]= 0.0862562;   fc[2]= 0.147269;   fll[2]= 0.730564;  
       fmcbb[2]= 0.0296323;   fmccc[2]= 0.0711745;   fmcc[2]= 0.152713;   fmcll[2]= 0.74648;  
      Fbb[3]= 1.19939;   Fcc[3]= 1.19939;   Fc[3]= 0.954396;   Fll[3]= 0.968575; 
         cafac[3]= 0.829742;  
         winpretag[3]= 114925 ;   wintag[3]= 13252.8 ;   woutpretag[3]= 95358.3;   wouttag[3]= 11651.5;  
       fbb[3]= 0.0578061;   fcc[3]= 0.116658;   fc[3]= 0.144674;   fll[3]= 0.680861;  
       fmcbb[3]= 0.0481963;   fmccc[3]= 0.097265;   fmcc[3]= 0.151587;   fmcll[3]= 0.702952;  
      Fbb[4]= 1.18786;   Fcc[4]= 1.18786;   Fc[4]= 0.945224;   Fll[4]= 0.959267; 
         cafac[4]= 0.887784;  
         winpretag[4]= 25025.5 ;   wintag[4]= 3747.58 ;   woutpretag[4]= 22217.2;   wouttag[4]= 3563.64;  
       fbb[4]= 0.0822379;   fcc[4]= 0.139829;   fc[4]= 0.13476;   fll[4]= 0.643173;  
       fmcbb[4]= 0.0692319;   fmccc[4]= 0.117715;   fmcc[4]= 0.142569;   fmcll[4]= 0.670484;  
      Fbb[5]= 1.17855;   Fcc[5]= 1.17855;   Fc[5]= 0.937813;   Fll[5]= 0.951746; 
         cafac[5]= 0.810527;  
         winpretag[5]= 6870.25 ;   wintag[5]= 1383.91 ;   woutpretag[5]= 5568.52;   wouttag[5]= 1194.4;  
       fbb[5]= 0.103897;   fcc[5]= 0.156038;   fc[5]= 0.119022;   fll[5]= 0.621044;  
       fmcbb[5]= 0.0881568;   fmccc[5]= 0.132398;   fmcc[5]= 0.126914;   fmcll[5]= 0.652531;  
      Fbb[6]= 1.19642;   Fcc[6]= 1.19642;   Fc[6]= 0.952033;   Fll[6]= 0.966178; 
         cafac[6]= 0.839141;  
         winpretag[6]= 146821 ;   wintag[6]= 18384.3 ;   woutpretag[6]= 123204;   wouttag[6]= 16402.3;  
       fbb[6]= 0.0641899;   fcc[6]= 0.122507;   fc[6]= 0.141754;   fll[6]= 0.67155;  
       fmcbb[6]= 0.0536517;   fmccc[6]= 0.102395;   fmcc[6]= 0.148896;   fmcll[6]= 0.695058;  
      Fbb[7]= 1.18584;   Fcc[7]= 1.18584;   Fc[7]= 0.943618;   Fll[7]= 0.957637; 
         cafac[7]= 0.870192;  
         winpretag[7]= 31895.8 ;   wintag[7]= 5131.49 ;   woutpretag[7]= 27755.4;   wouttag[7]= 4777.35;  
       fbb[7]= 0.0869321;   fcc[7]= 0.143342;   fc[7]= 0.131349;   fll[7]= 0.638377;  
       fmcbb[7]= 0.0733083;   fmccc[7]= 0.120877;   fmcc[7]= 0.139197;   fmcll[7]= 0.666617;  
 }  
 else if ( idsys==62   || sysname=="musms_down" ) { 
      Fbb[1]= 1.22535;   Fcc[1]= 1.22535;   Fc[1]= 0.977124;   Fll[1]= 0.990927; 
         cafac[1]= 0.999373;  
         winpretag[1]= 2.2435e+06 ;   wintag[1]= 79991.7 ;   woutpretag[1]= 2.24209e+06;   wouttag[1]= 81323.4;  
       fbb[1]= 0.014593;   fcc[1]= 0.0423984;   fc[1]= 0.129546;   fll[1]= 0.813463;  
       fmcbb[1]= 0.0119092;   fmccc[1]= 0.0346009;   fmcc[1]= 0.132579;   fmcll[1]= 0.820911;  
      Fbb[2]= 1.2103;   Fcc[2]= 1.2103;   Fc[2]= 0.965117;   Fll[2]= 0.97875; 
         cafac[2]= 0.918553;  
         winpretag[2]= 514664 ;   wintag[2]= 38653.1 ;   woutpretag[2]= 472746;   wouttag[2]= 37149.6;  
       fbb[2]= 0.03585;   fcc[2]= 0.0861075;   fc[2]= 0.14738;   fll[2]= 0.730662;  
       fmcbb[2]= 0.0296209;   fmccc[2]= 0.0711458;   fmcc[2]= 0.152707;   fmcll[2]= 0.746526;  
      Fbb[3]= 1.19788;   Fcc[3]= 1.19788;   Fc[3]= 0.955217;   Fll[3]= 0.96871; 
         cafac[3]= 0.828903;  
         winpretag[3]= 114915 ;   wintag[3]= 13245.4 ;   woutpretag[3]= 95253.6;   wouttag[3]= 11629.3;  
       fbb[3]= 0.0577591;   fcc[3]= 0.116478;   fc[3]= 0.144701;   fll[3]= 0.681062;  
       fmcbb[3]= 0.0482177;   fmccc[3]= 0.0972364;   fmcc[3]= 0.151485;   fmcll[3]= 0.703061;  
      Fbb[4]= 1.18646;   Fcc[4]= 1.18646;   Fc[4]= 0.946109;   Fll[4]= 0.959474; 
         cafac[4]= 0.889329;  
         winpretag[4]= 25013.4 ;   wintag[4]= 3748.68 ;   woutpretag[4]= 22245.1;   wouttag[4]= 3569.47;  
       fbb[4]= 0.0821038;   fcc[4]= 0.13968;   fc[4]= 0.134777;   fll[4]= 0.643439;  
       fmcbb[4]= 0.0692006;   fmccc[4]= 0.117728;   fmcc[4]= 0.142454;   fmcll[4]= 0.670617;  
      Fbb[5]= 1.17722;   Fcc[5]= 1.17722;   Fc[5]= 0.938739;   Fll[5]= 0.951999; 
         cafac[5]= 0.812537;  
         winpretag[5]= 6868.92 ;   wintag[5]= 1383.54 ;   woutpretag[5]= 5581.25;   wouttag[5]= 1196.38;  
       fbb[5]= 0.103816;   fcc[5]= 0.155906;   fc[5]= 0.119443;   fll[5]= 0.620835;  
       fmcbb[5]= 0.0881874;   fmccc[5]= 0.132436;   fmcc[5]= 0.127238;   fmcll[5]= 0.652139;  
      Fbb[6]= 1.19494;   Fcc[6]= 1.19494;   Fc[6]= 0.952871;   Fll[6]= 0.966331; 
         cafac[6]= 0.838883;  
         winpretag[6]= 146797 ;   wintag[6]= 18377.6 ;   woutpretag[6]= 123146;   wouttag[6]= 16385.2;  
       fbb[6]= 0.0641245;   fcc[6]= 0.122332;   fc[6]= 0.141798;   fll[6]= 0.671745;  
       fmcbb[6]= 0.0536633;   fmccc[6]= 0.102375;   fmcc[6]= 0.148811;   fmcll[6]= 0.69515;  
      Fbb[7]= 1.18446;   Fcc[7]= 1.18446;   Fc[7]= 0.944511;   Fll[7]= 0.957853; 
         cafac[7]= 0.871856;  
         winpretag[7]= 31882.3 ;   wintag[7]= 5132.22 ;   woutpretag[7]= 27796.8;   wouttag[7]= 4785.05;  
       fbb[7]= 0.0868103;   fcc[7]= 0.143197;   fc[7]= 0.131453;   fll[7]= 0.63854;  
       fmcbb[7]= 0.0732913;   fmccc[7]= 0.120897;   fmcc[7]= 0.139176;   fmcll[7]= 0.666636;  
 }  
 else if ( idsys==59   || sysname=="musid_up" ) { 
      Fbb[1]= 1.22353;   Fcc[1]= 1.22353;   Fc[1]= 0.977691;   Fll[1]= 0.990938; 
         cafac[1]= 0.999297;  
         winpretag[1]= 2.24368e+06 ;   wintag[1]= 80016.4 ;   woutpretag[1]= 2.2421e+06;   wouttag[1]= 81347;  
       fbb[1]= 0.0145685;   fcc[1]= 0.042342;   fc[1]= 0.129617;   fll[1]= 0.813472;  
       fmcbb[1]= 0.0119069;   fmccc[1]= 0.0346065;   fmcc[1]= 0.132575;   fmcll[1]= 0.820912;  
      Fbb[2]= 1.20859;   Fcc[2]= 1.20859;   Fc[2]= 0.965757;   Fll[2]= 0.978842; 
         cafac[2]= 0.918689;  
         winpretag[2]= 514695 ;   wintag[2]= 38654.6 ;   woutpretag[2]= 472844;   wouttag[2]= 37147.5;  
       fbb[2]= 0.035799;   fcc[2]= 0.0860135;   fc[2]= 0.147501;   fll[2]= 0.730686;  
       fmcbb[2]= 0.0296204;   fmccc[2]= 0.0711683;   fmcc[2]= 0.152731;   fmcll[2]= 0.74648;  
      Fbb[3]= 1.19629;   Fcc[3]= 1.19629;   Fc[3]= 0.95593;   Fll[3]= 0.968882; 
         cafac[3]= 0.828902;  
         winpretag[3]= 114941 ;   wintag[3]= 13251.7 ;   woutpretag[3]= 95274.9;   wouttag[3]= 11630.9;  
       fbb[3]= 0.0576766;   fcc[3]= 0.116342;   fc[3]= 0.144833;   fll[3]= 0.681149;  
       fmcbb[3]= 0.0482128;   fmccc[3]= 0.0972517;   fmcc[3]= 0.15151;   fmcll[3]= 0.703025;  
      Fbb[4]= 1.18498;   Fcc[4]= 1.18498;   Fc[4]= 0.946892;   Fll[4]= 0.959722; 
         cafac[4]= 0.889221;  
         winpretag[4]= 25016.8 ;   wintag[4]= 3746.7 ;   woutpretag[4]= 22245.4;   wouttag[4]= 3565.33;  
       fbb[4]= 0.0819808;   fcc[4]= 0.139515;   fc[4]= 0.134892;   fll[4]= 0.643612;  
       fmcbb[4]= 0.0691831;   fmccc[4]= 0.117736;   fmcc[4]= 0.142457;   fmcll[4]= 0.670624;  
      Fbb[5]= 1.17584;   Fcc[5]= 1.17584;   Fc[5]= 0.939589;   Fll[5]= 0.95232; 
         cafac[5]= 0.811836;  
         winpretag[5]= 6871.04 ;   wintag[5]= 1384.19 ;   woutpretag[5]= 5578.16;   wouttag[5]= 1195.48;  
       fbb[5]= 0.103555;   fcc[5]= 0.155781;   fc[5]= 0.119514;   fll[5]= 0.621149;  
       fmcbb[5]= 0.0880686;   fmccc[5]= 0.132485;   fmcc[5]= 0.127199;   fmcll[5]= 0.652248;  
      Fbb[6]= 1.19338;   Fcc[6]= 1.19338;   Fc[6]= 0.953603;   Fll[6]= 0.966524; 
         cafac[6]= 0.838826;  
         winpretag[6]= 146829 ;   wintag[6]= 18382.6 ;   woutpretag[6]= 123164;   wouttag[6]= 16382.3;  
       fbb[6]= 0.0640259;   fcc[6]= 0.122191;   fc[6]= 0.141925;   fll[6]= 0.671858;  
       fmcbb[6]= 0.0536508;   fmccc[6]= 0.102391;   fmcc[6]= 0.14883;   fmcll[6]= 0.695128;  
      Fbb[7]= 1.183;   Fcc[7]= 1.183;   Fc[7]= 0.945309;   Fll[7]= 0.958117; 
         cafac[7]= 0.871605;  
         winpretag[7]= 31887.8 ;   wintag[7]= 5130.89 ;   woutpretag[7]= 27793.6;   wouttag[7]= 4780.17;  
       fbb[7]= 0.0866579;   fcc[7]= 0.143042;   fc[7]= 0.131558;   fll[7]= 0.638742;  
       fmcbb[7]= 0.0732525;   fmccc[7]= 0.120914;   fmcc[7]= 0.139169;   fmcll[7]= 0.666664;  
 }  
 else if ( idsys==60   || sysname=="musid_down" ) { 
      Fbb[1]= 1.22571;   Fcc[1]= 1.22571;   Fc[1]= 0.976803;   Fll[1]= 0.990958; 
         cafac[1]= 0.999363;  
         winpretag[1]= 2.24366e+06 ;   wintag[1]= 79999.3 ;   woutpretag[1]= 2.24223e+06;   wouttag[1]= 81320.1;  
       fbb[1]= 0.0145957;   fcc[1]= 0.0424124;   fc[1]= 0.129495;   fll[1]= 0.813497;  
       fmcbb[1]= 0.0119079;   fmccc[1]= 0.0346024;   fmcc[1]= 0.13257;   fmcll[1]= 0.82092;  
      Fbb[2]= 1.21063;   Fcc[2]= 1.21063;   Fc[2]= 0.964789;   Fll[2]= 0.97877; 
         cafac[2]= 0.918524;  
         winpretag[2]= 514672 ;   wintag[2]= 38650.8 ;   woutpretag[2]= 472739;   wouttag[2]= 37145.6;  
       fbb[2]= 0.0358522;   fcc[2]= 0.0861447;   fc[2]= 0.147327;   fll[2]= 0.730676;  
       fmcbb[2]= 0.0296144;   fmccc[2]= 0.0711567;   fmcc[2]= 0.152704;   fmcll[2]= 0.746525;  
      Fbb[3]= 1.19819;   Fcc[3]= 1.19819;   Fc[3]= 0.954873;   Fll[3]= 0.96871; 
         cafac[3]= 0.829405;  
         winpretag[3]= 114910 ;   wintag[3]= 13247.2 ;   woutpretag[3]= 95307.2;   wouttag[3]= 11638.7;  
       fbb[3]= 0.0577694;   fcc[3]= 0.116552;   fc[3]= 0.144679;   fll[3]= 0.681;  
       fmcbb[3]= 0.0482138;   fmccc[3]= 0.0972732;   fmcc[3]= 0.151516;   fmcll[3]= 0.702997;  
      Fbb[4]= 1.18674;   Fcc[4]= 1.18674;   Fc[4]= 0.94575;   Fll[4]= 0.959454; 
         cafac[4]= 0.888699;  
         winpretag[4]= 25019.1 ;   wintag[4]= 3748.2 ;   woutpretag[4]= 22234.5;   wouttag[4]= 3566.36;  
       fbb[4]= 0.0821849;   fcc[4]= 0.139707;   fc[4]= 0.134697;   fll[4]= 0.643411;  
       fmcbb[4]= 0.0692525;   fmccc[4]= 0.117723;   fmcc[4]= 0.142423;   fmcll[4]= 0.670601;  
      Fbb[5]= 1.17751;   Fcc[5]= 1.17751;   Fc[5]= 0.938392;   Fll[5]= 0.951991; 
         cafac[5]= 0.811633;  
         winpretag[5]= 6868.57 ;   wintag[5]= 1384.19 ;   woutpretag[5]= 5574.75;   wouttag[5]= 1195.77;  
       fbb[5]= 0.103739;   fcc[5]= 0.155967;   fc[5]= 0.119405;   fll[5]= 0.620889;  
       fmcbb[5]= 0.0881003;   fmccc[5]= 0.132455;   fmcc[5]= 0.127244;   fmcll[5]= 0.6522;  
      Fbb[6]= 1.19524;   Fcc[6]= 1.19524;   Fc[6]= 0.952524;   Fll[6]= 0.966327; 
         cafac[6]= 0.839111;  
         winpretag[6]= 146798 ;   wintag[6]= 18379.6 ;   woutpretag[6]= 123180;   wouttag[6]= 16392.3;  
       fbb[6]= 0.0641437;   fcc[6]= 0.122399;   fc[6]= 0.141765;   fll[6]= 0.671693;  
       fmcbb[6]= 0.0536658;   fmccc[6]= 0.102405;   fmcc[6]= 0.148831;   fmcll[6]= 0.695099;  
      Fbb[7]= 1.18474;   Fcc[7]= 1.18474;   Fc[7]= 0.944155;   Fll[7]= 0.957837; 
         cafac[7]= 0.871155;  
         winpretag[7]= 31887.7 ;   wintag[7]= 5132.39 ;   woutpretag[7]= 27779.1;   wouttag[7]= 4781.4;  
       fbb[7]= 0.0868562;   fcc[7]= 0.143231;   fc[7]= 0.131383;   fll[7]= 0.63853;  
       fmcbb[7]= 0.0733123;   fmccc[7]= 0.120897;   fmcc[7]= 0.139154;   fmcll[7]= 0.666637;  
 }  
 else if ( idsys==57   || sysname=="mscale_one" ) { 
      Fbb[1]= 1.23206;   Fcc[1]= 1.23206;   Fc[1]= 0.957495;   Fll[1]= 0.993716; 
         cafac[1]= 0.997237;  
         winpretag[1]= 2.24829e+06 ;   wintag[1]= 80180.5 ;   woutpretag[1]= 2.24208e+06;   wouttag[1]= 80584.7;  
       fbb[1]= 0.0146814;   fcc[1]= 0.0426288;   fc[1]= 0.126972;   fll[1]= 0.815718;  
       fmcbb[1]= 0.0119161;   fmccc[1]= 0.0345995;   fmcc[1]= 0.132608;   fmcll[1]= 0.820876;  
      Fbb[2]= 1.21722;   Fcc[2]= 1.21722;   Fc[2]= 0.945959;   Fll[2]= 0.981744; 
         cafac[2]= 0.917974;  
         winpretag[2]= 515996 ;   wintag[2]= 38798.8 ;   woutpretag[2]= 473671;   wouttag[2]= 37085.2;  
       fbb[2]= 0.0360385;   fcc[2]= 0.0865713;   fc[2]= 0.144416;   fll[2]= 0.732974;  
       fmcbb[2]= 0.0296072;   fmccc[2]= 0.0711221;   fmcc[2]= 0.152666;   fmcll[2]= 0.746604;  
      Fbb[3]= 1.20446;   Fcc[3]= 1.20446;   Fc[3]= 0.93604;   Fll[3]= 0.97145; 
         cafac[3]= 0.831427;  
         winpretag[3]= 115262 ;   wintag[3]= 13261.9 ;   woutpretag[3]= 95831.5;   wouttag[3]= 11652.4;  
       fbb[3]= 0.0581283;   fcc[3]= 0.117144;   fc[3]= 0.14162;   fll[3]= 0.683108;  
       fmcbb[3]= 0.048261;   fmccc[3]= 0.0972587;   fmcc[3]= 0.151297;   fmcll[3]= 0.703183;  
      Fbb[4]= 1.1926;   Fcc[4]= 1.1926;   Fc[4]= 0.926823;   Fll[4]= 0.961884; 
         cafac[4]= 0.888202;  
         winpretag[4]= 25062.3 ;   wintag[4]= 3759.76 ;   woutpretag[4]= 22260.4;   wouttag[4]= 3570.58;  
       fbb[4]= 0.0827271;   fcc[4]= 0.140131;   fc[4]= 0.132081;   fll[4]= 0.645061;  
       fmcbb[4]= 0.0693673;   fmccc[4]= 0.117501;   fmcc[4]= 0.142509;   fmcll[4]= 0.670622;  
      Fbb[5]= 1.18285;   Fcc[5]= 1.18285;   Fc[5]= 0.919252;   Fll[5]= 0.954027; 
         cafac[5]= 0.824848;  
         winpretag[5]= 6895.28 ;   wintag[5]= 1387.06 ;   woutpretag[5]= 5687.56;   wouttag[5]= 1216.87;  
       fbb[5]= 0.104433;   fcc[5]= 0.156039;   fc[5]= 0.116734;   fll[5]= 0.622794;  
       fmcbb[5]= 0.0882892;   fmccc[5]= 0.131917;   fmcc[5]= 0.126989;   fmcll[5]= 0.652805;  
      Fbb[6]= 1.2014;   Fcc[6]= 1.2014;   Fc[6]= 0.933661;   Fll[6]= 0.968981; 
         cafac[6]= 0.841297;  
         winpretag[6]= 147219 ;   wintag[6]= 18408.7 ;   woutpretag[6]= 123855;   wouttag[6]= 16428.4;  
       fbb[6]= 0.0645496;   fcc[6]= 0.122936;   fc[6]= 0.1388;   fll[6]= 0.673714;  
       fmcbb[6]= 0.0537289;   fmccc[6]= 0.102328;   fmcc[6]= 0.148662;   fmcll[6]= 0.695281;  
      Fbb[7]= 1.19048;   Fcc[7]= 1.19048;   Fc[7]= 0.925179;   Fll[7]= 0.960178; 
         cafac[7]= 0.873918;  
         winpretag[7]= 31957.6 ;   wintag[7]= 5146.83 ;   woutpretag[7]= 27928.3;   wouttag[7]= 4804.5;  
       fbb[7]= 0.0874407;   fcc[7]= 0.143586;   fc[7]= 0.128748;   fll[7]= 0.640225;  
       fmcbb[7]= 0.0734499;   fmccc[7]= 0.120612;   fmcc[7]= 0.13916;   fmcll[7]= 0.666778;  
 }  
 else if ( idsys==55   || sysname=="metpileup_up" ) { 
      Fbb[1]= 1.26135;   Fcc[1]= 1.26135;   Fc[1]= 0.945856;   Fll[1]= 0.993947; 
         cafac[1]= 0.993671;  
         winpretag[1]= 2.24722e+06 ;   wintag[1]= 80103.4 ;   woutpretag[1]= 2.23299e+06;   wouttag[1]= 80026;  
       fbb[1]= 0.0150223;   fcc[1]= 0.0436176;   fc[1]= 0.125458;   fll[1]= 0.815902;  
       fmcbb[1]= 0.0119097;   fmccc[1]= 0.0345799;   fmcc[1]= 0.132639;   fmcll[1]= 0.820871;  
      Fbb[2]= 1.2445;   Fcc[2]= 1.2445;   Fc[2]= 0.933219;   Fll[2]= 0.980667; 
         cafac[2]= 0.916114;  
         winpretag[2]= 515407 ;   wintag[2]= 38774.8 ;   woutpretag[2]= 472172;   wouttag[2]= 37094.3;  
       fbb[2]= 0.0368609;   fcc[2]= 0.0885188;   fc[2]= 0.142552;   fll[2]= 0.732068;  
       fmcbb[2]= 0.029619;   fmccc[2]= 0.0711278;   fmcc[2]= 0.152753;   fmcll[2]= 0.7465;  
      Fbb[3]= 1.22989;   Fcc[3]= 1.22989;   Fc[3]= 0.922264;   Fll[3]= 0.969155; 
         cafac[3]= 0.83083;  
         winpretag[3]= 115133 ;   wintag[3]= 13263 ;   woutpretag[3]= 95655.6;   wouttag[3]= 11706;  
       fbb[3]= 0.0593415;   fcc[3]= 0.119675;   fc[3]= 0.139774;   fll[3]= 0.681209;  
       fmcbb[3]= 0.0482493;   fmccc[3]= 0.0973053;   fmcc[3]= 0.151556;   fmcll[3]= 0.70289;  
      Fbb[4]= 1.21636;   Fcc[4]= 1.21636;   Fc[4]= 0.912117;   Fll[4]= 0.958493; 
         cafac[4]= 0.887179;  
         winpretag[4]= 25065.4 ;   wintag[4]= 3747.62 ;   woutpretag[4]= 22237.5;   wouttag[4]= 3579.13;  
       fbb[4]= 0.0840825;   fcc[4]= 0.142829;   fc[4]= 0.129773;   fll[4]= 0.643316;  
       fmcbb[4]= 0.0691262;   fmccc[4]= 0.117423;   fmcc[4]= 0.142277;   fmcll[4]= 0.671174;  
      Fbb[5]= 1.20495;   Fcc[5]= 1.20495;   Fc[5]= 0.903562;   Fll[5]= 0.949502; 
         cafac[5]= 0.818634;  
         winpretag[5]= 6894.03 ;   wintag[5]= 1384.72 ;   woutpretag[5]= 5643.68;   wouttag[5]= 1213.79;  
       fbb[5]= 0.106377;   fcc[5]= 0.159429;   fc[5]= 0.115125;   fll[5]= 0.619069;  
       fmcbb[5]= 0.0882828;   fmccc[5]= 0.132312;   fmcc[5]= 0.127413;   fmcll[5]= 0.651993;  
      Fbb[6]= 1.22638;   Fcc[6]= 1.22638;   Fc[6]= 0.919628;   Fll[6]= 0.966386; 
         cafac[6]= 0.840294;  
         winpretag[6]= 147092 ;   wintag[6]= 18395.3 ;   woutpretag[6]= 123601;   wouttag[6]= 16491.5;  
       fbb[6]= 0.065836;   fcc[6]= 0.125549;   fc[6]= 0.13688;   fll[6]= 0.671734;  
       fmcbb[6]= 0.0536832;   fmccc[6]= 0.102374;   fmcc[6]= 0.148843;   fmcll[6]= 0.6951;  
      Fbb[7]= 1.21388;   Fcc[7]= 1.21388;   Fc[7]= 0.910258;   Fll[7]= 0.956539; 
         cafac[7]= 0.871631;  
         winpretag[7]= 31959.5 ;   wintag[7]= 5132.34 ;   woutpretag[7]= 27856.8;   wouttag[7]= 4811.18;  
       fbb[7]= 0.0889273;   fcc[7]= 0.146436;   fc[7]= 0.12659;   fll[7]= 0.638046;  
       fmcbb[7]= 0.0732585;   fmccc[7]= 0.120635;   fmcc[7]= 0.13907;   fmcll[7]= 0.667036;  
 }  
 else if ( idsys==56   || sysname=="metpileup_down" ) { 
      Fbb[1]= 1.2402;   Fcc[1]= 1.2402;   Fc[1]= 0.959486;   Fll[1]= 0.99294; 
         cafac[1]= 0.99888;  
         winpretag[1]= 2.24118e+06 ;   wintag[1]= 79860.6 ;   woutpretag[1]= 2.23867e+06;   wouttag[1]= 80565.9;  
       fbb[1]= 0.0147855;   fcc[1]= 0.0428788;   fc[1]= 0.127232;   fll[1]= 0.815103;  
       fmcbb[1]= 0.0119219;   fmccc[1]= 0.0345741;   fmcc[1]= 0.132605;   fmcll[1]= 0.820899;  
      Fbb[2]= 1.22459;   Fcc[2]= 1.22459;   Fc[2]= 0.947409;   Fll[2]= 0.980442; 
         cafac[2]= 0.918964;  
         winpretag[2]= 514854 ;   wintag[2]= 38689.3 ;   woutpretag[2]= 473133;   wouttag[2]= 37116.3;  
       fbb[2]= 0.0362623;   fcc[2]= 0.0871266;   fc[2]= 0.144615;   fll[2]= 0.731996;  
       fmcbb[2]= 0.0296117;   fmccc[2]= 0.0711476;   fmcc[2]= 0.152642;   fmcll[2]= 0.746598;  
      Fbb[3]= 1.21129;   Fcc[3]= 1.21129;   Fc[3]= 0.937122;   Fll[3]= 0.969796; 
         cafac[3]= 0.830283;  
         winpretag[3]= 114998 ;   wintag[3]= 13248.3 ;   woutpretag[3]= 95480.6;   wouttag[3]= 11653.7;  
       fbb[3]= 0.058418;   fcc[3]= 0.117882;   fc[3]= 0.141842;   fll[3]= 0.681858;  
       fmcbb[3]= 0.0482278;   fmccc[3]= 0.0973189;   fmcc[3]= 0.151359;   fmcll[3]= 0.703094;  
      Fbb[4]= 1.19891;   Fcc[4]= 1.19891;   Fc[4]= 0.92754;   Fll[4]= 0.95988; 
         cafac[4]= 0.887215;  
         winpretag[4]= 25027.9 ;   wintag[4]= 3753.89 ;   woutpretag[4]= 22205.1;   wouttag[4]= 3570.33;  
       fbb[4]= 0.083383;   fcc[4]= 0.14092;   fc[4]= 0.131928;   fll[4]= 0.643769;  
       fmcbb[4]= 0.0695491;   fmccc[4]= 0.11754;   fmcc[4]= 0.142234;   fmcll[4]= 0.670676;  
      Fbb[5]= 1.18886;   Fcc[5]= 1.18886;   Fc[5]= 0.919765;   Fll[5]= 0.951834; 
         cafac[5]= 0.828101;  
         winpretag[5]= 6884.24 ;   wintag[5]= 1390.31 ;   woutpretag[5]= 5700.85;   wouttag[5]= 1227.14;  
       fbb[5]= 0.105149;   fcc[5]= 0.156909;   fc[5]= 0.117049;   fll[5]= 0.620893;  
       fmcbb[5]= 0.0884453;   fmccc[5]= 0.131983;   fmcc[5]= 0.127259;   fmcll[5]= 0.652313;  
      Fbb[6]= 1.2081;   Fcc[6]= 1.2081;   Fc[6]= 0.934651;   Fll[6]= 0.967239; 
         cafac[6]= 0.84041;  
         winpretag[6]= 146910 ;   wintag[6]= 18392.5 ;   woutpretag[6]= 123464;   wouttag[6]= 16438.1;  
       fbb[6]= 0.064929;   fcc[6]= 0.123695;   fc[6]= 0.138959;   fll[6]= 0.672416;  
       fmcbb[6]= 0.0537447;   fmccc[6]= 0.102388;   fmcc[6]= 0.148675;   fmcll[6]= 0.695192;  
      Fbb[7]= 1.19673;   Fcc[7]= 1.19673;   Fc[7]= 0.925852;   Fll[7]= 0.958133; 
         cafac[7]= 0.873911;  
         winpretag[7]= 31912.1 ;   wintag[7]= 5144.19 ;   woutpretag[7]= 27888.4;   wouttag[7]= 4813.99;  
       fbb[7]= 0.0881096;   fcc[7]= 0.144392;   fc[7]= 0.128697;   fll[7]= 0.638801;  
       fmcbb[7]= 0.0736255;   fmccc[7]= 0.120656;   fmcc[7]= 0.139004;   fmcll[7]= 0.666715;  
 }  
 else if ( idsys==53   || sysname=="met_up" ) { 
      Fbb[1]= 1.23422;   Fcc[1]= 1.23422;   Fc[1]= 0.960591;   Fll[1]= 0.993093; 
         cafac[1]= 0.997516;  
         winpretag[1]= 2.24409e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.23851e+06;   wouttag[1]= 80609.8;  
       fbb[1]= 0.0147176;   fcc[1]= 0.0426987;   fc[1]= 0.12739;   fll[1]= 0.815194;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345956;   fmcc[1]= 0.132616;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21909;   Fcc[2]= 1.21909;   Fc[2]= 0.948809;   Fll[2]= 0.980912; 
         cafac[2]= 0.918693;  
         winpretag[2]= 515038 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473162;   wouttag[2]= 37107.7;  
       fbb[2]= 0.0360866;   fcc[2]= 0.086712;   fc[2]= 0.144922;   fll[2]= 0.732279;  
       fmcbb[2]= 0.0296014;   fmccc[2]= 0.0711288;   fmcc[2]= 0.152741;   fmcll[2]= 0.746528;  
      Fbb[3]= 1.20615;   Fcc[3]= 1.20615;   Fc[3]= 0.938744;   Fll[3]= 0.970507; 
         cafac[3]= 0.832245;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95764.3;   wouttag[3]= 11675.5;  
       fbb[3]= 0.0582375;   fcc[3]= 0.117333;   fc[3]= 0.142113;   fll[3]= 0.682317;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.19419;   Fcc[4]= 1.19419;   Fc[4]= 0.929434;   Fll[4]= 0.960882; 
         cafac[4]= 0.88847;  
         winpretag[4]= 25033.2 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22241.3;   wouttag[4]= 3569.26;  
       fbb[4]= 0.082788;   fcc[4]= 0.140374;   fc[4]= 0.13244;   fll[4]= 0.644398;  
       fmcbb[4]= 0.0693256;   fmccc[4]= 0.117547;   fmcc[4]= 0.142495;   fmcll[4]= 0.670632;  
      Fbb[5]= 1.18435;   Fcc[5]= 1.18435;   Fc[5]= 0.921777;   Fll[5]= 0.952966; 
         cafac[5]= 0.82427;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5679.55;   wouttag[5]= 1217.28;  
       fbb[5]= 0.10432;   fcc[5]= 0.156693;   fc[5]= 0.117041;   fll[5]= 0.621947;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.20306;   Fcc[6]= 1.20306;   Fc[6]= 0.936339;   Fll[6]= 0.968021; 
         cafac[6]= 0.841947;  
         winpretag[6]= 146991 ;   wintag[6]= 18399 ;   woutpretag[6]= 123759;   wouttag[6]= 16451.8;  
       fbb[6]= 0.0646439;   fcc[6]= 0.12316;   fc[6]= 0.139259;   fll[6]= 0.672937;  
       fmcbb[6]= 0.0537327;   fmccc[6]= 0.102372;   fmcc[6]= 0.148727;   fmcll[6]= 0.695168;  
      Fbb[7]= 1.19205;   Fcc[7]= 1.19205;   Fc[7]= 0.927771;   Fll[7]= 0.959162; 
         cafac[7]= 0.87399;  
         winpretag[7]= 31923.6 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27900.9;   wouttag[7]= 4803.9;  
       fbb[7]= 0.0874656;   fcc[7]= 0.143919;   fc[7]= 0.129095;   fll[7]= 0.639521;  
       fmcbb[7]= 0.0733739;   fmccc[7]= 0.120732;   fmcc[7]= 0.139145;   fmcll[7]= 0.666749;  
 }  
 else if ( idsys==54   || sysname=="met_down" ) { 
      Fbb[1]= 1.23869;   Fcc[1]= 1.23869;   Fc[1]= 0.962953;   Fll[1]= 0.992461; 
         cafac[1]= 1.0006;  
         winpretag[1]= 2.23889e+06 ;   wintag[1]= 79744.1 ;   woutpretag[1]= 2.24024e+06;   wouttag[1]= 80718.7;  
       fbb[1]= 0.0147818;   fcc[1]= 0.042831;   fc[1]= 0.127721;   fll[1]= 0.814666;  
       fmcbb[1]= 0.0119334;   fmccc[1]= 0.0345776;   fmcc[1]= 0.132635;   fmcll[1]= 0.820854;  
      Fbb[2]= 1.22307;   Fcc[2]= 1.22307;   Fc[2]= 0.950804;   Fll[2]= 0.97994; 
         cafac[2]= 0.919406;  
         winpretag[2]= 514508 ;   wintag[2]= 38652 ;   woutpretag[2]= 473042;   wouttag[2]= 37128.8;  
       fbb[2]= 0.0362422;   fcc[2]= 0.087038;   fc[2]= 0.145097;   fll[2]= 0.731623;  
       fmcbb[2]= 0.0296323;   fmccc[2]= 0.0711638;   fmcc[2]= 0.152605;   fmcll[2]= 0.746599;  
      Fbb[3]= 1.20985;   Fcc[3]= 1.20985;   Fc[3]= 0.940532;   Fll[3]= 0.969353; 
         cafac[3]= 0.832086;  
         winpretag[3]= 114946 ;   wintag[3]= 13251.5 ;   woutpretag[3]= 95645.2;   wouttag[3]= 11685;  
       fbb[3]= 0.05838;   fcc[3]= 0.117758;   fc[3]= 0.142497;   fll[3]= 0.681365;  
       fmcbb[3]= 0.0482539;   fmccc[3]= 0.0973324;   fmcc[3]= 0.151507;   fmcll[3]= 0.702907;  
      Fbb[4]= 1.19756;   Fcc[4]= 1.19756;   Fc[4]= 0.930975;   Fll[4]= 0.959504; 
         cafac[4]= 0.886989;  
         winpretag[4]= 25014.7 ;   wintag[4]= 3751.64 ;   woutpretag[4]= 22187.7;   wouttag[4]= 3568.09;  
       fbb[4]= 0.0833457;   fcc[4]= 0.14081;   fc[4]= 0.132557;   fll[4]= 0.643288;  
       fmcbb[4]= 0.0695964;   fmccc[4]= 0.117581;   fmcc[4]= 0.142385;   fmcll[4]= 0.670438;  
      Fbb[5]= 1.18773;   Fcc[5]= 1.18773;   Fc[5]= 0.923338;   Fll[5]= 0.951632; 
         cafac[5]= 0.825273;  
         winpretag[5]= 6879.09 ;   wintag[5]= 1388.24 ;   woutpretag[5]= 5677.13;   wouttag[5]= 1221.02;  
       fbb[5]= 0.105061;   fcc[5]= 0.156324;   fc[5]= 0.117186;   fll[5]= 0.62143;  
       fmcbb[5]= 0.088455;   fmccc[5]= 0.131615;   fmcc[5]= 0.126916;   fmcll[5]= 0.653014;  
      Fbb[6]= 1.20669;   Fcc[6]= 1.20669;   Fc[6]= 0.938073;   Fll[6]= 0.966819; 
         cafac[6]= 0.841609;  
         winpretag[6]= 146840 ;   wintag[6]= 18391.4 ;   woutpretag[6]= 123582;   wouttag[6]= 16464.2;  
       fbb[6]= 0.0648872;   fcc[6]= 0.12355;   fc[6]= 0.139586;   fll[6]= 0.671976;  
       fmcbb[6]= 0.0537729;   fmccc[6]= 0.102388;   fmcc[6]= 0.148801;   fmcll[6]= 0.695038;  
      Fbb[7]= 1.19543;   Fcc[7]= 1.19543;   Fc[7]= 0.929317;   Fll[7]= 0.957795; 
         cafac[7]= 0.873075;  
         winpretag[7]= 31893.7 ;   wintag[7]= 5139.89 ;   woutpretag[7]= 27845.6;   wouttag[7]= 4805.97;  
       fbb[7]= 0.0880598;   fcc[7]= 0.144178;   fc[7]= 0.12922;   fll[7]= 0.638543;  
       fmcbb[7]= 0.0736639;   fmccc[7]= 0.120608;   fmcc[7]= 0.139048;   fmcll[7]= 0.66668;  
 }  
 else if ( idsys==98   || sysname=="ifsr_nom" ) { 
      Fbb[1]= 1.18973;   Fcc[1]= 1.18973;   Fc[1]= 0.970615;   Fll[1]= 0.993995; 
         cafac[1]= 0.999281;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.24247e+06;   wouttag[1]= 80711.4;  
       fbb[1]= 0.0141871;   fcc[1]= 0.0411591;   fc[1]= 0.12872;   fll[1]= 0.815934;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.17779;   Fcc[2]= 1.17779;   Fc[2]= 0.960872;   Fll[2]= 0.984017; 
         cafac[2]= 0.920651;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 474173;   wouttag[2]= 36929.9;  
       fbb[2]= 0.0348639;   fcc[2]= 0.083774;   fc[2]= 0.146763;   fll[2]= 0.734599;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.16761;   Fcc[3]= 1.16761;   Fc[3]= 0.952567;   Fll[3]= 0.975512; 
         cafac[3]= 0.836232;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 96223.1;   wouttag[3]= 11617.9;  
       fbb[3]= 0.0563763;   fcc[3]= 0.113583;   fc[3]= 0.144205;   fll[3]= 0.685836;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.15818;   Fcc[4]= 1.15818;   Fc[4]= 0.944875;   Fll[4]= 0.967634; 
         cafac[4]= 0.895133;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22406.9;   wouttag[4]= 3553.81;  
       fbb[4]= 0.0802954;   fcc[4]= 0.136147;   fc[4]= 0.134647;   fll[4]= 0.64891;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.15043;   Fcc[5]= 1.15043;   Fc[5]= 0.938552;   Fll[5]= 0.961159; 
         cafac[5]= 0.838747;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5779.3;   wouttag[5]= 1225.22;  
       fbb[5]= 0.101331;   fcc[5]= 0.152204;   fc[5]= 0.119171;   fll[5]= 0.627294;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.16517;   Fcc[6]= 1.16517;   Fc[6]= 0.950584;   Fll[6]= 0.973481; 
         cafac[6]= 0.847023;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 124504;   wouttag[6]= 16379.6;  
       fbb[6]= 0.0626086;   fcc[6]= 0.119282;   fc[6]= 0.141379;   fll[6]= 0.67673;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.1565;   Fcc[7]= 1.1565;   Fc[7]= 0.943503;   Fll[7]= 0.966229; 
         cafac[7]= 0.882425;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 28169.1;   wouttag[7]= 4794.1;  
       fbb[7]= 0.0848599;   fcc[7]= 0.139631;   fc[7]= 0.131289;   fll[7]= 0.64422;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==47   || sysname=="ifsr_up" ) { 
      Fbb[1]= 1.2041;   Fcc[1]= 1.2041;   Fc[1]= 0.948881;   Fll[1]= 0.996692; 
         cafac[1]= 0.995469;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.23392e+06;   wouttag[1]= 79641.6;  
       fbb[1]= 0.0143584;   fcc[1]= 0.0416562;   fc[1]= 0.125838;   fll[1]= 0.818148;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.19184;   Fcc[2]= 1.19184;   Fc[2]= 0.939224;   Fll[2]= 0.986549; 
         cafac[2]= 0.916403;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 471985;   wouttag[2]= 36623.3;  
       fbb[2]= 0.03528;   fcc[2]= 0.0847739;   fc[2]= 0.143456;   fll[2]= 0.73649;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.1809;   Fcc[3]= 1.1809;   Fc[3]= 0.9306;   Fll[3]= 0.97749; 
         cafac[3]= 0.831397;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95666.7;   wouttag[3]= 11539.9;  
       fbb[3]= 0.0570181;   fcc[3]= 0.114876;   fc[3]= 0.14088;   fll[3]= 0.687226;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.17057;   Fcc[4]= 1.17057;   Fc[4]= 0.922462;   Fll[4]= 0.968942; 
         cafac[4]= 0.889244;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22259.5;   wouttag[4]= 3532.87;  
       fbb[4]= 0.0811548;   fcc[4]= 0.137604;   fc[4]= 0.131453;   fll[4]= 0.649788;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.16189;   Fcc[5]= 1.16189;   Fc[5]= 0.915617;   Fll[5]= 0.961752; 
         cafac[5]= 0.831365;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5728.44;   wouttag[5]= 1216.14;  
       fbb[5]= 0.102341;   fcc[5]= 0.15372;   fc[5]= 0.116259;   fll[5]= 0.627681;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.17823;   Fcc[6]= 1.17823;   Fc[6]= 0.928493;   Fll[6]= 0.975277; 
         cafac[6]= 0.841844;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123743;   wouttag[6]= 16273.7;  
       fbb[6]= 0.0633098;   fcc[6]= 0.120618;   fc[6]= 0.138094;   fll[6]= 0.677978;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.16869;   Fcc[7]= 1.16869;   Fc[7]= 0.920976;   Fll[7]= 0.967381; 
         cafac[7]= 0.876195;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27970.2;   wouttag[7]= 4764.63;  
       fbb[7]= 0.0857545;   fcc[7]= 0.141103;   fc[7]= 0.128154;   fll[7]= 0.644988;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==98   || sysname=="ifsr_nom" ) { 
      Fbb[1]= 1.2484;   Fcc[1]= 1.2484;   Fc[1]= 0.994767;   Fll[1]= 0.986768; 
         cafac[1]= 1.00356;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.25207e+06;   wouttag[1]= 82689.9;  
       fbb[1]= 0.0148867;   fcc[1]= 0.0431888;   fc[1]= 0.131923;   fll[1]= 0.810002;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.23074;   Fcc[2]= 1.23074;   Fc[2]= 0.9807;   Fll[2]= 0.972815; 
         cafac[2]= 0.925304;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 476569;   wouttag[2]= 37901.1;  
       fbb[2]= 0.0364315;   fcc[2]= 0.0875408;   fc[2]= 0.149791;   fll[2]= 0.726236;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.21669;   Fcc[3]= 1.21669;   Fc[3]= 0.969499;   Fll[3]= 0.961704; 
         cafac[3]= 0.838181;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 96447.4;   wouttag[3]= 11890.5;  
       fbb[3]= 0.0587461;   fcc[3]= 0.118357;   fc[3]= 0.146769;   fll[3]= 0.676128;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.20408;   Fcc[4]= 1.20408;   Fc[4]= 0.959457;   Fll[4]= 0.951743; 
         cafac[4]= 0.893713;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22371.4;   wouttag[4]= 3626.06;  
       fbb[4]= 0.0834781;   fcc[4]= 0.141544;   fc[4]= 0.136725;   fll[4]= 0.638253;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.19413;   Fcc[5]= 1.19413;   Fc[5]= 0.951527;   Fll[5]= 0.943876; 
         cafac[5]= 0.824172;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5678.87;   wouttag[5]= 1226.5;  
       fbb[5]= 0.105181;   fcc[5]= 0.157986;   fc[5]= 0.120818;   fll[5]= 0.616015;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.21345;   Fcc[6]= 1.21345;   Fc[6]= 0.96692;   Fll[6]= 0.959145; 
         cafac[6]= 0.847425;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 124563;   wouttag[6]= 16736.1;  
       fbb[6]= 0.0652025;   fcc[6]= 0.124224;   fc[6]= 0.143809;   fll[6]= 0.666764;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.20192;   Fcc[7]= 1.20192;   Fc[7]= 0.957734;   Fll[7]= 0.950034; 
         cafac[7]= 0.878025;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 28028.6;   wouttag[7]= 4871.17;  
       fbb[7]= 0.0881932;   fcc[7]= 0.145116;   fc[7]= 0.133269;   fll[7]= 0.633422;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==48   || sysname=="ifsr_down" ) { 
      Fbb[1]= 1.26303;   Fcc[1]= 1.26303;   Fc[1]= 0.973166;   Fll[1]= 0.989429; 
         cafac[1]= 0.999737;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.24349e+06;   wouttag[1]= 81618.3;  
       fbb[1]= 0.0150612;   fcc[1]= 0.043695;   fc[1]= 0.129058;   fll[1]= 0.812186;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.24497;   Fcc[2]= 1.24497;   Fc[2]= 0.959253;   Fll[2]= 0.975283; 
         cafac[2]= 0.921056;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 474381;   wouttag[2]= 37594.5;  
       fbb[2]= 0.0368527;   fcc[2]= 0.088553;   fc[2]= 0.146516;   fll[2]= 0.728079;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.23007;   Fcc[3]= 1.23007;   Fc[3]= 0.947773;   Fll[3]= 0.963611; 
         cafac[3]= 0.833362;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95892.9;   wouttag[3]= 11812.4;  
       fbb[3]= 0.0593924;   fcc[3]= 0.11966;   fc[3]= 0.143479;   fll[3]= 0.677469;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.21651;   Fcc[4]= 1.21651;   Fc[4]= 0.937321;   Fll[4]= 0.952984; 
         cafac[4]= 0.887861;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22224.9;   wouttag[4]= 3604.98;  
       fbb[4]= 0.0843395;   fcc[4]= 0.143004;   fc[4]= 0.133571;   fll[4]= 0.639086;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.20557;   Fcc[5]= 1.20557;   Fc[5]= 0.928894;   Fll[5]= 0.944417; 
         cafac[5]= 0.81689;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5628.7;   wouttag[5]= 1217.37;  
       fbb[5]= 0.106188;   fcc[5]= 0.1595;   fc[5]= 0.117945;   fll[5]= 0.616367;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.22658;   Fcc[6]= 1.22658;   Fc[6]= 0.945078;   Fll[6]= 0.960871; 
         cafac[6]= 0.842271;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123805;   wouttag[6]= 16629.9;  
       fbb[6]= 0.0659079;   fcc[6]= 0.125568;   fc[6]= 0.14056;   fll[6]= 0.667964;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.21413;   Fcc[7]= 1.21413;   Fc[7]= 0.935489;   Fll[7]= 0.951122; 
         cafac[7]= 0.871847;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27831.4;   wouttag[7]= 4841.54;  
       fbb[7]= 0.089089;   fcc[7]= 0.14659;   fc[7]= 0.130174;   fll[7]= 0.634147;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==72   || sysname=="iqopt3_one" ) { 
      Fbb[1]= 1.23347;   Fcc[1]= 1.23347;   Fc[1]= 0.96116;   Fll[1]= 0.993073; 
         cafac[1]= 1.00267;  
         winpretag[1]= 2.23319e+06 ;   wintag[1]= 79525.9 ;   woutpretag[1]= 2.23916e+06;   wouttag[1]= 80530.5;  
       fbb[1]= 0.0146766;   fcc[1]= 0.0425806;   fc[1]= 0.127451;   fll[1]= 0.815292;  
       fmcbb[1]= 0.0118986;   fmccc[1]= 0.0345211;   fmcc[1]= 0.132601;   fmcll[1]= 0.820979;  
      Fbb[2]= 1.2183;   Fcc[2]= 1.2183;   Fc[2]= 0.949344;   Fll[2]= 0.980865; 
         cafac[2]= 0.917148;  
         winpretag[2]= 515607 ;   wintag[2]= 38799.3 ;   woutpretag[2]= 472888;   wouttag[2]= 37108.9;  
       fbb[2]= 0.0361169;   fcc[2]= 0.0867566;   fc[2]= 0.144917;   fll[2]= 0.73221;  
       fmcbb[2]= 0.0296453;   fmccc[2]= 0.071211;   fmcc[2]= 0.15265;   fmcll[2]= 0.746494;  
      Fbb[3]= 1.20531;   Fcc[3]= 1.20531;   Fc[3]= 0.939224;   Fll[3]= 0.970408; 
         cafac[3]= 0.82596;  
         winpretag[3]= 115670 ;   wintag[3]= 13349.1 ;   woutpretag[3]= 95539.1;   wouttag[3]= 11666.4;  
       fbb[3]= 0.0584206;   fcc[3]= 0.117556;   fc[3]= 0.141704;   fll[3]= 0.682319;  
       fmcbb[3]= 0.0484691;   fmccc[3]= 0.0975316;   fmcc[3]= 0.150873;   fmcll[3]= 0.703126;  
      Fbb[4]= 1.19329;   Fcc[4]= 1.19329;   Fc[4]= 0.929855;   Fll[4]= 0.960728; 
         cafac[4]= 0.885964;  
         winpretag[4]= 24982.7 ;   wintag[4]= 3756.81 ;   woutpretag[4]= 22133.8;   wouttag[4]= 3562.89;  
       fbb[4]= 0.083071;   fcc[4]= 0.140863;   fc[4]= 0.131641;   fll[4]= 0.644425;  
       fmcbb[4]= 0.0696151;   fmccc[4]= 0.118046;   fmcc[4]= 0.141572;   fmcll[4]= 0.670767;  
      Fbb[5]= 1.1829;   Fcc[5]= 1.1829;   Fc[5]= 0.921757;   Fll[5]= 0.952361; 
         cafac[5]= 0.731815;  
         winpretag[5]= 7534.4 ;   wintag[5]= 1507.79 ;   woutpretag[5]= 5513.78;   wouttag[5]= 1176.53;  
       fbb[5]= 0.10549;   fcc[5]= 0.158393;   fc[5]= 0.11414;   fll[5]= 0.621977;  
       fmcbb[5]= 0.0891792;   fmccc[5]= 0.133902;   fmcc[5]= 0.123829;   fmcll[5]= 0.65309;  
      Fbb[6]= 1.20211;   Fcc[6]= 1.20211;   Fc[6]= 0.93673;   Fll[6]= 0.967832; 
         cafac[6]= 0.831066;  
         winpretag[6]= 148187 ;   wintag[6]= 18613.7 ;   woutpretag[6]= 123154;   wouttag[6]= 16431.9;  
       fbb[6]= 0.0650392;   fcc[6]= 0.123625;   fc[6]= 0.138571;   fll[6]= 0.672766;  
       fmcbb[6]= 0.054104;   fmccc[6]= 0.102839;   fmcc[6]= 0.14793;   fmcll[6]= 0.695127;  
      Fbb[7]= 1.19087;   Fcc[7]= 1.19087;   Fc[7]= 0.927966;   Fll[7]= 0.958776; 
         cafac[7]= 0.848035;  
         winpretag[7]= 32517.1 ;   wintag[7]= 5264.6 ;   woutpretag[7]= 27575.6;   wouttag[7]= 4775.99;  
       fbb[7]= 0.0883006;   fcc[7]= 0.144952;   fc[7]= 0.127559;   fll[7]= 0.639189;  
       fmcbb[7]= 0.0741482;   fmccc[7]= 0.12172;   fmcc[7]= 0.137461;   fmcll[7]= 0.666671;  
 }  
 else if ( idsys==71   || sysname=="ptjmin_one" ) { 
      Fbb[1]= 1.22979;   Fcc[1]= 1.22979;   Fc[1]= 0.960002;   Fll[1]= 0.993439; 
         cafac[1]= 0.997432;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.23832e+06;   wouttag[1]= 80532.3;  
       fbb[1]= 0.0146648;   fcc[1]= 0.0425451;   fc[1]= 0.127312;   fll[1]= 0.815478;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.21493;   Fcc[2]= 1.21493;   Fc[2]= 0.948402;   Fll[2]= 0.981436; 
         cafac[2]= 0.91819;  
         winpretag[2]= 514785 ;   wintag[2]= 38798.8 ;   woutpretag[2]= 472670;   wouttag[2]= 37108.9;  
       fbb[2]= 0.0361195;   fcc[2]= 0.0867003;   fc[2]= 0.144718;   fll[2]= 0.732462;  
       fmcbb[2]= 0.0297296;   fmccc[2]= 0.0713622;   fmcc[2]= 0.152592;   fmcll[2]= 0.746317;  
      Fbb[3]= 1.20233;   Fcc[3]= 1.20233;   Fc[3]= 0.938562;   Fll[3]= 0.971253; 
         cafac[3]= 0.83322;  
         winpretag[3]= 114807 ;   wintag[3]= 13242 ;   woutpretag[3]= 95659.5;   wouttag[3]= 11659.8;  
       fbb[3]= 0.058165;   fcc[3]= 0.117124;   fc[3]= 0.141864;   fll[3]= 0.682847;  
       fmcbb[3]= 0.0483771;   fmccc[3]= 0.0974145;   fmcc[3]= 0.151151;   fmcll[3]= 0.703058;  
      Fbb[4]= 1.1905;   Fcc[4]= 1.1905;   Fc[4]= 0.929327;   Fll[4]= 0.961697; 
         cafac[4]= 0.886728;  
         winpretag[4]= 24966.8 ;   wintag[4]= 3754.3 ;   woutpretag[4]= 22138.8;   wouttag[4]= 3559.17;  
       fbb[4]= 0.0827238;   fcc[4]= 0.140481;   fc[4]= 0.131897;   fll[4]= 0.644899;  
       fmcbb[4]= 0.0694868;   fmccc[4]= 0.118002;   fmcc[4]= 0.141927;   fmcll[4]= 0.670585;  
      Fbb[5]= 1.18108;   Fcc[5]= 1.18108;   Fc[5]= 0.921978;   Fll[5]= 0.954091; 
         cafac[5]= 0.827463;  
         winpretag[5]= 6795.16 ;   wintag[5]= 1367.63 ;   woutpretag[5]= 5622.74;   wouttag[5]= 1203.57;  
       fbb[5]= 0.104451;   fcc[5]= 0.155665;   fc[5]= 0.117209;   fll[5]= 0.622676;  
       fmcbb[5]= 0.0884366;   fmccc[5]= 0.131798;   fmcc[5]= 0.127127;   fmcll[5]= 0.652638;  
      Fbb[6]= 1.1993;   Fcc[6]= 1.1993;   Fc[6]= 0.936196;   Fll[6]= 0.968805; 
         cafac[6]= 0.842562;  
         winpretag[6]= 146569 ;   wintag[6]= 18363.9 ;   woutpretag[6]= 123493;   wouttag[6]= 16412.2;  
       fbb[6]= 0.0645583;   fcc[6]= 0.122946;   fc[6]= 0.138993;   fll[6]= 0.673502;  
       fmcbb[6]= 0.0538301;   fmccc[6]= 0.102515;   fmcc[6]= 0.148466;   fmcll[6]= 0.695189;  
      Fbb[7]= 1.18847;   Fcc[7]= 1.18847;   Fc[7]= 0.927745;   Fll[7]= 0.960059; 
         cafac[7]= 0.873415;  
         winpretag[7]= 31762 ;   wintag[7]= 5121.93 ;   woutpretag[7]= 27741.4;   wouttag[7]= 4778.12;  
       fbb[7]= 0.0874012;   fcc[7]= 0.14375;   fc[7]= 0.128735;   fll[7]= 0.640115;  
       fmcbb[7]= 0.0735409;   fmccc[7]= 0.120953;   fmcc[7]= 0.138761;   fmcll[7]= 0.666745;  
 }  
 else if ( idsys==69   || sysname=="powhe_one" ) { 
      Fbb[1]= 1.26145;   Fcc[1]= 1.26145;   Fc[1]= 0.96057;   Fll[1]= 0.991553; 
         cafac[1]= 0.997415;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.23828e+06;   wouttag[1]= 80883.5;  
       fbb[1]= 0.0150424;   fcc[1]= 0.0436405;   fc[1]= 0.127388;   fll[1]= 0.81393;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.24403;   Fcc[2]= 1.24403;   Fc[2]= 0.9473;   Fll[2]= 0.977856; 
         cafac[2]= 0.91825;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 472936;   wouttag[2]= 37318.2;  
       fbb[2]= 0.0368247;   fcc[2]= 0.0884857;   fc[2]= 0.14469;   fll[2]= 0.73;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.22931;   Fcc[3]= 1.22931;   Fc[3]= 0.936091;   Fll[3]= 0.966285; 
         cafac[3]= 0.83189;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 95723.5;   wouttag[3]= 11754.7;  
       fbb[3]= 0.0593554;   fcc[3]= 0.119585;   fc[3]= 0.141711;   fll[3]= 0.679349;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.21577;   Fcc[4]= 1.21577;   Fc[4]= 0.925781;   Fll[4]= 0.955643; 
         cafac[4]= 0.887966;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22227.5;   wouttag[4]= 3596.24;  
       fbb[4]= 0.0842881;   fcc[4]= 0.142917;   fc[4]= 0.131926;   fll[4]= 0.640869;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.20471;   Fcc[5]= 1.20471;   Fc[5]= 0.917361;   Fll[5]= 0.946951; 
         cafac[5]= 0.823272;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5672.68;   wouttag[5]= 1224.62;  
       fbb[5]= 0.106113;   fcc[5]= 0.159386;   fc[5]= 0.11648;   fll[5]= 0.618021;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.22581;   Fcc[6]= 1.22581;   Fc[6]= 0.933427;   Fll[6]= 0.963535; 
         cafac[6]= 0.841507;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 123693;   wouttag[6]= 16566.6;  
       fbb[6]= 0.0658667;   fcc[6]= 0.12549;   fc[6]= 0.138827;   fll[6]= 0.669816;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.21336;   Fcc[7]= 1.21336;   Fc[7]= 0.923951;   Fll[7]= 0.953753; 
         cafac[7]= 0.873363;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 27879.8;   wouttag[7]= 4838.57;  
       fbb[7]= 0.0890328;   fcc[7]= 0.146497;   fc[7]= 0.128568;   fll[7]= 0.635902;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==70   || sysname=="powpy_one" ) { 
      Fbb[1]= 1.24788;   Fcc[1]= 1.24788;   Fc[1]= 0.96507;   Fll[1]= 0.991595; 
         cafac[1]= 0.99717;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80031.8 ;   woutpretag[1]= 2.23773e+06;   wouttag[1]= 80911.6;  
       fbb[1]= 0.0148805;   fcc[1]= 0.0431709;   fc[1]= 0.127984;   fll[1]= 0.813964;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.23143;   Fcc[2]= 1.23143;   Fc[2]= 0.952348;   Fll[2]= 0.978523; 
         cafac[2]= 0.918957;  
         winpretag[2]= 515040 ;   wintag[2]= 38733.9 ;   woutpretag[2]= 473300;   wouttag[2]= 37286.1;  
       fbb[2]= 0.0364518;   fcc[2]= 0.0875895;   fc[2]= 0.145461;   fll[2]= 0.730498;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.21758;   Fcc[3]= 1.21758;   Fc[3]= 0.941637;   Fll[3]= 0.967519; 
         cafac[3]= 0.843527;  
         winpretag[3]= 115067 ;   wintag[3]= 13258.7 ;   woutpretag[3]= 97062.5;   wouttag[3]= 11888.2;  
       fbb[3]= 0.0587892;   fcc[3]= 0.118444;   fc[3]= 0.142551;   fll[3]= 0.680216;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.20485;   Fcc[4]= 1.20485;   Fc[4]= 0.931795;   Fll[4]= 0.957406; 
         cafac[4]= 0.912308;  
         winpretag[4]= 25032 ;   wintag[4]= 3753.02 ;   woutpretag[4]= 22836.9;   wouttag[4]= 3682.75;  
       fbb[4]= 0.0835314;   fcc[4]= 0.141634;   fc[4]= 0.132783;   fll[4]= 0.642051;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.19448;   Fcc[5]= 1.19448;   Fc[5]= 0.92377;   Fll[5]= 0.94916; 
         cafac[5]= 0.866466;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1387.25 ;   woutpretag[5]= 5970.3;   wouttag[5]= 1284.92;  
       fbb[5]= 0.105211;   fcc[5]= 0.158032;   fc[5]= 0.117294;   fll[5]= 0.619463;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.21429;   Fcc[6]= 1.21429;   Fc[6]= 0.939097;   Fll[6]= 0.964908; 
         cafac[6]= 0.857169;  
         winpretag[6]= 146990 ;   wintag[6]= 18399 ;   woutpretag[6]= 125995;   wouttag[6]= 16827.4;  
       fbb[6]= 0.0652479;   fcc[6]= 0.124311;   fc[6]= 0.139671;   fll[6]= 0.670771;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.2026;   Fcc[7]= 1.2026;   Fc[7]= 0.930051;   Fll[7]= 0.955614; 
         cafac[7]= 0.901943;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5140.26 ;   woutpretag[7]= 28792.1;   wouttag[7]= 4980.71;  
       fbb[7]= 0.0882428;   fcc[7]= 0.145198;   fc[7]= 0.129417;   fll[7]= 0.637143;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==15   || sysname=="btagb_up" ) { 
      Fbb[1]= 1.13993;   Fcc[1]= 1.13993;   Fc[1]= 0.933567;   Fll[1]= 1.0028; 
         cafac[1]= 0.993122;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 80499.3 ;   woutpretag[1]= 2.22865e+06;   wouttag[1]= 78662.9;  
       fbb[1]= 0.0135933;   fcc[1]= 0.0394364;   fc[1]= 0.123807;   fll[1]= 0.823164;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.13309;   Fcc[2]= 1.13309;   Fc[2]= 0.927961;   Fll[2]= 0.996782; 
         cafac[2]= 0.914768;  
         winpretag[2]= 515040 ;   wintag[2]= 39163.6 ;   woutpretag[2]= 471142;   wouttag[2]= 36273.1;  
       fbb[2]= 0.0335407;   fcc[2]= 0.0805947;   fc[2]= 0.141736;   fll[2]= 0.744128;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.1261;   Fcc[3]= 1.1261;   Fc[3]= 0.92224;   Fll[3]= 0.990636; 
         cafac[3]= 0.829061;  
         winpretag[3]= 115067 ;   wintag[3]= 13441.8 ;   woutpretag[3]= 95397.9;   wouttag[3]= 11435.9;  
       fbb[3]= 0.0543722;   fcc[3]= 0.109545;   fc[3]= 0.139614;   fll[3]= 0.696468;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.11916;   Fcc[4]= 1.11916;   Fc[4]= 0.916553;   Fll[4]= 0.984527; 
         cafac[4]= 0.8854;  
         winpretag[4]= 25032 ;   wintag[4]= 3813.38 ;   woutpretag[4]= 22163.3;   wouttag[4]= 3498.11;  
       fbb[4]= 0.0775901;   fcc[4]= 0.13156;   fc[4]= 0.130611;   fll[4]= 0.660239;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.11296;   Fcc[5]= 1.11296;   Fc[5]= 0.911479;   Fll[5]= 0.979077; 
         cafac[5]= 0.822526;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1409.9 ;   woutpretag[5]= 5667.54;   wouttag[5]= 1200.32;  
       fbb[5]= 0.0980312;   fcc[5]= 0.147247;   fc[5]= 0.115733;   fll[5]= 0.638988;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.12429;   Fcc[6]= 1.12429;   Fc[6]= 0.920757;   Fll[6]= 0.989043; 
         cafac[6]= 0.838939;  
         winpretag[6]= 146990 ;   wintag[6]= 18665.1 ;   woutpretag[6]= 123315;   wouttag[6]= 16120.1;  
       fbb[6]= 0.0604117;   fcc[6]= 0.115097;   fc[6]= 0.136943;   fll[6]= 0.687548;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.11781;   Fcc[7]= 1.11781;   Fc[7]= 0.915453;   Fll[7]= 0.983345; 
         cafac[7]= 0.871257;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5223.29 ;   woutpretag[7]= 27812.6;   wouttag[7]= 4715.08;  
       fbb[7]= 0.0820215;   fcc[7]= 0.134961;   fc[7]= 0.127386;   fll[7]= 0.655632;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==16   || sysname=="btagb_down" ) { 
      Fbb[1]= 1.33731;   Fcc[1]= 1.33731;   Fc[1]= 0.989835;   Fll[1]= 0.982526; 
         cafac[1]= 1.00232;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 79564.4 ;   woutpretag[1]= 2.24929e+06;   wouttag[1]= 82723;  
       fbb[1]= 0.015947;   fcc[1]= 0.0462649;   fc[1]= 0.131269;   fll[1]= 0.806519;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.31189;   Fcc[2]= 1.31189;   Fc[2]= 0.971016;   Fll[2]= 0.963847; 
         cafac[2]= 0.92279;  
         winpretag[2]= 515040 ;   wintag[2]= 38299.4 ;   woutpretag[2]= 475274;   wouttag[2]= 37973;  
       fbb[2]= 0.0388335;   fcc[2]= 0.0933126;   fc[2]= 0.148312;   fll[2]= 0.719542;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.29175;   Fcc[3]= 1.29175;   Fc[3]= 0.956107;   Fll[3]= 0.949047; 
         cafac[3]= 0.835738;  
         winpretag[3]= 115067 ;   wintag[3]= 13071.8 ;   woutpretag[3]= 96166.3;   wouttag[3]= 11918.7;  
       fbb[3]= 0.0623702;   fcc[3]= 0.125659;   fc[3]= 0.144741;   fll[3]= 0.66723;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.27379;   Fcc[4]= 1.27379;   Fc[4]= 0.942816;   Fll[4]= 0.935854; 
         cafac[4]= 0.891702;  
         winpretag[4]= 25032 ;   wintag[4]= 3691.16 ;   woutpretag[4]= 22321;   wouttag[4]= 3638.78;  
       fbb[4]= 0.0883106;   fcc[4]= 0.149738;   fc[4]= 0.134354;   fll[4]= 0.627598;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.25966;   Fcc[5]= 1.25966;   Fc[5]= 0.932361;   Fll[5]= 0.925477; 
         cafac[5]= 0.825628;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1363.29 ;   woutpretag[5]= 5688.91;   wouttag[5]= 1231.62;  
       fbb[5]= 0.110953;   fcc[5]= 0.166656;   fc[5]= 0.118385;   fll[5]= 0.604006;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.28712;   Fcc[6]= 1.28712;   Fc[6]= 0.952682;   Fll[6]= 0.945648; 
         cafac[6]= 0.845192;  
         winpretag[6]= 146990 ;   wintag[6]= 18126.2 ;   woutpretag[6]= 124235;   wouttag[6]= 16783.7;  
       fbb[6]= 0.069161;   fcc[6]= 0.131766;   fc[6]= 0.141691;   fll[6]= 0.657382;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.27071;   Fcc[7]= 1.27071;   Fc[7]= 0.940539;   Fll[7]= 0.933595; 
         cafac[7]= 0.876755;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5054.45 ;   woutpretag[7]= 27988.1;   wouttag[7]= 4888.51;  
       fbb[7]= 0.0932408;   fcc[7]= 0.153421;   fc[7]= 0.130876;   fll[7]= 0.622461;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==85   || sysname=="btagc_up" ) { 
      Fbb[1]= 1.16533;   Fcc[1]= 1.16533;   Fc[1]= 0.835306;   Fll[1]= 1.01724; 
         cafac[1]= 0.976144;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 88322.5 ;   woutpretag[1]= 2.19055e+06;   wouttag[1]= 80119;  
       fbb[1]= 0.0138961;   fcc[1]= 0.040315;   fc[1]= 0.110776;   fll[1]= 0.835013;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.16026;   Fcc[2]= 1.16026;   Fc[2]= 0.831675;   Fll[2]= 1.01282; 
         cafac[2]= 0.895943;  
         winpretag[2]= 515040 ;   wintag[2]= 42011.9 ;   woutpretag[2]= 461447;   wouttag[2]= 36857.5;  
       fbb[2]= 0.0343451;   fcc[2]= 0.0825275;   fc[2]= 0.127029;   fll[2]= 0.756098;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.15236;   Fcc[3]= 1.15236;   Fc[3]= 0.826012;   Fll[3]= 1.00592; 
         cafac[3]= 0.812145;  
         winpretag[3]= 115067 ;   wintag[3]= 14256.6 ;   woutpretag[3]= 93451.4;   wouttag[3]= 11634.4;  
       fbb[3]= 0.0556402;   fcc[3]= 0.1121;   fc[3]= 0.125047;   fll[3]= 0.707213;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.14361;   Fcc[4]= 1.14361;   Fc[4]= 0.819741;   Fll[4]= 0.998283; 
         cafac[4]= 0.86776;  
         winpretag[4]= 25032 ;   wintag[4]= 4013.24 ;   woutpretag[4]= 21721.7;   wouttag[4]= 3556.45;  
       fbb[4]= 0.0792857;   fcc[4]= 0.134435;   fc[4]= 0.116815;   fll[4]= 0.669464;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.13494;   Fcc[5]= 1.13494;   Fc[5]= 0.813525;   Fll[5]= 0.990712; 
         cafac[5]= 0.806579;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1469.02 ;   woutpretag[5]= 5557.65;   wouttag[5]= 1215.95;  
       fbb[5]= 0.0999672;   fcc[5]= 0.150155;   fc[5]= 0.103296;   fll[5]= 0.646582;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.15004;   Fcc[6]= 1.15004;   Fc[6]= 0.824345;   Fll[6]= 1.00389; 
         cafac[6]= 0.821907;  
         winpretag[6]= 146990 ;   wintag[6]= 19738.8 ;   woutpretag[6]= 120812;   wouttag[6]= 16393.5;  
       fbb[6]= 0.0617951;   fcc[6]= 0.117732;   fc[6]= 0.122604;   fll[6]= 0.697869;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.14173;   Fcc[7]= 1.14173;   Fc[7]= 0.818392;   Fll[7]= 0.996639; 
         cafac[7]= 0.853989;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5482.26 ;   woutpretag[7]= 27261.3;   wouttag[7]= 4789.38;  
       fbb[7]= 0.0837765;   fcc[7]= 0.137849;   fc[7]= 0.11388;   fll[7]= 0.664495;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==86   || sysname=="btagc_down" ) { 
      Fbb[1]= 1.29965;   Fcc[1]= 1.29965;   Fc[1]= 1.12459;   Fll[1]= 0.962889; 
         cafac[1]= 1.02706;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 71741.2 ;   woutpretag[1]= 2.3048e+06;   wouttag[1]= 81104.8;  
       fbb[1]= 0.0154979;   fcc[1]= 0.044962;   fc[1]= 0.14914;   fll[1]= 0.7904;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.27229;   Fcc[2]= 1.27229;   Fc[2]= 1.10091;   Fll[2]= 0.942614; 
         cafac[2]= 0.949963;  
         winpretag[2]= 515040 ;   wintag[2]= 35445.1 ;   woutpretag[2]= 489269;   wouttag[2]= 37361.3;  
       fbb[2]= 0.0376612;   fcc[2]= 0.0904958;   fc[2]= 0.168152;   fll[2]= 0.703691;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.25402;   Fcc[3]= 1.25402;   Fc[3]= 1.08511;   Fll[3]= 0.929081; 
         cafac[3]= 0.859935;  
         winpretag[3]= 115067 ;   wintag[3]= 12255.1 ;   woutpretag[3]= 98950.5;   wouttag[3]= 11712.1;  
       fbb[3]= 0.0605487;   fcc[3]= 0.121989;   fc[3]= 0.16427;   fll[3]= 0.653192;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.2391;   Fcc[4]= 1.2391;   Fc[4]= 1.0722;   Fll[4]= 0.918027; 
         cafac[4]= 0.916768;  
         winpretag[4]= 25032 ;   wintag[4]= 3490.06 ;   woutpretag[4]= 22948.5;   wouttag[4]= 3578.49;  
       fbb[4]= 0.0859059;   fcc[4]= 0.14566;   fc[4]= 0.152791;   fll[4]= 0.615643;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.22883;   Fcc[5]= 1.22883;   Fc[5]= 1.0633;   Fll[5]= 0.910414; 
         cafac[5]= 0.848133;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1304.42 ;   woutpretag[5]= 5843.98;   wouttag[5]= 1217.01;  
       fbb[5]= 0.108237;   fcc[5]= 0.162576;   fc[5]= 0.135011;   fll[5]= 0.594176;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.25026;   Fcc[6]= 1.25026;   Fc[6]= 1.08185;   Fll[6]= 0.926291; 
         cafac[6]= 0.869504;  
         winpretag[6]= 146990 ;   wintag[6]= 17049.6 ;   woutpretag[6]= 127808;   wouttag[6]= 16500.8;  
       fbb[6]= 0.0671803;   fcc[6]= 0.127992;   fc[6]= 0.160902;   fll[6]= 0.643926;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.23687;   Fcc[7]= 1.23687;   Fc[7]= 1.07026;   Fll[7]= 0.916373; 
         cafac[7]= 0.901255;  
         winpretag[7]= 31922.4 ;   wintag[7]= 4794.48 ;   woutpretag[7]= 28770.2;   wouttag[7]= 4813.29;  
       fbb[7]= 0.0907575;   fcc[7]= 0.149335;   fc[7]= 0.148928;   fll[7]= 0.610979;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==17   || sysname=="bmtag_up" ) { 
      Fbb[1]= 0.955108;   Fcc[1]= 0.955108;   Fc[1]= 0.93716;   Fll[1]= 1.0127; 
         cafac[1]= 0.994431;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 85157.6 ;   woutpretag[1]= 2.23159e+06;   wouttag[1]= 81624.7;  
       fbb[1]= 0.0113893;   fcc[1]= 0.0330423;   fc[1]= 0.124283;   fll[1]= 0.831285;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 0.959562;   Fcc[2]= 0.959562;   Fc[2]= 0.941531;   Fll[2]= 1.01742; 
         cafac[2]= 0.918058;  
         winpretag[2]= 515040 ;   wintag[2]= 41438.4 ;   woutpretag[2]= 472837;   wouttag[2]= 36915.6;  
       fbb[2]= 0.0284042;   fcc[2]= 0.0682521;   fc[2]= 0.143809;   fll[2]= 0.759535;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 0.961958;   Fcc[3]= 0.961958;   Fc[3]= 0.943882;   Fll[3]= 1.01996; 
         cafac[3]= 0.832464;  
         winpretag[3]= 115067 ;   wintag[3]= 14212.9 ;   woutpretag[3]= 95789.5;   wouttag[3]= 11524.9;  
       fbb[3]= 0.0464468;   fcc[3]= 0.0935777;   fc[3]= 0.14289;   fll[3]= 0.717085;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 0.963616;   Fcc[4]= 0.963616;   Fc[4]= 0.945509;   Fll[4]= 1.02172; 
         cafac[4]= 0.88954;  
         winpretag[4]= 25032 ;   wintag[4]= 4013.63 ;   woutpretag[4]= 22266.9;   wouttag[4]= 3480.69;  
       fbb[4]= 0.0668067;   fcc[4]= 0.113276;   fc[4]= 0.134737;   fll[4]= 0.68518;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 0.964352;   Fcc[5]= 0.964352;   Fc[5]= 0.946231;   Fll[5]= 1.0225; 
         cafac[5]= 0.828692;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1517.74 ;   woutpretag[5]= 5710.02;   wouttag[5]= 1233.36;  
       fbb[5]= 0.0849415;   fcc[5]= 0.127586;   fc[5]= 0.120146;   fll[5]= 0.667327;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 0.962352;   Fcc[6]= 0.962352;   Fc[6]= 0.944269;   Fll[6]= 1.02038; 
         cafac[6]= 0.842787;  
         winpretag[6]= 146990 ;   wintag[6]= 19744.3 ;   woutpretag[6]= 123881;   wouttag[6]= 16217.7;  
       fbb[6]= 0.0517103;   fcc[6]= 0.0985188;   fc[6]= 0.14044;   fll[6]= 0.709331;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 0.963775;   Fcc[7]= 0.963775;   Fc[7]= 0.945665;   Fll[7]= 1.02189; 
         cafac[7]= 0.875928;  
         winpretag[7]= 31922.4 ;   wintag[7]= 5531.37 ;   woutpretag[7]= 27961.7;   wouttag[7]= 4730.87;  
       fbb[7]= 0.0707187;   fcc[7]= 0.116363;   fc[7]= 0.13159;   fll[7]= 0.681329;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  
 else if ( idsys==18   || sysname=="bmtag_down" ) { 
      Fbb[1]= 1.50341;   Fcc[1]= 1.50341;   Fc[1]= 0.986526;   Fll[1]= 0.973647; 
         cafac[1]= 1.00111;  
         winpretag[1]= 2.24408e+06 ;   wintag[1]= 74906.1 ;   woutpretag[1]= 2.24657e+06;   wouttag[1]= 79705;  
       fbb[1]= 0.0179276;   fcc[1]= 0.0520111;   fc[1]= 0.13083;   fll[1]= 0.799231;  
       fmcbb[1]= 0.0119246;   fmccc[1]= 0.0345954;   fmcc[1]= 0.132617;   fmcll[1]= 0.820863;  
      Fbb[2]= 1.46107;   Fcc[2]= 1.46107;   Fc[2]= 0.958744;   Fll[2]= 0.946228; 
         cafac[2]= 0.919809;  
         winpretag[2]= 515040 ;   wintag[2]= 36081.5 ;   woutpretag[2]= 473739;   wouttag[2]= 37301.3;  
       fbb[2]= 0.0432495;   fcc[2]= 0.103924;   fc[2]= 0.146438;   fll[2]= 0.706389;  
       fmcbb[2]= 0.0296012;   fmccc[2]= 0.0711284;   fmcc[2]= 0.152739;   fmcll[2]= 0.746531;  
      Fbb[3]= 1.42813;   Fcc[3]= 1.42813;   Fc[3]= 0.93713;   Fll[3]= 0.924896; 
         cafac[3]= 0.832705;  
         winpretag[3]= 115067 ;   wintag[3]= 12317.5 ;   woutpretag[3]= 95817.3;   wouttag[3]= 11800.1;  
       fbb[3]= 0.0689555;   fcc[3]= 0.138927;   fc[3]= 0.141868;   fll[3]= 0.65025;  
       fmcbb[3]= 0.0482836;   fmccc[3]= 0.0972784;   fmcc[3]= 0.151386;   fmcll[3]= 0.703052;  
      Fbb[4]= 1.39919;   Fcc[4]= 1.39919;   Fc[4]= 0.918138;   Fll[4]= 0.906152; 
         cafac[4]= 0.888092;  
         winpretag[4]= 25032 ;   wintag[4]= 3498.79 ;   woutpretag[4]= 22230.7;   wouttag[4]= 3642.79;  
       fbb[4]= 0.0970047;   fcc[4]= 0.164479;   fc[4]= 0.130837;   fll[4]= 0.607679;  
       fmcbb[4]= 0.0693291;   fmccc[4]= 0.117553;   fmcc[4]= 0.142503;   fmcll[4]= 0.670615;  
      Fbb[5]= 1.37671;   Fcc[5]= 1.37671;   Fc[5]= 0.903384;   Fll[5]= 0.891591; 
         cafac[5]= 0.820536;  
         winpretag[5]= 6890.4 ;   wintag[5]= 1258.73 ;   woutpretag[5]= 5653.82;   wouttag[5]= 1199.55;  
       fbb[5]= 0.121262;   fcc[5]= 0.182141;   fc[5]= 0.114705;   fll[5]= 0.581891;  
       fmcbb[5]= 0.0880814;   fmccc[5]= 0.132302;   fmcc[5]= 0.126973;   fmcll[5]= 0.652643;  
      Fbb[6]= 1.42064;   Fcc[6]= 1.42064;   Fc[6]= 0.932213;   Fll[6]= 0.920044; 
         cafac[6]= 0.841801;  
         winpretag[6]= 146990 ;   wintag[6]= 17075 ;   woutpretag[6]= 123736;   wouttag[6]= 16642.1;  
       fbb[6]= 0.0763356;   fcc[6]= 0.145435;   fc[6]= 0.138647;   fll[6]= 0.639582;  
       fmcbb[6]= 0.0537332;   fmccc[6]= 0.102373;   fmcc[6]= 0.148729;   fmcll[6]= 0.695165;  
      Fbb[7]= 1.39428;   Fcc[7]= 1.39428;   Fc[7]= 0.914913;   Fll[7]= 0.902969; 
         cafac[7]= 0.872742;  
         winpretag[7]= 31922.4 ;   wintag[7]= 4757.53 ;   woutpretag[7]= 27860;   wouttag[7]= 4859.4;  
       fbb[7]= 0.102307;   fcc[7]= 0.16834;   fc[7]= 0.127311;   fll[7]= 0.602042;  
       fmcbb[7]= 0.0733768;   fmccc[7]= 0.120737;   fmcc[7]= 0.13915;   fmcll[7]= 0.666736;  
 }  

                else { cout << " WARNING: Ffactors request for unknown variation. Set to Nominal! " << idsys << "  " << sysname <<  endl; }

     return;


}
