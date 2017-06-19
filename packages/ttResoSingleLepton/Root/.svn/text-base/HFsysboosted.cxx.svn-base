/**
 * ########################################## WARNING #################################################################
 * #                                                                                                                  #
 * #  This script is depcrecated and currently only kept for histogrical reasons, HFsys for boosted are implementd    #
 * #  in the HFsys.cxx/.h files aswell, please use those for boosted Wjet SFs !! (deprecate warning added 2014/11/04) #
 * #                                                                                                                  #
 * ########################################## WARNING #################################################################
 */

#include "ttResoSingleLepton/HFsysboosted.h"
#include <iomanip>
#include<iostream>
#include<cmath>
#include<cstdlib>

using namespace std;

// FROM ELECTRON FILE ON TWIKI [iwatson]

void SetWflavors_boosted_elec(int idsys, TString sysname, int ijet, double Wjets_in[5],double Wjets_out[5], double& canorm, int _mode) {
	//-
	// This function corrects the Wjet ligth and heavy flavor components and return the normalization based from Charge Assymmetry.
	//
	// NOTE: this function shoud always be called with PRETAGGED event counts. NEVER WITH TAGGED counts.
	//
	// INPUTS:
	// for idsys there are two input modes: the individul systematics (long list), or the combined systematics (4 checks).
	// Long list:
	// input idsys>0..999 or sysname are passed to the function GetFFactors_boosted_elec (see below) that returns the relevant HF factors.
	// This allows to use the full list of systematics.
	// Short list:
	// most analysers will use the total sum of all systematics. In that case only four checks up and two down are needed.
	// Complete set: idsys=2000 (up), -2000 (down) ---> checks Fbb + Fcc + Fc versus Fl
	//                     2001 (up), -2001 (down) ---> checks Fbb+Fcc versus Fc
	//                     2002 (up), -2002 (down) ---> checks Fll  OBSOLETE, not needed anymore, return nominal
	//                     2003 (up), -2003 (down) ---> check Fbb+xx% in 1,3,4,5 jetbin. --> see remark *
	//                     2004 (up), -2004 (down) ---> check Fc+xx% in 1,3,4,5 jetbin.  --> see remark *
	//                     2005 (up), -2005 (down) ---> check canorm. Please note that canorm is returned and NOT yet applied to Wjets_out.
	// * please note that the xx% per jetbin change should be applied UNCORRELATED to all jetbins. So, for example, when you use simultaneoulsy 
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
	//      please note that you have to do the xx% change per jetbin per flavor, which is NOT available in the long list. For these check you can use idsys=+/-2003,2004.
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
	bool use8TeV=true; //set to true if using 8TeV SFs derived for ttbar resonance resolved selection
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
	  if (use8TeV && idsys==0) GetFFactors8TeV_boosted_elec( idsys, sysname,  Hbb,Hcc, Hc,  Hll, cafac,
							 fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);	
	  else GetFFactors_boosted_elec( idsys, sysname,  Hbb,Hcc, Hc,  Hll, cafac,
				 fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
	} else if (idsys*idsys<1999*1999) {
	  cout << " FATAL  HF KFACTORS probabbly called with old defintion of systematics ids " << endl;
	} else  {                  // go for the reduced errors.
	  if (use8TeV && (abs(idsys)==2000|| abs(idsys)==2003 || abs(idsys)==2004 || abs(idsys)==2005)) 
	    GetFFactors8TeV_boosted_elec( idsys, "",  Hbb,Hcc, Hc,  Hll,cafac, 
				  fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
	  else 	GetFFactors_boosted_elec( 0, "",  Hbb,Hcc, Hc,  Hll,cafac,
				  fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
	}

	//CHECK 8TeV SFs 
//         cout<<"use8TeV "<<use8TeV<<endl;
// 	for (int i=1;i<=7;i++)
// 	cout<<Hbb[i]<<setw(10)<<Hcc[i]<<setw(10)<<Hc[i]<<setw(10)<<Hll[i]<<setw(10)<<cafac[i]<<endl;

	// errors on fraction, not on K factors
	double sigma_hf[2]={ 0.24 ,  -0.25}; // error on ttoal HF
	double sigma_as[2] ={ 0.27 ,  -0.28}; // error on fbb (and fcc) wrt to fc
	double sigma_canorm_up[9]=  {0,  0.088, 0.099, 0.12, 0.12, 0.25, 0.11,0.13,0};
	double sigma_canorm_down[9]={0, -0.083, -0.07,-0.10,-0.12,-0.18, -0.094,-0.11,0};
	double sigma_wbbwcc[9]=  {0,  0.25, 0.0, 0.11, 0.15, 0.29, 0.11,0.15,0};
	double sigma_wc[9]=      {0,  0.25, 0.0, 0.13, 0.15, 0.38, 0.13,0.15,0};

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
	//   cout << " fraction in " << fbb_in <<  " " << fcc_in << " " << fc_in << " " << fll_in << endl;
	double fbb_out=Wjets_out[0]/nw_in;
	double fcc_out=Wjets_out[1]/nw_in;
	double fc_out= Wjets_out[2]/nw_in;
	double fll_out=Wjets_out[3]/nw_in;
	double fhf_out=fbb_out+fcc_out+fc_out;
	if (!use8TeV) { //for 7TeV only
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
		fbb_out*=1+sigma_wbbwcc[ijet];
		fcc_out*=1+sigma_wbbwcc[ijet];
		double rest=fc_out+fll_out;
		double rest_new=1-fbb_out-fcc_out;
		fc_out*=rest_new/rest;
		fll_out*=rest_new/rest;
		if (ijet==1) cname="WbbWccjet1_up";
		if (ijet==3 || ijet==6) cname="WbbWccjet3_up";
		if (ijet==4 || ijet==7) cname="WbbWccjet4_up";
		if (ijet==5) cname="WbbWccjet5_up"; 
		GetFFactors_boosted_elec( -1, cname,  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		canorm=cafac[ijet];
	} else if (ijet!=2 && idsys==-2003) {
		fbb_out*=1-sigma_wbbwcc[ijet];
		fcc_out*=1-sigma_wbbwcc[ijet];
		double rest=fc_out+fll_out;
		double rest_new=1-fbb_out-fcc_out;
		fc_out*=rest_new/rest;
		fll_out*=rest_new/rest;
		if (ijet==1) cname="WbbWccjet1_down";
		if (ijet==3 || ijet==6) cname="WbbWccjet3_down";
		if (ijet==4 || ijet==7) cname="WbbWccjet4_down";
		if (ijet==5) cname="WbbWccjet5_down"; 
		GetFFactors_boosted_elec( -1, cname,  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		canorm=cafac[ijet];
	}
	if (ijet!=2 && idsys==2004) { // change c. at the cost of fbb+fcc+fl
		fc_out*=1+sigma_wc[ijet];
		double rest=fbb_out+fcc_out+fll_out;
		double rest_new=1-fc_out; 
		fbb_out*=rest_new/rest;
		fcc_out*=rest_new/rest;
		fll_out*=rest_new/rest;
		if (ijet==1) cname="Wcjet1_up";
		if (ijet==3 || ijet==6) cname="Wcjet3_up";
		if (ijet==4 || ijet==7) cname="Wcjet4_up";
		if (ijet==5) cname="Wcjet5_up"; // 
		GetFFactors_boosted_elec( -1, cname,  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		canorm=cafac[ijet];
	} else if (ijet!=2 && idsys==-2004) {
		fc_out*=1-sigma_wc[ijet];
		double rest=fbb_out+fcc_out+fll_out;
		double rest_new=1-fc_out; 
		fbb_out*=rest_new/rest;
		fcc_out*=rest_new/rest;
		fll_out*=rest_new/rest;
		if (ijet==1) cname="Wcjet1_down";
		if (ijet==3 || ijet==6) cname="Wcjet3_down";
		if (ijet==4 || ijet==7) cname="Wcjet4_down";
		if (ijet==5) cname="Wcjet5_down"; 
		GetFFactors_boosted_elec( -1, cname,  Hbb,Hcc, Hc,  Hll,cafac,
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
	}// end if (!use8TeV)
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
		GetFFactors_boosted_elec( 0, "",  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		double nwnom=winpretag[ijet];	
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
		GetFFactors_boosted_elec( 0, "",  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		double nwprenom=winpretag[ijet];

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



void GetFFactors8TeV_boosted_elec(int idsys, TString sysname, double Fbb[9],double Fcc[9],double Fc[9], double Fll[9], double cafac[9],
		      double fbb[9],double fcc[9],double fc[9], double fll[9],
		      double winpretag[9],double wintag[9],double woutpretag[9], double wouttag[9]) {
  //8TeV nominal with option1
  Fbb[1]= 1.55801;   Fcc[1]= 1.55801;   Fc[1]= 1.0035;     Fll[1]= 0.967429; cafac[1]= 0.654; // 2013-02-22 
  Fbb[2]= 1.50299;   Fcc[2]= 1.50299;   Fc[2]= 0.968063;   Fll[2]= 0.933262; cafac[2]= 0.654; // 2013-02-22 
  Fbb[3]= 1.45578;   Fcc[3]= 1.45578;   Fc[3]= 0.937655;   Fll[3]= 0.903947; cafac[3]= 0.654; // 2013-02-22 
  Fbb[4]= 1.41531;   Fcc[4]= 1.41531;   Fc[4]= 0.911586;   Fll[4]= 0.878816; cafac[4]= 0.654; // 2013-02-22 
  Fbb[5]= 1.36036;   Fcc[5]= 1.36036;   Fc[5]= 0.876194;   Fll[5]= 0.844696; cafac[5]= 0.654; // 2013-02-22 
  Fbb[6]= 1.44281;   Fcc[6]= 1.44281;   Fc[6]= 0.929302;   Fll[6]= 0.895895; cafac[6]= 0.654; // 2013-02-22 
  Fbb[7]= 1.40189;   Fcc[7]= 1.40189;   Fc[7]= 0.902948;   Fll[7]= 0.870488; cafac[7]= 0.654; // 2013-02-22 
  //var3: BB/CC/C correlated variation,  
 if (idsys==2000) { 
   Fbb[2]= 1.596; Fcc[2]= 1.596; Fc[2]= 1.207;	Fll[2]= 0.795;	cafac[2]= 0.676;
   Fbb[3]= 1.596; Fcc[3]= 1.596; Fc[3]= 1.207;	Fll[3]= 0.795;	cafac[3]= 0.676;
   Fbb[7]= 1.531; Fcc[7]= 1.531; Fc[7]= 1.148;  Fll[7]=	0.751;  cafac[7]= 0.676;
  }
 else if (idsys==-2000) {
   Fbb[2]= 1.271; Fcc[2]= 1.271; Fc[2]= 0.583;  Fll[2]= 1.047; cafac[2]= 0.640;
   Fbb[3]= 1.271; Fcc[3]= 1.271; Fc[3]= 0.583;  Fll[3]= 1.047; cafac[3]= 0.640;
   Fbb[7]= 1.224; Fcc[7]= 1.224; Fc[7]=	0.565;  Fll[7]= 1.035; cafac[7]= 0.640;
 }
 //var1: BB/CC variation
 if (idsys==2003) { 
   Fbb[2]= 1.714; Fcc[2]= 1.714; Fc[2]= 0.886;	Fll[2]= 0.854;  cafac[2]= 0.649;
   Fbb[3]= 1.714; Fcc[3]= 1.714; Fc[3]= 0.886;	Fll[3]= 0.854;  cafac[3]= 0.649;
   Fbb[7]= 1.633; Fcc[7]= 1.633; Fc[7]=	0.831;  Fll[7]= 0.801; 	cafac[7]= 0.649;
  }
 else if (idsys==-2003) {
   Fbb[2]= 1.166; Fcc[2]= 1.166; Fc[2]=	0.996;  Fll[2]=	0.960; cafac[2]= 0.674;
   Fbb[3]= 1.166; Fcc[3]= 1.166; Fc[3]=	0.996;  Fll[3]=	0.960; cafac[3]= 0.674;
   Fbb[7]= 1.127; Fcc[7]= 1.127; Fc[7]=	0.989;  Fll[7]=	0.954; cafac[7]= 0.674;
 }
 //var2: C variation
 if (idsys==2004) { 
   Fbb[2]= 1.350; Fcc[2]= 1.350; Fc[2]=	1.272; Fll[2]= 0.838; cafac[2]=	0.69;
   Fbb[3]= 1.350; Fcc[3]= 1.350; Fc[3]=	1.272; Fll[3]= 0.838; cafac[3]=	0.69;
   Fbb[7]= 1.308; Fcc[7]= 1.308; Fc[7]=	1.242; Fll[7]= 0.812; cafac[7]=	0.69;
  }
 else if (idsys==-2004) {
   Fbb[2]= 1.579; Fcc[2]= 1.579; Fc[2]=	0.546; Fll[2]= 0.981; cafac[2]=	0.63;
   Fbb[3]= 1.579; Fcc[3]= 1.579; Fc[3]=	0.546; Fll[3]= 0.981; cafac[3]=	0.63;
   Fbb[7]= 1.511; Fcc[7]= 1.511; Fc[7]=	0.512; Fll[7]= 0.938; cafac[7]=	0.63;
 }
 //variation on the overall normalisation only
 if (idsys==2005) { 
   cafac[2]*=1+0.1345;  // +-10% (from resolved) +- 9.04% (stat)
   cafac[3]*=1+0.1345;  // +-10% (from resolved) +- 9.04% (stat)
   cafac[7]*=1+0.1345;  // +-10% (from resolved) +- 9.04% (stat)
  }
 else if (idsys==-2005) {
   cafac[2]*=1-0.1345;
   cafac[3]*=1-0.1345;
   cafac[7]*=1-0.1345;
 }
}

void GetFFactors8TeV_boosted_muon(int idsys, TString sysname, double Fbb[9],double Fcc[9],double Fc[9], double Fll[9], double cafac[9],
		      double fbb[9],double fcc[9],double fc[9], double fll[9],
		      double winpretag[9],double wintag[9],double woutpretag[9], double wouttag[9]) {
  //8TeV nominal with option1
  Fbb[1]= 1.73981;   Fcc[1]= 1.73981;   Fc[1]= 0.941812;   Fll[1]= 0.969285; cafac[1]= 0.809; // 2013-02-22
  Fbb[2]= 1.66549;   Fcc[2]= 1.66549;   Fc[2]= 0.901581;   Fll[2]= 0.92788;  cafac[2]= 0.809; // 2013-02-22
  Fbb[3]= 1.6002;    Fcc[3]= 1.6002;    Fc[3]= 0.866238;   Fll[3]= 0.891507; cafac[3]= 0.809; // 2013-02-22
  Fbb[3]= 1.6002;    Fcc[3]= 1.6002;    Fc[3]= 0.866238;   Fll[3]= 0.891507; cafac[3]= 0.809; // 2013-02-22
  Fbb[4]= 1.54469;   Fcc[4]= 1.54469;   Fc[4]= 0.836187;   Fll[4]= 0.860579; cafac[4]= 0.809; // 2013-02-22
  Fbb[5]= 1.46185;   Fcc[5]= 1.46185;   Fc[5]= 0.791343;   Fll[5]= 0.814426; cafac[5]= 0.809; // 2013-02-22
  Fbb[6]= 1.58193;   Fcc[6]= 1.58193;   Fc[6]= 0.856346;   Fll[6]= 0.881326; cafac[6]= 0.809; // 2013-02-22
  Fbb[7]= 1.52438;   Fcc[7]= 1.52438;   Fc[7]= 0.825194;   Fll[7]= 0.849265; cafac[7]= 0.809; // 2013-02-22
  Fbb[7]= 1.52438;   Fcc[7]= 1.52438;   Fc[7]= 0.825194;   Fll[7]= 0.849265; cafac[7]= 0.809; // 2013-02-22
  //var3: BB/CC/C correlated variation,  
 if (idsys==2000) { 
   Fbb[2]=1.764; Fcc[2]=1.764; Fc[2]= 1.063; Fll[2]= 0.801; cafac[2]= 0.845;
   Fbb[3]=1.764; Fcc[3]=1.764; Fc[3]= 1.063; Fll[3]= 0.801; cafac[3]= 0.845;
   Fbb[7]=1.674; Fcc[7]=1.674; Fc[7]= 1.001; Fll[7]= 0.746; cafac[7]= 0.845;
  }
 else if (idsys==-2000) {
   Fbb[2]= 1.393; Fcc[2]= 1.393; Fc[2]=	0.618; Fll[2]= 1.006; cafac[2]=	0.805;
   Fbb[3]= 1.393; Fcc[3]= 1.393; Fc[3]=	0.618; Fll[3]= 1.006; cafac[3]=	0.805;
   Fbb[7]= 1.326; Fcc[7]= 1.326; Fc[7]=	0.593; Fll[7]= 0.985; cafac[7]=	0.805;
 }
 //var1: BB/CC variation
 if (idsys==2003) { 
   Fbb[2]= 1.856; Fcc[2]= 1.856; Fc[2]=	0.818; Fll[2]= 0.843; cafac[2]=	0.818;
   Fbb[3]= 1.856; Fcc[3]= 1.856; Fc[3]=	0.818; Fll[3]= 0.843; cafac[3]=	0.818;
   Fbb[7]= 1.753; Fcc[7]= 1.753; Fc[7]=	0.759; Fll[7]= 0.781; cafac[7]=	0.818;
  }
 else if (idsys==-2003) {
   Fbb[2]= 1.311; Fcc[2]= 1.311; Fc[2]=	0.919; Fll[2]= 0.947; cafac[2]=	0.840;
   Fbb[3]= 1.311; Fcc[3]= 1.311; Fc[3]=	0.919; Fll[3]= 0.947; cafac[3]=	0.840;
   Fbb[7]= 1.253; Fcc[7]= 1.253; Fc[7]=	0.904; Fll[7]= 0.930; cafac[7]=	0.840;
 }
 //var2: C variation
 if (idsys==2004) { 
   Fbb[2]= 1.515; Fcc[2]= 1.515; Fc[2]=	1.121; Fll[2]= 0.845; cafac[2]=	0.858;
   Fbb[3]= 1.515; Fcc[3]= 1.515; Fc[3]=	1.121; Fll[3]= 0.845; cafac[3]=	0.858;
   Fbb[7]= 1.451; Fcc[7]= 1.451; Fc[7]=	1.085; Fll[7]= 0.808; cafac[7]=	0.858;
  }
 else if (idsys==-2004) {
   Fbb[2]= 1.694; Fcc[2]=1.694; Fc[2]= 0.580; Fll[2]= 0.944; cafac[2]= 0.797;
   Fbb[3]= 1.694; Fcc[3]=1.694; Fc[3]= 0.580; Fll[3]= 0.944; cafac[3]= 0.797;
   Fbb[7]= 1.606; Fcc[7]=1.606; Fc[7]= 0.538; Fll[7]= 0.895; cafac[7]= 0.797;
 }
 //variation on the overall normalisation only
 if (idsys==2005) { 
   cafac[2]*=1+0.1154;
   cafac[3]*=1+0.1154;
   cafac[7]*=1+0.1154;
  }
 else if (idsys==-2005) {
   cafac[2]*=1-0.1154; // +-9% (from resolved) +- 7.23% (stat)
   cafac[3]*=1-0.1154; // +-9% (from resolved) +- 7.23% (stat)
   cafac[7]*=1-0.1154; // +-9% (from resolved) +- 7.23% (stat)
 }
}


void GetFFactors_boosted_elec(int idsys, TString sysname, double Fbb[9],double Fcc[9],double Fc[9], double Fll[9], double cafac[9],
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
		//      GetFFactors_boosted_elec(idsys, "" , Fbb,Fcc,Fc,Fll); 
		//      .....
		//--
		// If you want to read this function by 'eye', i have two hints:
		// The HF factors are Fbb[2], Fcc[2] , Fc[2], Fll[2] for EACH jetbin, defined on the PRETAGGED sample.
		// The normalization factors are cafac[1], cafac[2], cafac[3], cafac[4], cafac[5] for jetbins 1,2,3,4,5 resp.
		// Ignore the other quantities which are used for additional checks and studies.
		double fmcbb[9];double fmccc[9];double fmcc[9];double fmcll[9];	


           Fbb[1]= 1.43394;   Fcc[1]= 1.43394;   Fc[1]= 0.740266;   Fll[1]= 1.0196; 
              cafac[1]= 0.980017;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06117e+06;   wouttag[1]= 36033;  
            fbb[1]= 0.0172893;   fcc[1]= 0.0503027;   fc[1]= 0.103702;   fll[1]= 0.828706;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.877209;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17583.5;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38474;   Fcc[3]= 1.38474;   Fc[3]= 0.714865;   Fll[3]= 0.984614; 
              cafac[3]= 0.807815;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46741.7;   wouttag[3]= 5659.07;  
            fbb[3]= 0.0697545;   fcc[3]= 0.135859;   fc[3]= 0.116675;   fll[3]= 0.677712;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35777;   Fcc[4]= 1.35777;   Fc[4]= 0.700942;   Fll[4]= 0.965438; 
              cafac[4]= 0.828149;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10564.8;   wouttag[4]= 1770.22;  
            fbb[4]= 0.0970496;   fcc[4]= 0.160392;   fc[4]= 0.105543;   fll[4]= 0.637016;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32667;   Fcc[5]= 1.32667;   Fc[5]= 0.68489;   Fll[5]= 0.943328; 
              cafac[5]= 0.817465;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2754.74;   wouttag[5]= 592.298;  
            fbb[5]= 0.123928;   fcc[5]= 0.192951;   fc[5]= 0.0924638;   fll[5]= 0.590657;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37727;   Fcc[6]= 1.37727;   Fc[6]= 0.711013;   Fll[6]= 0.979308; 
              cafac[6]= 0.811543;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60045.2;   wouttag[6]= 8024.64;  
            fbb[6]= 0.0770898;   fcc[6]= 0.142849;   fc[6]= 0.113583;   fll[6]= 0.666478;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.35115;   Fcc[7]= 1.35115;   Fc[7]= 0.697526;   Fll[7]= 0.960733; 
              cafac[7]= 0.825146;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13307.2;   wouttag[7]= 2364.1;  
            fbb[7]= 0.10277;   fcc[7]= 0.167321;   fc[7]= 0.102759;   fll[7]= 0.62715;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
       if ( idsys==0   || sysname=="norminal" ) { 
// OK
      }  
      else if ( idsys==1   || sysname=="tt_up" ) { 
           Fbb[1]= 1.43477;   Fcc[1]= 1.43477;   Fc[1]= 0.686119;   Fll[1]= 1.02889; 
              cafac[1]= 0.970342;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.0507e+06;   wouttag[1]= 34553.7;  
            fbb[1]= 0.0172993;   fcc[1]= 0.0503318;   fc[1]= 0.0961167;   fll[1]= 0.836252;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41383;   Fcc[2]= 1.41383;   Fc[2]= 0.676105;   Fll[2]= 1.01387; 
              cafac[2]= 0.867114;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 222948;   wouttag[2]= 17053.6;  
            fbb[2]= 0.0425567;   fcc[2]= 0.102102;   fc[2]= 0.109678;   fll[2]= 0.745663;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38865;   Fcc[3]= 1.38865;   Fc[3]= 0.664066;   Fll[3]= 0.995815; 
              cafac[3]= 0.798319;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46192.3;   wouttag[3]= 5522.7;  
            fbb[3]= 0.0699518;   fcc[3]= 0.136243;   fc[3]= 0.108384;   fll[3]= 0.685421;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.36096;   Fcc[4]= 1.36096;   Fc[4]= 0.650823;   Fll[4]= 0.975957; 
              cafac[4]= 0.819298;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10451.9;   wouttag[4]= 1736.95;  
            fbb[4]= 0.0972779;   fcc[4]= 0.160769;   fc[4]= 0.0979962;   fll[4]= 0.643957;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32903;   Fcc[5]= 1.32903;   Fc[5]= 0.635556;   Fll[5]= 0.953064; 
              cafac[5]= 0.809238;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2727.01;   wouttag[5]= 583.013;  
            fbb[5]= 0.124148;   fcc[5]= 0.193295;   fc[5]= 0.0858035;   fll[5]= 0.596753;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.38098;   Fcc[6]= 1.38098;   Fc[6]= 0.6604;   Fll[6]= 0.990318; 
              cafac[6]= 0.802229;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 59356.1;   wouttag[6]= 7845.46;  
            fbb[6]= 0.0772976;   fcc[6]= 0.143234;   fc[6]= 0.105498;   fll[6]= 0.673971;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.35416;   Fcc[7]= 1.35416;   Fc[7]= 0.647573;   Fll[7]= 0.971083; 
              cafac[7]= 0.816437;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13166.7;   wouttag[7]= 2321.56;  
            fbb[7]= 0.102999;   fcc[7]= 0.167694;   fc[7]= 0.0954003;   fll[7]= 0.633907;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==2   || sysname=="tt_down" ) { 
           Fbb[1]= 1.43307;   Fcc[1]= 1.43307;   Fc[1]= 0.797092;   Fll[1]= 1.00986; 
              cafac[1]= 0.99038;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.07239e+06;   wouttag[1]= 37617.5;  
            fbb[1]= 0.0172788;   fcc[1]= 0.0502721;   fc[1]= 0.111663;   fll[1]= 0.820786;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.40684;   Fcc[2]= 1.40684;   Fc[2]= 0.782504;   Fll[2]= 0.991374; 
              cafac[2]= 0.888019;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 228324;   wouttag[2]= 18150.9;  
            fbb[2]= 0.0423463;   fcc[2]= 0.101597;   fc[2]= 0.126938;   fll[2]= 0.729118;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38064;   Fcc[3]= 1.38064;   Fc[3]= 0.767933;   Fll[3]= 0.972914; 
              cafac[3]= 0.817979;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 47329.9;   wouttag[3]= 5805.04;  
            fbb[3]= 0.0695484;   fcc[3]= 0.135457;   fc[3]= 0.125336;   fll[3]= 0.669658;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35442;   Fcc[4]= 1.35442;   Fc[4]= 0.75335;   Fll[4]= 0.954438; 
              cafac[4]= 0.837612;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10685.6;   wouttag[4]= 1805.79;  
            fbb[4]= 0.0968107;   fcc[4]= 0.159997;   fc[4]= 0.113434;   fll[4]= 0.629758;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.3242;   Fcc[5]= 1.3242;   Fc[5]= 0.736536;   Fll[5]= 0.933137; 
              cafac[5]= 0.826259;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2784.37;   wouttag[5]= 602.221;  
            fbb[5]= 0.123696;   fcc[5]= 0.192591;   fc[5]= 0.0994363;   fll[5]= 0.584276;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37339;   Fcc[6]= 1.37339;   Fc[6]= 0.7639;   Fll[6]= 0.967805; 
              cafac[6]= 0.82151;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60782.7;   wouttag[6]= 8216.36;  
            fbb[6]= 0.0768726;   fcc[6]= 0.142446;   fc[6]= 0.122032;   fll[6]= 0.658649;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.34799;   Fcc[7]= 1.34799;   Fc[7]= 0.749774;   Fll[7]= 0.949907; 
              cafac[7]= 0.834455;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13457.3;   wouttag[7]= 2409.57;  
            fbb[7]= 0.10253;   fcc[7]= 0.16693;   fc[7]= 0.110456;   fll[7]= 0.620084;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==19   || sysname=="Kbb_up" ) { 
           Fbb[1]= 1.58228;   Fcc[1]= 1.58228;   Fc[1]= 0.734715;   Fll[1]= 1.01195; 
              cafac[1]= 0.972667;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06117e+06;   wouttag[1]= 36662.9;  
            fbb[1]= 0.0190779;   fcc[1]= 0.0555065;   fc[1]= 0.102924;   fll[1]= 0.822491;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.54318;   Fcc[2]= 1.54318;   Fc[2]= 0.716562;   Fll[2]= 0.986951; 
              cafac[2]= 0.863283;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 18160.4;  
            fbb[2]= 0.0464503;   fcc[2]= 0.111444;   fc[2]= 0.116241;   fll[2]= 0.725865;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.50494;   Fcc[3]= 1.50494;   Fc[3]= 0.698803;   Fll[3]= 0.962492; 
              cafac[3]= 0.789664;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46741.7;   wouttag[3]= 5882.43;  
            fbb[3]= 0.0758097;   fcc[3]= 0.147652;   fc[3]= 0.114054;   fll[3]= 0.662485;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.46732;   Fcc[4]= 1.46732;   Fc[4]= 0.681335;   Fll[4]= 0.938431; 
              cafac[4]= 0.804983;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10564.8;   wouttag[4]= 1843.97;  
            fbb[4]= 0.10488;   fcc[4]= 0.173333;   fc[4]= 0.10259;   fll[4]= 0.619196;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.42452;   Fcc[5]= 1.42452;   Fc[5]= 0.661459;   Fll[5]= 0.911056; 
              cafac[5]= 0.789499;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2754.74;   wouttag[5]= 617.354;  
            fbb[5]= 0.133067;   fcc[5]= 0.207182;   fc[5]= 0.0893005;   fll[5]= 0.57045;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.49449;   Fcc[6]= 1.49449;   Fc[6]= 0.693951;   Fll[6]= 0.955808; 
              cafac[6]= 0.792069;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60045.2;   wouttag[6]= 8351.98;  
            fbb[6]= 0.0836508;   fcc[6]= 0.155007;   fc[6]= 0.110858;   fll[6]= 0.650485;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.45816;   Fcc[7]= 1.45816;   Fc[7]= 0.677083;   Fll[7]= 0.932576; 
              cafac[7]= 0.800963;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13307.2;   wouttag[7]= 2463.69;  
            fbb[7]= 0.110909;   fcc[7]= 0.180573;   fc[7]= 0.0997477;   fll[7]= 0.60877;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==20   || sysname=="Kbb_down" ) { 
           Fbb[1]= 1.28334;   Fcc[1]= 1.28334;   Fc[1]= 0.745902;   Fll[1]= 1.02736; 
              cafac[1]= 0.987478;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06117e+06;   wouttag[1]= 35393.4;  
            fbb[1]= 0.0154735;   fcc[1]= 0.0450197;   fc[1]= 0.104492;   fll[1]= 0.835015;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.27329;   Fcc[2]= 1.27329;   Fc[2]= 0.74006;   Fll[2]= 1.01932; 
              cafac[2]= 0.891592;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 16987.7;  
            fbb[2]= 0.0383264;   fcc[2]= 0.0919525;   fc[2]= 0.120053;   fll[2]= 0.749668;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.25887;   Fcc[3]= 1.25887;   Fc[3]= 0.731683;   Fll[3]= 1.00778; 
              cafac[3]= 0.826819;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46741.7;   wouttag[3]= 5425.21;  
            fbb[3]= 0.0634145;   fcc[3]= 0.12351;   fc[3]= 0.11942;   fll[3]= 0.693655;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.24172;   Fcc[4]= 1.24172;   Fc[4]= 0.721712;   Fll[4]= 0.994045; 
              cafac[4]= 0.852689;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10564.8;   wouttag[4]= 1692.1;  
            fbb[4]= 0.088755;   fcc[4]= 0.146683;   fc[4]= 0.10867;   fll[4]= 0.655891;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.22164;   Fcc[5]= 1.22164;   Fc[5]= 0.710042;   Fll[5]= 0.977971; 
              cafac[5]= 0.847486;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2754.74;   wouttag[5]= 565.4;  
            fbb[5]= 0.114116;   fcc[5]= 0.177676;   fc[5]= 0.0958594;   fll[5]= 0.612348;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.25415;   Fcc[6]= 1.25415;   Fc[6]= 0.728935;   Fll[6]= 1.00399; 
              cafac[6]= 0.831999;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60045.2;   wouttag[6]= 7680.79;  
            fbb[6]= 0.0701981;   fcc[6]= 0.130078;   fc[6]= 0.116446;   fll[6]= 0.683277;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.23747;   Fcc[7]= 1.23747;   Fc[7]= 0.719242;   Fll[7]= 0.990643; 
              cafac[7]= 0.850835;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13307.2;   wouttag[7]= 2258.3;  
            fbb[7]= 0.0941231;   fcc[7]= 0.153243;   fc[7]= 0.105959;   fll[7]= 0.646675;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==21   || sysname=="Kc_up" ) { 
           Fbb[1]= 1.41802;   Fcc[1]= 1.41802;   Fc[1]= 0.811289;   Fll[1]= 1.00828; 
              cafac[1]= 0.969138;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06117e+06;   wouttag[1]= 37446.7;  
            fbb[1]= 0.0170974;   fcc[1]= 0.0497443;   fc[1]= 0.113652;   fll[1]= 0.819507;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.39261;   Fcc[2]= 1.39261;   Fc[2]= 0.796749;   Fll[2]= 0.990212; 
              cafac[2]= 0.866136;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17955.1;  
            fbb[2]= 0.0419179;   fcc[2]= 0.100569;   fc[2]= 0.129249;   fll[2]= 0.728264;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.36746;   Fcc[3]= 1.36746;   Fc[3]= 0.782364;   Fll[3]= 0.972334; 
              cafac[3]= 0.79774;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46741.7;   wouttag[3]= 5728.07;  
            fbb[3]= 0.0688846;   fcc[3]= 0.134164;   fc[3]= 0.127692;   fll[3]= 0.669259;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.34243;   Fcc[4]= 1.34243;   Fc[4]= 0.768041;   Fll[4]= 0.954533; 
              cafac[4]= 0.818795;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10564.8;   wouttag[4]= 1781.53;  
            fbb[4]= 0.0959534;   fcc[4]= 0.15858;   fc[4]= 0.115646;   fll[4]= 0.62982;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.31352;   Fcc[5]= 1.31352;   Fc[5]= 0.751504;   Fll[5]= 0.93398; 
              cafac[5]= 0.809365;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2754.74;   wouttag[5]= 594.087;  
            fbb[5]= 0.1227;   fcc[5]= 0.191039;   fc[5]= 0.101457;   fll[5]= 0.584804;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.36055;   Fcc[6]= 1.36055;   Fc[6]= 0.778406;   Fll[6]= 0.967414; 
              cafac[6]= 0.801687;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60045.2;   wouttag[6]= 8105.65;  
            fbb[6]= 0.0761535;   fcc[6]= 0.141114;   fc[6]= 0.124349;   fll[6]= 0.658383;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.33628;   Fcc[7]= 1.33628;   Fc[7]= 0.764525;   Fll[7]= 0.950164; 
              cafac[7]= 0.816069;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13307.2;   wouttag[7]= 2377;  
            fbb[7]= 0.101639;   fcc[7]= 0.16548;   fc[7]= 0.11263;   fll[7]= 0.620251;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==22   || sysname=="Kc_down" ) { 
           Fbb[1]= 1.45022;   Fcc[1]= 1.45022;   Fc[1]= 0.667631;   Fll[1]= 1.03118; 
              cafac[1]= 0.991142;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06117e+06;   wouttag[1]= 34587.1;  
            fbb[1]= 0.0174856;   fcc[1]= 0.0508738;   fc[1]= 0.0935269;   fll[1]= 0.838114;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.42868;   Fcc[2]= 1.42868;   Fc[2]= 0.657715;   Fll[2]= 1.01586; 
              cafac[2]= 0.88857;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17202.2;  
            fbb[2]= 0.0430037;   fcc[2]= 0.103174;   fc[2]= 0.106695;   fll[2]= 0.747127;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.40245;   Fcc[3]= 1.40245;   Fc[3]= 0.645639;   Fll[3]= 0.997208; 
              cafac[3]= 0.818147;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46741.7;   wouttag[3]= 5588.31;  
            fbb[3]= 0.0706468;   fcc[3]= 0.137596;   fc[3]= 0.105377;   fll[3]= 0.68638;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.37346;   Fcc[4]= 1.37346;   Fc[4]= 0.632293;   Fll[4]= 0.976595; 
              cafac[4]= 0.83772;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10564.8;   wouttag[4]= 1758.65;  
            fbb[4]= 0.0981711;   fcc[4]= 0.162245;   fc[4]= 0.0952061;   fll[4]= 0.644377;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.34008;   Fcc[5]= 1.34008;   Fc[5]= 0.616929;   Fll[5]= 0.952865; 
              cafac[5]= 0.82573;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2754.74;   wouttag[5]= 590.472;  
            fbb[5]= 0.12518;   fcc[5]= 0.194902;   fc[5]= 0.0832887;   fll[5]= 0.596629;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.39442;   Fcc[6]= 1.39442;   Fc[6]= 0.641943;   Fll[6]= 0.991498; 
              cafac[6]= 0.821645;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60045.2;   wouttag[6]= 7941.61;  
            fbb[6]= 0.0780494;   fcc[6]= 0.144627;   fc[6]= 0.102549;   fll[6]= 0.674774;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.36635;   Fcc[7]= 1.36635;   Fc[7]= 0.62902;   Fll[7]= 0.971539; 
              cafac[7]= 0.834427;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13307.2;   wouttag[7]= 2350.9;  
            fbb[7]= 0.103926;   fcc[7]= 0.169203;   fc[7]= 0.092667;   fll[7]= 0.634204;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==23   || sysname=="Kll_up" ) { 
           Fbb[1]= 1.42437;   Fcc[1]= 1.42437;   Fc[1]= 0.735326;   Fll[1]= 1.02101; 
              cafac[1]= 0.981368;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06117e+06;   wouttag[1]= 35880.9;  
            fbb[1]= 0.0171739;   fcc[1]= 0.049967;   fc[1]= 0.10301;   fll[1]= 0.829849;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.40203;   Fcc[2]= 1.40203;   Fc[2]= 0.723793;   Fll[2]= 1.00499; 
              cafac[2]= 0.879065;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17513.2;  
            fbb[2]= 0.0422015;   fcc[2]= 0.10125;   fc[2]= 0.117414;   fll[2]= 0.739135;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.37717;   Fcc[3]= 1.37717;   Fc[3]= 0.710959;   Fll[3]= 0.987173; 
              cafac[3]= 0.809914;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46741.7;   wouttag[3]= 5637.39;  
            fbb[3]= 0.0693734;   fcc[3]= 0.135116;   fc[3]= 0.116038;   fll[3]= 0.679473;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35079;   Fcc[4]= 1.35079;   Fc[4]= 0.697341;   Fll[4]= 0.968264; 
              cafac[4]= 0.830574;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10564.8;   wouttag[4]= 1763.89;  
            fbb[4]= 0.096551;   fcc[4]= 0.159568;   fc[4]= 0.105;   fll[4]= 0.638881;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32035;   Fcc[5]= 1.32035;   Fc[5]= 0.681626;   Fll[5]= 0.946444; 
              cafac[5]= 0.820165;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2754.74;   wouttag[5]= 590.29;  
            fbb[5]= 0.123337;   fcc[5]= 0.192032;   fc[5]= 0.0920231;   fll[5]= 0.592608;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.36987;   Fcc[6]= 1.36987;   Fc[6]= 0.707192;   Fll[6]= 0.981942; 
              cafac[6]= 0.813726;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60045.2;   wouttag[6]= 7994.33;  
            fbb[6]= 0.0766755;   fcc[6]= 0.142081;   fc[6]= 0.112973;   fll[6]= 0.66827;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.34431;   Fcc[7]= 1.34431;   Fc[7]= 0.693998;   Fll[7]= 0.963622; 
              cafac[7]= 0.827627;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13307.2;   wouttag[7]= 2355.72;  
            fbb[7]= 0.10225;   fcc[7]= 0.166474;   fc[7]= 0.10224;   fll[7]= 0.629036;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==24   || sysname=="Kll_down" ) { 
           Fbb[1]= 1.44364;   Fcc[1]= 1.44364;   Fc[1]= 0.745274;   Fll[1]= 1.01817; 
              cafac[1]= 0.978646;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06117e+06;   wouttag[1]= 36187.1;  
            fbb[1]= 0.0174063;   fcc[1]= 0.050643;   fc[1]= 0.104404;   fll[1]= 0.827547;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.4189;   Fcc[2]= 1.4189;   Fc[2]= 0.732501;   Fll[2]= 1.00073; 
              cafac[2]= 0.875332;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17654.6;  
            fbb[2]= 0.0427093;   fcc[2]= 0.102468;   fc[2]= 0.118827;   fll[2]= 0.735996;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.39239;   Fcc[3]= 1.39239;   Fc[3]= 0.718815;   Fll[3]= 0.982027; 
              cafac[3]= 0.805692;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46741.7;   wouttag[3]= 5681;  
            fbb[3]= 0.0701399;   fcc[3]= 0.136609;   fc[3]= 0.11732;   fll[3]= 0.675931;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.36481;   Fcc[4]= 1.36481;   Fc[4]= 0.704581;   Fll[4]= 0.962582; 
              cafac[4]= 0.8257;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10564.8;   wouttag[4]= 1776.61;  
            fbb[4]= 0.0975534;   fcc[4]= 0.161224;   fc[4]= 0.106091;   fll[4]= 0.635131;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.33305;   Fcc[5]= 1.33305;   Fc[5]= 0.688185;   Fll[5]= 0.940182; 
              cafac[5]= 0.814739;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2754.74;   wouttag[5]= 594.324;  
            fbb[5]= 0.124524;   fcc[5]= 0.19388;   fc[5]= 0.0929087;   fll[5]= 0.588688;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.38476;   Fcc[6]= 1.38476;   Fc[6]= 0.714876;   Fll[6]= 0.976646; 
              cafac[6]= 0.809337;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60045.2;   wouttag[6]= 8055.27;  
            fbb[6]= 0.0775086;   fcc[6]= 0.143625;   fc[6]= 0.1142;   fll[6]= 0.664666;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.35805;   Fcc[7]= 1.35805;   Fc[7]= 0.701091;   Fll[7]= 0.957814; 
              cafac[7]= 0.822639;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13307.2;   wouttag[7]= 2372.56;  
            fbb[7]= 0.103295;   fcc[7]= 0.168176;   fc[7]= 0.103285;   fll[7]= 0.625245;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==33   || sysname=="WbbWccjet1_up" ) { 
           Fbb[1]= 1.42628;   Fcc[1]= 1.42628;   Fc[1]= 0.736315;   Fll[1]= 1.01416; 
              cafac[1]= 0.97563;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 41132 ;   woutpretag[1]= 1.05642e+06;   wouttag[1]= 37346.8;  
            fbb[1]= 0.0214963;   fcc[1]= 0.0625428;   fc[1]= 0.101873;   fll[1]= 0.814088;  
            fmcbb[1]= 0.0150715;   fmccc[1]= 0.0438501;   fmcc[1]= 0.138355;   fmcll[1]= 0.802723;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.877209;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17583.5;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38474;   Fcc[3]= 1.38474;   Fc[3]= 0.714865;   Fll[3]= 0.984614; 
              cafac[3]= 0.807815;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46741.7;   wouttag[3]= 5659.07;  
            fbb[3]= 0.0697545;   fcc[3]= 0.135859;   fc[3]= 0.116675;   fll[3]= 0.677712;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35777;   Fcc[4]= 1.35777;   Fc[4]= 0.700942;   Fll[4]= 0.965438; 
              cafac[4]= 0.828149;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10564.8;   wouttag[4]= 1770.22;  
            fbb[4]= 0.0970496;   fcc[4]= 0.160392;   fc[4]= 0.105543;   fll[4]= 0.637016;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32667;   Fcc[5]= 1.32667;   Fc[5]= 0.68489;   Fll[5]= 0.943328; 
              cafac[5]= 0.817465;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2754.74;   wouttag[5]= 592.298;  
            fbb[5]= 0.123928;   fcc[5]= 0.192951;   fc[5]= 0.0924638;   fll[5]= 0.590657;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37727;   Fcc[6]= 1.37727;   Fc[6]= 0.711013;   Fll[6]= 0.979308; 
              cafac[6]= 0.811543;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60045.2;   wouttag[6]= 8024.64;  
            fbb[6]= 0.0770898;   fcc[6]= 0.142849;   fc[6]= 0.113583;   fll[6]= 0.666478;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.35115;   Fcc[7]= 1.35115;   Fc[7]= 0.697526;   Fll[7]= 0.960733; 
              cafac[7]= 0.825146;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13307.2;   wouttag[7]= 2364.1;  
            fbb[7]= 0.10277;   fcc[7]= 0.167321;   fc[7]= 0.102759;   fll[7]= 0.62715;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==34   || sysname=="WbbWccjet1_down" ) { 
           Fbb[1]= 1.44168;   Fcc[1]= 1.44168;   Fc[1]= 0.744261;   Fll[1]= 1.0251; 
              cafac[1]= 0.984491;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 39094.6 ;   woutpretag[1]= 1.06602e+06;   wouttag[1]= 34692.9;  
            fbb[1]= 0.013037;   fcc[1]= 0.0379306;   fc[1]= 0.105551;   fll[1]= 0.843481;  
            fmcbb[1]= 0.00904292;   fmccc[1]= 0.0263101;   fmcc[1]= 0.14182;   fmcll[1]= 0.822827;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.877209;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17583.5;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38474;   Fcc[3]= 1.38474;   Fc[3]= 0.714865;   Fll[3]= 0.984614; 
              cafac[3]= 0.807815;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46741.7;   wouttag[3]= 5659.07;  
            fbb[3]= 0.0697545;   fcc[3]= 0.135859;   fc[3]= 0.116675;   fll[3]= 0.677712;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35777;   Fcc[4]= 1.35777;   Fc[4]= 0.700942;   Fll[4]= 0.965438; 
              cafac[4]= 0.828149;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10564.8;   wouttag[4]= 1770.22;  
            fbb[4]= 0.0970496;   fcc[4]= 0.160392;   fc[4]= 0.105543;   fll[4]= 0.637016;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32667;   Fcc[5]= 1.32667;   Fc[5]= 0.68489;   Fll[5]= 0.943328; 
              cafac[5]= 0.817465;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2754.74;   wouttag[5]= 592.298;  
            fbb[5]= 0.123928;   fcc[5]= 0.192951;   fc[5]= 0.0924638;   fll[5]= 0.590657;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37727;   Fcc[6]= 1.37727;   Fc[6]= 0.711013;   Fll[6]= 0.979308; 
              cafac[6]= 0.811543;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60045.2;   wouttag[6]= 8024.64;  
            fbb[6]= 0.0770898;   fcc[6]= 0.142849;   fc[6]= 0.113583;   fll[6]= 0.666478;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.35115;   Fcc[7]= 1.35115;   Fc[7]= 0.697526;   Fll[7]= 0.960733; 
              cafac[7]= 0.825146;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13307.2;   wouttag[7]= 2364.1;  
            fbb[7]= 0.10277;   fcc[7]= 0.167321;   fc[7]= 0.102759;   fll[7]= 0.62715;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==35   || sysname=="WbbWccjet3_up" ) { 
           Fbb[1]= 1.43394;   Fcc[1]= 1.43394;   Fc[1]= 0.740266;   Fll[1]= 1.0196; 
              cafac[1]= 0.980017;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06117e+06;   wouttag[1]= 36033;  
            fbb[1]= 0.0172893;   fcc[1]= 0.0503027;   fc[1]= 0.103702;   fll[1]= 0.828706;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.877209;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17583.5;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.37459;   Fcc[3]= 1.37459;   Fc[3]= 0.709628;   Fll[3]= 0.977401; 
              cafac[3]= 0.805638;  
              winpretag[3]= 57862 ;   wintag[3]= 6813.51 ;   woutpretag[3]= 46615.8;   wouttag[3]= 5905.23;  
            fbb[3]= 0.0768603;   fcc[3]= 0.149698;   fc[3]= 0.113599;   fll[3]= 0.659843;  
            fmcbb[3]= 0.0559151;   fmccc[3]= 0.108904;   fmcc[3]= 0.160082;   fmcll[3]= 0.675099;  
           Fbb[4]= 1.35777;   Fcc[4]= 1.35777;   Fc[4]= 0.700942;   Fll[4]= 0.965438; 
              cafac[4]= 0.828149;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10564.8;   wouttag[4]= 1770.22;  
            fbb[4]= 0.0970496;   fcc[4]= 0.160392;   fc[4]= 0.105543;   fll[4]= 0.637016;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32667;   Fcc[5]= 1.32667;   Fc[5]= 0.68489;   Fll[5]= 0.943328; 
              cafac[5]= 0.817465;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2754.74;   wouttag[5]= 592.298;  
            fbb[5]= 0.123928;   fcc[5]= 0.192951;   fc[5]= 0.0924638;   fll[5]= 0.590657;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.36941;   Fcc[6]= 1.36941;   Fc[6]= 0.706955;   Fll[6]= 0.973719; 
              cafac[6]= 0.809987;  
              winpretag[6]= 73989 ;   wintag[6]= 9401.41 ;   woutpretag[6]= 59930.1;   wouttag[6]= 8266.77;  
            fbb[6]= 0.082584;   fcc[6]= 0.153591;   fc[6]= 0.111204;   fll[6]= 0.65262;  
            fmcbb[6]= 0.0603061;   fmccc[6]= 0.112159;   fmcc[6]= 0.1573;   fmcll[6]= 0.670235;  
           Fbb[7]= 1.35115;   Fcc[7]= 1.35115;   Fc[7]= 0.697526;   Fll[7]= 0.960733; 
              cafac[7]= 0.825146;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13307.2;   wouttag[7]= 2364.1;  
            fbb[7]= 0.10277;   fcc[7]= 0.167321;   fc[7]= 0.102759;   fll[7]= 0.62715;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==36   || sysname=="WbbWccjet3_down" ) { 
           Fbb[1]= 1.43394;   Fcc[1]= 1.43394;   Fc[1]= 0.740266;   Fll[1]= 1.0196; 
              cafac[1]= 0.980017;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06117e+06;   wouttag[1]= 36033;  
            fbb[1]= 0.0172893;   fcc[1]= 0.0503027;   fc[1]= 0.103702;   fll[1]= 0.828706;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.877209;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17583.5;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.39503;   Fcc[3]= 1.39503;   Fc[3]= 0.72018;   Fll[3]= 0.991935; 
              cafac[3]= 0.810035;  
              winpretag[3]= 57862 ;   wintag[3]= 6324.69 ;   woutpretag[3]= 46870.2;   wouttag[3]= 5407.9;  
            fbb[3]= 0.0625431;   fcc[3]= 0.121813;   fc[3]= 0.119797;   fll[3]= 0.695847;  
            fmcbb[3]= 0.0448328;   fmccc[3]= 0.0873194;   fmcc[3]= 0.166343;   fmcll[3]= 0.701504;  
           Fbb[4]= 1.35777;   Fcc[4]= 1.35777;   Fc[4]= 0.700942;   Fll[4]= 0.965438; 
              cafac[4]= 0.828149;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10564.8;   wouttag[4]= 1770.22;  
            fbb[4]= 0.0970496;   fcc[4]= 0.160392;   fc[4]= 0.105543;   fll[4]= 0.637016;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32667;   Fcc[5]= 1.32667;   Fc[5]= 0.68489;   Fll[5]= 0.943328; 
              cafac[5]= 0.817465;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2754.74;   wouttag[5]= 592.298;  
            fbb[5]= 0.123928;   fcc[5]= 0.192951;   fc[5]= 0.0924638;   fll[5]= 0.590657;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.38522;   Fcc[6]= 1.38522;   Fc[6]= 0.715118;   Fll[6]= 0.984962; 
              cafac[6]= 0.813123;  
              winpretag[6]= 73989 ;   wintag[6]= 8912.59 ;   woutpretag[6]= 60162.2;   wouttag[6]= 7778.76;  
            fbb[6]= 0.0715322;   fcc[6]= 0.131982;   fc[6]= 0.11599;   fll[6]= 0.680496;  
            fmcbb[6]= 0.0516394;   fmccc[6]= 0.0952788;   fmcc[6]= 0.162197;   fmcll[6]= 0.690885;  
           Fbb[7]= 1.35115;   Fcc[7]= 1.35115;   Fc[7]= 0.697526;   Fll[7]= 0.960733; 
              cafac[7]= 0.825146;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13307.2;   wouttag[7]= 2364.1;  
            fbb[7]= 0.10277;   fcc[7]= 0.167321;   fc[7]= 0.102759;   fll[7]= 0.62715;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==37   || sysname=="WbbWccjet4_up" ) { 
           Fbb[1]= 1.43394;   Fcc[1]= 1.43394;   Fc[1]= 0.740266;   Fll[1]= 1.0196; 
              cafac[1]= 0.980017;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06117e+06;   wouttag[1]= 36033;  
            fbb[1]= 0.0172893;   fcc[1]= 0.0503027;   fc[1]= 0.103702;   fll[1]= 0.828706;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.877209;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17583.5;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38474;   Fcc[3]= 1.38474;   Fc[3]= 0.714865;   Fll[3]= 0.984614; 
              cafac[3]= 0.807815;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46741.7;   wouttag[3]= 5659.07;  
            fbb[3]= 0.0697545;   fcc[3]= 0.135859;   fc[3]= 0.116675;   fll[3]= 0.677712;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.34093;   Fcc[4]= 1.34093;   Fc[4]= 0.692251;   Fll[4]= 0.953466; 
              cafac[4]= 0.823794;  
              winpretag[4]= 12757.2 ;   wintag[4]= 2060.72 ;   woutpretag[4]= 10509.3;   wouttag[4]= 1884.33;  
            fbb[4]= 0.110223;   fcc[4]= 0.182163;   fc[4]= 0.100576;   fll[4]= 0.607038;  
            fmcbb[4]= 0.082199;   fmccc[4]= 0.135849;   fmcc[4]= 0.145288;   fmcll[4]= 0.636664;  
           Fbb[5]= 1.32667;   Fcc[5]= 1.32667;   Fc[5]= 0.68489;   Fll[5]= 0.943328; 
              cafac[5]= 0.817465;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2754.74;   wouttag[5]= 592.298;  
            fbb[5]= 0.123928;   fcc[5]= 0.192951;   fc[5]= 0.0924638;   fll[5]= 0.590657;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37426;   Fcc[6]= 1.37426;   Fc[6]= 0.709455;   Fll[6]= 0.977162; 
              cafac[6]= 0.810669;  
              winpretag[6]= 73989 ;   wintag[6]= 9275.58 ;   woutpretag[6]= 59980.6;   wouttag[6]= 8144.91;  
            fbb[6]= 0.0794614;   fcc[6]= 0.146735;   fc[6]= 0.112688;   fll[6]= 0.661116;  
            fmcbb[6]= 0.0578214;   fmccc[6]= 0.106774;   fmcc[6]= 0.158837;   fmcll[6]= 0.676567;  
           Fbb[7]= 1.33792;   Fcc[7]= 1.33792;   Fc[7]= 0.6907;   Fll[7]= 0.95133; 
              cafac[7]= 0.822096;  
              winpretag[7]= 16127 ;   wintag[7]= 2706.48 ;   woutpretag[7]= 13258;   wouttag[7]= 2476.93;  
            fbb[7]= 0.113111;   fcc[7]= 0.184437;   fc[7]= 0.0988664;   fll[7]= 0.603586;  
            fmcbb[7]= 0.0845421;   fmccc[7]= 0.137853;   fmcc[7]= 0.14314;   fmcll[7]= 0.634465;  
      }  
      else if ( idsys==38   || sysname=="WbbWccjet4_down" ) { 
           Fbb[1]= 1.43394;   Fcc[1]= 1.43394;   Fc[1]= 0.740266;   Fll[1]= 1.0196; 
              cafac[1]= 0.980017;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06117e+06;   wouttag[1]= 36033;  
            fbb[1]= 0.0172893;   fcc[1]= 0.0503027;   fc[1]= 0.103702;   fll[1]= 0.828706;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.877209;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17583.5;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38474;   Fcc[3]= 1.38474;   Fc[3]= 0.714865;   Fll[3]= 0.984614; 
              cafac[3]= 0.807815;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46741.7;   wouttag[3]= 5659.07;  
            fbb[3]= 0.0697545;   fcc[3]= 0.135859;   fc[3]= 0.116675;   fll[3]= 0.677712;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.37503;   Fcc[4]= 1.37503;   Fc[4]= 0.709855;   Fll[4]= 0.977714; 
              cafac[4]= 0.832663;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1823.56 ;   woutpretag[4]= 10622.4;   wouttag[4]= 1651.95;  
            fbb[4]= 0.0835411;   fcc[4]= 0.138067;   fc[4]= 0.110636;   fll[4]= 0.667756;  
            fmcbb[4]= 0.0607558;   fmccc[4]= 0.10041;   fmcc[4]= 0.155857;   fmcll[4]= 0.682977;  
           Fbb[5]= 1.32667;   Fcc[5]= 1.32667;   Fc[5]= 0.68489;   Fll[5]= 0.943328; 
              cafac[5]= 0.817465;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2754.74;   wouttag[5]= 592.298;  
            fbb[5]= 0.123928;   fcc[5]= 0.192951;   fc[5]= 0.0924638;   fll[5]= 0.590657;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.3803;   Fcc[6]= 1.3803;   Fc[6]= 0.712578;   Fll[6]= 0.981464; 
              cafac[6]= 0.812423;  
              winpretag[6]= 73989 ;   wintag[6]= 9038.42 ;   woutpretag[6]= 60110.4;   wouttag[6]= 7903.57;  
            fbb[6]= 0.0747078;   fcc[6]= 0.138946;   fc[6]= 0.114483;   fll[6]= 0.671863;  
            fmcbb[6]= 0.0541242;   fmccc[6]= 0.100664;   fmcc[6]= 0.16066;   fmcll[6]= 0.684553;  
           Fbb[7]= 1.36464;   Fcc[7]= 1.36464;   Fc[7]= 0.704489;   Fll[7]= 0.970323; 
              cafac[7]= 0.82828;  
              winpretag[7]= 16127 ;   wintag[7]= 2469.32 ;   woutpretag[7]= 13357.7;   wouttag[7]= 2248.14;  
            fbb[7]= 0.0922216;   fcc[7]= 0.149863;   fc[7]= 0.10673;   fll[7]= 0.651185;  
            fmcbb[7]= 0.0675796;   fmccc[7]= 0.109819;   fmcc[7]= 0.1515;   fmcll[7]= 0.671101;  
      }  
      else if ( idsys==73   || sysname=="WbbWccjet5_up" ) { 
           Fbb[1]= 1.43394;   Fcc[1]= 1.43394;   Fc[1]= 0.740266;   Fll[1]= 1.0196; 
              cafac[1]= 0.980017;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06117e+06;   wouttag[1]= 36033;  
            fbb[1]= 0.0172893;   fcc[1]= 0.0503027;   fc[1]= 0.103702;   fll[1]= 0.828706;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.877209;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17583.5;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38474;   Fcc[3]= 1.38474;   Fc[3]= 0.714865;   Fll[3]= 0.984614; 
              cafac[3]= 0.807815;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46741.7;   wouttag[3]= 5659.07;  
            fbb[3]= 0.0697545;   fcc[3]= 0.135859;   fc[3]= 0.116675;   fll[3]= 0.677712;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35777;   Fcc[4]= 1.35777;   Fc[4]= 0.700942;   Fll[4]= 0.965438; 
              cafac[4]= 0.828149;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10564.8;   wouttag[4]= 1770.22;  
            fbb[4]= 0.0970496;   fcc[4]= 0.160392;   fc[4]= 0.105543;   fll[4]= 0.637016;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.28837;   Fcc[5]= 1.28837;   Fc[5]= 0.665117;   Fll[5]= 0.916094; 
              cafac[5]= 0.795047;  
              winpretag[5]= 3369.85 ;   wintag[5]= 734.467 ;   woutpretag[5]= 2679.19;   wouttag[5]= 659.573;  
            fbb[5]= 0.155251;   fcc[5]= 0.241721;   fc[5]= 0.0816227;   fll[5]= 0.521405;  
            fmcbb[5]= 0.120502;   fmccc[5]= 0.187618;   fmcc[5]= 0.122719;   fmcll[5]= 0.569161;  
           Fbb[6]= 1.37534;   Fcc[6]= 1.37534;   Fc[6]= 0.710015;   Fll[6]= 0.977934; 
              cafac[6]= 0.80995;  
              winpretag[6]= 73989 ;   wintag[6]= 9245.7 ;   woutpretag[6]= 59927.4;   wouttag[6]= 8106.56;  
            fbb[6]= 0.0786785;   fcc[6]= 0.14529;   fc[6]= 0.113027;   fll[6]= 0.663005;  
            fmcbb[6]= 0.0572066;   fmccc[6]= 0.10564;   fmcc[6]= 0.159189;   fmcll[6]= 0.677965;  
           Fbb[7]= 1.34265;   Fcc[7]= 1.34265;   Fc[7]= 0.693141;   Fll[7]= 0.954693; 
              cafac[7]= 0.81887;  
              winpretag[7]= 16127 ;   wintag[7]= 2676.61 ;   woutpretag[7]= 13205.9;   wouttag[7]= 2438.86;  
            fbb[7]= 0.109724;   fcc[7]= 0.178102;   fc[7]= 0.100334;   fll[7]= 0.61184;  
            fmcbb[7]= 0.0817215;   fmccc[7]= 0.132649;   fmcc[7]= 0.144753;   fmcll[7]= 0.640877;  
      }  
      else if ( idsys==74   || sysname=="WbbWccjet5_down" ) { 
           Fbb[1]= 1.43394;   Fcc[1]= 1.43394;   Fc[1]= 0.740266;   Fll[1]= 1.0196; 
              cafac[1]= 0.980017;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06117e+06;   wouttag[1]= 36033;  
            fbb[1]= 0.0172893;   fcc[1]= 0.0503027;   fc[1]= 0.103702;   fll[1]= 0.828706;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.877209;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17583.5;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38474;   Fcc[3]= 1.38474;   Fc[3]= 0.714865;   Fll[3]= 0.984614; 
              cafac[3]= 0.807815;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46741.7;   wouttag[3]= 5659.07;  
            fbb[3]= 0.0697545;   fcc[3]= 0.135859;   fc[3]= 0.116675;   fll[3]= 0.677712;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35777;   Fcc[4]= 1.35777;   Fc[4]= 0.700942;   Fll[4]= 0.965438; 
              cafac[4]= 0.828149;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10564.8;   wouttag[4]= 1770.22;  
            fbb[4]= 0.0970496;   fcc[4]= 0.160392;   fc[4]= 0.105543;   fll[4]= 0.637016;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.36732;   Fcc[5]= 1.36732;   Fc[5]= 0.705874;   Fll[5]= 0.972231; 
              cafac[5]= 0.842682;  
              winpretag[5]= 3369.85 ;   wintag[5]= 557.055 ;   woutpretag[5]= 2839.72;   wouttag[5]= 516.622;  
            fbb[5]= 0.0906845;   fcc[5]= 0.141193;   fc[5]= 0.103969;   fll[5]= 0.664153;  
            fmcbb[5]= 0.0663228;   fmccc[5]= 0.103263;   fmcc[5]= 0.147291;   fmcll[5]= 0.683123;  
           Fbb[6]= 1.37921;   Fcc[6]= 1.37921;   Fc[6]= 0.712014;   Fll[6]= 0.980687; 
              cafac[6]= 0.813147;  
              winpretag[6]= 73989 ;   wintag[6]= 9068.29 ;   woutpretag[6]= 60163.9;   wouttag[6]= 7942.16;  
            fbb[6]= 0.0754966;   fcc[6]= 0.140401;   fc[6]= 0.114142;   fll[6]= 0.669961;  
            fmcbb[6]= 0.054739;   fmccc[6]= 0.101798;   fmcc[6]= 0.160308;   fmcll[6]= 0.683155;  
           Fbb[7]= 1.35975;   Fcc[7]= 1.35975;   Fc[7]= 0.701967;   Fll[7]= 0.966849; 
              cafac[7]= 0.8316;  
              winpretag[7]= 16127 ;   wintag[7]= 2499.19 ;   woutpretag[7]= 13411.2;   wouttag[7]= 2287.2;  
            fbb[7]= 0.0957269;   fcc[7]= 0.156402;   fc[7]= 0.105216;   fll[7]= 0.642655;  
            fmcbb[7]= 0.0704003;   fmccc[7]= 0.115023;   fmcc[7]= 0.149887;   fmcll[7]= 0.66469;  
      }  
      else if ( idsys==41   || sysname=="Wcjet1_up" ) { 
           Fbb[1]= 1.44927;   Fcc[1]= 1.44927;   Fc[1]= 0.748181;   Fll[1]= 1.0305; 
              cafac[1]= 1.01677;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 45269.8 ;   woutpretag[1]= 1.10097e+06;   wouttag[1]= 41410.5;  
            fbb[1]= 0.0167625;   fcc[1]= 0.0487699;   fc[1]= 0.131014;   fll[1]= 0.803454;  
            fmcbb[1]= 0.0115662;   fmccc[1]= 0.0336514;   fmcc[1]= 0.175109;   fmcll[1]= 0.779673;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.877209;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17583.5;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38474;   Fcc[3]= 1.38474;   Fc[3]= 0.714865;   Fll[3]= 0.984614; 
              cafac[3]= 0.807815;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46741.7;   wouttag[3]= 5659.07;  
            fbb[3]= 0.0697545;   fcc[3]= 0.135859;   fc[3]= 0.116675;   fll[3]= 0.677712;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35777;   Fcc[4]= 1.35777;   Fc[4]= 0.700942;   Fll[4]= 0.965438; 
              cafac[4]= 0.828149;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10564.8;   wouttag[4]= 1770.22;  
            fbb[4]= 0.0970496;   fcc[4]= 0.160392;   fc[4]= 0.105543;   fll[4]= 0.637016;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32667;   Fcc[5]= 1.32667;   Fc[5]= 0.68489;   Fll[5]= 0.943328; 
              cafac[5]= 0.817465;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2754.74;   wouttag[5]= 592.298;  
            fbb[5]= 0.123928;   fcc[5]= 0.192951;   fc[5]= 0.0924638;   fll[5]= 0.590657;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37727;   Fcc[6]= 1.37727;   Fc[6]= 0.711013;   Fll[6]= 0.979308; 
              cafac[6]= 0.811543;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60045.2;   wouttag[6]= 8024.64;  
            fbb[6]= 0.0770898;   fcc[6]= 0.142849;   fc[6]= 0.113583;   fll[6]= 0.666478;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.35115;   Fcc[7]= 1.35115;   Fc[7]= 0.697526;   Fll[7]= 0.960733; 
              cafac[7]= 0.825146;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13307.2;   wouttag[7]= 2364.1;  
            fbb[7]= 0.10277;   fcc[7]= 0.167321;   fc[7]= 0.102759;   fll[7]= 0.62715;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==42   || sysname=="Wcjet1_down" ) { 
           Fbb[1]= 1.41893;   Fcc[1]= 1.41893;   Fc[1]= 0.732518;   Fll[1]= 1.00893; 
              cafac[1]= 0.946522;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 34956.9 ;   woutpretag[1]= 1.0249e+06;   wouttag[1]= 31131.8;  
            fbb[1]= 0.0178051;   fcc[1]= 0.0518034;   fc[1]= 0.0769625;   fll[1]= 0.853429;  
            fmcbb[1]= 0.0125483;   fmccc[1]= 0.0365088;   fmcc[1]= 0.105066;   fmcll[1]= 0.845877;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.877209;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17583.5;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38474;   Fcc[3]= 1.38474;   Fc[3]= 0.714865;   Fll[3]= 0.984614; 
              cafac[3]= 0.807815;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46741.7;   wouttag[3]= 5659.07;  
            fbb[3]= 0.0697545;   fcc[3]= 0.135859;   fc[3]= 0.116675;   fll[3]= 0.677712;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35777;   Fcc[4]= 1.35777;   Fc[4]= 0.700942;   Fll[4]= 0.965438; 
              cafac[4]= 0.828149;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10564.8;   wouttag[4]= 1770.22;  
            fbb[4]= 0.0970496;   fcc[4]= 0.160392;   fc[4]= 0.105543;   fll[4]= 0.637016;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32667;   Fcc[5]= 1.32667;   Fc[5]= 0.68489;   Fll[5]= 0.943328; 
              cafac[5]= 0.817465;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2754.74;   wouttag[5]= 592.298;  
            fbb[5]= 0.123928;   fcc[5]= 0.192951;   fc[5]= 0.0924638;   fll[5]= 0.590657;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37727;   Fcc[6]= 1.37727;   Fc[6]= 0.711013;   Fll[6]= 0.979308; 
              cafac[6]= 0.811543;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60045.2;   wouttag[6]= 8024.64;  
            fbb[6]= 0.0770898;   fcc[6]= 0.142849;   fc[6]= 0.113583;   fll[6]= 0.666478;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.35115;   Fcc[7]= 1.35115;   Fc[7]= 0.697526;   Fll[7]= 0.960733; 
              cafac[7]= 0.825146;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13307.2;   wouttag[7]= 2364.1;  
            fbb[7]= 0.10277;   fcc[7]= 0.167321;   fc[7]= 0.102759;   fll[7]= 0.62715;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==43   || sysname=="Wcjet3_up" ) { 
           Fbb[1]= 1.43394;   Fcc[1]= 1.43394;   Fc[1]= 0.740266;   Fll[1]= 1.0196; 
              cafac[1]= 0.980017;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06117e+06;   wouttag[1]= 36033;  
            fbb[1]= 0.0172893;   fcc[1]= 0.0503027;   fc[1]= 0.103702;   fll[1]= 0.828706;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.877209;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17583.5;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.39482;   Fcc[3]= 1.39482;   Fc[3]= 0.720071;   Fll[3]= 0.991785; 
              cafac[3]= 0.826765;  
              winpretag[3]= 57862 ;   wintag[3]= 6753.81 ;   woutpretag[3]= 47838.3;   wouttag[3]= 5895.21;  
            fbb[3]= 0.068481;   fcc[3]= 0.133378;   fc[3]= 0.132803;   fll[3]= 0.665338;  
            fmcbb[3]= 0.0490966;   fmccc[3]= 0.0956239;   fmcc[3]= 0.18443;   fmcll[3]= 0.670849;  
           Fbb[4]= 1.35777;   Fcc[4]= 1.35777;   Fc[4]= 0.700942;   Fll[4]= 0.965438; 
              cafac[4]= 0.828149;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10564.8;   wouttag[4]= 1770.22;  
            fbb[4]= 0.0970496;   fcc[4]= 0.160392;   fc[4]= 0.105543;   fll[4]= 0.637016;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32667;   Fcc[5]= 1.32667;   Fc[5]= 0.68489;   Fll[5]= 0.943328; 
              cafac[5]= 0.817465;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2754.74;   wouttag[5]= 592.298;  
            fbb[5]= 0.123928;   fcc[5]= 0.192951;   fc[5]= 0.0924638;   fll[5]= 0.590657;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.38506;   Fcc[6]= 1.38506;   Fc[6]= 0.715034;   Fll[6]= 0.984847; 
              cafac[6]= 0.825788;  
              winpretag[6]= 73989 ;   wintag[6]= 9341.71 ;   woutpretag[6]= 61099.2;   wouttag[6]= 8272.37;  
            fbb[6]= 0.0761423;   fcc[6]= 0.140962;   fc[6]= 0.12609;   fll[6]= 0.656805;  
            fmcbb[6]= 0.0549739;   fmccc[6]= 0.101773;   fmcc[6]= 0.176341;   fmcll[6]= 0.666911;  
           Fbb[7]= 1.35115;   Fcc[7]= 1.35115;   Fc[7]= 0.697526;   Fll[7]= 0.960733; 
              cafac[7]= 0.825146;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13307.2;   wouttag[7]= 2364.1;  
            fbb[7]= 0.10277;   fcc[7]= 0.167321;   fc[7]= 0.102759;   fll[7]= 0.62715;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==44   || sysname=="Wcjet3_down" ) { 
           Fbb[1]= 1.43394;   Fcc[1]= 1.43394;   Fc[1]= 0.740266;   Fll[1]= 1.0196; 
              cafac[1]= 0.980017;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06117e+06;   wouttag[1]= 36033;  
            fbb[1]= 0.0172893;   fcc[1]= 0.0503027;   fc[1]= 0.103702;   fll[1]= 0.828706;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.877209;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17583.5;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.3748;   Fcc[3]= 1.3748;   Fc[3]= 0.709734;   Fll[3]= 0.977547; 
              cafac[3]= 0.789967;  
              winpretag[3]= 57862 ;   wintag[3]= 6384.39 ;   woutpretag[3]= 45709.1;   wouttag[3]= 5436.69;  
            fbb[3]= 0.0710099;   fcc[3]= 0.138304;   fc[3]= 0.100779;   fll[3]= 0.689908;  
            fmcbb[3]= 0.0516512;   fmccc[3]= 0.100599;   fmcc[3]= 0.141995;   fmcll[3]= 0.705754;  
           Fbb[4]= 1.35777;   Fcc[4]= 1.35777;   Fc[4]= 0.700942;   Fll[4]= 0.965438; 
              cafac[4]= 0.828149;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10564.8;   wouttag[4]= 1770.22;  
            fbb[4]= 0.0970496;   fcc[4]= 0.160392;   fc[4]= 0.105543;   fll[4]= 0.637016;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32667;   Fcc[5]= 1.32667;   Fc[5]= 0.68489;   Fll[5]= 0.943328; 
              cafac[5]= 0.817465;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2754.74;   wouttag[5]= 592.298;  
            fbb[5]= 0.123928;   fcc[5]= 0.192951;   fc[5]= 0.0924638;   fll[5]= 0.590657;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.36957;   Fcc[6]= 1.36957;   Fc[6]= 0.707037;   Fll[6]= 0.973832; 
              cafac[6]= 0.797933;  
              winpretag[6]= 73989 ;   wintag[6]= 8972.29 ;   woutpretag[6]= 59038.2;   wouttag[6]= 7787.94;  
            fbb[6]= 0.0780267;   fcc[6]= 0.144715;   fc[6]= 0.101216;   fll[6]= 0.676042;  
            fmcbb[6]= 0.0569717;   fmccc[6]= 0.105664;   fmcc[6]= 0.143156;   fmcll[6]= 0.694209;  
           Fbb[7]= 1.35115;   Fcc[7]= 1.35115;   Fc[7]= 0.697526;   Fll[7]= 0.960733; 
              cafac[7]= 0.825146;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13307.2;   wouttag[7]= 2364.1;  
            fbb[7]= 0.10277;   fcc[7]= 0.167321;   fc[7]= 0.102759;   fll[7]= 0.62715;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==45   || sysname=="Wcjet4_up" ) { 
           Fbb[1]= 1.43394;   Fcc[1]= 1.43394;   Fc[1]= 0.740266;   Fll[1]= 1.0196; 
              cafac[1]= 0.980017;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06117e+06;   wouttag[1]= 36033;  
            fbb[1]= 0.0172893;   fcc[1]= 0.0503027;   fc[1]= 0.103702;   fll[1]= 0.828706;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.877209;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17583.5;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38474;   Fcc[3]= 1.38474;   Fc[3]= 0.714865;   Fll[3]= 0.984614; 
              cafac[3]= 0.807815;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46741.7;   wouttag[3]= 5659.07;  
            fbb[3]= 0.0697545;   fcc[3]= 0.135859;   fc[3]= 0.116675;   fll[3]= 0.677712;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.36938;   Fcc[4]= 1.36938;   Fc[4]= 0.706939;   Fll[4]= 0.973697; 
              cafac[4]= 0.849801;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1981.97 ;   woutpretag[4]= 10841.1;   wouttag[4]= 1837.1;  
            fbb[4]= 0.0951037;   fcc[4]= 0.157176;   fc[4]= 0.123477;   fll[4]= 0.624243;  
            fmcbb[4]= 0.0694502;   fmccc[4]= 0.114779;   fmcc[4]= 0.174664;   fmcll[4]= 0.641107;  
           Fbb[5]= 1.32667;   Fcc[5]= 1.32667;   Fc[5]= 0.68489;   Fll[5]= 0.943328; 
              cafac[5]= 0.817465;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2754.74;   wouttag[5]= 592.298;  
            fbb[5]= 0.123928;   fcc[5]= 0.192951;   fc[5]= 0.0924638;   fll[5]= 0.590657;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37932;   Fcc[6]= 1.37932;   Fc[6]= 0.712069;   Fll[6]= 0.980763; 
              cafac[6]= 0.815374;  
              winpretag[6]= 73989 ;   wintag[6]= 9196.82 ;   woutpretag[6]= 60328.7;   wouttag[6]= 8079.38;  
            fbb[6]= 0.0767222;   fcc[6]= 0.142264;   fc[6]= 0.11671;   fll[6]= 0.664304;  
            fmcbb[6]= 0.0556233;   fmccc[6]= 0.103141;   fmcc[6]= 0.163902;   fmcll[6]= 0.677333;  
           Fbb[7]= 1.36023;   Fcc[7]= 1.36023;   Fc[7]= 0.702215;   Fll[7]= 0.96719; 
              cafac[7]= 0.841242;  
              winpretag[7]= 16127 ;   wintag[7]= 2627.73 ;   woutpretag[7]= 13566.7;   wouttag[7]= 2431.39;  
            fbb[7]= 0.101279;   fcc[7]= 0.164841;   fc[7]= 0.116833;   fll[7]= 0.617048;  
            fmcbb[7]= 0.0744572;   fmccc[7]= 0.121186;   fmcc[7]= 0.166377;   fmcll[7]= 0.63798;  
      }  
      else if ( idsys==46   || sysname=="Wcjet4_down" ) { 
           Fbb[1]= 1.43394;   Fcc[1]= 1.43394;   Fc[1]= 0.740266;   Fll[1]= 1.0196; 
              cafac[1]= 0.980017;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06117e+06;   wouttag[1]= 36033;  
            fbb[1]= 0.0172893;   fcc[1]= 0.0503027;   fc[1]= 0.103702;   fll[1]= 0.828706;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.877209;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17583.5;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38474;   Fcc[3]= 1.38474;   Fc[3]= 0.714865;   Fll[3]= 0.984614; 
              cafac[3]= 0.807815;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46741.7;   wouttag[3]= 5659.07;  
            fbb[3]= 0.0697545;   fcc[3]= 0.135859;   fc[3]= 0.116675;   fll[3]= 0.677712;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.34635;   Fcc[4]= 1.34635;   Fc[4]= 0.695047;   Fll[4]= 0.957318; 
              cafac[4]= 0.807911;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1902.31 ;   woutpretag[4]= 10306.7;   wouttag[4]= 1707.71;  
            fbb[4]= 0.0989627;   fcc[4]= 0.163554;   fc[4]= 0.0879103;   fll[4]= 0.649573;  
            fmcbb[4]= 0.0735047;   fmccc[4]= 0.12148;   fmcc[4]= 0.126481;   fmcll[4]= 0.678535;  
           Fbb[5]= 1.32667;   Fcc[5]= 1.32667;   Fc[5]= 0.68489;   Fll[5]= 0.943328; 
              cafac[5]= 0.817465;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2754.74;   wouttag[5]= 592.298;  
            fbb[5]= 0.123928;   fcc[5]= 0.192951;   fc[5]= 0.0924638;   fll[5]= 0.590657;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37523;   Fcc[6]= 1.37523;   Fc[6]= 0.70996;   Fll[6]= 0.977858; 
              cafac[6]= 0.807759;  
              winpretag[6]= 73989 ;   wintag[6]= 9117.17 ;   woutpretag[6]= 59765.3;   wouttag[6]= 7970.57;  
            fbb[6]= 0.0774563;   fcc[6]= 0.143432;   fc[6]= 0.110466;   fll[6]= 0.668646;  
            fmcbb[6]= 0.0563223;   fmccc[6]= 0.104296;   fmcc[6]= 0.155595;   fmcll[6]= 0.683787;  
           Fbb[7]= 1.34219;   Fcc[7]= 1.34219;   Fc[7]= 0.6929;   Fll[7]= 0.95436; 
              cafac[7]= 0.809856;  
              winpretag[7]= 16127 ;   wintag[7]= 2548.07 ;   woutpretag[7]= 13060.6;   wouttag[7]= 2300.18;  
            fbb[7]= 0.10424;   fcc[7]= 0.169768;   fc[7]= 0.0888729;   fll[7]= 0.637118;  
            fmcbb[7]= 0.0776645;   fmccc[7]= 0.126486;   fmcc[7]= 0.128262;   fmcll[7]= 0.667587;  
      }  
      else if ( idsys==75   || sysname=="Wcjet5_up" ) { 
           Fbb[1]= 1.43394;   Fcc[1]= 1.43394;   Fc[1]= 0.740266;   Fll[1]= 1.0196; 
              cafac[1]= 0.980017;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06117e+06;   wouttag[1]= 36033;  
            fbb[1]= 0.0172893;   fcc[1]= 0.0503027;   fc[1]= 0.103702;   fll[1]= 0.828706;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.877209;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17583.5;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38474;   Fcc[3]= 1.38474;   Fc[3]= 0.714865;   Fll[3]= 0.984614; 
              cafac[3]= 0.807815;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46741.7;   wouttag[3]= 5659.07;  
            fbb[3]= 0.0697545;   fcc[3]= 0.135859;   fc[3]= 0.116675;   fll[3]= 0.677712;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35777;   Fcc[4]= 1.35777;   Fc[4]= 0.700942;   Fll[4]= 0.965438; 
              cafac[4]= 0.828149;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10564.8;   wouttag[4]= 1770.22;  
            fbb[4]= 0.0970496;   fcc[4]= 0.160392;   fc[4]= 0.105543;   fll[4]= 0.637016;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.35194;   Fcc[5]= 1.35194;   Fc[5]= 0.697934;   Fll[5]= 0.961293; 
              cafac[5]= 0.868538;  
              winpretag[5]= 3369.85 ;   wintag[5]= 663.538 ;   woutpretag[5]= 2926.85;   wouttag[5]= 637.246;  
            fbb[5]= 0.118798;   fcc[5]= 0.184964;   fc[5]= 0.13003;   fll[5]= 0.566208;  
            fmcbb[5]= 0.0878722;   fmccc[5]= 0.136814;   fmcc[5]= 0.186307;   fmcll[5]= 0.589006;  
           Fbb[6]= 1.37849;   Fcc[6]= 1.37849;   Fc[6]= 0.711642;   Fll[6]= 0.980174; 
              cafac[6]= 0.814399;  
              winpretag[6]= 73989 ;   wintag[6]= 9174.77 ;   woutpretag[6]= 60256.6;   wouttag[6]= 8056.14;  
            fbb[6]= 0.0768102;   fcc[6]= 0.142434;   fc[6]= 0.115347;   fll[6]= 0.66541;  
            fmcbb[6]= 0.0557205;   fmccc[6]= 0.103326;   fmcc[6]= 0.162085;   fmcll[6]= 0.678869;  
           Fbb[7]= 1.35654;   Fcc[7]= 1.35654;   Fc[7]= 0.700312;   Fll[7]= 0.964569; 
              cafac[7]= 0.837142;  
              winpretag[7]= 16127 ;   wintag[7]= 2605.68 ;   woutpretag[7]= 13500.6;   wouttag[7]= 2404.14;  
            fbb[7]= 0.10161;   fcc[7]= 0.165544;   fc[7]= 0.110677;   fll[7]= 0.62217;  
            fmcbb[7]= 0.0749032;   fmccc[7]= 0.122034;   fmcc[7]= 0.15804;   fmcll[7]= 0.645023;  
      }  
      else if ( idsys==76   || sysname=="Wcjet5_down" ) { 
           Fbb[1]= 1.43394;   Fcc[1]= 1.43394;   Fc[1]= 0.740266;   Fll[1]= 1.0196; 
              cafac[1]= 0.980017;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06117e+06;   wouttag[1]= 36033;  
            fbb[1]= 0.0172893;   fcc[1]= 0.0503027;   fc[1]= 0.103702;   fll[1]= 0.828706;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.877209;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225544;   wouttag[2]= 17583.5;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38474;   Fcc[3]= 1.38474;   Fc[3]= 0.714865;   Fll[3]= 0.984614; 
              cafac[3]= 0.807815;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46741.7;   wouttag[3]= 5659.07;  
            fbb[3]= 0.0697545;   fcc[3]= 0.135859;   fc[3]= 0.116675;   fll[3]= 0.677712;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35777;   Fcc[4]= 1.35777;   Fc[4]= 0.700942;   Fll[4]= 0.965438; 
              cafac[4]= 0.828149;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10564.8;   wouttag[4]= 1770.22;  
            fbb[4]= 0.0970496;   fcc[4]= 0.160392;   fc[4]= 0.105543;   fll[4]= 0.637016;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.30233;   Fcc[5]= 1.30233;   Fc[5]= 0.672325;   Fll[5]= 0.926022; 
              cafac[5]= 0.773642;  
              winpretag[5]= 3369.85 ;   wintag[5]= 627.984 ;   woutpretag[5]= 2607.06;   wouttag[5]= 553.729;  
            fbb[5]= 0.128869;   fcc[5]= 0.200645;   fc[5]= 0.0562758;   fll[5]= 0.61421;  
            fmcbb[5]= 0.0989526;   fmccc[5]= 0.154066;   fmcc[5]= 0.0837033;   fmcll[5]= 0.663278;  
           Fbb[6]= 1.37606;   Fcc[6]= 1.37606;   Fc[6]= 0.710385;   Fll[6]= 0.978444; 
              cafac[6]= 0.808712;  
              winpretag[6]= 73989 ;   wintag[6]= 9139.22 ;   woutpretag[6]= 59835.8;   wouttag[6]= 7993.41;  
            fbb[6]= 0.077369;   fcc[6]= 0.143263;   fc[6]= 0.111823;   fll[6]= 0.667544;  
            fmcbb[6]= 0.0562251;   fmccc[6]= 0.104112;   fmcc[6]= 0.157412;   fmcll[6]= 0.682251;  
           Fbb[7]= 1.3458;   Fcc[7]= 1.3458;   Fc[7]= 0.694763;   Fll[7]= 0.956927; 
              cafac[7]= 0.813579;  
              winpretag[7]= 16127 ;   wintag[7]= 2570.12 ;   woutpretag[7]= 13120.6;   wouttag[7]= 2325.49;  
            fbb[7]= 0.10392;   fcc[7]= 0.169084;   fc[7]= 0.0949045;   fll[7]= 0.632091;  
            fmcbb[7]= 0.0772186;   fmccc[7]= 0.125639;   fmcc[7]= 0.1366;   fmcll[7]= 0.660543;  
      }  
      else if ( idsys==3   || sysname=="wt_up" ) { 
           Fbb[1]= 1.42768;   Fcc[1]= 1.42768;   Fc[1]= 0.661809;   Fll[1]= 1.03349; 
              cafac[1]= 0.966097;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.0461e+06;   wouttag[1]= 33863.3;  
            fbb[1]= 0.0172139;   fcc[1]= 0.0500832;   fc[1]= 0.0927112;   fll[1]= 0.839992;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.40863;   Fcc[2]= 1.40863;   Fc[2]= 0.652977;   Fll[2]= 1.01969; 
              cafac[2]= 0.862712;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 221817;   wouttag[2]= 16789.6;  
            fbb[2]= 0.0424002;   fcc[2]= 0.101727;   fc[2]= 0.105926;   fll[2]= 0.749947;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38428;   Fcc[3]= 1.38428;   Fc[3]= 0.641688;   Fll[3]= 1.00207; 
              cafac[3]= 0.793385;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 45906.8;   wouttag[3]= 5445.61;  
            fbb[3]= 0.0697314;   fcc[3]= 0.135814;   fc[3]= 0.104732;   fll[3]= 0.689723;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35677;   Fcc[4]= 1.35677;   Fc[4]= 0.628938;   Fll[4]= 0.982155; 
              cafac[4]= 0.813422;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10377;   wouttag[4]= 1714.25;  
            fbb[4]= 0.0969786;   fcc[4]= 0.160274;   fc[4]= 0.0947009;   fll[4]= 0.648046;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32504;   Fcc[5]= 1.32504;   Fc[5]= 0.614229;   Fll[5]= 0.959185; 
              cafac[5]= 0.801336;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2700.39;   wouttag[5]= 574.543;  
            fbb[5]= 0.123775;   fcc[5]= 0.192714;   fc[5]= 0.0829242;   fll[5]= 0.600586;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37666;   Fcc[6]= 1.37666;   Fc[6]= 0.638158;   Fll[6]= 0.996553; 
              cafac[6]= 0.796967;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 58966.8;   wouttag[6]= 7737.98;  
            fbb[6]= 0.0770556;   fcc[6]= 0.142786;   fc[6]= 0.101945;   fll[6]= 0.678214;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.35002;   Fcc[7]= 1.35002;   Fc[7]= 0.625807;   Fll[7]= 0.977265; 
              cafac[7]= 0.810095;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13064.4;   wouttag[7]= 2290.51;  
            fbb[7]= 0.102683;   fcc[7]= 0.167181;   fc[7]= 0.0921937;   fll[7]= 0.637942;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==4   || sysname=="wt_down" ) { 
           Fbb[1]= 1.44002;   Fcc[1]= 1.44002;   Fc[1]= 0.816565;   Fll[1]= 1.0061; 
              cafac[1]= 0.993944;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.07625e+06;   wouttag[1]= 38203.7;  
            fbb[1]= 0.0173626;   fcc[1]= 0.0505161;   fc[1]= 0.114391;   fll[1]= 0.817731;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41213;   Fcc[2]= 1.41213;   Fc[2]= 0.800751;   Fll[2]= 0.986613; 
              cafac[2]= 0.891692;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 229268;   wouttag[2]= 18376.6;  
            fbb[2]= 0.0425057;   fcc[2]= 0.101979;   fc[2]= 0.129898;   fll[2]= 0.725617;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38518;   Fcc[3]= 1.38518;   Fc[3]= 0.785465;   Fll[3]= 0.967778; 
              cafac[3]= 0.822242;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 47576.5;   wouttag[3]= 5872.5;  
            fbb[3]= 0.0697767;   fcc[3]= 0.135902;   fc[3]= 0.128198;   fll[3]= 0.666124;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35872;   Fcc[4]= 1.35872;   Fc[4]= 0.770465;   Fll[4]= 0.949297; 
              cafac[4]= 0.842881;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10752.8;   wouttag[4]= 1826.21;  
            fbb[4]= 0.097118;   fcc[4]= 0.160505;   fc[4]= 0.116011;   fll[4]= 0.626366;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32824;   Fcc[5]= 1.32824;   Fc[5]= 0.753182;   Fll[5]= 0.928003; 
              cafac[5]= 0.833652;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2809.29;   wouttag[5]= 610.127;  
            fbb[5]= 0.124075;   fcc[5]= 0.19318;   fc[5]= 0.101684;   fll[5]= 0.581062;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37786;   Fcc[6]= 1.37786;   Fc[6]= 0.781317;   Fll[6]= 0.962668; 
              cafac[6]= 0.826122;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 61123.9;   wouttag[6]= 8311.36;  
            fbb[6]= 0.0771227;   fcc[6]= 0.14291;   fc[6]= 0.124814;   fll[6]= 0.655153;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.35224;   Fcc[7]= 1.35224;   Fc[7]= 0.766789;   Fll[7]= 0.944767; 
              cafac[7]= 0.840214;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13550.2;   wouttag[7]= 2437.78;  
            fbb[7]= 0.102853;   fcc[7]= 0.167456;   fc[7]= 0.112963;   fll[7]= 0.616728;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==5   || sysname=="st_up" ) { 
           Fbb[1]= 1.40954;   Fcc[1]= 1.40954;   Fc[1]= 0.741613;   Fll[1]= 1.02078; 
              cafac[1]= 0.980206;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06138e+06;   wouttag[1]= 35945.4;  
            fbb[1]= 0.0169952;   fcc[1]= 0.0494469;   fc[1]= 0.103891;   fll[1]= 0.829667;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.38834;   Fcc[2]= 1.38834;   Fc[2]= 0.730457;   Fll[2]= 1.00543; 
              cafac[2]= 0.876936;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225474;   wouttag[2]= 17484.9;  
            fbb[2]= 0.0417895;   fcc[2]= 0.100261;   fc[2]= 0.118495;   fll[2]= 0.739454;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.36459;   Fcc[3]= 1.36459;   Fc[3]= 0.717961;   Fll[3]= 0.988227; 
              cafac[3]= 0.806637;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46673.6;   wouttag[3]= 5614.07;  
            fbb[3]= 0.0687396;   fcc[3]= 0.133882;   fc[3]= 0.11718;   fll[3]= 0.680198;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.33928;   Fcc[4]= 1.33928;   Fc[4]= 0.704648;   Fll[4]= 0.969903; 
              cafac[4]= 0.826393;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10542.4;   wouttag[4]= 1754.19;  
            fbb[4]= 0.0957286;   fcc[4]= 0.158209;   fc[4]= 0.106101;   fll[4]= 0.639962;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.31004;   Fcc[5]= 1.31004;   Fc[5]= 0.689264;   Fll[5]= 0.948727; 
              cafac[5]= 0.815896;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2749.45;   wouttag[5]= 586.943;  
            fbb[5]= 0.122374;   fcc[5]= 0.190533;   fc[5]= 0.0930543;   fll[5]= 0.594038;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.35759;   Fcc[6]= 1.35759;   Fc[6]= 0.714279;   Fll[6]= 0.98316; 
              cafac[6]= 0.810265;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 59950.7;   wouttag[6]= 7957.9;  
            fbb[6]= 0.0759881;   fcc[6]= 0.140808;   fc[6]= 0.114105;   fll[6]= 0.669099;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.33307;   Fcc[7]= 1.33307;   Fc[7]= 0.701377;   Fll[7]= 0.9654; 
              cafac[7]= 0.823463;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13280;   wouttag[7]= 2342.65;  
            fbb[7]= 0.101394;   fcc[7]= 0.165082;   fc[7]= 0.103327;   fll[7]= 0.630197;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==6   || sysname=="st_down" ) { 
           Fbb[1]= 1.45611;   Fcc[1]= 1.45611;   Fc[1]= 0.739043;   Fll[1]= 1.01853; 
              cafac[1]= 0.979845;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06099e+06;   wouttag[1]= 36112.5;  
            fbb[1]= 0.0175566;   fcc[1]= 0.0510804;   fc[1]= 0.103531;   fll[1]= 0.827832;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.43042;   Fcc[2]= 1.43042;   Fc[2]= 0.726004;   Fll[2]= 1.00056; 
              cafac[2]= 0.877458;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225608;   wouttag[2]= 17672.9;  
            fbb[2]= 0.043056;   fcc[2]= 0.1033;   fc[2]= 0.117773;   fll[2]= 0.735871;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.40296;   Fcc[3]= 1.40296;   Fc[3]= 0.712066;   Fll[3]= 0.981347; 
              cafac[3]= 0.808882;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46803.5;   wouttag[3]= 5699.89;  
            fbb[3]= 0.0706724;   fcc[3]= 0.137646;   fc[3]= 0.116218;   fll[3]= 0.675463;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.37445;   Fcc[4]= 1.37445;   Fc[4]= 0.697598;   Fll[4]= 0.961407; 
              cafac[4]= 0.829742;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10585.2;   wouttag[4]= 1784.76;  
            fbb[4]= 0.0982421;   fcc[4]= 0.162363;   fc[4]= 0.105039;   fll[4]= 0.634356;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.34165;   Fcc[5]= 1.34165;   Fc[5]= 0.68095;   Fll[5]= 0.938464; 
              cafac[5]= 0.818892;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2759.55;   wouttag[5]= 597.145;  
            fbb[5]= 0.125327;   fcc[5]= 0.19513;   fc[5]= 0.0919319;   fll[5]= 0.587612;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.39506;   Fcc[6]= 1.39506;   Fc[6]= 0.708061;   Fll[6]= 0.975827; 
              cafac[6]= 0.812703;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60131.1;   wouttag[6]= 8085.15;  
            fbb[6]= 0.0780856;   fcc[6]= 0.144694;   fc[6]= 0.113112;   fll[6]= 0.664108;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.36746;   Fcc[7]= 1.36746;   Fc[7]= 0.694052;   Fll[7]= 0.95652; 
              cafac[7]= 0.826673;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13331.8;   wouttag[7]= 2383.54;  
            fbb[7]= 0.104011;   fcc[7]= 0.169341;   fc[7]= 0.102248;   fll[7]= 0.624401;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==7   || sysname=="z_up" ) { 
           Fbb[1]= 1.43175;   Fcc[1]= 1.43175;   Fc[1]= 0.677275;   Fll[1]= 1.03058; 
              cafac[1]= 0.966801;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.04686e+06;   wouttag[1]= 34229.1;  
            fbb[1]= 0.017263;   fcc[1]= 0.050226;   fc[1]= 0.0948778;   fll[1]= 0.837633;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41154;   Fcc[2]= 1.41154;   Fc[2]= 0.667715;   Fll[2]= 1.01604; 
              cafac[2]= 0.86333;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 221975;   wouttag[2]= 16912.8;  
            fbb[2]= 0.0424879;   fcc[2]= 0.101937;   fc[2]= 0.108317;   fll[2]= 0.747258;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.3867;   Fcc[3]= 1.3867;   Fc[3]= 0.655964;   Fll[3]= 0.998157; 
              cafac[3]= 0.796577;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46091.5;   wouttag[3]= 5494.29;  
            fbb[3]= 0.0698536;   fcc[3]= 0.136052;   fc[3]= 0.107062;   fll[3]= 0.687033;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35911;   Fcc[4]= 1.35911;   Fc[4]= 0.642911;   Fll[4]= 0.978295; 
              cafac[4]= 0.813652;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10379.9;   wouttag[4]= 1721.04;  
            fbb[4]= 0.0971456;   fcc[4]= 0.16055;   fc[4]= 0.0968048;   fll[4]= 0.645499;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32729;   Fcc[5]= 1.32729;   Fc[5]= 0.627859;   Fll[5]= 0.95539; 
              cafac[5]= 0.811297;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2733.95;   wouttag[5]= 583.398;  
            fbb[5]= 0.123985;   fcc[5]= 0.193041;   fc[5]= 0.0847643;   fll[5]= 0.59821;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37906;   Fcc[6]= 1.37906;   Fc[6]= 0.652351;   Fll[6]= 0.992658; 
              cafac[6]= 0.799996;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 59190.9;   wouttag[6]= 7802.18;  
            fbb[6]= 0.07719;   fcc[6]= 0.143035;   fc[6]= 0.104212;   fll[6]= 0.675563;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.35233;   Fcc[7]= 1.35233;   Fc[7]= 0.639707;   Fll[7]= 0.973418; 
              cafac[7]= 0.812604;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13104.9;   wouttag[7]= 2305.61;  
            fbb[7]= 0.10286;   fcc[7]= 0.167468;   fc[7]= 0.0942414;   fll[7]= 0.635431;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==8   || sysname=="z_down" ) { 
           Fbb[1]= 1.43607;   Fcc[1]= 1.43607;   Fc[1]= 0.801577;   Fll[1]= 1.00891; 
              cafac[1]= 0.993232;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.07548e+06;   wouttag[1]= 37836.7;  
            fbb[1]= 0.017315;   fcc[1]= 0.0503773;   fc[1]= 0.112291;   fll[1]= 0.820017;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.40931;   Fcc[2]= 1.40931;   Fc[2]= 0.786646;   Fll[2]= 0.990116; 
              cafac[2]= 0.891089;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 229113;   wouttag[2]= 18254.2;  
            fbb[2]= 0.0424208;   fcc[2]= 0.101776;   fc[2]= 0.12761;   fll[2]= 0.728193;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38283;   Fcc[3]= 1.38283;   Fc[3]= 0.771863;   Fll[3]= 0.97151; 
              cafac[3]= 0.818997;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 47388.8;   wouttag[3]= 5823.02;  
            fbb[3]= 0.0696586;   fcc[3]= 0.135672;   fc[3]= 0.125978;   fll[3]= 0.668692;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35646;   Fcc[4]= 1.35646;   Fc[4]= 0.757146;   Fll[4]= 0.952986; 
              cafac[4]= 0.842678;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10750.2;   wouttag[4]= 1819.53;  
            fbb[4]= 0.0969566;   fcc[4]= 0.160238;   fc[4]= 0.114005;   fll[4]= 0.6288;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32607;   Fcc[5]= 1.32607;   Fc[5]= 0.740182;   Fll[5]= 0.931634; 
              cafac[5]= 0.823482;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2775.01;   wouttag[5]= 601.017;  
            fbb[5]= 0.123872;   fcc[5]= 0.192864;   fc[5]= 0.0999285;   fll[5]= 0.583335;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37554;   Fcc[6]= 1.37554;   Fc[6]= 0.767793;   Fll[6]= 0.966387; 
              cafac[6]= 0.823045;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60896.3;   wouttag[6]= 8246.19;  
            fbb[6]= 0.0769928;   fcc[6]= 0.142669;   fc[6]= 0.122654;   fll[6]= 0.657684;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.35;   Fcc[7]= 1.35;   Fc[7]= 0.753537;   Fll[7]= 0.948444; 
              cafac[7]= 0.837674;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13509.2;   wouttag[7]= 2422.51;  
            fbb[7]= 0.102682;   fcc[7]= 0.167179;   fc[7]= 0.111011;   fll[7]= 0.619128;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==9   || sysname=="di_up" ) { 
           Fbb[1]= 1.43263;   Fcc[1]= 1.43263;   Fc[1]= 0.739075;   Fll[1]= 1.01988; 
              cafac[1]= 0.97977;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06091e+06;   wouttag[1]= 35992.2;  
            fbb[1]= 0.0172735;   fcc[1]= 0.0502566;   fc[1]= 0.103535;   fll[1]= 0.828935;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.40929;   Fcc[2]= 1.40929;   Fc[2]= 0.727035;   Fll[2]= 1.00327; 
              cafac[2]= 0.876809;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225441;   wouttag[2]= 17562.7;  
            fbb[2]= 0.04242;   fcc[2]= 0.101774;   fc[2]= 0.11794;   fll[2]= 0.737866;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38373;   Fcc[3]= 1.38373;   Fc[3]= 0.713849;   Fll[3]= 0.985072; 
              cafac[3]= 0.807345;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46714.6;   wouttag[3]= 5652.13;  
            fbb[3]= 0.0697038;   fcc[3]= 0.13576;   fc[3]= 0.116509;   fll[3]= 0.678027;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35683;   Fcc[4]= 1.35683;   Fc[4]= 0.699975;   Fll[4]= 0.965926; 
              cafac[4]= 0.827713;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10559.3;   wouttag[4]= 1768.28;  
            fbb[4]= 0.096983;   fcc[4]= 0.160282;   fc[4]= 0.105397;   fll[4]= 0.637338;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32582;   Fcc[5]= 1.32582;   Fc[5]= 0.683976;   Fll[5]= 0.943849; 
              cafac[5]= 0.817249;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2754.01;   wouttag[5]= 591.833;  
            fbb[5]= 0.123848;   fcc[5]= 0.192828;   fc[5]= 0.0923404;   fll[5]= 0.590984;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37629;   Fcc[6]= 1.37629;   Fc[6]= 0.71001;   Fll[6]= 0.979775; 
              cafac[6]= 0.811095;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60012.1;   wouttag[6]= 8015.2;  
            fbb[6]= 0.0770346;   fcc[6]= 0.142747;   fc[6]= 0.113423;   fll[6]= 0.666796;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.35023;   Fcc[7]= 1.35023;   Fc[7]= 0.69657;   Fll[7]= 0.961228; 
              cafac[7]= 0.824762;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13301;   wouttag[7]= 2361.68;  
            fbb[7]= 0.1027;   fcc[7]= 0.167208;   fc[7]= 0.102619;   fll[7]= 0.627474;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==10   || sysname=="di_down" ) { 
           Fbb[1]= 1.43525;   Fcc[1]= 1.43525;   Fc[1]= 0.741458;   Fll[1]= 1.01932; 
              cafac[1]= 0.980263;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06144e+06;   wouttag[1]= 36073.7;  
            fbb[1]= 0.0173051;   fcc[1]= 0.0503487;   fc[1]= 0.103869;   fll[1]= 0.828477;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41153;   Fcc[2]= 1.41153;   Fc[2]= 0.729207;   Fll[2]= 1.00248; 
              cafac[2]= 0.87761;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 225647;   wouttag[2]= 17604.3;  
            fbb[2]= 0.0424877;   fcc[2]= 0.101936;   fc[2]= 0.118292;   fll[2]= 0.737284;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38574;   Fcc[3]= 1.38574;   Fc[3]= 0.715881;   Fll[3]= 0.984157; 
              cafac[3]= 0.808284;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46768.9;   wouttag[3]= 5666.01;  
            fbb[3]= 0.0698052;   fcc[3]= 0.135957;   fc[3]= 0.116841;   fll[3]= 0.677397;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.3587;   Fcc[4]= 1.3587;   Fc[4]= 0.70191;   Fll[4]= 0.96495; 
              cafac[4]= 0.828586;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10570.4;   wouttag[4]= 1772.16;  
            fbb[4]= 0.0971161;   fcc[4]= 0.160502;   fc[4]= 0.105688;   fll[4]= 0.636694;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32752;   Fcc[5]= 1.32752;   Fc[5]= 0.685803;   Fll[5]= 0.942808; 
              cafac[5]= 0.817681;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2755.47;   wouttag[5]= 592.762;  
            fbb[5]= 0.124007;   fcc[5]= 0.193075;   fc[5]= 0.0925871;   fll[5]= 0.590332;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37826;   Fcc[6]= 1.37826;   Fc[6]= 0.712015;   Fll[6]= 0.978842; 
              cafac[6]= 0.811991;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60078.4;   wouttag[6]= 8034.07;  
            fbb[6]= 0.0771449;   fcc[6]= 0.142951;   fc[6]= 0.113743;   fll[6]= 0.666161;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.35206;   Fcc[7]= 1.35206;   Fc[7]= 0.698482;   Fll[7]= 0.960238; 
              cafac[7]= 0.825529;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13313.3;   wouttag[7]= 2366.52;  
            fbb[7]= 0.102839;   fcc[7]= 0.167434;   fc[7]= 0.1029;   fll[7]= 0.626827;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==11   || sysname=="qcd_up" ) { 
           Fbb[1]= 1.42509;   Fcc[1]= 1.42509;   Fc[1]= 0.495375;   Fll[1]= 1.06232; 
              cafac[1]= 0.937798;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.01546e+06;   wouttag[1]= 29509;  
            fbb[1]= 0.0171826;   fcc[1]= 0.0499922;   fc[1]= 0.0693959;   fll[1]= 0.863429;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41452;   Fcc[2]= 1.41452;   Fc[2]= 0.491702;   Fll[2]= 1.05445; 
              cafac[2]= 0.833153;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 214217;   wouttag[2]= 15216.6;  
            fbb[2]= 0.0425776;   fcc[2]= 0.102152;   fc[2]= 0.0797641;   fll[2]= 0.775506;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.39216;   Fcc[3]= 1.39216;   Fc[3]= 0.483928;   Fll[3]= 1.03777; 
              cafac[3]= 0.76628;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 44338.5;   wouttag[3]= 5042.97;  
            fbb[3]= 0.0701284;   fcc[3]= 0.136587;   fc[3]= 0.0789831;   fll[3]= 0.714302;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.36276;   Fcc[4]= 1.36276;   Fc[4]= 0.473708;   Fll[4]= 1.01586; 
              cafac[4]= 0.789406;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10070.6;   wouttag[4]= 1618.09;  
            fbb[4]= 0.0974064;   fcc[4]= 0.160982;   fc[4]= 0.0713275;   fll[4]= 0.670285;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32884;   Fcc[5]= 1.32884;   Fc[5]= 0.461919;   Fll[5]= 0.990576; 
              cafac[5]= 0.781778;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2634.48;   wouttag[5]= 549.687;  
            fbb[5]= 0.12413;   fcc[5]= 0.193267;   fc[5]= 0.0623615;   fll[5]= 0.620241;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.38401;   Fcc[6]= 1.38401;   Fc[6]= 0.481094;   Fll[6]= 1.0317; 
              cafac[6]= 0.770824;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 57032.5;   wouttag[6]= 7212.29;  
            fbb[6]= 0.0774666;   fcc[6]= 0.143547;   fc[6]= 0.0768541;   fll[6]= 0.702132;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.35553;   Fcc[7]= 1.35553;   Fc[7]= 0.471195;   Fll[7]= 1.01047; 
              cafac[7]= 0.787117;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 12693.9;   wouttag[7]= 2169.39;  
            fbb[7]= 0.103103;   fcc[7]= 0.167863;   fc[7]= 0.0694164;   fll[7]= 0.659618;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==12   || sysname=="qcd_down" ) { 
           Fbb[1]= 1.44208;   Fcc[1]= 1.44208;   Fc[1]= 0.965625;   Fll[1]= 0.980286; 
              cafac[1]= 1.02237;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.10704e+06;   wouttag[1]= 42578;  
            fbb[1]= 0.0173875;   fcc[1]= 0.0505884;   fc[1]= 0.135272;   fll[1]= 0.796752;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.40669;   Fcc[2]= 1.40669;   Fc[2]= 0.941928;   Fll[2]= 0.95623; 
              cafac[2]= 0.921266;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 236872;   wouttag[2]= 19950.4;  
            fbb[2]= 0.0423419;   fcc[2]= 0.101587;   fc[2]= 0.1528;   fll[2]= 0.703271;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.37805;   Fcc[3]= 1.37805;   Fc[3]= 0.922751;   Fll[3]= 0.936761; 
              cafac[3]= 0.849251;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 49139.4;   wouttag[3]= 6273.73;  
            fbb[3]= 0.069418;   fcc[3]= 0.135203;   fc[3]= 0.150605;   fll[3]= 0.644774;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35326;   Fcc[4]= 1.35326;   Fc[4]= 0.906147;   Fll[4]= 0.919905; 
              cafac[4]= 0.866556;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 11054.8;   wouttag[4]= 1921.03;  
            fbb[4]= 0.0967273;   fcc[4]= 0.159859;   fc[4]= 0.136441;   fll[4]= 0.606973;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.3247;   Fcc[5]= 1.3247;   Fc[5]= 0.887027;   Fll[5]= 0.900495; 
              cafac[5]= 0.852756;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2873.66;   wouttag[5]= 634.434;  
            fbb[5]= 0.123744;   fcc[5]= 0.192665;   fc[5]= 0.119753;   fll[5]= 0.563838;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37121;   Fcc[6]= 1.37121;   Fc[6]= 0.918166;   Fll[6]= 0.932107; 
              cafac[6]= 0.852098;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 63045.9;   wouttag[6]= 8833.72;  
            fbb[6]= 0.0767502;   fcc[6]= 0.14222;   fc[6]= 0.146676;   fll[6]= 0.634354;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.34719;   Fcc[7]= 1.34719;   Fc[7]= 0.902084;   Fll[7]= 0.915781; 
              cafac[7]= 0.862822;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13914.8;   wouttag[7]= 2557;  
            fbb[7]= 0.102468;   fcc[7]= 0.166831;   fc[7]= 0.132895;   fll[7]= 0.597806;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==13   || sysname=="jes_up" ) { 
           Fbb[1]= 1.46822;   Fcc[1]= 1.46822;   Fc[1]= 0.770737;   Fll[1]= 1.01288; 
              cafac[1]= 0.897635;  
              winpretag[1]= 1.16655e+06 ;   wintag[1]= 42032.9 ;   woutpretag[1]= 1.04714e+06;   wouttag[1]= 35299.8;  
            fbb[1]= 0.016423;   fcc[1]= 0.0483802;   fc[1]= 0.104963;   fll[1]= 0.830234;  
            fmcbb[1]= 0.0111856;   fmccc[1]= 0.0329515;   fmcc[1]= 0.136185;   fmcll[1]= 0.819678;  
           Fbb[2]= 1.44201;   Fcc[2]= 1.44201;   Fc[2]= 0.756978;   Fll[2]= 0.994797; 
              cafac[2]= 0.802066;  
              winpretag[2]= 286191 ;   wintag[2]= 21796.4 ;   woutpretag[2]= 229544;   wouttag[2]= 17923.8;  
            fbb[2]= 0.0407415;   fcc[2]= 0.100268;   fc[2]= 0.122636;   fll[2]= 0.736355;  
            fmcbb[2]= 0.0282532;   fmccc[2]= 0.0695331;   fmcc[2]= 0.162007;   fmcll[2]= 0.740207;  
           Fbb[3]= 1.41413;   Fcc[3]= 1.41413;   Fc[3]= 0.74234;   Fll[3]= 0.975561; 
              cafac[3]= 0.716057;  
              winpretag[3]= 67786.3 ;   wintag[3]= 7582.04 ;   woutpretag[3]= 48538.8;   wouttag[3]= 5869.58;  
            fbb[3]= 0.067748;   fcc[3]= 0.134747;   fc[3]= 0.122103;   fll[3]= 0.675402;  
            fmcbb[3]= 0.047908;   fmccc[3]= 0.0952862;   fmcc[3]= 0.164484;   fmcll[3]= 0.692322;  
           Fbb[4]= 1.38473;   Fcc[4]= 1.38473;   Fc[4]= 0.72691;   Fll[4]= 0.955283; 
              cafac[4]= 0.675252;  
              winpretag[4]= 15553.4 ;   wintag[4]= 2354.07 ;   woutpretag[4]= 10502.5;   wouttag[4]= 1763.59;  
            fbb[4]= 0.0967822;   fcc[4]= 0.159883;   fc[4]= 0.111033;   fll[4]= 0.632302;  
            fmcbb[4]= 0.0698923;   fmccc[4]= 0.115461;   fmcc[4]= 0.152746;   fmcll[4]= 0.6619;  
           Fbb[5]= 1.35349;   Fcc[5]= 1.35349;   Fc[5]= 0.710506;   Fll[5]= 0.933725; 
              cafac[5]= 0.619059;  
              winpretag[5]= 4385.19 ;   wintag[5]= 843.918 ;   woutpretag[5]= 2714.7;   wouttag[5]= 589.686;  
            fbb[5]= 0.122065;   fcc[5]= 0.19292;   fc[5]= 0.0999856;   fll[5]= 0.585029;  
            fmcbb[5]= 0.0901855;   fmccc[5]= 0.142536;   fmcc[5]= 0.140724;   fmcll[5]= 0.626554;  
           Fbb[6]= 1.40569;   Fcc[6]= 1.40569;   Fc[6]= 0.73791;   Fll[6]= 0.969739; 
              cafac[6]= 0.701766;  
              winpretag[6]= 87724.9 ;   wintag[6]= 10780 ;   woutpretag[6]= 61562.4;   wouttag[6]= 8272.94;  
            fbb[6]= 0.0757935;   fcc[6]= 0.142291;   fc[6]= 0.118962;   fll[6]= 0.662953;  
            fmcbb[6]= 0.0539191;   fmccc[6]= 0.101225;   fmcc[6]= 0.161215;   fmcll[6]= 0.683641;  
           Fbb[7]= 1.37774;   Fcc[7]= 1.37774;   Fc[7]= 0.723238;   Fll[7]= 0.950456; 
              cafac[7]= 0.661161;  
              winpretag[7]= 19938.6 ;   wintag[7]= 3197.99 ;   woutpretag[7]= 13182.6;   wouttag[7]= 2359.14;  
            fbb[7]= 0.102442;   fcc[7]= 0.167279;   fc[7]= 0.10856;   fll[7]= 0.621719;  
            fmcbb[7]= 0.0743554;   fmccc[7]= 0.121416;   fmcc[7]= 0.150102;   fmcll[7]= 0.654126;  
      }  
      else if ( idsys==14   || sysname=="jes_down" ) { 
           Fbb[1]= 1.42224;   Fcc[1]= 1.42224;   Fc[1]= 0.659054;   Fll[1]= 1.03496; 
              cafac[1]= 1.03684;  
              winpretag[1]= 1.01238e+06 ;   wintag[1]= 38575.5 ;   woutpretag[1]= 1.04968e+06;   wouttag[1]= 34782.5;  
            fbb[1]= 0.0182461;   fcc[1]= 0.0527068;   fc[1]= 0.0951646;   fll[1]= 0.833883;  
            fmcbb[1]= 0.0128291;   fmccc[1]= 0.037059;   fmcc[1]= 0.144396;   fmcll[1]= 0.805716;  
           Fbb[2]= 1.4013;   Fcc[2]= 1.4013;   Fc[2]= 0.649352;   Fll[2]= 1.01972; 
              cafac[2]= 0.949768;  
              winpretag[2]= 232271 ;   wintag[2]= 18443.6 ;   woutpretag[2]= 220603;   wouttag[2]= 17209;  
            fbb[2]= 0.0445163;   fcc[2]= 0.105361;   fc[2]= 0.106132;   fll[2]= 0.743991;  
            fmcbb[2]= 0.0317678;   fmccc[2]= 0.0751879;   fmcc[2]= 0.163442;   fmcll[2]= 0.729602;  
           Fbb[3]= 1.37559;   Fcc[3]= 1.37559;   Fc[3]= 0.637435;   Fll[3]= 1.00101; 
              cafac[3]= 0.91133;  
              winpretag[3]= 49741.1 ;   wintag[3]= 5747.9 ;   woutpretag[3]= 45330.6;   wouttag[3]= 5486.35;  
            fbb[3]= 0.0726474;   fcc[3]= 0.138336;   fc[3]= 0.102494;   fll[3]= 0.686522;  
            fmcbb[3]= 0.052812;   fmccc[3]= 0.100565;   fmcc[3]= 0.160792;   fmcll[3]= 0.685831;  
           Fbb[4]= 1.34707;   Fcc[4]= 1.34707;   Fc[4]= 0.624219;   Fll[4]= 0.980254; 
              cafac[4]= 0.970491;  
              winpretag[4]= 10579.3 ;   wintag[4]= 1643.51 ;   woutpretag[4]= 10267.1;   wouttag[4]= 1742.44;  
            fbb[4]= 0.101845;   fcc[4]= 0.16475;   fc[4]= 0.0926583;   fll[4]= 0.640746;  
            fmcbb[4]= 0.0756051;   fmccc[4]= 0.122303;   fmcc[4]= 0.148439;   fmcll[4]= 0.653653;  
           Fbb[5]= 1.31608;   Fcc[5]= 1.31608;   Fc[5]= 0.609861;   Fll[5]= 0.957708; 
              cafac[5]= 1.15519;  
              winpretag[5]= 2580.82 ;   wintag[5]= 508.876 ;   woutpretag[5]= 2981.33;   wouttag[5]= 646.921;  
            fbb[5]= 0.123476;   fcc[5]= 0.200978;   fc[5]= 0.0807516;   fll[5]= 0.594794;  
            fmcbb[5]= 0.0938207;   fmccc[5]= 0.152709;   fmcc[5]= 0.13241;   fmcll[5]= 0.62106;  
           Fbb[6]= 1.36818;   Fcc[6]= 1.36818;   Fc[6]= 0.634001;   Fll[6]= 0.995616; 
              cafac[6]= 0.932818;  
              winpretag[6]= 62901.2 ;   wintag[6]= 7900.29 ;   woutpretag[6]= 58675.4;   wouttag[6]= 7829.57;  
            fbb[6]= 0.0798031;   fcc[6]= 0.14552;   fc[6]= 0.0998867;   fll[6]= 0.67479;  
            fmcbb[6]= 0.0583281;   fmccc[6]= 0.106361;   fmcc[6]= 0.15755;   fmcll[6]= 0.677761;  
           Fbb[7]= 1.34088;   Fcc[7]= 1.34088;   Fc[7]= 0.62135;   Fll[7]= 0.975749; 
              cafac[7]= 1.00836;  
              winpretag[7]= 13160.1 ;   wintag[7]= 2152.39 ;   woutpretag[7]= 13270.2;   wouttag[7]= 2377.45;  
            fbb[7]= 0.106167;   fcc[7]= 0.171989;   fc[7]= 0.0902793;   fll[7]= 0.631565;  
            fmcbb[7]= 0.0791774;   fmccc[7]= 0.128266;   fmcc[7]= 0.145295;   fmcll[7]= 0.647261;  
      }  
      else if ( idsys==79   || sysname=="leff_up" ) { 
           Fbb[1]= 1.41873;   Fcc[1]= 1.41873;   Fc[1]= 0.722812;   Fll[1]= 1.02348; 
              cafac[1]= 0.952807;  
              winpretag[1]= 1.10982e+06 ;   wintag[1]= 41109.4 ;   woutpretag[1]= 1.05744e+06;   wouttag[1]= 35461.2;  
            fbb[1]= 0.017108;   fcc[1]= 0.0497751;   fc[1]= 0.101235;   fll[1]= 0.831882;  
            fmcbb[1]= 0.0120587;   fmccc[1]= 0.0350843;   fmcc[1]= 0.140057;   fmcll[1]= 0.8128;  
           Fbb[2]= 1.39754;   Fcc[2]= 1.39754;   Fc[2]= 0.712018;   Fll[2]= 1.00819; 
              cafac[2]= 0.85211;  
              winpretag[2]= 263551 ;   wintag[2]= 20309.2 ;   woutpretag[2]= 224574;   wouttag[2]= 17337.5;  
            fbb[2]= 0.0420715;   fcc[2]= 0.100938;   fc[2]= 0.115477;   fll[2]= 0.741514;  
            fmcbb[2]= 0.0301039;   fmccc[2]= 0.0722252;   fmcc[2]= 0.162183;   fmcll[2]= 0.735488;  
           Fbb[3]= 1.37326;   Fcc[3]= 1.37326;   Fc[3]= 0.699646;   Fll[3]= 0.990674; 
              cafac[3]= 0.784176;  
              winpretag[3]= 59317 ;   wintag[3]= 6733.68 ;   woutpretag[3]= 46515;   wouttag[3]= 5584.25;  
            fbb[3]= 0.0691777;   fcc[3]= 0.134739;   fc[3]= 0.114155;   fll[3]= 0.681929;  
            fmcbb[3]= 0.0503749;   fmccc[3]= 0.098116;   fmcc[3]= 0.163161;   fmcll[3]= 0.688349;  
           Fbb[4]= 1.34711;   Fcc[4]= 1.34711;   Fc[4]= 0.686324;   Fll[4]= 0.971811; 
              cafac[4]= 0.803713;  
              winpretag[4]= 13079.2 ;   wintag[4]= 1991.04 ;   woutpretag[4]= 10511.9;   wouttag[4]= 1748.61;  
            fbb[4]= 0.0962964;   fcc[4]= 0.159141;   fc[4]= 0.103314;   fll[4]= 0.641248;  
            fmcbb[4]= 0.0714836;   fmccc[4]= 0.118135;   fmcc[4]= 0.150533;   fmcll[4]= 0.659848;  
           Fbb[5]= 1.31692;   Fcc[5]= 1.31692;   Fc[5]= 0.670941;   Fll[5]= 0.950029; 
              cafac[5]= 0.793792;  
              winpretag[5]= 3455.39 ;   wintag[5]= 662.194 ;   woutpretag[5]= 2742.86;   wouttag[5]= 585.998;  
            fbb[5]= 0.123049;   fcc[5]= 0.191518;   fc[5]= 0.090551;   fll[5]= 0.594882;  
            fmcbb[5]= 0.0934369;   fmccc[5]= 0.145429;   fmcc[5]= 0.134961;   fmcll[5]= 0.626173;  
           Fbb[6]= 1.36602;   Fcc[6]= 1.36602;   Fc[6]= 0.69596;   Fll[6]= 0.985455; 
              cafac[6]= 0.787784;  
              winpretag[6]= 75851.5 ;   wintag[6]= 9386.92 ;   woutpretag[6]= 59754.6;   wouttag[6]= 7921.66;  
            fbb[6]= 0.0764651;   fcc[6]= 0.141688;   fc[6]= 0.111144;   fll[6]= 0.670703;  
            fmcbb[6]= 0.0559763;   fmccc[6]= 0.103723;   fmcc[6]= 0.159699;   fmcll[6]= 0.680602;  
           Fbb[7]= 1.34069;   Fcc[7]= 1.34069;   Fc[7]= 0.683052;   Fll[7]= 0.967177; 
              cafac[7]= 0.800916;  
              winpretag[7]= 16534.6 ;   wintag[7]= 2653.24 ;   woutpretag[7]= 13242.8;   wouttag[7]= 2336.15;  
            fbb[7]= 0.101988;   fcc[7]= 0.166029;   fc[7]= 0.100599;   fll[7]= 0.631384;  
            fmcbb[7]= 0.0760714;   fmccc[7]= 0.123839;   fmcc[7]= 0.147279;   fmcll[7]= 0.652811;  
      }  
      else if ( idsys==80   || sysname=="leff_down" ) { 
           Fbb[1]= 1.44908;   Fcc[1]= 1.44908;   Fc[1]= 0.757611;   Fll[1]= 1.01575; 
              cafac[1]= 1.00863;  
              winpretag[1]= 1.05581e+06 ;   wintag[1]= 39117.2 ;   woutpretag[1]= 1.06492e+06;   wouttag[1]= 36605.9;  
            fbb[1]= 0.0174697;   fcc[1]= 0.0508274;   fc[1]= 0.106157;   fll[1]= 0.825546;  
            fmcbb[1]= 0.0120557;   fmccc[1]= 0.0350757;   fmcc[1]= 0.14012;   fmcll[1]= 0.812748;  
           Fbb[2]= 1.42318;   Fcc[2]= 1.42318;   Fc[2]= 0.744074;   Fll[2]= 0.997597; 
              cafac[2]= 0.903597;  
              winpretag[2]= 250680 ;   wintag[2]= 19318.7 ;   woutpretag[2]= 226514;   wouttag[2]= 17829.5;  
            fbb[2]= 0.042833;   fcc[2]= 0.102765;   fc[2]= 0.120734;   fll[2]= 0.733668;  
            fmcbb[2]= 0.0300966;   fmccc[2]= 0.0722077;   fmcc[2]= 0.16226;   fmcll[2]= 0.735436;  
           Fbb[3]= 1.39611;   Fcc[3]= 1.39611;   Fc[3]= 0.729917;   Fll[3]= 0.978615; 
              cafac[3]= 0.832666;  
              winpretag[3]= 56406.9 ;   wintag[3]= 6404.51 ;   woutpretag[3]= 46968.1;   wouttag[3]= 5733.86;  
            fbb[3]= 0.0703259;   fcc[3]= 0.136968;   fc[3]= 0.119172;   fll[3]= 0.673535;  
            fmcbb[3]= 0.0503729;   fmccc[3]= 0.0981072;   fmcc[3]= 0.163267;   fmcll[3]= 0.688253;  
           Fbb[4]= 1.36831;   Fcc[4]= 1.36831;   Fc[4]= 0.715386;   Fll[4]= 0.959134; 
              cafac[4]= 0.853836;  
              winpretag[4]= 12435.2 ;   wintag[4]= 1893.24 ;   woutpretag[4]= 10617.6;   wouttag[4]= 1791.8;  
            fbb[4]= 0.0977946;   fcc[4]= 0.16163;   fc[4]= 0.107747;   fll[4]= 0.632828;  
            fmcbb[4]= 0.0714709;   fmccc[4]= 0.118123;   fmcc[4]= 0.150614;   fmcll[4]= 0.659791;  
           Fbb[5]= 1.33632;   Fcc[5]= 1.33632;   Fc[5]= 0.698658;   Fll[5]= 0.936706; 
              cafac[5]= 0.842364;  
              winpretag[5]= 3284.31 ;   wintag[5]= 629.328 ;   woutpretag[5]= 2766.59;   wouttag[5]= 598.576;  
            fbb[5]= 0.124794;   fcc[5]= 0.19437;   fc[5]= 0.0943549;   fll[5]= 0.586481;  
            fmcbb[5]= 0.0933867;   fmccc[5]= 0.145452;   fmcc[5]= 0.135052;   fmcll[5]= 0.62611;  
           Fbb[6]= 1.38841;   Fcc[6]= 1.38841;   Fc[6]= 0.725896;   Fll[6]= 0.973225; 
              cafac[6]= 0.83652;  
              winpretag[6]= 72126.4 ;   wintag[6]= 8927.08 ;   woutpretag[6]= 60335.2;   wouttag[6]= 8127.52;  
            fbb[6]= 0.0777083;   fcc[6]= 0.143998;   fc[6]= 0.115999;   fll[6]= 0.662295;  
            fmcbb[6]= 0.0559691;   fmccc[6]= 0.103714;   fmcc[6]= 0.159801;   fmcll[6]= 0.680516;  
           Fbb[7]= 1.3615;   Fcc[7]= 1.3615;   Fc[7]= 0.711825;   Fll[7]= 0.95436; 
              cafac[7]= 0.850618;  
              winpretag[7]= 15719.5 ;   wintag[7]= 2522.57 ;   woutpretag[7]= 13371.3;   wouttag[7]= 2391.99;  
            fbb[7]= 0.103542;   fcc[7]= 0.168599;   fc[7]= 0.104897;   fll[7]= 0.622962;  
            fmcbb[7]= 0.0760498;   fmccc[7]= 0.123833;   fmcc[7]= 0.147363;   fmcll[7]= 0.652754;  
      }  
      else if ( idsys==81   || sysname=="ca_up" ) { 
           Fbb[1]= 1.43394;   Fcc[1]= 1.43394;   Fc[1]= 0.740266;   Fll[1]= 1.0196; 
              cafac[1]= 0.985946;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06759e+06;   wouttag[1]= 36251;  
            fbb[1]= 0.0172893;   fcc[1]= 0.0503027;   fc[1]= 0.103702;   fll[1]= 0.828706;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.890703;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 229014;   wouttag[2]= 17854;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38474;   Fcc[3]= 1.38474;   Fc[3]= 0.714865;   Fll[3]= 0.984614; 
              cafac[3]= 0.831298;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 48100.6;   wouttag[3]= 5823.59;  
            fbb[3]= 0.0697545;   fcc[3]= 0.135859;   fc[3]= 0.116675;   fll[3]= 0.677712;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35777;   Fcc[4]= 1.35777;   Fc[4]= 0.700942;   Fll[4]= 0.965438; 
              cafac[4]= 0.885466;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 11296;   wouttag[4]= 1892.74;  
            fbb[4]= 0.0970496;   fcc[4]= 0.160392;   fc[4]= 0.105543;   fll[4]= 0.637016;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32667;   Fcc[5]= 1.32667;   Fc[5]= 0.68489;   Fll[5]= 0.943328; 
              cafac[5]= 0.945586;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 3186.49;   wouttag[5]= 685.128;  
            fbb[5]= 0.123928;   fcc[5]= 0.192951;   fc[5]= 0.0924638;   fll[5]= 0.590657;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37727;   Fcc[6]= 1.37727;   Fc[6]= 0.711013;   Fll[6]= 0.979308; 
              cafac[6]= 0.833424;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 61664.2;   wouttag[6]= 8241;  
            fbb[6]= 0.0770898;   fcc[6]= 0.142849;   fc[6]= 0.113583;   fll[6]= 0.666478;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.35115;   Fcc[7]= 1.35115;   Fc[7]= 0.697526;   Fll[7]= 0.960733; 
              cafac[7]= 0.878296;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 14164.3;   wouttag[7]= 2516.38;  
            fbb[7]= 0.10277;   fcc[7]= 0.167321;   fc[7]= 0.102759;   fll[7]= 0.62715;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==82   || sysname=="ca_down" ) { 
           Fbb[1]= 1.43394;   Fcc[1]= 1.43394;   Fc[1]= 0.740266;   Fll[1]= 1.0196; 
              cafac[1]= 0.974088;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.05475e+06;   wouttag[1]= 35815;  
            fbb[1]= 0.0172893;   fcc[1]= 0.0503027;   fc[1]= 0.103702;   fll[1]= 0.828706;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.863716;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 222075;   wouttag[2]= 17313;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38474;   Fcc[3]= 1.38474;   Fc[3]= 0.714865;   Fll[3]= 0.984614; 
              cafac[3]= 0.784331;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 45382.9;   wouttag[3]= 5494.56;  
            fbb[3]= 0.0697545;   fcc[3]= 0.135859;   fc[3]= 0.116675;   fll[3]= 0.677712;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35777;   Fcc[4]= 1.35777;   Fc[4]= 0.700942;   Fll[4]= 0.965438; 
              cafac[4]= 0.770832;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 9833.64;   wouttag[4]= 1647.7;  
            fbb[4]= 0.0970496;   fcc[4]= 0.160392;   fc[4]= 0.105543;   fll[4]= 0.637016;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32667;   Fcc[5]= 1.32667;   Fc[5]= 0.68489;   Fll[5]= 0.943328; 
              cafac[5]= 0.689344;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2322.99;   wouttag[5]= 499.467;  
            fbb[5]= 0.123928;   fcc[5]= 0.192951;   fc[5]= 0.0924638;   fll[5]= 0.590657;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37727;   Fcc[6]= 1.37727;   Fc[6]= 0.711013;   Fll[6]= 0.979308; 
              cafac[6]= 0.789662;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 58426.3;   wouttag[6]= 7808.27;  
            fbb[6]= 0.0770898;   fcc[6]= 0.142849;   fc[6]= 0.113583;   fll[6]= 0.666478;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.35115;   Fcc[7]= 1.35115;   Fc[7]= 0.697526;   Fll[7]= 0.960733; 
              cafac[7]= 0.771996;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 12450;   wouttag[7]= 2211.82;  
            fbb[7]= 0.10277;   fcc[7]= 0.167321;   fc[7]= 0.102759;   fll[7]= 0.62715;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==83   || sysname=="chmis_up" ) { 
           Fbb[1]= 1.43394;   Fcc[1]= 1.43394;   Fc[1]= 0.740266;   Fll[1]= 1.0196; 
              cafac[1]= 0.985897;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06754e+06;   wouttag[1]= 36249.2;  
            fbb[1]= 0.0172893;   fcc[1]= 0.0503027;   fc[1]= 0.103702;   fll[1]= 0.828706;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.881595;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 226672;   wouttag[2]= 17671.4;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38474;   Fcc[3]= 1.38474;   Fc[3]= 0.714865;   Fll[3]= 0.984614; 
              cafac[3]= 0.813469;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 47068.9;   wouttag[3]= 5698.69;  
            fbb[3]= 0.0697545;   fcc[3]= 0.135859;   fc[3]= 0.116675;   fll[3]= 0.677712;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35777;   Fcc[4]= 1.35777;   Fc[4]= 0.700942;   Fll[4]= 0.965438; 
              cafac[4]= 0.833946;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10638.8;   wouttag[4]= 1782.61;  
            fbb[4]= 0.0970496;   fcc[4]= 0.160392;   fc[4]= 0.105543;   fll[4]= 0.637016;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32667;   Fcc[5]= 1.32667;   Fc[5]= 0.68489;   Fll[5]= 0.943328; 
              cafac[5]= 0.824005;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2776.78;   wouttag[5]= 597.036;  
            fbb[5]= 0.123928;   fcc[5]= 0.192951;   fc[5]= 0.0924638;   fll[5]= 0.590657;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37727;   Fcc[6]= 1.37727;   Fc[6]= 0.711013;   Fll[6]= 0.979308; 
              cafac[6]= 0.818035;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60525.6;   wouttag[6]= 8088.83;  
            fbb[6]= 0.0770898;   fcc[6]= 0.142849;   fc[6]= 0.113583;   fll[6]= 0.666478;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.35115;   Fcc[7]= 1.35115;   Fc[7]= 0.697526;   Fll[7]= 0.960733; 
              cafac[7]= 0.831747;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13413.6;   wouttag[7]= 2383.01;  
            fbb[7]= 0.10277;   fcc[7]= 0.167321;   fc[7]= 0.102759;   fll[7]= 0.62715;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==84   || sysname=="chmis_down" ) { 
           Fbb[1]= 1.43394;   Fcc[1]= 1.43394;   Fc[1]= 0.740266;   Fll[1]= 1.0196; 
              cafac[1]= 0.974137;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.05481e+06;   wouttag[1]= 35816.8;  
            fbb[1]= 0.0172893;   fcc[1]= 0.0503027;   fc[1]= 0.103702;   fll[1]= 0.828706;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41041;   Fcc[2]= 1.41041;   Fc[2]= 0.728121;   Fll[2]= 1.00287; 
              cafac[2]= 0.872823;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 224416;   wouttag[2]= 17495.6;  
            fbb[2]= 0.0424539;   fcc[2]= 0.101855;   fc[2]= 0.118116;   fll[2]= 0.737575;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38474;   Fcc[3]= 1.38474;   Fc[3]= 0.714865;   Fll[3]= 0.984614; 
              cafac[3]= 0.80216;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46414.5;   wouttag[3]= 5619.46;  
            fbb[3]= 0.0697545;   fcc[3]= 0.135859;   fc[3]= 0.116675;   fll[3]= 0.677712;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35777;   Fcc[4]= 1.35777;   Fc[4]= 0.700942;   Fll[4]= 0.965438; 
              cafac[4]= 0.822352;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10490.9;   wouttag[4]= 1757.83;  
            fbb[4]= 0.0970496;   fcc[4]= 0.160392;   fc[4]= 0.105543;   fll[4]= 0.637016;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32667;   Fcc[5]= 1.32667;   Fc[5]= 0.68489;   Fll[5]= 0.943328; 
              cafac[5]= 0.810925;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2732.7;   wouttag[5]= 587.559;  
            fbb[5]= 0.123928;   fcc[5]= 0.192951;   fc[5]= 0.0924638;   fll[5]= 0.590657;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37727;   Fcc[6]= 1.37727;   Fc[6]= 0.711013;   Fll[6]= 0.979308; 
              cafac[6]= 0.805051;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 59564.9;   wouttag[6]= 7960.44;  
            fbb[6]= 0.0770898;   fcc[6]= 0.142849;   fc[6]= 0.113583;   fll[6]= 0.666478;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.35115;   Fcc[7]= 1.35115;   Fc[7]= 0.697526;   Fll[7]= 0.960733; 
              cafac[7]= 0.818545;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13200.7;   wouttag[7]= 2345.19;  
            fbb[7]= 0.10277;   fcc[7]= 0.167321;   fc[7]= 0.102759;   fll[7]= 0.62715;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==67   || sysname=="eer_up" ) { 
           Fbb[1]= 1.3817;   Fcc[1]= 1.3817;   Fc[1]= 0.751586;   Fll[1]= 1.02076; 
              cafac[1]= 0.978959;  
              winpretag[1]= 1.08199e+06 ;   wintag[1]= 40132.6 ;   woutpretag[1]= 1.05922e+06;   wouttag[1]= 35958.1;  
            fbb[1]= 0.0165963;   fcc[1]= 0.0484196;   fc[1]= 0.105385;   fll[1]= 0.829599;  
            fmcbb[1]= 0.0120115;   fmccc[1]= 0.0350434;   fmcc[1]= 0.140216;   fmcll[1]= 0.812729;  
           Fbb[2]= 1.36253;   Fcc[2]= 1.36253;   Fc[2]= 0.741157;   Fll[2]= 1.00659; 
              cafac[2]= 0.888716;  
              winpretag[2]= 256893 ;   wintag[2]= 19903 ;   woutpretag[2]= 228305;   wouttag[2]= 17748;  
            fbb[2]= 0.0410714;   fcc[2]= 0.0986623;   fc[2]= 0.120334;   fll[2]= 0.739932;  
            fmcbb[2]= 0.0301434;   fmccc[2]= 0.072411;   fmcc[2]= 0.16236;   fmcll[2]= 0.735086;  
           Fbb[3]= 1.34088;   Fcc[3]= 1.34088;   Fc[3]= 0.729377;   Fll[3]= 0.990595; 
              cafac[3]= 0.807736;  
              winpretag[3]= 57991.9 ;   wintag[3]= 6600.27 ;   woutpretag[3]= 46842.1;   wouttag[3]= 5615.98;  
            fbb[3]= 0.0672671;   fcc[3]= 0.131644;   fc[3]= 0.118828;   fll[3]= 0.682261;  
            fmcbb[3]= 0.0501666;   fmccc[3]= 0.0981779;   fmcc[3]= 0.162917;   fmcll[3]= 0.688738;  
           Fbb[4]= 1.31773;   Fcc[4]= 1.31773;   Fc[4]= 0.716785;   Fll[4]= 0.973494; 
              cafac[4]= 0.828741;  
              winpretag[4]= 12735.8 ;   wintag[4]= 1940.55 ;   woutpretag[4]= 10554.7;   wouttag[4]= 1745.74;  
            fbb[4]= 0.0941572;   fcc[4]= 0.155574;   fc[4]= 0.108147;   fll[4]= 0.642122;  
            fmcbb[4]= 0.0714543;   fmccc[4]= 0.118062;   fmcc[4]= 0.150878;   fmcll[4]= 0.659605;  
           Fbb[5]= 1.29045;   Fcc[5]= 1.29045;   Fc[5]= 0.701949;   Fll[5]= 0.953344; 
              cafac[5]= 0.843325;  
              winpretag[5]= 3354.31 ;   wintag[5]= 644.635 ;   woutpretag[5]= 2828.77;   wouttag[5]= 600.664;  
            fbb[5]= 0.120615;   fcc[5]= 0.188737;   fc[5]= 0.0953719;   fll[5]= 0.595276;  
            fmcbb[5]= 0.0934671;   fmccc[5]= 0.146256;   fmcc[5]= 0.135867;   fmcll[5]= 0.624409;  
           Fbb[6]= 1.33448;   Fcc[6]= 1.33448;   Fc[6]= 0.7259;   Fll[6]= 0.985873; 
              cafac[6]= 0.813058;  
              winpretag[6]= 74082 ;   wintag[6]= 9185.45 ;   woutpretag[6]= 60233;   wouttag[6]= 7959.39;  
            fbb[6]= 0.0744466;   fcc[6]= 0.138484;   fc[6]= 0.11587;   fll[6]= 0.671199;  
            fmcbb[6]= 0.0557868;   fmccc[6]= 0.103773;   fmcc[6]= 0.159623;   fmcll[6]= 0.680817;  
           Fbb[7]= 1.31195;   Fcc[7]= 1.31195;   Fc[7]= 0.713641;   Fll[7]= 0.969223; 
              cafac[7]= 0.831712;  
              winpretag[7]= 16090.1 ;   wintag[7]= 2585.18 ;   woutpretag[7]= 13382.4;   wouttag[7]= 2346.57;  
            fbb[7]= 0.0997647;   fcc[7]= 0.162603;   fc[7]= 0.105439;   fll[7]= 0.632193;  
            fmcbb[7]= 0.0760433;   fmccc[7]= 0.12394;   fmcc[7]= 0.147749;   fmcll[7]= 0.652268;  
      }  
      else if ( idsys==68   || sysname=="eer_down" ) { 
           Fbb[1]= 1.37384;   Fcc[1]= 1.37384;   Fc[1]= 0.743392;   Fll[1]= 1.0226; 
              cafac[1]= 0.977546;  
              winpretag[1]= 1.08204e+06 ;   wintag[1]= 40055.3 ;   woutpretag[1]= 1.05774e+06;   wouttag[1]= 35629.6;  
            fbb[1]= 0.0164774;   fcc[1]= 0.0481385;   fc[1]= 0.104159;   fll[1]= 0.831225;  
            fmcbb[1]= 0.0119937;   fmccc[1]= 0.0350394;   fmcc[1]= 0.140114;   fmcll[1]= 0.812853;  
           Fbb[2]= 1.35572;   Fcc[2]= 1.35572;   Fc[2]= 0.733591;   Fll[2]= 1.00912; 
              cafac[2]= 0.889023;  
              winpretag[2]= 256846 ;   wintag[2]= 19944.7 ;   woutpretag[2]= 228342;   wouttag[2]= 17714.1;  
            fbb[2]= 0.0409202;   fcc[2]= 0.0982389;   fc[2]= 0.119005;   fll[2]= 0.741836;  
            fmcbb[2]= 0.0301833;   fmccc[2]= 0.0724624;   fmcc[2]= 0.162222;   fmcll[2]= 0.735132;  
           Fbb[3]= 1.33466;   Fcc[3]= 1.33466;   Fc[3]= 0.722191;   Fll[3]= 0.993438; 
              cafac[3]= 0.80485;  
              winpretag[3]= 57979.8 ;   wintag[3]= 6602.87 ;   woutpretag[3]= 46665.1;   wouttag[3]= 5573.59;  
            fbb[3]= 0.0670745;   fcc[3]= 0.131266;   fc[3]= 0.117538;   fll[3]= 0.684122;  
            fmcbb[3]= 0.050256;   fmccc[3]= 0.0983516;   fmcc[3]= 0.162752;   fmcll[3]= 0.68864;  
           Fbb[4]= 1.31204;   Fcc[4]= 1.31204;   Fc[4]= 0.709953;   Fll[4]= 0.976604; 
              cafac[4]= 0.836009;  
              winpretag[4]= 12708.5 ;   wintag[4]= 1938.38 ;   woutpretag[4]= 10624.4;   wouttag[4]= 1752.73;  
            fbb[4]= 0.0939118;   fcc[4]= 0.15497;   fc[4]= 0.107119;   fll[4]= 0.643999;  
            fmcbb[4]= 0.071577;   fmccc[4]= 0.118114;   fmcc[4]= 0.150882;   fmcll[4]= 0.659427;  
           Fbb[5]= 1.2853;   Fcc[5]= 1.2853;   Fc[5]= 0.695482;   Fll[5]= 0.956698; 
              cafac[5]= 0.84586;  
              winpretag[5]= 3355.05 ;   wintag[5]= 647.361 ;   woutpretag[5]= 2837.91;   wouttag[5]= 602.392;  
            fbb[5]= 0.11901;   fcc[5]= 0.188727;   fc[5]= 0.0941824;   fll[5]= 0.59808;  
            fmcbb[5]= 0.0925935;   fmccc[5]= 0.146835;   fmcc[5]= 0.13542;   fmcll[5]= 0.625151;  
           Fbb[6]= 1.32841;   Fcc[6]= 1.32841;   Fc[6]= 0.718814;   Fll[6]= 0.988792; 
              cafac[6]= 0.812288;  
              winpretag[6]= 74043.4 ;   wintag[6]= 9188.6 ;   woutpretag[6]= 60144.5;   wouttag[6]= 7920.92;  
            fbb[6]= 0.0741705;   fcc[6]= 0.138076;   fc[6]= 0.114634;   fll[6]= 0.67312;  
            fmcbb[6]= 0.0558339;   fmccc[6]= 0.10394;   fmcc[6]= 0.159476;   fmcll[6]= 0.68075;  
           Fbb[7]= 1.30636;   Fcc[7]= 1.30636;   Fc[7]= 0.706881;   Fll[7]= 0.972378; 
              cafac[7]= 0.83789;  
              winpretag[7]= 16063.5 ;   wintag[7]= 2585.74 ;   woutpretag[7]= 13459.5;   wouttag[7]= 2355.57;  
            fbb[7]= 0.0992398;   fcc[7]= 0.162136;   fc[7]= 0.104373;   fll[7]= 0.634251;  
            fmcbb[7]= 0.0759665;   fmccc[7]= 0.124113;   fmcc[7]= 0.147653;   fmcll[7]= 0.652268;  
      }  
      else if ( idsys==65   || sysname=="ees_up" ) { 
           Fbb[1]= 1.41554;   Fcc[1]= 1.41554;   Fc[1]= 0.735009;   Fll[1]= 1.02163; 
              cafac[1]= 0.971324;  
              winpretag[1]= 1.08741e+06 ;   wintag[1]= 40278.9 ;   woutpretag[1]= 1.05623e+06;   wouttag[1]= 35645.8;  
            fbb[1]= 0.0169811;   fcc[1]= 0.0496253;   fc[1]= 0.103009;   fll[1]= 0.830385;  
            fmcbb[1]= 0.0119962;   fmccc[1]= 0.0350575;   fmcc[1]= 0.140147;   fmcll[1]= 0.8128;  
           Fbb[2]= 1.39393;   Fcc[2]= 1.39393;   Fc[2]= 0.723789;   Fll[2]= 1.00604; 
              cafac[2]= 0.881464;  
              winpretag[2]= 257892 ;   wintag[2]= 19975 ;   woutpretag[2]= 227323;   wouttag[2]= 17706.9;  
            fbb[2]= 0.042044;   fcc[2]= 0.10094;   fc[2]= 0.117519;   fll[2]= 0.739497;  
            fmcbb[2]= 0.0301621;   fmccc[2]= 0.0724136;   fmcc[2]= 0.162366;   fmcll[2]= 0.735058;  
           Fbb[3]= 1.36968;   Fcc[3]= 1.36968;   Fc[3]= 0.711194;   Fll[3]= 0.988534; 
              cafac[3]= 0.801559;  
              winpretag[3]= 58185.6 ;   wintag[3]= 6632.66 ;   woutpretag[3]= 46639.2;   wouttag[3]= 5628.57;  
            fbb[3]= 0.0687336;   fcc[3]= 0.134639;   fc[3]= 0.115721;   fll[3]= 0.680906;  
            fmcbb[3]= 0.0501823;   fmccc[3]= 0.0983;   fmcc[3]= 0.162713;   fmcll[3]= 0.688804;  
           Fbb[4]= 1.34389;   Fcc[4]= 1.34389;   Fc[4]= 0.697805;   Fll[4]= 0.969924; 
              cafac[4]= 0.831613;  
              winpretag[4]= 12774.6 ;   wintag[4]= 1948.34 ;   woutpretag[4]= 10623.5;   wouttag[4]= 1773.32;  
            fbb[4]= 0.0963376;   fcc[4]= 0.158852;   fc[4]= 0.104974;   fll[4]= 0.639837;  
            fmcbb[4]= 0.0716855;   fmccc[4]= 0.118203;   fmcc[4]= 0.150434;   fmcll[4]= 0.659678;  
           Fbb[5]= 1.31452;   Fcc[5]= 1.31452;   Fc[5]= 0.682554;   Fll[5]= 0.948724; 
              cafac[5]= 0.835406;  
              winpretag[5]= 3358.16 ;   wintag[5]= 645.644 ;   woutpretag[5]= 2805.42;   wouttag[5]= 600.371;  
            fbb[5]= 0.122829;   fcc[5]= 0.191361;   fc[5]= 0.0927136;   fll[5]= 0.593096;  
            fmcbb[5]= 0.0934402;   fmccc[5]= 0.145575;   fmcc[5]= 0.135833;   fmcll[5]= 0.625151;  
           Fbb[6]= 1.3626;   Fcc[6]= 1.3626;   Fc[6]= 0.707519;   Fll[6]= 0.983426; 
              cafac[6]= 0.808382;  
              winpretag[6]= 74318.4 ;   wintag[6]= 9226.64 ;   woutpretag[6]= 60077.7;   wouttag[6]= 7997.13;  
            fbb[6]= 0.0760783;   fcc[6]= 0.141516;   fc[6]= 0.11277;   fll[6]= 0.669636;  
            fmcbb[6]= 0.0558332;   fmccc[6]= 0.103857;   fmcc[6]= 0.159388;   fmcll[6]= 0.680922;  
           Fbb[7]= 1.33767;   Fcc[7]= 1.33767;   Fc[7]= 0.694575;   Fll[7]= 0.965433; 
              cafac[7]= 0.832052;  
              winpretag[7]= 16132.8 ;   wintag[7]= 2593.98 ;   woutpretag[7]= 13423.3;   wouttag[7]= 2374.54;  
            fbb[7]= 0.101949;   fcc[7]= 0.165738;   fc[7]= 0.102377;   fll[7]= 0.629936;  
            fmcbb[7]= 0.0762139;   fmccc[7]= 0.1239;   fmcc[7]= 0.147395;   fmcll[7]= 0.652491;  
      }  
      else if ( idsys==66   || sysname=="ees_down" ) { 
           Fbb[1]= 1.39672;   Fcc[1]= 1.39672;   Fc[1]= 0.751293;   Fll[1]= 1.01988; 
              cafac[1]= 0.98804;  
              winpretag[1]= 1.0713e+06 ;   wintag[1]= 39655.5 ;   woutpretag[1]= 1.05849e+06;   wouttag[1]= 35934.3;  
            fbb[1]= 0.0167992;   fcc[1]= 0.048976;   fc[1]= 0.105246;   fll[1]= 0.828979;  
            fmcbb[1]= 0.0120275;   fmccc[1]= 0.0350649;   fmcc[1]= 0.140086;   fmcll[1]= 0.812822;  
           Fbb[2]= 1.37616;   Fcc[2]= 1.37616;   Fc[2]= 0.740234;   Fll[2]= 1.00487; 
              cafac[2]= 0.894938;  
              winpretag[2]= 254763 ;   wintag[2]= 19735.8 ;   woutpretag[2]= 227997;   wouttag[2]= 17778.8;  
            fbb[2]= 0.0414418;   fcc[2]= 0.0995927;   fc[2]= 0.12005;   fll[2]= 0.738916;  
            fmcbb[2]= 0.030114;   fmccc[2]= 0.0723697;   fmcc[2]= 0.162179;   fmcll[2]= 0.735338;  
           Fbb[3]= 1.35319;   Fcc[3]= 1.35319;   Fc[3]= 0.727877;   Fll[3]= 0.988091; 
              cafac[3]= 0.81222;  
              winpretag[3]= 57483.9 ;   wintag[3]= 6559.84 ;   woutpretag[3]= 46689.6;   wouttag[3]= 5635.9;  
            fbb[3]= 0.0681296;   fcc[3]= 0.13311;   fc[3]= 0.118567;   fll[3]= 0.680193;  
            fmcbb[3]= 0.0503473;   fmccc[3]= 0.0983677;   fmcc[3]= 0.162894;   fmcll[3]= 0.688391;  
           Fbb[4]= 1.32959;   Fcc[4]= 1.32959;   Fc[4]= 0.71518;   Fll[4]= 0.970855; 
              cafac[4]= 0.8374;  
              winpretag[4]= 12636.9 ;   wintag[4]= 1915.49 ;   woutpretag[4]= 10582.1;   wouttag[4]= 1749.83;  
            fbb[4]= 0.0945694;   fcc[4]= 0.15638;   fc[4]= 0.10787;   fll[4]= 0.641181;  
            fmcbb[4]= 0.0711269;   fmccc[4]= 0.117616;   fmcc[4]= 0.150829;   fmcll[4]= 0.660429;  
           Fbb[5]= 1.30141;   Fcc[5]= 1.30141;   Fc[5]= 0.700023;   Fll[5]= 0.950279; 
              cafac[5]= 0.844388;  
              winpretag[5]= 3327.22 ;   wintag[5]= 640.212 ;   woutpretag[5]= 2809.46;   wouttag[5]= 599.892;  
            fbb[5]= 0.120984;   fcc[5]= 0.189196;   fc[5]= 0.095015;   fll[5]= 0.594806;  
            fmcbb[5]= 0.0929636;   fmccc[5]= 0.145378;   fmcc[5]= 0.135731;   fmcll[5]= 0.625928;  
           Fbb[6]= 1.34665;   Fcc[6]= 1.34665;   Fc[6]= 0.724359;   Fll[6]= 0.983315; 
              cafac[6]= 0.8181;  
              winpretag[6]= 73448 ;   wintag[6]= 9115.54 ;   woutpretag[6]= 60087.8;   wouttag[6]= 7982.12;  
            fbb[6]= 0.0752146;   fcc[6]= 0.139794;   fc[6]= 0.115599;   fll[6]= 0.669393;  
            fmcbb[6]= 0.055853;   fmccc[6]= 0.103809;   fmcc[6]= 0.159587;   fmcll[6]= 0.680751;  
           Fbb[7]= 1.32361;   Fcc[7]= 1.32361;   Fc[7]= 0.711967;   Fll[7]= 0.966494; 
              cafac[7]= 0.838592;  
              winpretag[7]= 15964.1 ;   wintag[7]= 2555.7 ;   woutpretag[7]= 13387.4;   wouttag[7]= 2350.38;  
            fbb[7]= 0.100169;   fcc[7]= 0.163336;   fc[7]= 0.105145;   fll[7]= 0.63135;  
            fmcbb[7]= 0.075678;   fmccc[7]= 0.123402;   fmcc[7]= 0.147682;   fmcll[7]= 0.653238;  
      }  
      else if ( idsys==63   || sysname=="jer_one" ) { 
           Fbb[1]= 1.41568;   Fcc[1]= 1.41568;   Fc[1]= 0.662682;   Fll[1]= 1.03342; 
              cafac[1]= 0.897295;  
              winpretag[1]= 1.1556e+06 ;   wintag[1]= 42132.6 ;   woutpretag[1]= 1.03692e+06;   wouttag[1]= 33082;  
            fbb[1]= 0.0162653;   fcc[1]= 0.0471894;   fc[1]= 0.0903574;   fll[1]= 0.846188;  
            fmcbb[1]= 0.0114894;   fmccc[1]= 0.0333333;   fmcc[1]= 0.136351;   fmcll[1]= 0.818826;  
           Fbb[2]= 1.40054;   Fcc[2]= 1.40054;   Fc[2]= 0.655594;   Fll[2]= 1.02236; 
              cafac[2]= 0.849241;  
              winpretag[2]= 271754 ;   wintag[2]= 20779.4 ;   woutpretag[2]= 230785;   wouttag[2]= 17239.7;  
            fbb[2]= 0.0404039;   fcc[2]= 0.0984776;   fc[2]= 0.107005;   fll[2]= 0.754114;  
            fmcbb[2]= 0.0288488;   fmccc[2]= 0.0703141;   fmcc[2]= 0.163218;   fmcll[2]= 0.737619;  
           Fbb[3]= 1.37763;   Fcc[3]= 1.37763;   Fc[3]= 0.64487;   Fll[3]= 1.00564; 
              cafac[3]= 0.768509;  
              winpretag[3]= 62149.6 ;   wintag[3]= 7067.68 ;   woutpretag[3]= 47762.5;   wouttag[3]= 5640.14;  
            fbb[3]= 0.067248;   fcc[3]= 0.132759;   fc[3]= 0.106613;   fll[3]= 0.69338;  
            fmcbb[3]= 0.0488143;   fmccc[3]= 0.0963677;   fmcc[3]= 0.165326;   fmcll[3]= 0.689492;  
           Fbb[4]= 1.35132;   Fcc[4]= 1.35132;   Fc[4]= 0.632555;   Fll[4]= 0.986434; 
              cafac[4]= 0.773679;  
              winpretag[4]= 13947.1 ;   wintag[4]= 2131.61 ;   woutpretag[4]= 10790.6;   wouttag[4]= 1775.26;  
            fbb[4]= 0.0969002;   fcc[4]= 0.158461;   fc[4]= 0.0990041;   fll[4]= 0.645634;  
            fmcbb[4]= 0.0717077;   fmccc[4]= 0.117264;   fmcc[4]= 0.156515;   fmcll[4]= 0.654514;  
           Fbb[5]= 1.32148;   Fcc[5]= 1.32148;   Fc[5]= 0.618585;   Fll[5]= 0.964649; 
              cafac[5]= 0.746767;  
              winpretag[5]= 3789.55 ;   wintag[5]= 725.951 ;   woutpretag[5]= 2829.91;   wouttag[5]= 591.816;  
            fbb[5]= 0.123263;   fcc[5]= 0.187874;   fc[5]= 0.0869842;   fll[5]= 0.601879;  
            fmcbb[5]= 0.0932764;   fmccc[5]= 0.14217;   fmcc[5]= 0.140618;   fmcll[5]= 0.623936;  
           Fbb[6]= 1.37021;   Fcc[6]= 1.37021;   Fc[6]= 0.641397;   Fll[6]= 1.00022; 
              cafac[6]= 0.767786;  
              winpretag[6]= 79886.3 ;   wintag[6]= 9925.24 ;   woutpretag[6]= 61335.6;   wouttag[6]= 8021.77;  
            fbb[6]= 0.0752524;   fcc[6]= 0.14002;   fc[6]= 0.104301;   fll[6]= 0.680427;  
            fmcbb[6]= 0.0549204;   fmccc[6]= 0.102189;   fmcc[6]= 0.162615;   fmcll[6]= 0.680276;  
           Fbb[7]= 1.34483;   Fcc[7]= 1.34483;   Fc[7]= 0.629518;   Fll[7]= 0.981697; 
              cafac[7]= 0.76699;  
              winpretag[7]= 17736.7 ;   wintag[7]= 2857.56 ;   woutpretag[7]= 13603.8;   wouttag[7]= 2370.05;  
            fbb[7]= 0.102632;   fcc[7]= 0.164857;   fc[7]= 0.0963906;   fll[7]= 0.636121;  
            fmcbb[7]= 0.076316;   fmccc[7]= 0.122585;   fmcc[7]= 0.153118;   fmcll[7]= 0.64798;  
      }  
      else if ( idsys==64   || sysname=="jef_one" ) { 
           Fbb[1]= 1.39584;   Fcc[1]= 1.39584;   Fc[1]= 0.731257;   Fll[1]= 1.02343; 
              cafac[1]= 0.975231;  
              winpretag[1]= 1.0815e+06 ;   wintag[1]= 40097.6 ;   woutpretag[1]= 1.05472e+06;   wouttag[1]= 35440.5;  
            fbb[1]= 0.0167825;   fcc[1]= 0.048938;   fc[1]= 0.102522;   fll[1]= 0.831758;  
            fmcbb[1]= 0.0120232;   fmccc[1]= 0.0350598;   fmcc[1]= 0.140199;   fmcll[1]= 0.812718;  
           Fbb[2]= 1.3763;   Fcc[2]= 1.3763;   Fc[2]= 0.721018;   Fll[2]= 1.0091; 
              cafac[2]= 0.886742;  
              winpretag[2]= 256688 ;   wintag[2]= 19898.4 ;   woutpretag[2]= 227616;   wouttag[2]= 17632.3;  
            fbb[2]= 0.0414552;   fcc[2]= 0.0997551;   fc[2]= 0.117066;   fll[2]= 0.741724;  
            fmcbb[2]= 0.0301208;   fmccc[2]= 0.0724807;   fmcc[2]= 0.162362;   fmcll[2]= 0.735037;  
           Fbb[3]= 1.35371;   Fcc[3]= 1.35371;   Fc[3]= 0.709184;   Fll[3]= 0.992536; 
              cafac[3]= 0.803039;  
              winpretag[3]= 57945.2 ;   wintag[3]= 6604.06 ;   woutpretag[3]= 46532.2;   wouttag[3]= 5582.58;  
            fbb[3]= 0.0680261;   fcc[3]= 0.132973;   fc[3]= 0.115539;   fll[3]= 0.683463;  
            fmcbb[3]= 0.0502516;   fmccc[3]= 0.0982282;   fmcc[3]= 0.162918;   fmcll[3]= 0.688603;  
           Fbb[4]= 1.32935;   Fcc[4]= 1.32935;   Fc[4]= 0.696421;   Fll[4]= 0.974673; 
              cafac[4]= 0.82301;  
              winpretag[4]= 12721 ;   wintag[4]= 1931.43 ;   woutpretag[4]= 10469.5;   wouttag[4]= 1728.82;  
            fbb[4]= 0.0952734;   fcc[4]= 0.1568;   fc[4]= 0.104938;   fll[4]= 0.642989;  
            fmcbb[4]= 0.0716692;   fmccc[4]= 0.117953;   fmcc[4]= 0.150682;   fmcll[4]= 0.659696;  
           Fbb[5]= 1.3004;   Fcc[5]= 1.3004;   Fc[5]= 0.681253;   Fll[5]= 0.953445; 
              cafac[5]= 0.842125;  
              winpretag[5]= 3344.24 ;   wintag[5]= 643.042 ;   woutpretag[5]= 2816.27;   wouttag[5]= 599.948;  
            fbb[5]= 0.121963;   fcc[5]= 0.190405;   fc[5]= 0.0920693;   fll[5]= 0.595563;  
            fmcbb[5]= 0.0937891;   fmccc[5]= 0.146421;   fmcc[5]= 0.135147;   fmcll[5]= 0.624643;  
           Fbb[6]= 1.34697;   Fcc[6]= 1.34697;   Fc[6]= 0.705654;   Fll[6]= 0.987595; 
              cafac[6]= 0.808344;  
              winpretag[6]= 74010.4 ;   wintag[6]= 9178.53 ;   woutpretag[6]= 59825.8;   wouttag[6]= 7908.52;  
            fbb[6]= 0.075296;   fcc[6]= 0.13981;   fc[6]= 0.112594;   fll[6]= 0.6723;  
            fmcbb[6]= 0.0559002;   fmccc[6]= 0.103796;   fmcc[6]= 0.15956;   fmcll[6]= 0.680744;  
           Fbb[7]= 1.32322;   Fcc[7]= 1.32322;   Fc[7]= 0.693208;   Fll[7]= 0.970177; 
              cafac[7]= 0.827007;  
              winpretag[7]= 16065.2 ;   wintag[7]= 2574.47 ;   woutpretag[7]= 13286;   wouttag[7]= 2328.72;  
            fbb[7]= 0.100927;   fcc[7]= 0.163918;   fc[7]= 0.102212;   fll[7]= 0.632943;  
            fmcbb[7]= 0.0762738;   fmccc[7]= 0.123879;   fmcc[7]= 0.147448;   fmcll[7]= 0.6524;  
      }  
      else if ( idsys==61   || sysname=="musms_up" ) { 
           Fbb[1]= 1.41115;   Fcc[1]= 1.41115;   Fc[1]= 0.727225;   Fll[1]= 1.02325; 
              cafac[1]= 0.974177;  
              winpretag[1]= 1.08193e+06 ;   wintag[1]= 40096.9 ;   woutpretag[1]= 1.05399e+06;   wouttag[1]= 35388.6;  
            fbb[1]= 0.0169434;   fcc[1]= 0.0494488;   fc[1]= 0.101947;   fll[1]= 0.831661;  
            fmcbb[1]= 0.0120068;   fmccc[1]= 0.0350416;   fmcc[1]= 0.140186;   fmcll[1]= 0.812765;  
           Fbb[2]= 1.39033;   Fcc[2]= 1.39033;   Fc[2]= 0.716496;   Fll[2]= 1.00815; 
              cafac[2]= 0.884354;  
              winpretag[2]= 256891 ;   wintag[2]= 19908.9 ;   woutpretag[2]= 227183;   wouttag[2]= 17635.7;  
            fbb[2]= 0.0418891;   fcc[2]= 0.100704;   fc[2]= 0.116318;   fll[2]= 0.741089;  
            fmcbb[2]= 0.030129;   fmccc[2]= 0.0724322;   fmcc[2]= 0.162342;   fmcll[2]= 0.735097;  
           Fbb[3]= 1.36659;   Fcc[3]= 1.36659;   Fc[3]= 0.704265;   Fll[3]= 0.990941; 
              cafac[3]= 0.800842;  
              winpretag[3]= 57992.6 ;   wintag[3]= 6605.33 ;   woutpretag[3]= 46442.9;   wouttag[3]= 5586.76;  
            fbb[3]= 0.0686555;   fcc[3]= 0.134122;   fc[3]= 0.114678;   fll[3]= 0.682544;  
            fmcbb[3]= 0.0502385;   fmccc[3]= 0.0981433;   fmcc[3]= 0.162834;   fmcll[3]= 0.688784;  
           Fbb[4]= 1.34116;   Fcc[4]= 1.34116;   Fc[4]= 0.69116;   Fll[4]= 0.972502; 
              cafac[4]= 0.823928;  
              winpretag[4]= 12735.3 ;   wintag[4]= 1936.54 ;   woutpretag[4]= 10493;   wouttag[4]= 1741.33;  
            fbb[4]= 0.0960273;   fcc[4]= 0.158334;   fc[4]= 0.104214;   fll[4]= 0.641425;  
            fmcbb[4]= 0.0716001;   fmccc[4]= 0.118057;   fmcc[4]= 0.150781;   fmcll[4]= 0.659562;  
           Fbb[5]= 1.31132;   Fcc[5]= 1.31132;   Fc[5]= 0.675783;   Fll[5]= 0.950865; 
              cafac[5]= 0.837226;  
              winpretag[5]= 3353.98 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2808.04;   wouttag[5]= 601.495;  
            fbb[5]= 0.122879;   fcc[5]= 0.191622;   fc[5]= 0.0916713;   fll[5]= 0.593828;  
            fmcbb[5]= 0.0937061;   fmccc[5]= 0.146129;   fmcc[5]= 0.135652;   fmcll[5]= 0.624513;  
           Fbb[6]= 1.35957;   Fcc[6]= 1.35957;   Fc[6]= 0.700644;   Fll[6]= 0.985847; 
              cafac[6]= 0.806553;  
              winpretag[6]= 74081.9 ;   wintag[6]= 9187.63 ;   woutpretag[6]= 59751;   wouttag[6]= 7926.46;  
            fbb[6]= 0.0759707;   fcc[6]= 0.14104;   fc[6]= 0.111775;   fll[6]= 0.671214;  
            fmcbb[6]= 0.0558787;   fmccc[6]= 0.103739;   fmcc[6]= 0.159531;   fmcll[6]= 0.680851;  
           Fbb[7]= 1.33483;   Fcc[7]= 1.33483;   Fc[7]= 0.687897;   Fll[7]= 0.967911; 
              cafac[7]= 0.826574;  
              winpretag[7]= 16089.3 ;   wintag[7]= 2582.3 ;   woutpretag[7]= 13299;   wouttag[7]= 2343.17;  
            fbb[7]= 0.101725;   fcc[7]= 0.165397;   fc[7]= 0.101552;   fll[7]= 0.631325;  
            fmcbb[7]= 0.0762083;   fmccc[7]= 0.123909;   fmcc[7]= 0.147627;   fmcll[7]= 0.652256;  
      }  
      else if ( idsys==62   || sysname=="musms_down" ) { 
           Fbb[1]= 1.4099;   Fcc[1]= 1.4099;   Fc[1]= 0.728011;   Fll[1]= 1.02319; 
              cafac[1]= 0.974329;  
              winpretag[1]= 1.08193e+06 ;   wintag[1]= 40097.6 ;   woutpretag[1]= 1.05416e+06;   wouttag[1]= 35404.9;  
            fbb[1]= 0.0169283;   fcc[1]= 0.049405;   fc[1]= 0.102058;   fll[1]= 0.831609;  
            fmcbb[1]= 0.0120068;   fmccc[1]= 0.0350416;   fmcc[1]= 0.140187;   fmcll[1]= 0.812764;  
           Fbb[2]= 1.38916;   Fcc[2]= 1.38916;   Fc[2]= 0.717302;   Fll[2]= 1.00813; 
              cafac[2]= 0.884477;  
              winpretag[2]= 256892 ;   wintag[2]= 19908.1 ;   woutpretag[2]= 227215;   wouttag[2]= 17637.2;  
            fbb[2]= 0.0418572;   fcc[2]= 0.100622;   fc[2]= 0.116448;   fll[2]= 0.741073;  
            fmcbb[2]= 0.0301314;   fmccc[2]= 0.0724342;   fmcc[2]= 0.162341;   fmcll[2]= 0.735093;  
           Fbb[3]= 1.36551;   Fcc[3]= 1.36551;   Fc[3]= 0.705094;   Fll[3]= 0.990977; 
              cafac[3]= 0.800997;  
              winpretag[3]= 57992.6 ;   wintag[3]= 6605.33 ;   woutpretag[3]= 46451.9;   wouttag[3]= 5586.91;  
            fbb[3]= 0.0686014;   fcc[3]= 0.134016;   fc[3]= 0.114813;   fll[3]= 0.682569;  
            fmcbb[3]= 0.0502385;   fmccc[3]= 0.0981433;   fmcc[3]= 0.162834;   fmcll[3]= 0.688784;  
           Fbb[4]= 1.34018;   Fcc[4]= 1.34018;   Fc[4]= 0.692012;   Fll[4]= 0.97259; 
              cafac[4]= 0.824207;  
              winpretag[4]= 12735.3 ;   wintag[4]= 1936.54 ;   woutpretag[4]= 10496.5;   wouttag[4]= 1741.49;  
            fbb[4]= 0.0959569;   fcc[4]= 0.158218;   fc[4]= 0.104342;   fll[4]= 0.641483;  
            fmcbb[4]= 0.0716001;   fmccc[4]= 0.118057;   fmcc[4]= 0.150781;   fmcll[4]= 0.659562;  
           Fbb[5]= 1.31044;   Fcc[5]= 1.31044;   Fc[5]= 0.676659;   Fll[5]= 0.951012; 
              cafac[5]= 0.837245;  
              winpretag[5]= 3353.98 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2808.1;   wouttag[5]= 601.334;  
            fbb[5]= 0.122797;   fcc[5]= 0.191493;   fc[5]= 0.0917902;   fll[5]= 0.59392;  
            fmcbb[5]= 0.0937061;   fmccc[5]= 0.146129;   fmcc[5]= 0.135652;   fmcll[5]= 0.624513;  
           Fbb[6]= 1.35851;   Fcc[6]= 1.35851;   Fc[6]= 0.70148;   Fll[6]= 0.985897; 
              cafac[6]= 0.806725;  
              winpretag[6]= 74081.9 ;   wintag[6]= 9187.63 ;   woutpretag[6]= 59763.7;   wouttag[6]= 7926.56;  
            fbb[6]= 0.075912;   fcc[6]= 0.140931;   fc[6]= 0.111908;   fll[6]= 0.671249;  
            fmcbb[6]= 0.0558787;   fmccc[6]= 0.103739;   fmcc[6]= 0.159531;   fmcll[6]= 0.680851;  
           Fbb[7]= 1.33387;   Fcc[7]= 1.33387;   Fc[7]= 0.688754;   Fll[7]= 0.968012; 
              cafac[7]= 0.826793;  
              winpretag[7]= 16089.3 ;   wintag[7]= 2582.3 ;   woutpretag[7]= 13302.5;   wouttag[7]= 2343.17;  
            fbb[7]= 0.101652;   fcc[7]= 0.165278;   fc[7]= 0.101679;   fll[7]= 0.631391;  
            fmcbb[7]= 0.0762083;   fmccc[7]= 0.123909;   fmcc[7]= 0.147627;   fmcll[7]= 0.652256;  
      }  
      else if ( idsys==59   || sysname=="musid_up" ) { 
           Fbb[1]= 1.40937;   Fcc[1]= 1.40937;   Fc[1]= 0.728115;   Fll[1]= 1.0232; 
              cafac[1]= 0.974348;  
              winpretag[1]= 1.08193e+06 ;   wintag[1]= 40096.9 ;   woutpretag[1]= 1.05418e+06;   wouttag[1]= 35404.7;  
            fbb[1]= 0.0169221;   fcc[1]= 0.0493866;   fc[1]= 0.102072;   fll[1]= 0.831619;  
            fmcbb[1]= 0.0120068;   fmccc[1]= 0.0350416;   fmcc[1]= 0.140187;   fmcll[1]= 0.812765;  
           Fbb[2]= 1.38868;   Fcc[2]= 1.38868;   Fc[2]= 0.717424;   Fll[2]= 1.00818; 
              cafac[2]= 0.884521;  
              winpretag[2]= 256891 ;   wintag[2]= 19908.1 ;   woutpretag[2]= 227226;   wouttag[2]= 17636.5;  
            fbb[2]= 0.0418393;   fcc[2]= 0.100588;   fc[2]= 0.116468;   fll[2]= 0.741105;  
            fmcbb[2]= 0.0301289;   fmccc[2]= 0.0724344;   fmcc[2]= 0.162341;   fmcll[2]= 0.735095;  
           Fbb[3]= 1.36508;   Fcc[3]= 1.36508;   Fc[3]= 0.70523;   Fll[3]= 0.991039; 
              cafac[3]= 0.80101;  
              winpretag[3]= 57992.6 ;   wintag[3]= 6605.33 ;   woutpretag[3]= 46452.7;   wouttag[3]= 5586.31;  
            fbb[3]= 0.0685793;   fcc[3]= 0.133973;   fc[3]= 0.114836;   fll[3]= 0.682612;  
            fmcbb[3]= 0.0502385;   fmccc[3]= 0.0981433;   fmcc[3]= 0.162834;   fmcll[3]= 0.688784;  
           Fbb[4]= 1.33978;   Fcc[4]= 1.33978;   Fc[4]= 0.69216;   Fll[4]= 0.972672; 
              cafac[4]= 0.824104;  
              winpretag[4]= 12735.3 ;   wintag[4]= 1936.54 ;   woutpretag[4]= 10495.2;   wouttag[4]= 1741.03;  
            fbb[4]= 0.0959281;   fcc[4]= 0.15817;   fc[4]= 0.104364;   fll[4]= 0.641537;  
            fmcbb[4]= 0.0716001;   fmccc[4]= 0.118057;   fmcc[4]= 0.150781;   fmcll[4]= 0.659562;  
           Fbb[5]= 1.31008;   Fcc[5]= 1.31008;   Fc[5]= 0.67682;   Fll[5]= 0.951116; 
              cafac[5]= 0.837284;  
              winpretag[5]= 3353.98 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2808.23;   wouttag[5]= 601.274;  
            fbb[5]= 0.122763;   fcc[5]= 0.191441;   fc[5]= 0.0918121;   fll[5]= 0.593984;  
            fmcbb[5]= 0.0937061;   fmccc[5]= 0.146129;   fmcc[5]= 0.135652;   fmcll[5]= 0.624513;  
           Fbb[6]= 1.35809;   Fcc[6]= 1.35809;   Fc[6]= 0.701619;   Fll[6]= 0.985965; 
              cafac[6]= 0.806718;  
              winpretag[6]= 74081.9 ;   wintag[6]= 9187.63 ;   woutpretag[6]= 59763.3;   wouttag[6]= 7925.46;  
            fbb[6]= 0.0758881;   fcc[6]= 0.140887;   fc[6]= 0.11193;   fll[6]= 0.671295;  
            fmcbb[6]= 0.0558787;   fmccc[6]= 0.103739;   fmcc[6]= 0.159531;   fmcll[6]= 0.680851;  
           Fbb[7]= 1.33348;   Fcc[7]= 1.33348;   Fc[7]= 0.688905;   Fll[7]= 0.968098; 
              cafac[7]= 0.826724;  
              winpretag[7]= 16089.3 ;   wintag[7]= 2582.3 ;   woutpretag[7]= 13301.4;   wouttag[7]= 2342.64;  
            fbb[7]= 0.101622;   fcc[7]= 0.16523;   fc[7]= 0.101701;   fll[7]= 0.631447;  
            fmcbb[7]= 0.0762083;   fmccc[7]= 0.123909;   fmcc[7]= 0.147627;   fmcll[7]= 0.652256;  
      }  
      else if ( idsys==60   || sysname=="musid_down" ) { 
           Fbb[1]= 1.40973;   Fcc[1]= 1.40973;   Fc[1]= 0.727882;   Fll[1]= 1.02322; 
              cafac[1]= 0.974303;  
              winpretag[1]= 1.08193e+06 ;   wintag[1]= 40096.9 ;   woutpretag[1]= 1.05413e+06;   wouttag[1]= 35399.9;  
            fbb[1]= 0.0169264;   fcc[1]= 0.0493993;   fc[1]= 0.102039;   fll[1]= 0.831635;  
            fmcbb[1]= 0.0120068;   fmccc[1]= 0.0350416;   fmcc[1]= 0.140186;   fmcll[1]= 0.812765;  
           Fbb[2]= 1.38902;   Fcc[2]= 1.38902;   Fc[2]= 0.717187;   Fll[2]= 1.00818; 
              cafac[2]= 0.884491;  
              winpretag[2]= 256892 ;   wintag[2]= 19908.9 ;   woutpretag[2]= 227218;   wouttag[2]= 17636.5;  
            fbb[2]= 0.0418495;   fcc[2]= 0.100609;   fc[2]= 0.116431;   fll[2]= 0.74111;  
            fmcbb[2]= 0.0301289;   fmccc[2]= 0.072432;   fmcc[2]= 0.162344;   fmcll[2]= 0.735095;  
           Fbb[3]= 1.36539;   Fcc[3]= 1.36539;   Fc[3]= 0.704985;   Fll[3]= 0.99103; 
              cafac[3]= 0.80097;  
              winpretag[3]= 57992.6 ;   wintag[3]= 6605.33 ;   woutpretag[3]= 46450.3;   wouttag[3]= 5586.29;  
            fbb[3]= 0.0685949;   fcc[3]= 0.134003;   fc[3]= 0.114796;   fll[3]= 0.682606;  
            fmcbb[3]= 0.0502385;   fmccc[3]= 0.0981433;   fmcc[3]= 0.162834;   fmcll[3]= 0.688784;  
           Fbb[4]= 1.34006;   Fcc[4]= 1.34006;   Fc[4]= 0.691908;   Fll[4]= 0.972648; 
              cafac[4]= 0.824092;  
              winpretag[4]= 12735.3 ;   wintag[4]= 1936.54 ;   woutpretag[4]= 10495.1;   wouttag[4]= 1741.12;  
            fbb[4]= 0.0959484;   fcc[4]= 0.158204;   fc[4]= 0.104326;   fll[4]= 0.641521;  
            fmcbb[4]= 0.0716001;   fmccc[4]= 0.118057;   fmcc[4]= 0.150781;   fmcll[4]= 0.659562;  
           Fbb[5]= 1.31034;   Fcc[5]= 1.31034;   Fc[5]= 0.676562;   Fll[5]= 0.951075; 
              cafac[5]= 0.837232;  
              winpretag[5]= 3353.98 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2808.06;   wouttag[5]= 601.286;  
            fbb[5]= 0.122787;   fcc[5]= 0.191478;   fc[5]= 0.091777;   fll[5]= 0.593959;  
            fmcbb[5]= 0.0937061;   fmccc[5]= 0.146129;   fmcc[5]= 0.135652;   fmcll[5]= 0.624513;  
           Fbb[6]= 1.35839;   Fcc[6]= 1.35839;   Fc[6]= 0.701372;   Fll[6]= 0.985952; 
              cafac[6]= 0.806682;  
              winpretag[6]= 74081.9 ;   wintag[6]= 9187.63 ;   woutpretag[6]= 59760.6;   wouttag[6]= 7925.55;  
            fbb[6]= 0.075905;   fcc[6]= 0.140918;   fc[6]= 0.111891;   fll[6]= 0.671286;  
            fmcbb[6]= 0.0558787;   fmccc[6]= 0.103739;   fmcc[6]= 0.159531;   fmcll[6]= 0.680851;  
           Fbb[7]= 1.33375;   Fcc[7]= 1.33375;   Fc[7]= 0.688652;   Fll[7]= 0.96807; 
              cafac[7]= 0.826702;  
              winpretag[7]= 16089.3 ;   wintag[7]= 2582.3 ;   woutpretag[7]= 13301.1;   wouttag[7]= 2342.75;  
            fbb[7]= 0.101643;   fcc[7]= 0.165264;   fc[7]= 0.101664;   fll[7]= 0.631429;  
            fmcbb[7]= 0.0762083;   fmccc[7]= 0.123909;   fmcc[7]= 0.147627;   fmcll[7]= 0.652256;  
      }  
      else if ( idsys==57   || sysname=="mscale_one" ) { 
           Fbb[1]= 1.39034;   Fcc[1]= 1.39034;   Fc[1]= 0.745819;   Fll[1]= 1.02125; 
              cafac[1]= 0.97764;  
              winpretag[1]= 1.08194e+06 ;   wintag[1]= 40097.6 ;   woutpretag[1]= 1.05775e+06;   wouttag[1]= 35801.9;  
            fbb[1]= 0.0166933;   fcc[1]= 0.0487188;   fc[1]= 0.104556;   fll[1]= 0.830032;  
            fmcbb[1]= 0.0120067;   fmccc[1]= 0.0350409;   fmcc[1]= 0.14019;   fmcll[1]= 0.812763;  
           Fbb[2]= 1.37062;   Fcc[2]= 1.37062;   Fc[2]= 0.735238;   Fll[2]= 1.00676; 
              cafac[2]= 0.888109;  
              winpretag[2]= 256890 ;   wintag[2]= 19909.2 ;   woutpretag[2]= 228147;   wouttag[2]= 17739.5;  
            fbb[2]= 0.0412997;   fcc[2]= 0.0992801;   fc[2]= 0.119358;   fll[2]= 0.740062;  
            fmcbb[2]= 0.0301323;   fmccc[2]= 0.0724347;   fmcc[2]= 0.162339;   fmcll[2]= 0.735094;  
           Fbb[3]= 1.34832;   Fcc[3]= 1.34832;   Fc[3]= 0.72328;   Fll[3]= 0.990385; 
              cafac[3]= 0.804292;  
              winpretag[3]= 57992.9 ;   wintag[3]= 6605.99 ;   woutpretag[3]= 46643.2;   wouttag[3]= 5602.95;  
            fbb[3]= 0.0677374;   fcc[3]= 0.132328;   fc[3]= 0.117783;   fll[3]= 0.682152;  
            fmcbb[3]= 0.0502382;   fmccc[3]= 0.0981427;   fmcc[3]= 0.162845;   fmcll[3]= 0.688774;  
           Fbb[4]= 1.32449;   Fcc[4]= 1.32449;   Fc[4]= 0.710494;   Fll[4]= 0.972877; 
              cafac[4]= 0.827535;  
              winpretag[4]= 12735.3 ;   wintag[4]= 1936.54 ;   woutpretag[4]= 10538.9;   wouttag[4]= 1743.23;  
            fbb[4]= 0.0948334;   fcc[4]= 0.156365;   fc[4]= 0.107129;   fll[4]= 0.641672;  
            fmcbb[4]= 0.0716001;   fmccc[4]= 0.118057;   fmcc[4]= 0.150781;   fmcll[4]= 0.659562;  
           Fbb[5]= 1.29647;   Fcc[5]= 1.29647;   Fc[5]= 0.695463;   Fll[5]= 0.952295; 
              cafac[5]= 0.841145;  
              winpretag[5]= 3353.98 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2821.18;   wouttag[5]= 601.745;  
            fbb[5]= 0.121487;   fcc[5]= 0.189451;   fc[5]= 0.094341;   fll[5]= 0.594721;  
            fmcbb[5]= 0.0937061;   fmccc[5]= 0.146129;   fmcc[5]= 0.135652;   fmcll[5]= 0.624513;  
           Fbb[6]= 1.34174;   Fcc[6]= 1.34174;   Fc[6]= 0.71975;   Fll[6]= 0.985552; 
              cafac[6]= 0.810078;  
              winpretag[6]= 74082.2 ;   wintag[6]= 9188.3 ;   woutpretag[6]= 60012.4;   wouttag[6]= 7944.16;  
            fbb[6]= 0.0749746;   fcc[6]= 0.139191;   fc[6]= 0.114829;   fll[6]= 0.671006;  
            fmcbb[6]= 0.0558785;   fmccc[6]= 0.103739;   fmcc[6]= 0.15954;   fmcll[6]= 0.680843;  
           Fbb[7]= 1.31855;   Fcc[7]= 1.31855;   Fc[7]= 0.707307;   Fll[7]= 0.968513; 
              cafac[7]= 0.830274;  
              winpretag[7]= 16089.3 ;   wintag[7]= 2582.3 ;   woutpretag[7]= 13358.5;   wouttag[7]= 2345.24;  
            fbb[7]= 0.100484;   fcc[7]= 0.16338;   fc[7]= 0.104418;   fll[7]= 0.631718;  
            fmcbb[7]= 0.0762083;   fmccc[7]= 0.123909;   fmcc[7]= 0.147627;   fmcll[7]= 0.652256;  
      }  
      else if ( idsys==55   || sysname=="metpileup_up" ) { 
           Fbb[1]= 1.40167;   Fcc[1]= 1.40167;   Fc[1]= 0.733401;   Fll[1]= 1.02279; 
              cafac[1]= 0.973658;  
              winpretag[1]= 1.08491e+06 ;   wintag[1]= 40208.8 ;   woutpretag[1]= 1.05633e+06;   wouttag[1]= 35535.2;  
            fbb[1]= 0.0167882;   fcc[1]= 0.0491208;   fc[1]= 0.102909;   fll[1]= 0.831182;  
            fmcbb[1]= 0.0119773;   fmccc[1]= 0.0350445;   fmcc[1]= 0.140318;   fmcll[1]= 0.81266;  
           Fbb[2]= 1.38139;   Fcc[2]= 1.38139;   Fc[2]= 0.722792;   Fll[2]= 1.00799; 
              cafac[2]= 0.886889;  
              winpretag[2]= 257255 ;   wintag[2]= 19957.2 ;   woutpretag[2]= 228156;   wouttag[2]= 17726.4;  
            fbb[2]= 0.0416474;   fcc[2]= 0.100091;   fc[2]= 0.117358;   fll[2]= 0.740903;  
            fmcbb[2]= 0.0301489;   fmccc[2]= 0.072457;   fmcc[2]= 0.162367;   fmcll[2]= 0.735027;  
           Fbb[3]= 1.35839;   Fcc[3]= 1.35839;   Fc[3]= 0.710758;   Fll[3]= 0.991213; 
              cafac[3]= 0.805543;  
              winpretag[3]= 58045.5 ;   wintag[3]= 6610.97 ;   woutpretag[3]= 46758.1;   wouttag[3]= 5616.99;  
            fbb[3]= 0.0682405;   fcc[3]= 0.13319;   fc[3]= 0.115717;   fll[3]= 0.682853;  
            fmcbb[3]= 0.0502362;   fmccc[3]= 0.0980495;   fmcc[3]= 0.162807;   fmcll[3]= 0.688907;  
           Fbb[4]= 1.33363;   Fcc[4]= 1.33363;   Fc[4]= 0.697799;   Fll[4]= 0.973141; 
              cafac[4]= 0.827429;  
              winpretag[4]= 12741.8 ;   wintag[4]= 1938.78 ;   woutpretag[4]= 10543;   wouttag[4]= 1746.51;  
            fbb[4]= 0.0955197;   fcc[4]= 0.157444;   fc[4]= 0.10522;   fll[4]= 0.641816;  
            fmcbb[4]= 0.0716241;   fmccc[4]= 0.118057;   fmcc[4]= 0.150788;   fmcll[4]= 0.659531;  
           Fbb[5]= 1.30432;   Fcc[5]= 1.30432;   Fc[5]= 0.682465;   Fll[5]= 0.951755; 
              cafac[5]= 0.837255;  
              winpretag[5]= 3356.58 ;   wintag[5]= 643.443 ;   woutpretag[5]= 2810.31;   wouttag[5]= 597.991;  
            fbb[5]= 0.122794;   fcc[5]= 0.190515;   fc[5]= 0.0923612;   fll[5]= 0.594329;  
            fmcbb[5]= 0.0941444;   fmccc[5]= 0.146065;   fmcc[5]= 0.135335;   fmcll[5]= 0.624456;  
           Fbb[6]= 1.35154;   Fcc[6]= 1.35154;   Fc[6]= 0.707174;   Fll[6]= 0.986214; 
              cafac[6]= 0.810787;  
              winpretag[6]= 74143.9 ;   wintag[6]= 9193.2 ;   woutpretag[6]= 60114.9;   wouttag[6]= 7959.47;  
            fbb[6]= 0.0755506;   fcc[6]= 0.140103;   fc[6]= 0.112793;   fll[6]= 0.671554;  
            fmcbb[6]= 0.0558996;   fmccc[6]= 0.103662;   fmcc[6]= 0.159498;   fmcll[6]= 0.680941;  
           Fbb[7]= 1.32741;   Fcc[7]= 1.32741;   Fc[7]= 0.694545;   Fll[7]= 0.968603; 
              cafac[7]= 0.82927;  
              winpretag[7]= 16098.4 ;   wintag[7]= 2582.23 ;   woutpretag[7]= 13349.9;   wouttag[7]= 2345.01;  
            fbb[7]= 0.101307;   fcc[7]= 0.164462;   fc[7]= 0.102491;   fll[7]= 0.63174;  
            fmcbb[7]= 0.0763197;   fmccc[7]= 0.123897;   fmcc[7]= 0.147566;   fmcll[7]= 0.652218;  
      }  
      else if ( idsys==56   || sysname=="metpileup_down" ) { 
           Fbb[1]= 1.416;   Fcc[1]= 1.416;   Fc[1]= 0.740768;   Fll[1]= 1.02056; 
              cafac[1]= 0.975431;  
              winpretag[1]= 1.07894e+06 ;   wintag[1]= 39922.4 ;   woutpretag[1]= 1.05243e+06;   wouttag[1]= 35595.5;  
            fbb[1]= 0.0170435;   fcc[1]= 0.0496991;   fc[1]= 0.103782;   fll[1]= 0.829475;  
            fmcbb[1]= 0.0120364;   fmccc[1]= 0.0350982;   fmcc[1]= 0.140101;   fmcll[1]= 0.812765;  
           Fbb[2]= 1.39411;   Fcc[2]= 1.39411;   Fc[2]= 0.729319;   Fll[2]= 1.00479; 
              cafac[2]= 0.886436;  
              winpretag[2]= 256620 ;   wintag[2]= 19867.1 ;   woutpretag[2]= 227477;   wouttag[2]= 17749.2;  
            fbb[2]= 0.0420931;   fcc[2]= 0.100978;   fc[2]= 0.118454;   fll[2]= 0.738475;  
            fmcbb[2]= 0.0301934;   fmccc[2]= 0.0724316;   fmcc[2]= 0.162417;   fmcll[2]= 0.734958;  
           Fbb[3]= 1.36981;   Fcc[3]= 1.36981;   Fc[3]= 0.716607;   Fll[3]= 0.987272; 
              cafac[3]= 0.807956;  
              winpretag[3]= 57904.1 ;   wintag[3]= 6585.72 ;   woutpretag[3]= 46784;   wouttag[3]= 5646.65;  
            fbb[3]= 0.0688037;   fcc[3]= 0.134533;   fc[3]= 0.116645;   fll[3]= 0.680018;  
            fmcbb[3]= 0.0502285;   fmccc[3]= 0.0982129;   fmcc[3]= 0.162774;   fmcll[3]= 0.688784;  
           Fbb[4]= 1.34428;   Fcc[4]= 1.34428;   Fc[4]= 0.70325;   Fll[4]= 0.968871; 
              cafac[4]= 0.822448;  
              winpretag[4]= 12731.9 ;   wintag[4]= 1936.48 ;   woutpretag[4]= 10471.3;   wouttag[4]= 1743.94;  
            fbb[4]= 0.0961589;   fcc[4]= 0.158878;   fc[4]= 0.106151;   fll[4]= 0.638812;  
            fmcbb[4]= 0.0715317;   fmccc[4]= 0.118188;   fmcc[4]= 0.150944;   fmcll[4]= 0.659336;  
           Fbb[5]= 1.31465;   Fcc[5]= 1.31465;   Fc[5]= 0.68775;   Fll[5]= 0.947517; 
              cafac[5]= 0.836939;  
              winpretag[5]= 3350.35 ;   wintag[5]= 648.715 ;   woutpretag[5]= 2804.04;   wouttag[5]= 604.827;  
            fbb[5]= 0.122577;   fcc[5]= 0.192346;   fc[5]= 0.0938935;   fll[5]= 0.591183;  
            fmcbb[5]= 0.0932391;   fmccc[5]= 0.14631;   fmcc[5]= 0.136523;   fmcll[5]= 0.623929;  
           Fbb[6]= 1.36277;   Fcc[6]= 1.36277;   Fc[6]= 0.712922;   Fll[6]= 0.982196; 
              cafac[6]= 0.811706;  
              winpretag[6]= 73986.4 ;   wintag[6]= 9170.92 ;   woutpretag[6]= 60055.2;   wouttag[6]= 7996.59;  
            fbb[6]= 0.0761;   fcc[6]= 0.141494;   fc[6]= 0.113747;   fll[6]= 0.668659;  
            fmcbb[6]= 0.0558421;   fmccc[6]= 0.103828;   fmcc[6]= 0.15955;   fmcll[6]= 0.68078;  
           Fbb[7]= 1.338;   Fcc[7]= 1.338;   Fc[7]= 0.699964;   Fll[7]= 0.964343; 
              cafac[7]= 0.825373;  
              winpretag[7]= 16082.2 ;   wintag[7]= 2585.2 ;   woutpretag[7]= 13273.9;   wouttag[7]= 2349.02;  
            fbb[7]= 0.10176;   fcc[7]= 0.165974;   fc[7]= 0.103552;   fll[7]= 0.628713;  
            fmcbb[7]= 0.0760539;   fmccc[7]= 0.124046;   fmcc[7]= 0.14794;   fmcll[7]= 0.65196;  
      }  
      else if ( idsys==87   || sysname=="jvf_up" ) { 
           Fbb[1]= 1.43893;   Fcc[1]= 1.43893;   Fc[1]= 0.748364;   Fll[1]= 1.01794; 
              cafac[1]= 0.986271;  
              winpretag[1]= 1.07776e+06 ;   wintag[1]= 39919.3 ;   woutpretag[1]= 1.06296e+06;   wouttag[1]= 36276.8;  
            fbb[1]= 0.0173474;   fcc[1]= 0.0504839;   fc[1]= 0.104903;   fll[1]= 0.827265;  
            fmcbb[1]= 0.0120558;   fmccc[1]= 0.0350843;   fmcc[1]= 0.140177;   fmcll[1]= 0.812683;  
           Fbb[2]= 1.41448;   Fcc[2]= 1.41448;   Fc[2]= 0.735646;   Fll[2]= 1.00064; 
              cafac[2]= 0.885568;  
              winpretag[2]= 255169 ;   wintag[2]= 19654.5 ;   woutpretag[2]= 225969;   wouttag[2]= 17679;  
            fbb[2]= 0.0425904;   fcc[2]= 0.102181;   fc[2]= 0.119371;   fll[2]= 0.735857;  
            fmcbb[2]= 0.0301103;   fmccc[2]= 0.0722397;   fmcc[2]= 0.162266;   fmcll[2]= 0.735383;  
           Fbb[3]= 1.38831;   Fcc[3]= 1.38831;   Fc[3]= 0.722037;   Fll[3]= 0.982133; 
              cafac[3]= 0.818141;  
              winpretag[3]= 57244.3 ;   wintag[3]= 6494.35 ;   woutpretag[3]= 46833.9;   wouttag[3]= 5685.36;  
            fbb[3]= 0.0699516;   fcc[3]= 0.136256;   fc[3]= 0.117878;   fll[3]= 0.675914;  
            fmcbb[3]= 0.0503862;   fmccc[3]= 0.0981451;   fmcc[3]= 0.163258;   fmcll[3]= 0.688211;  
           Fbb[4]= 1.36109;   Fcc[4]= 1.36109;   Fc[4]= 0.707882;   Fll[4]= 0.962879; 
              cafac[4]= 0.842048;  
              winpretag[4]= 12582.6 ;   wintag[4]= 1914.4 ;   woutpretag[4]= 10595.1;   wouttag[4]= 1779.32;  
            fbb[4]= 0.0973336;   fcc[4]= 0.160785;   fc[4]= 0.106591;   fll[4]= 0.63529;  
            fmcbb[4]= 0.0715113;   fmccc[4]= 0.11813;   fmcc[4]= 0.150577;   fmcll[4]= 0.659782;  
           Fbb[5]= 1.32975;   Fcc[5]= 1.32975;   Fc[5]= 0.691582;   Fll[5]= 0.940707; 
              cafac[5]= 0.834664;  
              winpretag[5]= 3309.76 ;   wintag[5]= 633.525 ;   woutpretag[5]= 2762.53;   wouttag[5]= 594.834;  
            fbb[5]= 0.1242;   fcc[5]= 0.193439;   fc[5]= 0.0933807;   fll[5]= 0.588981;  
            fmcbb[5]= 0.0934006;   fmccc[5]= 0.14547;   fmcc[5]= 0.135025;   fmcll[5]= 0.626105;  
           Fbb[6]= 1.38081;   Fcc[6]= 1.38081;   Fc[6]= 0.718136;   Fll[6]= 0.976826; 
              cafac[6]= 0.822872;  
              winpretag[6]= 73136.6 ;   wintag[6]= 9042.28 ;   woutpretag[6]= 60182.1;   wouttag[6]= 8060.26;  
            fbb[6]= 0.0772799;   fcc[6]= 0.143224;   fc[6]= 0.114757;   fll[6]= 0.664739;  
            fmcbb[6]= 0.0559672;   fmccc[6]= 0.103725;   fmcc[6]= 0.159799;   fmcll[6]= 0.680509;  
           Fbb[7]= 1.35444;   Fcc[7]= 1.35444;   Fc[7]= 0.704425;   Fll[7]= 0.958176; 
              cafac[7]= 0.839811;  
              winpretag[7]= 15892.3 ;   wintag[7]= 2547.93 ;   woutpretag[7]= 13346.5;   wouttag[7]= 2375.54;  
            fbb[7]= 0.103033;   fcc[7]= 0.167712;   fc[7]= 0.103789;   fll[7]= 0.625467;  
            fmcbb[7]= 0.07607;   fmccc[7]= 0.123824;   fmcc[7]= 0.147338;   fmcll[7]= 0.652768;  
      }  
      else if ( idsys==88   || sysname=="jvf_down" ) { 
           Fbb[1]= 1.42998;   Fcc[1]= 1.42998;   Fc[1]= 0.729733;   Fll[1]= 1.02159; 
              cafac[1]= 0.972428;  
              winpretag[1]= 1.08882e+06 ;   wintag[1]= 40344.6 ;   woutpretag[1]= 1.0588e+06;   wouttag[1]= 35732.5;  
            fbb[1]= 0.0172448;   fcc[1]= 0.050152;   fc[1]= 0.102111;   fll[1]= 0.830492;  
            fmcbb[1]= 0.0120595;   fmccc[1]= 0.0350718;   fmcc[1]= 0.13993;   fmcll[1]= 0.812939;  
           Fbb[2]= 1.40743;   Fcc[2]= 1.40743;   Fc[2]= 0.718224;   Fll[2]= 1.00548; 
              cafac[2]= 0.868421;  
              winpretag[2]= 259110 ;   wintag[2]= 19984 ;   woutpretag[2]= 225017;   wouttag[2]= 17475.9;  
            fbb[2]= 0.0423271;   fcc[2]= 0.101559;   fc[2]= 0.116445;   fll[2]= 0.739669;  
            fmcbb[2]= 0.0300741;   fmccc[2]= 0.0721595;   fmcc[2]= 0.162129;   fmcll[2]= 0.735638;  
           Fbb[3]= 1.38219;   Fcc[3]= 1.38219;   Fc[3]= 0.705348;   Fll[3]= 0.987454; 
              cafac[3]= 0.797936;  
              winpretag[3]= 58455.8 ;   wintag[3]= 6645.61 ;   woutpretag[3]= 46644;   wouttag[3]= 5633.71;  
            fbb[3]= 0.069576;   fcc[3]= 0.135501;   fc[3]= 0.115067;   fll[3]= 0.679856;  
            fmcbb[3]= 0.0503374;   fmccc[3]= 0.0980331;   fmcc[3]= 0.163136;   fmcll[3]= 0.688494;  
           Fbb[4]= 1.35539;   Fcc[4]= 1.35539;   Fc[4]= 0.691669;   Fll[4]= 0.968304; 
              cafac[4]= 0.814974;  
              winpretag[4]= 12918.8 ;   wintag[4]= 1969.46 ;   woutpretag[4]= 10528.4;   wouttag[4]= 1760.93;  
            fbb[4]= 0.0967742;   fcc[4]= 0.160081;   fc[4]= 0.10416;   fll[4]= 0.638985;  
            fmcbb[4]= 0.0713996;   fmccc[4]= 0.118107;   fmcc[4]= 0.150593;   fmcll[4]= 0.659901;  
           Fbb[5]= 1.32442;   Fcc[5]= 1.32442;   Fc[5]= 0.675864;   Fll[5]= 0.946178; 
              cafac[5]= 0.801262;  
              winpretag[5]= 3426.17 ;   wintag[5]= 658.192 ;   woutpretag[5]= 2745.26;   wouttag[5]= 590.032;  
            fbb[5]= 0.123702;   fcc[5]= 0.192549;   fc[5]= 0.0912513;   fll[5]= 0.592497;  
            fmcbb[5]= 0.0934013;   fmccc[5]= 0.145384;   fmcc[5]= 0.135014;   fmcll[5]= 0.6262;  
           Fbb[6]= 1.37475;   Fcc[6]= 1.37475;   Fc[6]= 0.70155;   Fll[6]= 0.982137; 
              cafac[6]= 0.800732;  
              winpretag[6]= 74800.7 ;   wintag[6]= 9273.25 ;   woutpretag[6]= 59895.3;   wouttag[6]= 7989.93;  
            fbb[6]= 0.0769139;   fcc[6]= 0.142519;   fc[6]= 0.112024;   fll[6]= 0.668543;  
            fmcbb[6]= 0.0559475;   fmccc[6]= 0.103669;   fmcc[6]= 0.159681;   fmcll[6]= 0.680702;  
           Fbb[7]= 1.34878;   Fcc[7]= 1.34878;   Fc[7]= 0.688295;   Fll[7]= 0.963581; 
              cafac[7]= 0.811262;  
              winpretag[7]= 16344.9 ;   wintag[7]= 2627.65 ;   woutpretag[7]= 13260;   wouttag[7]= 2352.74;  
            fbb[7]= 0.102523;   fcc[7]= 0.167012;   fc[7]= 0.101405;   fll[7]= 0.629061;  
            fmcbb[7]= 0.0760115;   fmccc[7]= 0.123825;   fmcc[7]= 0.147327;   fmcll[7]= 0.652837;  
      }  
      else if ( idsys==89   || sysname=="pdfs_up" ) { 
           Fbb[1]= 1.30823;   Fcc[1]= 1.30823;   Fc[1]= 0.680027;   Fll[1]= 1.03928; 
              cafac[1]= 0.961799;  
              winpretag[1]= 1.11666e+06 ;   wintag[1]= 42227.4 ;   woutpretag[1]= 1.074e+06;   wouttag[1]= 35110.5;  
            fbb[1]= 0.0155929;   fcc[1]= 0.0453648;   fc[1]= 0.0980797;   fll[1]= 0.840963;  
            fmcbb[1]= 0.0119191;   fmccc[1]= 0.0346764;   fmcc[1]= 0.144229;   fmcll[1]= 0.809175;  
           Fbb[2]= 1.29831;   Fcc[2]= 1.29831;   Fc[2]= 0.674872;   Fll[2]= 1.0314; 
              cafac[2]= 0.864437;  
              winpretag[2]= 264867 ;   wintag[2]= 20652.5 ;   woutpretag[2]= 228961;   wouttag[2]= 17107.7;  
            fbb[2]= 0.0389309;   fcc[2]= 0.0934907;   fc[2]= 0.110974;   fll[2]= 0.756605;  
            fmcbb[2]= 0.0299858;   fmccc[2]= 0.0720094;   fmcc[2]= 0.164437;   fmcll[2]= 0.733568;  
           Fbb[3]= 1.28323;   Fcc[3]= 1.28323;   Fc[3]= 0.667034;   Fll[3]= 1.01943; 
              cafac[3]= 0.755656;  
              winpretag[3]= 60033.3 ;   wintag[3]= 6901.47 ;   woutpretag[3]= 45364.5;   wouttag[3]= 5269.8;  
            fbb[3]= 0.0642725;   fcc[3]= 0.125096;   fc[3]= 0.110462;   fll[3]= 0.700169;  
            fmcbb[3]= 0.0500864;   fmccc[3]= 0.0974853;   fmcc[3]= 0.165601;   fmcll[3]= 0.686827;  
           Fbb[4]= 1.26505;   Fcc[4]= 1.26505;   Fc[4]= 0.65758;   Fll[4]= 1.00498; 
              cafac[4]= 0.785432;  
              winpretag[4]= 13316.2 ;   wintag[4]= 2066.49 ;   woutpretag[4]= 10458.9;   wouttag[4]= 1696.99;  
            fbb[4]= 0.0904689;   fcc[4]= 0.148201;   fc[4]= 0.102298;   fll[4]= 0.659031;  
            fmcbb[4]= 0.0715143;   fmccc[4]= 0.117151;   fmcc[4]= 0.155568;   fmcll[4]= 0.655767;  
           Fbb[5]= 1.24535;   Fcc[5]= 1.24535;   Fc[5]= 0.647343;   Fll[5]= 0.989331; 
              cafac[5]= 0.797774;  
              winpretag[5]= 3569.1 ;   wintag[5]= 712.377 ;   woutpretag[5]= 2847.33;   wouttag[5]= 604.229;  
            fbb[5]= 0.116696;   fcc[5]= 0.181862;   fc[5]= 0.0959864;   fll[5]= 0.605456;  
            fmcbb[5]= 0.0937055;   fmccc[5]= 0.146032;   fmcc[5]= 0.148278;   fmcll[5]= 0.611985;  
           Fbb[6]= 1.27825;   Fcc[6]= 1.27825;   Fc[6]= 0.664443;   Fll[6]= 1.01547; 
              cafac[6]= 0.763117;  
              winpretag[6]= 76918.5 ;   wintag[6]= 9680.34 ;   woutpretag[6]= 58697.9;   wouttag[6]= 7560.39;  
            fbb[6]= 0.0713517;   fcc[6]= 0.131842;   fc[6]= 0.108344;   fll[6]= 0.688462;  
            fmcbb[6]= 0.0558199;   fmccc[6]= 0.103142;   fmcc[6]= 0.16306;   fmcll[6]= 0.677977;  
           Fbb[7]= 1.26083;   Fcc[7]= 1.26083;   Fc[7]= 0.65539;   Fll[7]= 1.00163; 
              cafac[7]= 0.788039;  
              winpretag[7]= 16885.3 ;   wintag[7]= 2778.87 ;   woutpretag[7]= 13306.3;   wouttag[7]= 2301.23;  
            fbb[7]= 0.0960816;   fcc[7]= 0.155405;   fc[7]= 0.100947;   fll[7]= 0.647566;  
            fmcbb[7]= 0.0762049;   fmccc[7]= 0.123256;   fmcc[7]= 0.154027;   fmcll[7]= 0.646513;  
      }  
      else if ( idsys==90   || sysname=="pdfs_down" ) { 
           Fbb[1]= 1.46025;   Fcc[1]= 1.46025;   Fc[1]= 0.86527;   Fll[1]= 0.995629; 
              cafac[1]= 1.00215;  
              winpretag[1]= 1.04721e+06 ;   wintag[1]= 37966.5 ;   woutpretag[1]= 1.04946e+06;   wouttag[1]= 37603.5;  
            fbb[1]= 0.0176696;   fcc[1]= 0.0517382;   fc[1]= 0.11757;   fll[1]= 0.813023;  
            fmcbb[1]= 0.0121004;   fmccc[1]= 0.035431;   fmcc[1]= 0.135876;   fmcll[1]= 0.816592;  
           Fbb[2]= 1.42786;   Fcc[2]= 1.42786;   Fc[2]= 0.846073;   Fll[2]= 0.97354; 
              cafac[2]= 0.923193;  
              winpretag[2]= 248915 ;   wintag[2]= 19163.6 ;   woutpretag[2]= 229796;   wouttag[2]= 18811.2;  
            fbb[2]= 0.0432374;   fcc[2]= 0.104065;   fc[2]= 0.135467;   fll[2]= 0.717231;  
            fmcbb[2]= 0.0302813;   fmccc[2]= 0.0728821;   fmcc[2]= 0.160112;   fmcll[2]= 0.736724;  
           Fbb[3]= 1.39853;   Fcc[3]= 1.39853;   Fc[3]= 0.828696;   Fll[3]= 0.953545; 
              cafac[3]= 0.866897;  
              winpretag[3]= 55952 ;   wintag[3]= 6309.19 ;   woutpretag[3]= 48504.6;   wouttag[3]= 6056.24;  
            fbb[3]= 0.0704883;   fcc[3]= 0.138243;   fc[3]= 0.13248;   fll[3]= 0.658789;  
            fmcbb[3]= 0.0504018;   fmccc[3]= 0.0988492;   fmcc[3]= 0.159865;   fmcll[3]= 0.690884;  
           Fbb[4]= 1.37077;   Fcc[4]= 1.37077;   Fc[4]= 0.812246;   Fll[4]= 0.934616; 
              cafac[4]= 0.882176;  
              winpretag[4]= 12154.5 ;   wintag[4]= 1806.59 ;   woutpretag[4]= 10722.4;   wouttag[4]= 1807.31;  
            fbb[4]= 0.098276;   fcc[4]= 0.16319;   fc[4]= 0.118211;   fll[4]= 0.620323;  
            fmcbb[4]= 0.0716941;   fmccc[4]= 0.11905;   fmcc[4]= 0.145536;   fmcll[4]= 0.663719;  
           Fbb[5]= 1.33808;   Fcc[5]= 1.33808;   Fc[5]= 0.792879;   Fll[5]= 0.912332; 
              cafac[5]= 0.902008;  
              winpretag[5]= 3138.86 ;   wintag[5]= 579.144 ;   woutpretag[5]= 2831.28;   wouttag[5]= 604.173;  
            fbb[5]= 0.125388;   fcc[5]= 0.195679;   fc[5]= 0.0961732;   fll[5]= 0.582761;  
            fmcbb[5]= 0.0937068;   fmccc[5]= 0.146238;   fmcc[5]= 0.121296;   fmcll[5]= 0.638759;  
           Fbb[6]= 1.39095;   Fcc[6]= 1.39095;   Fc[6]= 0.824208;   Fll[6]= 0.94838; 
              cafac[6]= 0.870992;  
              winpretag[6]= 71245.3 ;   wintag[6]= 8694.92 ;   woutpretag[6]= 62054.1;   wouttag[6]= 8469.01;  
            fbb[6]= 0.077813;   fcc[6]= 0.145192;   fc[6]= 0.128347;   fll[6]= 0.648648;  
            fmcbb[6]= 0.0559421;   fmccc[6]= 0.104383;   fmcc[6]= 0.155721;   fmcll[6]= 0.683953;  
           Fbb[7]= 1.36393;   Fcc[7]= 1.36393;   Fc[7]= 0.808194;   Fll[7]= 0.929954; 
              cafac[7]= 0.886201;  
              winpretag[7]= 15293.3 ;   wintag[7]= 2385.73 ;   woutpretag[7]= 13553;   wouttag[7]= 2411.55;  
            fbb[7]= 0.103948;   fcc[7]= 0.169987;   fc[7]= 0.113601;   fll[7]= 0.612465;  
            fmcbb[7]= 0.0762121;   fmccc[7]= 0.12463;   fmcc[7]= 0.140561;   fmcll[7]= 0.658597;  
      }  
      else if ( idsys==91   || sysname=="pdfa_up" ) { 
           Fbb[1]= 1.18961;   Fcc[1]= 1.18961;   Fc[1]= 0.978189;   Fll[1]= 0.992769; 
              cafac[1]= 0.920972;  
              winpretag[1]= 1.10168e+06 ;   wintag[1]= 40878.1 ;   woutpretag[1]= 1.01462e+06;   wouttag[1]= 38111.9;  
            fbb[1]= 0.014317;   fcc[1]= 0.0417263;   fc[1]= 0.137058;   fll[1]= 0.806899;  
            fmcbb[1]= 0.012035;   fmccc[1]= 0.0350754;   fmcc[1]= 0.140114;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.17714;   Fcc[2]= 1.17714;   Fc[2]= 0.967927;   Fll[2]= 0.982355; 
              cafac[2]= 0.84872;  
              winpretag[2]= 261596 ;   wintag[2]= 20308.7 ;   woutpretag[2]= 222021;   wouttag[2]= 17885.9;  
            fbb[2]= 0.0354968;   fcc[2]= 0.0853172;   fc[2]= 0.1574;   fll[2]= 0.721786;  
            fmcbb[2]= 0.0301552;   fmccc[2]= 0.0724786;   fmcc[2]= 0.162615;   fmcll[2]= 0.734751;  
           Fbb[3]= 1.16673;   Fcc[3]= 1.16673;   Fc[3]= 0.959371;   Fll[3]= 0.973671; 
              cafac[3]= 0.74727;  
              winpretag[3]= 59200.1 ;   wintag[3]= 6752.58 ;   woutpretag[3]= 44238.4;   wouttag[3]= 5311.71;  
            fbb[3]= 0.0586074;   fcc[3]= 0.114608;   fc[3]= 0.15654;   fll[3]= 0.670244;  
            fmcbb[3]= 0.0502322;   fmccc[3]= 0.0982304;   fmcc[3]= 0.163169;   fmcll[3]= 0.688368;  
           Fbb[4]= 1.15728;   Fcc[4]= 1.15728;   Fc[4]= 0.951598;   Fll[4]= 0.965782; 
              cafac[4]= 0.768097;  
              winpretag[4]= 13079.4 ;   wintag[4]= 1998.79 ;   woutpretag[4]= 10046.2;   wouttag[4]= 1630.05;  
            fbb[4]= 0.0828863;   fcc[4]= 0.136864;   fc[4]= 0.143843;   fll[4]= 0.636406;  
            fmcbb[4]= 0.0716218;   fmccc[4]= 0.118264;   fmcc[4]= 0.15116;   fmcll[4]= 0.658954;  
           Fbb[5]= 1.14621;   Fcc[5]= 1.14621;   Fc[5]= 0.942496;   Fll[5]= 0.956545; 
              cafac[5]= 0.762661;  
              winpretag[5]= 3485.13 ;   wintag[5]= 677.6 ;   woutpretag[5]= 2657.97;   wouttag[5]= 551.007;  
            fbb[5]= 0.10703;   fcc[5]= 0.167343;   fc[5]= 0.1305;   fll[5]= 0.595128;  
            fmcbb[5]= 0.0933773;   fmccc[5]= 0.145997;   fmcc[5]= 0.138462;   fmcll[5]= 0.622164;  
           Fbb[6]= 1.16413;   Fcc[6]= 1.16413;   Fc[6]= 0.957233;   Fll[6]= 0.971501; 
              cafac[6]= 0.751759;  
              winpretag[6]= 75764.6 ;   wintag[6]= 9428.97 ;   woutpretag[6]= 56956.7;   wouttag[6]= 7488.16;  
            fbb[6]= 0.0650858;   fcc[6]= 0.120937;   fc[6]= 0.153119;   fll[6]= 0.660859;  
            fmcbb[6]= 0.0559094;   fmccc[6]= 0.103886;   fmcc[6]= 0.15996;   fmcll[6]= 0.680245;  
           Fbb[7]= 1.15493;   Fcc[7]= 1.15493;   Fc[7]= 0.949668;   Fll[7]= 0.963824; 
              cafac[7]= 0.766629;  
              winpretag[7]= 16564.5 ;   wintag[7]= 2676.39 ;   woutpretag[7]= 12698.8;   wouttag[7]= 2181.73;  
            fbb[7]= 0.0880046;   fcc[7]= 0.143325;   fc[7]= 0.141015;   fll[7]= 0.627655;  
            fmcbb[7]= 0.0761991;   fmccc[7]= 0.124099;   fmcc[7]= 0.148488;   fmcll[7]= 0.651214;  
      }  
      else if ( idsys==92   || sysname=="pdfa_down" ) { 
           Fbb[1]= 1.55418;   Fcc[1]= 1.55418;   Fc[1]= 0.568017;   Fll[1]= 1.04251; 
              cafac[1]= 1.05324;  
              winpretag[1]= 1.06218e+06 ;   wintag[1]= 39315.8 ;   woutpretag[1]= 1.11873e+06;   wouttag[1]= 34724.9;  
            fbb[1]= 0.0186154;   fcc[1]= 0.0544066;   fc[1]= 0.0796711;   fll[1]= 0.847307;  
            fmcbb[1]= 0.0119776;   fmccc[1]= 0.0350065;   fmcc[1]= 0.140262;   fmcll[1]= 0.812754;  
           Fbb[2]= 1.52662;   Fcc[2]= 1.52662;   Fc[2]= 0.557942;   Fll[2]= 1.02402; 
              cafac[2]= 0.943655;  
              winpretag[2]= 252186 ;   wintag[2]= 19507.4 ;   woutpretag[2]= 237977;   wouttag[2]= 18033;  
            fbb[2]= 0.0459539;   fcc[2]= 0.110503;   fc[2]= 0.0904192;   fll[2]= 0.753124;  
            fmcbb[2]= 0.0301018;   fmccc[2]= 0.0723841;   fmcc[2]= 0.162058;   fmcll[2]= 0.735456;  
           Fbb[3]= 1.49255;   Fcc[3]= 1.49255;   Fc[3]= 0.545491;   Fll[3]= 1.00117; 
              cafac[3]= 0.880675;  
              winpretag[3]= 56785.2 ;   wintag[3]= 6458.08 ;   woutpretag[3]= 50009.3;   wouttag[3]= 6024.59;  
            fbb[3]= 0.0749933;   fcc[3]= 0.146348;   fc[3]= 0.0886339;   fll[3]= 0.690025;  
            fmcbb[3]= 0.0502451;   fmccc[3]= 0.0980524;   fmcc[3]= 0.162485;   fmcll[3]= 0.689218;  
           Fbb[4]= 1.45512;   Fcc[4]= 1.45512;   Fc[4]= 0.531812;   Fll[4]= 0.976066; 
              cafac[4]= 0.909381;  
              winpretag[4]= 12391.3 ;   wintag[4]= 1874.29 ;   woutpretag[4]= 11268.4;   wouttag[4]= 1893.75;  
            fbb[4]= 0.104154;   fcc[4]= 0.17147;   fc[4]= 0.0799743;   fll[4]= 0.644402;  
            fmcbb[4]= 0.0715772;   fmccc[4]= 0.117839;   fmcc[4]= 0.150381;   fmcll[4]= 0.660203;  
           Fbb[5]= 1.40961;   Fcc[5]= 1.40961;   Fc[5]= 0.51518;   Fll[5]= 0.945539; 
              cafac[5]= 0.950944;  
              winpretag[5]= 3222.83 ;   wintag[5]= 613.922 ;   woutpretag[5]= 3064.73;   wouttag[5]= 668.556;  
            fbb[5]= 0.132591;   fcc[5]= 0.206186;   fc[5]= 0.0683199;   fll[5]= 0.592904;  
            fmcbb[5]= 0.0940618;   fmccc[5]= 0.146271;   fmcc[5]= 0.132614;   fmcll[5]= 0.627054;  
           Fbb[6]= 1.48214;   Fcc[6]= 1.48214;   Fc[6]= 0.541688;   Fll[6]= 0.99419; 
              cafac[6]= 0.8889;  
              winpretag[6]= 72399.3 ;   wintag[6]= 8946.29 ;   woutpretag[6]= 64355.7;   wouttag[6]= 8581.01;  
            fbb[6]= 0.0827726;   fcc[6]= 0.153528;   fc[6]= 0.0861735;   fll[6]= 0.677526;  
            fmcbb[6]= 0.0558466;   fmccc[6]= 0.103585;   fmcc[6]= 0.159083;   fmcll[6]= 0.681485;  
           Fbb[7]= 1.44549;   Fcc[7]= 1.44549;   Fc[7]= 0.528292;   Fll[7]= 0.969604; 
              cafac[7]= 0.918212;  
              winpretag[7]= 15614.1 ;   wintag[7]= 2488.21 ;   woutpretag[7]= 14337;   wouttag[7]= 2561.46;  
            fbb[7]= 0.110173;   fcc[7]= 0.178818;   fc[7]= 0.0775075;   fll[7]= 0.633502;  
            fmcbb[7]= 0.0762181;   fmccc[7]= 0.123708;   fmcc[7]= 0.146713;   fmcll[7]= 0.653361;  
      }  
      else if ( idsys==53   || sysname=="met_up" ) { 
           Fbb[1]= 1.37815;   Fcc[1]= 1.37815;   Fc[1]= 0.738125;   Fll[1]= 1.02337; 
              cafac[1]= 0.973709;  
              winpretag[1]= 1.08787e+06 ;   wintag[1]= 40321.4 ;   woutpretag[1]= 1.05927e+06;   wouttag[1]= 35619.6;  
            fbb[1]= 0.0164764;   fcc[1]= 0.0483172;   fc[1]= 0.10364;   fll[1]= 0.831567;  
            fmcbb[1]= 0.0119555;   fmccc[1]= 0.0350596;   fmcc[1]= 0.140409;   fmcll[1]= 0.812575;  
           Fbb[2]= 1.35988;   Fcc[2]= 1.35988;   Fc[2]= 0.728342;   Fll[2]= 1.00981; 
              cafac[2]= 0.888314;  
              winpretag[2]= 257677 ;   wintag[2]= 19987.4 ;   woutpretag[2]= 228898;   wouttag[2]= 17705.7;  
            fbb[2]= 0.0409686;   fcc[2]= 0.0985378;   fc[2]= 0.118311;   fll[2]= 0.742183;  
            fmcbb[2]= 0.0301266;   fmccc[2]= 0.0724607;   fmcc[2]= 0.162439;   fmcll[2]= 0.734974;  
           Fbb[3]= 1.33859;   Fcc[3]= 1.33859;   Fc[3]= 0.71694;   Fll[3]= 0.994001; 
              cafac[3]= 0.805574;  
              winpretag[3]= 58102.9 ;   wintag[3]= 6607.9 ;   woutpretag[3]= 46806.2;   wouttag[3]= 5582.9;  
            fbb[3]= 0.0671265;   fcc[3]= 0.131161;   fc[3]= 0.116563;   fll[3]= 0.68515;  
            fmcbb[3]= 0.0501471;   fmccc[3]= 0.097984;   fmcc[3]= 0.162583;   fmcll[3]= 0.689285;  
           Fbb[4]= 1.31516;   Fcc[4]= 1.31516;   Fc[4]= 0.704392;   Fll[4]= 0.976604; 
              cafac[4]= 0.82432;  
              winpretag[4]= 12738 ;   wintag[4]= 1939.13 ;   woutpretag[4]= 10500.2;   wouttag[4]= 1728.85;  
            fbb[4]= 0.0942493;   fcc[4]= 0.155625;   fc[4]= 0.10591;   fll[4]= 0.644216;  
            fmcbb[4]= 0.0716635;   fmccc[4]= 0.118332;   fmcc[4]= 0.150356;   fmcll[4]= 0.659649;  
           Fbb[5]= 1.2874;   Fcc[5]= 1.2874;   Fc[5]= 0.689522;   Fll[5]= 0.955987; 
              cafac[5]= 0.833277;  
              winpretag[5]= 3357.31 ;   wintag[5]= 647.422 ;   woutpretag[5]= 2797.57;   wouttag[5]= 594.628;  
            fbb[5]= 0.121872;   fcc[5]= 0.18805;   fc[5]= 0.0925594;   fll[5]= 0.597519;  
            fmcbb[5]= 0.094665;   fmccc[5]= 0.14607;   fmcc[5]= 0.134237;   fmcll[5]= 0.625028;  
           Fbb[6]= 1.33212;   Fcc[6]= 1.33212;   Fc[6]= 0.713475;   Fll[6]= 0.989196; 
              cafac[6]= 0.810048;  
              winpretag[6]= 74198.2 ;   wintag[6]= 9194.45 ;   woutpretag[6]= 60104.1;   wouttag[6]= 7905.72;  
            fbb[6]= 0.074406;   fcc[6]= 0.138078;   fc[6]= 0.113586;   fll[6]= 0.673929;  
            fmcbb[6]= 0.0558553;   fmccc[6]= 0.103653;   fmcc[6]= 0.159202;   fmcll[6]= 0.68129;  
           Fbb[7]= 1.30927;   Fcc[7]= 1.30927;   Fc[7]= 0.701238;   Fll[7]= 0.97223; 
              cafac[7]= 0.825986;  
              winpretag[7]= 16095.3 ;   wintag[7]= 2586.55 ;   woutpretag[7]= 13294.5;   wouttag[7]= 2324.03;  
            fbb[7]= 0.100109;   fcc[7]= 0.162504;   fc[7]= 0.103078;   fll[7]= 0.63431;  
            fmcbb[7]= 0.0764614;   fmccc[7]= 0.124117;   fmcc[7]= 0.146994;   fmcll[7]= 0.652427;  
      }  
      else if ( idsys==54   || sysname=="met_down" ) { 
           Fbb[1]= 1.39067;   Fcc[1]= 1.39067;   Fc[1]= 0.754379;   Fll[1]= 1.01967; 
              cafac[1]= 0.979207;  
              winpretag[1]= 1.07589e+06 ;   wintag[1]= 39746.1 ;   woutpretag[1]= 1.05352e+06;   wouttag[1]= 35743.6;  
            fbb[1]= 0.0167732;   fcc[1]= 0.0488204;   fc[1]= 0.105682;   fll[1]= 0.828724;  
            fmcbb[1]= 0.0120613;   fmccc[1]= 0.0351058;   fmcc[1]= 0.140092;   fmcll[1]= 0.812741;  
           Fbb[2]= 1.37058;   Fcc[2]= 1.37058;   Fc[2]= 0.743484;   Fll[2]= 1.00494; 
              cafac[2]= 0.888258;  
              winpretag[2]= 256327 ;   wintag[2]= 19847.9 ;   woutpretag[2]= 227684;   wouttag[2]= 17749.3;  
            fbb[2]= 0.0413779;   fcc[2]= 0.0992419;   fc[2]= 0.12072;   fll[2]= 0.73866;  
            fmcbb[2]= 0.03019;   fmccc[2]= 0.0724086;   fmcc[2]= 0.162371;   fmcll[2]= 0.73503;  
           Fbb[3]= 1.34816;   Fcc[3]= 1.34816;   Fc[3]= 0.731322;   Fll[3]= 0.9885; 
              cafac[3]= 0.809685;  
              winpretag[3]= 57828.1 ;   wintag[3]= 6572.98 ;   woutpretag[3]= 46822.6;   wouttag[3]= 5625.29;  
            fbb[3]= 0.0677455;   fcc[3]= 0.132308;   fc[3]= 0.119064;   fll[3]= 0.680883;  
            fmcbb[3]= 0.0502503;   fmccc[3]= 0.0981394;   fmcc[3]= 0.162807;   fmcll[3]= 0.688804;  
           Fbb[4]= 1.32439;   Fcc[4]= 1.32439;   Fc[4]= 0.718426;   Fll[4]= 0.971069; 
              cafac[4]= 0.823964;  
              winpretag[4]= 12724.7 ;   wintag[4]= 1941.39 ;   woutpretag[4]= 10484.7;   wouttag[4]= 1741.41;  
            fbb[4]= 0.094563;   fcc[4]= 0.157017;   fc[4]= 0.108585;   fll[4]= 0.639834;  
            fmcbb[4]= 0.0714012;   fmccc[4]= 0.118559;   fmcc[4]= 0.151143;   fmcll[4]= 0.658897;  
           Fbb[5]= 1.29637;   Fcc[5]= 1.29637;   Fc[5]= 0.703226;   Fll[5]= 0.950523; 
              cafac[5]= 0.845124;  
              winpretag[5]= 3338.79 ;   wintag[5]= 650.537 ;   woutpretag[5]= 2821.69;   wouttag[5]= 607.942;  
            fbb[5]= 0.121408;   fcc[5]= 0.190592;   fc[5]= 0.0959978;   fll[5]= 0.592001;  
            fmcbb[5]= 0.0936527;   fmccc[5]= 0.14702;   fmcc[5]= 0.136511;   fmcll[5]= 0.622816;  
           Fbb[6]= 1.34159;   Fcc[6]= 1.34159;   Fc[6]= 0.727759;   Fll[6]= 0.983684; 
              cafac[6]= 0.813772;  
              winpretag[6]= 73891.6 ;   wintag[6]= 9164.91 ;   woutpretag[6]= 60130.9;   wouttag[6]= 7974.16;  
            fbb[6]= 0.074933;   fcc[6]= 0.139344;   fc[6]= 0.116158;   fll[6]= 0.669566;  
            fmcbb[6]= 0.0558538;   fmccc[6]= 0.103864;   fmcc[6]= 0.15961;   fmcll[6]= 0.680672;  
           Fbb[7]= 1.31846;   Fcc[7]= 1.31846;   Fc[7]= 0.715213;   Fll[7]= 0.966726; 
              cafac[7]= 0.828439;  
              winpretag[7]= 16063.5 ;   wintag[7]= 2591.92 ;   woutpretag[7]= 13307.6;   wouttag[7]= 2349.14;  
            fbb[7]= 0.100238;   fcc[7]= 0.164115;   fc[7]= 0.105924;   fll[7]= 0.629723;  
            fmcbb[7]= 0.0760262;   fmccc[7]= 0.124474;   fmcc[7]= 0.148102;   fmcll[7]= 0.651398;  
      }  
      else if ( idsys==98   || sysname=="ifsr_nom" ) { 
           Fbb[1]= 1.3378;   Fcc[1]= 1.3378;   Fc[1]= 0.753543;   Fll[1]= 1.02289; 
              cafac[1]= 0.982541;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06391e+06;   wouttag[1]= 35920.5;  
            fbb[1]= 0.0161301;   fcc[1]= 0.0469301;   fc[1]= 0.105562;   fll[1]= 0.831378;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.3227;   Fcc[2]= 1.3227;   Fc[2]= 0.745038;   Fll[2]= 1.01134; 
              cafac[2]= 0.879316;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 226086;   wouttag[2]= 17306;  
            fbb[2]= 0.0398137;   fcc[2]= 0.0955209;   fc[2]= 0.12086;   fll[2]= 0.743805;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.30429;   Fcc[3]= 1.30429;   Fc[3]= 0.734671;   Fll[3]= 0.997271; 
              cafac[3]= 0.812183;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46994.5;   wouttag[3]= 5553.52;  
            fbb[3]= 0.0657024;   fcc[3]= 0.127966;   fc[3]= 0.119908;   fll[3]= 0.686423;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.28382;   Fcc[4]= 1.28382;   Fc[4]= 0.723141;   Fll[4]= 0.98162; 
              cafac[4]= 0.837737;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10687.2;   wouttag[4]= 1743.41;  
            fbb[4]= 0.0917645;   fcc[4]= 0.151657;   fc[4]= 0.108885;   fll[4]= 0.647693;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.26;   Fcc[5]= 1.26;   Fc[5]= 0.709723;   Fll[5]= 0.963405; 
              cafac[5]= 0.837415;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2821.97;   wouttag[5]= 590.001;  
            fbb[5]= 0.1177;   fcc[5]= 0.183255;   fc[5]= 0.0958164;   fll[5]= 0.603229;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.29864;   Fcc[6]= 1.29864;   Fc[6]= 0.731489;   Fll[6]= 0.992952; 
              cafac[6]= 0.817807;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60508.7;   wouttag[6]= 7883.19;  
            fbb[6]= 0.0726888;   fcc[6]= 0.134694;   fc[6]= 0.116854;   fll[6]= 0.675763;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.27877;   Fcc[7]= 1.27877;   Fc[7]= 0.720296;   Fll[7]= 0.977757; 
              cafac[7]= 0.837265;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13502.6;   wouttag[7]= 2334.26;  
            fbb[7]= 0.0972646;   fcc[7]= 0.158358;   fc[7]= 0.106114;   fll[7]= 0.638263;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==47   || sysname=="ifsr_up" ) { 
           Fbb[1]= 1.39644;   Fcc[1]= 1.39644;   Fc[1]= 0.716036;   Fll[1]= 1.02595; 
              cafac[1]= 0.97564;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.05643e+06;   wouttag[1]= 35176.8;  
            fbb[1]= 0.0168372;   fcc[1]= 0.0489874;   fc[1]= 0.100308;   fll[1]= 0.833868;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.37773;   Fcc[2]= 1.37773;   Fc[2]= 0.70644;   Fll[2]= 1.0122; 
              cafac[2]= 0.872506;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 224335;   wouttag[2]= 17184.7;  
            fbb[2]= 0.0414701;   fcc[2]= 0.0994949;   fc[2]= 0.114599;   fll[2]= 0.744436;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.35527;   Fcc[3]= 1.35527;   Fc[3]= 0.694923;   Fll[3]= 0.9957; 
              cafac[3]= 0.803791;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46508.9;   wouttag[3]= 5539.56;  
            fbb[3]= 0.0682701;   fcc[3]= 0.132968;   fc[3]= 0.11342;   fll[3]= 0.685342;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.33054;   Fcc[4]= 1.33054;   Fc[4]= 0.682242;   Fll[4]= 0.97753; 
              cafac[4]= 0.831046;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10601.8;   wouttag[4]= 1750.07;  
            fbb[4]= 0.0951033;   fcc[4]= 0.157175;   fc[4]= 0.102727;   fll[4]= 0.644995;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.30192;   Fcc[5]= 1.30192;   Fc[5]= 0.667568;   Fll[5]= 0.956505; 
              cafac[5]= 0.819861;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2762.81;   wouttag[5]= 585.781;  
            fbb[5]= 0.121615;   fcc[5]= 0.189351;   fc[5]= 0.0901252;   fll[5]= 0.598908;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.34843;   Fcc[6]= 1.34843;   Fc[6]= 0.691417;   Fll[6]= 0.990676; 
              cafac[6]= 0.809156;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 59868.6;   wouttag[6]= 7874.07;  
            fbb[6]= 0.0754753;   fcc[6]= 0.139857;   fc[6]= 0.110453;   fll[6]= 0.674214;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.32445;   Fcc[7]= 1.32445;   Fc[7]= 0.679123;   Fll[7]= 0.973061; 
              cafac[7]= 0.827963;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13352.6;   wouttag[7]= 2337.41;  
            fbb[7]= 0.100739;   fcc[7]= 0.164015;   fc[7]= 0.100048;   fll[7]= 0.635198;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==98   || sysname=="ifsr_nom" ) { 
           Fbb[1]= 1.41215;   Fcc[1]= 1.41215;   Fc[1]= 0.801244;   Fll[1]= 1.01035; 
              cafac[1]= 0.99132;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.07341e+06;   wouttag[1]= 37635;  
            fbb[1]= 0.0170266;   fcc[1]= 0.0495382;   fc[1]= 0.112244;   fll[1]= 0.821191;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.3878;   Fcc[2]= 1.3878;   Fc[2]= 0.787431;   Fll[2]= 0.992936; 
              cafac[2]= 0.888724;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 228505;   wouttag[2]= 18103.6;  
            fbb[2]= 0.0417732;   fcc[2]= 0.100222;   fc[2]= 0.127737;   fll[2]= 0.730267;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.36323;   Fcc[3]= 1.36323;   Fc[3]= 0.773488;   Fll[3]= 0.975354; 
              cafac[3]= 0.820221;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 47459.6;   wouttag[3]= 5793.17;  
            fbb[3]= 0.0686711;   fcc[3]= 0.133748;   fc[3]= 0.126243;   fll[3]= 0.671338;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.33846;   Fcc[4]= 1.33846;   Fc[4]= 0.759437;   Fll[4]= 0.957636; 
              cafac[4]= 0.831727;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10610.5;   wouttag[4]= 1783.4;  
            fbb[4]= 0.0956699;   fcc[4]= 0.158112;   fc[4]= 0.11435;   fll[4]= 0.631868;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.30986;   Fcc[5]= 1.30986;   Fc[5]= 0.743206;   Fll[5]= 0.937169; 
              cafac[5]= 0.832447;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2805.22;   wouttag[5]= 603.259;  
            fbb[5]= 0.122357;   fcc[5]= 0.190506;   fc[5]= 0.100337;   fll[5]= 0.586801;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.35638;   Fcc[6]= 1.35638;   Fc[6]= 0.769604;   Fll[6]= 0.970457; 
              cafac[6]= 0.822523;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60857.7;   wouttag[6]= 8184.79;  
            fbb[6]= 0.0759205;   fcc[6]= 0.140682;   fc[6]= 0.122943;   fll[6]= 0.660454;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.33238;   Fcc[7]= 1.33238;   Fc[7]= 0.755987;   Fll[7]= 0.953286; 
              cafac[7]= 0.831431;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13408.5;   wouttag[7]= 2387.55;  
            fbb[7]= 0.101342;   fcc[7]= 0.164997;   fc[7]= 0.111372;   fll[7]= 0.622289;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==48   || sysname=="ifsr_down" ) { 
           Fbb[1]= 1.47126;   Fcc[1]= 1.47126;   Fc[1]= 0.764385;   Fll[1]= 1.01328; 
              cafac[1]= 0.984412;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06593e+06;   wouttag[1]= 36892.8;  
            fbb[1]= 0.0177393;   fcc[1]= 0.0516119;   fc[1]= 0.107081;   fll[1]= 0.823568;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.44275;   Fcc[2]= 1.44275;   Fc[2]= 0.749572;   Fll[2]= 0.993642; 
              cafac[2]= 0.881913;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 226754;   wouttag[2]= 17982.3;  
            fbb[2]= 0.0434271;   fcc[2]= 0.10419;   fc[2]= 0.121596;   fll[2]= 0.730787;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.41377;   Fcc[3]= 1.41377;   Fc[3]= 0.73452;   Fll[3]= 0.973689; 
              cafac[3]= 0.81182;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46973.5;   wouttag[3]= 5778.03;  
            fbb[3]= 0.0712174;   fcc[3]= 0.138708;   fc[3]= 0.119883;   fll[3]= 0.670192;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.38453;   Fcc[4]= 1.38453;   Fc[4]= 0.719328;   Fll[4]= 0.95355; 
              cafac[4]= 0.825167;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10526.8;   wouttag[4]= 1789.56;  
            fbb[4]= 0.0989629;   fcc[4]= 0.163554;   fc[4]= 0.108311;   fll[4]= 0.629172;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.35094;   Fcc[5]= 1.35094;   Fc[5]= 0.701872;   Fll[5]= 0.930411; 
              cafac[5]= 0.815014;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2746.48;   wouttag[5]= 598.563;  
            fbb[5]= 0.126194;   fcc[5]= 0.19648;   fc[5]= 0.0947564;   fll[5]= 0.582569;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.40568;   Fcc[6]= 1.40568;   Fc[6]= 0.730313;   Fll[6]= 0.968113; 
              cafac[6]= 0.813894;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60219.2;   wouttag[6]= 8173.65;  
            fbb[6]= 0.0786797;   fcc[6]= 0.145795;   fc[6]= 0.116666;   fll[6]= 0.658859;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.37738;   Fcc[7]= 1.37738;   Fc[7]= 0.715609;   Fll[7]= 0.94862; 
              cafac[7]= 0.822251;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13260.5;   wouttag[7]= 2389.72;  
            fbb[7]= 0.104764;   fcc[7]= 0.170569;   fc[7]= 0.105423;   fll[7]= 0.619244;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==72   || sysname=="iqopt3_one" ) { 
           Fbb[1]= 1.38386;   Fcc[1]= 1.38386;   Fc[1]= 0.73248;   Fll[1]= 1.02397; 
              cafac[1]= 0.98123;  
              winpretag[1]= 1.0761e+06 ;   wintag[1]= 39813.6 ;   woutpretag[1]= 1.0559e+06;   wouttag[1]= 35373;  
            fbb[1]= 0.0165743;   fcc[1]= 0.0483832;   fc[1]= 0.102679;   fll[1]= 0.832364;  
            fmcbb[1]= 0.0119769;   fmccc[1]= 0.0349626;   fmcc[1]= 0.14018;   fmcll[1]= 0.812881;  
           Fbb[2]= 1.36523;   Fcc[2]= 1.36523;   Fc[2]= 0.722622;   Fll[2]= 1.01019; 
              cafac[2]= 0.883748;  
              winpretag[2]= 257295 ;   wintag[2]= 19954.8 ;   woutpretag[2]= 227384;   wouttag[2]= 17583.5;  
            fbb[2]= 0.0412049;   fcc[2]= 0.0990061;   fc[2]= 0.117229;   fll[2]= 0.74256;  
            fmcbb[2]= 0.0301816;   fmccc[2]= 0.0725196;   fmcc[2]= 0.162228;   fmcll[2]= 0.735071;  
           Fbb[3]= 1.34322;   Fcc[3]= 1.34322;   Fc[3]= 0.710971;   Fll[3]= 0.9939; 
              cafac[3]= 0.795176;  
              winpretag[3]= 58394.8 ;   wintag[3]= 6661.81 ;   woutpretag[3]= 46434.1;   wouttag[3]= 5559.73;  
            fbb[3]= 0.0677397;   fcc[3]= 0.132221;   fc[3]= 0.115347;   fll[3]= 0.684692;  
            fmcbb[3]= 0.0504308;   fmccc[3]= 0.098436;   fmcc[3]= 0.162238;   fmcll[3]= 0.688895;  
           Fbb[4]= 1.31947;   Fcc[4]= 1.31947;   Fc[4]= 0.698397;   Fll[4]= 0.976323; 
              cafac[4]= 0.818224;  
              winpretag[4]= 12742.8 ;   wintag[4]= 1942.21 ;   woutpretag[4]= 10426.5;   wouttag[4]= 1722.19;  
            fbb[4]= 0.0949527;   fcc[4]= 0.156103;   fc[4]= 0.104569;   fll[4]= 0.644375;  
            fmcbb[4]= 0.071963;   fmccc[4]= 0.118308;   fmcc[4]= 0.149727;   fmcll[4]= 0.660002;  
           Fbb[5]= 1.29081;   Fcc[5]= 1.29081;   Fc[5]= 0.68323;   Fll[5]= 0.955119; 
              cafac[5]= 0.741412;  
              winpretag[5]= 3704.94 ;   wintag[5]= 713.409 ;   woutpretag[5]= 2746.89;   wouttag[5]= 583.264;  
            fbb[5]= 0.120932;   fcc[5]= 0.190669;   fc[5]= 0.0908522;   fll[5]= 0.597547;  
            fmcbb[5]= 0.0936871;   fmccc[5]= 0.147712;   fmcc[5]= 0.132975;   fmcll[5]= 0.625626;  
           Fbb[6]= 1.33644;   Fcc[6]= 1.33644;   Fc[6]= 0.707381;   Fll[6]= 0.988881; 
              cafac[6]= 0.795629;  
              winpretag[6]= 74842.5 ;   wintag[6]= 9317.43 ;   woutpretag[6]= 59546.8;   wouttag[6]= 7879.02;  
            fbb[6]= 0.0751589;   fcc[6]= 0.139335;   fc[6]= 0.112233;   fll[6]= 0.673273;  
            fmcbb[6]= 0.0562382;   fmccc[6]= 0.104259;   fmcc[6]= 0.158659;   fmcll[6]= 0.680844;  
           Fbb[7]= 1.3129;   Fcc[7]= 1.3129;   Fc[7]= 0.694922;   Fll[7]= 0.971465; 
              cafac[7]= 0.798544;  
              winpretag[7]= 16447.7 ;   wintag[7]= 2655.62 ;   woutpretag[7]= 13134.2;   wouttag[7]= 2311.36;  
            fbb[7]= 0.100905;   fcc[7]= 0.164023;   fc[7]= 0.101426;   fll[7]= 0.633646;  
            fmcbb[7]= 0.0768565;   fmccc[7]= 0.124931;   fmcc[7]= 0.145954;   fmcll[7]= 0.652258;  
      }  
      else if ( idsys==71   || sysname=="ptjmin_one" ) { 
           Fbb[1]= 1.38055;   Fcc[1]= 1.38055;   Fc[1]= 0.732292;   Fll[1]= 1.02415; 
              cafac[1]= 0.975324;  
              winpretag[1]= 1.08193e+06 ;   wintag[1]= 40096.9 ;   woutpretag[1]= 1.05523e+06;   wouttag[1]= 35385.4;  
            fbb[1]= 0.016576;   fcc[1]= 0.0483768;   fc[1]= 0.102658;   fll[1]= 0.832389;  
            fmcbb[1]= 0.0120068;   fmccc[1]= 0.0350416;   fmcc[1]= 0.140187;   fmcll[1]= 0.812765;  
           Fbb[2]= 1.36216;   Fcc[2]= 1.36216;   Fc[2]= 0.722533;   Fll[2]= 1.0105; 
              cafac[2]= 0.88404;  
              winpretag[2]= 257025 ;   wintag[2]= 19964.9 ;   woutpretag[2]= 227220;   wouttag[2]= 17583.5;  
            fbb[2]= 0.0412267;   fcc[2]= 0.0990095;   fc[2]= 0.11718;   fll[2]= 0.742584;  
            fmcbb[2]= 0.0302658;   fmccc[2]= 0.0726859;   fmcc[2]= 0.162179;   fmcll[2]= 0.73487;  
           Fbb[3]= 1.34077;   Fcc[3]= 1.34077;   Fc[3]= 0.711191;   Fll[3]= 0.994635; 
              cafac[3]= 0.803258;  
              winpretag[3]= 57904.8 ;   wintag[3]= 6602.08 ;   woutpretag[3]= 46512.5;   wouttag[3]= 5559.83;  
            fbb[3]= 0.0675103;   fcc[3]= 0.131765;   fc[3]= 0.115621;   fll[3]= 0.685103;  
            fmcbb[3]= 0.0503518;   fmccc[3]= 0.0982754;   fmcc[3]= 0.162574;   fmcll[3]= 0.688799;  
           Fbb[4]= 1.31722;   Fcc[4]= 1.31722;   Fc[4]= 0.698697;   Fll[4]= 0.977162; 
              cafac[4]= 0.820802;  
              winpretag[4]= 12713.7 ;   wintag[4]= 1937.63 ;   woutpretag[4]= 10435.5;   wouttag[4]= 1721.11;  
            fbb[4]= 0.094628;   fcc[4]= 0.155657;   fc[4]= 0.104821;   fll[4]= 0.644894;  
            fmcbb[4]= 0.0718392;   fmccc[4]= 0.118171;   fmcc[4]= 0.150023;   fmcll[4]= 0.659967;  
           Fbb[5]= 1.29022;   Fcc[5]= 1.29022;   Fc[5]= 0.684377;   Fll[5]= 0.957135; 
              cafac[5]= 0.827874;  
              winpretag[5]= 3311.49 ;   wintag[5]= 638.742 ;   woutpretag[5]= 2741.5;   wouttag[5]= 582.546;  
            fbb[5]= 0.12111;   fcc[5]= 0.188173;   fc[5]= 0.0927873;   fll[5]= 0.597929;  
            fmcbb[5]= 0.0938679;   fmccc[5]= 0.145846;   fmcc[5]= 0.135579;   fmcll[5]= 0.624707;  
           Fbb[6]= 1.33433;   Fcc[6]= 1.33433;   Fc[6]= 0.707773;   Fll[6]= 0.989854; 
              cafac[6]= 0.807341;  
              winpretag[6]= 73930 ;   wintag[6]= 9178.45 ;   woutpretag[6]= 59686.7;   wouttag[6]= 7863.63;  
            fbb[6]= 0.0747173;   fcc[6]= 0.13854;   fc[6]= 0.112682;   fll[6]= 0.674061;  
            fmcbb[6]= 0.0559962;   fmccc[6]= 0.103828;   fmcc[6]= 0.159207;   fmcll[6]= 0.68097;  
           Fbb[7]= 1.31155;   Fcc[7]= 1.31155;   Fc[7]= 0.695689;   Fll[7]= 0.972955; 
              cafac[7]= 0.82202;  
              winpretag[7]= 16025.2 ;   wintag[7]= 2576.37 ;   woutpretag[7]= 13173.1;   wouttag[7]= 2304.23;  
            fbb[7]= 0.100191;   fcc[7]= 0.162487;   fc[7]= 0.102293;   fll[7]= 0.635029;  
            fmcbb[7]= 0.0763913;   fmccc[7]= 0.123889;   fmcc[7]= 0.147039;   fmcll[7]= 0.652681;  
      }  
      else if ( idsys==69   || sysname=="powhe_one" ) { 
           Fbb[1]= 1.43718;   Fcc[1]= 1.43718;   Fc[1]= 0.758781;   Fll[1]= 1.01622; 
              cafac[1]= 0.983348;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06478e+06;   wouttag[1]= 36562.6;  
            fbb[1]= 0.0173284;   fcc[1]= 0.0504164;   fc[1]= 0.106296;   fll[1]= 0.825959;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.41242;   Fcc[2]= 1.41242;   Fc[2]= 0.745709;   Fll[2]= 0.998714; 
              cafac[2]= 0.880666;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 226433;   wouttag[2]= 17780.9;  
            fbb[2]= 0.0425142;   fcc[2]= 0.102;   fc[2]= 0.120969;   fll[2]= 0.734517;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.38629;   Fcc[3]= 1.38629;   Fc[3]= 0.731913;   Fll[3]= 0.980237; 
              cafac[3]= 0.811081;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 46930.7;   wouttag[3]= 5711.71;  
            fbb[3]= 0.0698327;   fcc[3]= 0.136011;   fc[3]= 0.119457;   fll[3]= 0.674699;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.35932;   Fcc[4]= 1.35932;   Fc[4]= 0.717677;   Fll[4]= 0.961171; 
              cafac[4]= 0.83117;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 10603.4;   wouttag[4]= 1783.47;  
            fbb[4]= 0.097161;   fcc[4]= 0.160576;   fc[4]= 0.108063;   fll[4]= 0.634201;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.32825;   Fcc[5]= 1.32825;   Fc[5]= 0.701269;   Fll[5]= 0.939196; 
              cafac[5]= 0.820165;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2763.84;   wouttag[5]= 596.025;  
            fbb[5]= 0.124075;   fcc[5]= 0.19318;   fc[5]= 0.094675;   fll[5]= 0.58807;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.37883;   Fcc[6]= 1.37883;   Fc[6]= 0.727974;   Fll[6]= 0.974962; 
              cafac[6]= 0.814732;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 60281.2;   wouttag[6]= 8094.45;  
            fbb[6]= 0.0771768;   fcc[6]= 0.14301;   fc[6]= 0.116293;   fll[6]= 0.66352;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.35271;   Fcc[7]= 1.35271;   Fc[7]= 0.714185;   Fll[7]= 0.956495; 
              cafac[7]= 0.828088;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13354.6;   wouttag[7]= 2381.08;  
            fbb[7]= 0.102888;   fcc[7]= 0.167514;   fc[7]= 0.105214;   fll[7]= 0.624384;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==70   || sysname=="powpy_one" ) { 
           Fbb[1]= 1.3452;   Fcc[1]= 1.3452;   Fc[1]= 0.783947;   Fll[1]= 1.01722; 
              cafac[1]= 0.987303;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40113.3 ;   woutpretag[1]= 1.06906e+06;   wouttag[1]= 36776.5;  
            fbb[1]= 0.0162194;   fcc[1]= 0.0471898;   fc[1]= 0.109821;   fll[1]= 0.82677;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.32802;   Fcc[2]= 1.32802;   Fc[2]= 0.773937;   Fll[2]= 1.00423; 
              cafac[2]= 0.88519;  
              winpretag[2]= 257116 ;   wintag[2]= 19814 ;   woutpretag[2]= 227596;   wouttag[2]= 17643.2;  
            fbb[2]= 0.0399739;   fcc[2]= 0.0959054;   fc[2]= 0.125548;   fll[2]= 0.738573;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.30876;   Fcc[3]= 1.30876;   Fc[3]= 0.762709;   Fll[3]= 0.98966; 
              cafac[3]= 0.825836;  
              winpretag[3]= 57862 ;   wintag[3]= 6569.1 ;   woutpretag[3]= 47784.5;   wouttag[3]= 5700.74;  
            fbb[3]= 0.0659272;   fcc[3]= 0.128404;   fc[3]= 0.124484;   fll[3]= 0.681185;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.28816;   Fcc[4]= 1.28816;   Fc[4]= 0.750705;   Fll[4]= 0.974084; 
              cafac[4]= 0.864373;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1942.14 ;   woutpretag[4]= 11027;   wouttag[4]= 1811.83;  
            fbb[4]= 0.0920743;   fcc[4]= 0.152169;   fc[4]= 0.113036;   fll[4]= 0.642721;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.26421;   Fcc[5]= 1.26421;   Fc[5]= 0.736748;   Fll[5]= 0.955974; 
              cafac[5]= 0.874888;  
              winpretag[5]= 3369.85 ;   wintag[5]= 645.761 ;   woutpretag[5]= 2948.24;   wouttag[5]= 619.992;  
            fbb[5]= 0.118093;   fcc[5]= 0.183867;   fc[5]= 0.0994649;   fll[5]= 0.598575;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.30307;   Fcc[6]= 1.30307;   Fc[6]= 0.759396;   Fll[6]= 0.985362; 
              cafac[6]= 0.835111;  
              winpretag[6]= 73989 ;   wintag[6]= 9157 ;   woutpretag[6]= 61789;   wouttag[6]= 8120.49;  
            fbb[6]= 0.0729366;   fcc[6]= 0.135153;   fc[6]= 0.121312;   fll[6]= 0.670598;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.28308;   Fcc[7]= 1.28308;   Fc[7]= 0.747745;   Fll[7]= 0.970243; 
              cafac[7]= 0.866423;  
              winpretag[7]= 16127 ;   wintag[7]= 2587.9 ;   woutpretag[7]= 13972.8;   wouttag[7]= 2432.13;  
            fbb[7]= 0.0975922;   fcc[7]= 0.158892;   fc[7]= 0.110158;   fll[7]= 0.633359;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==15   || sysname=="btagb_up" ) { 
           Fbb[1]= 1.33533;   Fcc[1]= 1.33533;   Fc[1]= 0.704012;   Fll[1]= 1.03157; 
              cafac[1]= 0.974093;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 40355.5 ;   woutpretag[1]= 1.05476e+06;   wouttag[1]= 34857.3;  
            fbb[1]= 0.0161003;   fcc[1]= 0.0468434;   fc[1]= 0.0986233;   fll[1]= 0.838433;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.32274;   Fcc[2]= 1.32274;   Fc[2]= 0.697379;   Fll[2]= 1.02185; 
              cafac[2]= 0.871493;  
              winpretag[2]= 257116 ;   wintag[2]= 20036.7 ;   woutpretag[2]= 224074;   wouttag[2]= 17086;  
            fbb[2]= 0.039815;   fcc[2]= 0.0955242;   fc[2]= 0.113129;   fll[2]= 0.751532;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.30504;   Fcc[3]= 1.30504;   Fc[3]= 0.688042;   Fll[3]= 1.00817; 
              cafac[3]= 0.801973;  
              winpretag[3]= 57862 ;   wintag[3]= 6670.06 ;   woutpretag[3]= 46403.7;   wouttag[3]= 5517.8;  
            fbb[3]= 0.0657397;   fcc[3]= 0.128039;   fc[3]= 0.112297;   fll[3]= 0.693924;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.28416;   Fcc[4]= 1.28416;   Fc[4]= 0.677038;   Fll[4]= 0.992044; 
              cafac[4]= 0.823288;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1975.24 ;   woutpretag[4]= 10502.8;   wouttag[4]= 1732.94;  
            fbb[4]= 0.0917886;   fcc[4]= 0.151697;   fc[4]= 0.101943;   fll[4]= 0.654571;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.25986;   Fcc[5]= 1.25986;   Fc[5]= 0.664225;   Fll[5]= 0.97327; 
              cafac[5]= 0.815935;  
              winpretag[5]= 3369.85 ;   wintag[5]= 660.795 ;   woutpretag[5]= 2749.58;   wouttag[5]= 586.432;  
            fbb[5]= 0.117687;   fcc[5]= 0.183234;   fc[5]= 0.0896739;   fll[5]= 0.609405;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.29927;   Fcc[6]= 1.29927;   Fc[6]= 0.685004;   Fll[6]= 1.00372; 
              cafac[6]= 0.806208;  
              winpretag[6]= 73989 ;   wintag[6]= 9306.09 ;   woutpretag[6]= 59650.5;   wouttag[6]= 7836.97;  
            fbb[6]= 0.0727239;   fcc[6]= 0.134759;   fc[6]= 0.109428;   fll[6]= 0.683089;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.27901;   Fcc[7]= 1.27901;   Fc[7]= 0.67432;   Fll[7]= 0.988061; 
              cafac[7]= 0.821172;  
              winpretag[7]= 16127 ;   wintag[7]= 2636.04 ;   woutpretag[7]= 13243.1;   wouttag[7]= 2320.71;  
            fbb[7]= 0.0972824;   fcc[7]= 0.158387;   fc[7]= 0.0993406;   fll[7]= 0.64499;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==16   || sysname=="btagb_down" ) { 
           Fbb[1]= 1.54246;   Fcc[1]= 1.54246;   Fc[1]= 0.77812;   Fll[1]= 1.00678; 
              cafac[1]= 0.986246;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 39871.1 ;   woutpretag[1]= 1.06792e+06;   wouttag[1]= 37282.1;  
            fbb[1]= 0.0185977;   fcc[1]= 0.0541095;   fc[1]= 0.109005;   fll[1]= 0.818288;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.50557;   Fcc[2]= 1.50557;   Fc[2]= 0.759514;   Fll[2]= 0.982709; 
              cafac[2]= 0.883107;  
              winpretag[2]= 257116 ;   wintag[2]= 19589.7 ;   woutpretag[2]= 227061;   wouttag[2]= 18097.3;  
            fbb[2]= 0.0453183;   fcc[2]= 0.108728;   fc[2]= 0.123209;   fll[2]= 0.722746;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.4704;   Fcc[3]= 1.4704;   Fc[3]= 0.741772;   Fll[3]= 0.959753; 
              cafac[3]= 0.813818;  
              winpretag[3]= 57862 ;   wintag[3]= 6466.42 ;   woutpretag[3]= 47089.1;   wouttag[3]= 5800.76;  
            fbb[3]= 0.07407;   fcc[3]= 0.144264;   fc[3]= 0.121066;   fll[3]= 0.6606;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.43629;   Fcc[4]= 1.43629;   Fc[4]= 0.72456;   Fll[4]= 0.937484; 
              cafac[4]= 0.833058;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1907.92 ;   woutpretag[4]= 10627.5;   wouttag[4]= 1806.06;  
            fbb[4]= 0.102662;   fcc[4]= 0.169667;   fc[4]= 0.109099;   fll[4]= 0.618571;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.39733;   Fcc[5]= 1.39733;   Fc[5]= 0.704909;   Fll[5]= 0.912057; 
              cafac[5]= 0.818782;  
              winpretag[5]= 3369.85 ;   wintag[5]= 630.274 ;   woutpretag[5]= 2759.18;   wouttag[5]= 596.906;  
            fbb[5]= 0.130528;   fcc[5]= 0.203228;   fc[5]= 0.0951664;   fll[5]= 0.571077;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.46094;   Fcc[6]= 1.46094;   Fc[6]= 0.736998;   Fll[6]= 0.953576; 
              cafac[6]= 0.816982;  
              winpretag[6]= 73989 ;   wintag[6]= 9004.61 ;   woutpretag[6]= 60447.7;   wouttag[6]= 8210.15;  
            fbb[6]= 0.0817729;   fcc[6]= 0.151527;   fc[6]= 0.117734;   fll[6]= 0.648966;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.42797;   Fcc[7]= 1.42797;   Fc[7]= 0.720364;   Fll[7]= 0.932055; 
              cafac[7]= 0.829087;  
              winpretag[7]= 16127 ;   wintag[7]= 2538.2 ;   woutpretag[7]= 13370.7;   wouttag[7]= 2404.71;  
            fbb[7]= 0.108613;   fcc[7]= 0.176834;   fc[7]= 0.106124;   fll[7]= 0.60843;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==85   || sysname=="btagc_up" ) { 
           Fbb[1]= 1.33771;   Fcc[1]= 1.33771;   Fc[1]= 0.639227;   Fll[1]= 1.0426; 
              cafac[1]= 0.96267;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 44385.4 ;   woutpretag[1]= 1.04239e+06;   wouttag[1]= 35806.8;  
            fbb[1]= 0.0161291;   fcc[1]= 0.046927;   fc[1]= 0.0895477;   fll[1]= 0.847396;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.32794;   Fcc[2]= 1.32794;   Fc[2]= 0.634558;   Fll[2]= 1.03498; 
              cafac[2]= 0.859483;  
              winpretag[2]= 257116 ;   wintag[2]= 21595.5 ;   woutpretag[2]= 220986;   wouttag[2]= 17457.8;  
            fbb[2]= 0.0399715;   fcc[2]= 0.0958995;   fc[2]= 0.102938;   fll[2]= 0.761191;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.31073;   Fcc[3]= 1.31073;   Fc[3]= 0.626336;   Fll[3]= 1.02157; 
              cafac[3]= 0.790641;  
              winpretag[3]= 57862 ;   wintag[3]= 7095.09 ;   woutpretag[3]= 45748;   wouttag[3]= 5618.55;  
            fbb[3]= 0.0660268;   fcc[3]= 0.128598;   fc[3]= 0.102226;   fll[3]= 0.703149;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.28897;   Fcc[4]= 1.28897;   Fc[4]= 0.615935;   Fll[4]= 1.00461; 
              cafac[4]= 0.812654;  
              winpretag[4]= 12757.2 ;   wintag[4]= 2072.96 ;   woutpretag[4]= 10367.2;   wouttag[4]= 1757.28;  
            fbb[4]= 0.0921321;   fcc[4]= 0.152265;   fc[4]= 0.0927429;   fll[4]= 0.66286;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.26363;   Fcc[5]= 1.26363;   Fc[5]= 0.603826;   Fll[5]= 0.984856; 
              cafac[5]= 0.805898;  
              winpretag[5]= 3369.85 ;   wintag[5]= 688.563 ;   woutpretag[5]= 2715.76;   wouttag[5]= 592.679;  
            fbb[5]= 0.118038;   fcc[5]= 0.183782;   fc[5]= 0.0815197;   fll[5]= 0.61666;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.30472;   Fcc[6]= 1.30472;   Fc[6]= 0.623462;   Fll[6]= 1.01688; 
              cafac[6]= 0.795064;  
              winpretag[6]= 73989 ;   wintag[6]= 9856.62 ;   woutpretag[6]= 58826;   wouttag[6]= 7968.03;  
            fbb[6]= 0.0730288;   fcc[6]= 0.135324;   fc[6]= 0.0995971;   fll[6]= 0.69205;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.28359;   Fcc[7]= 1.28359;   Fc[7]= 0.613365;   Fll[7]= 1.00041; 
              cafac[7]= 0.81067;  
              winpretag[7]= 16127 ;   wintag[7]= 2761.52 ;   woutpretag[7]= 13073.7;   wouttag[7]= 2351.28;  
            fbb[7]= 0.0976309;   fcc[7]= 0.158955;   fc[7]= 0.0903607;   fll[7]= 0.653054;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==86   || sysname=="btagc_down" ) { 
           Fbb[1]= 1.53665;   Fcc[1]= 1.53665;   Fc[1]= 0.872676;   Fll[1]= 0.990822; 
              cafac[1]= 1.00386;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 35841.2 ;   woutpretag[1]= 1.08699e+06;   wouttag[1]= 36245;  
            fbb[1]= 0.0185278;   fcc[1]= 0.053906;   fc[1]= 0.122251;   fll[1]= 0.805315;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.49552;   Fcc[2]= 1.49552;   Fc[2]= 0.849317;   Fll[2]= 0.964299; 
              cafac[2]= 0.901361;  
              winpretag[2]= 257116 ;   wintag[2]= 18026.2 ;   woutpretag[2]= 231754;   wouttag[2]= 17709.6;  
            fbb[2]= 0.0450157;   fcc[2]= 0.108002;   fc[2]= 0.137776;   fll[2]= 0.709206;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.45988;   Fcc[3]= 1.45988;   Fc[3]= 0.829077;   Fll[3]= 0.94132; 
              cafac[3]= 0.830891;  
              winpretag[3]= 57862 ;   wintag[3]= 6039.83 ;   woutpretag[3]= 48077;   wouttag[3]= 5698.6;  
            fbb[3]= 0.0735401;   fcc[3]= 0.143232;   fc[3]= 0.135316;   fll[3]= 0.647912;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.42742;   Fcc[4]= 1.42742;   Fc[4]= 0.810641;   Fll[4]= 0.920388; 
              cafac[4]= 0.848874;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1809.67 ;   woutpretag[4]= 10829.2;   wouttag[4]= 1782.42;  
            fbb[4]= 0.102028;   fcc[4]= 0.16862;   fc[4]= 0.12206;   fll[4]= 0.607291;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.39033;   Fcc[5]= 1.39033;   Fc[5]= 0.789577;   Fll[5]= 0.896472; 
              cafac[5]= 0.833434;  
              winpretag[5]= 3369.85 ;   wintag[5]= 602.02 ;   woutpretag[5]= 2808.55;   wouttag[5]= 590.454;  
            fbb[5]= 0.129874;   fcc[5]= 0.20221;   fc[5]= 0.106597;   fll[5]= 0.561319;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.45089;   Fcc[6]= 1.45089;   Fc[6]= 0.823969;   Fll[6]= 0.93552; 
              cafac[6]= 0.833702;  
              winpretag[6]= 73989 ;   wintag[6]= 8451.52 ;   woutpretag[6]= 61684.8;   wouttag[6]= 8078.35;  
            fbb[6]= 0.0812103;   fcc[6]= 0.150484;   fc[6]= 0.131628;   fll[6]= 0.636678;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.41951;   Fcc[7]= 1.41951;   Fc[7]= 0.806147;   Fll[7]= 0.915286; 
              cafac[7]= 0.844636;  
              winpretag[7]= 16127 ;   wintag[7]= 2411.69 ;   woutpretag[7]= 13621.5;   wouttag[7]= 2374.63;  
            fbb[7]= 0.107969;   fcc[7]= 0.175786;   fc[7]= 0.118761;   fll[7]= 0.597483;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==17   || sysname=="bmtag_up" ) { 
           Fbb[1]= 1.21879;   Fcc[1]= 1.21879;   Fc[1]= 0.708452;   Fll[1]= 1.03756; 
              cafac[1]= 0.975563;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 42541 ;   woutpretag[1]= 1.05635e+06;   wouttag[1]= 36560.1;  
            fbb[1]= 0.0146952;   fcc[1]= 0.0427554;   fc[1]= 0.0992453;   fll[1]= 0.843304;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.21549;   Fcc[2]= 1.21549;   Fc[2]= 0.706532;   Fll[2]= 1.03475; 
              cafac[2]= 0.873603;  
              winpretag[2]= 257116 ;   wintag[2]= 21020.6 ;   woutpretag[2]= 224617;   wouttag[2]= 17484.1;  
            fbb[2]= 0.0365867;   fcc[2]= 0.0877787;   fc[2]= 0.114614;   fll[2]= 0.761021;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.20582;   Fcc[3]= 1.20582;   Fc[3]= 0.700912;   Fll[3]= 1.02652; 
              cafac[3]= 0.803416;  
              winpretag[3]= 57862 ;   wintag[3]= 6967.84 ;   woutpretag[3]= 46487.2;   wouttag[3]= 5561.86;  
            fbb[3]= 0.060742;   fcc[3]= 0.118305;   fc[3]= 0.114398;   fll[3]= 0.706555;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.19213;   Fcc[4]= 1.19213;   Fc[4]= 0.692951;   Fll[4]= 1.01486; 
              cafac[4]= 0.825355;  
              winpretag[4]= 12757.2 ;   wintag[4]= 2067.94 ;   woutpretag[4]= 10529.2;   wouttag[4]= 1744.65;  
            fbb[4]= 0.0852101;   fcc[4]= 0.140825;   fc[4]= 0.104339;   fll[4]= 0.669626;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.17597;   Fcc[5]= 1.17597;   Fc[5]= 0.683558;   Fll[5]= 1.0011; 
              cafac[5]= 0.821593;  
              winpretag[5]= 3369.85 ;   wintag[5]= 688.38 ;   woutpretag[5]= 2768.65;   wouttag[5]= 588.842;  
            fbb[5]= 0.10985;   fcc[5]= 0.171033;   fc[5]= 0.092284;   fll[5]= 0.626833;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.20205;   Fcc[6]= 1.20205;   Fc[6]= 0.69872;   Fll[6]= 1.02331; 
              cafac[6]= 0.808105;  
              winpretag[6]= 73989 ;   wintag[6]= 9724.16 ;   woutpretag[6]= 59790.9;   wouttag[6]= 7891.26;  
            fbb[6]= 0.0672821;   fcc[6]= 0.124675;   fc[6]= 0.11162;   fll[6]= 0.696423;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.18871;   Fcc[7]= 1.18871;   Fc[7]= 0.690967;   Fll[7]= 1.01195; 
              cafac[7]= 0.824203;  
              winpretag[7]= 16127 ;   wintag[7]= 2756.32 ;   woutpretag[7]= 13292;   wouttag[7]= 2334.34;  
            fbb[7]= 0.0904145;   fcc[7]= 0.147205;   fc[7]= 0.101793;   fll[7]= 0.660587;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
      else if ( idsys==18   || sysname=="bmtag_down" ) { 
           Fbb[1]= 1.6519;   Fcc[1]= 1.6519;   Fc[1]= 0.771882;   Fll[1]= 1.00151; 
              cafac[1]= 0.984458;  
              winpretag[1]= 1.08281e+06 ;   wintag[1]= 37685.7 ;   woutpretag[1]= 1.06598e+06;   wouttag[1]= 35540.8;  
            fbb[1]= 0.0199173;   fcc[1]= 0.0579486;   fc[1]= 0.108131;   fll[1]= 0.814003;  
            fmcbb[1]= 0.0120572;   fmccc[1]= 0.0350801;   fmcc[1]= 0.140088;   fmcll[1]= 0.812775;  
           Fbb[2]= 1.60253;   Fcc[2]= 1.60253;   Fc[2]= 0.748815;   Fll[2]= 0.971581; 
              cafac[2]= 0.880677;  
              winpretag[2]= 257116 ;   wintag[2]= 18608.5 ;   woutpretag[2]= 226436;   wouttag[2]= 17682.8;  
            fbb[2]= 0.0482366;   fcc[2]= 0.115729;   fc[2]= 0.121473;   fll[2]= 0.714561;  
            fmcbb[2]= 0.0301003;   fmccc[2]= 0.0722167;   fmcc[2]= 0.16222;   fmcll[2]= 0.735463;  
           Fbb[3]= 1.5575;   Fcc[3]= 1.5575;   Fc[3]= 0.727775;   Fll[3]= 0.944282; 
              cafac[3]= 0.812;  
              winpretag[3]= 57862 ;   wintag[3]= 6172.49 ;   woutpretag[3]= 46983.9;   wouttag[3]= 5752.02;  
            fbb[3]= 0.0784575;   fcc[3]= 0.152809;   fc[3]= 0.118782;   fll[3]= 0.649951;  
            fmcbb[3]= 0.0503739;   fmccc[3]= 0.0981117;   fmcc[3]= 0.163213;   fmcll[3]= 0.688302;  
           Fbb[4]= 1.51515;   Fcc[4]= 1.51515;   Fc[4]= 0.707985;   Fll[4]= 0.918605; 
              cafac[4]= 0.830723;  
              winpretag[4]= 12757.2 ;   wintag[4]= 1816.28 ;   woutpretag[4]= 10597.7;   wouttag[4]= 1791.95;  
            fbb[4]= 0.108299;   fcc[4]= 0.178983;   fc[4]= 0.106603;   fll[4]= 0.606114;  
            fmcbb[4]= 0.0714774;   fmccc[4]= 0.118129;   fmcc[4]= 0.150573;   fmcll[4]= 0.659821;  
           Fbb[5]= 1.46724;   Fcc[5]= 1.46724;   Fc[5]= 0.685596;   Fll[5]= 0.889555; 
              cafac[5]= 0.813564;  
              winpretag[5]= 3369.85 ;   wintag[5]= 605.543 ;   woutpretag[5]= 2741.59;   wouttag[5]= 596.034;  
            fbb[5]= 0.137058;   fcc[5]= 0.213395;   fc[5]= 0.0925591;   fll[5]= 0.556988;  
            fmcbb[5]= 0.0934124;   fmccc[5]= 0.14544;   fmcc[5]= 0.135005;   fmcll[5]= 0.626142;  
           Fbb[6]= 1.54572;   Fcc[6]= 1.54572;   Fc[6]= 0.72227;   Fll[6]= 0.93714; 
              cafac[6]= 0.81477;  
              winpretag[6]= 73989 ;   wintag[6]= 8594.31 ;   woutpretag[6]= 60284;   wouttag[6]= 8149.84;  
            fbb[6]= 0.0865183;   fcc[6]= 0.16032;   fc[6]= 0.115382;   fll[6]= 0.63778;  
            fmcbb[6]= 0.0559728;   fmccc[6]= 0.103719;   fmcc[6]= 0.159749;   fmcll[6]= 0.68056;  
           Fbb[7]= 1.50488;   Fcc[7]= 1.50488;   Fc[7]= 0.703187;   Fll[7]= 0.912379; 
              cafac[7]= 0.825943;  
              winpretag[7]= 16127 ;   wintag[7]= 2421.83 ;   woutpretag[7]= 13320;   wouttag[7]= 2390.2;  
            fbb[7]= 0.114463;   fcc[7]= 0.186358;   fc[7]= 0.103593;   fll[7]= 0.595586;  
            fmcbb[7]= 0.0760609;   fmccc[7]= 0.123836;   fmcc[7]= 0.14732;   fmcll[7]= 0.652783;  
      }  
       else { cout << " WARNING: Ffactors request for unknown variation. Set to Nominal! " << idsys << "  " << sysname <<  endl; }

     return;


}

// FROM MUON FILE ON TWIKI [iwatson]

void SetWflavors_boosted_muon(int idsys, TString sysname, int ijet, double Wjets_in[5],double Wjets_out[5], double& canorm, int _mode) {
	//-
	// This function corrects the Wjet ligth and heavy flavor components and return the normalization based from Charge Assymmetry.
	//
	// NOTE: this function shoud always be called with PRETAGGED event counts. NEVER WITH TAGGED counts.
	//
	// INPUTS:
	// for idsys there are two input modes: the individul systematics (long list), or the combined systematics (4 checks).
	// Long list:
	// input idsys>0..999 or sysname are passed to the function GetFFactors_boosted_muon (see below) that returns the relevant HF factors.
	// This allows to use the full list of systematics.
	// Short list:
	// most analysers will use the total sum of all systematics. In that case only four checks up and two down are needed.
	// Complete set: idsys=2000 (up), -2000 (down) ---> checks Fbb + Fcc + Fc versus Fl
	//                     2001 (up), -2001 (down) ---> checks Fbb+Fcc versus Fc
	//                     2002 (up), -2002 (down) ---> checks Fll  OBSOLETE, not needed anymore, return nominal
	//                     2003 (up), -2003 (down) ---> check Fbb+xx% in 1,3,4,5 jetbin. --> see remark *
	//                     2004 (up), -2004 (down) ---> check Fc+xx% in 1,3,4,5 jetbin.  --> see remark *
	//                     2005 (up), -2005 (down) ---> check canorm. Please note that canorm is returned and NOT yet applied to Wjets_out.
	// * please note that the xx% per jetbin change should be applied UNCORRELATED to all jetbins. So, for example, when you use simultaneoulsy 
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
	//      please note that you have to do the xx% change per jetbin per flavor, which is NOT available in the long list. For these check you can use idsys=+/-2003,2004.
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
        bool use8TeV=true;//set to true if using 8TeV SFs derived for ttbar resonance resolved selection
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
	  if (use8TeV && idsys==0)   GetFFactors8TeV_boosted_muon( idsys, sysname,  Hbb,Hcc, Hc,  Hll, cafac,
							   fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
	  else 	GetFFactors_boosted_muon( idsys, sysname,  Hbb,Hcc, Hc,  Hll, cafac,
				  fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
	} else if (idsys*idsys<1999*1999) {
		cout << " FATAL  HF KFACTORS probabbly called with old defintion of systematics ids " << endl;
	} else  {                  // go for the reduced errors.
	  if (use8TeV && (abs(idsys)==2000|| abs(idsys)==2003 || abs(idsys)==2004 || abs(idsys)==2005)) 
	    GetFFactors8TeV_boosted_muon( idsys, "",  Hbb,Hcc, Hc,  Hll,cafac,
				  fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
	  else GetFFactors_boosted_muon( 0, "",  Hbb,Hcc, Hc,  Hll,cafac,
				 fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
	}
	//CHECK 8TeV SFs 
// 	cout<<"use8TeV "<<use8TeV<<endl;
// 	for (int i=1;i<=7;i++)
// 	cout<<Hbb[i]<<setw(10)<<Hcc[i]<<setw(10)<<Hc[i]<<setw(10)<<Hll[i]<<setw(10)<<cafac[i]<<endl;

	// errors on fraction, not on K factors
	double sigma_hf[2]={ 0.21 ,  - 0.19}; // error on ttoal HF
	double sigma_as[2] ={ 0.26 ,  - 0.28}; // error on fbb (and fcc) wrt to fc
	double sigma_canorm_up[9]=  {0,  0.089, 0.071, 0.084, 0.11, 0.18, 0.087,0.11,0};
	double sigma_canorm_down[9]={0, -0.073, -0.058,-0.068,-0.096,-0.16, -0.069,-0.098,0};
	double sigma_wbbwcc[9]=  {0,  0.25, 0.0, 0.11, 0.15, 0.29, 0.11,0.15,0};
	double sigma_wc[9]=      {0,  0.25, 0.0, 0.13, 0.15, 0.38, 0.13,0.15,0};

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
	//   cout << " fraction in " << fbb_in <<  " " << fcc_in << " " << fc_in << " " << fll_in << endl;
	double fbb_out=Wjets_out[0]/nw_in;
	double fcc_out=Wjets_out[1]/nw_in;
	double fc_out= Wjets_out[2]/nw_in;
	double fll_out=Wjets_out[3]/nw_in;
	double fhf_out=fbb_out+fcc_out+fc_out;
	if (!use8TeV) { //7TeV only
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
		fbb_out*=1+sigma_wbbwcc[ijet];
		fcc_out*=1+sigma_wbbwcc[ijet];
		double rest=fc_out+fll_out;
		double rest_new=1-fbb_out-fcc_out;
		fc_out*=rest_new/rest;
		fll_out*=rest_new/rest;
		if (ijet==1) cname="WbbWccjet1_up";
		if (ijet==3 || ijet==6) cname="WbbWccjet3_up";
		if (ijet==4 || ijet==7) cname="WbbWccjet4_up";
		if (ijet==5) cname="WbbWccjet5_up"; 
		GetFFactors_boosted_muon( -1, cname,  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		canorm=cafac[ijet];
	} else if (ijet!=2 && idsys==-2003) {
		fbb_out*=1-sigma_wbbwcc[ijet];
		fcc_out*=1-sigma_wbbwcc[ijet];
		double rest=fc_out+fll_out;
		double rest_new=1-fbb_out-fcc_out;
		fc_out*=rest_new/rest;
		fll_out*=rest_new/rest;
		if (ijet==1) cname="WbbWccjet1_down";
		if (ijet==3 || ijet==6) cname="WbbWccjet3_down";
		if (ijet==4 || ijet==7) cname="WbbWccjet4_down";
		if (ijet==5) cname="WbbWccjet5_down"; 
		GetFFactors_boosted_muon( -1, cname,  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		canorm=cafac[ijet];
	}
	if (ijet!=2 && idsys==2004) { // change c. at the cost of fbb+fcc+fl
		fc_out*=1+sigma_wc[ijet];
		double rest=fbb_out+fcc_out+fll_out;
		double rest_new=1-fc_out; 
		fbb_out*=rest_new/rest;
		fcc_out*=rest_new/rest;
		fll_out*=rest_new/rest;
		if (ijet==1) cname="Wcjet1_up";
		if (ijet==3 || ijet==6) cname="Wcjet3_up";
		if (ijet==4 || ijet==7) cname="Wcjet4_up";
		if (ijet==5) cname="Wcjet5_up"; // 
		GetFFactors_boosted_muon( -1, cname,  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		canorm=cafac[ijet];
	} else if (ijet!=2 && idsys==-2004) {
		fc_out*=1-sigma_wc[ijet];
		double rest=fbb_out+fcc_out+fll_out;
		double rest_new=1-fc_out; 
		fbb_out*=rest_new/rest;
		fcc_out*=rest_new/rest;
		fll_out*=rest_new/rest;
		if (ijet==1) cname="Wcjet1_down";
		if (ijet==3 || ijet==6) cname="Wcjet3_down";
		if (ijet==4 || ijet==7) cname="Wcjet4_down";
		if (ijet==5) cname="Wcjet5_down"; 
		GetFFactors_boosted_muon( -1, cname,  Hbb,Hcc, Hc,  Hll,cafac,
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
	} //end if (!use8TeV)
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
		GetFFactors_boosted_muon( 0, "",  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		double nwnom=winpretag[ijet];	
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
		GetFFactors_boosted_muon( 0, "",  Hbb,Hcc, Hc,  Hll,cafac,
			fbb, fcc, fc,  fll, winpretag, wintag, woutpretag, wouttag);
		double nwprenom=winpretag[ijet];

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





void GetFFactors_boosted_muon(int idsys, TString sysname, double Fbb[9],double Fcc[9],double Fc[9], double Fll[9], double cafac[9],
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
		//      GetFFactors_boosted_muon(idsys, "" , Fbb,Fcc,Fc,Fll); 
		//      .....
		//--
		// If you want to read this function by 'eye', i have two hints:
		// The HF factors are Fbb[2], Fcc[2] , Fc[2], Fll[2] for EACH jetbin, defined on the PRETAGGED sample.
		// The normalization factors are cafac[1], cafac[2], cafac[3], cafac[4], cafac[5] for jetbins 1,2,3,4,5 resp.
		// Ignore the other quantities which are used for additional checks and studies.
		double fmcbb[9];double fmccc[9];double fmcc[9];double fmcll[9];	




          Fbb[1]= 1.25694;   Fcc[1]= 1.25694;   Fc[1]= 0.989716;   Fll[1]= 0.987144; 
             cafac[1]= 1.05204;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29297e+06;   wouttag[1]= 83670.4;  
           fbb[1]= 0.0148818;   fcc[1]= 0.0434135;   fc[1]= 0.131966;   fll[1]= 0.809738;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.973572;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 37902.7;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2232;   Fcc[3]= 1.2232;   Fc[3]= 0.963152;   Fll[3]= 0.960649; 
             cafac[3]= 0.886748;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97025.1;   wouttag[3]= 11805.4;  
           fbb[3]= 0.0600158;   fcc[3]= 0.121519;   fc[3]= 0.148382;   fll[3]= 0.670083;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21025;   Fcc[4]= 1.21025;   Fc[4]= 0.952949;   Fll[4]= 0.950473; 
             cafac[4]= 0.951006;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22713.7;   wouttag[4]= 3661.08;  
           fbb[4]= 0.0853226;   fcc[4]= 0.143756;   fc[4]= 0.137189;   fll[4]= 0.633732;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19645;   Fcc[5]= 1.19645;   Fc[5]= 0.942087;   Fll[5]= 0.939639; 
             cafac[5]= 0.903508;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5631.03;   wouttag[5]= 1179.03;  
           fbb[5]= 0.112189;   fcc[5]= 0.167499;   fc[5]= 0.125899;   fll[5]= 0.594413;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.21975;   Fcc[6]= 1.21975;   Fc[6]= 0.960432;   Fll[6]= 0.957937; 
             cafac[6]= 0.899307;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 125483;   wouttag[6]= 16624.4;  
           fbb[6]= 0.0667573;   fcc[6]= 0.127449;   fc[6]= 0.145427;   fll[6]= 0.660367;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.20737;   Fcc[7]= 1.20737;   Fc[7]= 0.950681;   Fll[7]= 0.94821; 
             cafac[7]= 0.940464;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28323.3;   wouttag[7]= 4850.32;  
           fbb[7]= 0.0909332;   fcc[7]= 0.148714;   fc[7]= 0.134831;   fll[7]= 0.625521;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     if ( idsys==0   || sysname=="norminal" ) { 
// OK
     }  
     else if ( idsys==1   || sysname=="tt_up" ) { 
          Fbb[1]= 1.25474;   Fcc[1]= 1.25474;   Fc[1]= 0.949581;   Fll[1]= 0.993793; 
             cafac[1]= 1.04453;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.2766e+06;   wouttag[1]= 81322.4;  
           fbb[1]= 0.0148557;   fcc[1]= 0.0433374;   fc[1]= 0.126615;   fll[1]= 0.815192;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23798;   Fcc[2]= 1.23798;   Fc[2]= 0.936895;   Fll[2]= 0.980517; 
             cafac[2]= 0.965333;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 477540;   wouttag[2]= 37075.2;  
           fbb[2]= 0.0370273;   fcc[2]= 0.0891109;   fc[2]= 0.144967;   fll[2]= 0.728895;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.22329;   Fcc[3]= 1.22329;   Fc[3]= 0.92578;   Fll[3]= 0.968885; 
             cafac[3]= 0.879387;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 96219.6;   wouttag[3]= 11594.9;  
           fbb[3]= 0.06002;   fcc[3]= 0.121528;   fc[3]= 0.142624;   fll[3]= 0.675828;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21018;   Fcc[4]= 1.21018;   Fc[4]= 0.915858;   Fll[4]= 0.958501; 
             cafac[4]= 0.943324;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22530.3;   wouttag[4]= 3604.73;  
           fbb[4]= 0.0853178;   fcc[4]= 0.143748;   fc[4]= 0.131849;   fll[4]= 0.639085;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19626;   Fcc[5]= 1.19626;   Fc[5]= 0.905325;   Fll[5]= 0.947476; 
             cafac[5]= 0.896418;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5586.84;   wouttag[5]= 1162.75;  
           fbb[5]= 0.112171;   fcc[5]= 0.167472;   fc[5]= 0.120986;   fll[5]= 0.599371;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.2198;   Fcc[6]= 1.2198;   Fc[6]= 0.923137;   Fll[6]= 0.966118; 
             cafac[6]= 0.891898;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 124449;   wouttag[6]= 16341.2;  
           fbb[6]= 0.0667599;   fcc[6]= 0.127454;   fc[6]= 0.13978;   fll[6]= 0.666007;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.20727;   Fcc[7]= 1.20727;   Fc[7]= 0.913658;   Fll[7]= 0.956198; 
             cafac[7]= 0.932914;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28095.9;   wouttag[7]= 4777.62;  
           fbb[7]= 0.0909261;   fcc[7]= 0.148703;   fc[7]= 0.12958;   fll[7]= 0.630791;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==2   || sysname=="tt_down" ) { 
          Fbb[1]= 1.25927;   Fcc[1]= 1.25927;   Fc[1]= 1.03207;   Fll[1]= 0.980128; 
             cafac[1]= 1.06008;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.3105e+06;   wouttag[1]= 86185.2;  
           fbb[1]= 0.0149093;   fcc[1]= 0.0434938;   fc[1]= 0.137614;   fll[1]= 0.803983;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.2387;   Fcc[2]= 1.2387;   Fc[2]= 1.01521;   Fll[2]= 0.964117; 
             cafac[2]= 0.982393;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 485980;   wouttag[2]= 38788.9;  
           fbb[2]= 0.0370488;   fcc[2]= 0.0891627;   fc[2]= 0.157085;   fll[2]= 0.716703;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.22311;   Fcc[3]= 1.22311;   Fc[3]= 1.00244;   Fll[3]= 0.95199; 
             cafac[3]= 0.894622;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97886.6;   wouttag[3]= 12030.5;  
           fbb[3]= 0.0600114;   fcc[3]= 0.12151;   fc[3]= 0.154435;   fll[3]= 0.664044;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21032;   Fcc[4]= 1.21032;   Fc[4]= 0.991955;   Fll[4]= 0.94203; 
             cafac[4]= 0.95922;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22909.9;   wouttag[4]= 3721.33;  
           fbb[4]= 0.0853277;   fcc[4]= 0.143765;   fc[4]= 0.142804;   fll[4]= 0.628103;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19665;   Fcc[5]= 1.19665;   Fc[5]= 0.980757;   Fll[5]= 0.931395; 
             cafac[5]= 0.911089;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5678.28;   wouttag[5]= 1196.43;  
           fbb[5]= 0.112208;   fcc[5]= 0.167527;   fc[5]= 0.131067;   fll[5]= 0.589198;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.2197;   Fcc[6]= 1.2197;   Fc[6]= 0.999646;   Fll[6]= 0.949334; 
             cafac[6]= 0.907232;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 126589;   wouttag[6]= 16927.3;  
           fbb[6]= 0.0667547;   fcc[6]= 0.127444;   fc[6]= 0.151365;   fll[6]= 0.654436;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.20747;   Fcc[7]= 1.20747;   Fc[7]= 0.989617;   Fll[7]= 0.93981; 
             cafac[7]= 0.948537;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28566.4;   wouttag[7]= 4928.07;  
           fbb[7]= 0.0909407;   fcc[7]= 0.148727;   fc[7]= 0.140353;   fll[7]= 0.619979;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==19   || sysname=="Kbb_up" ) { 
          Fbb[1]= 1.35577;   Fcc[1]= 1.35577;   Fc[1]= 0.984899;   Fll[1]= 0.982339; 
             cafac[1]= 1.04691;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29297e+06;   wouttag[1]= 84509.5;  
           fbb[1]= 0.016052;   fcc[1]= 0.0468271;   fc[1]= 0.131324;   fll[1]= 0.805797;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.32817;   Fcc[2]= 1.32817;   Fc[2]= 0.964844;   Fll[2]= 0.962337; 
             cafac[2]= 0.963373;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 38693.6;  
           fbb[2]= 0.0397249;   fcc[2]= 0.095603;   fc[2]= 0.149292;   fll[2]= 0.71538;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.30595;   Fcc[3]= 1.30595;   Fc[3]= 0.948701;   Fll[3]= 0.946236; 
             cafac[3]= 0.873444;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97025.1;   wouttag[3]= 12104.5;  
           fbb[3]= 0.0640755;   fcc[3]= 0.129739;   fc[3]= 0.146155;   fll[3]= 0.66003;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.28705;   Fcc[4]= 1.28705;   Fc[4]= 0.934978;   Fll[4]= 0.932548; 
             cafac[4]= 0.933071;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22713.7;   wouttag[4]= 3765.07;  
           fbb[4]= 0.0907376;   fcc[4]= 0.15288;   fc[4]= 0.134602;   fll[4]= 0.621781;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.26711;   Fcc[5]= 1.26711;   Fc[5]= 0.920486;   Fll[5]= 0.918094; 
             cafac[5]= 0.882791;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5631.03;   wouttag[5]= 1212.56;  
           fbb[5]= 0.118814;   fcc[5]= 0.17739;   fc[5]= 0.123012;   fll[5]= 0.580783;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.3009;   Fcc[6]= 1.3009;   Fc[6]= 0.945033;   Fll[6]= 0.942577; 
             cafac[6]= 0.884888;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 125483;   wouttag[6]= 17066.3;  
           fbb[6]= 0.0711985;   fcc[6]= 0.135928;   fc[6]= 0.143095;   fll[6]= 0.649778;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.28287;   Fcc[7]= 1.28287;   Fc[7]= 0.931942;   Fll[7]= 0.92952; 
             cafac[7]= 0.921926;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28323.3;   wouttag[7]= 4989.11;  
           fbb[7]= 0.0966202;   fcc[7]= 0.158015;   fc[7]= 0.132173;   fll[7]= 0.613191;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==20   || sysname=="Kbb_down" ) { 
          Fbb[1]= 1.15714;   Fcc[1]= 1.15714;   Fc[1]= 0.994581;   Fll[1]= 0.991996; 
             cafac[1]= 1.05721;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29297e+06;   wouttag[1]= 82823.1;  
           fbb[1]= 0.0137001;   fcc[1]= 0.0399663;   fc[1]= 0.132615;   fll[1]= 0.813718;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.14656;   Fcc[2]= 1.14656;   Fc[2]= 0.985492;   Fll[2]= 0.982931; 
             cafac[2]= 0.983989;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 37095;  
           fbb[2]= 0.0342931;   fcc[2]= 0.0825308;   fc[2]= 0.152486;   fll[2]= 0.73069;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.1379;   Fcc[3]= 1.1379;   Fc[3]= 0.978049;   Fll[3]= 0.975508; 
             cafac[3]= 0.900464;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97025.1;   wouttag[3]= 11496.9;  
           fbb[3]= 0.0558305;   fcc[3]= 0.113045;   fc[3]= 0.150677;   fll[3]= 0.680448;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.13043;   Fcc[4]= 1.13043;   Fc[4]= 0.971625;   Fll[4]= 0.9691; 
             cafac[4]= 0.969644;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22713.7;   wouttag[4]= 3553.01;  
           fbb[4]= 0.0796954;   fcc[4]= 0.134275;   fc[4]= 0.139877;   fll[4]= 0.646152;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.1224;   Fcc[5]= 1.1224;   Fc[5]= 0.964727;   Fll[5]= 0.96222; 
             cafac[5]= 0.925221;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5631.03;   wouttag[5]= 1143.88;  
           fbb[5]= 0.105246;   fcc[5]= 0.157132;   fc[5]= 0.128924;   fll[5]= 0.608697;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.13592;   Fcc[6]= 1.13592;   Fc[6]= 0.976342;   Fll[6]= 0.973805; 
             cafac[6]= 0.914204;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 125483;   wouttag[6]= 16167.9;  
           fbb[6]= 0.0621691;   fcc[6]= 0.118689;   fc[6]= 0.147836;   fll[6]= 0.671306;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.12876;   Fcc[7]= 1.12876;   Fc[7]= 0.970189;   Fll[7]= 0.967668; 
             cafac[7]= 0.959762;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28323.3;   wouttag[7]= 4705.84;  
           fbb[7]= 0.0850128;   fcc[7]= 0.139032;   fc[7]= 0.137598;   fll[7]= 0.638357;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==21   || sysname=="Kc_up" ) { 
          Fbb[1]= 1.24791;   Fcc[1]= 1.24791;   Fc[1]= 1.0365;   Fll[1]= 0.98005; 
             cafac[1]= 1.04447;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29297e+06;   wouttag[1]= 85603;  
           fbb[1]= 0.0147749;   fcc[1]= 0.0431015;   fc[1]= 0.138205;   fll[1]= 0.803919;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.22816;   Fcc[2]= 1.22816;   Fc[2]= 1.0201;   Fll[2]= 0.964543; 
             cafac[2]= 0.965581;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 38396.9;  
           fbb[2]= 0.0367338;   fcc[2]= 0.0884045;   fc[2]= 0.157841;   fll[2]= 0.71702;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.21333;   Fcc[3]= 1.21333;   Fc[3]= 1.00778;   Fll[3]= 0.952893; 
             cafac[3]= 0.879589;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97025.1;   wouttag[3]= 11900.2;  
           fbb[3]= 0.0595313;   fcc[3]= 0.120538;   fc[3]= 0.155257;   fll[3]= 0.664674;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.20121;   Fcc[4]= 1.20121;   Fc[4]= 0.997712;   Fll[4]= 0.943374; 
             cafac[4]= 0.943903;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22713.7;   wouttag[4]= 3679.75;  
           fbb[4]= 0.0846853;   fcc[4]= 0.142683;   fc[4]= 0.143633;   fll[4]= 0.628999;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.18825;   Fcc[5]= 1.18825;   Fc[5]= 0.986946;   Fll[5]= 0.933195; 
             cafac[5]= 0.897312;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5631.03;   wouttag[5]= 1183.17;  
           fbb[5]= 0.11142;   fcc[5]= 0.16635;   fc[5]= 0.131894;   fll[5]= 0.590336;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.2101;   Fcc[6]= 1.2101;   Fc[6]= 1.0051;   Fll[6]= 0.950356; 
             cafac[6]= 0.892191;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 125483;   wouttag[6]= 16741.5;  
           fbb[6]= 0.066229;   fcc[6]= 0.126441;   fc[6]= 0.15219;   fll[6]= 0.655141;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.1985;   Fcc[7]= 1.1985;   Fc[7]= 0.995465;   Fll[7]= 0.941249; 
             cafac[7]= 0.933559;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28323.3;   wouttag[7]= 4872.96;  
           fbb[7]= 0.0902657;   fcc[7]= 0.147623;   fc[7]= 0.141183;   fll[7]= 0.620929;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==22   || sysname=="Kc_down" ) { 
          Fbb[1]= 1.2661;   Fcc[1]= 1.2661;   Fc[1]= 0.94225;   Fll[1]= 0.994342; 
             cafac[1]= 1.05971;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29297e+06;   wouttag[1]= 81709.7;  
           fbb[1]= 0.0149903;   fcc[1]= 0.04373;   fc[1]= 0.125637;   fll[1]= 0.815642;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.24866;   Fcc[2]= 1.24866;   Fc[2]= 0.929266;   Fll[2]= 0.98064; 
             cafac[2]= 0.981696;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 37400.3;  
           fbb[2]= 0.0373468;   fcc[2]= 0.0898799;   fc[2]= 0.143786;   fll[2]= 0.728987;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.23324;   Fcc[3]= 1.23324;   Fc[3]= 0.917792;   Fll[3]= 0.968532; 
             cafac[3]= 0.894025;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97025.1;   wouttag[3]= 11708.9;  
           fbb[3]= 0.0605083;   fcc[3]= 0.122516;   fc[3]= 0.141393;   fll[3]= 0.675582;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21942;   Fcc[4]= 1.21942;   Fc[4]= 0.907508;   Fll[4]= 0.957679; 
             cafac[4]= 0.958216;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22713.7;   wouttag[4]= 3642.12;  
           fbb[4]= 0.0859695;   fcc[4]= 0.144846;   fc[4]= 0.130647;   fll[4]= 0.638537;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.20477;   Fcc[5]= 1.20477;   Fc[5]= 0.896604;   Fll[5]= 0.946173; 
             cafac[5]= 0.909791;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5631.03;   wouttag[5]= 1174.83;  
           fbb[5]= 0.112969;   fcc[5]= 0.168664;   fc[5]= 0.119821;   fll[5]= 0.598546;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.22956;   Fcc[6]= 1.22956;   Fc[6]= 0.915051;   Fll[6]= 0.965639; 
             cafac[6]= 0.906539;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 125483;   wouttag[6]= 16505.5;  
           fbb[6]= 0.0672941;   fcc[6]= 0.128474;   fc[6]= 0.138555;   fll[6]= 0.665677;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.21636;   Fcc[7]= 1.21636;   Fc[7]= 0.90523;   Fll[7]= 0.955275; 
             cafac[7]= 0.947471;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28323.3;   wouttag[7]= 4827.34;  
           fbb[7]= 0.0916107;   fcc[7]= 0.149823;   fc[7]= 0.128385;   fll[7]= 0.630182;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==23   || sysname=="Kll_up" ) { 
          Fbb[1]= 1.25135;   Fcc[1]= 1.25135;   Fc[1]= 0.985315;   Fll[1]= 0.988176; 
             cafac[1]= 1.05313;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29297e+06;   wouttag[1]= 83420.3;  
           fbb[1]= 0.0148156;   fcc[1]= 0.0432204;   fc[1]= 0.13138;   fll[1]= 0.810584;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23341;   Fcc[2]= 1.23341;   Fc[2]= 0.971186;   Fll[2]= 0.974006; 
             cafac[2]= 0.975054;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 37800.5;  
           fbb[2]= 0.0368906;   fcc[2]= 0.088782;   fc[2]= 0.150273;   fll[2]= 0.724054;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2187;   Fcc[3]= 1.2187;   Fc[3]= 0.959605;   Fll[3]= 0.962391; 
             cafac[3]= 0.888356;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97025.1;   wouttag[3]= 11775.9;  
           fbb[3]= 0.0597948;   fcc[3]= 0.121072;   fc[3]= 0.147835;   fll[3]= 0.671298;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.20603;   Fcc[4]= 1.20603;   Fc[4]= 0.949629;   Fll[4]= 0.952387; 
             cafac[4]= 0.952921;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22713.7;   wouttag[4]= 3652.25;  
           fbb[4]= 0.0850254;   fcc[4]= 0.143255;   fc[4]= 0.136711;   fll[4]= 0.635008;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19254;   Fcc[5]= 1.19254;   Fc[5]= 0.939008;   Fll[5]= 0.941735; 
             cafac[5]= 0.905523;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5631.03;   wouttag[5]= 1176.36;  
           fbb[5]= 0.111823;   fcc[5]= 0.166952;   fc[5]= 0.125487;   fll[5]= 0.595738;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.21532;   Fcc[6]= 1.21532;   Fc[6]= 0.956946;   Fll[6]= 0.959725; 
             cafac[6]= 0.900986;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 125483;   wouttag[6]= 16583.2;  
           fbb[6]= 0.066515;   fcc[6]= 0.126987;   fc[6]= 0.144899;   fll[6]= 0.661599;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.20321;   Fcc[7]= 1.20321;   Fc[7]= 0.947412;   Fll[7]= 0.950162; 
             cafac[7]= 0.9424;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28323.3;   wouttag[7]= 4838.76;  
           fbb[7]= 0.0906205;   fcc[7]= 0.148203;   fc[7]= 0.134367;   fll[7]= 0.626809;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==24   || sysname=="Kll_down" ) { 
          Fbb[1]= 1.26258;   Fcc[1]= 1.26258;   Fc[1]= 0.994157;   Fll[1]= 0.986104; 
             cafac[1]= 1.05093;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29297e+06;   wouttag[1]= 83922.7;  
           fbb[1]= 0.0149486;   fcc[1]= 0.0436083;   fc[1]= 0.132559;   fll[1]= 0.808885;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.24328;   Fcc[2]= 1.24328;   Fc[2]= 0.978963;   Fll[2]= 0.971033; 
             cafac[2]= 0.972078;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 38005.8;  
           fbb[2]= 0.0371861;   fcc[2]= 0.089493;   fc[2]= 0.151476;   fll[2]= 0.721845;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.22774;   Fcc[3]= 1.22774;   Fc[3]= 0.966725;   Fll[3]= 0.958894; 
             cafac[3]= 0.885128;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97025.1;   wouttag[3]= 11835.1;  
           fbb[3]= 0.0602385;   fcc[3]= 0.12197;   fc[3]= 0.148932;   fll[3]= 0.668859;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21449;   Fcc[4]= 1.21449;   Fc[4]= 0.956292;   Fll[4]= 0.948546; 
             cafac[4]= 0.949078;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22713.7;   wouttag[4]= 3669.97;  
           fbb[4]= 0.0856219;   fcc[4]= 0.144261;   fc[4]= 0.13767;   fll[4]= 0.632447;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.20039;   Fcc[5]= 1.20039;   Fc[5]= 0.945187;   Fll[5]= 0.93753; 
             cafac[5]= 0.90148;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5631.03;   wouttag[5]= 1181.71;  
           fbb[5]= 0.112558;   fcc[5]= 0.16805;   fc[5]= 0.126313;   fll[5]= 0.593078;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.22421;   Fcc[6]= 1.22421;   Fc[6]= 0.963944;   Fll[6]= 0.956135; 
             cafac[6]= 0.897616;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 125483;   wouttag[6]= 16665.9;  
           fbb[6]= 0.0670014;   fcc[6]= 0.127915;   fc[6]= 0.145959;   fll[6]= 0.659125;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.21155;   Fcc[7]= 1.21155;   Fc[7]= 0.953973;   Fll[7]= 0.946245; 
             cafac[7]= 0.938514;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28323.3;   wouttag[7]= 4861.96;  
           fbb[7]= 0.0912481;   fcc[7]= 0.149229;   fc[7]= 0.135298;   fll[7]= 0.624225;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==33   || sysname=="WbbWccjet1_up" ) { 
          Fbb[1]= 1.25303;   Fcc[1]= 1.25303;   Fc[1]= 0.986634;   Fll[1]= 0.98407; 
             cafac[1]= 1.04795;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 79384.7 ;   woutpretag[1]= 2.28407e+06;   wouttag[1]= 85961.8;  
           fbb[1]= 0.0185443;   fcc[1]= 0.0540979;   fc[1]= 0.129956;   fll[1]= 0.797402;  
           fmcbb[1]= 0.0147996;   fmccc[1]= 0.0431738;   fmcc[1]= 0.131716;   fmcll[1]= 0.81031;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.973572;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 37902.7;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2232;   Fcc[3]= 1.2232;   Fc[3]= 0.963152;   Fll[3]= 0.960649; 
             cafac[3]= 0.886748;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97025.1;   wouttag[3]= 11805.4;  
           fbb[3]= 0.0600158;   fcc[3]= 0.121519;   fc[3]= 0.148382;   fll[3]= 0.670083;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21025;   Fcc[4]= 1.21025;   Fc[4]= 0.952949;   Fll[4]= 0.950473; 
             cafac[4]= 0.951006;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22713.7;   wouttag[4]= 3661.08;  
           fbb[4]= 0.0853226;   fcc[4]= 0.143756;   fc[4]= 0.137189;   fll[4]= 0.633732;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19645;   Fcc[5]= 1.19645;   Fc[5]= 0.942087;   Fll[5]= 0.939639; 
             cafac[5]= 0.903508;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5631.03;   wouttag[5]= 1179.03;  
           fbb[5]= 0.112189;   fcc[5]= 0.167499;   fc[5]= 0.125899;   fll[5]= 0.594413;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.21975;   Fcc[6]= 1.21975;   Fc[6]= 0.960432;   Fll[6]= 0.957937; 
             cafac[6]= 0.899307;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 125483;   wouttag[6]= 16624.4;  
           fbb[6]= 0.0667573;   fcc[6]= 0.127449;   fc[6]= 0.145427;   fll[6]= 0.660367;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.20737;   Fcc[7]= 1.20737;   Fc[7]= 0.950681;   Fll[7]= 0.94821; 
             cafac[7]= 0.940464;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28323.3;   wouttag[7]= 4850.32;  
           fbb[7]= 0.0909332;   fcc[7]= 0.148714;   fc[7]= 0.134831;   fll[7]= 0.625521;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==34   || sysname=="WbbWccjet1_down" ) { 
          Fbb[1]= 1.26088;   Fcc[1]= 1.26088;   Fc[1]= 0.992818;   Fll[1]= 0.990238; 
             cafac[1]= 1.05618;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 75347.4 ;   woutpretag[1]= 2.30199e+06;   wouttag[1]= 81346.5;  
           fbb[1]= 0.0111963;   fcc[1]= 0.0326622;   fc[1]= 0.13399;   fll[1]= 0.822152;  
           fmcbb[1]= 0.00887978;   fmccc[1]= 0.0259043;   fmcc[1]= 0.134959;   fmcll[1]= 0.830257;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.973572;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 37902.7;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2232;   Fcc[3]= 1.2232;   Fc[3]= 0.963152;   Fll[3]= 0.960649; 
             cafac[3]= 0.886748;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97025.1;   wouttag[3]= 11805.4;  
           fbb[3]= 0.0600158;   fcc[3]= 0.121519;   fc[3]= 0.148382;   fll[3]= 0.670083;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21025;   Fcc[4]= 1.21025;   Fc[4]= 0.952949;   Fll[4]= 0.950473; 
             cafac[4]= 0.951006;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22713.7;   wouttag[4]= 3661.08;  
           fbb[4]= 0.0853226;   fcc[4]= 0.143756;   fc[4]= 0.137189;   fll[4]= 0.633732;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19645;   Fcc[5]= 1.19645;   Fc[5]= 0.942087;   Fll[5]= 0.939639; 
             cafac[5]= 0.903508;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5631.03;   wouttag[5]= 1179.03;  
           fbb[5]= 0.112189;   fcc[5]= 0.167499;   fc[5]= 0.125899;   fll[5]= 0.594413;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.21975;   Fcc[6]= 1.21975;   Fc[6]= 0.960432;   Fll[6]= 0.957937; 
             cafac[6]= 0.899307;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 125483;   wouttag[6]= 16624.4;  
           fbb[6]= 0.0667573;   fcc[6]= 0.127449;   fc[6]= 0.145427;   fll[6]= 0.660367;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.20737;   Fcc[7]= 1.20737;   Fc[7]= 0.950681;   Fll[7]= 0.94821; 
             cafac[7]= 0.940464;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28323.3;   wouttag[7]= 4850.32;  
           fbb[7]= 0.0909332;   fcc[7]= 0.148714;   fc[7]= 0.134831;   fll[7]= 0.625521;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==35   || sysname=="WbbWccjet3_up" ) { 
          Fbb[1]= 1.25694;   Fcc[1]= 1.25694;   Fc[1]= 0.989716;   Fll[1]= 0.987144; 
             cafac[1]= 1.05204;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29297e+06;   wouttag[1]= 83670.4;  
           fbb[1]= 0.0148818;   fcc[1]= 0.0434135;   fc[1]= 0.131966;   fll[1]= 0.809738;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.973572;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 37902.7;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.21799;   Fcc[3]= 1.21799;   Fc[3]= 0.959048;   Fll[3]= 0.956556; 
             cafac[3]= 0.882926;  
             winpretag[3]= 109417 ;   wintag[3]= 12844.3 ;   woutpretag[3]= 96606.9;   wouttag[3]= 12218;  
           fbb[3]= 0.0663337;   fcc[3]= 0.134312;   fc[3]= 0.144917;   fll[3]= 0.654438;  
           fmcbb[3]= 0.0544615;   fmccc[3]= 0.110273;   fmcc[3]= 0.151105;   fmcll[3]= 0.68416;  
          Fbb[4]= 1.21025;   Fcc[4]= 1.21025;   Fc[4]= 0.952949;   Fll[4]= 0.950473; 
             cafac[4]= 0.951006;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22713.7;   wouttag[4]= 3661.08;  
           fbb[4]= 0.0853226;   fcc[4]= 0.143756;   fc[4]= 0.137189;   fll[4]= 0.633732;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19645;   Fcc[5]= 1.19645;   Fc[5]= 0.942087;   Fll[5]= 0.939639; 
             cafac[5]= 0.903508;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5631.03;   wouttag[5]= 1179.03;  
           fbb[5]= 0.112189;   fcc[5]= 0.167499;   fc[5]= 0.125899;   fll[5]= 0.594413;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.21568;   Fcc[6]= 1.21568;   Fc[6]= 0.95723;   Fll[6]= 0.954742; 
             cafac[6]= 0.896449;  
             winpretag[6]= 139533 ;   wintag[6]= 17594.4 ;   woutpretag[6]= 125084;   wouttag[6]= 17036.8;  
           fbb[6]= 0.0716797;   fcc[6]= 0.137442;   fc[6]= 0.142725;   fll[6]= 0.648153;  
           fmcbb[6]= 0.0589626;   fmccc[6]= 0.113057;   fmcc[6]= 0.149102;   fmcll[6]= 0.678878;  
          Fbb[7]= 1.20737;   Fcc[7]= 1.20737;   Fc[7]= 0.950681;   Fll[7]= 0.94821; 
             cafac[7]= 0.940464;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28323.3;   wouttag[7]= 4850.32;  
           fbb[7]= 0.0909332;   fcc[7]= 0.148714;   fc[7]= 0.134831;   fll[7]= 0.625521;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==36   || sysname=="WbbWccjet3_down" ) { 
          Fbb[1]= 1.25694;   Fcc[1]= 1.25694;   Fc[1]= 0.989716;   Fll[1]= 0.987144; 
             cafac[1]= 1.05204;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29297e+06;   wouttag[1]= 83670.4;  
           fbb[1]= 0.0148818;   fcc[1]= 0.0434135;   fc[1]= 0.131966;   fll[1]= 0.809738;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.973572;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 37902.7;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.22846;   Fcc[3]= 1.22846;   Fc[3]= 0.967291;   Fll[3]= 0.964777; 
             cafac[3]= 0.890636;  
             winpretag[3]= 109417 ;   wintag[3]= 11947.1 ;   woutpretag[3]= 97450.5;   wouttag[3]= 11385.5;  
           fbb[3]= 0.0536436;   fcc[3]= 0.108617;   fc[3]= 0.151876;   fll[3]= 0.685864;  
           fmcbb[3]= 0.0436673;   fmccc[3]= 0.0884171;   fmcc[3]= 0.157012;   fmcll[3]= 0.710904;  
          Fbb[4]= 1.21025;   Fcc[4]= 1.21025;   Fc[4]= 0.952949;   Fll[4]= 0.950473; 
             cafac[4]= 0.951006;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22713.7;   wouttag[4]= 3661.08;  
           fbb[4]= 0.0853226;   fcc[4]= 0.143756;   fc[4]= 0.137189;   fll[4]= 0.633732;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19645;   Fcc[5]= 1.19645;   Fc[5]= 0.942087;   Fll[5]= 0.939639; 
             cafac[5]= 0.903508;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5631.03;   wouttag[5]= 1179.03;  
           fbb[5]= 0.112189;   fcc[5]= 0.167499;   fc[5]= 0.125899;   fll[5]= 0.594413;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.22384;   Fcc[6]= 1.22384;   Fc[6]= 0.963657;   Fll[6]= 0.961152; 
             cafac[6]= 0.902204;  
             winpretag[6]= 139533 ;   wintag[6]= 16697.1 ;   woutpretag[6]= 125887;   wouttag[6]= 16206.5;  
           fbb[6]= 0.0618019;   fcc[6]= 0.11739;   fc[6]= 0.148147;   fll[6]= 0.672662;  
           fmcbb[6]= 0.0504982;   fmccc[6]= 0.0959187;   fmcc[6]= 0.153734;   fmcll[6]= 0.699849;  
          Fbb[7]= 1.20737;   Fcc[7]= 1.20737;   Fc[7]= 0.950681;   Fll[7]= 0.94821; 
             cafac[7]= 0.940464;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28323.3;   wouttag[7]= 4850.32;  
           fbb[7]= 0.0909332;   fcc[7]= 0.148714;   fc[7]= 0.134831;   fll[7]= 0.625521;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==37   || sysname=="WbbWccjet4_up" ) { 
          Fbb[1]= 1.25694;   Fcc[1]= 1.25694;   Fc[1]= 0.989716;   Fll[1]= 0.987144; 
             cafac[1]= 1.05204;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29297e+06;   wouttag[1]= 83670.4;  
           fbb[1]= 0.0148818;   fcc[1]= 0.0434135;   fc[1]= 0.131966;   fll[1]= 0.809738;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.973572;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 37902.7;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2232;   Fcc[3]= 1.2232;   Fc[3]= 0.963152;   Fll[3]= 0.960649; 
             cafac[3]= 0.886748;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97025.1;   wouttag[3]= 11805.4;  
           fbb[3]= 0.0600158;   fcc[3]= 0.121519;   fc[3]= 0.148382;   fll[3]= 0.670083;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.2014;   Fcc[4]= 1.2014;   Fc[4]= 0.945984;   Fll[4]= 0.943526; 
             cafac[4]= 0.94532;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3762.34 ;   woutpretag[4]= 22577.9;   wouttag[4]= 3869.82;  
           fbb[4]= 0.0974038;   fcc[4]= 0.164111;   fc[4]= 0.131417;   fll[4]= 0.607068;  
           fmcbb[4]= 0.0810752;   fmccc[4]= 0.1366;   fmcc[4]= 0.138921;   fmcll[4]= 0.643404;  
          Fbb[5]= 1.19645;   Fcc[5]= 1.19645;   Fc[5]= 0.942087;   Fll[5]= 0.939639; 
             cafac[5]= 0.903508;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5631.03;   wouttag[5]= 1179.03;  
           fbb[5]= 0.112189;   fcc[5]= 0.167499;   fc[5]= 0.125899;   fll[5]= 0.594413;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.2182;   Fcc[6]= 1.2182;   Fc[6]= 0.959214;   Fll[6]= 0.956721; 
             cafac[6]= 0.898185;  
             winpretag[6]= 139533 ;   wintag[6]= 17359.3 ;   woutpretag[6]= 125326;   wouttag[6]= 16830.4;  
           fbb[6]= 0.0688778;   fcc[6]= 0.131003;   fc[6]= 0.144415;   fll[6]= 0.655705;  
           fmcbb[6]= 0.0565405;   fmccc[6]= 0.107538;   fmcc[6]= 0.150555;   fmcll[6]= 0.685367;  
          Fbb[7]= 1.20037;   Fcc[7]= 1.20037;   Fc[7]= 0.945175;   Fll[7]= 0.942719; 
             cafac[7]= 0.936156;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4963.63 ;   woutpretag[7]= 28193.6;   wouttag[7]= 5054.66;  
           fbb[7]= 0.100474;   fcc[7]= 0.164815;   fc[7]= 0.130271;   fll[7]= 0.604441;  
           fmcbb[7]= 0.083702;   fmccc[7]= 0.137303;   fmcc[7]= 0.137827;   fmcll[7]= 0.641168;  
     }  
     else if ( idsys==38   || sysname=="WbbWccjet4_down" ) { 
          Fbb[1]= 1.25694;   Fcc[1]= 1.25694;   Fc[1]= 0.989716;   Fll[1]= 0.987144; 
             cafac[1]= 1.05204;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29297e+06;   wouttag[1]= 83670.4;  
           fbb[1]= 0.0148818;   fcc[1]= 0.0434135;   fc[1]= 0.131966;   fll[1]= 0.809738;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.973572;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 37902.7;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2232;   Fcc[3]= 1.2232;   Fc[3]= 0.963152;   Fll[3]= 0.960649; 
             cafac[3]= 0.886748;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97025.1;   wouttag[3]= 11805.4;  
           fbb[3]= 0.0600158;   fcc[3]= 0.121519;   fc[3]= 0.148382;   fll[3]= 0.670083;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21922;   Fcc[4]= 1.21922;   Fc[4]= 0.960018;   Fll[4]= 0.957523; 
             cafac[4]= 0.956846;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3335.12 ;   woutpretag[4]= 22853.2;   wouttag[4]= 3446.66;  
           fbb[4]= 0.0730622;   fcc[4]= 0.123099;   fc[4]= 0.143047;   fll[4]= 0.660792;  
           fmcbb[4]= 0.0599252;   fmccc[4]= 0.100965;   fmcc[4]= 0.149004;   fmcll[4]= 0.690106;  
          Fbb[5]= 1.19645;   Fcc[5]= 1.19645;   Fc[5]= 0.942087;   Fll[5]= 0.939639; 
             cafac[5]= 0.903508;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5631.03;   wouttag[5]= 1179.03;  
           fbb[5]= 0.112189;   fcc[5]= 0.167499;   fc[5]= 0.125899;   fll[5]= 0.594413;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.2213;   Fcc[6]= 1.2213;   Fc[6]= 0.961654;   Fll[6]= 0.959155; 
             cafac[6]= 0.900436;  
             winpretag[6]= 139533 ;   wintag[6]= 16932.1 ;   woutpretag[6]= 125641;   wouttag[6]= 16417.3;  
           fbb[6]= 0.0646315;   fcc[6]= 0.123887;   fc[6]= 0.146442;   fll[6]= 0.66504;  
           fmcbb[6]= 0.0529202;   fmccc[6]= 0.101438;   fmcc[6]= 0.152281;   fmcll[6]= 0.693361;  
          Fbb[7]= 1.21444;   Fcc[7]= 1.21444;   Fc[7]= 0.956252;   Fll[7]= 0.953767; 
             cafac[7]= 0.944862;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4536.41 ;   woutpretag[7]= 28455.8;   wouttag[7]= 4641.67;  
           fbb[7]= 0.0812811;   fcc[7]= 0.132426;   fc[7]= 0.139445;   fll[7]= 0.646849;  
           fmcbb[7]= 0.0669288;   fmccc[7]= 0.109043;   fmcc[7]= 0.145824;   fmcll[7]= 0.678204;  
     }  
     else if ( idsys==73   || sysname=="WbbWccjet5_up" ) { 
          Fbb[1]= 1.25694;   Fcc[1]= 1.25694;   Fc[1]= 0.989716;   Fll[1]= 0.987144; 
             cafac[1]= 1.05204;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29297e+06;   wouttag[1]= 83670.4;  
           fbb[1]= 0.0148818;   fcc[1]= 0.0434135;   fc[1]= 0.131966;   fll[1]= 0.809738;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.973572;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 37902.7;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2232;   Fcc[3]= 1.2232;   Fc[3]= 0.963152;   Fll[3]= 0.960649; 
             cafac[3]= 0.886748;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97025.1;   wouttag[3]= 11805.4;  
           fbb[3]= 0.0600158;   fcc[3]= 0.121519;   fc[3]= 0.148382;   fll[3]= 0.670083;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21025;   Fcc[4]= 1.21025;   Fc[4]= 0.952949;   Fll[4]= 0.950473; 
             cafac[4]= 0.951006;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22713.7;   wouttag[4]= 3661.08;  
           fbb[4]= 0.0853226;   fcc[4]= 0.143756;   fc[4]= 0.137189;   fll[4]= 0.633732;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.17601;   Fcc[5]= 1.17601;   Fc[5]= 0.925993;   Fll[5]= 0.923586; 
             cafac[5]= 0.888432;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1353.68 ;   woutpretag[5]= 5537.07;   wouttag[5]= 1308.99;  
           fbb[5]= 0.142252;   fcc[5]= 0.212382;   fc[5]= 0.1128;   fll[5]= 0.532566;  
           fmcbb[5]= 0.120961;   fmccc[5]= 0.180596;   fmcc[5]= 0.121815;   fmcll[5]= 0.576629;  
          Fbb[6]= 1.21879;   Fcc[6]= 1.21879;   Fc[6]= 0.959673;   Fll[6]= 0.957179; 
             cafac[6]= 0.898359;  
             winpretag[6]= 139533 ;   wintag[6]= 17298.1 ;   woutpretag[6]= 125351;   wouttag[6]= 16771.3;  
           fbb[6]= 0.0681849;   fcc[6]= 0.129559;   fc[6]= 0.144805;   fll[6]= 0.657452;  
           fmcbb[6]= 0.055945;   fmccc[6]= 0.106301;   fmcc[6]= 0.15089;   fmcll[6]= 0.686864;  
          Fbb[7]= 1.203;   Fcc[7]= 1.203;   Fc[7]= 0.947243;   Fll[7]= 0.944781; 
             cafac[7]= 0.936701;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4902.41 ;   woutpretag[7]= 28210;   wouttag[7]= 4996.19;  
           fbb[7]= 0.0973741;   fcc[7]= 0.158284;   fc[7]= 0.132026;   fll[7]= 0.612316;  
           fmcbb[7]= 0.0809428;   fmccc[7]= 0.131574;   fmcc[7]= 0.139379;   fmcll[7]= 0.648104;  
     }  
     else if ( idsys==74   || sysname=="WbbWccjet5_down" ) { 
          Fbb[1]= 1.25694;   Fcc[1]= 1.25694;   Fc[1]= 0.989716;   Fll[1]= 0.987144; 
             cafac[1]= 1.05204;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29297e+06;   wouttag[1]= 83670.4;  
           fbb[1]= 0.0148818;   fcc[1]= 0.0434135;   fc[1]= 0.131966;   fll[1]= 0.809738;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.973572;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 37902.7;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2232;   Fcc[3]= 1.2232;   Fc[3]= 0.963152;   Fll[3]= 0.960649; 
             cafac[3]= 0.886748;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97025.1;   wouttag[3]= 11805.4;  
           fbb[3]= 0.0600158;   fcc[3]= 0.121519;   fc[3]= 0.148382;   fll[3]= 0.670083;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21025;   Fcc[4]= 1.21025;   Fc[4]= 0.952949;   Fll[4]= 0.950473; 
             cafac[4]= 0.951006;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22713.7;   wouttag[4]= 3661.08;  
           fbb[4]= 0.0853226;   fcc[4]= 0.143756;   fc[4]= 0.137189;   fll[4]= 0.633732;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.21761;   Fcc[5]= 1.21761;   Fc[5]= 0.958751;   Fll[5]= 0.95626; 
             cafac[5]= 0.919667;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1048.9 ;   woutpretag[5]= 5731.74;   wouttag[5]= 1039.74;  
           fbb[5]= 0.0810633;   fcc[5]= 0.121028;   fc[5]= 0.139462;   fll[5]= 0.658447;  
           fmcbb[5]= 0.0665755;   fmccc[5]= 0.0993975;   fmcc[5]= 0.145462;   fmcll[5]= 0.688565;  
          Fbb[6]= 1.22072;   Fcc[6]= 1.22072;   Fc[6]= 0.961193;   Fll[6]= 0.958695; 
             cafac[6]= 0.90026;  
             winpretag[6]= 139533 ;   wintag[6]= 16993.3 ;   woutpretag[6]= 125616;   wouttag[6]= 16476.9;  
           fbb[6]= 0.0653275;   fcc[6]= 0.125336;   fc[6]= 0.14605;   fll[6]= 0.663286;  
           fmcbb[6]= 0.0535158;   fmccc[6]= 0.102675;   fmcc[6]= 0.151946;   fmcll[6]= 0.691863;  
          Fbb[7]= 1.21176;   Fcc[7]= 1.21176;   Fc[7]= 0.954144;   Fll[7]= 0.951665; 
             cafac[7]= 0.944285;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4597.63 ;   woutpretag[7]= 28438.4;   wouttag[7]= 4702.2;  
           fbb[7]= 0.0844454;   fcc[7]= 0.139075;   fc[7]= 0.137657;   fll[7]= 0.638822;  
           fmcbb[7]= 0.069688;   fmccc[7]= 0.114771;   fmcc[7]= 0.144273;   fmcll[7]= 0.671268;  
     }  
     else if ( idsys==41   || sysname=="Wcjet1_up" ) { 
          Fbb[1]= 1.25744;   Fcc[1]= 1.25744;   Fc[1]= 0.990108;   Fll[1]= 0.987535; 
             cafac[1]= 1.10119;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 87281.1 ;   woutpretag[1]= 2.40011e+06;   wouttag[1]= 98299.6;  
           fbb[1]= 0.0143151;   fcc[1]= 0.0417602;   fc[1]= 0.165023;   fll[1]= 0.778901;  
           fmcbb[1]= 0.0113843;   fmccc[1]= 0.0332106;   fmcc[1]= 0.166672;   fmcll[1]= 0.788733;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.973572;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 37902.7;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2232;   Fcc[3]= 1.2232;   Fc[3]= 0.963152;   Fll[3]= 0.960649; 
             cafac[3]= 0.886748;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97025.1;   wouttag[3]= 11805.4;  
           fbb[3]= 0.0600158;   fcc[3]= 0.121519;   fc[3]= 0.148382;   fll[3]= 0.670083;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21025;   Fcc[4]= 1.21025;   Fc[4]= 0.952949;   Fll[4]= 0.950473; 
             cafac[4]= 0.951006;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22713.7;   wouttag[4]= 3661.08;  
           fbb[4]= 0.0853226;   fcc[4]= 0.143756;   fc[4]= 0.137189;   fll[4]= 0.633732;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19645;   Fcc[5]= 1.19645;   Fc[5]= 0.942087;   Fll[5]= 0.939639; 
             cafac[5]= 0.903508;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5631.03;   wouttag[5]= 1179.03;  
           fbb[5]= 0.112189;   fcc[5]= 0.167499;   fc[5]= 0.125899;   fll[5]= 0.594413;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.21975;   Fcc[6]= 1.21975;   Fc[6]= 0.960432;   Fll[6]= 0.957937; 
             cafac[6]= 0.899307;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 125483;   wouttag[6]= 16624.4;  
           fbb[6]= 0.0667573;   fcc[6]= 0.127449;   fc[6]= 0.145427;   fll[6]= 0.660367;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.20737;   Fcc[7]= 1.20737;   Fc[7]= 0.950681;   Fll[7]= 0.94821; 
             cafac[7]= 0.940464;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28323.3;   wouttag[7]= 4850.32;  
           fbb[7]= 0.0909332;   fcc[7]= 0.148714;   fc[7]= 0.134831;   fll[7]= 0.625521;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==42   || sysname=="Wcjet1_down" ) { 
          Fbb[1]= 1.25644;   Fcc[1]= 1.25644;   Fc[1]= 0.989325;   Fll[1]= 0.986754; 
             cafac[1]= 1.00711;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 67451.1 ;   woutpretag[1]= 2.19506e+06;   wouttag[1]= 70301.6;  
           fbb[1]= 0.0154481;   fcc[1]= 0.0450655;   fc[1]= 0.0989357;   fll[1]= 0.840551;  
           fmcbb[1]= 0.0122951;   fmccc[1]= 0.0358675;   fmcc[1]= 0.100003;   fmcll[1]= 0.851834;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.973572;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 37902.7;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2232;   Fcc[3]= 1.2232;   Fc[3]= 0.963152;   Fll[3]= 0.960649; 
             cafac[3]= 0.886748;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97025.1;   wouttag[3]= 11805.4;  
           fbb[3]= 0.0600158;   fcc[3]= 0.121519;   fc[3]= 0.148382;   fll[3]= 0.670083;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21025;   Fcc[4]= 1.21025;   Fc[4]= 0.952949;   Fll[4]= 0.950473; 
             cafac[4]= 0.951006;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22713.7;   wouttag[4]= 3661.08;  
           fbb[4]= 0.0853226;   fcc[4]= 0.143756;   fc[4]= 0.137189;   fll[4]= 0.633732;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19645;   Fcc[5]= 1.19645;   Fc[5]= 0.942087;   Fll[5]= 0.939639; 
             cafac[5]= 0.903508;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5631.03;   wouttag[5]= 1179.03;  
           fbb[5]= 0.112189;   fcc[5]= 0.167499;   fc[5]= 0.125899;   fll[5]= 0.594413;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.21975;   Fcc[6]= 1.21975;   Fc[6]= 0.960432;   Fll[6]= 0.957937; 
             cafac[6]= 0.899307;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 125483;   wouttag[6]= 16624.4;  
           fbb[6]= 0.0667573;   fcc[6]= 0.127449;   fc[6]= 0.145427;   fll[6]= 0.660367;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.20737;   Fcc[7]= 1.20737;   Fc[7]= 0.950681;   Fll[7]= 0.94821; 
             cafac[7]= 0.940464;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28323.3;   wouttag[7]= 4850.32;  
           fbb[7]= 0.0909332;   fcc[7]= 0.148714;   fc[7]= 0.134831;   fll[7]= 0.625521;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==43   || sysname=="Wcjet3_up" ) { 
          Fbb[1]= 1.25694;   Fcc[1]= 1.25694;   Fc[1]= 0.989716;   Fll[1]= 0.987144; 
             cafac[1]= 1.05204;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29297e+06;   wouttag[1]= 83670.4;  
           fbb[1]= 0.0148818;   fcc[1]= 0.0434135;   fc[1]= 0.131966;   fll[1]= 0.809738;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.973572;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 37902.7;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.22427;   Fcc[3]= 1.22427;   Fc[3]= 0.963993;   Fll[3]= 0.961488; 
             cafac[3]= 0.912393;  
             winpretag[3]= 109417 ;   wintag[3]= 12731.2 ;   woutpretag[3]= 99831.1;   wouttag[3]= 12422.8;  
           fbb[3]= 0.0586461;   fcc[3]= 0.118746;   fc[3]= 0.167818;   fll[3]= 0.654791;  
           fmcbb[3]= 0.0479028;   fmccc[3]= 0.0969931;   fmcc[3]= 0.174086;   fmcll[3]= 0.681018;  
          Fbb[4]= 1.21025;   Fcc[4]= 1.21025;   Fc[4]= 0.952949;   Fll[4]= 0.950473; 
             cafac[4]= 0.951006;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22713.7;   wouttag[4]= 3661.08;  
           fbb[4]= 0.0853226;   fcc[4]= 0.143756;   fc[4]= 0.137189;   fll[4]= 0.633732;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19645;   Fcc[5]= 1.19645;   Fc[5]= 0.942087;   Fll[5]= 0.939639; 
             cafac[5]= 0.903508;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5631.03;   wouttag[5]= 1179.03;  
           fbb[5]= 0.112189;   fcc[5]= 0.167499;   fc[5]= 0.125899;   fll[5]= 0.594413;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.22058;   Fcc[6]= 1.22058;   Fc[6]= 0.961088;   Fll[6]= 0.958591; 
             cafac[6]= 0.918882;  
             winpretag[6]= 139533 ;   wintag[6]= 17481.3 ;   woutpretag[6]= 128214;   wouttag[6]= 17264.3;  
           fbb[6]= 0.0656911;   fcc[6]= 0.125285;   fc[6]= 0.16062;   fll[6]= 0.648404;  
           fmcbb[6]= 0.0538195;   fmccc[6]= 0.102644;   fmcc[6]= 0.167123;   fmcll[6]= 0.676414;  
          Fbb[7]= 1.20737;   Fcc[7]= 1.20737;   Fc[7]= 0.950681;   Fll[7]= 0.94821; 
             cafac[7]= 0.940464;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28323.3;   wouttag[7]= 4850.32;  
           fbb[7]= 0.0909332;   fcc[7]= 0.148714;   fc[7]= 0.134831;   fll[7]= 0.625521;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==44   || sysname=="Wcjet3_down" ) { 
          Fbb[1]= 1.25694;   Fcc[1]= 1.25694;   Fc[1]= 0.989716;   Fll[1]= 0.987144; 
             cafac[1]= 1.05204;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29297e+06;   wouttag[1]= 83670.4;  
           fbb[1]= 0.0148818;   fcc[1]= 0.0434135;   fc[1]= 0.131966;   fll[1]= 0.809738;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.973572;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 37902.7;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.22214;   Fcc[3]= 1.22214;   Fc[3]= 0.962312;   Fll[3]= 0.959811; 
             cafac[3]= 0.862546;  
             winpretag[3]= 109417 ;   wintag[3]= 12060.2 ;   woutpretag[3]= 94377;   wouttag[3]= 11222.7;  
           fbb[3]= 0.0613831;   fcc[3]= 0.124288;   fc[3]= 0.128979;   fll[3]= 0.68535;  
           fmcbb[3]= 0.050226;   fmccc[3]= 0.101697;   fmcc[3]= 0.134031;   fmcll[3]= 0.714046;  
          Fbb[4]= 1.21025;   Fcc[4]= 1.21025;   Fc[4]= 0.952949;   Fll[4]= 0.950473; 
             cafac[4]= 0.951006;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22713.7;   wouttag[4]= 3661.08;  
           fbb[4]= 0.0853226;   fcc[4]= 0.143756;   fc[4]= 0.137189;   fll[4]= 0.633732;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19645;   Fcc[5]= 1.19645;   Fc[5]= 0.942087;   Fll[5]= 0.939639; 
             cafac[5]= 0.903508;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5631.03;   wouttag[5]= 1179.03;  
           fbb[5]= 0.112189;   fcc[5]= 0.167499;   fc[5]= 0.125899;   fll[5]= 0.594413;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.21892;   Fcc[6]= 1.21892;   Fc[6]= 0.959778;   Fll[6]= 0.957284; 
             cafac[6]= 0.880574;  
             winpretag[6]= 139533 ;   wintag[6]= 16810.2 ;   woutpretag[6]= 122869;   wouttag[6]= 16012;  
           fbb[6]= 0.0678221;   fcc[6]= 0.12961;   fc[6]= 0.130255;   fll[6]= 0.672313;  
           fmcbb[6]= 0.0556412;   fmccc[6]= 0.106332;   fmcc[6]= 0.135713;   fmcll[6]= 0.702313;  
          Fbb[7]= 1.20737;   Fcc[7]= 1.20737;   Fc[7]= 0.950681;   Fll[7]= 0.94821; 
             cafac[7]= 0.940464;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28323.3;   wouttag[7]= 4850.32;  
           fbb[7]= 0.0909332;   fcc[7]= 0.148714;   fc[7]= 0.134831;   fll[7]= 0.625521;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==45   || sysname=="Wcjet4_up" ) { 
          Fbb[1]= 1.25694;   Fcc[1]= 1.25694;   Fc[1]= 0.989716;   Fll[1]= 0.987144; 
             cafac[1]= 1.05204;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29297e+06;   wouttag[1]= 83670.4;  
           fbb[1]= 0.0148818;   fcc[1]= 0.0434135;   fc[1]= 0.131966;   fll[1]= 0.809738;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.973572;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 37902.7;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2232;   Fcc[3]= 1.2232;   Fc[3]= 0.963152;   Fll[3]= 0.960649; 
             cafac[3]= 0.886748;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97025.1;   wouttag[3]= 11805.4;  
           fbb[3]= 0.0600158;   fcc[3]= 0.121519;   fc[3]= 0.148382;   fll[3]= 0.670083;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21178;   Fcc[4]= 1.21178;   Fc[4]= 0.954157;   Fll[4]= 0.951678; 
             cafac[4]= 0.983737;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3627.58 ;   woutpretag[4]= 23495.5;   wouttag[4]= 3853.49;  
           fbb[4]= 0.083132;   fcc[4]= 0.140065;   fc[4]= 0.159341;   fll[4]= 0.617462;  
           fmcbb[4]= 0.0686032;   fmccc[4]= 0.115586;   fmcc[4]= 0.166996;   fmcll[4]= 0.648814;  
          Fbb[5]= 1.19645;   Fcc[5]= 1.19645;   Fc[5]= 0.942087;   Fll[5]= 0.939639; 
             cafac[5]= 0.903508;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5631.03;   wouttag[5]= 1179.03;  
           fbb[5]= 0.112189;   fcc[5]= 0.167499;   fc[5]= 0.125899;   fll[5]= 0.594413;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.22002;   Fcc[6]= 1.22002;   Fc[6]= 0.960642;   Fll[6]= 0.958146; 
             cafac[6]= 0.904976;  
             winpretag[6]= 139533 ;   wintag[6]= 17224.6 ;   woutpretag[6]= 126274;   wouttag[6]= 16789.9;  
           fbb[6]= 0.0663758;   fcc[6]= 0.12681;   fc[6]= 0.149246;   fll[6]= 0.657568;  
           fmcbb[6]= 0.0544056;   fmccc[6]= 0.103941;   fmcc[6]= 0.155361;   fmcll[6]= 0.686293;  
          Fbb[7]= 1.20858;   Fcc[7]= 1.20858;   Fc[7]= 0.951634;   Fll[7]= 0.949161; 
             cafac[7]= 0.965461;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4828.87 ;   woutpretag[7]= 29076.1;   wouttag[7]= 5044.54;  
           fbb[7]= 0.0892062;   fcc[7]= 0.1458;   fc[7]= 0.15235;   fll[7]= 0.612644;  
           fmcbb[7]= 0.073811;   fmccc[7]= 0.120638;   fmcc[7]= 0.160093;   fmcll[7]= 0.645458;  
     }  
     else if ( idsys==46   || sysname=="Wcjet4_down" ) { 
          Fbb[1]= 1.25694;   Fcc[1]= 1.25694;   Fc[1]= 0.989716;   Fll[1]= 0.987144; 
             cafac[1]= 1.05204;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29297e+06;   wouttag[1]= 83670.4;  
           fbb[1]= 0.0148818;   fcc[1]= 0.0434135;   fc[1]= 0.131966;   fll[1]= 0.809738;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.973572;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 37902.7;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2232;   Fcc[3]= 1.2232;   Fc[3]= 0.963152;   Fll[3]= 0.960649; 
             cafac[3]= 0.886748;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97025.1;   wouttag[3]= 11805.4;  
           fbb[3]= 0.0600158;   fcc[3]= 0.121519;   fc[3]= 0.148382;   fll[3]= 0.670083;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.20872;   Fcc[4]= 1.20872;   Fc[4]= 0.951744;   Fll[4]= 0.949271; 
             cafac[4]= 0.920458;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3469.88 ;   woutpretag[4]= 21984.1;   wouttag[4]= 3481.5;  
           fbb[4]= 0.0875076;   fcc[4]= 0.147438;   fc[4]= 0.115093;   fll[4]= 0.649962;  
           fmcbb[4]= 0.0723972;   fmccc[4]= 0.121979;   fmcc[4]= 0.120928;   fmcll[4]= 0.684696;  
          Fbb[5]= 1.19645;   Fcc[5]= 1.19645;   Fc[5]= 0.942087;   Fll[5]= 0.939639; 
             cafac[5]= 0.903508;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5631.03;   wouttag[5]= 1179.03;  
           fbb[5]= 0.112189;   fcc[5]= 0.167499;   fc[5]= 0.125899;   fll[5]= 0.594413;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.21948;   Fcc[6]= 1.21948;   Fc[6]= 0.960223;   Fll[6]= 0.957727; 
             cafac[6]= 0.893712;  
             winpretag[6]= 139533 ;   wintag[6]= 17066.9 ;   woutpretag[6]= 124702;   wouttag[6]= 16461;  
           fbb[6]= 0.0671387;   fcc[6]= 0.128089;   fc[6]= 0.141609;   fll[6]= 0.663164;  
           fmcbb[6]= 0.0550551;   fmccc[6]= 0.105035;   fmcc[6]= 0.147475;   fmcll[6]= 0.692435;  
          Fbb[7]= 1.20616;   Fcc[7]= 1.20616;   Fc[7]= 0.94973;   Fll[7]= 0.947262; 
             cafac[7]= 0.916774;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4671.17 ;   woutpretag[7]= 27609.9;   wouttag[7]= 4666.27;  
           fbb[7]= 0.0926568;   fcc[7]= 0.151623;   fc[7]= 0.117347;   fll[7]= 0.638373;  
           fmcbb[7]= 0.0768198;   fmccc[7]= 0.125707;   fmcc[7]= 0.123559;   fmcll[7]= 0.673914;  
     }  
     else if ( idsys==75   || sysname=="Wcjet5_up" ) { 
          Fbb[1]= 1.25694;   Fcc[1]= 1.25694;   Fc[1]= 0.989716;   Fll[1]= 0.987144; 
             cafac[1]= 1.05204;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29297e+06;   wouttag[1]= 83670.4;  
           fbb[1]= 0.0148818;   fcc[1]= 0.0434135;   fc[1]= 0.131966;   fll[1]= 0.809738;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.973572;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 37902.7;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2232;   Fcc[3]= 1.2232;   Fc[3]= 0.963152;   Fll[3]= 0.960649; 
             cafac[3]= 0.886748;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97025.1;   wouttag[3]= 11805.4;  
           fbb[3]= 0.0600158;   fcc[3]= 0.121519;   fc[3]= 0.148382;   fll[3]= 0.670083;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21025;   Fcc[4]= 1.21025;   Fc[4]= 0.952949;   Fll[4]= 0.950473; 
             cafac[4]= 0.951006;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22713.7;   wouttag[4]= 3661.08;  
           fbb[4]= 0.0853226;   fcc[4]= 0.143756;   fc[4]= 0.137189;   fll[4]= 0.633732;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.20053;   Fcc[5]= 1.20053;   Fc[5]= 0.945296;   Fll[5]= 0.94284; 
             cafac[5]= 0.978992;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1246.56 ;   woutpretag[5]= 6101.48;   wouttag[5]= 1313.81;  
           fbb[5]= 0.105973;   fcc[5]= 0.158218;   fc[5]= 0.174332;   fll[5]= 0.561477;  
           fmcbb[5]= 0.088272;   fmccc[5]= 0.131791;   fmcc[5]= 0.184421;   fmcll[5]= 0.595517;  
          Fbb[6]= 1.21994;   Fcc[6]= 1.21994;   Fc[6]= 0.960581;   Fll[6]= 0.958085; 
             cafac[6]= 0.903034;  
             winpretag[6]= 139533 ;   wintag[6]= 17191 ;   woutpretag[6]= 126003;   wouttag[6]= 16725.8;  
           fbb[6]= 0.0664682;   fcc[6]= 0.127022;   fc[6]= 0.147628;   fll[6]= 0.658882;  
           fmcbb[6]= 0.0544849;   fmccc[6]= 0.104121;   fmcc[6]= 0.153686;   fmcll[6]= 0.687707;  
          Fbb[7]= 1.20822;   Fcc[7]= 1.20822;   Fc[7]= 0.951355;   Fll[7]= 0.948883; 
             cafac[7]= 0.956773;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4795.3 ;   woutpretag[7]= 28814.5;   wouttag[7]= 4969.34;  
           fbb[7]= 0.0896235;   fcc[7]= 0.146768;   fc[7]= 0.144925;   fll[7]= 0.618684;  
           fmcbb[7]= 0.074178;   fmccc[7]= 0.121475;   fmcc[7]= 0.152335;   fmcll[7]= 0.652013;  
     }  
     else if ( idsys==76   || sysname=="Wcjet5_down" ) { 
          Fbb[1]= 1.25694;   Fcc[1]= 1.25694;   Fc[1]= 0.989716;   Fll[1]= 0.987144; 
             cafac[1]= 1.05204;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29297e+06;   wouttag[1]= 83670.4;  
           fbb[1]= 0.0148818;   fcc[1]= 0.0434135;   fc[1]= 0.131966;   fll[1]= 0.809738;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.973572;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481616;   wouttag[2]= 37902.7;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2232;   Fcc[3]= 1.2232;   Fc[3]= 0.963152;   Fll[3]= 0.960649; 
             cafac[3]= 0.886748;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97025.1;   wouttag[3]= 11805.4;  
           fbb[3]= 0.0600158;   fcc[3]= 0.121519;   fc[3]= 0.148382;   fll[3]= 0.670083;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21025;   Fcc[4]= 1.21025;   Fc[4]= 0.952949;   Fll[4]= 0.950473; 
             cafac[4]= 0.951006;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22713.7;   wouttag[4]= 3661.08;  
           fbb[4]= 0.0853226;   fcc[4]= 0.143756;   fc[4]= 0.137189;   fll[4]= 0.633732;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.1924;   Fcc[5]= 1.1924;   Fc[5]= 0.9389;   Fll[5]= 0.93646; 
             cafac[5]= 0.839238;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1156.01 ;   woutpretag[5]= 5230.47;   wouttag[5]= 1064.27;  
           fbb[5]= 0.118364;   fcc[5]= 0.176717;   fc[5]= 0.0777932;   fll[5]= 0.627126;  
           fmcbb[5]= 0.0992647;   fmccc[5]= 0.148203;   fmcc[5]= 0.0828557;   fmcll[5]= 0.669677;  
          Fbb[6]= 1.21956;   Fcc[6]= 1.21956;   Fc[6]= 0.960284;   Fll[6]= 0.957789; 
             cafac[6]= 0.895613;  
             winpretag[6]= 139533 ;   wintag[6]= 17100.4 ;   woutpretag[6]= 124968;   wouttag[6]= 16523.9;  
           fbb[6]= 0.0670464;   fcc[6]= 0.127876;   fc[6]= 0.143226;   fll[6]= 0.661851;  
           fmcbb[6]= 0.0549759;   fmccc[6]= 0.104854;   fmcc[6]= 0.14915;   fmcll[6]= 0.69102;  
          Fbb[7]= 1.20651;   Fcc[7]= 1.20651;   Fc[7]= 0.950007;   Fll[7]= 0.947539; 
             cafac[7]= 0.924723;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4704.74 ;   woutpretag[7]= 27849.2;   wouttag[7]= 4735.45;  
           fbb[7]= 0.0922411;   fcc[7]= 0.150658;   fc[7]= 0.124752;   fll[7]= 0.632349;  
           fmcbb[7]= 0.0764528;   fmccc[7]= 0.124871;   fmcc[7]= 0.131317;   fmcll[7]= 0.66736;  
     }  
     else if ( idsys==3   || sysname=="wt_up" ) { 
          Fbb[1]= 1.25556;   Fcc[1]= 1.25556;   Fc[1]= 0.927871;   Fll[1]= 0.997275; 
             cafac[1]= 1.04048;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.26779e+06;   wouttag[1]= 80085.5;  
           fbb[1]= 0.0148655;   fcc[1]= 0.0433659;   fc[1]= 0.12372;   fll[1]= 0.818048;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23963;   Fcc[2]= 1.23963;   Fc[2]= 0.916096;   Fll[2]= 0.98462; 
             cafac[2]= 0.960888;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 475341;   wouttag[2]= 36648.2;  
           fbb[2]= 0.0370767;   fcc[2]= 0.0892299;   fc[2]= 0.141749;   fll[2]= 0.731945;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.22504;   Fcc[3]= 1.22504;   Fc[3]= 0.905315;   Fll[3]= 0.973032; 
             cafac[3]= 0.8753;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 95772.5;   wouttag[3]= 11486.7;  
           fbb[3]= 0.0601059;   fcc[3]= 0.121702;   fc[3]= 0.139471;   fll[3]= 0.678721;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21173;   Fcc[4]= 1.21173;   Fc[4]= 0.895479;   Fll[4]= 0.96246; 
             cafac[4]= 0.940929;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22473.1;   wouttag[4]= 3583.3;  
           fbb[4]= 0.0854272;   fcc[4]= 0.143932;   fc[4]= 0.128915;   fll[4]= 0.641725;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19762;   Fcc[5]= 1.19762;   Fc[5]= 0.885053;   Fll[5]= 0.951255; 
             cafac[5]= 0.891395;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5555.54;   wouttag[5]= 1153.15;  
           fbb[5]= 0.112299;   fcc[5]= 0.167663;   fc[5]= 0.118277;   fll[5]= 0.601761;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.22149;   Fcc[6]= 1.22149;   Fc[6]= 0.902695;   Fll[6]= 0.970216; 
             cafac[6]= 0.888074;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 123916;   wouttag[6]= 16201.6;  
           fbb[6]= 0.0668528;   fcc[6]= 0.127631;   fc[6]= 0.136684;   fll[6]= 0.668831;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.20878;   Fcc[7]= 1.20878;   Fc[7]= 0.893301;   Fll[7]= 0.96012; 
             cafac[7]= 0.929943;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28006.4;   wouttag[7]= 4747.09;  
           fbb[7]= 0.09104;   fcc[7]= 0.148889;   fc[7]= 0.126693;   fll[7]= 0.633378;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==4   || sysname=="wt_down" ) { 
          Fbb[1]= 1.25829;   Fcc[1]= 1.25829;   Fc[1]= 1.05017;   Fll[1]= 0.977242; 
             cafac[1]= 1.06358;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.31812e+06;   wouttag[1]= 87252.3;  
           fbb[1]= 0.0148978;   fcc[1]= 0.0434601;   fc[1]= 0.140027;   fll[1]= 0.801615;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23706;   Fcc[2]= 1.23706;   Fc[2]= 1.03245;   Fll[2]= 0.960753; 
             cafac[2]= 0.986243;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 487884;   wouttag[2]= 39156;  
           fbb[2]= 0.0369999;   fcc[2]= 0.0890449;   fc[2]= 0.159752;   fll[2]= 0.714203;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.22142;   Fcc[3]= 1.22142;   Fc[3]= 1.0194;   Fll[3]= 0.948606; 
             cafac[3]= 0.898172;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 98275;   wouttag[3]= 12123.4;  
           fbb[3]= 0.0599282;   fcc[3]= 0.121342;   fc[3]= 0.157046;   fll[3]= 0.661683;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.2088;   Fcc[4]= 1.2088;   Fc[4]= 1.00887;   Fll[4]= 0.938809; 
             cafac[4]= 0.96101;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22952.7;   wouttag[4]= 3738.33;  
           fbb[4]= 0.0852209;   fcc[4]= 0.143585;   fc[4]= 0.145239;   fll[4]= 0.625955;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19531;   Fcc[5]= 1.19531;   Fc[5]= 0.997607;   Fll[5]= 0.928331; 
             cafac[5]= 0.915619;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5706.51;   wouttag[5]= 1204.9;  
           fbb[5]= 0.112082;   fcc[5]= 0.16734;   fc[5]= 0.133318;   fll[5]= 0.587259;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.21805;   Fcc[6]= 1.21805;   Fc[6]= 1.01659;   Fll[6]= 0.945994; 
             cafac[6]= 0.910508;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 127046;   wouttag[6]= 17046;  
           fbb[6]= 0.0666646;   fcc[6]= 0.127272;   fc[6]= 0.15393;   fll[6]= 0.652134;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.20599;   Fcc[7]= 1.20599;   Fc[7]= 1.00652;   Fll[7]= 0.936621; 
             cafac[7]= 0.950928;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28638.4;   wouttag[7]= 4953.03;  
           fbb[7]= 0.0908294;   fcc[7]= 0.148545;   fc[7]= 0.14275;   fll[7]= 0.617876;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==5   || sysname=="st_up" ) { 
          Fbb[1]= 1.23763;   Fcc[1]= 1.23763;   Fc[1]= 0.990436;   Fll[1]= 0.988119; 
             cafac[1]= 1.05208;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29308e+06;   wouttag[1]= 83500.8;  
           fbb[1]= 0.0146532;   fcc[1]= 0.0427466;   fc[1]= 0.132062;   fll[1]= 0.810538;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.22066;   Fcc[2]= 1.22066;   Fc[2]= 0.976857;   Fll[2]= 0.974572; 
             cafac[2]= 0.973331;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481497;   wouttag[2]= 37735.2;  
           fbb[2]= 0.0365095;   fcc[2]= 0.0878648;   fc[2]= 0.15115;   fll[2]= 0.724475;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.20685;   Fcc[3]= 1.20685;   Fc[3]= 0.965802;   Fll[3]= 0.963543; 
             cafac[3]= 0.885771;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 96918.2;   wouttag[3]= 11732.7;  
           fbb[3]= 0.0592135;   fcc[3]= 0.119895;   fc[3]= 0.14879;   fll[3]= 0.672102;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.19499;   Fcc[4]= 1.19499;   Fc[4]= 0.956314;   Fll[4]= 0.954076; 
             cafac[4]= 0.949409;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22675.6;   wouttag[4]= 3634.17;  
           fbb[4]= 0.0842473;   fcc[4]= 0.141945;   fc[4]= 0.137673;   fll[4]= 0.636135;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.18235;   Fcc[5]= 1.18235;   Fc[5]= 0.946196;   Fll[5]= 0.943982; 
             cafac[5]= 0.901524;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5618.67;   wouttag[5]= 1169.72;  
           fbb[5]= 0.110867;   fcc[5]= 0.165525;   fc[5]= 0.126448;   fll[5]= 0.59716;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.20369;   Fcc[6]= 1.20369;   Fc[6]= 0.963275;   Fll[6]= 0.961021; 
             cafac[6]= 0.898185;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 125326;   wouttag[6]= 16515.5;  
           fbb[6]= 0.0658785;   fcc[6]= 0.125771;   fc[6]= 0.145857;   fll[6]= 0.662493;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.19236;   Fcc[7]= 1.19236;   Fc[7]= 0.954202;   Fll[7]= 0.95197; 
             cafac[7]= 0.938794;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28273;   wouttag[7]= 4813.98;  
           fbb[7]= 0.0898027;   fcc[7]= 0.146866;   fc[7]= 0.13533;   fll[7]= 0.628001;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==6   || sysname=="st_down" ) { 
          Fbb[1]= 1.27448;   Fcc[1]= 1.27448;   Fc[1]= 0.989062;   Fll[1]= 0.986259; 
             cafac[1]= 1.05199;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29287e+06;   wouttag[1]= 83824.5;  
           fbb[1]= 0.0150895;   fcc[1]= 0.0440195;   fc[1]= 0.131879;   fll[1]= 0.809012;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.25434;   Fcc[2]= 1.25434;   Fc[2]= 0.973429;   Fll[2]= 0.97067; 
             cafac[2]= 0.97379;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481724;   wouttag[2]= 38054.8;  
           fbb[2]= 0.0375167;   fcc[2]= 0.0902888;   fc[2]= 0.15062;   fll[2]= 0.721575;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.23801;   Fcc[3]= 1.23801;   Fc[3]= 0.960753;   Fll[3]= 0.958029; 
             cafac[3]= 0.887634;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97122;   wouttag[3]= 11871.3;  
           fbb[3]= 0.060742;   fcc[3]= 0.12299;   fc[3]= 0.148012;   fll[3]= 0.668256;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.22403;   Fcc[4]= 1.22403;   Fc[4]= 0.949909;   Fll[4]= 0.947216; 
             cafac[4]= 0.952454;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22748.3;   wouttag[4]= 3685.48;  
           fbb[4]= 0.0862945;   fcc[4]= 0.145394;   fc[4]= 0.136751;   fll[4]= 0.631561;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.20918;   Fcc[5]= 1.20918;   Fc[5]= 0.93838;   Fll[5]= 0.93572; 
             cafac[5]= 0.905308;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5642.25;   wouttag[5]= 1187.46;  
           fbb[5]= 0.113382;   fcc[5]= 0.169281;   fc[5]= 0.125403;   fll[5]= 0.591933;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.23428;   Fcc[6]= 1.23428;   Fc[6]= 0.957861;   Fll[6]= 0.955146; 
             cafac[6]= 0.900326;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 125625;   wouttag[6]= 16723.2;  
           fbb[6]= 0.0675525;   fcc[6]= 0.128967;   fc[6]= 0.145037;   fll[6]= 0.658443;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.22093;   Fcc[7]= 1.22093;   Fc[7]= 0.9475;   Fll[7]= 0.944814; 
             cafac[7]= 0.941978;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28368.9;   wouttag[7]= 4883.27;  
           fbb[7]= 0.0919547;   fcc[7]= 0.150385;   fc[7]= 0.13438;   fll[7]= 0.62328;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==7   || sysname=="z_up" ) { 
          Fbb[1]= 1.25584;   Fcc[1]= 1.25584;   Fc[1]= 0.937667;   Fll[1]= 0.995667; 
             cafac[1]= 1.03627;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.2586e+06;   wouttag[1]= 80182.1;  
           fbb[1]= 0.0148688;   fcc[1]= 0.0433757;   fc[1]= 0.125026;   fll[1]= 0.816729;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23948;   Fcc[2]= 1.23948;   Fc[2]= 0.925449;   Fll[2]= 0.982693; 
             cafac[2]= 0.958484;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 474152;   wouttag[2]= 36677.5;  
           fbb[2]= 0.0370723;   fcc[2]= 0.0892192;   fc[2]= 0.143196;   fll[2]= 0.730513;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2248;   Fcc[3]= 1.2248;   Fc[3]= 0.91449;   Fll[3]= 0.971056; 
             cafac[3]= 0.871965;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 95407.5;   wouttag[3]= 11469.4;  
           fbb[3]= 0.0600942;   fcc[3]= 0.121678;   fc[3]= 0.140885;   fll[3]= 0.677343;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21155;   Fcc[4]= 1.21155;   Fc[4]= 0.904591;   Fll[4]= 0.960545; 
             cafac[4]= 0.939721;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22444.2;   wouttag[4]= 3584.95;  
           fbb[4]= 0.0854142;   fcc[4]= 0.14391;   fc[4]= 0.130227;   fll[4]= 0.640448;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19748;   Fcc[5]= 1.19748;   Fc[5]= 0.894092;   Fll[5]= 0.949397; 
             cafac[5]= 0.894448;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5574.57;   wouttag[5]= 1158.73;  
           fbb[5]= 0.112286;   fcc[5]= 0.167644;   fc[5]= 0.119485;   fll[5]= 0.600586;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.22127;   Fcc[6]= 1.22127;   Fc[6]= 0.911853;   Fll[6]= 0.968256; 
             cafac[6]= 0.885472;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 123553;   wouttag[6]= 16188.5;  
           fbb[6]= 0.0668406;   fcc[6]= 0.127608;   fc[6]= 0.138071;   fll[6]= 0.66748;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.20861;   Fcc[7]= 1.20861;   Fc[7]= 0.902398;   Fll[7]= 0.958217; 
             cafac[7]= 0.929668;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 27998.2;   wouttag[7]= 4753.55;  
           fbb[7]= 0.0910268;   fcc[7]= 0.148868;   fc[7]= 0.127983;   fll[7]= 0.632122;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==8   || sysname=="z_down" ) { 
          Fbb[1]= 1.25801;   Fcc[1]= 1.25801;   Fc[1]= 1.04036;   Fll[1]= 0.978852; 
             cafac[1]= 1.06783;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.32739e+06;   wouttag[1]= 87166.5;  
           fbb[1]= 0.0148945;   fcc[1]= 0.0434505;   fc[1]= 0.138719;   fll[1]= 0.802936;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23721;   Fcc[2]= 1.23721;   Fc[2]= 1.02316;   Fll[2]= 0.962667; 
             cafac[2]= 0.98866;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 489080;   wouttag[2]= 39128;  
           fbb[2]= 0.0370044;   fcc[2]= 0.0890558;   fc[2]= 0.158314;   fll[2]= 0.715626;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.22166;   Fcc[3]= 1.22166;   Fc[3]= 1.01029;   Fll[3]= 0.950566; 
             cafac[3]= 0.901547;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 98644.3;   wouttag[3]= 12141.8;  
           fbb[3]= 0.0599399;   fcc[3]= 0.121366;   fc[3]= 0.155644;   fll[3]= 0.66305;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.20899;   Fcc[4]= 1.20899;   Fc[4]= 0.999818;   Fll[4]= 0.94071; 
             cafac[4]= 0.962216;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22981.5;   wouttag[4]= 3736.67;  
           fbb[4]= 0.0852341;   fcc[4]= 0.143607;   fc[4]= 0.143936;   fll[4]= 0.627223;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19545;   Fcc[5]= 1.19545;   Fc[5]= 0.988623;   Fll[5]= 0.930177; 
             cafac[5]= 0.912468;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5686.87;   wouttag[5]= 1199.1;  
           fbb[5]= 0.112096;   fcc[5]= 0.167359;   fc[5]= 0.132118;   fll[5]= 0.588427;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.21828;   Fcc[6]= 1.21828;   Fc[6]= 1.0075;   Fll[6]= 0.947938; 
             cafac[6]= 0.913135;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 127412;   wouttag[6]= 17060;  
           fbb[6]= 0.0666769;   fcc[6]= 0.127296;   fc[6]= 0.152554;   fll[6]= 0.653474;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.20616;   Fcc[7]= 1.20616;   Fc[7]= 0.997481;   Fll[7]= 0.938511; 
             cafac[7]= 0.951178;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28646;   wouttag[7]= 4946.35;  
           fbb[7]= 0.0908427;   fcc[7]= 0.148566;   fc[7]= 0.141468;   fll[7]= 0.619122;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==9   || sysname=="di_up" ) { 
          Fbb[1]= 1.25551;   Fcc[1]= 1.25551;   Fc[1]= 0.988862;   Fll[1]= 0.987364; 
             cafac[1]= 1.05184;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29254e+06;   wouttag[1]= 83602.5;  
           fbb[1]= 0.0148649;   fcc[1]= 0.0433641;   fc[1]= 0.131852;   fll[1]= 0.809919;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23706;   Fcc[2]= 1.23706;   Fc[2]= 0.974327;   Fll[2]= 0.972852; 
             cafac[2]= 0.973196;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481430;   wouttag[2]= 37865.3;  
           fbb[2]= 0.0369998;   fcc[2]= 0.0890447;   fc[2]= 0.150759;   fll[2]= 0.723197;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.22204;   Fcc[3]= 1.22204;   Fc[3]= 0.962499;   Fll[3]= 0.961041; 
             cafac[3]= 0.886337;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 96980.1;   wouttag[3]= 11793.1;  
           fbb[3]= 0.0599586;   fcc[3]= 0.121403;   fc[3]= 0.148281;   fll[3]= 0.670357;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.20916;   Fcc[4]= 1.20916;   Fc[4]= 0.952353;   Fll[4]= 0.950911; 
             cafac[4]= 0.950613;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22704.3;   wouttag[4]= 3657.47;  
           fbb[4]= 0.0852458;   fcc[4]= 0.143627;   fc[4]= 0.137103;   fll[4]= 0.634025;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19544;   Fcc[5]= 1.19544;   Fc[5]= 0.941551;   Fll[5]= 0.940125; 
             cafac[5]= 0.903204;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5629.14;   wouttag[5]= 1177.99;  
           fbb[5]= 0.112095;   fcc[5]= 0.167358;   fc[5]= 0.125827;   fll[5]= 0.59472;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.2186;   Fcc[6]= 1.2186;   Fc[6]= 0.959795;   Fll[6]= 0.958341; 
             cafac[6]= 0.898907;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 125427;   wouttag[6]= 16607.4;  
           fbb[6]= 0.0666947;   fcc[6]= 0.12733;   fc[6]= 0.14533;   fll[6]= 0.660646;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.20629;   Fcc[7]= 1.20629;   Fc[7]= 0.950098;   Fll[7]= 0.948659; 
             cafac[7]= 0.940091;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28312.1;   wouttag[7]= 4845.65;  
           fbb[7]= 0.0908524;   fcc[7]= 0.148582;   fc[7]= 0.134748;   fll[7]= 0.625817;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==10   || sysname=="di_down" ) { 
          Fbb[1]= 1.25837;   Fcc[1]= 1.25837;   Fc[1]= 0.990571;   Fll[1]= 0.986924; 
             cafac[1]= 1.05223;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.2934e+06;   wouttag[1]= 83738.3;  
           fbb[1]= 0.0148988;   fcc[1]= 0.043463;   fc[1]= 0.13208;   fll[1]= 0.809558;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.2396;   Fcc[2]= 1.2396;   Fc[2]= 0.975791;   Fll[2]= 0.972199; 
             cafac[2]= 0.973948;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 481802;   wouttag[2]= 37940.2;  
           fbb[2]= 0.0370758;   fcc[2]= 0.0892275;   fc[2]= 0.150985;   fll[2]= 0.722711;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.22437;   Fcc[3]= 1.22437;   Fc[3]= 0.963805;   Fll[3]= 0.960256; 
             cafac[3]= 0.88716;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97070.1;   wouttag[3]= 11817.7;  
           fbb[3]= 0.060073;   fcc[3]= 0.121635;   fc[3]= 0.148482;   fll[3]= 0.66981;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21134;   Fcc[4]= 1.21134;   Fc[4]= 0.953545;   Fll[4]= 0.950035; 
             cafac[4]= 0.951399;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22723.1;   wouttag[4]= 3664.69;  
           fbb[4]= 0.0853995;   fcc[4]= 0.143886;   fc[4]= 0.137275;   fll[4]= 0.63344;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19746;   Fcc[5]= 1.19746;   Fc[5]= 0.942623;   Fll[5]= 0.939153; 
             cafac[5]= 0.903812;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5632.93;   wouttag[5]= 1180.06;  
           fbb[5]= 0.112284;   fcc[5]= 0.16764;   fc[5]= 0.12597;   fll[5]= 0.594105;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.2209;   Fcc[6]= 1.2209;   Fc[6]= 0.96107;   Fll[6]= 0.957532; 
             cafac[6]= 0.899708;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 125539;   wouttag[6]= 16641.4;  
           fbb[6]= 0.06682;   fcc[6]= 0.127569;   fc[6]= 0.145523;   fll[6]= 0.660088;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.20844;   Fcc[7]= 1.20844;   Fc[7]= 0.951264;   Fll[7]= 0.947762; 
             cafac[7]= 0.940836;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28334.5;   wouttag[7]= 4855;  
           fbb[7]= 0.0910141;   fcc[7]= 0.148847;   fc[7]= 0.134914;   fll[7]= 0.625225;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==11   || sysname=="qcd_up" ) { 
          Fbb[1]= 1.25424;   Fcc[1]= 1.25424;   Fc[1]= 0.809856;   Fll[1]= 1.01653; 
             cafac[1]= 1.01918;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.22135e+06;   wouttag[1]= 73475;  
           fbb[1]= 0.0148498;   fcc[1]= 0.0433203;   fc[1]= 0.107984;   fll[1]= 0.833846;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.24333;   Fcc[2]= 1.24333;   Fc[2]= 0.802813;   Fll[2]= 1.00769; 
             cafac[2]= 0.937437;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 463741;   wouttag[2]= 34340;  
           fbb[2]= 0.0371875;   fcc[2]= 0.0894964;   fc[2]= 0.12422;   fll[2]= 0.749096;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.22968;   Fcc[3]= 1.22968;   Fc[3]= 0.794;   Fll[3]= 0.99663; 
             cafac[3]= 0.854404;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 93486;   wouttag[3]= 10904.5;  
           fbb[3]= 0.0603337;   fcc[3]= 0.122163;   fc[3]= 0.122322;   fll[3]= 0.695181;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21562;   Fcc[4]= 1.21562;   Fc[4]= 0.784917;   Fll[4]= 0.985229; 
             cafac[4]= 0.917253;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 21907.6;   wouttag[4]= 3421.55;  
           fbb[4]= 0.0857012;   fcc[4]= 0.144394;   fc[4]= 0.112998;   fll[4]= 0.656906;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.20083;   Fcc[5]= 1.20083;   Fc[5]= 0.775367;   Fll[5]= 0.973243; 
             cafac[5]= 0.872277;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5436.39;   wouttag[5]= 1109.97;  
           fbb[5]= 0.1126;   fcc[5]= 0.168112;   fc[5]= 0.103619;   fll[5]= 0.61567;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.22594;   Fcc[6]= 1.22594;   Fc[6]= 0.791582;   Fll[6]= 0.993595; 
             cafac[6]= 0.86674;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 120939;   wouttag[6]= 15414.8;  
           fbb[6]= 0.067096;   fcc[6]= 0.128096;   fc[6]= 0.11986;   fll[6]= 0.684948;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.21253;   Fcc[7]= 1.21253;   Fc[7]= 0.782921;   Fll[7]= 0.982724; 
             cafac[7]= 0.907269;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 27323.6;   wouttag[7]= 4541.44;  
           fbb[7]= 0.0913219;   fcc[7]= 0.14935;   fc[7]= 0.111038;   fll[7]= 0.64829;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==12   || sysname=="qcd_down" ) { 
          Fbb[1]= 1.25948;   Fcc[1]= 1.25948;   Fc[1]= 1.15871;   Fll[1]= 0.95953; 
             cafac[1]= 1.0849;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.3646e+06;   wouttag[1]= 93867.9;  
           fbb[1]= 0.0149119;   fcc[1]= 0.0435012;   fc[1]= 0.1545;   fll[1]= 0.787087;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23368;   Fcc[2]= 1.23368;   Fc[2]= 1.13498;   Fll[2]= 0.939876; 
             cafac[2]= 1.00971;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 499491;   wouttag[2]= 41465.5;  
           fbb[2]= 0.0368988;   fcc[2]= 0.0888017;   fc[2]= 0.175616;   fll[2]= 0.698683;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2172;   Fcc[3]= 1.2172;   Fc[3]= 1.11982;   Fll[3]= 0.927323; 
             cafac[3]= 0.918969;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 100551;   wouttag[3]= 12702.8;  
           fbb[3]= 0.0597214;   fcc[3]= 0.120923;   fc[3]= 0.172518;   fll[3]= 0.646838;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.20527;   Fcc[4]= 1.20527;   Fc[4]= 1.10884;   Fll[4]= 0.918228; 
             cafac[4]= 0.984619;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 23516.5;   wouttag[4]= 3899.62;  
           fbb[4]= 0.0849715;   fcc[4]= 0.143165;   fc[4]= 0.159631;   fll[4]= 0.612233;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19239;   Fcc[5]= 1.19239;   Fc[5]= 1.09699;   Fll[5]= 0.908417; 
             cafac[5]= 0.934599;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5824.8;   wouttag[5]= 1247.77;  
           fbb[5]= 0.111808;   fcc[5]= 0.16693;   fc[5]= 0.1466;   fll[5]= 0.574662;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.21402;   Fcc[6]= 1.21402;   Fc[6]= 1.11689;   Fll[6]= 0.924895; 
             cafac[6]= 0.931748;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 130010;   wouttag[6]= 17829.2;  
           fbb[6]= 0.0664436;   fcc[6]= 0.12685;   fc[6]= 0.169117;   fll[6]= 0.637589;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.20258;   Fcc[7]= 1.20258;   Fc[7]= 1.10636;   Fll[7]= 0.916181; 
             cafac[7]= 0.973519;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 29318.8;   wouttag[7]= 5157.9;  
           fbb[7]= 0.0905727;   fcc[7]= 0.148125;   fc[7]= 0.156911;   fll[7]= 0.604392;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==13   || sysname=="jes_up" ) { 
          Fbb[1]= 1.42282;   Fcc[1]= 1.42282;   Fc[1]= 0.994306;   Fll[1]= 0.978602; 
             cafac[1]= 0.972798;  
             winpretag[1]= 2.32497e+06 ;   wintag[1]= 80928.8 ;   woutpretag[1]= 2.26173e+06;   wouttag[1]= 82624.6;  
           fbb[1]= 0.0157652;   fcc[1]= 0.046246;   fc[1]= 0.128992;   fll[1]= 0.808996;  
           fmcbb[1]= 0.0110803;   fmccc[1]= 0.0325031;   fmcc[1]= 0.129731;   fmcll[1]= 0.826686;  
          Fbb[2]= 1.38924;   Fcc[2]= 1.38924;   Fc[2]= 0.970841;   Fll[2]= 0.955508; 
             cafac[2]= 0.875904;  
             winpretag[2]= 548078 ;   wintag[2]= 39905.2 ;   woutpretag[2]= 480064;   wouttag[2]= 38341.8;  
           fbb[2]= 0.0391714;   fcc[2]= 0.0957872;   fc[2]= 0.149251;   fll[2]= 0.71579;  
           fmcbb[2]= 0.0281963;   fmccc[2]= 0.0689494;   fmcc[2]= 0.153734;   fmcll[2]= 0.74912;  
          Fbb[3]= 1.36221;   Fcc[3]= 1.36221;   Fc[3]= 0.951955;   Fll[3]= 0.93692; 
             cafac[3]= 0.764027;  
             winpretag[3]= 126797 ;   wintag[3]= 14018.6 ;   woutpretag[3]= 96876.6;   wouttag[3]= 11967.7;  
           fbb[3]= 0.0638567;   fcc[3]= 0.130686;   fc[3]= 0.148329;   fll[3]= 0.657128;  
           fmcbb[3]= 0.0468772;   fmccc[3]= 0.0959364;   fmcc[3]= 0.155816;   fmcll[3]= 0.701371;  
          Fbb[4]= 1.33967;   Fcc[4]= 1.33967;   Fc[4]= 0.936198;   Fll[4]= 0.921412; 
             cafac[4]= 0.782178;  
             winpretag[4]= 29029.5 ;   wintag[4]= 4218.47 ;   woutpretag[4]= 22706.3;   wouttag[4]= 3741.65;  
           fbb[4]= 0.0880951;   fcc[4]= 0.156663;   fc[4]= 0.137587;   fll[4]= 0.617654;  
           fmcbb[4]= 0.065759;   fmccc[4]= 0.116942;   fmcc[4]= 0.146964;   fmcll[4]= 0.670335;  
          Fbb[5]= 1.31425;   Fcc[5]= 1.31425;   Fc[5]= 0.918434;   Fll[5]= 0.903928; 
             cafac[5]= 0.671983;  
             winpretag[5]= 8017.11 ;   wintag[5]= 1528.76 ;   woutpretag[5]= 5387.36;   wouttag[5]= 1168.06;  
           fbb[5]= 0.121649;   fcc[5]= 0.179817;   fc[5]= 0.123581;   fll[5]= 0.574953;  
           fmcbb[5]= 0.0925617;   fmccc[5]= 0.136822;   fmcc[5]= 0.134556;   fmcll[5]= 0.636061;  
          Fbb[6]= 1.35575;   Fcc[6]= 1.35575;   Fc[6]= 0.947437;   Fll[6]= 0.932474; 
             cafac[6]= 0.76165;  
             winpretag[6]= 163844 ;   wintag[6]= 19765.9 ;   woutpretag[6]= 124792;   wouttag[6]= 16926.8;  
           fbb[6]= 0.0711199;   fcc[6]= 0.137824;   fc[6]= 0.145154;   fll[6]= 0.645902;  
           fmcbb[6]= 0.052458;   fmccc[6]= 0.101659;   fmcc[6]= 0.153207;   fmcll[6]= 0.692676;  
          Fbb[7]= 1.33408;   Fcc[7]= 1.33408;   Fc[7]= 0.932295;   Fll[7]= 0.917571; 
             cafac[7]= 0.755926;  
             winpretag[7]= 37046.6 ;   wintag[7]= 5747.23 ;   woutpretag[7]= 28004.5;   wouttag[7]= 4934.8;  
           fbb[7]= 0.0954659;   fcc[7]= 0.161749;   fc[7]= 0.134511;   fll[7]= 0.608274;  
           fmcbb[7]= 0.0715593;   fmccc[7]= 0.121244;   fmcc[7]= 0.144279;   fmcll[7]= 0.662918;  
     }  
     else if ( idsys==14   || sysname=="jes_down" ) { 
          Fbb[1]= 1.23702;   Fcc[1]= 1.23702;   Fc[1]= 0.909894;   Fll[1]= 1.00087; 
             cafac[1]= 1.12126;  
             winpretag[1]= 2.04224e+06 ;   wintag[1]= 74389.9 ;   woutpretag[1]= 2.28989e+06;   wouttag[1]= 82000.4;  
           fbb[1]= 0.0155365;   fcc[1]= 0.0449792;   fc[1]= 0.124275;   fll[1]= 0.81521;  
           fmcbb[1]= 0.0125596;   fmccc[1]= 0.0363609;   fmcc[1]= 0.136582;   fmcll[1]= 0.814498;  
          Fbb[2]= 1.22254;   Fcc[2]= 1.22254;   Fc[2]= 0.899242;   Fll[2]= 0.989157; 
             cafac[2]= 1.06241;  
             winpretag[2]= 447997 ;   wintag[2]= 34025.2 ;   woutpretag[2]= 475955;   wouttag[2]= 37203.9;  
           fbb[2]= 0.0386842;   fcc[2]= 0.0912733;   fc[2]= 0.139679;   fll[2]= 0.730363;  
           fmcbb[2]= 0.0316425;   fmccc[2]= 0.0746588;   fmcc[2]= 0.15533;   fmcll[2]= 0.738369;  
          Fbb[3]= 1.20902;   Fcc[3]= 1.20902;   Fc[3]= 0.889299;   Fll[3]= 0.97822; 
             cafac[3]= 1.01912;  
             winpretag[3]= 94590.2 ;   wintag[3]= 10926.4 ;   woutpretag[3]= 96399;   wouttag[3]= 11696.3;  
           fbb[3]= 0.0619305;   fcc[3]= 0.123918;   fc[3]= 0.137;   fll[3]= 0.677152;  
           fmcbb[3]= 0.0512237;   fmccc[3]= 0.102494;   fmcc[3]= 0.154053;   fmcll[3]= 0.692229;  
          Fbb[4]= 1.197;   Fcc[4]= 1.197;   Fc[4]= 0.880454;   Fll[4]= 0.96849; 
             cafac[4]= 1.11522;  
             winpretag[4]= 19673.5 ;   wintag[4]= 2981.83 ;   woutpretag[4]= 21940.1;   wouttag[4]= 3545.73;  
           fbb[4]= 0.0870048;   fcc[4]= 0.143424;   fc[4]= 0.124804;   fll[4]= 0.644767;  
           fmcbb[4]= 0.0726859;   fmccc[4]= 0.11982;   fmcc[4]= 0.14175;   fmcll[4]= 0.665744;  
          Fbb[5]= 1.18215;   Fcc[5]= 1.18215;   Fc[5]= 0.869531;   Fll[5]= 0.956475; 
             cafac[5]= 1.1024;  
             winpretag[5]= 4932.84 ;   wintag[5]= 970.11 ;   woutpretag[5]= 5437.95;   wouttag[5]= 1143.73;  
           fbb[5]= 0.116535;   fcc[5]= 0.170996;   fc[5]= 0.113659;   fll[5]= 0.59881;  
           fmcbb[5]= 0.098579;   fmccc[5]= 0.144648;   fmcc[5]= 0.130713;   fmcll[5]= 0.626059;  
          Fbb[6]= 1.20589;   Fcc[6]= 1.20589;   Fc[6]= 0.886994;   Fll[6]= 0.975684; 
             cafac[6]= 1.04061;  
             winpretag[6]= 119196 ;   wintag[6]= 14878.3 ;   woutpretag[6]= 124037;   wouttag[6]= 16346.4;  
           fbb[6]= 0.0684049;   fcc[6]= 0.129148;   fc[6]= 0.133986;   fll[6]= 0.66846;  
           fmcbb[6]= 0.0567258;   fmccc[6]= 0.107098;   fmcc[6]= 0.151057;   fmcll[6]= 0.685119;  
          Fbb[7]= 1.19399;   Fcc[7]= 1.19399;   Fc[7]= 0.878242;   Fll[7]= 0.966057; 
             cafac[7]= 1.11231;  
             winpretag[7]= 24606.3 ;   wintag[7]= 3951.94 ;   woutpretag[7]= 27369.8;   wouttag[7]= 4693.18;  
           fbb[7]= 0.092984;   fcc[7]= 0.149007;   fc[7]= 0.122548;   fll[7]= 0.635461;  
           fmcbb[7]= 0.0778767;   fmccc[7]= 0.124798;   fmcc[7]= 0.139537;   fmcll[7]= 0.657788;  
     }  
     else if ( idsys==79   || sysname=="leff_up" ) { 
          Fbb[1]= 1.25046;   Fcc[1]= 1.25046;   Fc[1]= 0.978765;   Fll[1]= 0.989292; 
             cafac[1]= 1.03797;  
             winpretag[1]= 2.20771e+06 ;   wintag[1]= 78410 ;   woutpretag[1]= 2.29153e+06;   wouttag[1]= 83120.6;  
           fbb[1]= 0.0148072;   fcc[1]= 0.043189;   fc[1]= 0.130595;   fll[1]= 0.811409;  
           fmcbb[1]= 0.0118414;   fmccc[1]= 0.0345385;   fmcc[1]= 0.133429;   fmcll[1]= 0.820191;  
          Fbb[2]= 1.23286;   Fcc[2]= 1.23286;   Fc[2]= 0.964994;   Fll[2]= 0.975372; 
             cafac[2]= 0.959931;  
             winpretag[2]= 501187 ;   wintag[2]= 37242 ;   woutpretag[2]= 481105;   wouttag[2]= 37683.8;  
           fbb[2]= 0.0368707;   fcc[2]= 0.0887405;   fc[2]= 0.149408;   fll[2]= 0.724981;  
           fmcbb[2]= 0.0299065;   fmccc[2]= 0.0719792;   fmcc[2]= 0.154828;   fmcll[2]= 0.743286;  
          Fbb[3]= 1.21826;   Fcc[3]= 1.21826;   Fc[3]= 0.953566;   Fll[3]= 0.963822; 
             cafac[3]= 0.873851;  
             winpretag[3]= 110861 ;   wintag[3]= 12561.2 ;   woutpretag[3]= 96876;   wouttag[3]= 11739.2;  
           fbb[3]= 0.0597704;   fcc[3]= 0.121019;   fc[3]= 0.146988;   fll[3]= 0.672223;  
           fmcbb[3]= 0.049062;   fmccc[3]= 0.0993371;   fmcc[3]= 0.154145;   fmcll[3]= 0.697456;  
          Fbb[4]= 1.2056;   Fcc[4]= 1.2056;   Fc[4]= 0.943653;   Fll[4]= 0.953802; 
             cafac[4]= 0.937149;  
             winpretag[4]= 24202.4 ;   wintag[4]= 3597.38 ;   woutpretag[4]= 22681.2;   wouttag[4]= 3643.28;  
           fbb[4]= 0.085;   fcc[4]= 0.143196;   fc[4]= 0.135914;   fll[4]= 0.63589;  
           fmcbb[4]= 0.0705044;   fmccc[4]= 0.118776;   fmcc[4]= 0.14403;   fmcll[4]= 0.66669;  
          Fbb[5]= 1.1921;   Fcc[5]= 1.1921;   Fc[5]= 0.933089;   Fll[5]= 0.943124; 
             cafac[5]= 0.889886;  
             winpretag[5]= 6315.54 ;   wintag[5]= 1217.56 ;   woutpretag[5]= 5620.1;   wouttag[5]= 1172.95;  
           fbb[5]= 0.111782;   fcc[5]= 0.166964;   fc[5]= 0.124747;   fll[5]= 0.596508;  
           fmcbb[5]= 0.0937687;   fmccc[5]= 0.140059;   fmcc[5]= 0.133692;   fmcll[5]= 0.632481;  
          Fbb[6]= 1.21489;   Fcc[6]= 1.21489;   Fc[6]= 0.950924;   Fll[6]= 0.961151; 
             cafac[6]= 0.886203;  
             winpretag[6]= 141379 ;   wintag[6]= 17376.1 ;   woutpretag[6]= 125290;   wouttag[6]= 16534.3;  
           fbb[6]= 0.0664905;   fcc[6]= 0.126936;   fc[6]= 0.144065;   fll[6]= 0.662508;  
           fmcbb[6]= 0.0547298;   fmccc[6]= 0.104484;   fmcc[6]= 0.1515;   fmcll[6]= 0.689286;  
          Fbb[7]= 1.20278;   Fcc[7]= 1.20278;   Fc[7]= 0.941447;   Fll[7]= 0.951572; 
             cafac[7]= 0.926661;  
             winpretag[7]= 30517.9 ;   wintag[7]= 4814.94 ;   woutpretag[7]= 28279.7;   wouttag[7]= 4826.45;  
           fbb[7]= 0.090592;   fcc[7]= 0.148158;   fc[7]= 0.133582;   fll[7]= 0.627667;  
           fmcbb[7]= 0.0753189;   fmccc[7]= 0.12318;   fmcc[7]= 0.14189;   fmcll[7]= 0.659611;  
     }  
     else if ( idsys==80   || sysname=="leff_down" ) { 
          Fbb[1]= 1.2634;   Fcc[1]= 1.2634;   Fc[1]= 1.00077;   Fll[1]= 0.984984; 
             cafac[1]= 1.06645;  
             winpretag[1]= 2.1514e+06 ;   wintag[1]= 76322.2 ;   woutpretag[1]= 2.29437e+06;   wouttag[1]= 84221.5;  
           fbb[1]= 0.0149561;   fcc[1]= 0.0436374;   fc[1]= 0.133347;   fll[1]= 0.80806;  
           fmcbb[1]= 0.0118379;   fmccc[1]= 0.0345396;   fmcc[1]= 0.133244;   fmcll[1]= 0.820378;  
          Fbb[2]= 1.24376;   Fcc[2]= 1.24376;   Fc[2]= 0.985209;   Fll[2]= 0.969668; 
             cafac[2]= 0.987559;  
             winpretag[2]= 488192 ;   wintag[2]= 36252.8 ;   woutpretag[2]= 482119;   wouttag[2]= 38121.7;  
           fbb[2]= 0.0372041;   fcc[2]= 0.0895297;   fc[2]= 0.152344;   fll[2]= 0.720922;  
           fmcbb[2]= 0.0299126;   fmccc[2]= 0.0719831;   fmcc[2]= 0.154632;   fmcll[2]= 0.743473;  
          Fbb[3]= 1.22811;   Fcc[3]= 1.22811;   Fc[3]= 0.972813;   Fll[3]= 0.957469; 
             cafac[3]= 0.899978;  
             winpretag[3]= 107972 ;   wintag[3]= 12230.2 ;   woutpretag[3]= 97172.8;   wouttag[3]= 11871.5;  
           fbb[3]= 0.0602597;   fcc[3]= 0.122017;   fc[3]= 0.149783;   fll[3]= 0.66794;  
           fmcbb[3]= 0.0490669;   fmccc[3]= 0.0993533;   fmcc[3]= 0.153969;   fmcll[3]= 0.697611;  
          Fbb[4]= 1.21486;   Fcc[4]= 1.21486;   Fc[4]= 0.962319;   Fll[4]= 0.94714; 
             cafac[4]= 0.965221;  
             winpretag[4]= 23565.4 ;   wintag[4]= 3500.08 ;   woutpretag[4]= 22745.9;   wouttag[4]= 3678.81;  
           fbb[4]= 0.0856428;   fcc[4]= 0.144313;   fc[4]= 0.138471;   fll[4]= 0.631573;  
           fmcbb[4]= 0.0704958;   fmccc[4]= 0.11879;   fmcc[4]= 0.143893;   fmcll[4]= 0.666821;  
          Fbb[5]= 1.20077;   Fcc[5]= 1.20077;   Fc[5]= 0.951157;   Fll[5]= 0.936154; 
             cafac[5]= 0.917484;  
             winpretag[5]= 6149.28 ;   wintag[5]= 1185.01 ;   woutpretag[5]= 5641.87;   wouttag[5]= 1185.09;  
           fbb[5]= 0.112594;   fcc[5]= 0.168028;   fc[5]= 0.127058;   fll[5]= 0.59232;  
           fmcbb[5]= 0.0937679;   fmccc[5]= 0.139933;   fmcc[5]= 0.133583;   fmcll[5]= 0.632716;  
          Fbb[6]= 1.22458;   Fcc[6]= 1.22458;   Fc[6]= 0.970016;   Fll[6]= 0.954716; 
             cafac[6]= 0.91275;  
             winpretag[6]= 137687 ;   wintag[6]= 16915.3 ;   woutpretag[6]= 125674;   wouttag[6]= 16714.4;  
           fbb[6]= 0.0670225;   fcc[6]= 0.127959;   fc[6]= 0.146796;   fll[6]= 0.658222;  
           fmcbb[6]= 0.054731;   fmccc[6]= 0.104492;   fmcc[6]= 0.151334;   fmcll[6]= 0.689443;  
          Fbb[7]= 1.21192;   Fcc[7]= 1.21192;   Fc[7]= 0.959988;   Fll[7]= 0.944845; 
             cafac[7]= 0.954623;  
             winpretag[7]= 29714.7 ;   wintag[7]= 4685.1 ;   woutpretag[7]= 28366.4;   wouttag[7]= 4874.11;  
           fbb[7]= 0.0912719;   fcc[7]= 0.149266;   fc[7]= 0.136087;   fll[7]= 0.623374;  
           fmcbb[7]= 0.0753119;   fmccc[7]= 0.123165;   fmcc[7]= 0.14176;   fmcll[7]= 0.659763;  
     }  
     else if ( idsys==81   || sysname=="ca_up" ) { 
          Fbb[1]= 1.25694;   Fcc[1]= 1.25694;   Fc[1]= 0.989716;   Fll[1]= 0.987144; 
             cafac[1]= 1.056;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.30161e+06;   wouttag[1]= 83985.8;  
           fbb[1]= 0.0148818;   fcc[1]= 0.0434135;   fc[1]= 0.131966;   fll[1]= 0.809738;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.9845;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 487022;   wouttag[2]= 38328.2;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2232;   Fcc[3]= 1.2232;   Fc[3]= 0.963152;   Fll[3]= 0.960649; 
             cafac[3]= 0.903182;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 98823.2;   wouttag[3]= 12024.1;  
           fbb[3]= 0.0600158;   fcc[3]= 0.121519;   fc[3]= 0.148382;   fll[3]= 0.670083;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21025;   Fcc[4]= 1.21025;   Fc[4]= 0.952949;   Fll[4]= 0.950473; 
             cafac[4]= 0.988593;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 23611.5;   wouttag[4]= 3805.78;  
           fbb[4]= 0.0853226;   fcc[4]= 0.143756;   fc[4]= 0.137189;   fll[4]= 0.633732;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19645;   Fcc[5]= 1.19645;   Fc[5]= 0.942087;   Fll[5]= 0.939639; 
             cafac[5]= 0.993761;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 6193.52;   wouttag[5]= 1296.8;  
           fbb[5]= 0.112189;   fcc[5]= 0.167499;   fc[5]= 0.125899;   fll[5]= 0.594413;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.21975;   Fcc[6]= 1.21975;   Fc[6]= 0.960432;   Fll[6]= 0.957937; 
             cafac[6]= 0.9144;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 127589;   wouttag[6]= 16903.4;  
           fbb[6]= 0.0667573;   fcc[6]= 0.127449;   fc[6]= 0.145427;   fll[6]= 0.660367;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.20737;   Fcc[7]= 1.20737;   Fc[7]= 0.950681;   Fll[7]= 0.94821; 
             cafac[7]= 0.975851;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 29389;   wouttag[7]= 5032.83;  
           fbb[7]= 0.0909332;   fcc[7]= 0.148714;   fc[7]= 0.134831;   fll[7]= 0.625521;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==82   || sysname=="ca_down" ) { 
          Fbb[1]= 1.25694;   Fcc[1]= 1.25694;   Fc[1]= 0.989716;   Fll[1]= 0.987144; 
             cafac[1]= 1.04807;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.28433e+06;   wouttag[1]= 83355.1;  
           fbb[1]= 0.0148818;   fcc[1]= 0.0434135;   fc[1]= 0.131966;   fll[1]= 0.809738;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.962643;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 476210;   wouttag[2]= 37477.3;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2232;   Fcc[3]= 1.2232;   Fc[3]= 0.963152;   Fll[3]= 0.960649; 
             cafac[3]= 0.870314;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 95226.9;   wouttag[3]= 11586.6;  
           fbb[3]= 0.0600158;   fcc[3]= 0.121519;   fc[3]= 0.148382;   fll[3]= 0.670083;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21025;   Fcc[4]= 1.21025;   Fc[4]= 0.952949;   Fll[4]= 0.950473; 
             cafac[4]= 0.913419;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 21816;   wouttag[4]= 3516.38;  
           fbb[4]= 0.0853226;   fcc[4]= 0.143756;   fc[4]= 0.137189;   fll[4]= 0.633732;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19645;   Fcc[5]= 1.19645;   Fc[5]= 0.942087;   Fll[5]= 0.939639; 
             cafac[5]= 0.813256;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5068.54;   wouttag[5]= 1061.25;  
           fbb[5]= 0.112189;   fcc[5]= 0.167499;   fc[5]= 0.125899;   fll[5]= 0.594413;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.21975;   Fcc[6]= 1.21975;   Fc[6]= 0.960432;   Fll[6]= 0.957937; 
             cafac[6]= 0.884215;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 123377;   wouttag[6]= 16345.4;  
           fbb[6]= 0.0667573;   fcc[6]= 0.127449;   fc[6]= 0.145427;   fll[6]= 0.660367;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.20737;   Fcc[7]= 1.20737;   Fc[7]= 0.950681;   Fll[7]= 0.94821; 
             cafac[7]= 0.905076;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 27257.5;   wouttag[7]= 4667.81;  
           fbb[7]= 0.0909332;   fcc[7]= 0.148714;   fc[7]= 0.134831;   fll[7]= 0.625521;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==83   || sysname=="chmis_up" ) { 
          Fbb[1]= 1.25694;   Fcc[1]= 1.25694;   Fc[1]= 0.989716;   Fll[1]= 0.987144; 
             cafac[1]= 1.05835;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.30673e+06;   wouttag[1]= 84172.4;  
           fbb[1]= 0.0148818;   fcc[1]= 0.0434135;   fc[1]= 0.131966;   fll[1]= 0.809738;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.97844;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 484024;   wouttag[2]= 38092.2;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2232;   Fcc[3]= 1.2232;   Fc[3]= 0.963152;   Fll[3]= 0.960649; 
             cafac[3]= 0.892955;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97704.3;   wouttag[3]= 11888;  
           fbb[3]= 0.0600158;   fcc[3]= 0.121519;   fc[3]= 0.148382;   fll[3]= 0.670083;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21025;   Fcc[4]= 1.21025;   Fc[4]= 0.952949;   Fll[4]= 0.950473; 
             cafac[4]= 0.957663;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22872.7;   wouttag[4]= 3686.71;  
           fbb[4]= 0.0853226;   fcc[4]= 0.143756;   fc[4]= 0.137189;   fll[4]= 0.633732;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19645;   Fcc[5]= 1.19645;   Fc[5]= 0.942087;   Fll[5]= 0.939639; 
             cafac[5]= 0.910736;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5676.08;   wouttag[5]= 1188.46;  
           fbb[5]= 0.112189;   fcc[5]= 0.167499;   fc[5]= 0.125899;   fll[5]= 0.594413;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.21975;   Fcc[6]= 1.21975;   Fc[6]= 0.960432;   Fll[6]= 0.957937; 
             cafac[6]= 0.906502;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 126487;   wouttag[6]= 16757.4;  
           fbb[6]= 0.0667573;   fcc[6]= 0.127449;   fc[6]= 0.145427;   fll[6]= 0.660367;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.20737;   Fcc[7]= 1.20737;   Fc[7]= 0.950681;   Fll[7]= 0.94821; 
             cafac[7]= 0.947987;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28549.9;   wouttag[7]= 4889.12;  
           fbb[7]= 0.0909332;   fcc[7]= 0.148714;   fc[7]= 0.134831;   fll[7]= 0.625521;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==84   || sysname=="chmis_down" ) { 
          Fbb[1]= 1.25694;   Fcc[1]= 1.25694;   Fc[1]= 0.989716;   Fll[1]= 0.987144; 
             cafac[1]= 1.04572;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.27921e+06;   wouttag[1]= 83168.4;  
           fbb[1]= 0.0148818;   fcc[1]= 0.0434135;   fc[1]= 0.131966;   fll[1]= 0.809738;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23833;   Fcc[2]= 1.23833;   Fc[2]= 0.975059;   Fll[2]= 0.972525; 
             cafac[2]= 0.968704;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 479208;   wouttag[2]= 37713.2;  
           fbb[2]= 0.0370378;   fcc[2]= 0.0891361;   fc[2]= 0.150872;   fll[2]= 0.722954;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2232;   Fcc[3]= 1.2232;   Fc[3]= 0.963152;   Fll[3]= 0.960649; 
             cafac[3]= 0.880541;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 96345.9;   wouttag[3]= 11722.7;  
           fbb[3]= 0.0600158;   fcc[3]= 0.121519;   fc[3]= 0.148382;   fll[3]= 0.670083;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.21025;   Fcc[4]= 1.21025;   Fc[4]= 0.952949;   Fll[4]= 0.950473; 
             cafac[4]= 0.944349;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22554.7;   wouttag[4]= 3635.45;  
           fbb[4]= 0.0853226;   fcc[4]= 0.143756;   fc[4]= 0.137189;   fll[4]= 0.633732;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.19645;   Fcc[5]= 1.19645;   Fc[5]= 0.942087;   Fll[5]= 0.939639; 
             cafac[5]= 0.89628;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5585.98;   wouttag[5]= 1169.59;  
           fbb[5]= 0.112189;   fcc[5]= 0.167499;   fc[5]= 0.125899;   fll[5]= 0.594413;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.21975;   Fcc[6]= 1.21975;   Fc[6]= 0.960432;   Fll[6]= 0.957937; 
             cafac[6]= 0.892113;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 124479;   wouttag[6]= 16491.4;  
           fbb[6]= 0.0667573;   fcc[6]= 0.127449;   fc[6]= 0.145427;   fll[6]= 0.660367;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.20737;   Fcc[7]= 1.20737;   Fc[7]= 0.950681;   Fll[7]= 0.94821; 
             cafac[7]= 0.93294;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28096.7;   wouttag[7]= 4811.52;  
           fbb[7]= 0.0909332;   fcc[7]= 0.148714;   fc[7]= 0.134831;   fll[7]= 0.625521;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==67   || sysname=="eer_up" ) { 
          Fbb[1]= 1.22738;   Fcc[1]= 1.22738;   Fc[1]= 1.01174;   Fll[1]= 0.985236; 
             cafac[1]= 1.05632;  
             winpretag[1]= 2.17956e+06 ;   wintag[1]= 77367 ;   woutpretag[1]= 2.30231e+06;   wouttag[1]= 84656;  
           fbb[1]= 0.0145317;   fcc[1]= 0.0423926;   fc[1]= 0.134903;   fll[1]= 0.808173;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.0345392;   fmcc[1]= 0.133338;   fmcll[1]= 0.820283;  
          Fbb[2]= 1.21042;   Fcc[2]= 1.21042;   Fc[2]= 0.99776;   Fll[2]= 0.971625; 
             cafac[2]= 0.978637;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 484121;   wouttag[2]= 38114.1;  
           fbb[2]= 0.0362031;   fcc[2]= 0.0871274;   fc[2]= 0.154385;   fll[2]= 0.722285;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.19714;   Fcc[3]= 1.19714;   Fc[3]= 0.986815;   Fll[3]= 0.960967; 
             cafac[3]= 0.891366;  
             winpretag[3]= 109416 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97529.8;   wouttag[3]= 11830.6;  
           fbb[3]= 0.0587374;   fcc[3]= 0.118931;   fc[3]= 0.152028;   fll[3]= 0.670304;  
           fmcbb[3]= 0.0490647;   fmccc[3]= 0.0993456;   fmcc[3]= 0.154059;   fmcll[3]= 0.697531;  
          Fbb[4]= 1.186;   Fcc[4]= 1.186;   Fc[4]= 0.977633;   Fll[4]= 0.952026; 
             cafac[4]= 0.955808;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22828.4;   wouttag[4]= 3660.43;  
           fbb[4]= 0.0836135;   fcc[4]= 0.140877;   fc[4]= 0.140742;   fll[4]= 0.634768;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.17409;   Fcc[5]= 1.17409;   Fc[5]= 0.967816;   Fll[5]= 0.942466; 
             cafac[5]= 0.908274;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5660.74;   wouttag[5]= 1178.18;  
           fbb[5]= 0.110093;   fcc[5]= 0.164369;   fc[5]= 0.129337;   fll[5]= 0.596201;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.19418;   Fcc[6]= 1.19418;   Fc[6]= 0.984369;   Fll[6]= 0.958585; 
             cafac[6]= 0.904004;  
             winpretag[6]= 139532 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 126138;   wouttag[6]= 16647.2;  
           fbb[6]= 0.0653579;   fcc[6]= 0.124777;   fc[6]= 0.149052;   fll[6]= 0.660813;  
           fmcbb[6]= 0.0547306;   fmccc[6]= 0.104488;   fmcc[6]= 0.151419;   fmcll[6]= 0.689362;  
          Fbb[7]= 1.18352;   Fcc[7]= 1.18352;   Fc[7]= 0.975585;   Fll[7]= 0.950031; 
             cafac[7]= 0.945275;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28468.2;   wouttag[7]= 4848.53;  
           fbb[7]= 0.0891372;   fcc[7]= 0.145777;   fc[7]= 0.138363;   fll[7]= 0.626722;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==68   || sysname=="eer_down" ) { 
          Fbb[1]= 1.22734;   Fcc[1]= 1.22734;   Fc[1]= 1.0116;   Fll[1]= 0.98526; 
             cafac[1]= 1.0563;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.30225e+06;   wouttag[1]= 84646.9;  
           fbb[1]= 0.0145314;   fcc[1]= 0.0423914;   fc[1]= 0.134885;   fll[1]= 0.808193;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.0345391;   fmcc[1]= 0.133337;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.2104;   Fcc[2]= 1.2104;   Fc[2]= 0.997634;   Fll[2]= 0.971655; 
             cafac[2]= 0.978595;  
             winpretag[2]= 494691 ;   wintag[2]= 36747 ;   woutpretag[2]= 484102;   wouttag[2]= 38110;  
           fbb[2]= 0.0362013;   fcc[2]= 0.0871255;   fc[2]= 0.154366;   fll[2]= 0.722307;  
           fmcbb[2]= 0.0299086;   fmccc[2]= 0.0719809;   fmcc[2]= 0.154732;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.19712;   Fcc[3]= 1.19712;   Fc[3]= 0.986693;   Fll[3]= 0.960999; 
             cafac[3]= 0.891368;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97530.5;   wouttag[3]= 11830.2;  
           fbb[3]= 0.0587361;   fcc[3]= 0.118928;   fc[3]= 0.152008;   fll[3]= 0.670327;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.18598;   Fcc[4]= 1.18598;   Fc[4]= 0.977509;   Fll[4]= 0.952053; 
             cafac[4]= 0.955868;  
             winpretag[4]= 23884.4 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22830.4;   wouttag[4]= 3660.52;  
           fbb[4]= 0.0836099;   fcc[4]= 0.140897;   fc[4]= 0.140721;   fll[4]= 0.634772;  
           fmcbb[4]= 0.0704986;   fmccc[4]= 0.118802;   fmcc[4]= 0.143959;   fmcll[4]= 0.66674;  
          Fbb[5]= 1.17408;   Fcc[5]= 1.17408;   Fc[5]= 0.967698;   Fll[5]= 0.942498; 
             cafac[5]= 0.908354;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5661.24;   wouttag[5]= 1178.25;  
           fbb[5]= 0.110091;   fcc[5]= 0.164367;   fc[5]= 0.129321;   fll[5]= 0.596221;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.19415;   Fcc[6]= 1.19415;   Fc[6]= 0.984247;   Fll[6]= 0.958616; 
             cafac[6]= 0.90402;  
             winpretag[6]= 139534 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 126141;   wouttag[6]= 16646.9;  
           fbb[6]= 0.0653563;   fcc[6]= 0.124779;   fc[6]= 0.149032;   fll[6]= 0.660833;  
           fmcbb[6]= 0.0547301;   fmccc[6]= 0.104491;   fmcc[6]= 0.151418;   fmcll[6]= 0.689361;  
          Fbb[7]= 1.1835;   Fcc[7]= 1.1835;   Fc[7]= 0.975462;   Fll[7]= 0.95006; 
             cafac[7]= 0.945339;  
             winpretag[7]= 30116.8 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28470.6;   wouttag[7]= 4848.69;  
           fbb[7]= 0.0891339;   fcc[7]= 0.145793;   fc[7]= 0.138343;   fll[7]= 0.62673;  
           fmcbb[7]= 0.0753141;   fmccc[7]= 0.123188;   fmcc[7]= 0.141823;   fmcll[7]= 0.659674;  
     }  
     else if ( idsys==65   || sysname=="ees_up" ) { 
          Fbb[1]= 1.22745;   Fcc[1]= 1.22745;   Fc[1]= 1.01188;   Fll[1]= 0.98521; 
             cafac[1]= 1.05636;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.30239e+06;   wouttag[1]= 84664.9;  
           fbb[1]= 0.0145326;   fcc[1]= 0.0423949;   fc[1]= 0.134921;   fll[1]= 0.808152;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133337;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.21048;   Fcc[2]= 1.21048;   Fc[2]= 0.997889;   Fll[2]= 0.971591; 
             cafac[2]= 0.978653;  
             winpretag[2]= 494692 ;   wintag[2]= 36747 ;   woutpretag[2]= 484131;   wouttag[2]= 38116.5;  
           fbb[2]= 0.0362037;   fcc[2]= 0.0871314;   fc[2]= 0.154404;   fll[2]= 0.722261;  
           fmcbb[2]= 0.0299086;   fmccc[2]= 0.0719808;   fmcc[2]= 0.154731;   fmcll[2]= 0.74338;  
          Fbb[3]= 1.1972;   Fcc[3]= 1.1972;   Fc[3]= 0.986937;   Fll[3]= 0.960928; 
             cafac[3]= 0.891477;  
             winpretag[3]= 109419 ;   wintag[3]= 12397.4 ;   woutpretag[3]= 97544.9;   wouttag[3]= 11834.1;  
           fbb[3]= 0.0587382;   fcc[3]= 0.11894;   fc[3]= 0.15205;   fll[3]= 0.670271;  
           fmcbb[3]= 0.0490632;   fmccc[3]= 0.0993492;   fmcc[3]= 0.154062;   fmcll[3]= 0.697525;  
          Fbb[4]= 1.18605;   Fcc[4]= 1.18605;   Fc[4]= 0.977747;   Fll[4]= 0.95198; 
             cafac[4]= 0.955955;  
             winpretag[4]= 23883.8 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22831.9;   wouttag[4]= 3661.13;  
           fbb[4]= 0.0836168;   fcc[4]= 0.140909;   fc[4]= 0.140759;   fll[4]= 0.634716;  
           fmcbb[4]= 0.0705004;   fmccc[4]= 0.118805;   fmcc[4]= 0.143963;   fmcll[4]= 0.666732;  
          Fbb[5]= 1.17414;   Fcc[5]= 1.17414;   Fc[5]= 0.967932;   Fll[5]= 0.942424; 
             cafac[5]= 0.908519;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5662.26;   wouttag[5]= 1178.55;  
           fbb[5]= 0.110097;   fcc[5]= 0.164376;   fc[5]= 0.129353;   fll[5]= 0.596174;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.19423;   Fcc[6]= 1.19423;   Fc[6]= 0.98449;   Fll[6]= 0.958545; 
             cafac[6]= 0.904128;  
             winpretag[6]= 139536 ;   wintag[6]= 17147.4 ;   woutpretag[6]= 126158;   wouttag[6]= 16651.8;  
           fbb[6]= 0.0653592;   fcc[6]= 0.124791;   fc[6]= 0.149073;   fll[6]= 0.660777;  
           fmcbb[6]= 0.0547293;   fmccc[6]= 0.104495;   fmcc[6]= 0.151421;   fmcll[6]= 0.689354;  
          Fbb[7]= 1.18356;   Fcc[7]= 1.18356;   Fc[7]= 0.9757;   Fll[7]= 0.949987; 
             cafac[7]= 0.945444;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28473.2;   wouttag[7]= 4849.59;  
           fbb[7]= 0.0891408;   fcc[7]= 0.145804;   fc[7]= 0.13838;   fll[7]= 0.626676;  
           fmcbb[7]= 0.0753156;   fmccc[7]= 0.123191;   fmcc[7]= 0.141826;   fmcll[7]= 0.659668;  
     }  
     else if ( idsys==66   || sysname=="ees_down" ) { 
          Fbb[1]= 1.22727;   Fcc[1]= 1.22727;   Fc[1]= 1.01152;   Fll[1]= 0.985277; 
             cafac[1]= 1.05628;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77364.9 ;   woutpretag[1]= 2.3022e+06;   wouttag[1]= 84639.6;  
           fbb[1]= 0.0145302;   fcc[1]= 0.042389;   fc[1]= 0.134873;   fll[1]= 0.808207;  
           fmcbb[1]= 0.0118395;   fmccc[1]= 0.0345394;   fmcc[1]= 0.133337;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.21033;   Fcc[2]= 1.21033;   Fc[2]= 0.997564;   Fll[2]= 0.971679; 
             cafac[2]= 0.978576;  
             winpretag[2]= 494692 ;   wintag[2]= 36747 ;   woutpretag[2]= 484094;   wouttag[2]= 38107.6;  
           fbb[2]= 0.0361992;   fcc[2]= 0.0871219;   fc[2]= 0.154355;   fll[2]= 0.722324;  
           fmcbb[2]= 0.0299086;   fmccc[2]= 0.0719821;   fmcc[2]= 0.154732;   fmcll[2]= 0.743377;  
          Fbb[3]= 1.19706;   Fcc[3]= 1.19706;   Fc[3]= 0.986627;   Fll[3]= 0.961025; 
             cafac[3]= 0.891326;  
             winpretag[3]= 109417 ;   wintag[3]= 12396.4 ;   woutpretag[3]= 97526.6;   wouttag[3]= 11829.7;  
           fbb[3]= 0.0587325;   fcc[3]= 0.118929;   fc[3]= 0.151997;   fll[3]= 0.670342;  
           fmcbb[3]= 0.0490641;   fmccc[3]= 0.099351;   fmcc[3]= 0.154057;   fmcll[3]= 0.697528;  
          Fbb[4]= 1.18592;   Fcc[4]= 1.18592;   Fc[4]= 0.977447;   Fll[4]= 0.952083; 
             cafac[4]= 0.955909;  
             winpretag[4]= 23884.4 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22831.4;   wouttag[4]= 3660.55;  
           fbb[4]= 0.0836057;   fcc[4]= 0.14089;   fc[4]= 0.140712;   fll[4]= 0.634792;  
           fmcbb[4]= 0.0704986;   fmccc[4]= 0.118802;   fmcc[4]= 0.143959;   fmcll[4]= 0.66674;  
          Fbb[5]= 1.17402;   Fcc[5]= 1.17402;   Fc[5]= 0.967639;   Fll[5]= 0.94253; 
             cafac[5]= 0.908369;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5661.32;   wouttag[5]= 1178.23;  
           fbb[5]= 0.110086;   fcc[5]= 0.164359;   fc[5]= 0.129314;   fll[5]= 0.596242;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.19409;   Fcc[6]= 1.19409;   Fc[6]= 0.984182;   Fll[6]= 0.958644; 
             cafac[6]= 0.903997;  
             winpretag[6]= 139534 ;   wintag[6]= 17146.4 ;   woutpretag[6]= 126139;   wouttag[6]= 16646.3;  
           fbb[6]= 0.0653524;   fcc[6]= 0.124778;   fc[6]= 0.149022;   fll[6]= 0.660848;  
           fmcbb[6]= 0.0547299;   fmccc[6]= 0.104496;   fmcc[6]= 0.151417;   fmcll[6]= 0.689357;  
          Fbb[7]= 1.18344;   Fcc[7]= 1.18344;   Fc[7]= 0.975401;   Fll[7]= 0.950091; 
             cafac[7]= 0.945374;  
             winpretag[7]= 30116.8 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28471.7;   wouttag[7]= 4848.69;  
           fbb[7]= 0.0891295;   fcc[7]= 0.145785;   fc[7]= 0.138335;   fll[7]= 0.626751;  
           fmcbb[7]= 0.0753141;   fmccc[7]= 0.123188;   fmcc[7]= 0.141823;   fmcll[7]= 0.659674;  
     }  
     else if ( idsys==63   || sysname=="jer_one" ) { 
          Fbb[1]= 1.30931;   Fcc[1]= 1.30931;   Fc[1]= 0.94837;   Fll[1]= 0.991443; 
             cafac[1]= 0.985712;  
             winpretag[1]= 2.27557e+06 ;   wintag[1]= 78950.9 ;   woutpretag[1]= 2.24305e+06;   wouttag[1]= 78758.2;  
           fbb[1]= 0.0148219;   fcc[1]= 0.0433548;   fc[1]= 0.122561;   fll[1]= 0.819263;  
           fmcbb[1]= 0.0113204;   fmccc[1]= 0.0331128;   fmcc[1]= 0.129233;   fmcll[1]= 0.826334;  
          Fbb[2]= 1.28836;   Fcc[2]= 1.28836;   Fc[2]= 0.933196;   Fll[2]= 0.975579; 
             cafac[2]= 0.939908;  
             winpretag[2]= 512839 ;   wintag[2]= 37957.8 ;   woutpretag[2]= 482021;   wouttag[2]= 37656.7;  
           fbb[2]= 0.0373835;   fcc[2]= 0.0902103;   fc[2]= 0.14434;   fll[2]= 0.728066;  
           fmcbb[2]= 0.0290164;   fmccc[2]= 0.0700196;   fmcc[2]= 0.154673;   fmcll[2]= 0.746291;  
          Fbb[3]= 1.27019;   Fcc[3]= 1.27019;   Fc[3]= 0.920035;   Fll[3]= 0.96182; 
             cafac[3]= 0.852751;  
             winpretag[3]= 115457 ;   wintag[3]= 12750.1 ;   woutpretag[3]= 98456.4;   wouttag[3]= 11722.6;  
           fbb[3]= 0.0616053;   fcc[3]= 0.122461;   fc[3]= 0.143259;   fll[3]= 0.672674;  
           fmcbb[3]= 0.048501;   fmccc[3]= 0.0964121;   fmcc[3]= 0.155711;   fmcll[3]= 0.699376;  
          Fbb[4]= 1.25402;   Fcc[4]= 1.25402;   Fc[4]= 0.908325;   Fll[4]= 0.949579; 
             cafac[4]= 0.876206;  
             winpretag[4]= 25732 ;   wintag[4]= 3792.04 ;   woutpretag[4]= 22546.5;   wouttag[4]= 3627.91;  
           fbb[4]= 0.0841441;   fcc[4]= 0.148316;   fc[4]= 0.132406;   fll[4]= 0.635134;  
           fmcbb[4]= 0.0670994;   fmccc[4]= 0.118272;   fmcc[4]= 0.145769;   fmcll[4]= 0.668859;  
          Fbb[5]= 1.23657;   Fcc[5]= 1.23657;   Fc[5]= 0.895683;   Fll[5]= 0.936362; 
             cafac[5]= 0.792337;  
             winpretag[5]= 6899.23 ;   wintag[5]= 1275.91 ;   woutpretag[5]= 5466.52;   wouttag[5]= 1109.49;  
           fbb[5]= 0.112551;   fcc[5]= 0.171705;   fc[5]= 0.118274;   fll[5]= 0.59747;  
           fmcbb[5]= 0.0910189;   fmccc[5]= 0.138856;   fmcc[5]= 0.132049;   fmcll[5]= 0.638075;  
          Fbb[6]= 1.26575;   Fcc[6]= 1.26575;   Fc[6]= 0.91682;   Fll[6]= 0.958459; 
             cafac[6]= 0.853429;  
             winpretag[6]= 148089 ;   wintag[6]= 17818.1 ;   woutpretag[6]= 126383;   wouttag[6]= 16480.8;  
           fbb[6]= 0.0679878;   fcc[6]= 0.129344;   fc[6]= 0.140164;   fll[6]= 0.662504;  
           fmcbb[6]= 0.0537135;   fmccc[6]= 0.102188;   fmcc[6]= 0.152881;   fmcll[6]= 0.691218;  
          Fbb[7]= 1.25029;   Fcc[7]= 1.25029;   Fc[7]= 0.905622;   Fll[7]= 0.946753; 
             cafac[7]= 0.856945;  
             winpretag[7]= 32631.2 ;   wintag[7]= 5067.95 ;   woutpretag[7]= 27963.2;   wouttag[7]= 4750.88;  
           fbb[7]= 0.0902168;   fcc[7]= 0.153316;   fc[7]= 0.129385;   fll[7]= 0.627083;  
           fmcbb[7]= 0.0721567;   fmccc[7]= 0.122624;   fmcc[7]= 0.142868;   fmcll[7]= 0.662351;  
     }  
     else if ( idsys==64   || sysname=="jef_one" ) { 
          Fbb[1]= 1.27613;   Fcc[1]= 1.27613;   Fc[1]= 0.975122;   Fll[1]= 0.988426; 
             cafac[1]= 1.04987;  
             winpretag[1]= 2.17846e+06 ;   wintag[1]= 77310.8 ;   woutpretag[1]= 2.2871e+06;   wouttag[1]= 83011.1;  
           fbb[1]= 0.0151118;   fcc[1]= 0.0440934;   fc[1]= 0.130019;   fll[1]= 0.810775;  
           fmcbb[1]= 0.0118419;   fmccc[1]= 0.0345525;   fmcc[1]= 0.133337;   fmcll[1]= 0.820269;  
          Fbb[2]= 1.25642;   Fcc[2]= 1.25642;   Fc[2]= 0.960062;   Fll[2]= 0.97316; 
             cafac[2]= 0.971863;  
             winpretag[2]= 494285 ;   wintag[2]= 36723.7 ;   woutpretag[2]= 480378;   wouttag[2]= 37800;  
           fbb[2]= 0.0376035;   fcc[2]= 0.090441;   fc[2]= 0.148605;   fll[2]= 0.72335;  
           fmcbb[2]= 0.0299292;   fmccc[2]= 0.0719832;   fmcc[2]= 0.154787;   fmcll[2]= 0.7433;  
          Fbb[3]= 1.24008;   Fcc[3]= 1.24008;   Fc[3]= 0.947577;   Fll[3]= 0.960505; 
             cafac[3]= 0.884199;  
             winpretag[3]= 109277 ;   wintag[3]= 12378.2 ;   woutpretag[3]= 96622.7;   wouttag[3]= 11776.8;  
           fbb[3]= 0.0608315;   fcc[3]= 0.123185;   fc[3]= 0.14594;   fll[3]= 0.670043;  
           fmcbb[3]= 0.0490546;   fmccc[3]= 0.0993368;   fmcc[3]= 0.154013;   fmcll[3]= 0.697595;  
          Fbb[4]= 1.22586;   Fcc[4]= 1.22586;   Fc[4]= 0.936715;   Fll[4]= 0.949495; 
             cafac[4]= 0.947992;  
             winpretag[4]= 23858.1 ;   wintag[4]= 3550.38 ;   woutpretag[4]= 22617.3;   wouttag[4]= 3662.91;  
           fbb[4]= 0.0865892;   fcc[4]= 0.145609;   fc[4]= 0.135125;   fll[4]= 0.632677;  
           fmcbb[4]= 0.0706352;   fmccc[4]= 0.11878;   fmcc[4]= 0.144254;   fmcll[4]= 0.66633;  
          Fbb[5]= 1.21089;   Fcc[5]= 1.21089;   Fc[5]= 0.925277;   Fll[5]= 0.9379; 
             cafac[5]= 0.90345;  
             winpretag[5]= 6220.52 ;   wintag[5]= 1196.14 ;   woutpretag[5]= 5619.93;   wouttag[5]= 1178.43;  
           fbb[5]= 0.113126;   fcc[5]= 0.169801;   fc[5]= 0.123554;   fll[5]= 0.593519;  
           fmcbb[5]= 0.0934233;   fmccc[5]= 0.140228;   fmcc[5]= 0.133532;   fmcll[5]= 0.632817;  
          Fbb[6]= 1.23629;   Fcc[6]= 1.23629;   Fc[6]= 0.944685;   Fll[6]= 0.957574; 
             cafac[6]= 0.896784;  
             winpretag[6]= 139356 ;   wintag[6]= 17124.7 ;   woutpretag[6]= 124972;   wouttag[6]= 16596.8;  
           fbb[6]= 0.0676621;   fcc[6]= 0.129182;   fc[6]= 0.143052;   fll[6]= 0.660104;  
           fmcbb[6]= 0.0547297;   fmccc[6]= 0.104491;   fmcc[6]= 0.151428;   fmcll[6]= 0.689351;  
          Fbb[7]= 1.22274;   Fcc[7]= 1.22274;   Fc[7]= 0.934326;   Fll[7]= 0.947073; 
             cafac[7]= 0.938113;  
             winpretag[7]= 30078.7 ;   wintag[7]= 4746.52 ;   woutpretag[7]= 28217.2;   wouttag[7]= 4851.1;  
           fbb[7]= 0.0921308;   fcc[7]= 0.150661;   fc[7]= 0.132709;   fll[7]= 0.624499;  
           fmcbb[7]= 0.075348;   fmccc[7]= 0.123216;   fmcc[7]= 0.142037;   fmcll[7]= 0.659399;  
     }  
     else if ( idsys==61   || sysname=="musms_up" ) { 
          Fbb[1]= 1.27113;   Fcc[1]= 1.27113;   Fc[1]= 0.977042;   Fll[1]= 0.988403; 
             cafac[1]= 1.04927;  
             winpretag[1]= 2.17952e+06 ;   wintag[1]= 77387.8 ;   woutpretag[1]= 2.28691e+06;   wouttag[1]= 83074.8;  
           fbb[1]= 0.0150471;   fcc[1]= 0.0439061;   fc[1]= 0.130302;   fll[1]= 0.810745;  
           fmcbb[1]= 0.0118376;   fmccc[1]= 0.034541;   fmcc[1]= 0.133363;   fmcll[1]= 0.820258;  
          Fbb[2]= 1.2518;   Fcc[2]= 1.2518;   Fc[2]= 0.962183;   Fll[2]= 0.973371; 
             cafac[2]= 0.971365;  
             winpretag[2]= 494788 ;   wintag[2]= 36767.1 ;   woutpretag[2]= 480620;   wouttag[2]= 37803.9;  
           fbb[2]= 0.0374362;   fcc[2]= 0.0900703;   fc[2]= 0.148922;   fll[2]= 0.723572;  
           fmcbb[2]= 0.0299058;   fmccc[2]= 0.0719526;   fmcc[2]= 0.154775;   fmcll[2]= 0.743367;  
          Fbb[3]= 1.23581;   Fcc[3]= 1.23581;   Fc[3]= 0.949889;   Fll[3]= 0.960934; 
             cafac[3]= 0.883774;  
             winpretag[3]= 109446 ;   wintag[3]= 12402.2 ;   woutpretag[3]= 96726;   wouttag[3]= 11783;  
           fbb[3]= 0.0606056;   fcc[3]= 0.122681;   fc[3]= 0.146309;   fll[3]= 0.670404;  
           fmcbb[3]= 0.0490414;   fmccc[3]= 0.0992721;   fmcc[3]= 0.154027;   fmcll[3]= 0.697659;  
          Fbb[4]= 1.22191;   Fcc[4]= 1.22191;   Fc[4]= 0.939208;   Fll[4]= 0.950128; 
             cafac[4]= 0.948255;  
             winpretag[4]= 23895.3 ;   wintag[4]= 3549.37 ;   woutpretag[4]= 22658.9;   wouttag[4]= 3659.26;  
           fbb[4]= 0.0861879;   fcc[4]= 0.145103;   fc[4]= 0.13526;   fll[4]= 0.633449;  
           fmcbb[4]= 0.0705354;   fmccc[4]= 0.118751;   fmcc[4]= 0.144015;   fmcll[4]= 0.666699;  
          Fbb[5]= 1.20711;   Fcc[5]= 1.20711;   Fc[5]= 0.92783;   Fll[5]= 0.938619; 
             cafac[5]= 0.9026;  
             winpretag[5]= 6234.47 ;   wintag[5]= 1200.38 ;   woutpretag[5]= 5627.23;   wouttag[5]= 1180.03;  
           fbb[5]= 0.113243;   fcc[5]= 0.169203;   fc[5]= 0.123957;   fll[5]= 0.593597;  
           fmcbb[5]= 0.0938132;   fmccc[5]= 0.140172;   fmcc[5]= 0.133599;   fmcll[5]= 0.632415;  
          Fbb[6]= 1.2321;   Fcc[6]= 1.2321;   Fc[6]= 0.94704;   Fll[6]= 0.958051; 
             cafac[6]= 0.896468;  
             winpretag[6]= 139576 ;   wintag[6]= 17151.9 ;   woutpretag[6]= 125126;   wouttag[6]= 16600.9;  
           fbb[6]= 0.0674216;   fcc[6]= 0.128673;   fc[6]= 0.143383;   fll[6]= 0.660523;  
           fmcbb[6]= 0.054721;   fmccc[6]= 0.104434;   fmcc[6]= 0.151401;   fmcll[6]= 0.689444;  
          Fbb[7]= 1.21882;   Fcc[7]= 1.21882;   Fc[7]= 0.936831;   Fll[7]= 0.947724; 
             cafac[7]= 0.938123;  
             winpretag[7]= 30129.8 ;   wintag[7]= 4749.75 ;   woutpretag[7]= 28265.5;   wouttag[7]= 4849.37;  
           fbb[7]= 0.0918403;   fcc[7]= 0.150138;   fc[7]= 0.132899;   fll[7]= 0.625123;  
           fmcbb[7]= 0.0753521;   fmccc[7]= 0.123183;   fmcc[7]= 0.14186;   fmcll[7]= 0.659605;  
     }  
     else if ( idsys==62   || sysname=="musms_down" ) { 
          Fbb[1]= 1.27373;   Fcc[1]= 1.27373;   Fc[1]= 0.979292;   Fll[1]= 0.987891; 
             cafac[1]= 1.04992;  
             winpretag[1]= 2.17928e+06 ;   wintag[1]= 77355.9 ;   woutpretag[1]= 2.28808e+06;   wouttag[1]= 83218;  
           fbb[1]= 0.0150766;   fcc[1]= 0.0439922;   fc[1]= 0.130561;   fll[1]= 0.81037;  
           fmcbb[1]= 0.0118365;   fmccc[1]= 0.034538;   fmcc[1]= 0.133322;   fmcll[1]= 0.820304;  
          Fbb[2]= 1.25405;   Fcc[2]= 1.25405;   Fc[2]= 0.964159;   Fll[2]= 0.972625; 
             cafac[2]= 0.970898;  
             winpretag[2]= 494557 ;   wintag[2]= 36745.7 ;   woutpretag[2]= 480165;   wouttag[2]= 37815.6;  
           fbb[2]= 0.0375418;   fcc[2]= 0.0902792;   fc[2]= 0.149179;   fll[2]= 0.723;  
           fmcbb[2]= 0.0299364;   fmccc[2]= 0.07199;   fmcc[2]= 0.154725;   fmcll[2]= 0.743349;  
          Fbb[3]= 1.23787;   Fcc[3]= 1.23787;   Fc[3]= 0.951722;   Fll[3]= 0.960078; 
             cafac[3]= 0.884143;  
             winpretag[3]= 109406 ;   wintag[3]= 12386.6 ;   woutpretag[3]= 96730.1;   wouttag[3]= 11788;  
           fbb[3]= 0.0607292;   fcc[3]= 0.122904;   fc[3]= 0.14667;   fll[3]= 0.669697;  
           fmcbb[3]= 0.0490593;   fmccc[3]= 0.0992861;   fmcc[3]= 0.15411;   fmcll[3]= 0.697545;  
          Fbb[4]= 1.22386;   Fcc[4]= 1.22386;   Fc[4]= 0.940945;   Fll[4]= 0.949207; 
             cafac[4]= 0.947345;  
             winpretag[4]= 23882.6 ;   wintag[4]= 3544 ;   woutpretag[4]= 22625.1;   wouttag[4]= 3653.89;  
           fbb[4]= 0.0862822;   fcc[4]= 0.145355;   fc[4]= 0.135423;   fll[4]= 0.63294;  
           fmcbb[4]= 0.0705003;   fmccc[4]= 0.118768;   fmcc[4]= 0.143922;   fmcll[4]= 0.66681;  
          Fbb[5]= 1.2089;   Fcc[5]= 1.2089;   Fc[5]= 0.929442;   Fll[5]= 0.937603; 
             cafac[5]= 0.900454;  
             winpretag[5]= 6229.76 ;   wintag[5]= 1200.45 ;   woutpretag[5]= 5609.61;   wouttag[5]= 1178.52;  
           fbb[5]= 0.113496;   fcc[5]= 0.169402;   fc[5]= 0.12403;   fll[5]= 0.593072;  
           fmcbb[5]= 0.0938844;   fmccc[5]= 0.14013;   fmcc[5]= 0.133445;   fmcll[5]= 0.632541;  
          Fbb[6]= 1.23413;   Fcc[6]= 1.23413;   Fc[6]= 0.948846;   Fll[6]= 0.957177; 
             cafac[6]= 0.896464;  
             winpretag[6]= 139518 ;   wintag[6]= 17131 ;   woutpretag[6]= 125073;   wouttag[6]= 16600.6;  
           fbb[6]= 0.0675454;   fcc[6]= 0.128899;   fc[6]= 0.143696;   fll[6]= 0.659859;  
           fmcbb[6]= 0.0547311;   fmccc[6]= 0.104445;   fmcc[6]= 0.151443;   fmcll[6]= 0.689381;  
          Fbb[7]= 1.22073;   Fcc[7]= 1.22073;   Fc[7]= 0.938542;   Fll[7]= 0.946782; 
             cafac[7]= 0.936915;  
             winpretag[7]= 30112.4 ;   wintag[7]= 4744.45 ;   woutpretag[7]= 28212.7;   wouttag[7]= 4842.69;  
           fbb[7]= 0.0919675;   fcc[7]= 0.150378;   fc[7]= 0.133043;   fll[7]= 0.624611;  
           fmcbb[7]= 0.0753381;   fmccc[7]= 0.123187;   fmcc[7]= 0.141755;   fmcll[7]= 0.65972;  
     }  
     else if ( idsys==59   || sysname=="musid_up" ) { 
          Fbb[1]= 1.27226;   Fcc[1]= 1.27226;   Fc[1]= 0.979936;   Fll[1]= 0.987866; 
             cafac[1]= 1.05009;  
             winpretag[1]= 2.17946e+06 ;   wintag[1]= 77367.8 ;   woutpretag[1]= 2.28864e+06;   wouttag[1]= 83257.4;  
           fbb[1]= 0.0150644;   fcc[1]= 0.0439472;   fc[1]= 0.130664;   fll[1]= 0.810324;  
           fmcbb[1]= 0.0118407;   fmccc[1]= 0.0345426;   fmcc[1]= 0.133339;   fmcll[1]= 0.820277;  
          Fbb[2]= 1.2527;   Fcc[2]= 1.2527;   Fc[2]= 0.964869;   Fll[2]= 0.972677; 
             cafac[2]= 0.971096;  
             winpretag[2]= 494664 ;   wintag[2]= 36741.8 ;   woutpretag[2]= 480366;   wouttag[2]= 37813.1;  
           fbb[2]= 0.0374677;   fcc[2]= 0.0901693;   fc[2]= 0.14929;   fll[2]= 0.723073;  
           fmcbb[2]= 0.0299096;   fmccc[2]= 0.07198;   fmcc[2]= 0.154725;   fmcll[2]= 0.743385;  
          Fbb[3]= 1.23659;   Fcc[3]= 1.23659;   Fc[3]= 0.952463;   Fll[3]= 0.960171; 
             cafac[3]= 0.884363;  
             winpretag[3]= 109418 ;   wintag[3]= 12397.3 ;   woutpretag[3]= 96765.4;   wouttag[3]= 11797.6;  
           fbb[3]= 0.0606719;   fcc[3]= 0.122819;   fc[3]= 0.146751;   fll[3]= 0.669757;  
           fmcbb[3]= 0.0490638;   fmccc[3]= 0.0993207;   fmcc[3]= 0.154076;   fmcll[3]= 0.69754;  
          Fbb[4]= 1.22268;   Fcc[4]= 1.22268;   Fc[4]= 0.941748;   Fll[4]= 0.949369; 
             cafac[4]= 0.947918;  
             winpretag[4]= 23886.2 ;   wintag[4]= 3548.44 ;   woutpretag[4]= 22642.1;   wouttag[4]= 3659.56;  
           fbb[4]= 0.086159;   fcc[4]= 0.145253;   fc[4]= 0.135594;   fll[4]= 0.632994;  
           fmcbb[4]= 0.0704673;   fmccc[4]= 0.118799;   fmcc[4]= 0.143982;   fmcll[4]= 0.666752;  
          Fbb[5]= 1.20789;   Fcc[5]= 1.20789;   Fc[5]= 0.930359;   Fll[5]= 0.937888; 
             cafac[5]= 0.901586;  
             winpretag[5]= 6232.4 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5619.05;   wouttag[5]= 1180.38;  
           fbb[5]= 0.113262;   fcc[5]= 0.169101;   fc[5]= 0.124332;   fll[5]= 0.593305;  
           fmcbb[5]= 0.0937684;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.23288;   Fcc[6]= 1.23288;   Fc[6]= 0.949606;   Fll[6]= 0.95729; 
             cafac[6]= 0.896802;  
             winpretag[6]= 139537 ;   wintag[6]= 17147 ;   woutpretag[6]= 125137;   wouttag[6]= 16617.1;  
           fbb[6]= 0.0674688;   fcc[6]= 0.128801;   fc[6]= 0.143803;   fll[6]= 0.659926;  
           fmcbb[6]= 0.0547244;   fmccc[6]= 0.104472;   fmcc[6]= 0.151435;   fmcll[6]= 0.689369;  
          Fbb[7]= 1.21959;   Fcc[7]= 1.21959;   Fc[7]= 0.939368;   Fll[7]= 0.94697; 
             cafac[7]= 0.937632;  
             winpretag[7]= 30118.6 ;   wintag[7]= 4749.72 ;   woutpretag[7]= 28240.1;   wouttag[7]= 4850.15;  
           fbb[7]= 0.0918218;   fcc[7]= 0.150236;   fc[7]= 0.133241;   fll[7]= 0.624701;  
           fmcbb[7]= 0.075289;   fmccc[7]= 0.123185;   fmcc[7]= 0.141841;   fmcll[7]= 0.659685;  
     }  
     else if ( idsys==60   || sysname=="musid_down" ) { 
          Fbb[1]= 1.27202;   Fcc[1]= 1.27202;   Fc[1]= 0.979613;   Fll[1]= 0.987934; 
             cafac[1]= 1.04985;  
             winpretag[1]= 2.17948e+06 ;   wintag[1]= 77360 ;   woutpretag[1]= 2.28813e+06;   wouttag[1]= 83212.5;  
           fbb[1]= 0.0150596;   fcc[1]= 0.0439349;   fc[1]= 0.13062;   fll[1]= 0.810386;  
           fmcbb[1]= 0.0118391;   fmccc[1]= 0.0345396;   fmcc[1]= 0.133338;   fmcll[1]= 0.820283;  
          Fbb[2]= 1.25249;   Fcc[2]= 1.25249;   Fc[2]= 0.964575;   Fll[2]= 0.972768; 
             cafac[2]= 0.971049;  
             winpretag[2]= 494702 ;   wintag[2]= 36749.7 ;   woutpretag[2]= 480379;   wouttag[2]= 37813.6;  
           fbb[2]= 0.037465;   fcc[2]= 0.0901451;   fc[2]= 0.149252;   fll[2]= 0.723138;  
           fmcbb[2]= 0.0299124;   fmccc[2]= 0.0719727;   fmcc[2]= 0.154733;   fmcll[2]= 0.743382;  
          Fbb[3]= 1.2364;   Fcc[3]= 1.2364;   Fc[3]= 0.952183;   Fll[3]= 0.960272; 
             cafac[3]= 0.884069;  
             winpretag[3]= 109424 ;   wintag[3]= 12395.1 ;   woutpretag[3]= 96738.7;   wouttag[3]= 11790.4;  
           fbb[3]= 0.0606583;   fcc[3]= 0.122808;   fc[3]= 0.146676;   fll[3]= 0.669857;  
           fmcbb[3]= 0.0490604;   fmccc[3]= 0.0993274;   fmcc[3]= 0.154042;   fmcll[3]= 0.69757;  
          Fbb[4]= 1.22249;   Fcc[4]= 1.22249;   Fc[4]= 0.941473;   Fll[4]= 0.949471; 
             cafac[4]= 0.948979;  
             winpretag[4]= 23882.8 ;   wintag[4]= 3547.65 ;   woutpretag[4]= 22664.3;   wouttag[4]= 3662.61;  
           fbb[4]= 0.08619;   fcc[4]= 0.145216;   fc[4]= 0.135512;   fll[4]= 0.633082;  
           fmcbb[4]= 0.0705035;   fmccc[4]= 0.118787;   fmcc[4]= 0.143936;   fmcll[4]= 0.666774;  
          Fbb[5]= 1.20773;   Fcc[5]= 1.20773;   Fc[5]= 0.930106;   Fll[5]= 0.938007; 
             cafac[5]= 0.902858;  
             winpretag[5]= 6232.63 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5627.18;   wouttag[5]= 1181.91;  
           fbb[5]= 0.113243;   fcc[5]= 0.169071;   fc[5]= 0.124413;   fll[5]= 0.593272;  
           fmcbb[5]= 0.093765;   fmccc[5]= 0.139991;   fmcc[5]= 0.133763;   fmcll[5]= 0.632482;  
          Fbb[6]= 1.23269;   Fcc[6]= 1.23269;   Fc[6]= 0.949329;   Fll[6]= 0.957393; 
             cafac[6]= 0.896834;  
             winpretag[6]= 139540 ;   wintag[6]= 17144.1 ;   woutpretag[6]= 125144;   wouttag[6]= 16613.3;  
           fbb[6]= 0.0674618;   fcc[6]= 0.128785;   fc[6]= 0.143735;   fll[6]= 0.660019;  
           fmcbb[6]= 0.0547272;   fmccc[6]= 0.104474;   fmcc[6]= 0.151407;   fmcll[6]= 0.689392;  
          Fbb[7]= 1.21941;   Fcc[7]= 1.21941;   Fc[7]= 0.939098;   Fll[7]= 0.947075; 
             cafac[7]= 0.938742;  
             winpretag[7]= 30115.4 ;   wintag[7]= 4748.94 ;   woutpretag[7]= 28270.6;   wouttag[7]= 4854.73;  
           fbb[7]= 0.0918429;   fcc[7]= 0.150201;   fc[7]= 0.133193;   fll[7]= 0.624764;  
           fmcbb[7]= 0.0753176;   fmccc[7]= 0.123175;   fmcc[7]= 0.14183;   fmcll[7]= 0.659677;  
     }  
     else if ( idsys==57   || sysname=="mscale_one" ) { 
          Fbb[1]= 1.23893;   Fcc[1]= 1.23893;   Fc[1]= 1.00116;   Fll[1]= 0.9863; 
             cafac[1]= 1.05195;  
             winpretag[1]= 2.18839e+06 ;   wintag[1]= 77706.6 ;   woutpretag[1]= 2.30207e+06;   wouttag[1]= 84336.8;  
           fbb[1]= 0.0146692;   fcc[1]= 0.0427948;   fc[1]= 0.133565;   fll[1]= 0.808971;  
           fmcbb[1]= 0.0118402;   fmccc[1]= 0.0345416;   fmcc[1]= 0.13341;   fmcll[1]= 0.820208;  
          Fbb[2]= 1.22142;   Fcc[2]= 1.22142;   Fc[2]= 0.987007;   Fll[2]= 0.972358; 
             cafac[2]= 0.97481;  
             winpretag[2]= 496446 ;   wintag[2]= 36893 ;   woutpretag[2]= 483941;   wouttag[2]= 38086.1;  
           fbb[2]= 0.0365576;   fcc[2]= 0.0878774;   fc[2]= 0.152814;   fll[2]= 0.722751;  
           fmcbb[2]= 0.0299304;   fmccc[2]= 0.0719468;   fmcc[2]= 0.154826;   fmcll[2]= 0.743297;  
          Fbb[3]= 1.20748;   Fcc[3]= 1.20748;   Fc[3]= 0.975742;   Fll[3]= 0.96126; 
             cafac[3]= 0.888133;  
             winpretag[3]= 109781 ;   wintag[3]= 12442.8 ;   woutpretag[3]= 97499.9;   wouttag[3]= 11840.6;  
           fbb[3]= 0.0592229;   fcc[3]= 0.119816;   fc[3]= 0.150375;   fll[3]= 0.670587;  
           fmcbb[3]= 0.0490467;   fmccc[3]= 0.0992281;   fmcc[3]= 0.154113;   fmcll[3]= 0.697612;  
          Fbb[4]= 1.19558;   Fcc[4]= 1.19558;   Fc[4]= 0.966128;   Fll[4]= 0.951789; 
             cafac[4]= 0.952512;  
             winpretag[4]= 23945.4 ;   wintag[4]= 3558.11 ;   woutpretag[4]= 22808.3;   wouttag[4]= 3663.52;  
           fbb[4]= 0.0842886;   fcc[4]= 0.142021;   fc[4]= 0.139053;   fll[4]= 0.634638;  
           fmcbb[4]= 0.0705;   fmccc[4]= 0.118788;   fmcc[4]= 0.143928;   fmcll[4]= 0.666784;  
          Fbb[5]= 1.18287;   Fcc[5]= 1.18287;   Fc[5]= 0.955858;   Fll[5]= 0.941672; 
             cafac[5]= 0.904038;  
             winpretag[5]= 6241.23 ;   wintag[5]= 1199.65 ;   woutpretag[5]= 5642.31;   wouttag[5]= 1173.33;  
           fbb[5]= 0.110797;   fcc[5]= 0.165923;   fc[5]= 0.128124;   fll[5]= 0.595156;  
           fmcbb[5]= 0.0936672;   fmccc[5]= 0.140271;   fmcc[5]= 0.134041;   fmcll[5]= 0.632021;  
          Fbb[6]= 1.20431;   Fcc[6]= 1.20431;   Fc[6]= 0.973182;   Fll[6]= 0.958739; 
             cafac[6]= 0.90069;  
             winpretag[6]= 139967 ;   wintag[6]= 17200.6 ;   woutpretag[6]= 126067;   wouttag[6]= 16656.1;  
           fbb[6]= 0.0658838;   fcc[6]= 0.125736;   fc[6]= 0.147413;   fll[6]= 0.660967;  
           fmcbb[6]= 0.0547065;   fmccc[6]= 0.104405;   fmcc[6]= 0.151476;   fmcll[6]= 0.689413;  
          Fbb[7]= 1.19293;   Fcc[7]= 1.19293;   Fc[7]= 0.963987;   Fll[7]= 0.949679; 
             cafac[7]= 0.941764;  
             winpretag[7]= 30186.7 ;   wintag[7]= 4757.75 ;   woutpretag[7]= 28428.7;   wouttag[7]= 4846.83;  
           fbb[7]= 0.0898158;   fcc[7]= 0.147005;   fc[7]= 0.136774;   fll[7]= 0.626405;  
           fmcbb[7]= 0.0752899;   fmccc[7]= 0.12323;   fmcc[7]= 0.141884;   fmcll[7]= 0.659597;  
     }  
     else if ( idsys==55   || sysname=="metpileup_up" ) { 
          Fbb[1]= 1.22326;   Fcc[1]= 1.22326;   Fc[1]= 1.01066;   Fll[1]= 0.985645; 
             cafac[1]= 1.05518;  
             winpretag[1]= 2.18241e+06 ;   wintag[1]= 77559.1 ;   woutpretag[1]= 2.30284e+06;   wouttag[1]= 84680.6;  
           fbb[1]= 0.0144787;   fcc[1]= 0.0422373;   fc[1]= 0.134831;   fll[1]= 0.808453;  
           fmcbb[1]= 0.0118362;   fmccc[1]= 0.0345284;   fmcc[1]= 0.133408;   fmcll[1]= 0.820228;  
          Fbb[2]= 1.20668;   Fcc[2]= 1.20668;   Fc[2]= 0.996967;   Fll[2]= 0.972287; 
             cafac[2]= 0.978432;  
             winpretag[2]= 494926 ;   wintag[2]= 36775.1 ;   woutpretag[2]= 484251;   wouttag[2]= 38087.3;  
           fbb[2]= 0.0361349;   fcc[2]= 0.0868706;   fc[2]= 0.154289;   fll[2]= 0.722705;  
           fmcbb[2]= 0.0299456;   fmccc[2]= 0.0719911;   fmcc[2]= 0.154759;   fmcll[2]= 0.743304;  
          Fbb[3]= 1.19375;   Fcc[3]= 1.19375;   Fc[3]= 0.98628;   Fll[3]= 0.961865; 
             cafac[3]= 0.890461;  
             winpretag[3]= 109494 ;   wintag[3]= 12398.4 ;   woutpretag[3]= 97500.1;   wouttag[3]= 11804.5;  
           fbb[3]= 0.0585711;   fcc[3]= 0.118364;   fc[3]= 0.152125;   fll[3]= 0.67094;  
           fmcbb[3]= 0.0490648;   fmccc[3]= 0.0991532;   fmcc[3]= 0.154241;   fmcll[3]= 0.697541;  
          Fbb[4]= 1.18279;   Fcc[4]= 1.18279;   Fc[4]= 0.977225;   Fll[4]= 0.953033; 
             cafac[4]= 0.954268;  
             winpretag[4]= 23894.9 ;   wintag[4]= 3544.14 ;   woutpretag[4]= 22802.1;   wouttag[4]= 3645.18;  
           fbb[4]= 0.0832768;   fcc[4]= 0.140613;   fc[4]= 0.140415;   fll[4]= 0.635695;  
           fmcbb[4]= 0.0704071;   fmccc[4]= 0.118883;   fmcc[4]= 0.143688;   fmcll[4]= 0.667023;  
          Fbb[5]= 1.17107;   Fcc[5]= 1.17107;   Fc[5]= 0.967545;   Fll[5]= 0.943594; 
             cafac[5]= 0.908396;  
             winpretag[5]= 6236.64 ;   wintag[5]= 1199.86 ;   woutpretag[5]= 5665.33;   wouttag[5]= 1175.49;  
           fbb[5]= 0.109848;   fcc[5]= 0.163983;   fc[5]= 0.129873;   fll[5]= 0.596297;  
           fmcbb[5]= 0.0938009;   fmccc[5]= 0.140028;   fmcc[5]= 0.134229;   fmcll[5]= 0.631942;  
          Fbb[6]= 1.19083;   Fcc[6]= 1.19083;   Fc[6]= 0.983869;   Fll[6]= 0.959513; 
             cafac[6]= 0.903043;  
             winpretag[6]= 139626 ;   wintag[6]= 17142.4 ;   woutpretag[6]= 126088;   wouttag[6]= 16603.3;  
           fbb[6]= 0.0651569;   fcc[6]= 0.12427;   fc[6]= 0.149097;   fll[6]= 0.661477;  
           fmcbb[6]= 0.0547154;   fmccc[6]= 0.104355;   fmcc[6]= 0.151541;   fmcll[6]= 0.689388;  
          Fbb[7]= 1.18035;   Fcc[7]= 1.18035;   Fc[7]= 0.975206;   Fll[7]= 0.951064; 
             cafac[7]= 0.944092;  
             winpretag[7]= 30131.5 ;   wintag[7]= 4744 ;   woutpretag[7]= 28446.9;   wouttag[7]= 4830.21;  
           fbb[7]= 0.08882;   fcc[7]= 0.145488;   fc[7]= 0.138216;   fll[7]= 0.627476;  
           fmcbb[7]= 0.0752491;   fmccc[7]= 0.123259;   fmcc[7]= 0.14173;   fmcll[7]= 0.659762;  
     }  
     else if ( idsys==56   || sysname=="metpileup_down" ) { 
          Fbb[1]= 1.22052;   Fcc[1]= 1.22052;   Fc[1]= 1.01758;   Fll[1]= 0.984675; 
             cafac[1]= 1.0602;  
             winpretag[1]= 2.17614e+06 ;   wintag[1]= 77215.4 ;   woutpretag[1]= 2.30715e+06;   wouttag[1]= 84984.8;  
           fbb[1]= 0.014468;   fcc[1]= 0.0421379;   fc[1]= 0.135681;   fll[1]= 0.807714;  
           fmcbb[1]= 0.011854;   fmccc[1]= 0.0345247;   fmcc[1]= 0.133337;   fmcll[1]= 0.820285;  
          Fbb[2]= 1.20392;   Fcc[2]= 1.20392;   Fc[2]= 1.00374;   Fll[2]= 0.971287; 
             cafac[2]= 0.979846;  
             winpretag[2]= 494317 ;   wintag[2]= 36710.4 ;   woutpretag[2]= 484354;   wouttag[2]= 38136.3;  
           fbb[2]= 0.035991;   fcc[2]= 0.0866043;   fc[2]= 0.155358;   fll[2]= 0.722047;  
           fmcbb[2]= 0.0298947;   fmccc[2]= 0.0719352;   fmcc[2]= 0.154778;   fmcll[2]= 0.743392;  
          Fbb[3]= 1.19107;   Fcc[3]= 1.19107;   Fc[3]= 0.993026;   Fll[3]= 0.960915; 
             cafac[3]= 0.891705;  
             winpretag[3]= 109325 ;   wintag[3]= 12397.5 ;   woutpretag[3]= 97485.7;   wouttag[3]= 11828.7;  
           fbb[3]= 0.0583561;   fcc[3]= 0.118302;   fc[3]= 0.153042;   fll[3]= 0.6703;  
           fmcbb[3]= 0.0489948;   fmccc[3]= 0.0993246;   fmcc[3]= 0.154117;   fmcll[3]= 0.697564;  
          Fbb[4]= 1.18037;   Fcc[4]= 1.18037;   Fc[4]= 0.984112;   Fll[4]= 0.95229; 
             cafac[4]= 0.954028;  
             winpretag[4]= 23883.2 ;   wintag[4]= 3553.44 ;   woutpretag[4]= 22785.3;   wouttag[4]= 3653.77;  
           fbb[4]= 0.0831036;   fcc[4]= 0.140128;   fc[4]= 0.141481;   fll[4]= 0.635288;  
           fmcbb[4]= 0.0704044;   fmccc[4]= 0.118715;   fmcc[4]= 0.143765;   fmcll[4]= 0.667116;  
          Fbb[5]= 1.16888;   Fcc[5]= 1.16888;   Fc[5]= 0.974531;   Fll[5]= 0.943018; 
             cafac[5]= 0.901301;  
             winpretag[5]= 6229.42 ;   wintag[5]= 1197.78 ;   woutpretag[5]= 5614.58;   wouttag[5]= 1164.34;  
           fbb[5]= 0.109332;   fcc[5]= 0.163761;   fc[5]= 0.13024;   fll[5]= 0.596667;  
           fmcbb[5]= 0.0935356;   fmccc[5]= 0.1401;   fmcc[5]= 0.133644;   fmcll[5]= 0.63272;  
          Fbb[6]= 1.18822;   Fcc[6]= 1.18822;   Fc[6]= 0.990649;   Fll[6]= 0.958615; 
             cafac[6]= 0.903592;  
             winpretag[6]= 139438 ;   wintag[6]= 17148.7 ;   woutpretag[6]= 125995;   wouttag[6]= 16627.9;  
           fbb[6]= 0.0649381;   fcc[6]= 0.12413;   fc[6]= 0.150013;   fll[6]= 0.660919;  
           fmcbb[6]= 0.0546518;   fmccc[6]= 0.104468;   fmcc[6]= 0.151429;   fmcll[6]= 0.689452;  
          Fbb[7]= 1.17798;   Fcc[7]= 1.17798;   Fc[7]= 0.982115;   Fll[7]= 0.950357; 
             cafac[7]= 0.942325;  
             winpretag[7]= 30112.6 ;   wintag[7]= 4751.22 ;   woutpretag[7]= 28375.9;   wouttag[7]= 4828.43;  
           fbb[7]= 0.0885717;   fcc[7]= 0.145055;   fc[7]= 0.139137;   fll[7]= 0.627236;  
           fmcbb[7]= 0.0751895;   fmccc[7]= 0.123139;   fmcc[7]= 0.141671;   fmcll[7]= 0.660001;  
     }  
     else if ( idsys==87   || sysname=="jvf_up" ) { 
          Fbb[1]= 1.28345;   Fcc[1]= 1.28345;   Fc[1]= 0.996901;   Fll[1]= 0.984456; 
             cafac[1]= 1.06375;  
             winpretag[1]= 2.15996e+06 ;   wintag[1]= 76555.9 ;   woutpretag[1]= 2.29766e+06;   wouttag[1]= 84318.9;  
           fbb[1]= 0.0152099;   fcc[1]= 0.0443805;   fc[1]= 0.133133;   fll[1]= 0.807277;  
           fmcbb[1]= 0.0118508;   fmccc[1]= 0.034579;   fmcc[1]= 0.133547;   fmcll[1]= 0.820024;  
          Fbb[2]= 1.26212;   Fcc[2]= 1.26212;   Fc[2]= 0.980331;   Fll[2]= 0.968093; 
             cafac[2]= 0.987008;  
             winpretag[2]= 488512 ;   wintag[2]= 36142.2 ;   woutpretag[2]= 482165;   wouttag[2]= 38129;  
           fbb[2]= 0.0378211;   fcc[2]= 0.0910043;   fc[2]= 0.15183;   fll[2]= 0.719344;  
           fmcbb[2]= 0.0299664;   fmccc[2]= 0.0721043;   fmcc[2]= 0.154876;   fmcll[2]= 0.743053;  
          Fbb[3]= 1.24508;   Fcc[3]= 1.24508;   Fc[3]= 0.967093;   Fll[3]= 0.95502; 
             cafac[3]= 0.902454;  
             winpretag[3]= 107668 ;   wintag[3]= 12155.8 ;   woutpretag[3]= 97165.6;   wouttag[3]= 11892;  
           fbb[3]= 0.0611961;   fcc[3]= 0.123891;   fc[3]= 0.149148;   fll[3]= 0.665765;  
           fmcbb[3]= 0.0491505;   fmccc[3]= 0.0995048;   fmcc[3]= 0.154223;   fmcll[3]= 0.697122;  
          Fbb[4]= 1.23063;   Fcc[4]= 1.23063;   Fc[4]= 0.955871;   Fll[4]= 0.943938; 
             cafac[4]= 0.973448;  
             winpretag[4]= 23416.1 ;   wintag[4]= 3470.34 ;   woutpretag[4]= 22794.4;   wouttag[4]= 3699.33;  
           fbb[4]= 0.0868628;   fcc[4]= 0.146402;   fc[4]= 0.137768;   fll[4]= 0.628967;  
           fmcbb[4]= 0.0705841;   fmccc[4]= 0.118965;   fmcc[4]= 0.144129;   fmcll[4]= 0.666322;  
          Fbb[5]= 1.21512;   Fcc[5]= 1.21512;   Fc[5]= 0.943824;   Fll[5]= 0.932042; 
             cafac[5]= 0.931774;  
             winpretag[5]= 6071.02 ;   wintag[5]= 1163.52 ;   woutpretag[5]= 5656.82;   wouttag[5]= 1189.4;  
           fbb[5]= 0.114268;   fcc[5]= 0.170677;   fc[5]= 0.126287;   fll[5]= 0.588768;  
           fmcbb[5]= 0.0940385;   fmccc[5]= 0.140461;   fmcc[5]= 0.133804;   fmcll[5]= 0.631697;  
          Fbb[6]= 1.24123;   Fcc[6]= 1.24123;   Fc[6]= 0.964108;   Fll[6]= 0.952073; 
             cafac[6]= 0.916831;  
             winpretag[6]= 137155 ;   wintag[6]= 16789.6 ;   woutpretag[6]= 125748;   wouttag[6]= 16753.9;  
           fbb[6]= 0.0680155;   fcc[6]= 0.129883;   fc[6]= 0.146155;   fll[6]= 0.655947;  
           fmcbb[6]= 0.0547967;   fmccc[6]= 0.10464;   fmcc[6]= 0.151596;   fmcll[6]= 0.688968;  
          Fbb[7]= 1.2274;   Fcc[7]= 1.2274;   Fc[7]= 0.953366;   Fll[7]= 0.941464; 
             cafac[7]= 0.964205;  
             winpretag[7]= 29487.1 ;   wintag[7]= 4633.87 ;   woutpretag[7]= 28431.6;   wouttag[7]= 4897.85;  
           fbb[7]= 0.0925622;   fcc[7]= 0.15145;   fc[7]= 0.135381;   fll[7]= 0.620607;  
           fmcbb[7]= 0.0754131;   fmccc[7]= 0.123391;   fmcc[7]= 0.142003;   fmcll[7]= 0.659193;  
     }  
     else if ( idsys==88   || sysname=="jvf_down" ) { 
          Fbb[1]= 1.24477;   Fcc[1]= 1.24477;   Fc[1]= 0.979249;   Fll[1]= 0.989525; 
             cafac[1]= 1.04156;  
             winpretag[1]= 2.19698e+06 ;   wintag[1]= 78079.9 ;   woutpretag[1]= 2.28828e+06;   wouttag[1]= 83019.6;  
           fbb[1]= 0.0147514;   fcc[1]= 0.0430107;   fc[1]= 0.130523;   fll[1]= 0.811715;  
           fmcbb[1]= 0.0118507;   fmccc[1]= 0.0345531;   fmcc[1]= 0.133289;   fmcll[1]= 0.820307;  
          Fbb[2]= 1.22766;   Fcc[2]= 1.22766;   Fc[2]= 0.965794;   Fll[2]= 0.975929; 
             cafac[2]= 0.961839;  
             winpretag[2]= 499830 ;   wintag[2]= 37199.4 ;   woutpretag[2]= 480755;   wouttag[2]= 37666.4;  
           fbb[2]= 0.0366991;   fcc[2]= 0.0883353;   fc[2]= 0.149386;   fll[2]= 0.72558;  
           fmcbb[2]= 0.0298934;   fmccc[2]= 0.0719539;   fmcc[2]= 0.154677;   fmcll[2]= 0.743476;  
          Fbb[3]= 1.21345;   Fcc[3]= 1.21345;   Fc[3]= 0.954615;   Fll[3]= 0.964633; 
             cafac[3]= 0.873383;  
             winpretag[3]= 110801 ;   wintag[3]= 12572 ;   woutpretag[3]= 96771.4;   wouttag[3]= 11723.6;  
           fbb[3]= 0.059519;   fcc[3]= 0.120482;   fc[3]= 0.146954;   fll[3]= 0.673045;  
           fmcbb[3]= 0.0490492;   fmccc[3]= 0.0992882;   fmcc[3]= 0.153941;   fmcll[3]= 0.697721;  
          Fbb[4]= 1.20113;   Fcc[4]= 1.20113;   Fc[4]= 0.944921;   Fll[4]= 0.954837; 
             cafac[4]= 0.933901;  
             winpretag[4]= 24240.6 ;   wintag[4]= 3605.47 ;   woutpretag[4]= 22638.4;   wouttag[4]= 3632.52;  
           fbb[4]= 0.0846621;   fcc[4]= 0.142544;   fc[4]= 0.135918;   fll[4]= 0.636875;  
           fmcbb[4]= 0.0704853;   fmccc[4]= 0.118675;   fmcc[4]= 0.143841;   fmcll[4]= 0.666999;  
          Fbb[5]= 1.18809;   Fcc[5]= 1.18809;   Fc[5]= 0.934662;   Fll[5]= 0.944471; 
             cafac[5]= 0.881885;  
             winpretag[5]= 6351.01 ;   wintag[5]= 1226.12 ;   woutpretag[5]= 5600.86;   wouttag[5]= 1168.17;  
           fbb[5]= 0.11123;   fcc[5]= 0.165962;   fc[5]= 0.124792;   fll[5]= 0.598015;  
           fmcbb[5]= 0.093621;   fmccc[5]= 0.139688;   fmcc[5]= 0.133516;   fmcll[5]= 0.633175;  
          Fbb[6]= 1.21017;   Fcc[6]= 1.21017;   Fc[6]= 0.952027;   Fll[6]= 0.962018; 
             cafac[6]= 0.884852;  
             winpretag[6]= 141392 ;   wintag[6]= 17403.6 ;   woutpretag[6]= 125111;   wouttag[6]= 16506.9;  
           fbb[6]= 0.066228;   fcc[6]= 0.126373;   fc[6]= 0.144034;   fll[6]= 0.663364;  
           fmcbb[6]= 0.0547263;   fmccc[6]= 0.104427;   fmcc[6]= 0.151292;   fmcll[6]= 0.689555;  
          Fbb[7]= 1.1984;   Fcc[7]= 1.1984;   Fc[7]= 0.942773;   Fll[7]= 0.952666; 
             cafac[7]= 0.922331;  
             winpretag[7]= 30591.6 ;   wintag[7]= 4831.59 ;   woutpretag[7]= 28215.6;   wouttag[7]= 4811.71;  
           fbb[7]= 0.0902257;   fcc[7]= 0.147448;   fc[7]= 0.133588;   fll[7]= 0.628738;  
           fmcbb[7]= 0.0752884;   fmccc[7]= 0.123037;   fmcc[7]= 0.141697;   fmcll[7]= 0.659977;  
     }  
     else if ( idsys==89   || sysname=="pdfs_up" ) { 
          Fbb[1]= 1.15489;   Fcc[1]= 1.15489;   Fc[1]= 0.941073;   Fll[1]= 1.00121; 
             cafac[1]= 1.0382;  
             winpretag[1]= 2.2481e+06 ;   wintag[1]= 81568.4 ;   woutpretag[1]= 2.33397e+06;   wouttag[1]= 83652.3;  
           fbb[1]= 0.0135645;   fcc[1]= 0.0394877;   fc[1]= 0.129455;   fll[1]= 0.817493;  
           fmcbb[1]= 0.0117453;   fmccc[1]= 0.0341918;   fmcc[1]= 0.137561;   fmcll[1]= 0.816502;  
          Fbb[2]= 1.14645;   Fcc[2]= 1.14645;   Fc[2]= 0.934194;   Fll[2]= 0.993896; 
             cafac[2]= 0.95132;  
             winpretag[2]= 509470 ;   wintag[2]= 38052.1 ;   woutpretag[2]= 484670;   wouttag[2]= 36820.3;  
           fbb[2]= 0.034122;   fcc[2]= 0.082091;   fc[2]= 0.146455;   fll[2]= 0.737332;  
           fmcbb[2]= 0.0297633;   fmccc[2]= 0.0716048;   fmcc[2]= 0.156772;   fmcll[2]= 0.74186;  
          Fbb[3]= 1.13842;   Fcc[3]= 1.13842;   Fc[3]= 0.927657;   Fll[3]= 0.986941; 
             cafac[3]= 0.857534;  
             winpretag[3]= 113198 ;   wintag[3]= 12933.2 ;   woutpretag[3]= 97070.8;   wouttag[3]= 11437.8;  
           fbb[3]= 0.0557627;   fcc[3]= 0.112289;   fc[3]= 0.145569;   fll[3]= 0.686379;  
           fmcbb[3]= 0.0489823;   fmccc[3]= 0.0986355;   fmcc[3]= 0.156921;   fmcll[3]= 0.695461;  
          Fbb[4]= 1.13061;   Fcc[4]= 1.13061;   Fc[4]= 0.921293;   Fll[4]= 0.98017; 
             cafac[4]= 0.914157;  
             winpretag[4]= 24894.4 ;   wintag[4]= 3751.03 ;   woutpretag[4]= 22757.4;   wouttag[4]= 3570.08;  
           fbb[4]= 0.0797031;   fcc[4]= 0.134769;   fc[4]= 0.136269;   fll[4]= 0.649258;  
           fmcbb[4]= 0.0704954;   fmccc[4]= 0.1192;   fmcc[4]= 0.147911;   fmcll[4]= 0.662394;  
          Fbb[5]= 1.12308;   Fcc[5]= 1.12308;   Fc[5]= 0.915158;   Fll[5]= 0.973643; 
             cafac[5]= 0.846138;  
             winpretag[5]= 6595.3 ;   wintag[5]= 1348.17 ;   woutpretag[5]= 5580.54;   wouttag[5]= 1188.56;  
           fbb[5]= 0.1052;   fcc[5]= 0.156511;   fc[5]= 0.132489;   fll[5]= 0.6058;  
           fmcbb[5]= 0.0936707;   fmccc[5]= 0.139358;   fmcc[5]= 0.144772;   fmcll[5]= 0.622199;  
          Fbb[6]= 1.13637;   Fcc[6]= 1.13637;   Fc[6]= 0.92598;   Fll[6]= 0.985157; 
             cafac[6]= 0.86725;  
             winpretag[6]= 144687 ;   wintag[6]= 18032.4 ;   woutpretag[6]= 125480;   wouttag[6]= 16183.2;  
           fbb[6]= 0.0621829;   fcc[6]= 0.118216;   fc[6]= 0.143357;   fll[6]= 0.676243;  
           fmcbb[6]= 0.0547209;   fmccc[6]= 0.10403;   fmcc[6]= 0.154817;   fmcll[6]= 0.686432;  
          Fbb[7]= 1.12903;   Fcc[7]= 1.12903;   Fc[7]= 0.920001;   Fll[7]= 0.978796; 
             cafac[7]= 0.898876;  
             winpretag[7]= 31489.7 ;   wintag[7]= 5099.2 ;   woutpretag[7]= 28305.4;   wouttag[7]= 4774.8;  
           fbb[7]= 0.0850716;   fcc[7]= 0.139347;   fc[7]= 0.135473;   fll[7]= 0.640108;  
           fmcbb[7]= 0.0753493;   fmccc[7]= 0.123422;   fmcc[7]= 0.147253;   fmcll[7]= 0.653975;  
     }  
     else if ( idsys==90   || sysname=="pdfs_down" ) { 
          Fbb[1]= 1.35314;   Fcc[1]= 1.35314;   Fc[1]= 1.10491;   Fll[1]= 0.963532; 
             cafac[1]= 1.07768;  
             winpretag[1]= 2.11101e+06 ;   wintag[1]= 73163.8 ;   woutpretag[1]= 2.27499e+06;   wouttag[1]= 86870.1;  
           fbb[1]= 0.0161568;   fcc[1]= 0.0472365;   fc[1]= 0.142356;   fll[1]= 0.79425;  
           fmcbb[1]= 0.0119402;   fmccc[1]= 0.0349088;   fmcc[1]= 0.12884;   fmcll[1]= 0.824311;  
          Fbb[2]= 1.32012;   Fcc[2]= 1.32012;   Fc[2]= 1.07795;   Fll[2]= 0.940018; 
             cafac[2]= 1.01018;  
             winpretag[2]= 479909 ;   wintag[2]= 35442.6 ;   woutpretag[2]= 484795;   wouttag[2]= 40165.9;  
           fbb[2]= 0.039689;   fcc[2]= 0.0955508;   fc[2]= 0.164457;   fll[2]= 0.700303;  
           fmcbb[2]= 0.0300648;   fmccc[2]= 0.0723807;   fmcc[2]= 0.152565;   fmcll[2]= 0.74499;  
          Fbb[3]= 1.29731;   Fcc[3]= 1.29731;   Fc[3]= 1.05933;   Fll[3]= 0.923781; 
             cafac[3]= 0.929734;  
             winpretag[3]= 105636 ;   wintag[3]= 11858.2 ;   woutpretag[3]= 98213.2;   wouttag[3]= 12457.1;  
           fbb[3]= 0.0637661;   fcc[3]= 0.129868;   fc[3]= 0.159948;   fll[3]= 0.646417;  
           fmcbb[3]= 0.0491524;   fmccc[3]= 0.100105;   fmcc[3]= 0.150991;   fmcll[3]= 0.699751;  
          Fbb[4]= 1.28035;   Fcc[4]= 1.28035;   Fc[4]= 1.04548;   Fll[4]= 0.911703; 
             cafac[4]= 1.00024;  
             winpretag[4]= 22873.4 ;   wintag[4]= 3346.43 ;   woutpretag[4]= 22878.9;   wouttag[4]= 3806.88;  
           fbb[4]= 0.0902717;   fcc[4]= 0.151502;   fc[4]= 0.146016;   fll[4]= 0.61221;  
           fmcbb[4]= 0.0705054;   fmccc[4]= 0.118328;   fmcc[4]= 0.139665;   fmcll[4]= 0.671501;  
          Fbb[5]= 1.26219;   Fcc[5]= 1.26219;   Fc[5]= 1.03065;   Fll[5]= 0.898771; 
             cafac[5]= 0.982336;  
             winpretag[5]= 5869.52 ;   wintag[5]= 1054.41 ;   woutpretag[5]= 5765.84;   wouttag[5]= 1178.36;  
           fbb[5]= 0.118492;   fcc[5]= 0.177607;   fc[5]= 0.12484;   fll[5]= 0.579061;  
           fmcbb[5]= 0.093878;   fmccc[5]= 0.140714;   fmcc[5]= 0.121128;   fmcll[5]= 0.64428;  
          Fbb[6]= 1.29283;   Fcc[6]= 1.29283;   Fc[6]= 1.05566;   Fll[6]= 0.920586; 
             cafac[6]= 0.945265;  
             winpretag[6]= 134379 ;   wintag[6]= 16259.1 ;   woutpretag[6]= 127024;   wouttag[6]= 17415.5;  
           fbb[6]= 0.0707701;   fcc[6]= 0.135722;   fc[6]= 0.155983;   fll[6]= 0.637524;  
           fmcbb[6]= 0.0547406;   fmccc[6]= 0.104981;   fmcc[6]= 0.147758;   fmcll[6]= 0.69252;  
          Fbb[7]= 1.2766;   Fcc[7]= 1.2766;   Fc[7]= 1.04241;   Fll[7]= 0.909032; 
             cafac[7]= 0.996264;  
             winpretag[7]= 28742.9 ;   wintag[7]= 4400.84 ;   woutpretag[7]= 28635.5;   wouttag[7]= 4989.33;  
           fbb[7]= 0.0961003;   fcc[7]= 0.156894;   fc[7]= 0.141642;   fll[7]= 0.605363;  
           fmcbb[7]= 0.0752783;   fmccc[7]= 0.1229;   fmcc[7]= 0.135879;   fmcll[7]= 0.665943;  
     }  
     else if ( idsys==91   || sysname=="pdfa_up" ) { 
          Fbb[1]= 1.03974;   Fcc[1]= 1.03974;   Fc[1]= 1.23865;   Fll[1]= 0.958952; 
             cafac[1]= 0.997386;  
             winpretag[1]= 2.22128e+06 ;   wintag[1]= 78958.7 ;   woutpretag[1]= 2.21547e+06;   wouttag[1]= 89181.2;  
           fbb[1]= 0.0123354;   fcc[1]= 0.0359443;   fc[1]= 0.165171;   fll[1]= 0.786549;  
           fmcbb[1]= 0.011864;   fmccc[1]= 0.0345706;   fmcc[1]= 0.133348;   fmcll[1]= 0.820218;  
          Fbb[2]= 1.02891;   Fcc[2]= 1.02891;   Fc[2]= 1.22575;   Fll[2]= 0.948966; 
             cafac[2]= 0.935019;  
             winpretag[2]= 504376 ;   wintag[2]= 37540.1 ;   woutpretag[2]= 471601;   wouttag[2]= 38344.1;  
           fbb[2]= 0.0308329;   fcc[2]= 0.0741325;   fc[2]= 0.189892;   fll[2]= 0.705143;  
           fmcbb[2]= 0.0299666;   fmccc[2]= 0.0720497;   fmcc[2]= 0.154919;   fmcll[2]= 0.743065;  
          Fbb[3]= 1.02529;   Fcc[3]= 1.02529;   Fc[3]= 1.22144;   Fll[3]= 0.94563; 
             cafac[3]= 0.843804;  
             winpretag[3]= 111905 ;   wintag[3]= 12689.8 ;   woutpretag[3]= 94426;   wouttag[3]= 11467.1;  
           fbb[3]= 0.0504725;   fcc[3]= 0.102065;   fc[3]= 0.188297;   fll[3]= 0.659165;  
           fmcbb[3]= 0.0492276;   fmccc[3]= 0.0995473;   fmcc[3]= 0.15416;   fmcll[3]= 0.697065;  
          Fbb[4]= 1.02453;   Fcc[4]= 1.02453;   Fc[4]= 1.22054;   Fll[4]= 0.944931; 
             cafac[4]= 0.872373;  
             winpretag[4]= 24637.9 ;   wintag[4]= 3655.98 ;   woutpretag[4]= 21493.4;   wouttag[4]= 3375.71;  
           fbb[4]= 0.0721278;   fcc[4]= 0.122733;   fc[4]= 0.176827;   fll[4]= 0.628313;  
           fmcbb[4]= 0.0704006;   fmccc[4]= 0.119794;   fmcc[4]= 0.144876;   fmcll[4]= 0.664929;  
          Fbb[5]= 1.02329;   Fcc[5]= 1.02329;   Fc[5]= 1.21905;   Fll[5]= 0.94378; 
             cafac[5]= 0.805948;  
             winpretag[5]= 6499.03 ;   wintag[5]= 1265.95 ;   woutpretag[5]= 5237.88;   wouttag[5]= 1072.43;  
           fbb[5]= 0.0956924;   fcc[5]= 0.143383;   fc[5]= 0.166711;   fll[5]= 0.594214;  
           fmcbb[5]= 0.0935148;   fmccc[5]= 0.14012;   fmcc[5]= 0.136754;   fmcll[5]= 0.62961;  
          Fbb[6]= 1.02507;   Fcc[6]= 1.02507;   Fc[6]= 1.22118;   Fll[6]= 0.945425; 
             cafac[6]= 0.84712;  
             winpretag[6]= 143042 ;   wintag[6]= 17611.7 ;   woutpretag[6]= 121174;   wouttag[6]= 15918.5;  
           fbb[6]= 0.0562626;   fcc[6]= 0.107507;   fc[6]= 0.185338;   fll[6]= 0.650892;  
           fmcbb[6]= 0.0548866;   fmccc[6]= 0.104878;   fmcc[6]= 0.15177;   fmcll[6]= 0.688465;  
          Fbb[7]= 1.02427;   Fcc[7]= 1.02427;   Fc[7]= 1.22023;   Fll[7]= 0.944691; 
             cafac[7]= 0.857454;  
             winpretag[7]= 31136.9 ;   wintag[7]= 4921.94 ;   woutpretag[7]= 26698.5;   wouttag[7]= 4459.21;  
           fbb[7]= 0.0770511;   fcc[7]= 0.127047;   fc[7]= 0.174713;   fll[7]= 0.621189;  
           fmcbb[7]= 0.0752251;   fmccc[7]= 0.124036;   fmcc[7]= 0.143181;   fmcll[7]= 0.657557;  
     }  
     else if ( idsys==92   || sysname=="pdfa_down" ) { 
          Fbb[1]= 1.44874;   Fcc[1]= 1.44874;   Fc[1]= 0.804668;   Fll[1]= 1.00641; 
             cafac[1]= 1.12704;  
             winpretag[1]= 2.13783e+06 ;   wintag[1]= 75773.5 ;   woutpretag[1]= 2.40941e+06;   wouttag[1]= 81497.8;  
           fbb[1]= 0.0171162;   fcc[1]= 0.0499907;   fc[1]= 0.107284;   fll[1]= 0.825609;  
           fmcbb[1]= 0.0118145;   fmccc[1]= 0.0345063;   fmcc[1]= 0.133327;   fmcll[1]= 0.820352;  
          Fbb[2]= 1.42;   Fcc[2]= 1.42;   Fc[2]= 0.788702;   Fll[2]= 0.986439; 
             cafac[2]= 1.03083;  
             winpretag[2]= 485003 ;   wintag[2]= 35954.7 ;   woutpretag[2]= 499954;   wouttag[2]= 38642.1;  
           fbb[2]= 0.0423872;   fcc[2]= 0.102112;   fc[2]= 0.121883;   fll[2]= 0.733618;  
           fmcbb[2]= 0.0298502;   fmccc[2]= 0.0719098;   fmcc[2]= 0.154536;   fmcll[2]= 0.743704;  
          Fbb[3]= 1.39191;   Fcc[3]= 1.39191;   Fc[3]= 0.773104;   Fll[3]= 0.966931; 
             cafac[3]= 0.949885;  
             winpretag[3]= 106928 ;   wintag[3]= 12101.6 ;   woutpretag[3]= 101570;   wouttag[3]= 12470.5;  
           fbb[3]= 0.0680559;   fcc[3]= 0.137985;   fc[3]= 0.119021;   fll[3]= 0.674938;  
           fmcbb[3]= 0.0488937;   fmccc[3]= 0.0991335;   fmcc[3]= 0.153952;   fmcll[3]= 0.698021;  
          Fbb[4]= 1.36563;   Fcc[4]= 1.36563;   Fc[4]= 0.758507;   Fll[4]= 0.948673; 
             cafac[4]= 1.05987;  
             winpretag[4]= 23129.9 ;   wintag[4]= 3441.48 ;   woutpretag[4]= 24514.8;   wouttag[4]= 4060.7;  
           fbb[4]= 0.0964222;   fcc[4]= 0.160742;   fc[4]= 0.108458;   fll[4]= 0.634377;  
           fmcbb[4]= 0.0706063;   fmccc[4]= 0.117705;   fmcc[4]= 0.142989;   fmcll[4]= 0.668699;  
          Fbb[5]= 1.33698;   Fcc[5]= 1.33698;   Fc[5]= 0.742589;   Fll[5]= 0.928765; 
             cafac[5]= 1.0417;  
             winpretag[5]= 5965.79 ;   wintag[5]= 1136.62 ;   woutpretag[5]= 6214.54;   wouttag[5]= 1330.37;  
           fbb[5]= 0.125735;   fcc[5]= 0.186992;   fc[5]= 0.0967174;   fll[5]= 0.590556;  
           fmcbb[5]= 0.0940445;   fmccc[5]= 0.139862;   fmcc[5]= 0.130243;   fmcll[5]= 0.63585;  
          Fbb[6]= 1.38489;   Fcc[6]= 1.38489;   Fc[6]= 0.769201;   Fll[6]= 0.962049; 
             cafac[6]= 0.974083;  
             winpretag[6]= 136024 ;   wintag[6]= 16679.7 ;   woutpretag[6]= 132499;   wouttag[6]= 17796.8;  
           fbb[6]= 0.0755678;   fcc[6]= 0.144136;   fc[6]= 0.116186;   fll[6]= 0.66411;  
           fmcbb[6]= 0.054566;   fmccc[6]= 0.104078;   fmcc[6]= 0.151048;   fmcll[6]= 0.690308;  
          Fbb[7]= 1.35966;   Fcc[7]= 1.35966;   Fc[7]= 0.755188;   Fll[7]= 0.944522; 
             cafac[7]= 1.0558;  
             winpretag[7]= 29095.7 ;   wintag[7]= 4578.1 ;   woutpretag[7]= 30719.3;   wouttag[7]= 5398.65;  
           fbb[7]= 0.102535;   fcc[7]= 0.166216;   fc[7]= 0.10601;   fll[7]= 0.62524;  
           fmcbb[7]= 0.075412;   fmccc[7]= 0.122248;   fmcc[7]= 0.140376;   fmcll[7]= 0.661964;  
     }  
     else if ( idsys==53   || sysname=="met_up" ) { 
          Fbb[1]= 1.23769;   Fcc[1]= 1.23769;   Fc[1]= 1.00356;   Fll[1]= 0.985989; 
             cafac[1]= 1.05193;  
             winpretag[1]= 2.18462e+06 ;   wintag[1]= 77640.6 ;   woutpretag[1]= 2.29807e+06;   wouttag[1]= 84352.9;  
           fbb[1]= 0.0146411;   fcc[1]= 0.0427313;   fc[1]= 0.133886;   fll[1]= 0.808741;  
           fmcbb[1]= 0.0118294;   fmccc[1]= 0.0345252;   fmcc[1]= 0.133411;   fmcll[1]= 0.820234;  
          Fbb[2]= 1.22016;   Fcc[2]= 1.22016;   Fc[2]= 0.989352;   Fll[2]= 0.972031; 
             cafac[2]= 0.975838;  
             winpretag[2]= 495288 ;   wintag[2]= 36815 ;   woutpretag[2]= 483321;   wouttag[2]= 38063.6;  
           fbb[2]= 0.036531;   fcc[2]= 0.0878182;   fc[2]= 0.15315;   fll[2]= 0.722501;  
           fmcbb[2]= 0.0299394;   fmccc[2]= 0.0719724;   fmcc[2]= 0.154798;   fmcll[2]= 0.74329;  
          Fbb[3]= 1.2063;   Fcc[3]= 1.2063;   Fc[3]= 0.978108;   Fll[3]= 0.960984; 
             cafac[3]= 0.889955;  
             winpretag[3]= 109536 ;   wintag[3]= 12410.2 ;   woutpretag[3]= 97482.5;   wouttag[3]= 11835.9;  
           fbb[3]= 0.0592008;   fcc[3]= 0.119662;   fc[3]= 0.150915;   fll[3]= 0.670222;  
           fmcbb[3]= 0.0490764;   fmccc[3]= 0.0991981;   fmcc[3]= 0.154293;   fmcll[3]= 0.697433;  
          Fbb[4]= 1.19454;   Fcc[4]= 1.19454;   Fc[4]= 0.968574;   Fll[4]= 0.951617; 
             cafac[4]= 0.951588;  
             winpretag[4]= 23901.6 ;   wintag[4]= 3545.37 ;   woutpretag[4]= 22744.5;   wouttag[4]= 3647.62;  
           fbb[4]= 0.0840058;   fcc[4]= 0.141927;   fc[4]= 0.139221;   fll[4]= 0.634846;  
           fmcbb[4]= 0.0703248;   fmccc[4]= 0.118813;   fmcc[4]= 0.143738;   fmcll[4]= 0.667124;  
          Fbb[5]= 1.18205;   Fcc[5]= 1.18205;   Fc[5]= 0.958445;   Fll[5]= 0.941665; 
             cafac[5]= 0.909697;  
             winpretag[5]= 6233.78 ;   wintag[5]= 1197.99 ;   woutpretag[5]= 5670.85;   wouttag[5]= 1179.5;  
           fbb[5]= 0.110712;   fcc[5]= 0.165103;   fc[5]= 0.128236;   fll[5]= 0.595949;  
           fmcbb[5]= 0.0936612;   fmccc[5]= 0.139675;   fmcc[5]= 0.133796;   fmcll[5]= 0.632867;  
          Fbb[6]= 1.20317;   Fcc[6]= 1.20317;   Fc[6]= 0.975571;   Fll[6]= 0.958492; 
             cafac[6]= 0.902223;  
             winpretag[6]= 139672 ;   wintag[6]= 17153.6 ;   woutpretag[6]= 126015;   wouttag[6]= 16642.1;  
           fbb[6]= 0.0658164;   fcc[6]= 0.125564;   fc[6]= 0.147869;   fll[6]= 0.66075;  
           fmcbb[6]= 0.0547025;   fmccc[6]= 0.104361;   fmcc[6]= 0.151572;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.19193;   Fcc[7]= 1.19193;   Fc[7]= 0.966461;   Fll[7]= 0.949541; 
             cafac[7]= 0.942317;  
             winpretag[7]= 30135.4 ;   wintag[7]= 4743.36 ;   woutpretag[7]= 28397.1;   wouttag[7]= 4836.21;  
           fbb[7]= 0.0895764;   fcc[7]= 0.146761;   fc[7]= 0.13693;   fll[7]= 0.626732;  
           fmcbb[7]= 0.0751522;   fmccc[7]= 0.123129;   fmcc[7]= 0.141682;   fmcll[7]= 0.660037;  
     }  
     else if ( idsys==54   || sysname=="met_down" ) { 
          Fbb[1]= 1.21085;   Fcc[1]= 1.21085;   Fc[1]= 1.01979;   Fll[1]= 0.984858; 
             cafac[1]= 1.0622;  
             winpretag[1]= 2.1738e+06 ;   wintag[1]= 77083.5 ;   woutpretag[1]= 2.30901e+06;   wouttag[1]= 84994.4;  
           fbb[1]= 0.0143567;   fcc[1]= 0.04182;   fc[1]= 0.135946;   fll[1]= 0.807877;  
           fmcbb[1]= 0.0118568;   fmccc[1]= 0.0345379;   fmcc[1]= 0.133307;   fmcll[1]= 0.820298;  
          Fbb[2]= 1.19497;   Fcc[2]= 1.19497;   Fc[2]= 1.00642;   Fll[2]= 0.971943; 
             cafac[2]= 0.980914;  
             winpretag[2]= 494091 ;   wintag[2]= 36732.8 ;   woutpretag[2]= 484661;   wouttag[2]= 38143.1;  
           fbb[2]= 0.0357386;   fcc[2]= 0.0860052;   fc[2]= 0.15573;   fll[2]= 0.722526;  
           fmcbb[2]= 0.0299076;   fmccc[2]= 0.0719728;   fmcc[2]= 0.154737;   fmcll[2]= 0.743383;  
          Fbb[3]= 1.18276;   Fcc[3]= 1.18276;   Fc[3]= 0.996142;   Fll[3]= 0.962015; 
             cafac[3]= 0.893286;  
             winpretag[3]= 109326 ;   wintag[3]= 12399.6 ;   woutpretag[3]= 97659.2;   wouttag[3]= 11825.9;  
           fbb[3]= 0.0580088;   fcc[3]= 0.117309;   fc[3]= 0.153662;   fll[3]= 0.67102;  
           fmcbb[3]= 0.0490452;   fmccc[3]= 0.0991823;   fmcc[3]= 0.154257;   fmcll[3]= 0.697515;  
          Fbb[4]= 1.17264;   Fcc[4]= 1.17264;   Fc[4]= 0.987618;   Fll[4]= 0.953784; 
             cafac[4]= 0.958977;  
             winpretag[4]= 23882.3 ;   wintag[4]= 3551.61 ;   woutpretag[4]= 22902.6;   wouttag[4]= 3661.26;  
           fbb[4]= 0.0823377;   fcc[4]= 0.139222;   fc[4]= 0.142017;   fll[4]= 0.636424;  
           fmcbb[4]= 0.0702156;   fmccc[4]= 0.118725;   fmcc[4]= 0.143797;   fmcll[4]= 0.667262;  
          Fbb[5]= 1.16156;   Fcc[5]= 1.16156;   Fc[5]= 0.97829;   Fll[5]= 0.944775; 
             cafac[5]= 0.905853;  
             winpretag[5]= 6219.71 ;   wintag[5]= 1196 ;   woutpretag[5]= 5634.14;   wouttag[5]= 1165.33;  
           fbb[5]= 0.109061;   fcc[5]= 0.162804;   fc[5]= 0.130919;   fll[5]= 0.597217;  
           fmcbb[5]= 0.093891;   fmccc[5]= 0.140159;   fmcc[5]= 0.133824;   fmcll[5]= 0.632126;  
          Fbb[6]= 1.18006;   Fcc[6]= 1.18006;   Fc[6]= 0.993863;   Fll[6]= 0.959815; 
             cafac[6]= 0.905949;  
             winpretag[6]= 139428 ;   wintag[6]= 17147.2 ;   woutpretag[6]= 126314;   wouttag[6]= 16630.8;  
           fbb[6]= 0.064516;   fcc[6]= 0.123148;   fc[6]= 0.150624;   fll[6]= 0.661712;  
           fmcbb[6]= 0.054672;   fmccc[6]= 0.104358;   fmcc[6]= 0.151554;   fmcll[6]= 0.689416;  
          Fbb[7]= 1.17033;   Fcc[7]= 1.17033;   Fc[7]= 0.985676;   Fll[7]= 0.951908; 
             cafac[7]= 0.947206;  
             winpretag[7]= 30102 ;   wintag[7]= 4747.61 ;   woutpretag[7]= 28512.8;   wouttag[7]= 4836.94;  
           fbb[7]= 0.0879009;   fcc[7]= 0.144131;   fc[7]= 0.139706;   fll[7]= 0.628261;  
           fmcbb[7]= 0.0751074;   fmccc[7]= 0.123154;   fmcc[7]= 0.141737;   fmcll[7]= 0.660002;  
     }  
     else if ( idsys==98   || sysname=="ifsr_nom" ) { 
          Fbb[1]= 1.22154;   Fcc[1]= 1.22154;   Fc[1]= 1.01403;   Fll[1]= 0.985194; 
             cafac[1]= 1.05667;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.30306e+06;   wouttag[1]= 84720.3;  
           fbb[1]= 0.0144626;   fcc[1]= 0.0421907;   fc[1]= 0.135208;   fll[1]= 0.808139;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.20498;   Fcc[2]= 1.20498;   Fc[2]= 1.00029;   Fll[2]= 0.971845; 
             cafac[2]= 0.978765;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 484185;   wouttag[2]= 38096.5;  
           fbb[2]= 0.0360405;   fcc[2]= 0.0867361;   fc[2]= 0.154775;   fll[2]= 0.722448;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.19208;   Fcc[3]= 1.19208;   Fc[3]= 0.989573;   Fll[3]= 0.961436; 
             cafac[3]= 0.892509;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97655.4;   wouttag[3]= 11833.1;  
           fbb[3]= 0.0584886;   fcc[3]= 0.118427;   fc[3]= 0.152452;   fll[3]= 0.670632;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.18128;   Fcc[4]= 1.18128;   Fc[4]= 0.980607;   Fll[4]= 0.952725; 
             cafac[4]= 0.95839;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22890.1;   wouttag[4]= 3665.23;  
           fbb[4]= 0.0832803;   fcc[4]= 0.140315;   fc[4]= 0.14117;   fll[4]= 0.635234;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.16972;   Fcc[5]= 1.16972;   Fc[5]= 0.971015;   Fll[5]= 0.943405; 
             cafac[5]= 0.91837;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5723.66;   wouttag[5]= 1189.53;  
           fbb[5]= 0.109683;   fcc[5]= 0.163757;   fc[5]= 0.129765;   fll[5]= 0.596795;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.1892;   Fcc[6]= 1.1892;   Fc[6]= 0.987185;   Fll[6]= 0.959116; 
             cafac[6]= 0.90589;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 126402;   wouttag[6]= 16662.1;  
           fbb[6]= 0.0650854;   fcc[6]= 0.124257;   fc[6]= 0.149478;   fll[6]= 0.66118;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.17887;   Fcc[7]= 1.17887;   Fc[7]= 0.978606;   Fll[7]= 0.950781; 
             cafac[7]= 0.949507;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28595.6;   wouttag[7]= 4863.33;  
           fbb[7]= 0.0887869;   fcc[7]= 0.145204;   fc[7]= 0.138792;   fll[7]= 0.627217;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==47   || sysname=="ifsr_up" ) { 
          Fbb[1]= 1.22935;   Fcc[1]= 1.22935;   Fc[1]= 0.980206;   Fll[1]= 0.99025; 
             cafac[1]= 1.05025;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.28907e+06;   wouttag[1]= 82824.6;  
           fbb[1]= 0.0145552;   fcc[1]= 0.0424607;   fc[1]= 0.130698;   fll[1]= 0.812286;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.21351;   Fcc[2]= 1.21351;   Fc[2]= 0.967571;   Fll[2]= 0.977486; 
             cafac[2]= 0.971682;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 480681;   wouttag[2]= 37474.3;  
           fbb[2]= 0.0362955;   fcc[2]= 0.0873497;   fc[2]= 0.149713;   fll[2]= 0.726641;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.20032;   Fcc[3]= 1.20032;   Fc[3]= 0.957056;   Fll[3]= 0.966863; 
             cafac[3]= 0.886085;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 96952.5;   wouttag[3]= 11683.2;  
           fbb[3]= 0.0588931;   fcc[3]= 0.119246;   fc[3]= 0.147443;   fll[3]= 0.674418;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.18886;   Fcc[4]= 1.18886;   Fc[4]= 0.947917;   Fll[4]= 0.957631; 
             cafac[4]= 0.952351;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22745.9;   wouttag[4]= 3629.97;  
           fbb[4]= 0.0838148;   fcc[4]= 0.141216;   fc[4]= 0.136464;   fll[4]= 0.638505;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.17664;   Fcc[5]= 1.17664;   Fc[5]= 0.938174;   Fll[5]= 0.947787; 
             cafac[5]= 0.913049;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5690.5;   wouttag[5]= 1180.06;  
           fbb[5]= 0.110332;   fcc[5]= 0.164725;   fc[5]= 0.125376;   fll[5]= 0.599567;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.19727;   Fcc[6]= 1.19727;   Fc[6]= 0.954623;   Fll[6]= 0.964405; 
             cafac[6]= 0.899582;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 125521;   wouttag[6]= 16467.1;  
           fbb[6]= 0.065527;   fcc[6]= 0.1251;   fc[6]= 0.144547;   fll[6]= 0.664826;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.18631;   Fcc[7]= 1.18631;   Fc[7]= 0.945884;   Fll[7]= 0.955577; 
             cafac[7]= 0.943621;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28418.4;   wouttag[7]= 4818.58;  
           fbb[7]= 0.0893474;   fcc[7]= 0.146121;   fc[7]= 0.134151;   fll[7]= 0.630381;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==98   || sysname=="ifsr_nom" ) { 
          Fbb[1]= 1.27635;   Fcc[1]= 1.27635;   Fc[1]= 1.03292;   Fll[1]= 0.979024; 
             cafac[1]= 1.06027;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.31091e+06;   wouttag[1]= 86420.5;  
           fbb[1]= 0.0151116;   fcc[1]= 0.0440838;   fc[1]= 0.137728;   fll[1]= 0.803077;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.2542;   Fcc[2]= 1.2542;   Fc[2]= 1.015;   Fll[2]= 0.962036; 
             cafac[2]= 0.982546;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 486055;   wouttag[2]= 38953.4;  
           fbb[2]= 0.0375125;   fcc[2]= 0.0902787;   fc[2]= 0.157052;   fll[2]= 0.715157;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.23743;   Fcc[3]= 1.23743;   Fc[3]= 1.00143;   Fll[3]= 0.94917; 
             cafac[3]= 0.893796;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97796.2;   wouttag[3]= 12076.1;  
           fbb[3]= 0.0607136;   fcc[3]= 0.122932;   fc[3]= 0.154278;   fll[3]= 0.662076;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.22365;   Fcc[4]= 1.22365;   Fc[4]= 0.990281;   Fll[4]= 0.938606; 
             cafac[4]= 0.955606;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22823.6;   wouttag[4]= 3726.51;  
           fbb[4]= 0.0862679;   fcc[4]= 0.145349;   fc[4]= 0.142563;   fll[4]= 0.62582;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.20897;   Fcc[5]= 1.20897;   Fc[5]= 0.978397;   Fll[5]= 0.927343; 
             cafac[5]= 0.899053;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5603.26;   wouttag[5]= 1186.71;  
           fbb[5]= 0.113363;   fcc[5]= 0.169252;   fc[5]= 0.130751;   fll[5]= 0.586634;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.23375;   Fcc[6]= 1.23375;   Fc[6]= 0.998452;   Fll[6]= 0.946352; 
             cafac[6]= 0.90528;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 126316;   wouttag[6]= 16973.7;  
           fbb[6]= 0.0675237;   fcc[6]= 0.128912;   fc[6]= 0.151184;   fll[6]= 0.65238;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.22059;   Fcc[7]= 1.22059;   Fc[7]= 0.987798;   Fll[7]= 0.936253; 
             cafac[7]= 0.943067;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28401.7;   wouttag[7]= 4925.15;  
           fbb[7]= 0.091929;   fcc[7]= 0.150343;   fc[7]= 0.140095;   fll[7]= 0.617633;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==48   || sysname=="ifsr_down" ) { 
          Fbb[1]= 1.28453;   Fcc[1]= 1.28453;   Fc[1]= 0.999227;   Fll[1]= 0.984038; 
             cafac[1]= 1.05383;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29688e+06;   wouttag[1]= 84519.1;  
           fbb[1]= 0.0152084;   fcc[1]= 0.0443664;   fc[1]= 0.133235;   fll[1]= 0.807191;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.26305;   Fcc[2]= 1.26305;   Fc[2]= 0.982519;   Fll[2]= 0.967584; 
             cafac[2]= 0.975462;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 482551;   wouttag[2]= 38331.2;  
           fbb[2]= 0.0377772;   fcc[2]= 0.0909157;   fc[2]= 0.152026;   fll[2]= 0.719281;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.24594;   Fcc[3]= 1.24594;   Fc[3]= 0.969207;   Fll[3]= 0.954475; 
             cafac[3]= 0.887405;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97096.9;   wouttag[3]= 11926.9;  
           fbb[3]= 0.0611312;   fcc[3]= 0.123778;   fc[3]= 0.149314;   fll[3]= 0.665777;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.23145;   Fcc[4]= 1.23145;   Fc[4]= 0.957938;   Fll[4]= 0.943376; 
             cafac[4]= 0.949645;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22681.2;   wouttag[4]= 3691.71;  
           fbb[4]= 0.0868175;   fcc[4]= 0.146275;   fc[4]= 0.137907;   fll[4]= 0.629001;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.21605;   Fcc[5]= 1.21605;   Fc[5]= 0.945958;   Fll[5]= 0.931579; 
             cafac[5]= 0.893935;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5571.37;   wouttag[5]= 1177.59;  
           fbb[5]= 0.114027;   fcc[5]= 0.170243;   fc[5]= 0.126416;   fll[5]= 0.589314;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.24207;   Fcc[6]= 1.24207;   Fc[6]= 0.966201;   Fll[6]= 0.951514; 
             cafac[6]= 0.899023;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 125443;   wouttag[6]= 16780.3;  
           fbb[6]= 0.0679791;   fcc[6]= 0.129782;   fc[6]= 0.1463;   fll[6]= 0.655939;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.22823;   Fcc[7]= 1.22823;   Fc[7]= 0.955434;   Fll[7]= 0.940911; 
             cafac[7]= 0.937286;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28227.6;   wouttag[7]= 4881.21;  
           fbb[7]= 0.0925047;   fcc[7]= 0.151285;   fc[7]= 0.135505;   fll[7]= 0.620706;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==72   || sysname=="iqopt3_one" ) { 
          Fbb[1]= 1.25654;   Fcc[1]= 1.25654;   Fc[1]= 0.989957;   Fll[1]= 0.98716; 
             cafac[1]= 1.05747;  
             winpretag[1]= 2.16895e+06 ;   wintag[1]= 76874.1 ;   woutpretag[1]= 2.29361e+06;   wouttag[1]= 83573.3;  
           fbb[1]= 0.0148454;   fcc[1]= 0.0433079;   fc[1]= 0.131985;   fll[1]= 0.809861;  
           fmcbb[1]= 0.0118145;   fmccc[1]= 0.034466;   fmcc[1]= 0.133324;   fmcll[1]= 0.820395;  
          Fbb[2]= 1.23788;   Fcc[2]= 1.23788;   Fc[2]= 0.975256;   Fll[2]= 0.972501; 
             cafac[2]= 0.971877;  
             winpretag[2]= 495238 ;   wintag[2]= 36811.3 ;   woutpretag[2]= 481310;   wouttag[2]= 37902.8;  
           fbb[2]= 0.0370791;   fcc[2]= 0.0892045;   fc[2]= 0.150808;   fll[2]= 0.722908;  
           fmcbb[2]= 0.0299537;   fmccc[2]= 0.0720623;   fmcc[2]= 0.154635;   fmcll[2]= 0.743349;  
          Fbb[3]= 1.22269;   Fcc[3]= 1.22269;   Fc[3]= 0.963292;   Fll[3]= 0.96057; 
             cafac[3]= 0.879497;  
             winpretag[3]= 109990 ;   wintag[3]= 12480.2 ;   woutpretag[3]= 96735.8;   wouttag[3]= 11788.8;  
           fbb[3]= 0.0602051;   fcc[3]= 0.12177;   fc[3]= 0.147895;   fll[3]= 0.670131;  
           fmcbb[3]= 0.0492397;   fmccc[3]= 0.0995912;   fmcc[3]= 0.153531;   fmcll[3]= 0.697639;  
          Fbb[4]= 1.20968;   Fcc[4]= 1.20968;   Fc[4]= 0.953042;   Fll[4]= 0.950349; 
             cafac[4]= 0.948247;  
             winpretag[4]= 23843.9 ;   wintag[4]= 3554.12 ;   woutpretag[4]= 22609.9;   wouttag[4]= 3656.07;  
           fbb[4]= 0.0855787;   fcc[4]= 0.144226;   fc[4]= 0.13628;   fll[4]= 0.633916;  
           fmcbb[4]= 0.0707446;   fmccc[4]= 0.119226;   fmcc[4]= 0.142994;   fmcll[4]= 0.667035;  
          Fbb[5]= 1.19539;   Fcc[5]= 1.19539;   Fc[5]= 0.941778;   Fll[5]= 0.939117; 
             cafac[5]= 0.810657;  
             winpretag[5]= 6824.64 ;   wintag[5]= 1318.58 ;   woutpretag[5]= 5532.44;   wouttag[5]= 1162.07;  
           fbb[5]= 0.113228;   fcc[5]= 0.169145;   fc[5]= 0.122876;   fll[5]= 0.59475;  
           fmcbb[5]= 0.0947208;   fmccc[5]= 0.141498;   fmcc[5]= 0.130472;   fmcll[5]= 0.633308;  
          Fbb[6]= 1.21912;   Fcc[6]= 1.21912;   Fc[6]= 0.960476;   Fll[6]= 0.957762; 
             cafac[6]= 0.887899;  
             winpretag[6]= 140658 ;   wintag[6]= 17352.9 ;   woutpretag[6]= 124891;   wouttag[6]= 16614.8;  
           fbb[6]= 0.0671636;   fcc[6]= 0.12795;   fc[6]= 0.144672;   fll[6]= 0.660214;  
           fmcbb[6]= 0.0550918;   fmccc[6]= 0.104953;   fmcc[6]= 0.150626;   fmcll[6]= 0.689329;  
          Fbb[7]= 1.20647;   Fcc[7]= 1.20647;   Fc[7]= 0.950512;   Fll[7]= 0.947826; 
             cafac[7]= 0.915367;  
             winpretag[7]= 30668.5 ;   wintag[7]= 4872.7 ;   woutpretag[7]= 28072.9;   wouttag[7]= 4844.27;  
           fbb[7]= 0.0917885;   fcc[7]= 0.149823;   fc[7]= 0.133269;   fll[7]= 0.625119;  
           fmcbb[7]= 0.07608;   fmccc[7]= 0.124183;   fmcc[7]= 0.140208;   fmcll[7]= 0.65953;  
     }  
     else if ( idsys==71   || sysname=="ptjmin_one" ) { 
          Fbb[1]= 1.25366;   Fcc[1]= 1.25366;   Fc[1]= 0.988576;   Fll[1]= 0.987515; 
             cafac[1]= 1.05183;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29253e+06;   wouttag[1]= 83570.2;  
           fbb[1]= 0.014843;   fcc[1]= 0.0433004;   fc[1]= 0.131814;   fll[1]= 0.810042;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.23527;   Fcc[2]= 1.23527;   Fc[2]= 0.974069;   Fll[2]= 0.973023; 
             cafac[2]= 0.972841;  
             winpretag[2]= 494445 ;   wintag[2]= 36811.8 ;   woutpretag[2]= 481017;   wouttag[2]= 37902.8;  
           fbb[2]= 0.037106;   fcc[2]= 0.089204;   fc[2]= 0.150566;   fll[2]= 0.723124;  
           fmcbb[2]= 0.0300389;   fmccc[2]= 0.0722143;   fmcc[2]= 0.154574;   fmcll[2]= 0.743173;  
          Fbb[3]= 1.22042;   Fcc[3]= 1.22042;   Fc[3]= 0.962364;   Fll[3]= 0.96133; 
             cafac[3]= 0.887383;  
             winpretag[3]= 109166 ;   wintag[3]= 12379.5 ;   woutpretag[3]= 96872.1;   wouttag[3]= 11785.5;  
           fbb[3]= 0.0599895;   fcc[3]= 0.12141;   fc[3]= 0.148028;   fll[3]= 0.670573;  
           fmcbb[3]= 0.0491546;   fmccc[3]= 0.0994816;   fmcc[3]= 0.153817;   fmcll[3]= 0.697547;  
          Fbb[4]= 1.20754;   Fcc[4]= 1.20754;   Fc[4]= 0.952204;   Fll[4]= 0.951181; 
             cafac[4]= 0.948974;  
             winpretag[4]= 23822.4 ;   wintag[4]= 3550.46 ;   woutpretag[4]= 22606.8;   wouttag[4]= 3651.5;  
           fbb[4]= 0.0852975;   fcc[4]= 0.143967;   fc[4]= 0.136513;   fll[4]= 0.634223;  
           fmcbb[4]= 0.0706375;   fmccc[4]= 0.119223;   fmcc[4]= 0.143365;   fmcll[4]= 0.666774;  
          Fbb[5]= 1.19414;   Fcc[5]= 1.19414;   Fc[5]= 0.941636;   Fll[5]= 0.940625; 
             cafac[5]= 0.90939;  
             winpretag[5]= 6144.54 ;   wintag[5]= 1186.52 ;   woutpretag[5]= 5587.78;   wouttag[5]= 1170.86;  
           fbb[5]= 0.112463;   fcc[5]= 0.166577;   fc[5]= 0.126044;   fll[5]= 0.594916;  
           fmcbb[5]= 0.0941793;   fmccc[5]= 0.139495;   fmcc[5]= 0.133857;   fmcll[5]= 0.632469;  
          Fbb[6]= 1.21702;   Fcc[6]= 1.21702;   Fc[6]= 0.959678;   Fll[6]= 0.958647; 
             cafac[6]= 0.899738;  
             winpretag[6]= 139133 ;   wintag[6]= 17116.5 ;   woutpretag[6]= 125183;   wouttag[6]= 16586.1;  
           fbb[6]= 0.0667185;   fcc[6]= 0.127335;   fc[6]= 0.145051;   fll[6]= 0.660895;  
           fmcbb[6]= 0.0548214;   fmccc[6]= 0.104629;   fmcc[6]= 0.151146;   fmcll[6]= 0.689404;  
          Fbb[7]= 1.20477;   Fcc[7]= 1.20477;   Fc[7]= 0.950018;   Fll[7]= 0.948998; 
             cafac[7]= 0.940224;  
             winpretag[7]= 29966.9 ;   wintag[7]= 4736.98 ;   woutpretag[7]= 28175.6;   wouttag[7]= 4830.86;  
           fbb[7]= 0.0909172;   fcc[7]= 0.148644;   fc[7]= 0.134347;   fll[7]= 0.626092;  
           fmcbb[7]= 0.0754646;   fmccc[7]= 0.12338;   fmcc[7]= 0.141416;   fmcll[7]= 0.65974;  
     }  
     else if ( idsys==69   || sysname=="powhe_one" ) { 
          Fbb[1]= 1.28866;   Fcc[1]= 1.28866;   Fc[1]= 1.00134;   Fll[1]= 0.983461; 
             cafac[1]= 1.0541;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29748e+06;   wouttag[1]= 84677.1;  
           fbb[1]= 0.0152574;   fcc[1]= 0.0445091;   fc[1]= 0.133517;   fll[1]= 0.806717;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.26672;   Fcc[2]= 1.26672;   Fc[2]= 0.984292;   Fll[2]= 0.966712; 
             cafac[2]= 0.975475;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 482557;   wouttag[2]= 38393.3;  
           fbb[2]= 0.0378869;   fcc[2]= 0.0911796;   fc[2]= 0.152301;   fll[2]= 0.718633;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.2493;   Fcc[3]= 1.2493;   Fc[3]= 0.970755;   Fll[3]= 0.953418; 
             cafac[3]= 0.888389;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 97204.6;   wouttag[3]= 11958.8;  
           fbb[3]= 0.061296;   fcc[3]= 0.124111;   fc[3]= 0.149553;   fll[3]= 0.66504;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.23458;   Fcc[4]= 1.23458;   Fc[4]= 0.959322;   Fll[4]= 0.942189; 
             cafac[4]= 0.952716;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 22754.6;   wouttag[4]= 3709.43;  
           fbb[4]= 0.0870382;   fcc[4]= 0.146647;   fc[4]= 0.138106;   fll[4]= 0.628209;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.21894;   Fcc[5]= 1.21894;   Fc[5]= 0.94717;   Fll[5]= 0.930254; 
             cafac[5]= 0.904758;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5638.82;   wouttag[5]= 1193.62;  
           fbb[5]= 0.114298;   fcc[5]= 0.170648;   fc[5]= 0.126578;   fll[5]= 0.588476;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.24537;   Fcc[6]= 1.24537;   Fc[6]= 0.967705;   Fll[6]= 0.950422; 
             cafac[6]= 0.900906;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 125706;   wouttag[6]= 16841.9;  
           fbb[6]= 0.0681595;   fcc[6]= 0.130126;   fc[6]= 0.146528;   fll[6]= 0.655186;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.23131;   Fcc[7]= 1.23131;   Fc[7]= 0.956781;   Fll[7]= 0.939694; 
             cafac[7]= 0.942053;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 28371.2;   wouttag[7]= 4913.64;  
           fbb[7]= 0.0927368;   fcc[7]= 0.151664;   fc[7]= 0.135696;   fll[7]= 0.619903;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==70   || sysname=="powpy_one" ) { 
          Fbb[1]= 1.24542;   Fcc[1]= 1.24542;   Fc[1]= 1.00467;   Fll[1]= 0.985365; 
             cafac[1]= 1.0538;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77366.1 ;   woutpretag[1]= 2.29682e+06;   wouttag[1]= 84337.8;  
           fbb[1]= 0.0147454;   fcc[1]= 0.0430157;   fc[1]= 0.13396;   fll[1]= 0.808279;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.2272;   Fcc[2]= 1.2272;   Fc[2]= 0.989967;   Fll[2]= 0.970947; 
             cafac[2]= 0.976383;  
             winpretag[2]= 494690 ;   wintag[2]= 36747.4 ;   woutpretag[2]= 483007;   wouttag[2]= 38095.9;  
           fbb[2]= 0.036705;   fcc[2]= 0.0883352;   fc[2]= 0.153179;   fll[2]= 0.721781;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.21276;   Fcc[3]= 1.21276;   Fc[3]= 0.978317;   Fll[3]= 0.959522; 
             cafac[3]= 0.900693;  
             winpretag[3]= 109417 ;   wintag[3]= 12395.7 ;   woutpretag[3]= 98550.9;   wouttag[3]= 11993.9;  
           fbb[3]= 0.0595033;   fcc[3]= 0.120482;   fc[3]= 0.150718;   fll[3]= 0.669297;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.20056;   Fcc[4]= 1.20056;   Fc[4]= 0.968477;   Fll[4]= 0.94987; 
             cafac[4]= 0.97721;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3548.73 ;   woutpretag[4]= 23339.6;   wouttag[4]= 3758.36;  
           fbb[4]= 0.0846397;   fcc[4]= 0.142606;   fc[4]= 0.139424;   fll[4]= 0.633331;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.18754;   Fcc[5]= 1.18754;   Fc[5]= 0.957978;   Fll[5]= 0.939573; 
             cafac[5]= 0.950564;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1201.29 ;   woutpretag[5]= 5924.3;   wouttag[5]= 1238.61;  
           fbb[5]= 0.111354;   fcc[5]= 0.166252;   fc[5]= 0.128023;   fll[5]= 0.594371;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.20951;   Fcc[6]= 1.20951;   Fc[6]= 0.975695;   Fll[6]= 0.95695; 
             cafac[6]= 0.917282;  
             winpretag[6]= 139533 ;   wintag[6]= 17145.7 ;   woutpretag[6]= 127991;   wouttag[6]= 16953.6;  
           fbb[6]= 0.0661968;   fcc[6]= 0.126379;   fc[6]= 0.147738;   fll[6]= 0.659686;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.19784;   Fcc[7]= 1.19784;   Fc[7]= 0.966285;   Fll[7]= 0.947721; 
             cafac[7]= 0.971241;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4750.02 ;   woutpretag[7]= 29250.2;   wouttag[7]= 5003.48;  
           fbb[7]= 0.090216;   fcc[7]= 0.147542;   fc[7]= 0.137044;   fll[7]= 0.625198;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==15   || sysname=="btagb_up" ) { 
          Fbb[1]= 1.16533;   Fcc[1]= 1.16533;   Fc[1]= 0.964772;   Fll[1]= 0.996379; 
             cafac[1]= 1.04772;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 77820.2 ;   woutpretag[1]= 2.28356e+06;   wouttag[1]= 81816;  
           fbb[1]= 0.0137971;   fcc[1]= 0.0402492;   fc[1]= 0.12864;   fll[1]= 0.817313;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.15527;   Fcc[2]= 1.15527;   Fc[2]= 0.956449;   Fll[2]= 0.987783; 
             cafac[2]= 0.969858;  
             winpretag[2]= 494690 ;   wintag[2]= 37163.7 ;   woutpretag[2]= 479779;   wouttag[2]= 37117.1;  
           fbb[2]= 0.0345537;   fcc[2]= 0.0831578;   fc[2]= 0.147992;   fll[2]= 0.734296;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.14632;   Fcc[3]= 1.14632;   Fc[3]= 0.949034;   Fll[3]= 0.980125; 
             cafac[3]= 0.88358;  
             winpretag[3]= 109417 ;   wintag[3]= 12574.2 ;   woutpretag[3]= 96678.4;   wouttag[3]= 11581.2;  
           fbb[3]= 0.0562434;   fcc[3]= 0.113881;   fc[3]= 0.146207;   fll[3]= 0.683669;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.13823;   Fcc[4]= 1.13823;   Fc[4]= 0.942338;   Fll[4]= 0.973209; 
             cafac[4]= 0.947688;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3606 ;   woutpretag[4]= 22634.5;   wouttag[4]= 3592.85;  
           fbb[4]= 0.0802453;   fcc[4]= 0.135202;   fc[4]= 0.135661;   fll[4]= 0.648892;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.12958;   Fcc[5]= 1.12958;   Fc[5]= 0.935175;   Fll[5]= 0.965812; 
             cafac[5]= 0.901395;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1223.19 ;   woutpretag[5]= 5617.86;   wouttag[5]= 1161.61;  
           fbb[5]= 0.105918;   fcc[5]= 0.158137;   fc[5]= 0.124975;   fll[5]= 0.61097;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.14417;   Fcc[6]= 1.14417;   Fc[6]= 0.947255;   Fll[6]= 0.978288; 
             cafac[6]= 0.896268;  
             winpretag[6]= 139533 ;   wintag[6]= 17403.3 ;   woutpretag[6]= 125059;   wouttag[6]= 16311;  
           fbb[6]= 0.0626207;   fcc[6]= 0.119552;   fc[6]= 0.143432;   fll[6]= 0.674396;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.13643;   Fcc[7]= 1.13643;   Fc[7]= 0.940846;   Fll[7]= 0.971669; 
             cafac[7]= 0.937465;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4829.19 ;   woutpretag[7]= 28233;   wouttag[7]= 4763.89;  
           fbb[7]= 0.0855904;   fcc[7]= 0.139977;   fc[7]= 0.133436;   fll[7]= 0.640997;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==16   || sysname=="btagb_down" ) { 
          Fbb[1]= 1.35887;   Fcc[1]= 1.35887;   Fc[1]= 1.01574;   Fll[1]= 0.977152; 
             cafac[1]= 1.05655;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 76911.9 ;   woutpretag[1]= 2.30281e+06;   wouttag[1]= 85643.3;  
           fbb[1]= 0.0160886;   fcc[1]= 0.046934;   fc[1]= 0.135436;   fll[1]= 0.801542;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.3296;   Fcc[2]= 1.3296;   Fc[2]= 0.993856;   Fll[2]= 0.956103; 
             cafac[2]= 0.977323;  
             winpretag[2]= 494690 ;   wintag[2]= 36326.5 ;   woutpretag[2]= 483472;   wouttag[2]= 38714.6;  
           fbb[2]= 0.0397676;   fcc[2]= 0.0957059;   fc[2]= 0.15378;   fll[2]= 0.710746;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.30692;   Fcc[3]= 1.30692;   Fc[3]= 0.976907;   Fll[3]= 0.939798; 
             cafac[3]= 0.889901;  
             winpretag[3]= 109417 ;   wintag[3]= 12213.7 ;   woutpretag[3]= 97370;   wouttag[3]= 12031;  
           fbb[3]= 0.0641234;   fcc[3]= 0.129836;   fc[3]= 0.150501;   fll[3]= 0.65554;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.28808;   Fcc[4]= 1.28808;   Fc[4]= 0.96282;   Fll[4]= 0.926246; 
             cafac[4]= 0.954282;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3489.98 ;   woutpretag[4]= 22792;   wouttag[4]= 3728.13;  
           fbb[4]= 0.0908097;   fcc[4]= 0.153001;   fc[4]= 0.13861;   fll[4]= 0.617579;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.26814;   Fcc[5]= 1.26814;   Fc[5]= 0.947921;   Fll[5]= 0.911914; 
             cafac[5]= 0.905477;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1178.13 ;   woutpretag[5]= 5643.31;   wouttag[5]= 1194.61;  
           fbb[5]= 0.118912;   fcc[5]= 0.177536;   fc[5]= 0.126679;   fll[5]= 0.576874;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.30188;   Fcc[6]= 1.30188;   Fc[6]= 0.973141;   Fll[6]= 0.936175; 
             cafac[6]= 0.902309;  
             winpretag[6]= 139533 ;   wintag[6]= 16881.8 ;   woutpretag[6]= 125902;   wouttag[6]= 16936.7;  
           fbb[6]= 0.0712526;   fcc[6]= 0.136031;   fc[6]= 0.147351;   fll[6]= 0.645365;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.2839;   Fcc[7]= 1.2839;   Fc[7]= 0.959699;   Fll[7]= 0.923244; 
             cafac[7]= 0.943391;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4668.11 ;   woutpretag[7]= 28411.4;   wouttag[7]= 4933.73;  
           fbb[7]= 0.0966975;   fcc[7]= 0.158142;   fc[7]= 0.13611;   fll[7]= 0.609051;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==85   || sysname=="btagc_up" ) { 
          Fbb[1]= 1.19243;   Fcc[1]= 1.19243;   Fc[1]= 0.859941;   Fll[1]= 1.01189; 
             cafac[1]= 1.02836;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 85533.5 ;   woutpretag[1]= 2.24137e+06;   wouttag[1]= 83169.6;  
           fbb[1]= 0.014118;   fcc[1]= 0.0411854;   fc[1]= 0.114662;   fll[1]= 0.830034;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.18441;   Fcc[2]= 1.18441;   Fc[2]= 0.854157;   Fll[2]= 1.00508; 
             cafac[2]= 0.948137;  
             winpretag[2]= 494690 ;   wintag[2]= 39954.9 ;   woutpretag[2]= 469034;   wouttag[2]= 37663.3;  
           fbb[2]= 0.0354252;   fcc[2]= 0.0852552;   fc[2]= 0.132165;   fll[2]= 0.747155;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.17449;   Fcc[3]= 1.17449;   Fc[3]= 0.847005;   Fll[3]= 0.996665; 
             cafac[3]= 0.864036;  
             winpretag[3]= 109417 ;   wintag[3]= 13368.3 ;   woutpretag[3]= 94539.9;   wouttag[3]= 11775.4;  
           fbb[3]= 0.0576259;   fcc[3]= 0.11668;   fc[3]= 0.130488;   fll[3]= 0.695206;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.16427;   Fcc[4]= 1.16427;   Fc[4]= 0.839633;   Fll[4]= 0.987991; 
             cafac[4]= 0.927248;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3801.84 ;   woutpretag[4]= 22146.3;   wouttag[4]= 3652.44;  
           fbb[4]= 0.0820814;   fcc[4]= 0.138295;   fc[4]= 0.120876;   fll[4]= 0.658748;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.15346;   Fcc[5]= 1.15346;   Fc[5]= 0.831837;   Fll[5]= 0.978817; 
             cafac[5]= 0.882103;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1280.92 ;   woutpretag[5]= 5497.63;   wouttag[5]= 1178.41;  
           fbb[5]= 0.108158;   fcc[5]= 0.16148;   fc[5]= 0.111165;   fll[5]= 0.619196;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.17178;   Fcc[6]= 1.17178;   Fc[6]= 0.845047;   Fll[6]= 0.994361; 
             cafac[6]= 0.876525;  
             winpretag[6]= 139533 ;   wintag[6]= 18451.1 ;   woutpretag[6]= 122304;   wouttag[6]= 16582.4;  
           fbb[6]= 0.0641318;   fcc[6]= 0.122437;   fc[6]= 0.127955;   fll[6]= 0.685476;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.16202;   Fcc[7]= 1.16202;   Fc[7]= 0.838008;   Fll[7]= 0.986078; 
             cafac[7]= 0.917261;  
             winpretag[7]= 30116.3 ;   wintag[7]= 5082.76 ;   woutpretag[7]= 27624.5;   wouttag[7]= 4840.57;  
           fbb[7]= 0.0875179;   fcc[7]= 0.143129;   fc[7]= 0.118851;   fll[7]= 0.650502;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==86   || sysname=="btagc_down" ) { 
          Fbb[1]= 1.31769;   Fcc[1]= 1.31769;   Fc[1]= 1.15854;   Fll[1]= 0.956267; 
             cafac[1]= 1.08462;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 69198.7 ;   woutpretag[1]= 2.36398e+06;   wouttag[1]= 84126.3;  
           fbb[1]= 0.015601;   fcc[1]= 0.0455117;   fc[1]= 0.154477;   fll[1]= 0.784411;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.28631;   Fcc[2]= 1.28631;   Fc[2]= 1.13095;   Fll[2]= 0.933499; 
             cafac[2]= 1.00854;  
             winpretag[2]= 494690 ;   wintag[2]= 33529.2 ;   woutpretag[2]= 498916;   wouttag[2]= 38143.4;  
           fbb[2]= 0.0384731;   fcc[2]= 0.0925903;   fc[2]= 0.174994;   fll[2]= 0.693943;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.26571;   Fcc[3]= 1.26571;   Fc[3]= 1.11284;   Fll[3]= 0.918545; 
             cafac[3]= 0.917736;  
             winpretag[3]= 109417 ;   wintag[3]= 11417.7 ;   woutpretag[3]= 100416;   wouttag[3]= 11827.8;  
           fbb[3]= 0.0621013;   fcc[3]= 0.125742;   fc[3]= 0.171442;   fll[3]= 0.640715;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.25042;   Fcc[4]= 1.25042;   Fc[4]= 1.09939;   Fll[4]= 0.907449; 
             cafac[4]= 0.983257;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3293.16 ;   woutpretag[4]= 23484;   wouttag[4]= 3665.99;  
           fbb[4]= 0.0881547;   fcc[4]= 0.148528;   fc[4]= 0.158271;   fll[4]= 0.605046;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.23404;   Fcc[5]= 1.23404;   Fc[5]= 1.08499;   Fll[5]= 0.895561; 
             cafac[5]= 0.932625;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1120.61 ;   woutpretag[5]= 5812.5;   wouttag[5]= 1178.59;  
           fbb[5]= 0.115714;   fcc[5]= 0.172761;   fc[5]= 0.144996;   fll[5]= 0.566529;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.26162;   Fcc[6]= 1.26162;   Fc[6]= 1.10924;   Fll[6]= 0.915579; 
             cafac[6]= 0.930383;  
             winpretag[6]= 139533 ;   wintag[6]= 15831.4 ;   woutpretag[6]= 129819;   wouttag[6]= 16654;  
           fbb[6]= 0.069049;   fcc[6]= 0.131824;   fc[6]= 0.16796;   fll[6]= 0.631167;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.24699;   Fcc[7]= 1.24699;   Fc[7]= 1.09638;   Fll[7]= 0.904963; 
             cafac[7]= 0.971985;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4413.77 ;   woutpretag[7]= 29272.6;   wouttag[7]= 4855.26;  
           fbb[7]= 0.0939178;   fcc[7]= 0.153595;   fc[7]= 0.155495;   fll[7]= 0.596991;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==17   || sysname=="bmtag_up" ) { 
          Fbb[1]= 1.02197;   Fcc[1]= 1.02197;   Fc[1]= 0.96255;   Fll[1]= 1.00485; 
             cafac[1]= 1.04788;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 82195.6 ;   woutpretag[1]= 2.28391e+06;   wouttag[1]= 84763.5;  
           fbb[1]= 0.0120998;   fcc[1]= 0.0352979;   fc[1]= 0.128344;   fll[1]= 0.824258;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.02192;   Fcc[2]= 1.02192;   Fc[2]= 0.962506;   Fll[2]= 1.0048; 
             cafac[2]= 0.971686;  
             winpretag[2]= 494690 ;   wintag[2]= 39102.2 ;   woutpretag[2]= 480683;   wouttag[2]= 37728;  
           fbb[2]= 0.0305653;   fcc[2]= 0.0735592;   fc[2]= 0.14893;   fll[2]= 0.746946;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.02108;   Fcc[3]= 1.02108;   Fc[3]= 0.961713;   Fll[3]= 1.00397; 
             cafac[3]= 0.885477;  
             winpretag[3]= 109417 ;   wintag[3]= 13172.7 ;   woutpretag[3]= 96886;   wouttag[3]= 11635.9;  
           fbb[3]= 0.0500988;   fcc[3]= 0.101439;   fc[3]= 0.14816;   fll[3]= 0.700302;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.01993;   Fcc[4]= 1.01993;   Fc[4]= 0.960631;   Fll[4]= 1.00284; 
             cafac[4]= 0.949668;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3775.27 ;   woutpretag[4]= 22681.8;   wouttag[4]= 3587.33;  
           fbb[4]= 0.0719055;   fcc[4]= 0.12115;   fc[4]= 0.138295;   fll[4]= 0.66865;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.01871;   Fcc[5]= 1.01871;   Fc[5]= 0.959484;   Fll[5]= 1.00164; 
             cafac[5]= 0.90483;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1283.62 ;   woutpretag[5]= 5639.27;   wouttag[5]= 1163.64;  
           fbb[5]= 0.0955232;   fcc[5]= 0.142617;   fc[5]= 0.128224;   fll[5]= 0.633637;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.02078;   Fcc[6]= 1.02078;   Fc[6]= 0.961428;   Fll[6]= 1.00367; 
             cafac[6]= 0.898428;  
             winpretag[6]= 139533 ;   wintag[6]= 18231.6 ;   woutpretag[6]= 125360;   wouttag[6]= 16357;  
           fbb[6]= 0.0558676;   fcc[6]= 0.106659;   fc[6]= 0.145578;   fll[6]= 0.691896;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.01968;   Fcc[7]= 1.01968;   Fc[7]= 0.960393;   Fll[7]= 1.00259; 
             cafac[7]= 0.93985;  
             winpretag[7]= 30116.3 ;   wintag[7]= 5058.9 ;   woutpretag[7]= 28304.8;   wouttag[7]= 4759.19;  
           fbb[7]= 0.0767977;   fcc[7]= 0.125597;   fc[7]= 0.136209;   fll[7]= 0.661397;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
     }  
     else if ( idsys==18   || sysname=="bmtag_down" ) { 
          Fbb[1]= 1.48649;   Fcc[1]= 1.48649;   Fc[1]= 1.01833;   Fll[1]= 0.969515; 
             cafac[1]= 1.05652;  
             winpretag[1]= 2.17955e+06 ;   wintag[1]= 72536.5 ;   woutpretag[1]= 2.30275e+06;   wouttag[1]= 82666.8;  
           fbb[1]= 0.0175996;   fcc[1]= 0.0513418;   fc[1]= 0.135781;   fll[1]= 0.795277;  
           fmcbb[1]= 0.0118397;   fmccc[1]= 0.034539;   fmcc[1]= 0.133338;   fmcll[1]= 0.820284;  
          Fbb[2]= 1.44355;   Fcc[2]= 1.44355;   Fc[2]= 0.988916;   Fll[2]= 0.941512; 
             cafac[2]= 0.975794;  
             winpretag[2]= 494690 ;   wintag[2]= 34422.8 ;   woutpretag[2]= 482715;   wouttag[2]= 38077.2;  
           fbb[2]= 0.043176;   fcc[2]= 0.103908;   fc[2]= 0.153016;   fll[2]= 0.699899;  
           fmcbb[2]= 0.0299095;   fmccc[2]= 0.0719811;   fmcc[2]= 0.154731;   fmcll[2]= 0.743378;  
          Fbb[3]= 1.41065;   Fcc[3]= 1.41065;   Fc[3]= 0.966378;   Fll[3]= 0.920054; 
             cafac[3]= 0.888307;  
             winpretag[3]= 109417 ;   wintag[3]= 11619.4 ;   woutpretag[3]= 97195.7;   wouttag[3]= 11945.5;  
           fbb[3]= 0.0692128;   fcc[3]= 0.140141;   fc[3]= 0.148878;   fll[3]= 0.641767;  
           fmcbb[3]= 0.0490644;   fmccc[3]= 0.0993451;   fmcc[3]= 0.154058;   fmcll[3]= 0.697532;  
          Fbb[4]= 1.38356;   Fcc[4]= 1.38356;   Fc[4]= 0.947815;   Fll[4]= 0.902381; 
             cafac[4]= 0.952617;  
             winpretag[4]= 23883.9 ;   wintag[4]= 3325.1 ;   woutpretag[4]= 22752.2;   wouttag[4]= 3723.93;  
           fbb[4]= 0.0975409;   fcc[4]= 0.164342;   fc[4]= 0.13645;   fll[4]= 0.601667;  
           fmcbb[4]= 0.0705002;   fmccc[4]= 0.118783;   fmcc[4]= 0.143962;   fmcll[4]= 0.666755;  
          Fbb[5]= 1.35519;   Fcc[5]= 1.35519;   Fc[5]= 0.92838;   Fll[5]= 0.883878; 
             cafac[5]= 0.902685;  
             winpretag[5]= 6232.41 ;   wintag[5]= 1120.28 ;   woutpretag[5]= 5625.9;   wouttag[5]= 1189.82;  
           fbb[5]= 0.127073;   fcc[5]= 0.189721;   fc[5]= 0.124067;   fll[5]= 0.559138;  
           fmcbb[5]= 0.0937683;   fmccc[5]= 0.139997;   fmcc[5]= 0.133638;   fmcll[5]= 0.632597;  
          Fbb[6]= 1.40338;   Fcc[6]= 1.40338;   Fc[6]= 0.961397;   Fll[6]= 0.915312; 
             cafac[6]= 0.900499;  
             winpretag[6]= 139533 ;   wintag[6]= 16064.8 ;   woutpretag[6]= 125649;   wouttag[6]= 16846.8;  
           fbb[6]= 0.0768076;   fcc[6]= 0.146636;   fc[6]= 0.145573;   fll[6]= 0.630983;  
           fmcbb[6]= 0.0547304;   fmccc[6]= 0.104488;   fmcc[6]= 0.151418;   fmcll[6]= 0.689364;  
          Fbb[7]= 1.37759;   Fcc[7]= 1.37759;   Fc[7]= 0.943726;   Fll[7]= 0.898489; 
             cafac[7]= 0.941403;  
             winpretag[7]= 30116.3 ;   wintag[7]= 4445.38 ;   woutpretag[7]= 28351.6;   wouttag[7]= 4925.58;  
           fbb[7]= 0.103754;   fcc[7]= 0.169681;   fc[7]= 0.133845;   fll[7]= 0.592721;  
           fmcbb[7]= 0.0753154;   fmccc[7]= 0.123173;   fmcc[7]= 0.141826;   fmcll[7]= 0.659686;  
        }  
       else { cout << " WARNING: Ffactors request for unknown variation. Set to Nominal! " << idsys << "  " << sysname <<  endl; }

     return;


}
