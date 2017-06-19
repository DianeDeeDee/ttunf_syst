#ifndef muon_SF_EPS_h
#define muon_SF_EPS_h

/******************************************************************************
 * muon_SF_EPS.h                                                              *
 *                                                                            *
 * Simple functions which return data/MC scale factors and trigger effs       *
 * given a muon's (or two muons') eta and phi (defined b/w -pi and pi).       *
 * Also functions for corresponding uncertainties.                            *
 *                     .                                                      *
 * double mu_trigger_SF(double eta, double phi)                               * 
 * double mu_trigger_SF_err(double eta, double phi, bool isUp)                * 
 * double mu_trigger_eff(double eta, double phi)                              * 
 * double mu_trigger_eff_err(double eta, double phi, bool isUp)               * 
 *                                                                            * 
 * double mu_ID_SF()                                    	              *
 * double mu_ID_SF_err()                                          	      *
 * double mu_ID_SF_staterr()                                                  *
 * double mu_ID_SF_systerr()                                   	              *
 *                                                                            *
 * History                                                                    *
 *       08 Dec 2010 -- created by S. Caughron                                *   
 *       11 Feb 2011 -- Finalization for Moriond                              *   
 *       15 Mar 2011 -- removing reco SF (provided by MCP now)                *   
 *       09 May 2011 -- update for pLHC                                       *   
 *       14 Jun 2011 -- update for EPS                                        *
 *       23 Jun 2011 -- added trigger efficiencies for EPS                    *
 *       30 Jun 2011 -- updated trigger eff uncertianties for EPS             *
 *       13 Jul 2011 -- put new eff maps for LP (improved binning, B to G5)   *
 *****************************************************************************/

#include <iostream>
#include <cmath>
#include <cstdio>

// forward declarations
double mu_trigger_SF(double eta, double phi, double pt);
double mu_trigger_SF_err(double eta, double phi, double pt, bool isUp);
double mu_trigger_eff(double eta, double phi, double pt);
double mu_trigger_eff_err(double eta, double phi, double pt, bool isUp);

double mu_ID_SF();
double mu_ID_SF_err();
double mu_ID_SF_staterr();
double mu_ID_SF_systerr();
double check_phi_range(double phi);
void printWarning();

//these are deprecated
double mu_recoID_SF();
double mu_recoID_SF_err();
double mu_recoID_SF_staterr();
double mu_recoID_SF_systerr();


// For trigger scale factors 
double mu_trigger_SF(double eta, double phi, double pt)
{

    double mu_eta = eta;
    double mu_phi = check_phi_range(phi);
    double mu_pt = pt;
    int etaI=-1, phiI=-1;   

    const double pi = acos(-1.);

    const double EC_eta = 1.05;
    const double etabins[7] = {-1.05,-0.75,-0.45,-0.15,0.15,0.45,0.75};   // lower edges of eta bins
    const double phibins[3] = {-pi,0.,pi/2.};   

    const double EC_SF_pos = 1.0140;
    const double EC_SF_neg = 1.0097;
    const double trig_SFmatrix_lo[7][3] = {{1.0388,1.0660,1.0388},   // rows are eta, columns are phi
                              			   {1.0196,1.0660,1.0196},     // period B2 numbers
                              			   {0.9996,1.0660,0.9996},
                              			   {1.0315,1.0315,1.0315},
                              			   {1.0026,1.0026,1.0026},
                              			   {0.9070,1.0369,1.0369},
                              			   {0.9070,0.9946,0.9946}};
    const double trig_SFmatrix_mid[7][3] = {{1.0893,1.1525,1.0893},   // rows are eta, columns are phi
                              				{1.0767,1.1525,1.0767},     // period B2 numbers
                              				{1.0880,1.1525,1.0880},
                              				{1.1071,1.1071,1.1071},
                              				{1.0609,1.0609,1.0609},
                              				{0.9400,1.1205,1.1205},
                              				{0.9400,1.0684,1.0684}};
    const double trig_SFmatrix_hi[7][3] = {{1.3432,0.8027,1.3432},   // rows are eta, columns are phi
                              			   {1.1803,0.8027,1.1803},     // period B2 numbers
                              			   {0.9868,0.8027,0.9868},
                              			   {1.2126,1.2126,1.2126},
                              			   {0.8520,0.8520,0.8520},
                              			   {1.5726,1.2885,1.2885},
                              			   {1.5726,1.3099,1.3099}};
    
    if ( mu_eta > EC_eta )    // check if in endcap
        return EC_SF_pos;       
	else if ( mu_eta < -EC_eta)
        return EC_SF_neg;       
    else {                          // if in barrel, look up appropriate SF 
        
        for (int i=2; i>=0; i--){    // find phi index
            if ( mu_phi > phibins[i] ) {
                phiI = i;
                break;
            }
        }
        for (int i=6; i>=0; i--){    // find eta index
            if ( mu_eta > etabins[i] ) {
                etaI = i;
                break;
            }
        }

		if (mu_pt < 60.)
        	return trig_SFmatrix_lo[etaI][phiI];   
		else if (mu_pt < 120.)
        	return trig_SFmatrix_mid[etaI][phiI];   
		else 
        	return trig_SFmatrix_hi[etaI][phiI];   

    } //else

    return 0;
}

// Trigger scale factor uncertainties 
double mu_trigger_SF_err(double eta, double phi, double pt, bool isUp) 
{

    double mu_eta = eta;
    double mu_phi = check_phi_range(phi);
    double mu_pt = pt;
    int etaI=-1, phiI=-1;   

    const double pi = acos(-1.);

    const double EC_eta = 1.05;
    const double etabins[7] = {-1.05,-0.75,-0.45,-0.15,0.15,0.45,0.75};   // lower edges of eta bins
    const double phibins[3] = {-pi,0.,pi/2.};   
    
    const double EC_SF_pos_up = 0.07000;
    const double EC_SF_pos_down = 0.07000;
    const double EC_SF_neg_up = 0.07000;
    const double EC_SF_neg_down = 0.07000;
    
    const double trig_SFmatrix_up_lo[7][3] = {{0.00644,0.00543,0.00644},   // rows are eta, columns are phi
					      {0.00577,0.00543,0.00577},   // positive errors
					      {0.07000,0.00543,0.07000},
					      {0.00609,0.00609,0.00609},
					      {0.00459,0.00459,0.00459},
					      {0.07000,0.07000,0.07000},
					      {0.07000,0.07000,0.07000}};

    const double trig_SFmatrix_down_lo[7][3] = {{0.00640,0.00528,0.00640},   // rows are eta, columns are phi
						{0.00609,0.00536,0.00609},   // negative errors
						{0.07000,0.00536,0.07000},
						{0.00615,0.00615,0.00615},
						{0.00473,0.00473,0.00473},
						{0.07000,0.07000,0.07000},
						{0.07000,0.07000,0.07000}};
    
    const double trig_SFmatrix_up_mid[7][3] = {{0.03443,0.02918,0.03443},   // rows are eta, columns are phi
					       {0.03345,0.02918,0.03345},   // negative errors
					       {0.03101,0.02918,0.03101},
					       {0.03677,0.03677,0.03677},
					       {0.07492,0.07492,0.07492},
					       {0.03053,0.03296,0.03296},
					       {0.03053,0.07998,0.07998}};

    const double trig_SFmatrix_down_mid[7][3] = {{0.03426,0.03131,0.03426},   // rows are eta, columns are phi
						 {0.03308,0.03131,0.03308},   // positive errors
						 {0.03456,0.03131,0.03456},
						 {0.03503,0.03503,0.03503},
						 {0.07483,0.07483,0.07483},
						 {0.02877,0.03211,0.03211},
						 {0.02877,0.08131,0.08131}};
    
    const double trig_SFmatrix_up_hi[7][3] = {{0.1446,0.2084,0.1446},   // rows are eta, columns are phi
					      {0.1480,0.2084,0.1480},   // positive errors
					      {0.1622,0.2084,0.1621},
					      {0.1445,0.1445,0.1445},
					      {0.1512,0.1512,0.1512},
					      {0.1487,0.1704,0.1704},
					      {0.1487,0.1763,0.1763}};

    const double trig_SFmatrix_down_hi[7][3] = {{0.3481,0.1229,0.3481},   // rows are eta, columns are phi
						{0.1426,0.1229,0.1426},   // negative errors
						{0.1677,0.1229,0.1677},
						{0.2142,0.2142,0.2142},
						{0.1490,0.1490,0.1490},
						{0.2493,0.2549,0.2549},
						{0.2493,0.1763,0.1763}};
    
    
    if ( mu_eta > EC_eta && isUp)    // check if in endcap
        return EC_SF_pos_up;       
    else if ( mu_eta > EC_eta && !isUp)    
        return EC_SF_pos_down;       
	else if ( mu_eta < -EC_eta && isUp)
        return EC_SF_neg_up;       
	else if ( mu_eta < -EC_eta && !isUp)
        return EC_SF_neg_down;       
    else {                          // if in barrel, look up appropriate SF 
        
        for (int i=2; i>=0; i--){    // find phi index
            if ( mu_phi > phibins[i] ) {
                phiI = i;
                break;
            }
        }
        for (int i=6; i>=0; i--){    // find eta index
            if ( mu_eta > etabins[i] ) {
                etaI = i;
                break;
            }
        }

		if (mu_pt < 60.) {
			if (isUp)
        		return trig_SFmatrix_up_lo[etaI][phiI];   
			else
        		return trig_SFmatrix_down_lo[etaI][phiI];   
		} else if (mu_pt < 120.) {
			if (isUp)
        		return trig_SFmatrix_up_mid[etaI][phiI];   
			else
        		return trig_SFmatrix_down_mid[etaI][phiI];   
		} else {
			if (isUp)
        		return trig_SFmatrix_up_hi[etaI][phiI];   
			else
        		return trig_SFmatrix_down_hi[etaI][phiI];   
		}

    } //else

    return 0;
}

// For trigger efficiencies as measured in data
double mu_trigger_eff(double eta, double phi, double pt) //pt in GeV
{
  
  double mu_eta = eta;
  double mu_phi = check_phi_range(phi);
  double mu_pt = pt;
  int etaI=-1, phiI=-1;   
  
  const double etabins[10] = {-2.625, -2.4, -1.3125, -1.05 , -0.2625, 0., 0.2625, 1.05, 1.3125, 2.4} ;
  const double phibins[7]  = {-3.5  , -2.1, -1.8   , -1.325, -1.025 , 0., 1.7 } ;
  
  // central values for eff
  // rows: eta, columns: phi
  const double trig_effmatrix[10][7] = { {0.43969,0.40218,0.34967,0.38809,0.43872,0.39262,0.46143},
					 {0.89834,0.87006,0.90091,0.88300,0.89272,0.89046,0.90262}, 
					 {0.82306,0.79439,0.87077,0.81320,0.76661,0.82105,0.85080}, 
					 {0.76700,0.44346,0.60514,0.53356,0.76977,0.80560,0.75631}, 
					 {0.75339,0.72363,0.81760,0.70446,0.77366,0.80438,0.76362}, 
					 {0.79608,0.53052,0.47011,0.59173,0.84150,0.84608,0.76750}, 
					 {0.73785,0.43620,0.57768,0.49414,0.73457,0.76332,0.77936}, 
					 {0.82384,0.82609,0.88498,0.83331,0.78040,0.83383,0.86723}, 
					 {0.92032,0.87821,0.91550,0.89843,0.89815,0.88828,0.90471},
					 {0.44864,0.52562,0.44621,0.41199,0.47011,0.45515,0.47574} };

  for (int i=6; i>=0; i--){    // find phi index
    if ( mu_phi > phibins[i] ) {
      phiI = i;
      break;
    }
  }
  for (int i=9; i>=0; i--){    // find eta index
    if ( mu_eta > etabins[i] ) {
      etaI = i;
      break;
    }
  }
  
  return trig_effmatrix[etaI][phiI];   
}

// Trigger efficiency uncertainties
//pt should be in GeV, isUp=true is +1 sigma unc, false is -1 sigma unc
double mu_trigger_eff_err(double eta, double phi, double pt, bool isUp) 
{

  double mu_eta = eta;
  double mu_phi = check_phi_range(phi);
  double mu_pt = pt;
  int etaI=-1, phiI=-1;   
  
  const double etabins[10] = {-2.625, -2.4, -1.3125, -1.05 , -0.2625, 0., 0.2625, 1.05, 1.3125, 2.4} ;
  const double phibins[7]  = {-3.5  , -2.1, -1.8   , -1.325, -1.025 , 0., 1.7 } ;
  
  //uncertainties
  // rows: eta, columns: phi
  const double trig_effmatrix_err[10][7] = { {0.0206373,0.0236883,0.027266,0.0241146,0.0219176,0.0161734,0.0160637}, 
					     {0.0030955,0.0042006,0.004556,0.0039498,0.0036496,0.0023653,0.0024603}, 
					     {0.0073139,0.0091679,0.009217,0.0092056,0.0098387,0.0055790,0.0054113}, 
					     {0.0041850,0.0062754,0.007403,0.0059580,0.0048004,0.0030273,0.0038699}, 
					     {0.0081684,0.0102559,0.009914,0.0100281,0.0098609,0.0059393,0.0063520}, 
					     {0.0076389,0.0118545,0.014384,0.0121222,0.0075882,0.0054894,0.0063555}, 
					     {0.0044918,0.0059000,0.009171,0.0059351,0.0061848,0.0037100,0.0032804}, 
					     {0.0083015,0.0087119,0.009145,0.0097753,0.0086994,0.0052928,0.0050537}, 
					     {0.0028978,0.0040778,0.004154,0.0036926,0.0037117,0.0023811,0.0024050}, 
					     {0.0199521,0.0259806,0.031607,0.0244118,0.0236672,0.0167017,0.0166448} };

  //syst uncertainty from estimating true trigger eff from tag and probe
  const double trig_effmatrix_tpsyst[10][7] = { {0.00167,0.018063,0.012382,0.023539,0.028395,0.025499,0.026552}, 
						{0.03472,0.031546,0.027462,0.030010,0.030677,0.030553,0.028537}, 
						{0.02786,0.034864,0.027545,0.025706,0.030521,0.027495,0.029857}, 
						{0.02569,0.014712,0.015116,0.023172,0.023206,0.026210,0.026698}, 
						{0.17623,0.095552,0.083617,0.163108,0.158083,0.141751,0.166614}, 
						{0.18989,0.082965,0.106507,0.135236,0.109509,0.092677,0.172741}, 
						{0.02529,0.019854,0.017741,0.020040,0.025080,0.025643,0.025920}, 
						{0.02787,0.024241,0.031086,0.041334,0.028221,0.031177,0.034623}, 
						{0.03090,0.027656,0.028907,0.026535,0.029224,0.029647,0.029066}, 
						{0.03871,0.041957,0.045360,0.043003,0.031884,0.030473,0.042689} };
  
  for (int i=6; i>=0; i--){    // find phi index
    if ( mu_phi > phibins[i] ) {
      phiI = i;
      break;
    }
  }
  for (int i=9; i>=0; i--){    // find eta index
    if ( mu_eta > etabins[i] ) {
      etaI = i;
      break;
    }
  }

  //bias from T&P only overestimates eff, so only associated to -1 sigma variation
  if(isUp)
    return trig_effmatrix_err[etaI][phiI];
  else
    return sqrt( (trig_effmatrix_err[etaI][phiI] * trig_effmatrix_err[etaI][phiI]) + (trig_effmatrix_tpsyst[etaI][phiI] * trig_effmatrix_tpsyst[etaI][phiI]) );
  
}

// For reco+ID scale factors
double mu_recoID_SF()
{

	printWarning();   //this is not to be used anymore
    return -99.;

}

// For reco+ID scale factor total uncertainties. 
double mu_recoID_SF_err()
{

    printWarning();   //this is not to be used anymore
    return -99.;

}

// For reco+ID scale factor statistical uncertainties (symmetric)
double mu_recoID_SF_staterr()
{

    printWarning();   //this is not to be used anymore
	return -99.;    

}

// For reco+ID scale factor systematic uncertainties 
double mu_recoID_SF_systerr()
{

    printWarning();   //this is not to be used anymore
	return -99.;

}

//-----------------------------------------------------------------

// For ID scale factors
double mu_ID_SF()
{

	const double ID_SF = 1.0008;  // for EPS, rel16, just a number
    return ID_SF;

}

// For ID scale factor total uncertainties. 
double mu_ID_SF_err()
{

    double staterr, systerr, toterr;

    staterr = mu_ID_SF_staterr();
    systerr = mu_ID_SF_systerr();

    toterr = sqrt(pow(staterr,2) + pow(systerr,2));

    return toterr;

}

// For ID scale factor statistical uncertainties (symmetric)
double mu_ID_SF_staterr()
{

	const double SFerr = 0.0003;  //for EPS, in rel16, just a number
	return SFerr;    

}

// For ID scale factor systematic uncertainties 
double mu_ID_SF_systerr()
{

	const double SFerr = 0.0003;
	return SFerr;

}

// Helper function to deal with possible phi ambiguity
double check_phi_range(double phi)
{

	double newphi = phi;
	const double pi = acos(-1.);

	if (newphi > pi) {
		printf("<muon_SF>: WARNING: Muon phi %4.2f > pi! Using (phi-2*pi) \n", phi);
		newphi -= 2*pi;
	}
	if (newphi < -pi) {
		printf("<muon_SF>: WARNING: Muon phi %4.2f < -pi! Using (phi+2*pi) \n", phi);
		newphi += 2*pi;
	}

	return newphi;
}

void printWarning()
{

   printf("<muon_SF>: WARNING: Reco SFs no longer provided by this tool! Use MCP group tool instead.\n");
   printf("<muon_SF>: See https://twiki.cern.ch/twiki/bin/view/AtlasProtected/MCPAnalysisGuidelinesRel16\n");
   return;

}

#endif
