#ifndef electron_SF_pLHC_sgtop_h
#define electron_SF_pLHC_sgtop_h

/****************************************************************************
 * electron_SF_pLHC_sgtop.h                                                        *
 *                                                                          *
 * Simple functions which return data/MC scale factors                      *
 * given a electron's eta 													*
 * Also functions for corresponding uncertainties.                          *
 *                                                                          *
 * double ele_ID_SF(double et, double ET)                   	            *
 * double ele_ID_SF_err(double eta, double ET)               	      		*
 * double ele_reco_SF();                                    	            *
 * double ele_reco_SF_err();                                 	      		*
 * double ele_recoID_SF(double eta, double ET);             	            *
 * double ele_recoID_SF_err(double eta, double ET);          	      		*
 * double ele_trigger_SF();                                 	            *
 * double ele_trigger_SF_err();                              	      		*
 *                                                                          *
 * History                                                                  *
 *         26 Jan 2011 -- created by S. Caughron                            *   
 *         06 May 2011 -- updated for pLHC / single top                     *   
 *         14 May 2011 -- updated for pLHC / single top                     *   
 ***************************************************************************/

#include <iostream>
#include <cmath>

// forward declaration of functions
double ele_ID_SF(double eta, double ET);
double ele_ID_SF_err(double eta, double ET);
double ele_reco_SF();
double ele_reco_SF_err();
double ele_recoID_SF(double eta, double ET);
double ele_recoID_SF_err(double eta, double ET);
double ele_trigger_SF();
double ele_trigger_SF_err();


// For ID+iso scale factors
double ele_ID_SF(double eta, double ET)
{

    double ele_eta = eta;
    double ele_ET = ET;
    int etaI=-1; int ET_I=-1;   

    const double etabins[8] = {-2.47,-2.01,-1.37,-0.8,  // lower edges of eta bins
							   0.0,0.8,1.52,2.01};   
    const double ETbins[6] = {20.,25.,30.,35.,40.,45.}; // lower edges of eta bins

    const double SFmatrix[6][8] = {{ 0.888, 0.916, 0.940, 0.886, 0.898, 0.934, 0.950, 0.894 },
								   { 0.925, 0.955, 0.980, 0.924, 0.936, 0.974, 0.990, 0.932 },
								   { 0.976, 1.007, 1.033, 0.974, 0.987, 1.026, 1.044, 0.983 },
								   { 0.984, 1.016, 1.043, 0.983, 0.996, 1.036, 1.053, 0.992 },
								   { 0.980, 1.011, 1.038, 0.978, 0.991, 1.031, 1.048, 0.987 },
								   { 0.987, 1.019, 1.046, 0.986, 0.999, 1.039, 1.056, 0.995 }};

	const double etamax = 2.47;
	const double etacrack[2] = {1.37,1.52};

    if ( std::fabs(ele_eta) > etamax || (std::fabs(ele_eta) > etacrack[0] && std::fabs(ele_eta) < etacrack[1]) )  // check forward, crack regions

        return 1.0;       

    else {           
        
        for (int i=7; i>=0; i--){    // find eta index
            if ( ele_eta > etabins[i] ) {
                etaI = i;
                break;
            }
        }
        for (int i=5; i>=0; i--){    // find eta index
            if ( ele_ET > ETbins[i] ) {
                ET_I = i;
                break;
            }
        }

        return SFmatrix[ET_I][etaI];

    } //else

    return 0;
}

// For ID+iso scale factor uncertainties (symmetric)
double ele_ID_SF_err(double eta, double ET)
{

    double ele_eta = eta;
    double ele_ET = ET;
    int etaI=-1; int ET_I=-1;   

    const double etabins[8] = {-2.47,-2.01,-1.37,-0.8,  // lower edges of eta bins
							   0.0,0.8,1.52,2.01};   
    const double ETbins[6] = {20.,25.,30.,35.,40.,45.}; // lower edges of eta bins

    const double errmatrix[6][8] = {{ 0.047, 0.050, 0.049, 0.046, 0.045, 0.048, 0.049, 0.052 },
									{ 0.036, 0.040, 0.039, 0.035, 0.034, 0.037, 0.038, 0.043 },
									{ 0.035, 0.039, 0.038, 0.034, 0.032, 0.036, 0.037, 0.042 },
									{ 0.034, 0.039, 0.037, 0.033, 0.032, 0.036, 0.037, 0.042 },
									{ 0.035, 0.039, 0.038, 0.033, 0.032, 0.036, 0.037, 0.042 },
									{ 0.039, 0.043, 0.041, 0.037, 0.036, 0.040, 0.041, 0.045 }};

	const double etamax = 2.47;
	const double etacrack[2] = {1.37,1.52};

    if ( std::fabs(ele_eta) > etamax || (std::fabs(ele_eta) > etacrack[0] && std::fabs(ele_eta) < etacrack[1]) )  // check forward, crack regions

        return 1.0;       

    else {           
        
        for (int i=7; i>=0; i--){    // find eta index
            if ( ele_eta > etabins[i] ) {
                etaI = i;
                break;
            }
        }
        for (int i=5; i>=0; i--){    // find ET index
            if ( ele_ET > ETbins[i] ) {
                ET_I = i;
                break;
            }
        }
		
        return errmatrix[ET_I][etaI];

    } //else

    return 0;
}

double ele_reco_SF()
{

	const double SF = 1.013;  //for rel15 (and rel16), just a number
    return SF;

}

// For reco SF uncertainties (symmetric)
double ele_reco_SF_err()
{

	const double SFerr = 0.015;  //for rel15, just a number
    return SFerr;

}

//returns cumulative reco+ID SF
double ele_recoID_SF(double eta, double ET)
{

	double ele_eta = eta;
	double ele_ET = ET;
	double total_SF = ele_ID_SF(ele_eta, ele_ET) * ele_reco_SF();

	return total_SF;

}

//returns cumulative reco+ID SF uncertainty
double ele_recoID_SF_err(double eta, double ET)
{

	double ele_eta = eta;
	double ele_ET = ET;
	double ID_err = ele_ID_SF_err(ele_eta, ele_ET) / ele_ID_SF(ele_eta, ele_ET);  //need relative errors
	double reco_err = ele_reco_SF_err() / ele_reco_SF();

	double tot_rel_err = sqrt( pow(ID_err,2) + pow(reco_err,2) );
	double tot_abs_err = tot_rel_err * ele_recoID_SF(ele_eta, ele_ET);

	return tot_abs_err;

}

// For trigger SFs
double ele_trigger_SF()
{

	const double SF = 0.995;  //for rel15 (and rel16), just a number
    return SF;

}

// For trigger SF uncertainties (symmetric)
double ele_trigger_SF_err()
{

	const double SFerr = 0.010;  //for rel15 (and rel16) , just a number
    return SFerr;

}

#endif
