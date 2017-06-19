#ifndef electron_SF_R16_h
#define electron_SF_R16_h

/****************************************************************************
 * electron_SF_R16.h                                                        *
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
 ***************************************************************************/

#include <iostream>
#include <cmath>
#include "Rtypes.h"

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

	const double SFmatrix[6][8] = {{0.917,0.946,0.968,0.907,0.912,0.970,0.961,0.953},   //rows are ET
  							 	   {0.960,0.990,1.013,0.949,0.955,1.016,1.006,0.998},	  //columns are eta
  							 	   {0.998,1.029,1.053,0.987,0.993,1.056,1.046,1.037},
  							 	   {0.996,1.027,1.051,0.985,0.991,1.054,1.044,1.035},
  							 	   {0.998,1.029,1.053,0.987,0.993,1.056,1.046,1.037},
  							 	   {1.007,1.038,1.062,0.995,1.002,1.065,1.055,1.046}};


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

	const double errmatrix[6][8] = {{ 0.082,  0.083,  0.082,  0.081,  0.081,  0.082,  0.086,  0.085 },     //rows are ET
									{ 0.028,  0.032,  0.030,  0.027,  0.027,  0.029,  0.038,  0.037  },    //columns are eta
									{ 0.027,  0.032,  0.029,  0.026,  0.026,  0.028,  0.038,  0.036  },
									{ 0.025,  0.030,  0.027,  0.023,  0.023,  0.025,  0.036,  0.034  },
									{ 0.025,  0.030,  0.028,  0.024,  0.024,  0.026,  0.037,  0.035  },
									{ 0.034,  0.038,  0.035,  0.033,  0.033,  0.034,  0.043,  0.041  }};

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

	const double SF = 1.0;  //for rel15 (and rel16), just a number
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

	const double SFerr = 0.005;  //for rel15 (and rel16) , just a number
    return SFerr;

}

#endif
