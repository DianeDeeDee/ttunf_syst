#ifndef electron_SF_R15_h
#define electron_SF_R15_h

/****************************************************************************
 * electron_SF_R15.h                                                        *
 *                                                                          *
 * Simple functions which return data/MC scale factors                      *
 * given a electron's eta 													*
 * Also functions for corresponding uncertainties.                          *
 *                     .                                                    *
 * double ele_ID_SF(double eta)	                            	            *
 * double ele_ID_SF_err(double eta)                          	      		*
 * double ele_reco_SF(double eta);                          	            *
 * double ele_reco_SF_err(double eta);                       	      		*
 * double ele_recoID_SF(double eta);                        	            *
 * double ele_recoID_SF_err(double eta);                     	      		*
 *                                                                          *
 * History                                                                  *
 *         26 Jan 2011 -- created by S. Caughron                            *   
 ***************************************************************************/

#include <iostream>
#include <cmath>

// forward declaration of functions
double ele_ID_SF(double eta);
double ele_ID_SF_err(double eta);
double ele_reco_SF(double eta);
double ele_reco_SF_err(double eta);
double ele_recoID_SF(double eta);
double ele_recoID_SF_err(double eta);


// For ID+iso scale factors
double ele_ID_SF(double eta)
{

    double ele_eta = eta;
    int etaI=-1;   

    const double etabins[8] = {-2.47,-2.01,-1.37,-0.8,  // lower edges of eta bins
							   0.0,0.8,1.52,2.01};   

    const double SFvector[8] = {0.921,0.985,0.993,0.978,   
                                0.984,1.024,1.022,0.934};

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
		
        return SFvector[etaI];

    } //else

    return 0;
}

// For ID+iso uncertainties (symmetric)
double ele_ID_SF_err(double eta)
{

    double ele_eta = eta;
    int etaI=-1;   

    const double etabins[8] = {-2.47,-2.01,-1.37,-0.8,  // lower edges of eta bins
							   0.0,0.8,1.52,2.01};   

    const double errvector[8] = {0.039,0.031,0.029,0.027,   
                                0.027,0.032,0.049,0.039};

	const double etamax = 2.47;
	const double etacrack[2] = {1.37,1.52};

    if ( std::fabs(ele_eta) > etamax || (std::fabs(ele_eta) > etacrack[0] && std::fabs(ele_eta) < etacrack[1]) )  // check forward, crack regions

        return 0.0;       

    else {           
        
        for (int i=7; i>=0; i--){    // find eta index
            if ( ele_eta > etabins[i] ) {
                etaI = i;
                break;
            }
        }
		
        return errvector[etaI];

    } //else

    return 0;
}

//reco SFs (flat and equal to 1 in rel15)
double ele_reco_SF(double eta)
{

// Will be binned same as ID in rel16
/*
    double ele_eta = eta;
    int etaI=-1;   

    const double etabins[8] = {-2.47,-2.01,-1.37,-0.8,  // lower edges of eta bins
							   0.0,0.8,1.52,2.01};   

    const double SFvector[8] = {0.921,0.985,0.993,0.978,   
                                0.984,1.024,1.022,0.934};

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
		
        return SFvector[etaI];

    } //else
*/

	const double SF = 1.0;  //for rel15, just a number
    return SF;
}

// For reco SF uncertainties (symmetric)
double ele_reco_SF_err(double eta)
{

// Will be binned same as ID in rel16
/*
    double ele_eta = eta;
    int etaI=-1;   

    const double etabins[8] = {-2.47,-2.01,-1.37,-0.8,  // lower edges of eta bins
							   0.0,0.8,1.52,2.01};   

    const double errvector[8] = {0.039,0.031,0.029,0.027,   
                                0.027,0.032,0.049,0.039};

	const double etamax = 2.47;
	const double etacrack[2] = {1.37,1.52};

    if ( std::fabs(ele_eta) > etamax || (std::fabs(ele_eta) > etacrack[0] && std::fabs(ele_eta) < etacrack[1]) )  // check forward, crack regions

        return 0.0;       

    else {           
        
        for (int i=7; i>=0; i--){    // find eta index
            if ( ele_eta > etabins[i] ) {
                etaI = i;
                break;
            }
        }
		
        return errvector[etaI];

    } //else
*/

	const double SFerr = 0.015;  //for rel15, just a number
    return SFerr;
}

//returns cumulative reco+ID SF
double ele_recoID_SF(double eta)
{

	double ele_eta = eta;
	double total_SF = ele_ID_SF(ele_eta) * ele_reco_SF(ele_eta);

	return total_SF;

}

//returns cumulative reco+ID SF uncertainty
double ele_recoID_SF_err(double eta)
{

	double ele_eta = eta;
	double ID_err = ele_ID_SF_err(ele_eta) / ele_ID_SF(ele_eta);  //need relative errors
	double reco_err = ele_reco_SF_err(ele_eta) / ele_reco_SF(ele_eta);

	double tot_rel_err = sqrt( pow(ID_err,2) + pow(reco_err,2) );
	double tot_abs_err = tot_rel_err * ele_recoID_SF(ele_eta);

	return tot_abs_err;

}
