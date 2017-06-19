#ifndef TopElectronSFTool_h
#define TopElectronSFTool_h

/****************************************************************************
 * TopElectronSFTool.h                                                      *
 *                                                                          *
 * Simple functions which return data/MC scale factors                      *
 * given a electron's eta 												    *
 * Also functions for corresponding uncertainties.                          *
 *                                                                          *
 * double ele_ID_SF(double et, double ET)                   	            *
 * double ele_ID_SF_err(double eta, double ET)               	      	    *
 * double ele_reco_SF(double eta);                                          *
 * double ele_reco_SF_err(double eta);                                      *
 * double ele_recoID_SF(double eta, double ET);             	            *
 * double ele_recoID_SF_err(double eta, double ET);          	      	    *
 * double ele_trigger_SF(double eta, double ET, int set);                   *
 * double ele_trigger_SF_err(double eta, double ET, int set);          	    *
 * double ele_ID_SF_AFII(double eta, double ET);                            *
 * double ele_ID_SF_err_AFII(double eta, double ET);                        *
 * double ele_recoID_SF_AFII(double eta, double ET);                        *
 * double ele_recoID_SF_err_AFII(double eta, double ET);                    *
 * double ele_trigger_SF_AFII(double eta, double ET, int set);              *
 * double ele_trigger_SF_err_AFII(double eta, double ET, int set);          *
 *                                                                          *
 * set = 0, 1, 2                                                            *
 * correpond to trigger e20_medium, e22_medium, e22vh_medium1               *
 * History                                                                  *
 *         13 Jan 2012 -- created by H. Zhu                                 *   
 *         27 Jan 2012 -- updated by H. Zhu                                 *   
 *         10 Feb 2012 -- updated by H. Zhu                                 *
 *         20 Feb 2012 -- updated by H. Zhu                                 *
 ***************************************************************************/

#include <iostream>
#include <cmath>



class TopElectronSFTool {

  public:
  TopElectronSFTool(){}; 
  ~TopElectronSFTool(){};

  private:
	// declarations of functions
	double ele_ID_SF(double eta, double ET);
	double ele_ID_SF_err(double eta, double ET);
	double ele_reco_SF(double eta);
	double ele_reco_SF_err(double eta);
	double ele_recoID_SF(double eta, double ET);
	double ele_recoID_SF_err(double eta, double ET);
	double ele_trigger_SF(double eta, double ET, int set);
	double ele_trigger_SF_err(double eta, double ET, int set);
	
	double ele_ID_SF_AFII(double eta, double ET);
	double ele_ID_SF_err_AFII(double eta, double ET);
	double ele_recoID_SF_AFII(double eta, double ET);
	double ele_recoID_SF_err_AFII(double eta, double ET);
	double ele_trigger_SF_AFII(double eta, double ET, int set);
	double ele_trigger_SF_err_AFII(double eta, double ET, int set);
	
	static const double ele_ID_etamax;
	static const double ele_ID_etacrack[2];
	static const double ele_ID_EtMin;
	static const double ele_ID_etabins[18];
	static const double ele_ID_ETbins[8];
	static const double ele_ID_SFmatrix[8][18]; 
	static const double ele_ID_errmatrix[8][18];
    

	static const double ele_reco_etamax;
	static const double ele_reco_etacrack[2]; 
	static const double ele_reco_etabins[9]; 
	static const double ele_reco_SFmatrix[9];
	static const double ele_reco_errmatrix[9];


	static const double ele_trigger_etamax;
	static const double ele_trigger_etacrack[2];
	static const double ele_trigger_EtMin;
	static const double ele_trigger_etabins[18];
	static const double ele_trigger_ETbins[6];
    
	static const double ele_trigger_SFmatrix_e20_medium[6][18];
	
	static const double ele_trigger_SFmatrix_e22_medium[6][18];

	static const double ele_trigger_SFmatrix_e22vh_medium[6][18];

	static const double ele_trigger_errmatrix_e20_medium[6][18];
	
	static const double ele_trigger_errmatrix_e22_medium[6][18];
	
	static const double ele_trigger_errmatrix_e22vh_medium[6][18];

	
	static const double ele_ID_AFII_SFmatrix[8][18];

	static const double ele_ID_AFII_errmatrix[8][18];

	static const double ele_trigger_AFII_SFmatrix_e20_medium[6][18];
	
	static const double ele_trigger_AFII_SFmatrix_e22_medium[6][18];

	static const double ele_trigger_AFII_SFmatrix_e22vh_medium[6][18];

	static const double ele_trigger_AFII_errmatrix_e20_medium[6][18];

	static const double ele_trigger_AFII_errmatrix_e22_medium[6][18];

	static const double ele_trigger_AFII_errmatrix_e22vh_medium[6][18];
};



#endif
