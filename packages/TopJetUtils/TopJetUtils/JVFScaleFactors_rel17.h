#ifndef JVFScaleFactors_rel17_h
#define JVFScaleFactors_rel17_h 1

/************************************************************************
 * JVFScaleFactors_rel17.h                                              *
 *                                                                      *
 * Simple functions which return data/MC scale factors given:		*
 * pT 		= jet pT (in GeV!!!!)					*
 * JVF 		= jet JVF value						*
 * matched 	= 1 when the jet was matched to a truth jet within 	*
 * 		DeltaR=0.4. Otherwise matched =0. Truth jets have been	*
 *		defined with a pT>5 GeV and eta<5.			*
 * year		= "2011", "2012" or "2012mc12a" (default) data		*
 * ps 		= parton showering used : 0 = Pythia, 1 = Herwig, 2 = Other *
 * 	    								*
 * double JVF_SF(double pT, double JVF, bool matched, std::string year, std::string ps) *
 * 									*
 * Also functions for corresponding uncertainties                       *
 * (not yet defined for 2012).                                          *
 *                                                                      *
 *									*
 * Example of use: Getting the nominal SFs				*
 *									*
 * #include "TopJetUtils/JVFScaleFactors_rel17.h"			*
 *									*
 * JVFweight = 1.;							*
 *									*
 * for (unsigned int i=0; i<jets.size(); ++i) {				*	
 *     									*
 *	double pT = jets[i]->Pt(); //In GeV!				*	
 *	double JVF = jets[i]->getMoment("JVF");				*
 *	bool matched = jets[i]->getMoment("matched");			*
 *									*
 *	JVFweight *= topjetutils::JVF_SF(pT,JVF,matched,"2012mc12a");	*
 * }									*
 * 									*
 * eventWeight *= JVFweight;						*
 *  
 * Example of use: To get the systematic uncertainty up/down
 * 
 * #include "TopJetUtils/JVFScaleFactors_rel17.h"
 * 
 * JVFweight = 1.;
 * 
 * Double_t SF_SasS = 1.; //SF for hard-scatter jet selection efficiency
 * Double_t SF_SasB = 1.; //SF for hard-scatter jet selection inefficiency
 * Double_t SF_BasB = 1.; //SF for pile-up jet rejection efficiency
 * Double_t SF_BasS = 1.; //SF for pile-up jet rejection inefficiency
 * Double_t sigmaSF_SasS = 1.; //SF variation up(down) for hard-scatter jet selection efficiency
 * Double_t sigmaSF_SasB = 1.; //SF variation up(down) for hard-scatter jet selection inefficiency
 * Double_t sigmaSF_BasB = 1.; //SF variation up(down) for pile-up jet rejection efficiency
 * Double_t sigmaSF_BasS = 1.; //SF variation up(down) for pile-up jet rejection inefficiency
 * 
 * int variation = 1 if doing up variation, -1 if doing down variation
 * 
 * for (unsigned int i=0; i<jets.size(); ++i) {	
 * 
 * 	double pT = jets[i]->Pt(); //In GeV!				*	
 * 	double JVF = jets[i]->getMoment("JVF");				*
 * 	bool matched = jets[i]->getMoment("matched");			*
 * 									*
 * 	JVFweight *= topjetutils::JVF_SF(pT,JVF,matched, "2011");
 * 
 * 	if (doJVFSFSystematic){
 *		if (fabs(JVF)>=0.75){
 *			if (matched){
 *				SF_SasS *= topjetutils::JVFSasS_SF(pT,"2011");
 *				sigmaSF_SasS *= topjetutils::JVFSasS_SF_err(pT,variation,"2011");
 *			} else {
 *				SF_BasS *= topjetutils::JVFBasS_SF(pT,"2011");
 *				sigmaSF_BasS *= topjetutils::JVFBasS_SF_err(pT, variation,"2011");
 *				}		
 *		} else {
 *			if (matched){
 *				SF_SasB *= topjetutils::JVFSasB_SF(pT,"2011");
 *				sigmaSF_SasB *= topjetutils::JVFSasB_SF_err(pT, variation,"2011");
 *			} else {
 *				SF_BasB *= topjetutils::JVFBasB_SF(pT,"2011");
 *				sigmaSF_BasB *= topjetutils::JVFBasB_SF_err(pT, variation,"2011");
 *			}	
 *		}
 *	}//Out of doJVFSFSystematic condition
 * 
 * } //Out of loop over jets
 * 
 * Double_t sigmaSFHS_2 = pow( (sigmaSF_SasS*sigmaSF_SasB - SF_SasS*SF_SasB)/(SF_SasS*SF_SasB) ,2);
 * Double_t sigmaSFPU_2 = pow( (sigmaSF_BasB*sigmaSF_BasS - SF_BasB*SF_BasS)/(SF_BasB*SF_BasS) ,2);
 *	
 * JVFweight *= (1 + variation*sqrt( sigmaSFHS_2 + sigmaSFPU_2 ) );
 
 * For a more detailed information please see: https://twiki.cern.ch/twiki/bin/view/AtlasProtected/TopCommonScales#JVF  
 * 									*
 * History                                                              *
 *         13 Fev 2012 -- created by R. Camacho/ S. Calvet              *    
 *         04 Aug 2012 -- 2012 Preliminary SFs by R. Camacho/ S. Calvet *    
 ************************************************************************/

//Headers

#include <string>

//Declaration of functions

namespace topjetutils{

double JVF_SF(const double pT, const double JVF, const bool matched, const std::string year, const int ps);

double JVFSasS_PythiaSF(const double jet_pT); //Returns the SasS SF for Pythia parton showering (year==2012mc12a)
double JVFSasS_HerwigSF(const double jet_pT); //Returns the SasS SF for Herwig parton showering (year==2012mc12a)
double JVFSasB_PythiaSF(const double jet_pT); //Returns the SasB SF for Pythia parton showering (year==2012mc12a)
double JVFSasB_HerwigSF(const double jet_pT); //Returns the SasB SF for Herwig parton showering (year==2012mc12a)

double JVFSasS_SF( const double pT, const std::string year, const int ps ); //SFs hard-scatter selected as hard-scatter jet 
double JVFSasB_SF( const double pT, const std::string year, const int ps ); //SFs hard-scatter selected as pile-up jet
double JVFBasS_SF( const double pT, const std::string year, const int ps ); //SFs pile-up selected as hard-scatter jet
double JVFBasB_SF( const double pT, const std::string year, const int ps ); //SFs pile-up selected as pile-up jet 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//General tools
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double exp_er(double x, double * par, double *matrix); 
double explinear_er(double x, double * par, double *matrix); 
double erfPlus_er(double x, double * par, double *matrix);
double erfMinus_er(double x, double * par, double *matrix);
double exp_func(double x, double *par);
double explinear_func(double x, double *par);
double erfPlus_func(double x, double *par);
double erfMinus_func(double x, double *par);

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Selection systematic uncertainty functions
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double JVFSasS_SF_errSelection(double pT, int up, std::string year);
double JVFSasB_SF_errSelection(double pT, int up, std::string year);
double JVFBasB_SF_errSelection(double pT, int up, std::string year);
double JVFBasS_SF_errSelection(double pT, int up, std::string year);

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Fit systematic uncertainty functions
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double JVFSasS_SF_errFit(double pT, std::string year);
double JVFSasB_SF_errFit(double pT, std::string year);
double JVFBasB_SF_errFit(double pT, std::string year);
double JVFBasS_SF_errFit(double pT, std::string year);

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Showering dependant systematics
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double JVFSasS_Averaged_Sys_SF( const int pT, const std::string year, const int up);
double JVFSasB_Averaged_Sys_SF( const int pT, const std::string year, const int up);
double JVFSasS_Herwig_Sys_SF( const int pT, const std::string year, const int up);
double JVFSasB_Herwig_Sys_SF( const int pT, const std::string year, const int up);
double JVFSasS_Pythia_Sys_SF( const int pT, const std::string year, const int up);
double JVFSasB_Pythia_Sys_SF( const int pT, const std::string year, const int up);


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Total systematic uncertainty functions
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double JVFSasS_SF_err(double pT, int up, std::string year, const int ps);
double JVFSasB_SF_err(double pT, int up, std::string year, const int ps);
double JVFBasB_SF_err(double pT, int up, std::string year, const int ps);
double JVFBasS_SF_err(double pT, int up, std::string year, const int ps);

}

#endif
