#ifndef WJetSF_EPS_ttbar_h
#define WJetSF_EPS_ttbar_h

/*********************************************************************************************************************************
 * WJetSF_EPS_ttbar.h                        											 *
 *                                                         									 *
 * double WJet_SF(int HFOR) is a function which returns a data/MC flavor scale factor of a given event.          		 *
 *  The inputs of this function is the HFOR flag. The flavor fractions defined are bb, cc,and  c flavor jets that corresponds    *
 *  with the HFOR flag equals to 0, 1, and 2 respectively. If the HFOR flag is different than those vales the function returns 1.*
 *																 *
 *  This file also has defined several functions for the systematic uncertainties, like:					 *
 *    double WJet_SF_err_TtbarDown(int HFOR, int jetn)										 *
 *  The input of this function is also the HFOR flag. The output in this case is the flavor scale factor after the ttbar 
 *  cross section has been shifted down by 30%.				                                                 	 *
 *																 *
 *  The rest of the functions for the different systematic uncertainties are:							 *
 *    double WJet_SF_err_TtbarUp(int HFOR);										         *
 *    double WJet_SF_err_SgTopDown(int HFOR);										         *
 *    double WJet_SF_err_SgTopUp(int HFOR);										         *
 *    double WJet_SF_err_ZjetsDown(int HFOR);										         *
 *    double WJet_SF_err_ZjetsUp(int HFOR);										         *
 *    double WJet_SF_err_DibDown(int HFOR);										         *
 *    double WJet_SF_err_DibUp(int HFOR);										         *
 *    double WJet_SF_err_QCDDown(int HFOR);										         *
 *    double WJet_SF_err_QCDUp(int HFOR);										         *
 *    double WJet_SF_err_JESDown(int HFOR);										         *
 *    double WJet_SF_err_JESUp(int HFOR);										         *
 *    double WJet_SF_err_BtagSFDown(int HFOR);									                 *	         
 *    double WJet_SF_err_BtagSFUp(int HFOR);										         *
 *    double WJet_SF_err_LqtagSFDown(int HFOR);									                 *	         
 *    double WJet_SF_err_LqtagSFUp(int HFOR);										         *
 *																 *
 *  The WJet_SF_Staterr(int HFOR), WJet_SF_Systerr(int HFOR), and WJet_SF_Totalerr(int HFOR) 	                                 *
 *  functions are a little bit different since they return the value of the statistical and the systematic shift, 		 *
 *  they don't return the scale factor with those statistical, or the systematic uncert. applied as the previous functions.      *
 *  i.e., the scale factor will be expresed like :the nominal +/- its statistical uncertainty +/- its systematic uncertainties:	 *
 *  WJet_SF(int HFOR) +/- WJet_SF_Staterr(int HFOR) (stat) +/- WJet_SF_Systerr(int HFOR) (sys)	                                 *
 *																 *
 *  or showing only the total uncertainty:											 *
 *																 *
 *  WJet_SF(int HFOR) +/- WJet_SF_Totalerr(int HFOR)(sys)						                    	 *
 *																 *
 * History                                                    									 *
 *         14 Feb 2011 -- created by Jenny Holzbauer          									 *
 *         13 May 2011 -- updated by  Barbara Alvarez, including the pLHC files                                          	 *
 *         21 Jun 2011 -- updated by  Barbara Alvarez, with preliminary results for EPS                                      	 *
 *                                                            									 *
 *********************************************************************************************************************************/

#include <iostream>
#include <math.h>

// forward declaration of functions
double WJet_SF(int HFOR);
//systematic uncertainties
double WJet_SF_err_TtbarDown(int HFOR);
double WJet_SF_err_TtbarUp(int HFOR);
double WJet_SF_err_SgTopDown(int HFOR);
double WJet_SF_err_SgTopUp(int HFOR);
double WJet_SF_err_ZjetsDown(int HFOR);
double WJet_SF_err_ZjetsUp(int HFOR);
double WJet_SF_err_DibDown(int HFOR);
double WJet_SF_err_DibUp(int HFOR);
double WJet_SF_err_QCDDown(int HFOR);
double WJet_SF_err_QCDUp(int HFOR);
double WJet_SF_err_JESDown(int HFOR);
double WJet_SF_err_JESUp(int HFOR);
double WJet_SF_err_BtagSFDown(int HFOR);
double WJet_SF_err_BtagSFUp(int HFOR);
double WJet_SF_err_LqtagSFDown(int HFOR);
double WJet_SF_err_LqtagSFUp(int HFOR);

//statistical error
double WJet_SF_Staterr(int HFOR);

//total systematic uncertainty
double WJet_SF_Systerr(int HFOR);
//total uncertainty : stat+syst
double WJet_SF_Totalerr(int HFOR);

//fucntion that returns the flavor scale factor
double WJet_SF(int HFOR)
{

// HFOR is the heavy flavor flag: 0 for bb, 1 for cc, and 2 for c 
   
  //rel 16 numbers
  const double SFmatrix[2] = {1.6331, 1.1071}; //HOFR  bb/cc, and c in this order  

  //check HF flag 
  if(HFOR>2 || HFOR<0 ){
    return 1.0;       

  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor
    
    return SFmatrix[HFOR_index];
  }
}
//ttbar cross section shift down
double WJet_SF_err_TtbarDown(int HFOR)
{
  //30% shift down on the ttbar contribution
  
  const double SFmatrix_err[2] = {1.8676, 1.0455}; //HOFR  bb/cc, and c, in this order  
  
  //check HF flag 
  if(HFOR>2 || HFOR<0 ){
    return 1.0;       
    
  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor
    
    return SFmatrix_err[HFOR_index];
  }

}
//ttbar cross section shift up
double WJet_SF_err_TtbarUp(int HFOR)
{
  //30% shift up on the ttbar contribution
  
  const double SFmatrix_err[2] = {1.3969, 1.1693}; //HOFR  bb/cc, and c, in this order  

  
  
  //check HF flag 
  if(HFOR>2 || HFOR<0 ){
    return 1.0;       
    
  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor
    
    return SFmatrix_err[HFOR_index];
  }

}

//single top cross section shift down
double WJet_SF_err_SgTopDown(int HFOR)
{
  //14% shift on the single top cross section
  

  const double SFmatrix_err[2] = {1.7421, 1.0846}; //HOFR  bb/cc, and c, in this order  
  
  
  //check HF flag 
  if(HFOR>2 || HFOR<0 ){
    return 1.0;       
    
  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor
    
    return SFmatrix_err[HFOR_index];
  }

}
//single top cross section shift up
double WJet_SF_err_SgTopUp(int HFOR)
{
  
  const double SFmatrix_err[2] = {1.5236, 1.1298}; //HOFR  bb/cc, and c, in this order  

  
  
  //check HF flag 
  if(HFOR>2 || HFOR<0 ){
    return 1.0;       
    
  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor
    
    return SFmatrix_err[HFOR_index];
  }

}

//Z+jest cross section shift down 30% shift
double WJet_SF_err_ZjetsDown(int HFOR)
{
  
  const double SFmatrix_err[2] = {1.5735, 1.0926}; //HOFR  bb/cc, and c, in this order  
  
  //check HF flag 
  if(HFOR>2 || HFOR<0 ){
    return 1.0;       
    
  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor
    
    return SFmatrix_err[HFOR_index];
  }

}
//Z+jest cross section shift up
double WJet_SF_err_ZjetsUp(int HFOR)
{
  
  const double SFmatrix_err[2] = {1.6990, 1.1213}; //HOFR  bb/cc, and c, in this order  

  //check HF flag 
  if(HFOR>2 || HFOR<0 ){
    return 1.0;       
    
  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor
    
    return SFmatrix_err[HFOR_index];
  }

}

//Diboson cross section shift down 5% shiff
double WJet_SF_err_DibDown(int HFOR)
{
  
  const double SFmatrix_err[2] = {1.6342, 1.1072}; //HOFR  bb/cc, and c, in this order  
  
  
  //check HF flag 
  if(HFOR>2 || HFOR<0 ){
    return 1.0;       
    
  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor
    
    return SFmatrix_err[HFOR_index];
  }

}
//Diboson cross section shift up
double WJet_SF_err_DibUp(int HFOR)
{
  
  
  const double SFmatrix_err[2] = {1.6320, 1.1070}; //HOFR  bb/cc, and c, in this order  
  //check HF flag 
  if(HFOR>2 || HFOR<0 ){
    return 1.0;       
    
  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor
    
    return SFmatrix_err[HFOR_index];
  }

}

//QCD shift Down estimated by the data driven methods
double WJet_SF_err_QCDDown(int HFOR)
{
  
  const double SFmatrix_err[2] = {1.6167, 1.1688}; //HOFR  bb/cc, and c, in this order  
  //check HF flag 
  if(HFOR>2 || HFOR<0 ){
    return 1.0;       
    
  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor
    
    return SFmatrix_err[HFOR_index];
  }

}
//QCD shift Up estimated by the data driven methods 
double WJet_SF_err_QCDUp(int HFOR)
{
  
  const double SFmatrix_err[2] = {1.6479, 1.0443}; //HOFR  bb/cc, and c, in this order  
  
  //check HF flag 
  if(HFOR>2 || HFOR<0 ){
    return 1.0;       
    
  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor
    
    return SFmatrix_err[HFOR_index];
  }

}
//JES up, NEED to be UPDATED
double WJet_SF_err_JESUp(int HFOR)
{
  
  const double SFmatrix_err[2] = {2.2577, 1.2748}; //HOFR  bb/cc, and c, in this order  
  //check HF flag 
  if(HFOR>2 || HFOR<0 ){
    return 1.0;       
    
  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor
    
    return SFmatrix_err[HFOR_index];
  }

}
//JES down. NEED to be UPDATED
double WJet_SF_err_JESDown(int HFOR)
{
  
  const double SFmatrix_err[2] = {1.3043, 0.9816}; //HOFR  bb/cc, and c, in this order  
  
  //check HF flag 
  if(HFOR>2 || HFOR<0 ){
    return 1.0;       
    
  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor
    
    return SFmatrix_err[HFOR_index];
  }

}
//B/C tag scale factors shiftted up
double WJet_SF_err_BtagSFUp(int HFOR)
{
  
  const double SFmatrix_err[2] = {1.2683, 0.8911}; //HOFR  bb/cc, and c, in this order  
  
  //check HF flag 
  if(HFOR>2 || HFOR<0 ){
    return 1.0;       
    
  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor
    
    return SFmatrix_err[HFOR_index];
  }

}

//B/C tag scale factors shiftted down
 
double WJet_SF_err_BtagSFDown(int HFOR)
{
  
  const double SFmatrix_err[2] = {2.1477, 1.4730}; //HOFR  bb/cc, and c, in this order  
  //check HF flag 
  if(HFOR>2 || HFOR<0 ){
    return 1.0;       
    
  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor
    
    return SFmatrix_err[HFOR_index];
  }

}

//Lq tag scale factor shitd up 
double WJet_SF_err_LqtagSFUp(int HFOR)
{
  

  
  const double SFmatrix_err[2] = {1.4915, 1.0585}; //HOFR  bb/cc, and c, in this order  
  
  //check HF flag 
  if(HFOR>2 || HFOR<0 ){
    return 1.0;       
    
  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor
    
    return SFmatrix_err[HFOR_index];
  }

}

//Lq tag scale factor shitd down
double WJet_SF_err_LqtagSFDown(int HFOR)
{
  
  const double SFmatrix_err[2] = {1.7731, 1.1545}; //HOFR  bb/cc, and c, in this order  
  //check HF flag 
  if(HFOR>2 || HFOR<0 ){
    return 1.0;       
    
  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor
    
    return SFmatrix_err[HFOR_index];
  }

}

//Statistical ucnertainty (symmetric)
double WJet_SF_Staterr(int HFOR)
{
  const double SFmatrix_err[2] = {12.7754*1.6331, 9.0186*1.1071}; //HOFR  bb/cc, and c, in this order  

  //check HF flag 
  if(HFOR>2 || HFOR<0 ){
    return 0.0;       
  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor
    
    return SFmatrix_err[HFOR_index]/100;
  }

}
//Total systematic uncertanties up and down squared root (symmetric)
double WJet_SF_Systerr(int HFOR)
{
  
  
  const double SFmatrix_err[2] = {44.9723*1.6331, 31.6671*1.1071}; //HOFR  bb/cc, and c, in this order  
  
  //check HF flag 
  if(HFOR>2 || HFOR<0 ){
    return 0.0;       
  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor
    
    return SFmatrix_err[HFOR_index]/100;
  }

}
//Total uncertaintu stat + syst. (symmetric)
double WJet_SF_Totalerr(int HFOR)
{
  
  
  const double SFmatrix_err[2] = {46.7516*1.6331, 32.9263*1.1071}; //HOFR  bb/cc, and c, in this order  
  
  //check HF flag 
  if(HFOR>2 || HFOR<0 ){
    return 0.0;       
  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor
    
    return SFmatrix_err[HFOR_index]/100;
  }

}

#endif
