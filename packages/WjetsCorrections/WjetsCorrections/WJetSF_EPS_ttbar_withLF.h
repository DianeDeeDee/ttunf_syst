#ifndef WJetSF_EPS_ttbar_withLF_h
#define WJetSF_EPS_ttbar_withLF_h

/***********************************************************************************************************************
 * WJetSF_EPS_ttbar_withLF.h                        							       	       *
 *                                                         							       *
 * double WJet_SF(int HFOR) is a function which returns a data/MC flavor scale factor of a given event.                *
 *  The inputs of this function is the HFOR flag.                                                                      *
 *  The flavor fractions defined are bb, cc, c and light flavor jets that correspond                                   *
 *  to the HFOR flag equals to 0, 1, 2 and 3 respectively.                                                             *
 *  If the HFOR flag is different than those vales the function returns 1.                                             *
 *	                                                                                                               *
 * double WJet_SFLF(double Wbb_SF, double Wcc_SF, double Wc_SF) is a function returning the light flavor scale factor  *
 *  needed to keep the pretag normalisation unchanged while applying the scale factors given in input.                 *
 *  The flavor fractions before applying the other scale factors are defined in FractionMatrix[4] and can be           *
 *  adjusted to the needs of different analyses.                                                                       *
 *                                                                                                                     *
 *  This file also has defined several functions for the systematic uncertainties, like:			       *
 *    double WJet_SF_err_TtbarDown(int HFOR)	       								       *
 *  The input of this function is also the HFOR flag.                                                                  *
 *  The output in this case is the flavor scale factor after the ttbar cross section has been shifted down by 30%.     *
 *														       *
 *  The rest of the functions for the different systematic uncertainties are:					       *
 *    double WJet_SF_err_TtbarUp(int HFOR);									       *
 *    double WJet_SF_err_SgTopDown(int HFOR);									       *
 *    double WJet_SF_err_SgTopUp(int HFOR);									       *
 *    double WJet_SF_err_ZjetsDown(int HFOR);									       *
 *    double WJet_SF_err_ZjetsUp(int HFOR);									       *
 *    double WJet_SF_err_DibDown(int HFOR);									       *
 *    double WJet_SF_err_DibUp(int HFOR);									       *
 *    double WJet_SF_err_QCDDown(int HFOR);									       *
 *    double WJet_SF_err_QCDUp(int HFOR);									       *
 *    double WJet_SF_err_JESDown(int HFOR);									       *
 *    double WJet_SF_err_JESUp(int HFOR);									       *
 *    double WJet_SF_err_BtagSFDown(int HFOR);									       *
 *    double WJet_SF_err_BtagSFUp(int HFOR);									       *
 *    double WJet_SF_err_LqtagSFDown(int HFOR);									       *
 *    double WJet_SF_err_LqtagSFUp(int HFOR);									       *
 *														       *
 *  The WJet_SF_Staterr(int HFOR), WJet_SF_Systerr(int HFOR), and WJet_SF_Totalerr(int HFOR) 	                       *
 *  functions are a little bit different since they return the value of the statistical and the systematic shift,      *
 *  they don't return the scale factor with those statistical, or the systematic uncertainty applied.                  *
 *  i.e., the scale factor will be expresed like nominal +/- stat +/- syst 	                                       *
 *  WJet_SF(int HFOR) +/- WJet_SF_Staterr(int HFOR) (stat) +/- WJet_SF_Systerr(int HFOR) (sys)	                       *
 *  or showing only the total uncertainty:									       *
 *  WJet_SF(int HFOR) +/- WJet_SF_Totalerr(int HFOR)(sys)						               *
 *                                                                                                                     *
 *  The variations due to the statistical uncertainty alone, the extrapolation uncertainty alone (25%)                 *
 *  and the combination of the two for each flavor (they are uncorrelated) are given by the functions:                 *
 *    double WJet_SF_err_StatBBCCUp(int HFOR);                                                                         *
 *    double WJet_SF_err_StatBBCCDown(int HFOR);                                                                       *
 *    double WJet_SF_err_StatCUp(int HFOR);                                                                            *
 *    double WJet_SF_err_StatCDown(int HFOR);                                                                          *
 *    double WJet_SF_err_ExtrBBCCUp(int HFOR);                                                                         *
 *    double WJet_SF_err_ExtrBBCCDown(int HFOR);                                                                       *
 *    double WJet_SF_err_ExtrCUp(int HFOR);                                                                            *
 *    double WJet_SF_err_ExtrCDown(int HFOR);                                                                          * 
 *    double WJet_SF_err_StatExtrBBCCUp(int HFOR);                                                                     *
 *    double WJet_SF_err_StatExtrBBCCDown(int HFOR);                                                                   *
 *    double WJet_SF_err_StatExtrCUp(int HFOR);                                                                        *
 *    double WJet_SF_err_StatExtrCDown(int HFOR);                                                                      *
 *                                                                                                                     *
 *  The combination of the uncorrelated uncertainties on BtagSF and LqtagSF (sum in quadrature) is given in:           *
 *    double WJet_SF_err_BLqtagSFDown(int HFOR);                                                                       *
 *    double WJet_SF_err_BLqtagSFUp(int HFOR);                                                                         *
 *														       *
 * History                                                    							       *
 *         14 Feb 2011 -- created by Jenny Holzbauer          							       *
 *         13 May 2011 -- updated by Barbara Alvarez, including the pLHC files                                         *
 *         21 Jun 2011 -- updated by Barbara Alvarez, with preliminary results for EPS                                 *
 *         12 Aug 2011 -- updated by Lucia Masetti, including light flavor scale factors and more functions            *
 *                                                            							       *
 ***********************************************************************************************************************/

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

//new functions to help analysers

//LF scaling factor for a set of HF SFs
double WJet_SFLF(double Wbb_SF, double Wcc_SF, double Wc_SF);
//statistical variation, uncorrelated
double WJet_SF_err_StatBBCCUp(int HFOR);
double WJet_SF_err_StatBBCCDown(int HFOR);
double WJet_SF_err_StatCUp(int HFOR);
double WJet_SF_err_StatCDown(int HFOR);
//extrapolation uncertanty, uncorrelated
double WJet_SF_err_ExtrBBCCUp(int HFOR);
double WJet_SF_err_ExtrBBCCDown(int HFOR);
double WJet_SF_err_ExtrCUp(int HFOR);
double WJet_SF_err_ExtrCDown(int HFOR);
//statistical and extrapolation variation added in quadrature
double WJet_SF_err_StatExtrBBCCUp(int HFOR);
double WJet_SF_err_StatExtrBBCCDown(int HFOR);
double WJet_SF_err_StatExtrCUp(int HFOR);
double WJet_SF_err_StatExtrCDown(int HFOR);
//combined uncertainty on Btag and Lqtag efficiency
double WJet_SF_err_BLqtagSFDown(int HFOR);
double WJet_SF_err_BLqtagSFUp(int HFOR);

//LF scaling factor for a set of HF SFs
double WJet_SFLF(double Wbb_SF, double Wcc_SF, double Wc_SF)
{

  //Flavor factions in the considered pretag jet multiplicity bin without scaling in the HFOR order: bb, cc, c, light
  const double FractionMatrix[4] = {0.0621, 0.1102, 0.1124, 0.7154};
  double SFmatrix[4];
  SFmatrix[0] = Wbb_SF;
  SFmatrix[1] = Wcc_SF;
  SFmatrix[2] = Wc_SF;

  double CorrLFfrac = 1.;
  for (int i=0; i<3; i++){
    CorrLFfrac -= FractionMatrix[i]*SFmatrix[i];
  }
  SFmatrix[3] = CorrLFfrac/FractionMatrix[3];

  return SFmatrix[3];

}

//fucntion that returns the flavor scale factor
double WJet_SF(int HFOR)
{

// HFOR is the heavy flavor flag: 0 for bb, 1 for cc, and 2 for c 
   
  //rel 16 numbers
  const double SFmatrix[2] = {1.6331, 1.1071}; //HOFR  bb/cc, and c in this order  

  //check HF flag 
  if(HFOR>3 || HFOR<0 ){
    return 1.0;       

  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor

    if(HFOR_index>1){
      return WJet_SFLF(SFmatrix[0], SFmatrix[0], SFmatrix[1]);
    }
    else{
      return SFmatrix[HFOR_index];
    }
  }
}
//ttbar cross section shift down
double WJet_SF_err_TtbarDown(int HFOR)
{
  //30% shift down on the ttbar contribution
  
  const double SFmatrix_err[2] = {1.8676, 1.0455}; //HOFR  bb/cc, and c, in this order  
  
  //check HF flag 
  if(HFOR>3 || HFOR<0 ){
    return 1.0;       

  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor

    if(HFOR_index>1){
      return WJet_SFLF(SFmatrix_err[0], SFmatrix_err[0], SFmatrix_err[1]);
    }
    else{
      return SFmatrix_err[HFOR_index];
    }
  }
}
//ttbar cross section shift up
double WJet_SF_err_TtbarUp(int HFOR)
{
  //30% shift up on the ttbar contribution
  
  const double SFmatrix_err[2] = {1.3969, 1.1693}; //HOFR  bb/cc, and c, in this order  

  //check HF flag 
  if(HFOR>3 || HFOR<0 ){
    return 1.0;       

  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor

    if(HFOR_index>1){
      return WJet_SFLF(SFmatrix_err[0], SFmatrix_err[0], SFmatrix_err[1]);
    }
    else{
      return SFmatrix_err[HFOR_index];
    }
  }
}
//single top cross section shift down
double WJet_SF_err_SgTopDown(int HFOR)
{
  //14% shift on the single top cross section
  
  const double SFmatrix_err[2] = {1.7421, 1.0846}; //HOFR  bb/cc, and c, in this order  
  
  //check HF flag 
  if(HFOR>3 || HFOR<0 ){
    return 1.0;       

  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor

    if(HFOR_index>1){
      return WJet_SFLF(SFmatrix_err[0], SFmatrix_err[0], SFmatrix_err[1]);
    }
    else{
      return SFmatrix_err[HFOR_index];
    }
  }  
}
//single top cross section shift up
double WJet_SF_err_SgTopUp(int HFOR)
{
  
  const double SFmatrix_err[2] = {1.5236, 1.1298}; //HOFR  bb/cc, and c, in this order  

  //check HF flag 
  if(HFOR>3 || HFOR<0 ){
    return 1.0;       

  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor

    if(HFOR_index>1){
      return WJet_SFLF(SFmatrix_err[0], SFmatrix_err[0], SFmatrix_err[1]);
    }
    else{
      return SFmatrix_err[HFOR_index];
    }
  }  
}
//Z+jest cross section shift down 30% shift
double WJet_SF_err_ZjetsDown(int HFOR)
{
  
  const double SFmatrix_err[2] = {1.5735, 1.0926}; //HOFR  bb/cc, and c, in this order  
  
  //check HF flag 
  if(HFOR>3 || HFOR<0 ){
    return 1.0;       

  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor

    if(HFOR_index>1){
      return WJet_SFLF(SFmatrix_err[0], SFmatrix_err[0], SFmatrix_err[1]);
    }
    else{
      return SFmatrix_err[HFOR_index];
    }
  }  
}
//Z+jest cross section shift up
double WJet_SF_err_ZjetsUp(int HFOR)
{
  
  const double SFmatrix_err[2] = {1.6990, 1.1213}; //HOFR  bb/cc, and c, in this order  

  //check HF flag 
  if(HFOR>3 || HFOR<0 ){
    return 1.0;       

  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor

    if(HFOR_index>1){
      return WJet_SFLF(SFmatrix_err[0], SFmatrix_err[0], SFmatrix_err[1]);
    }
    else{
      return SFmatrix_err[HFOR_index];
    }
  }  
}

//Diboson cross section shift down 5% shiff
double WJet_SF_err_DibDown(int HFOR)
{
  
  const double SFmatrix_err[2] = {1.6342, 1.1072}; //HOFR  bb/cc, and c, in this order  
  
  //check HF flag 
  if(HFOR>3 || HFOR<0 ){
    return 1.0;       

  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor

    if(HFOR_index>1){
      return WJet_SFLF(SFmatrix_err[0], SFmatrix_err[0], SFmatrix_err[1]);
    }
    else{
      return SFmatrix_err[HFOR_index];
    }
  }  
}
//Diboson cross section shift up
double WJet_SF_err_DibUp(int HFOR)
{
  
  const double SFmatrix_err[2] = {1.6320, 1.1070}; //HOFR  bb/cc, and c, in this order  

  //check HF flag 
  if(HFOR>3 || HFOR<0 ){
    return 1.0;       

  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor

    if(HFOR_index>1){
      return WJet_SFLF(SFmatrix_err[0], SFmatrix_err[0], SFmatrix_err[1]);
    }
    else{
      return SFmatrix_err[HFOR_index];
    }
  }  
}
//QCD shift Down estimated by the data driven methods
double WJet_SF_err_QCDDown(int HFOR)
{
  
  const double SFmatrix_err[2] = {1.6167, 1.1688}; //HOFR  bb/cc, and c, in this order  

  //check HF flag 
  if(HFOR>3 || HFOR<0 ){
    return 1.0;       

  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor

    if(HFOR_index>1){
      return WJet_SFLF(SFmatrix_err[0], SFmatrix_err[0], SFmatrix_err[1]);
    }
    else{
      return SFmatrix_err[HFOR_index];
    }
  } 
}
//QCD shift Up estimated by the data driven methods 
double WJet_SF_err_QCDUp(int HFOR)
{
  
  const double SFmatrix_err[2] = {1.6479, 1.0443}; //HOFR  bb/cc, and c, in this order  
  
  //check HF flag 
  if(HFOR>3 || HFOR<0 ){
    return 1.0;       

  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor

    if(HFOR_index>1){
      return WJet_SFLF(SFmatrix_err[0], SFmatrix_err[0], SFmatrix_err[1]);
    }
    else{
      return SFmatrix_err[HFOR_index];
    }
  } 
}
//JES up, NEED to be UPDATED
double WJet_SF_err_JESUp(int HFOR)
{
  
  const double SFmatrix_err[2] = {2.2577, 1.2748}; //HOFR  bb/cc, and c, in this order  

  //check HF flag 
  if(HFOR>3 || HFOR<0 ){
    return 1.0;       

  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor

    if(HFOR_index>1){
      return WJet_SFLF(SFmatrix_err[0], SFmatrix_err[0], SFmatrix_err[1]);
    }
    else{
      return SFmatrix_err[HFOR_index];
    }
  } 
}
//JES down. NEED to be UPDATED
double WJet_SF_err_JESDown(int HFOR)
{
  
  const double SFmatrix_err[2] = {1.3043, 0.9816}; //HOFR  bb/cc, and c, in this order  
  
  //check HF flag 
  if(HFOR>3 || HFOR<0 ){
    return 1.0;       

  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor

    if(HFOR_index>1){
      return WJet_SFLF(SFmatrix_err[0], SFmatrix_err[0], SFmatrix_err[1]);
    }
    else{
      return SFmatrix_err[HFOR_index];
    }
  } 
}
//B/C tag scale factors shiftted up
double WJet_SF_err_BtagSFUp(int HFOR)
{
  
  const double SFmatrix_err[2] = {1.2683, 0.8911}; //HOFR  bb/cc, and c, in this order  
  
  //check HF flag 
  if(HFOR>3 || HFOR<0 ){
    return 1.0;       

  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor

    if(HFOR_index>1){
      return WJet_SFLF(SFmatrix_err[0], SFmatrix_err[0], SFmatrix_err[1]);
    }
    else{
      return SFmatrix_err[HFOR_index];
    }
  } 
}
//B/C tag scale factors shiftted down 
double WJet_SF_err_BtagSFDown(int HFOR)
{
  
  const double SFmatrix_err[2] = {2.1477, 1.4730}; //HOFR  bb/cc, and c, in this order  

  //check HF flag 
  if(HFOR>3 || HFOR<0 ){
    return 1.0;       

  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor

    if(HFOR_index>1){
      return WJet_SFLF(SFmatrix_err[0], SFmatrix_err[0], SFmatrix_err[1]);
    }
    else{
      return SFmatrix_err[HFOR_index];
    }
  } 
}
//Lq tag scale factor shitd up 
double WJet_SF_err_LqtagSFUp(int HFOR)
{
  
  const double SFmatrix_err[2] = {1.4915, 1.0585}; //HOFR  bb/cc, and c, in this order  
  
  //check HF flag 
  if(HFOR>3 || HFOR<0 ){
    return 1.0;       

  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor

    if(HFOR_index>1){
      return WJet_SFLF(SFmatrix_err[0], SFmatrix_err[0], SFmatrix_err[1]);
    }
    else{
      return SFmatrix_err[HFOR_index];
    }
  } 
}
//Lq tag scale factor shitd down
double WJet_SF_err_LqtagSFDown(int HFOR)
{
  
  const double SFmatrix_err[2] = {1.7731, 1.1545}; //HOFR  bb/cc, and c, in this order  

  //check HF flag 
  if(HFOR>3 || HFOR<0 ){
    return 1.0;       

  }else { 
    
    int HFOR_index = -999;
    if(HFOR==0)     HFOR_index = HFOR;
    else if(HFOR>0) HFOR_index = HFOR-1;//bb and cc events have the same scale factor

    if(HFOR_index>1){
      return WJet_SFLF(SFmatrix_err[0], SFmatrix_err[0], SFmatrix_err[1]);
    }
    else{
      return SFmatrix_err[HFOR_index];
    }
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
//statistical variation Wbb and Wcc SF
double WJet_SF_err_StatBBCCUp(int HFOR)
{  
  if (HFOR>3 || HFOR<0){
    return 1.0;
  }
  else {
    double SFmatrix[4]; 
    SFmatrix[0] = WJet_SF(0) + WJet_SF_Staterr(0);
    SFmatrix[1] = WJet_SF(1) + WJet_SF_Staterr(1);
    SFmatrix[2] = WJet_SF(2);
    SFmatrix[3] = WJet_SFLF(SFmatrix[0],SFmatrix[1],SFmatrix[2]);
    return SFmatrix[HFOR];
  }
}
double WJet_SF_err_StatBBCCDown(int HFOR)
{
  if (HFOR>3 || HFOR<0){
    return 1.0;
  }
  else {
    double SFmatrix[4]; 
    SFmatrix[0] = WJet_SF(0) - WJet_SF_Staterr(0);
    SFmatrix[1] = WJet_SF(1) - WJet_SF_Staterr(1);
    SFmatrix[2] = WJet_SF(2);
    SFmatrix[3] = WJet_SFLF(SFmatrix[0],SFmatrix[1],SFmatrix[2]);
    return SFmatrix[HFOR];
  }
}
double WJet_SF_err_StatCUp(int HFOR)
{  
  if (HFOR>3 || HFOR<0){
    return 1.0;
  }
  else {
    double SFmatrix[4]; 
    SFmatrix[0] = WJet_SF(0);
    SFmatrix[1] = WJet_SF(1);
    SFmatrix[2] = WJet_SF(2) + WJet_SF_Staterr(2);
    SFmatrix[3] = WJet_SFLF(SFmatrix[0],SFmatrix[1],SFmatrix[2]);
    return SFmatrix[HFOR];
  }
}
double WJet_SF_err_StatCDown(int HFOR)
{  
  if (HFOR>3 || HFOR<0){
    return 1.0;
  }
  else {
    double SFmatrix[4]; 
    SFmatrix[0] = WJet_SF(0);
    SFmatrix[1] = WJet_SF(1);
    SFmatrix[2] = WJet_SF(2) - WJet_SF_Staterr(2);
    SFmatrix[3] = WJet_SFLF(SFmatrix[0],SFmatrix[1],SFmatrix[2]);
    return SFmatrix[HFOR];
  }
}
//extrapolation uncertanty, uncorrelated
double WJet_SF_err_ExtrBBCCUp(int HFOR)
{  
  if (HFOR>3 || HFOR<0){
    return 1.0;
  }
  else {
    double SFmatrix[4]; 
    SFmatrix[0] = WJet_SF(0)*(1.+0.25);
    SFmatrix[1] = WJet_SF(1)*(1.+0.25);
    SFmatrix[2] = WJet_SF(2);
    SFmatrix[3] = WJet_SFLF(SFmatrix[0],SFmatrix[1],SFmatrix[2]);
    return SFmatrix[HFOR];
  }
}
double WJet_SF_err_ExtrBBCCDown(int HFOR)
{  
  if (HFOR>3 || HFOR<0){
    return 1.0;
  }
  else {
    double SFmatrix[4]; 
    SFmatrix[0] = WJet_SF(0)*(1.-0.25);
    SFmatrix[1] = WJet_SF(1)*(1.-0.25);
    SFmatrix[2] = WJet_SF(2);
    SFmatrix[3] = WJet_SFLF(SFmatrix[0],SFmatrix[1],SFmatrix[2]);
    return SFmatrix[HFOR];
  }
}
double WJet_SF_err_ExtrCUp(int HFOR)
{  
  if (HFOR>3 || HFOR<0){
    return 1.0;
  }
  else {
    double SFmatrix[4]; 
    SFmatrix[0] = WJet_SF(0);
    SFmatrix[1] = WJet_SF(1);
    SFmatrix[2] = WJet_SF(2)*(1.+0.25);
    SFmatrix[3] = WJet_SFLF(SFmatrix[0],SFmatrix[1],SFmatrix[2]);
    return SFmatrix[HFOR];
  }
}
double WJet_SF_err_ExtrCDown(int HFOR)
{  
  if (HFOR>3 || HFOR<0){
    return 1.0;
  }
  else {
    double SFmatrix[4]; 
    SFmatrix[0] = WJet_SF(0);
    SFmatrix[1] = WJet_SF(1);
    SFmatrix[2] = WJet_SF(2)*(1.-0.25);
    SFmatrix[3] = WJet_SFLF(SFmatrix[0],SFmatrix[1],SFmatrix[2]);
    return SFmatrix[HFOR];
  }
}
//Statistical and extrapolation variation added in quadrature
double WJet_SF_err_StatExtrBBCCUp(int HFOR)
{  
  if (HFOR>3 || HFOR<0){
    return 1.0;
  }
  else {
    double SFmatrix[4]; 
    double err_stat = WJet_SF_Staterr(0)/WJet_SF(0);
    double err_extr = 0.25;
    double err = err_stat*err_stat + err_extr*err_extr;
    err = sqrt(err);
    SFmatrix[0] = WJet_SF(0) * (1+err);
    SFmatrix[1] = WJet_SF(1) * (1+err);
    SFmatrix[2] = WJet_SF(2);
    SFmatrix[3] = WJet_SFLF(SFmatrix[0],SFmatrix[1],SFmatrix[2]);
    return SFmatrix[HFOR];
  }
}
double WJet_SF_err_StatExtrBBCCDown(int HFOR)
{
  if (HFOR>3 || HFOR<0){
    return 1.0;
  }
  else {
    double SFmatrix[4]; 
    double err_stat = WJet_SF_Staterr(0)/WJet_SF(0);
    double err_extr = 0.25;
    double err = err_stat*err_stat + err_extr*err_extr;
    err = sqrt(err);
    SFmatrix[0] = WJet_SF(0) * (1-err);
    SFmatrix[1] = WJet_SF(1) * (1-err);
    SFmatrix[2] = WJet_SF(2);
    SFmatrix[3] = WJet_SFLF(SFmatrix[0],SFmatrix[1],SFmatrix[2]);
    return SFmatrix[HFOR];
  }
}
double WJet_SF_err_StatExtrCUp(int HFOR)
{  
  if (HFOR>3 || HFOR<0){
    return 1.0;
  }
  else {
    double SFmatrix[4]; 
    double err_stat = WJet_SF_Staterr(2)/WJet_SF(2);
    double err_extr = 0.25;
    double err = err_stat*err_stat + err_extr*err_extr;
    err = sqrt(err);
    SFmatrix[0] = WJet_SF(0);
    SFmatrix[1] = WJet_SF(1);
    SFmatrix[2] = WJet_SF(2) * (1+err);
    SFmatrix[3] = WJet_SFLF(SFmatrix[0],SFmatrix[1],SFmatrix[2]);
    return SFmatrix[HFOR];
  }
}
double WJet_SF_err_StatExtrCDown(int HFOR)
{  
  if (HFOR>3 || HFOR<0){
    return 1.0;
  }
  else {
    double SFmatrix[4]; 
    double err_stat = WJet_SF_Staterr(2)/WJet_SF(2);
    double err_extr = 0.25;
    double err = err_stat*err_stat + err_extr*err_extr;
    err = sqrt(err);
    SFmatrix[0] = WJet_SF(0);
    SFmatrix[1] = WJet_SF(1);
    SFmatrix[2] = WJet_SF(2) * (1-err);
    SFmatrix[3] = WJet_SFLF(SFmatrix[0],SFmatrix[1],SFmatrix[2]);
    return SFmatrix[HFOR];
  }
}
//combined uncertainty on Btag and Lqtag efficiency
double WJet_SF_err_BLqtagSFUp(int HFOR)
{
  if(HFOR>3 || HFOR<0 ){
    return 1.0;       
  }
  else { 
    double SFmatrix[4];
    double err_Btag[2] = {WJet_SF_err_BtagSFUp(0)-WJet_SF(0),WJet_SF_err_BtagSFUp(2)-WJet_SF(2)};
    double err_Lqtag[2] = {WJet_SF_err_LqtagSFUp(0)-WJet_SF(0),WJet_SF_err_LqtagSFUp(2)-WJet_SF(2)};
    SFmatrix[0] = WJet_SF(0) - sqrt(err_Btag[0]*err_Btag[0] + err_Lqtag[0]*err_Lqtag[0]);
    SFmatrix[1] = SFmatrix[0];
    SFmatrix[2] = WJet_SF(2) - sqrt(err_Btag[1]*err_Btag[1] + err_Lqtag[1]*err_Lqtag[1]);
    SFmatrix[3] = WJet_SFLF(SFmatrix[0], SFmatrix[1], SFmatrix[2]);
    return SFmatrix[HFOR];
  }
}
double WJet_SF_err_BLqtagSFDown(int HFOR)
{
  if(HFOR>3 || HFOR<0 ){
    return 1.0;       
  }
  else { 
    double SFmatrix[4];
    double err_Btag[2] = {WJet_SF_err_BtagSFDown(0)-WJet_SF(0),WJet_SF_err_BtagSFDown(2)-WJet_SF(2)};
    double err_Lqtag[2] = {WJet_SF_err_LqtagSFDown(0)-WJet_SF(0),WJet_SF_err_LqtagSFDown(2)-WJet_SF(2)};
    SFmatrix[0] = WJet_SF(0) + sqrt(err_Btag[0]*err_Btag[0] + err_Lqtag[0]*err_Lqtag[0]);
    SFmatrix[1] = SFmatrix[0];
    SFmatrix[2] = WJet_SF(2) + sqrt(err_Btag[1]*err_Btag[1] + err_Lqtag[1]*err_Lqtag[1]);
    SFmatrix[3] = WJet_SFLF(SFmatrix[0], SFmatrix[1], SFmatrix[2]);
    return SFmatrix[HFOR];
  }
}
#endif
