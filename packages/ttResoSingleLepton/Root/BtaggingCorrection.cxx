#include "ttResoSingleLepton/BtaggingCorrection.h"

#include <iostream>
#include <TMath.h>

BtaggingCorrection::BtaggingCorrection(BtaggingCorrection::display t_display) : m_display(t_display){ }

BtaggingCorrection::~BtaggingCorrection(){ }

std::pair<double, double> BtaggingCorrection::GetBtaggingCorrection_Eigenvector_GeV(double jet_pT, double jet_eta, double jet_jfit_deltaR, int flavor, bool isBtagged, double BtagDataEffiency, double BtagMCEffiency, double SFin_err_up, double SFin_err_dw){
  if(m_display>=VERBOSE)std::cout << "GetBtaggingCorrection_GeV1: pT= " << jet_pT << "\teta= " << jet_eta << "\tflavor= " << flavor << std::endl;
  
  if(BtagMCEffiency<=0.) {
  	std::cerr << "GetBtaggingCorrection_GeV_v2 : WARNING! strange value for BtagMCEffiency" << BtagMCEffiency << std::endl;
  	return std::pair<double, double>(SFin_err_up, SFin_err_dw);
  }
  if(isBtagged)return std::pair<double, double>(SFin_err_up,  SFin_err_dw); // nothing to be done
  
  double k=1;
  const double D = BtagDataEffiency;
  const double M = BtagMCEffiency;
  const double SFin = ( 1 - D ) / ( 1 - M );
  double SFout_err_up = SFin_err_up;
  double SFout_err_dw = SFin_err_dw;
  
  if(m_display>=VERBOSE)std::cout << "GetBtaggingCorrection_GeV2: pT= " << jet_pT << "\teta= " << jet_eta << "\tflavor= " << flavor << std::endl;
  
  switch (flavor){
  	case 5:
  		if(m_display>=VERBOSE)std::cout << "GetBtaggingCorrection_GeV2.5: pT= " << jet_pT << "\teta= " << jet_eta << "\tflavor= " << flavor << std::endl;
		if(jet_jfit_deltaR<0.4){
		   if (TMath::Abs(jet_eta)<2.1) k = (      (jet_pT<100)*( (((1.72513+(-0.0110655*jet_pT))+(jet_pT*(jet_pT*6.77799e-05)))*TMath::Gaus(jet_jfit_deltaR,(((0.0804822+(jet_pT*-0.0036982))+(jet_pT*(jet_pT*5.25122e-05)))+(jet_pT*(jet_pT*(jet_pT*-2.78081e-07)))),((0.260559+(jet_pT*-0.00222124))+(jet_pT*(jet_pT*1.51114e-05))))>0.1)*((1.72513+(-0.0110655*jet_pT))+(jet_pT*(jet_pT*6.77799e-05)))*TMath::Gaus(jet_jfit_deltaR,(((0.0804822+(jet_pT*-0.0036982))+(jet_pT*(jet_pT*5.25122e-05)))+(jet_pT*(jet_pT*(jet_pT*-2.78081e-07)))),((0.260559+(jet_pT*-0.00222124))+(jet_pT*(jet_pT*1.51114e-05)))) + (((1.72513+(-0.0110655*jet_pT))+(jet_pT*(jet_pT*6.77799e-05)))*TMath::Gaus(jet_jfit_deltaR,(((0.0804822+(jet_pT*-0.0036982))+(jet_pT*(jet_pT*5.25122e-05)))+(jet_pT*(jet_pT*(jet_pT*-2.78081e-07)))),((0.260559+(jet_pT*-0.00222124))+(jet_pT*(jet_pT*1.51114e-05))))<=0.1) )   + (100<jet_pT&&jet_pT<600)*( (((((1.47049+(jet_pT*-0.00290711))+(jet_pT*(jet_pT*1.06591e-05)))+(jet_pT*(jet_pT*(jet_pT*-1.29257e-08))))+(jet_jfit_deltaR*((-2.14819+(jet_pT*-0.00486483))+(jet_pT*(jet_pT*1.29238e-05)))))+(jet_jfit_deltaR*(jet_jfit_deltaR*((-4.43455+(jet_pT*-0.00141062))+(jet_pT*(jet_pT*-2.37868e-06)))))>0.1)*((((1.47049+(jet_pT*-0.00290711))+(jet_pT*(jet_pT*1.06591e-05)))+(jet_pT*(jet_pT*(jet_pT*-1.29257e-08))))+(jet_jfit_deltaR*((-2.14819+(jet_pT*-0.00486483))+(jet_pT*(jet_pT*1.29238e-05)))))+(jet_jfit_deltaR*(jet_jfit_deltaR*((-4.43455+(jet_pT*-0.00141062))+(jet_pT*(jet_pT*-2.37868e-06))))) + (((((1.47049+(jet_pT*-0.00290711))+(jet_pT*(jet_pT*1.06591e-05)))+(jet_pT*(jet_pT*(jet_pT*-1.29257e-08))))+(jet_jfit_deltaR*((-2.14819+(jet_pT*-0.00486483))+(jet_pT*(jet_pT*1.29238e-05)))))+(jet_jfit_deltaR*(jet_jfit_deltaR*((-4.43455+(jet_pT*-0.00141062))+(jet_pT*(jet_pT*-2.37868e-06)))))<=0.1) )   + (600<jet_pT)*( (((((1.47049+(600*-0.00290711))+(600*(600*1.06591e-05)))+(600*(600*(600*-1.29257e-08))))+(jet_jfit_deltaR*((-2.14819+(600*-0.00486483))+(600*(600*1.29238e-05)))))+(jet_jfit_deltaR*(jet_jfit_deltaR*((-4.43455+(600*-0.00141062))+(600*(600*-2.37868e-06)))))>0.1)*((((1.47049+(600*-0.00290711))+(600*(600*1.06591e-05)))+(600*(600*(600*-1.29257e-08))))+(jet_jfit_deltaR*((-2.14819+(600*-0.00486483))+(600*(600*1.29238e-05)))))+(jet_jfit_deltaR*(jet_jfit_deltaR*((-4.43455+(600*-0.00141062))+(600*(600*-2.37868e-06))))) + (((((1.47049+(600*-0.00290711))+(600*(600*1.06591e-05)))+(600*(600*(600*-1.29257e-08))))+(jet_jfit_deltaR*((-2.14819+(600*-0.00486483))+(600*(600*1.29238e-05)))))+(jet_jfit_deltaR*(jet_jfit_deltaR*((-4.43455+(600*-0.00141062))+(600*(600*-2.37868e-06)))))<=0.1) ) );
		   else 			k = (      (jet_pT<80)*( ((1.9111)*TMath::Gaus(jet_jfit_deltaR,(0.0115839+(-0.00183693*jet_pT)),(0.187359))>0.1)*(1.9111)*TMath::Gaus(jet_jfit_deltaR,(0.0115839+(-0.00183693*jet_pT)),(0.187359)) + ((1.9111)*TMath::Gaus(jet_jfit_deltaR,(0.0115839+(-0.00183693*jet_pT)),(0.187359))<=0.1) )   + (80<jet_pT&&jet_pT<200)*( (((1.71334+(jet_pT*-0.00220366))+(jet_jfit_deltaR*(-9.05431+(jet_pT*0.0229291))))+(jet_jfit_deltaR*(jet_jfit_deltaR*(19.3483+(jet_pT*-0.113171))))>0.1)*((1.71334+(jet_pT*-0.00220366))+(jet_jfit_deltaR*(-9.05431+(jet_pT*0.0229291))))+(jet_jfit_deltaR*(jet_jfit_deltaR*(19.3483+(jet_pT*-0.113171)))) + (((1.71334+(jet_pT*-0.00220366))+(jet_jfit_deltaR*(-9.05431+(jet_pT*0.0229291))))+(jet_jfit_deltaR*(jet_jfit_deltaR*(19.3483+(jet_pT*-0.113171))))<=0.1) )   + (200<jet_pT)*( (((1.71334+(200*-0.00220366))+(jet_jfit_deltaR*(-9.05431+(200*0.0229291))))+(jet_jfit_deltaR*(jet_jfit_deltaR*(19.3483+(200*-0.113171))))>0.1)*((1.71334+(200*-0.00220366))+(jet_jfit_deltaR*(-9.05431+(200*0.0229291))))+(jet_jfit_deltaR*(jet_jfit_deltaR*(19.3483+(200*-0.113171)))) + (((1.71334+(200*-0.00220366))+(jet_jfit_deltaR*(-9.05431+(200*0.0229291))))+(jet_jfit_deltaR*(jet_jfit_deltaR*(19.3483+(200*-0.113171))))<=0.1) ) );

		}
		else{
			k = (      (jet_pT<200)*( ((0.0322498+(jet_pT*0.000650857))+(jet_eta*((-0.000229652+(jet_pT*-2.14801e-05))+(jet_pT*(jet_pT*1.00641e-08)))))+(jet_eta*(jet_eta*(-0.00381115+(jet_pT*-2.6435e-05)))))   + (200<jet_pT&&jet_pT<1000)*( (0.196118+(jet_pT*-0.000165288)))   + (1000<jet_pT)*((0.196118+(1000*-0.000165288))) );
		}
		;
		 break;
  	case 4:
  		if(m_display>=VERBOSE)std::cout << "GetBtaggingCorrection_GeV2.4: pT= " << jet_pT << "\teta= " << jet_eta << "\tflavor= " << flavor << std::endl;
		if(jet_jfit_deltaR<0.4){
		   if (TMath::Abs(jet_eta)<2.1) k = (      (jet_pT<90)*( (TMath::Exp((1.80712+(jet_pT*-0.00859894))+((-5.54601+(jet_pT*-0.0454658))*jet_jfit_deltaR))>0.1)*TMath::Exp((1.80712+(jet_pT*-0.00859894))+((-5.54601+(jet_pT*-0.0454658))*jet_jfit_deltaR)) + (TMath::Exp((1.80712+(jet_pT*-0.00859894))+((-5.54601+(jet_pT*-0.0454658))*jet_jfit_deltaR))<=0.1) )   + (90<jet_pT&&jet_pT<1200)*( (TMath::Exp((((1.25305+(jet_pT*-0.00233846))+(jet_pT*(jet_pT*3.22361e-06)))+(jet_pT*(jet_pT*(jet_pT*-1.95202e-09))))+(((-10.5329+(jet_pT*0.0150788))+(jet_pT*(jet_pT*-2.35955e-05)))*jet_jfit_deltaR))>0.1)*TMath::Exp((((1.25305+(jet_pT*-0.00233846))+(jet_pT*(jet_pT*3.22361e-06)))+(jet_pT*(jet_pT*(jet_pT*-1.95202e-09))))+(((-10.5329+(jet_pT*0.0150788))+(jet_pT*(jet_pT*-2.35955e-05)))*jet_jfit_deltaR)) + (TMath::Exp((((1.25305+(jet_pT*-0.00233846))+(jet_pT*(jet_pT*3.22361e-06)))+(jet_pT*(jet_pT*(jet_pT*-1.95202e-09))))+(((-10.5329+(jet_pT*0.0150788))+(jet_pT*(jet_pT*-2.35955e-05)))*jet_jfit_deltaR))<=0.1) )   + (1200<jet_pT)*( (TMath::Exp((((1.25305+(1200*-0.00233846))+(1200*(1200*3.22361e-06)))+(1200*(1200*(1200*-1.95202e-09))))+(((-10.5329+(1200*0.0150788))+(1200*(1200*-2.35955e-05)))*jet_jfit_deltaR))>0.1)*TMath::Exp((((1.25305+(1200*-0.00233846))+(1200*(1200*3.22361e-06)))+(1200*(1200*(1200*-1.95202e-09))))+(((-10.5329+(1200*0.0150788))+(1200*(1200*-2.35955e-05)))*jet_jfit_deltaR)) + (TMath::Exp((((1.25305+(1200*-0.00233846))+(1200*(1200*3.22361e-06)))+(1200*(1200*(1200*-1.95202e-09))))+(((-10.5329+(1200*0.0150788))+(1200*(1200*-2.35955e-05)))*jet_jfit_deltaR))<=0.1) ) );
		   else 			k = (      (jet_pT<80)*( (TMath::Exp((2.47119+(jet_pT*-0.0161987))+((-6.60286+(jet_pT*-0.121647))*jet_jfit_deltaR))>0.1)*TMath::Exp((2.47119+(jet_pT*-0.0161987))+((-6.60286+(jet_pT*-0.121647))*jet_jfit_deltaR)) + (TMath::Exp((2.47119+(jet_pT*-0.0161987))+((-6.60286+(jet_pT*-0.121647))*jet_jfit_deltaR))<=0.1) )   + (80<jet_pT&&jet_pT<400)*( (TMath::Exp((0.897403+(jet_pT*0.0014461))+((1.04715+(jet_pT*-0.107781))*jet_jfit_deltaR))>0.1)*TMath::Exp((0.897403+(jet_pT*0.0014461))+((1.04715+(jet_pT*-0.107781))*jet_jfit_deltaR)) + (TMath::Exp((0.897403+(jet_pT*0.0014461))+((1.04715+(jet_pT*-0.107781))*jet_jfit_deltaR))<=0.1) )   + (400<jet_pT)*( (TMath::Exp((0.897403+(400*0.0014461))+((1.04715+(400*-0.107781))*jet_jfit_deltaR))>0.1)*TMath::Exp((0.897403+(400*0.0014461))+((1.04715+(400*-0.107781))*jet_jfit_deltaR)) + (TMath::Exp((0.897403+(400*0.0014461))+((1.04715+(400*-0.107781))*jet_jfit_deltaR))<=0.1) ) );
		}
		else{
			k = (      (jet_pT<200)*( (0.0169573+(jet_pT*0.000376827)))   + (200<jet_pT&&jet_pT<1000)*( (0.118225+(jet_pT*-0.000115122)))   + (1000<jet_pT)*((0.118225+(1000*-0.000115122))) );
		}
		;
		 break;
	default: break;
  	
  }
  if(m_display>=VERBOSE)std::cout << "GetBtaggingCorrection_GeV3: pT= " << jet_pT << "\teta= " << jet_eta << "\tisBtagged= " << isBtagged << std::endl;
  
  if(!isBtagged){
  	    if( k * M < 1 ) {
	
		double SFout = ( 1 - k * D ) / ( 1 - k * M );
		if(SFout<=0) {
			if(m_display>=INFO)std::cout << "Bad btagging correction 2 : " << k <<"\t-->newSFineff= " << SFout <<"\twith : D= " << D <<"\tM=" << M << "\t Going to use k=1 --> SFineff=" << ( 1 - D ) / ( 1 - M )  <<std::endl;
			if(m_display>=INFO)std::cout << "                          : flavor=" << flavor << " pt= " << jet_pT << "\teta= "<< jet_eta << "\tDR= " << jet_jfit_deltaR << std::endl;
			SFout = ( 1 - D ) / ( 1 - M );
			k = 1;
		}
		
		// SFin_err is SFineff_err = ( M * SFeff_err ) / ( 1 - M)
		SFout_err_up = SFout + (SFin_err_up-SFin) * k * (1 - M) / ( 1 - k * M ) ;
		SFout_err_dw = SFout - (SFin-SFin_err_dw) * k * (1 - M) / ( 1 - k * M ) ;
		
		if(SFout_err_up<0) {
			if(m_display>=INFO)std::cout << "Bad btagging correction 3: " << k <<"\t-->SF_err= " << SFout_err_up <<" and newSFineff= " << SFout << " (was : " << (1-D)/(1-M) <<" )\twith : D= " << D <<"\tM=" << M << "\t Going to use old error --> SFineff_err=" << (SFin_err_up * M)/fabs(1 - M) << std::endl;
			if(m_display>=INFO)std::cout << "                         : flavor=" << flavor << " pt= " << jet_pT << "\teta= "<< jet_eta << "\tDR= " << jet_jfit_deltaR << std::endl;
			k=1.;
			SFout_err_up = SFin_err_up;
		}
		else{
			if(m_display>=DEBUG)std::cout << "Good btagging correction 2: " << k <<"\t-->newSFineff_err= " << SFout_err_up << " (was : " << SFin_err_up <<" )\t SF=" << SFout << " (was : " << SFin << ")\twith : D= " << D <<"\tM=" << M << std::endl;
			if(m_display>=DEBUG)std::cout << "                          : flavor=" << flavor << " pt= " << jet_pT << "\teta= "<< jet_eta << "\tDR= " << jet_jfit_deltaR << std::endl;
		}
		if(SFout_err_dw<0) {
			if(m_display>=INFO)std::cout << "Bad btagging correction 3: " << k <<"\t-->SF_err= " << SFout_err_dw <<" and newSFineff= " << SFout << " (was : " << (1-D)/(1-M) <<" )\twith : D= " << D <<"\tM=" << M << "\t Going to use old error --> SFineff_err=" << (SFin_err_dw * M)/fabs(1 - M) << std::endl;
			if(m_display>=INFO)std::cout << "                         : flavor=" << flavor << " pt= " << jet_pT << "\teta= "<< jet_eta << "\tDR= " << jet_jfit_deltaR << std::endl;
			k=1.;
			SFout_err_dw = SFin_err_dw;
		}
		else{
			if(m_display>=DEBUG)std::cout << "Good btagging correction 2: " << k <<"\t-->newSFineff_err= " << SFout_err_dw << " (was : " << SFin_err_dw <<" )\t SF=" << SFout << " (was : " << SFin << ")\twith : D= " << D <<"\tM=" << M << std::endl;
			if(m_display>=DEBUG)std::cout << "                          : flavor=" << flavor << " pt= " << jet_pT << "\teta= "<< jet_eta << "\tDR= " << jet_jfit_deltaR << std::endl;
		}

	   }
	   else {
			if(m_display>=DEBUG)std::cout << "Bad btagging correction 4: " << k  <<" )\twith : D= " << D <<"\tM=" << M << "\t Going to use old SF" << std::endl;
	   
	   }
  }
  
  if(m_display>=VERBOSE)std::cout << "GetBtaggingCorrection_GeV: SFup=" <<SFout_err_up  << "\tSFdw= " << SFout_err_dw << std::endl;
  
  
  return std::pair<double, double>(SFout_err_up, SFout_err_dw);

}

std::pair<double, double> BtaggingCorrection::GetBtaggingCorrection_GeV(double jet_pT, double jet_eta, double jet_jfit_deltaR, int flavor, bool isBtagged, double SFeff, double SFeff_err, double BtagDataEffiency, double BtagMCEffiency){
  if(m_display>=VERBOSE)std::cout << "GetBtaggingCorrection_GeV1: pT= " << jet_pT << "\teta= " << jet_eta << "\tflavor= " << flavor << std::endl;
  
  if(isBtagged)return std::pair<double, double>(SFeff, SFeff_err); // nothing to be done
  
  double k=1;
  const double D = BtagDataEffiency;
  const double M = BtagMCEffiency;
  
  double SFout = SFeff;
  double SFout_err = SFeff_err;
  if(m_display>=VERBOSE)std::cout << "GetBtaggingCorrection_GeV2: pT= " << jet_pT << "\teta= " << jet_eta << "\tflavor= " << flavor << std::endl;
  
  switch (flavor){
  	case 5:
  		if(m_display>=VERBOSE)std::cout << "GetBtaggingCorrection_GeV2.5: pT= " << jet_pT << "\teta= " << jet_eta << "\tflavor= " << flavor << std::endl;
		if(jet_jfit_deltaR<0.4){
		   if (TMath::Abs(jet_eta)<2.1) k = (      (jet_pT<100)*( (((1.72513+(-0.0110655*jet_pT))+(jet_pT*(jet_pT*6.77799e-05)))*TMath::Gaus(jet_jfit_deltaR,(((0.0804822+(jet_pT*-0.0036982))+(jet_pT*(jet_pT*5.25122e-05)))+(jet_pT*(jet_pT*(jet_pT*-2.78081e-07)))),((0.260559+(jet_pT*-0.00222124))+(jet_pT*(jet_pT*1.51114e-05))))>0.1)*((1.72513+(-0.0110655*jet_pT))+(jet_pT*(jet_pT*6.77799e-05)))*TMath::Gaus(jet_jfit_deltaR,(((0.0804822+(jet_pT*-0.0036982))+(jet_pT*(jet_pT*5.25122e-05)))+(jet_pT*(jet_pT*(jet_pT*-2.78081e-07)))),((0.260559+(jet_pT*-0.00222124))+(jet_pT*(jet_pT*1.51114e-05)))) + (((1.72513+(-0.0110655*jet_pT))+(jet_pT*(jet_pT*6.77799e-05)))*TMath::Gaus(jet_jfit_deltaR,(((0.0804822+(jet_pT*-0.0036982))+(jet_pT*(jet_pT*5.25122e-05)))+(jet_pT*(jet_pT*(jet_pT*-2.78081e-07)))),((0.260559+(jet_pT*-0.00222124))+(jet_pT*(jet_pT*1.51114e-05))))<=0.1) )   + (100<jet_pT&&jet_pT<600)*( (((((1.47049+(jet_pT*-0.00290711))+(jet_pT*(jet_pT*1.06591e-05)))+(jet_pT*(jet_pT*(jet_pT*-1.29257e-08))))+(jet_jfit_deltaR*((-2.14819+(jet_pT*-0.00486483))+(jet_pT*(jet_pT*1.29238e-05)))))+(jet_jfit_deltaR*(jet_jfit_deltaR*((-4.43455+(jet_pT*-0.00141062))+(jet_pT*(jet_pT*-2.37868e-06)))))>0.1)*((((1.47049+(jet_pT*-0.00290711))+(jet_pT*(jet_pT*1.06591e-05)))+(jet_pT*(jet_pT*(jet_pT*-1.29257e-08))))+(jet_jfit_deltaR*((-2.14819+(jet_pT*-0.00486483))+(jet_pT*(jet_pT*1.29238e-05)))))+(jet_jfit_deltaR*(jet_jfit_deltaR*((-4.43455+(jet_pT*-0.00141062))+(jet_pT*(jet_pT*-2.37868e-06))))) + (((((1.47049+(jet_pT*-0.00290711))+(jet_pT*(jet_pT*1.06591e-05)))+(jet_pT*(jet_pT*(jet_pT*-1.29257e-08))))+(jet_jfit_deltaR*((-2.14819+(jet_pT*-0.00486483))+(jet_pT*(jet_pT*1.29238e-05)))))+(jet_jfit_deltaR*(jet_jfit_deltaR*((-4.43455+(jet_pT*-0.00141062))+(jet_pT*(jet_pT*-2.37868e-06)))))<=0.1) )   + (600<jet_pT)*( (((((1.47049+(600*-0.00290711))+(600*(600*1.06591e-05)))+(600*(600*(600*-1.29257e-08))))+(jet_jfit_deltaR*((-2.14819+(600*-0.00486483))+(600*(600*1.29238e-05)))))+(jet_jfit_deltaR*(jet_jfit_deltaR*((-4.43455+(600*-0.00141062))+(600*(600*-2.37868e-06)))))>0.1)*((((1.47049+(600*-0.00290711))+(600*(600*1.06591e-05)))+(600*(600*(600*-1.29257e-08))))+(jet_jfit_deltaR*((-2.14819+(600*-0.00486483))+(600*(600*1.29238e-05)))))+(jet_jfit_deltaR*(jet_jfit_deltaR*((-4.43455+(600*-0.00141062))+(600*(600*-2.37868e-06))))) + (((((1.47049+(600*-0.00290711))+(600*(600*1.06591e-05)))+(600*(600*(600*-1.29257e-08))))+(jet_jfit_deltaR*((-2.14819+(600*-0.00486483))+(600*(600*1.29238e-05)))))+(jet_jfit_deltaR*(jet_jfit_deltaR*((-4.43455+(600*-0.00141062))+(600*(600*-2.37868e-06)))))<=0.1) ) );
		   else 			k = (      (jet_pT<80)*( ((1.9111)*TMath::Gaus(jet_jfit_deltaR,(0.0115839+(-0.00183693*jet_pT)),(0.187359))>0.1)*(1.9111)*TMath::Gaus(jet_jfit_deltaR,(0.0115839+(-0.00183693*jet_pT)),(0.187359)) + ((1.9111)*TMath::Gaus(jet_jfit_deltaR,(0.0115839+(-0.00183693*jet_pT)),(0.187359))<=0.1) )   + (80<jet_pT&&jet_pT<200)*( (((1.71334+(jet_pT*-0.00220366))+(jet_jfit_deltaR*(-9.05431+(jet_pT*0.0229291))))+(jet_jfit_deltaR*(jet_jfit_deltaR*(19.3483+(jet_pT*-0.113171))))>0.1)*((1.71334+(jet_pT*-0.00220366))+(jet_jfit_deltaR*(-9.05431+(jet_pT*0.0229291))))+(jet_jfit_deltaR*(jet_jfit_deltaR*(19.3483+(jet_pT*-0.113171)))) + (((1.71334+(jet_pT*-0.00220366))+(jet_jfit_deltaR*(-9.05431+(jet_pT*0.0229291))))+(jet_jfit_deltaR*(jet_jfit_deltaR*(19.3483+(jet_pT*-0.113171))))<=0.1) )   + (200<jet_pT)*( (((1.71334+(200*-0.00220366))+(jet_jfit_deltaR*(-9.05431+(200*0.0229291))))+(jet_jfit_deltaR*(jet_jfit_deltaR*(19.3483+(200*-0.113171))))>0.1)*((1.71334+(200*-0.00220366))+(jet_jfit_deltaR*(-9.05431+(200*0.0229291))))+(jet_jfit_deltaR*(jet_jfit_deltaR*(19.3483+(200*-0.113171)))) + (((1.71334+(200*-0.00220366))+(jet_jfit_deltaR*(-9.05431+(200*0.0229291))))+(jet_jfit_deltaR*(jet_jfit_deltaR*(19.3483+(200*-0.113171))))<=0.1) ) );

		}
		else{
			k = (      (jet_pT<200)*( ((0.0322498+(jet_pT*0.000650857))+(jet_eta*((-0.000229652+(jet_pT*-2.14801e-05))+(jet_pT*(jet_pT*1.00641e-08)))))+(jet_eta*(jet_eta*(-0.00381115+(jet_pT*-2.6435e-05)))))   + (200<jet_pT&&jet_pT<1000)*( (0.196118+(jet_pT*-0.000165288)))   + (1000<jet_pT)*((0.196118+(1000*-0.000165288))) );
		}
		;
		 break;
  	case 4:
  		if(m_display>=VERBOSE)std::cout << "GetBtaggingCorrection_GeV2.4: pT= " << jet_pT << "\teta= " << jet_eta << "\tflavor= " << flavor << std::endl;
		if(jet_jfit_deltaR<0.4){
		   if (TMath::Abs(jet_eta)<2.1) k = (      (jet_pT<90)*( (TMath::Exp((1.80712+(jet_pT*-0.00859894))+((-5.54601+(jet_pT*-0.0454658))*jet_jfit_deltaR))>0.1)*TMath::Exp((1.80712+(jet_pT*-0.00859894))+((-5.54601+(jet_pT*-0.0454658))*jet_jfit_deltaR)) + (TMath::Exp((1.80712+(jet_pT*-0.00859894))+((-5.54601+(jet_pT*-0.0454658))*jet_jfit_deltaR))<=0.1) )   + (90<jet_pT&&jet_pT<1200)*( (TMath::Exp((((1.25305+(jet_pT*-0.00233846))+(jet_pT*(jet_pT*3.22361e-06)))+(jet_pT*(jet_pT*(jet_pT*-1.95202e-09))))+(((-10.5329+(jet_pT*0.0150788))+(jet_pT*(jet_pT*-2.35955e-05)))*jet_jfit_deltaR))>0.1)*TMath::Exp((((1.25305+(jet_pT*-0.00233846))+(jet_pT*(jet_pT*3.22361e-06)))+(jet_pT*(jet_pT*(jet_pT*-1.95202e-09))))+(((-10.5329+(jet_pT*0.0150788))+(jet_pT*(jet_pT*-2.35955e-05)))*jet_jfit_deltaR)) + (TMath::Exp((((1.25305+(jet_pT*-0.00233846))+(jet_pT*(jet_pT*3.22361e-06)))+(jet_pT*(jet_pT*(jet_pT*-1.95202e-09))))+(((-10.5329+(jet_pT*0.0150788))+(jet_pT*(jet_pT*-2.35955e-05)))*jet_jfit_deltaR))<=0.1) )   + (1200<jet_pT)*( (TMath::Exp((((1.25305+(1200*-0.00233846))+(1200*(1200*3.22361e-06)))+(1200*(1200*(1200*-1.95202e-09))))+(((-10.5329+(1200*0.0150788))+(1200*(1200*-2.35955e-05)))*jet_jfit_deltaR))>0.1)*TMath::Exp((((1.25305+(1200*-0.00233846))+(1200*(1200*3.22361e-06)))+(1200*(1200*(1200*-1.95202e-09))))+(((-10.5329+(1200*0.0150788))+(1200*(1200*-2.35955e-05)))*jet_jfit_deltaR)) + (TMath::Exp((((1.25305+(1200*-0.00233846))+(1200*(1200*3.22361e-06)))+(1200*(1200*(1200*-1.95202e-09))))+(((-10.5329+(1200*0.0150788))+(1200*(1200*-2.35955e-05)))*jet_jfit_deltaR))<=0.1) ) );
		   else 			k = (      (jet_pT<80)*( (TMath::Exp((2.47119+(jet_pT*-0.0161987))+((-6.60286+(jet_pT*-0.121647))*jet_jfit_deltaR))>0.1)*TMath::Exp((2.47119+(jet_pT*-0.0161987))+((-6.60286+(jet_pT*-0.121647))*jet_jfit_deltaR)) + (TMath::Exp((2.47119+(jet_pT*-0.0161987))+((-6.60286+(jet_pT*-0.121647))*jet_jfit_deltaR))<=0.1) )   + (80<jet_pT&&jet_pT<400)*( (TMath::Exp((0.897403+(jet_pT*0.0014461))+((1.04715+(jet_pT*-0.107781))*jet_jfit_deltaR))>0.1)*TMath::Exp((0.897403+(jet_pT*0.0014461))+((1.04715+(jet_pT*-0.107781))*jet_jfit_deltaR)) + (TMath::Exp((0.897403+(jet_pT*0.0014461))+((1.04715+(jet_pT*-0.107781))*jet_jfit_deltaR))<=0.1) )   + (400<jet_pT)*( (TMath::Exp((0.897403+(400*0.0014461))+((1.04715+(400*-0.107781))*jet_jfit_deltaR))>0.1)*TMath::Exp((0.897403+(400*0.0014461))+((1.04715+(400*-0.107781))*jet_jfit_deltaR)) + (TMath::Exp((0.897403+(400*0.0014461))+((1.04715+(400*-0.107781))*jet_jfit_deltaR))<=0.1) ) );
		}
		else{
			k = (      (jet_pT<200)*( (0.0169573+(jet_pT*0.000376827)))   + (200<jet_pT&&jet_pT<1000)*( (0.118225+(jet_pT*-0.000115122)))   + (1000<jet_pT)*((0.118225+(1000*-0.000115122))) );
		}
		;
		 break;
	default: break;
  	
  }
  if(m_display>=VERBOSE)std::cout << "GetBtaggingCorrection_GeV3: pT= " << jet_pT << "\teta= " << jet_eta << "\tisBtagged= " << isBtagged << std::endl;
  
  if(!isBtagged){
  	    if( k * M < 1 ) {
	
		SFout = ( 1 - k * D ) / ( 1 - k * M );
		if(SFout<=0) {
			if(m_display>=INFO)std::cout << "Bad btagging correction 2 : " << k <<"\t-->newSFineff= " << SFout <<"\twith : D= " << D <<"\tM=" << M << "\t Going to use k=1 --> SFineff=" << ( 1 - D ) / ( 1 - M )  <<std::endl;
			if(m_display>=INFO)std::cout << "                          : flavor=" << flavor << " pt= " << jet_pT << "\teta= "<< jet_eta << "\tDR= " << jet_jfit_deltaR << std::endl;
			SFout = ( 1 - D ) / ( 1 - M );
			k = 1;
		}
		
		// SFeff_err is actually SFineff_err = ( M * SFeff_err ) / ( 1 - M)
		SFout_err = SFeff_err * k * (1 - M) / ( 1 - k * M ) ;
		
		if(SFout_err<0) {
			if(m_display>=INFO)std::cout << "Bad btagging correction 3: " << k <<"\t-->SF_err= " << SFout_err <<" and newSFineff= " << SFout << " (was : " << (1-D)/(1-M) <<" )\twith : D= " << D <<"\tM=" << M << "\t Going to use old error --> SFineff_err=" << (SFeff_err * M)/fabs(1 - M) << std::endl;
			if(m_display>=INFO)std::cout << "                         : flavor=" << flavor << " pt= " << jet_pT << "\teta= "<< jet_eta << "\tDR= " << jet_jfit_deltaR << std::endl;
			k=1.;
			SFout_err = SFeff_err;
		}
		else{
			if(m_display>=DEBUG)std::cout << "Good btagging correction 2: " << k <<"\t-->newSFineff_err= " << SFout_err << " (was : " << SFeff_err <<" )\t SF=" << SFout << " (was : " << SFeff << ")\twith : D= " << D <<"\tM=" << M << std::endl;
			if(m_display>=DEBUG)std::cout << "                          : flavor=" << flavor << " pt= " << jet_pT << "\teta= "<< jet_eta << "\tDR= " << jet_jfit_deltaR << std::endl;
		}

	   }
	   else {
			if(m_display>=DEBUG)std::cout << "Bad btagging correction 4: " << k  <<" )\twith : D= " << D <<"\tM=" << M << "\t Going to use old SF: SFineff=" << SFout << std::endl;
	   
	   }
  }
  
  if(m_display>=VERBOSE)std::cout << "GetBtaggingCorrection_GeV: SF=" << SFout << "\tSFerr= " << SFout_err << std::endl;
  
  
  return std::pair<double, double>(SFout, SFout_err);

}

std::pair<double, double> BtaggingCorrection::GetBtaggingCorrection_Eigenvector_MeV(double jet_pT, double jet_eta, double jet_jfit_deltaR, int flavor, bool isBtagged, double BtagDataEffiency, double BtagMCEffiency, double SFin_err_up, double SFin_err_dw){
 
 return GetBtaggingCorrection_Eigenvector_GeV(jet_pT/1000., jet_eta, jet_jfit_deltaR, flavor, isBtagged, BtagDataEffiency, BtagMCEffiency, SFin_err_up, SFin_err_dw);

}

std::pair<double, double> BtaggingCorrection::GetBtaggingCorrection_MeV(double jet_pT, double jet_eta, double jet_jfit_deltaR, int flavor, bool isBtagged, double SFeff, double SFeff_err, double BtagDataEffiency, double BtagMCEffiency){
 
 return GetBtaggingCorrection_GeV(jet_pT/1000., jet_eta, jet_jfit_deltaR, flavor, isBtagged, SFeff, SFeff_err, BtagDataEffiency, BtagMCEffiency);

}
