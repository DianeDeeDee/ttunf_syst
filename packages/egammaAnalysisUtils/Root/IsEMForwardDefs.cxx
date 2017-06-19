#include "egammaAnalysisUtils/IsEMForwardDefs.h"

#include <TMath.h>
#include <iostream>

// 2011 optimised cuts
// Define cuts for EMEC and FCAL - [4 bins in pile-up][for the 6 var's]
const unsigned int nvtx_bins_2011 = 4;
const double cuts_nvtx_2011[nvtx_bins_2011+1] = {-0.5,3.5,6.5,10.5,1e200};
const unsigned int nvariables     = 6;

// 2^0: <pho_{E}^{1}>: not used
// 2^1: <lambda^{2}>; 2^2: lateral; 2^3: longitudinal;
// 2^4: f_{max}; 2^5: <R^{2}>; 2^6: lambda_{Center}
// 2^7 - out of eta range
// 2^8 - out of pVtx range; should not happen

const double Cuts_EMEC_Loose_2011[nvtx_bins_2011][nvariables] =
  { { 4000,.60,.49,-.22,3500,255},
    { 4200,.64,.48,-.22,3500,255},
    { 4200,.69,.48,-.18,3900,260},
    { 4500,.69,.55,-.22,3900,255} };
const double Cuts_FCAL_Loose_2011[nvtx_bins_2011][nvariables] = 
  { { 6700,.38,.38,-.28, 900,250},
    { 7900,.40,.49,-.28,1100,250},
    { 8900,.42,.59,-.28,1400,246},
    {10000,.63,.73,-.23,1900,252} };

const double Cuts_EMEC_Medium_2011[nvtx_bins_2011][nvariables] =
  { { 3600,.54,.25,-.26,2600,255},
    { 4200,.58,.26,-.25,2700,255},
    { 4000,.62,.27,-.24,3000,250},
    { 4500,.64,.29,-.23,3300,255} };
const double Cuts_FCAL_Medium_2011[nvtx_bins_2011][nvariables] =
  { { 4300,.38,.28,-.30, 900,244},
    { 5500,.40,.33,-.29, 950,246},
    { 6600,.42,.40,-.29,1050,246},
    { 9700,.45,.52,-.26,1200,250} };

const double Cuts_EMEC_Tight_2011[nvtx_bins_2011][nvariables] = 
  { { 3000,.48,.25,-.26,2200,255},
    { 3000,.54,.25,-.41,2300,250},
    { 3200,.55,.27,-.26,2600,250},
    { 2800,.64,.24,-.39,3000,250} };
const double Cuts_FCAL_Tight_2011[nvtx_bins_2011][nvariables] = 
  { { 3000,.38,.22,-.30, 850,243},
    { 4000,.40,.27,-.30, 950,244},
    { 4500,.42,.31,-.29,1000,246},
    { 7500,.42,.41,-.28,1100,250} };

bool isForward_Loose(int vxp_n, float el_cl_eta,
		     float el_secondlambda,float el_lateral,float el_longitudinal,
		     float el_cellmaxfrac, float el_secondR,float el_centerlambda,
		     egammaForwardMenu::egForwardMenu menu)
{
  switch (menu) {
  case egammaForwardMenu::eg2012:
    return (identifyForwardElectron2012(el_cellmaxfrac, el_longitudinal, el_secondlambda,
					el_lateral, el_secondR, el_centerlambda, el_cl_eta, vxp_n) >= 1);
  case egammaForwardMenu::eg2011:
    return (Forward_IsEM(vxp_n, el_cl_eta,
			 el_secondlambda, el_lateral, el_longitudinal,
			 el_cellmaxfrac, el_secondR, el_centerlambda,
			 egammaForwardMenu::Loose, menu) == 0);
  default:
    // should bever be reached
    std::cout << "Undefined Forward IsEM menu\n";
    return false;
  }
}

bool isForward_Medium(int vxp_n, float el_cl_eta,
		      float el_secondlambda,float el_lateral,float el_longitudinal,
		      float el_cellmaxfrac, float el_secondR,float el_centerlambda,
		      egammaForwardMenu::egForwardMenu menu)
{
  switch (menu) {
  case egammaForwardMenu::eg2012:
    return (identifyForwardElectron2012(el_cellmaxfrac, el_longitudinal, el_secondlambda,
					el_lateral, el_secondR, el_centerlambda, el_cl_eta, vxp_n) >= 3);
  case egammaForwardMenu::eg2011:
    return (Forward_IsEM(vxp_n, el_cl_eta,
			 el_secondlambda, el_lateral, el_longitudinal,
			 el_cellmaxfrac, el_secondR, el_centerlambda,
			 egammaForwardMenu::Medium, menu) == 0);
  default:
    // should bever be reached
    std::cout << "Undefined Forward IsEM menu\n";
    return false;
  }
}

bool isForward_Tight(int vxp_n, float el_cl_eta,
		     float el_secondlambda,float el_lateral,float el_longitudinal,
		     float el_cellmaxfrac, float el_secondR,float el_centerlambda,
		     egammaForwardMenu::egForwardMenu menu)
{
  switch (menu) {
  case egammaForwardMenu::eg2012:
    return (identifyForwardElectron2012(el_cellmaxfrac, el_longitudinal, el_secondlambda,
					el_lateral, el_secondR, el_centerlambda, el_cl_eta, vxp_n) >= 7);
  case egammaForwardMenu::eg2011:
    return (Forward_IsEM(vxp_n, el_cl_eta,
			 el_secondlambda, el_lateral, el_longitudinal,
			 el_cellmaxfrac, el_secondR, el_centerlambda,
			 egammaForwardMenu::Tight, menu) == 0);
  default:
    // should bever be reached
    std::cout << "Undefined Forward IsEM menu\n";
    return false;
  }
}

unsigned int Forward_IsEM(int vxp_n, float el_cl_eta,
			  float el_secondlambda,float el_lateral,float el_longitudinal,
			  float el_cellmaxfrac, float el_secondR,float el_centerlambda,
			  egammaForwardMenu::egForwardCut cut,
			  egammaForwardMenu::egForwardMenu menu)
{
  // put variables in array to allow access via index
  double variables[nvariables] = {el_secondlambda,el_lateral,el_longitudinal,
				  -el_cellmaxfrac,el_secondR,el_centerlambda};

  int emec_fcal_bin = getFwdEtaBin(el_cl_eta, menu);
  if (emec_fcal_bin < 0)
    return 128;

  int nvtx_bin = getNvtxBin(vxp_n, menu);
  if (nvtx_bin < 0)
    return 256;

  if (menu == egammaForwardMenu::eg2011) {
    if (emec_fcal_bin == 0 && cut == egammaForwardMenu::Loose) {
      return Forward_ID(Cuts_EMEC_Loose_2011[nvtx_bin], variables);
    } else if (emec_fcal_bin == 1 && cut == egammaForwardMenu::Loose) {
      return Forward_ID(Cuts_FCAL_Loose_2011[nvtx_bin], variables);
      
    } else if (emec_fcal_bin == 0 && cut == egammaForwardMenu::Medium) {
      return Forward_ID(Cuts_EMEC_Medium_2011[nvtx_bin], variables);
    } else if (emec_fcal_bin == 1 && cut == egammaForwardMenu::Medium) {
      return Forward_ID(Cuts_FCAL_Medium_2011[nvtx_bin], variables);

    } else if (emec_fcal_bin == 0 && cut == egammaForwardMenu::Tight) {
      return Forward_ID(Cuts_EMEC_Tight_2011[nvtx_bin], variables);
    } else if (emec_fcal_bin == 1 && cut == egammaForwardMenu::Tight) {
      return Forward_ID(Cuts_FCAL_Tight_2011[nvtx_bin], variables);
    }
    
  }

  // non-2011 period selected or other strange configuration
  return 512;
}

int getNvtxBin(int vxp_n,
	       egammaForwardMenu::egForwardMenu menu)
{
  // get nvertex bin number, non-2011 defaults right now to bin number=0
  if (menu == egammaForwardMenu::eg2011) {
    for (unsigned int i=0; i<nvtx_bins_2011; i++){
      if ( vxp_n>cuts_nvtx_2011[i] && vxp_n<cuts_nvtx_2011[i+1] ) {
        return i;
      }
    }
    return -1;
  } else {
    return 0;
  }
}


int getFwdEtaBin(float el_cl_eta,
		 egammaForwardMenu::egForwardMenu menu)
{
  // Forward eta bin choice
  // coded as year dependendt, but actually it is not (yet)

  el_cl_eta = TMath::Abs(el_cl_eta);

  if (menu == egammaForwardMenu::eg2011) {
    if ( el_cl_eta>=2.5 && el_cl_eta<3.2 ) {
      // Electron in EMEC
      return 0;
    } else if ( el_cl_eta>=3.2 && el_cl_eta<4.9 ) {
      // Electron in FCAL
      return 1;
    }
  } else {
    if ( el_cl_eta>=2.5 && el_cl_eta<3.2 ) {
      // Electron in EMEC
      return 0;
    } else if ( el_cl_eta>=3.2 && el_cl_eta<4.9 ) {
      // Electron in FCAL
      return 1;
    }
  }
  return -1;
}

unsigned int Forward_ID (const double *cuts, const double *variables)
{
  // Returns IsEM style bit mask for the cuts defined by the various variables
  // 2^1: First variable; 2^2: Second variable; ...
  // bit = 1 for failed cut

  unsigned int BinOut=0;
  for (unsigned int j=0; j<nvariables; j++) {
    if ( variables[j] > cuts[j] ) {
      BinOut |= 0x0001 << (j+1);
    }
  }
  return BinOut;
}

int identifyForwardElectron2012(double frac, double lon, double sL, double lat, double sR, double cL, double eta, double npvtx)
{
	//Identify forward electron by binning in eta and npvtx, then switch to correct value for mva, then set return bits accordingly
	int switchVar=0;
	double looseCut,mediumCut,tightCut;
	double mva;

	eta = fabs(eta);
	if (eta <= 2.5) switchVar = 0;
	else if(eta <= 2.6) switchVar = 10;
	else if(eta <= 2.7) switchVar = 20;
	else if(eta <= 2.8) switchVar = 30;
	else if(eta <= 2.9) switchVar = 40;
	else if(eta <= 3.0) switchVar = 50;
	else if(eta <= 3.2) switchVar = 60;
	else if(eta <= 3.6) switchVar = 70;
	else if(eta <= 4.0) switchVar = 80;
	else if(eta <= 4.9) switchVar = 90;
	else switchVar = 0;

	if(npvtx < 6.5) switchVar += 1;
	else if(npvtx < 10.5) switchVar += 2;
	else if(npvtx < 13.5) switchVar += 3;
	else if(npvtx < 16.5) switchVar += 4;
	else if(npvtx < 19.5) switchVar += 5;
	else switchVar += 6;


		switch(switchVar)
		{
		case 11:
			//Determined at: Mon Jul 8 12:11:02 CEST 2013
			looseCut = 0.122708;
			mediumCut = 0.192463;
			tightCut = 0.242224;
			//These values correspond to background efficiencies of 4.24152% (@90.2337%), 2.72614% (@80.1935%), 1.90595% (@70.4089%)
			mva = (0.876963767733) +  frac*(-0.0624124249142) + lon*(-0.138468346605) + sL*(9.26747013451e-07) + lat*(-1.53764447032) + sR*(1.47216226637e-05) + cL*(-3.85816681824e-05);
		break;

		case 21:
			//Determined at: Mon Jul 8 12:11:08 CEST 2013
			looseCut = 0.0915239;
			mediumCut = 0.166847;
			tightCut = 0.216093;
			//These values correspond to background efficiencies of 2.7254% (@90.2383%), 1.62184% (@80.2599%), 1.2376% (@70.4305%)
			mva = (0.991170596184) +  frac*(-0.131513542599) + lon*(-0.0597599328482) + sL*(6.59456121774e-07) + lat*(-1.53298933175) + sR*(1.14948645643e-05) + cL*(-7.53782650298e-05);
		break;

		case 31:
			//Determined at: Mon Jul 8 12:11:15 CEST 2013
			looseCut = 0.0925174;
			mediumCut = 0.171629;
			tightCut = 0.226572;
			//These values correspond to background efficiencies of 3.468% (@90.2347%), 2.04659% (@80.2521%), 1.38012% (@70.4434%)
			mva = (1.16404329161) +  frac*(-0.198484763252) + lon*(-0.0678343028792) + sL*(6.70423459426e-07) + lat*(-1.72564963799) + sR*(1.09061380176e-05) + cL*(-7.56685755552e-05);
		break;

		case 41:
			//Determined at: Mon Jul 8 12:11:21 CEST 2013
			looseCut = 0.113614;
			mediumCut = 0.18373;
			tightCut = 0.228498;
			//These values correspond to background efficiencies of 2.94658% (@90.2376%), 1.88835% (@80.2505%), 1.40042% (@70.4399%)
			mva = (0.992087117841) +  frac*(-0.135503538624) + lon*(-0.0677394380253) + sL*(6.20758969948e-07) + lat*(-1.58326653881) + sR*(1.53105989283e-05) + cL*(-7.12785016075e-05);
		break;

		case 51:
			//Determined at: Mon Jul 8 12:11:28 CEST 2013
			looseCut = 0.134251;
			mediumCut = 0.206829;
			tightCut = 0.24935;
			//These values correspond to background efficiencies of 3.18457% (@90.2311%), 2.00131% (@80.2171%), 1.46812% (@70.4307%)
			mva = (0.968250456191) +  frac*(-0.14249618878) + lon*(-0.0757201749792) + sL*(6.45185849215e-07) + lat*(-1.66336383884) + sR*(2.98142037611e-05) + cL*(-0.000104492286648);
		break;

		case 61:
			//Determined at: Mon Jul 8 12:11:33 CEST 2013
			looseCut = 0.098628;
			mediumCut = 0.192883;
			tightCut = 0.242083;
			//These values correspond to background efficiencies of 9.99638% (@90.2339%), 6.17157% (@80.192%), 4.56081% (@70.3818%)
			mva = (0.811248002175) +  frac*(0.0349886957416) + lon*(-0.288372120803) + sL*(1.09378654827e-06) + lat*(-1.7687716983) + sR*(8.96818768695e-05) + cL*(-0.000245811841717);
		break;

		case 71:
			//Determined at: Mon Jul 8 12:11:40 CEST 2013
			looseCut = 0.162322;
			mediumCut = 0.219539;
			tightCut = 0.256869;
			//These values correspond to background efficiencies of 3.16202% (@83.4459%), 2.04276% (@73.5%), 1.49568% (@64.2027%)
			mva = (0.540305263197) +  frac*(0.102078795589) + lon*(-0.327337234705) + sL*(1.61159405884e-06) + lat*(-1.40519752653) + sR*(0.000244035914856) + cL*(-0.000364305302652);
		break;

		case 81:
			//Determined at: Mon Jul 8 12:11:49 CEST 2013
			looseCut = 0.18715;
			mediumCut = 0.23296;
			tightCut = 0.26023;
			//These values correspond to background efficiencies of 1.84597% (@83.4489%), 1.24496% (@73.4848%), 0.953035% (@64.1885%)
			mva = (0.359676121586) +  frac*(0.181579764304) + lon*(-0.382507620656) + sL*(1.92245115016e-06) + lat*(-1.14373597842) + sR*(0.000255923316224) + cL*(-0.000339202852779);
		break;

		case 91:
			//Determined at: Mon Jul 8 12:11:57 CEST 2013
			looseCut = 0.142903;
			mediumCut = 0.191653;
			tightCut = 0.220188;
			//These values correspond to background efficiencies of 3.09474% (@83.4304%), 2.02105% (@73.4844%), 1.37895% (@64.2071%)
			mva = (0.195589700186) +  frac*(0.261268926473) + lon*(-0.507940775141) + sL*(1.00936843551e-06) + lat*(-0.303950171783) + sR*(9.21168563385e-05) + cL*(-0.000272161204899);
		break;

		case 12:
			//Determined at: Mon Jul 8 12:12:04 CEST 2013
			looseCut = 0.0969765;
			mediumCut = 0.170095;
			tightCut = 0.220594;
			//These values correspond to background efficiencies of 5.12503% (@90.19%), 3.27247% (@80.3849%), 2.36716% (@70.6491%)
			mva = (0.846227872176) +  frac*(-0.0342412185352) + lon*(-0.133526274972) + sL*(8.99273137154e-07) + lat*(-1.4744228414) + sR*(1.39723818052e-05) + cL*(-4.72608944054e-05);
		break;

		case 22:
			//Determined at: Mon Jul 8 12:12:18 CEST 2013
			looseCut = 0.0754831;
			mediumCut = 0.149073;
			tightCut = 0.197828;
			//These values correspond to background efficiencies of 3.11684% (@90.198%), 1.84804% (@80.3711%), 1.2977% (@70.637%)
			mva = (0.934838608788) +  frac*(-0.0750042007021) + lon*(-0.0570312754878) + sL*(6.54627070407e-07) + lat*(-1.45318943577) + sR*(1.03420067983e-05) + cL*(-8.02708947512e-05);
		break;

		case 32:
			//Determined at: Mon Jul 8 12:12:39 CEST 2013
			looseCut = 0.0821713;
			mediumCut = 0.155119;
			tightCut = 0.206532;
			//These values correspond to background efficiencies of 3.513% (@90.1882%), 2.14329% (@80.392%), 1.51528% (@70.6658%)
			mva = (1.07270217508) +  frac*(-0.112454108194) + lon*(-0.0637364125313) + sL*(6.34877152127e-07) + lat*(-1.6056700155) + sR*(9.7848783804e-06) + cL*(-7.73222254584e-05);
		break;

		case 42:
			//Determined at: Mon Jul 8 12:12:56 CEST 2013
			looseCut = 0.0943547;
			mediumCut = 0.16392;
			tightCut = 0.210093;
			//These values correspond to background efficiencies of 3.43148% (@90.187%), 2.08638% (@80.3869%), 1.49305% (@70.6596%)
			mva = (0.948770382993) +  frac*(-0.0908296170906) + lon*(-0.0697194429395) + sL*(6.18339176576e-07) + lat*(-1.50775738803) + sR*(1.54646464184e-05) + cL*(-8.81829228469e-05);
		break;

		case 52:
			//Determined at: Mon Jul 8 12:13:12 CEST 2013
			looseCut = 0.104919;
			mediumCut = 0.182135;
			tightCut = 0.226944;
			//These values correspond to background efficiencies of 3.57174% (@90.191%), 2.00628% (@80.3681%), 1.39364% (@70.6572%)
			mva = (0.931868674702) +  frac*(-0.102921486525) + lon*(-0.0736289855894) + sL*(6.40719933092e-07) + lat*(-1.58661816213) + sR*(2.81934558899e-05) + cL*(-0.00011909699237);
		break;

		case 62:
			//Determined at: Mon Jul 8 12:13:25 CEST 2013
			looseCut = 0.0729059;
			mediumCut = 0.166212;
			tightCut = 0.219529;
			//These values correspond to background efficiencies of 11.2265% (@90.1931%), 6.9468% (@80.3832%), 4.95008% (@70.6659%)
			mva = (0.748339347461) +  frac*(0.0963144910611) + lon*(-0.228300765723) + sL*(9.51838791374e-07) + lat*(-1.68681203046) + sR*(8.73169274317e-05) + cL*(-0.000255595595366);
		break;

		case 72:
			//Determined at: Mon Jul 8 12:13:44 CEST 2013
			looseCut = 0.145965;
			mediumCut = 0.204785;
			tightCut = 0.241523;
			//These values correspond to background efficiencies of 3.56423% (@82.3599%), 2.33523% (@72.21%), 1.71417% (@62.9192%)
			mva = (0.51051352291) +  frac*(0.131688087264) + lon*(-0.318034158544) + sL*(1.52323832415e-06) + lat*(-1.29493852409) + sR*(0.000214110087352) + cL*(-0.000340884283335);
		break;

		case 82:
			//Determined at: Mon Jul 8 12:14:00 CEST 2013
			looseCut = 0.163874;
			mediumCut = 0.213235;
			tightCut = 0.242104;
			//These values correspond to background efficiencies of 2.41392% (@82.3519%), 1.53924% (@72.2326%), 1.10728% (@62.9114%)
			mva = (0.351870955981) +  frac*(0.191722077154) + lon*(-0.346064847547) + sL*(1.7076458452e-06) + lat*(-1.05105661567) + sR*(0.000218643059098) + cL*(-0.00033066021321);
		break;

		case 92:
			//Determined at: Mon Jul 8 12:14:25 CEST 2013
			looseCut = 0.102398;
			mediumCut = 0.150302;
			tightCut = 0.183735;
			//These values correspond to background efficiencies of 3.39192% (@82.3516%), 2.14639% (@72.2462%), 1.39048% (@62.9124%)
			mva = (0.173271300875) +  frac*(0.274069827082) + lon*(-0.442558407018) + sL*(8.77600383326e-07) + lat*(-0.247154040918) + sR*(7.66255932899e-05) + cL*(-0.000275851616272);
		break;

		case 13:
			//Determined at: Mon Jul 8 12:14:45 CEST 2013
			looseCut = 0.0812017;
			mediumCut = 0.15342;
			tightCut = 0.200732;
			//These values correspond to background efficiencies of 5.38491% (@89.9527%), 3.39502% (@80.1482%), 2.39686% (@70.2989%)
			mva = (0.807750030614) +  frac*(-0.0237481320981) + lon*(-0.128070060405) + sL*(8.66081718723e-07) + lat*(-1.39662906464) + sR*(1.28746189658e-05) + cL*(-4.63442294123e-05);
		break;

		case 23:
			//Determined at: Mon Jul 8 12:14:58 CEST 2013
			looseCut = 0.0631263;
			mediumCut = 0.132661;
			tightCut = 0.181534;
			//These values correspond to background efficiencies of 3.23542% (@89.9412%), 1.84412% (@80.1493%), 1.22941% (@70.3302%)
			mva = (0.888197814402) +  frac*(-0.0541079742581) + lon*(-0.0549678354625) + sL*(6.40985514567e-07) + lat*(-1.37121340058) + sR*(8.99256555319e-06) + cL*(-7.71409259602e-05);
		break;

		case 33:
			//Determined at: Mon Jul 8 12:15:15 CEST 2013
			looseCut = 0.0694034;
			mediumCut = 0.13596;
			tightCut = 0.185739;
			//These values correspond to background efficiencies of 3.48455% (@89.9339%), 2.10983% (@80.1559%), 1.40655% (@70.3183%)
			mva = (0.997921001403) +  frac*(-0.0656780321652) + lon*(-0.0567083168506) + sL*(6.0773327017e-07) + lat*(-1.49944401718) + sR*(8.88931739623e-06) + cL*(-7.82105505005e-05);
		break;

		case 43:
			//Determined at: Mon Jul 8 12:15:30 CEST 2013
			looseCut = 0.0776768;
			mediumCut = 0.146736;
			tightCut = 0.193339;
			//These values correspond to background efficiencies of 3.7377% (@89.9464%), 2.29329% (@80.1498%), 1.62702% (@70.304%)
			mva = (0.907075010543) +  frac*(-0.0797003650421) + lon*(-0.0718796704557) + sL*(5.63636711772e-07) + lat*(-1.4174661927) + sR*(1.33496106643e-05) + cL*(-7.84753480769e-05);
		break;

		case 53:
			//Determined at: Mon Jul 8 12:15:43 CEST 2013
			looseCut = 0.0888762;
			mediumCut = 0.162011;
			tightCut = 0.208104;
			//These values correspond to background efficiencies of 4.13902% (@89.9502%), 2.48571% (@80.1552%), 1.68783% (@70.3069%)
			mva = (0.874158661946) +  frac*(-0.06192628762) + lon*(-0.0714700029372) + sL*(6.19606677666e-07) + lat*(-1.48155342819) + sR*(2.50140950669e-05) + cL*(-0.000116716084872);
		break;

		case 63:
			//Determined at: Mon Jul 8 12:15:57 CEST 2013
			looseCut = 0.0536402;
			mediumCut = 0.146147;
			tightCut = 0.199172;
			//These values correspond to background efficiencies of 11.7856% (@89.9513%), 7.13184% (@80.1592%), 5.071% (@70.3106%)
			mva = (0.679211206296) +  frac*(0.150522221532) + lon*(-0.214788217326) + sL*(9.29056749766e-07) + lat*(-1.54632561927) + sR*(7.76222671825e-05) + cL*(-0.000248076069685);
		break;

		case 73:
			//Determined at: Mon Jul 8 12:16:13 CEST 2013
			looseCut = 0.131403;
			mediumCut = 0.185235;
			tightCut = 0.223506;
			//These values correspond to background efficiencies of 3.71591% (@82.0845%), 2.46343% (@72.2067%), 1.72362% (@62.4434%)
			mva = (0.470380008098) +  frac*(0.15770053768) + lon*(-0.288859041465) + sL*(1.42674751492e-06) + lat*(-1.19245264904) + sR*(0.000193195503845) + cL*(-0.000336760051633);
		break;

		case 83:
			//Determined at: Mon Jul 8 12:16:29 CEST 2013
			looseCut = 0.133179;
			mediumCut = 0.182709;
			tightCut = 0.215964;
			//These values correspond to background efficiencies of 2.64779% (@82.0814%), 1.63504% (@72.1874%), 1.12711% (@62.5109%)
			mva = (0.328120441138) +  frac*(0.19149617877) + lon*(-0.309813154135) + sL*(1.46499664818e-06) + lat*(-0.941973288656) + sR*(0.000186648813775) + cL*(-0.000301831216721);
		break;

		case 93:
			//Determined at: Mon Jul 8 12:16:49 CEST 2013
			looseCut = 0.0770711;
			mediumCut = 0.12391;
			tightCut = 0.156041;
			//These values correspond to background efficiencies of 3.74281% (@82.0841%), 2.3177% (@72.2048%), 1.51335% (@62.5155%)
			mva = (0.150282158834) +  frac*(0.274136209437) + lon*(-0.385059434183) + sL*(7.36083864861e-07) + lat*(-0.193976120709) + sR*(6.01768245261e-05) + cL*(-0.000252538925823);
		break;

		case 14:
			//Determined at: Mon Jul 8 12:17:08 CEST 2013
			looseCut = 0.0745222;
			mediumCut = 0.139701;
			tightCut = 0.186164;
			//These values correspond to background efficiencies of 5.57842% (@89.8576%), 3.39996% (@79.989%), 2.41361% (@70.2738%)
			mva = (0.786931327272) +  frac*(-0.0204759471869) + lon*(-0.123680508086) + sL*(8.4687738698e-07) + lat*(-1.34064260639) + sR*(1.17532544693e-05) + cL*(-4.70036500677e-05);
		break;

		case 24:
			//Determined at: Mon Jul 8 12:17:17 CEST 2013
			looseCut = 0.052924;
			mediumCut = 0.118817;
			tightCut = 0.163188;
			//These values correspond to background efficiencies of 3.40278% (@89.8946%), 1.93995% (@79.9985%), 1.36321% (@70.3115%)
			mva = (0.845791141118) +  frac*(-0.0404669797191) + lon*(-0.0499778243721) + sL*(6.25934399542e-07) + lat*(-1.30071463724) + sR*(8.25124351466e-06) + cL*(-7.62182419349e-05);
		break;

		case 34:
			//Determined at: Mon Jul 8 12:17:30 CEST 2013
			looseCut = 0.0551299;
			mediumCut = 0.120911;
			tightCut = 0.169049;
			//These values correspond to background efficiencies of 3.91232% (@89.8912%), 2.2369% (@79.9966%), 1.49576% (@70.3075%)
			mva = (0.952807601494) +  frac*(-0.057747621141) + lon*(-0.0553039142205) + sL*(5.80192660776e-07) + lat*(-1.41237587503) + sR*(7.28408923572e-06) + cL*(-7.03741143517e-05);
		break;

		case 44:
			//Determined at: Mon Jul 8 12:17:40 CEST 2013
			looseCut = 0.0613888;
			mediumCut = 0.125127;
			tightCut = 0.171309;
			//These values correspond to background efficiencies of 4.13007% (@89.893%), 2.41511% (@79.9943%), 1.63203% (@70.3135%)
			mva = (0.844434009779) +  frac*(-0.0465127704727) + lon*(-0.0639552395354) + sL*(5.39977599223e-07) + lat*(-1.32302720902) + sR*(1.21699002361e-05) + cL*(-7.84234766595e-05);
		break;

		case 54:
			//Determined at: Mon Jul 8 12:17:50 CEST 2013
			looseCut = 0.0724987;
			mediumCut = 0.147672;
			tightCut = 0.193411;
			//These values correspond to background efficiencies of 4.5882% (@89.8965%), 2.42131% (@80.0082%), 1.63517% (@70.2941%)
			mva = (0.832754994539) +  frac*(-0.0394728362206) + lon*(-0.0733909599271) + sL*(5.7680940883e-07) + lat*(-1.39433434981) + sR*(2.29595430567e-05) + cL*(-0.000114068539248);
		break;

		case 64:
			//Determined at: Mon Jul 8 12:17:59 CEST 2013
			looseCut = 0.0447151;
			mediumCut = 0.132454;
			tightCut = 0.183492;
			//These values correspond to background efficiencies of 12.413% (@89.8688%), 7.61585% (@79.9961%), 5.47877% (@70.2591%)
			mva = (0.642422297613) +  frac*(0.171882248548) + lon*(-0.182217561159) + sL*(8.7937909845e-07) + lat*(-1.48292649436) + sR*(7.81564255971e-05) + cL*(-0.000267546869005);
		break;

		case 74:
			//Determined at: Mon Jul 8 12:18:11 CEST 2013
			looseCut = 0.11414;
			mediumCut = 0.171161;
			tightCut = 0.20547;
			//These values correspond to background efficiencies of 3.82899% (@81.5532%), 2.39432% (@70.7695%), 1.73084% (@61.4947%)
			mva = (0.438334376211) +  frac*(0.168513599809) + lon*(-0.270930473098) + sL*(1.28097585099e-06) + lat*(-1.08808066499) + sR*(0.000171516758683) + cL*(-0.000314555912186);
		break;

		case 84:
			//Determined at: Mon Jul 8 12:18:23 CEST 2013
			looseCut = 0.113809;
			mediumCut = 0.167273;
			tightCut = 0.198005;
			//These values correspond to background efficiencies of 2.80588% (@81.5596%), 1.58382% (@70.7678%), 1.05735% (@61.4878%)
			mva = (0.304529064676) +  frac*(0.204988256347) + lon*(-0.270548052133) + sL*(1.34957966673e-06) + lat*(-0.847950100221) + sR*(0.00015753251253) + cL*(-0.000300184120103);
		break;

		case 94:
			//Determined at: Mon Jul 8 12:18:40 CEST 2013
			looseCut = 0.061232;
			mediumCut = 0.105595;
			tightCut = 0.13517;
			//These values correspond to background efficiencies of 3.79362% (@81.5727%), 2.22455% (@70.7449%), 1.46174% (@61.5124%)
			mva = (0.132841973725) +  frac*(0.264998624544) + lon*(-0.328056710248) + sL*(6.5018927045e-07) + lat*(-0.166253205988) + sR*(5.11373207107e-05) + cL*(-0.000241898124685);
		break;

		case 15:
			//Determined at: Mon Jul 8 12:18:55 CEST 2013
			looseCut = 0.058091;
			mediumCut = 0.122189;
			tightCut = 0.168602;
			//These values correspond to background efficiencies of 6.0423% (@89.5923%), 3.80029% (@79.9568%), 2.58388% (@70.1668%)
			mva = (0.733077216445) +  frac*(0.0218135436197) + lon*(-0.102931660529) + sL*(7.9844234061e-07) + lat*(-1.27428334884) + sR*(1.19662822014e-05) + cL*(-5.757069163e-05);
		break;

		case 25:
			//Determined at: Mon Jul 8 12:19:01 CEST 2013
			looseCut = 0.0437049;
			mediumCut = 0.102874;
			tightCut = 0.147955;
			//These values correspond to background efficiencies of 4.01987% (@89.5852%), 2.28524% (@79.9563%), 1.40758% (@70.1747%)
			mva = (0.800826391422) +  frac*(-0.0241532807936) + lon*(-0.0452774784568) + sL*(5.85314861197e-07) + lat*(-1.22253102962) + sR*(7.26084313236e-06) + cL*(-7.03914445488e-05);
		break;

		case 35:
			//Determined at: Mon Jul 8 12:19:08 CEST 2013
			looseCut = 0.0476752;
			mediumCut = 0.106326;
			tightCut = 0.146797;
			//These values correspond to background efficiencies of 3.67317% (@89.6066%), 2.27261% (@79.8932%), 1.63311% (@70.1554%)
			mva = (0.900823904355) +  frac*(-0.0541833715351) + lon*(-0.0465740596588) + sL*(5.45978970342e-07) + lat*(-1.32106789472) + sR*(5.94923812309e-06) + cL*(-5.79818576736e-05);
		break;

		case 45:
			//Determined at: Mon Jul 8 12:19:14 CEST 2013
			looseCut = 0.0520327;
			mediumCut = 0.11554;
			tightCut = 0.160059;
			//These values correspond to background efficiencies of 4.58415% (@89.6271%), 2.67457% (@79.937%), 1.75429% (@70.1681%)
			mva = (0.803985836488) +  frac*(-0.0366210844606) + lon*(-0.0750122778923) + sL*(5.25991886577e-07) + lat*(-1.2295198758) + sR*(9.75059199108e-06) + cL*(-6.51039797594e-05);
		break;

		case 55:
			//Determined at: Mon Jul 8 12:19:20 CEST 2013
			looseCut = 0.0505588;
			mediumCut = 0.119957;
			tightCut = 0.16535;
			//These values correspond to background efficiencies of 4.83476% (@89.6181%), 2.66558% (@79.9523%), 1.72039% (@70.1372%)
			mva = (0.784951582811) +  frac*(-0.0435680568678) + lon*(-0.0576349363955) + sL*(5.53592070028e-07) + lat*(-1.2929825158) + sR*(1.91845434315e-05) + cL*(-0.000102713521199);
		break;

		case 65:
			//Determined at: Mon Jul 8 12:19:26 CEST 2013
			looseCut = 0.0305308;
			mediumCut = 0.116786;
			tightCut = 0.166101;
			//These values correspond to background efficiencies of 12.6828% (@89.6252%), 7.49156% (@79.9432%), 5.35433% (@70.1723%)
			mva = (0.60238845275) +  frac*(0.196340873155) + lon*(-0.154496538166) + sL*(8.19341379858e-07) + lat*(-1.39467944503) + sR*(7.31447141015e-05) + cL*(-0.00028018274984);
		break;

		case 75:
			//Determined at: Mon Jul 8 12:19:33 CEST 2013
			looseCut = 0.113454;
			mediumCut = 0.167967;
			tightCut = 0.19874;
			//These values correspond to background efficiencies of 3.88851% (@79.1963%), 2.4583% (@67.6732%), 1.73635% (@58.6146%)
			mva = (0.416913906087) +  frac*(0.172276403216) + lon*(-0.250061610152) + sL*(1.15886614236e-06) + lat*(-1.00676360021) + sR*(0.000152419312401) + cL*(-0.000299620158854);
		break;

		case 85:
			//Determined at: Mon Jul 8 12:19:40 CEST 2013
			looseCut = 0.104009;
			mediumCut = 0.151615;
			tightCut = 0.179783;
			//These values correspond to background efficiencies of 2.63733% (@79.2668%), 1.48187% (@67.678%), 0.995528% (@58.7021%)
			mva = (0.272929925991) +  frac*(0.212747407791) + lon*(-0.244713108674) + sL*(1.17172326102e-06) + lat*(-0.728429036623) + sR*(0.000129601118052) + cL*(-0.000272187694872);
		break;

		case 95:
			//Determined at: Mon Jul 8 12:19:48 CEST 2013
			looseCut = 0.0538611;
			mediumCut = 0.0904982;
			tightCut = 0.112329;
			//These values correspond to background efficiencies of 3.25467% (@79.2527%), 1.8244% (@67.6868%), 1.20718% (@58.7189%)
			mva = (0.107209842608) +  frac*(0.254976438662) + lon*(-0.26902856323) + sL*(5.80582207725e-07) + lat*(-0.12334172567) + sR*(3.92563115434e-05) + cL*(-0.00022686082384);
		break;

		case 16:
			//Determined at: Mon Jul 8 12:19:58 CEST 2013
			looseCut = 0.0306351;
			mediumCut = 0.0943776;
			tightCut = 0.135055;
			//These values correspond to background efficiencies of 7.8791% (@90.9677%), 5.11496% (@80%), 3.35831% (@69.5484%)
			mva = (0.661088517729) +  frac*(0.0423356105143) + lon*(-0.112242512469) + sL*(8.23027954869e-07) + lat*(-1.13186600927) + sR*(9.72477392539e-06) + cL*(-5.60114733163e-05);
		break;

		case 26:
			//Determined at: Mon Jul 8 12:20:02 CEST 2013
			looseCut = 0.0183159;
			mediumCut = 0.0885929;
			tightCut = 0.124867;
			//These values correspond to background efficiencies of 4.94929% (@90.9414%), 2.19067% (@80.0178%), 1.4334% (@69.627%)
			mva = (0.757542214267) +  frac*(-0.0236122086563) + lon*(-0.0412087574698) + sL*(5.07784257842e-07) + lat*(-1.14046973612) + sR*(6.61535715577e-06) + cL*(-6.04744817954e-05);
		break;

		case 36:
			//Determined at: Mon Jul 8 12:20:07 CEST 2013
			looseCut = 0.0173461;
			mediumCut = 0.0808189;
			tightCut = 0.124762;
			//These values correspond to background efficiencies of 5.07223% (@90.839%), 2.75772% (@80.0386%), 1.80565% (@69.6239%)
			mva = (0.818419090942) +  frac*(-0.00681473831964) + lon*(-0.0227622508601) + sL*(4.65201698761e-07) + lat*(-1.21095755207) + sR*(5.23007180628e-06) + cL*(-6.09188871245e-05);
		break;

		case 46:
			//Determined at: Mon Jul 8 12:20:12 CEST 2013
			looseCut = 0.0211886;
			mediumCut = 0.0903187;
			tightCut = 0.130861;
			//These values correspond to background efficiencies of 5.59956% (@90.9389%), 2.78136% (@80.0218%), 1.84196% (@69.4323%)
			mva = (0.739747542159) +  frac*(-0.0247689172186) + lon*(-0.066249845201) + sL*(4.51497220147e-07) + lat*(-1.11886724066) + sR*(8.43345696161e-06) + cL*(-5.41242479526e-05);
		break;

		case 56:
			//Determined at: Mon Jul 8 12:20:16 CEST 2013
			looseCut = 0.036909;
			mediumCut = 0.108834;
			tightCut = 0.156308;
			//These values correspond to background efficiencies of 5.44188% (@90.9091%), 2.45973% (@80%), 1.34959% (@69.5758%)
			mva = (0.738212893356) +  frac*(-0.0250700940633) + lon*(-0.0517040746474) + sL*(4.84810321125e-07) + lat*(-1.19181373403) + sR*(1.64744148292e-05) + cL*(-9.48723786628e-05);
		break;

		case 66:
			//Determined at: Mon Jul 8 12:20:21 CEST 2013
			looseCut = 0.0110357;
			mediumCut = 0.10597;
			tightCut = 0.15126;
			//These values correspond to background efficiencies of 13.7378% (@90.9295%), 7.28379% (@80.06%), 5.23695% (@69.5652%)
			mva = (0.5660998571) +  frac*(0.167767773713) + lon*(-0.156156577524) + sL*(7.94104383263e-07) + lat*(-1.25248761069) + sR*(6.63023870846e-05) + cL*(-0.00028041203656);
		break;

		case 76:
			//Determined at: Mon Jul 8 12:20:26 CEST 2013
			looseCut = 0.0868362;
			mediumCut = 0.142746;
			tightCut = 0.169258;
			//These values correspond to background efficiencies of 4.01888% (@82.5655%), 2.35984% (@70.822%), 1.84497% (@60.4336%)
			mva = (0.375790197155) +  frac*(0.160964088344) + lon*(-0.236096105888) + sL*(9.80187824948e-07) + lat*(-0.867426358236) + sR*(0.000128007622043) + cL*(-0.00024992334284);
		break;

		case 86:
			//Determined at: Mon Jul 8 12:20:30 CEST 2013
			looseCut = 0.0641773;
			mediumCut = 0.110752;
			tightCut = 0.146215;
			//These values correspond to background efficiencies of 3.81348% (@82.5509%), 2.07569% (@70.7706%), 1.22611% (@60.5846%)
			mva = (0.227441216978) +  frac*(0.223361738275) + lon*(-0.212647181148) + sL*(9.37120102103e-07) + lat*(-0.561088901978) + sR*(9.30413996175e-05) + cL*(-0.000248903172209);
		break;

		case 96:
			//Determined at: Mon Jul 8 12:20:35 CEST 2013
			looseCut = 0.0332412;
			mediumCut = 0.0624907;
			tightCut = 0.0886686;
			//These values correspond to background efficiencies of 3.48803% (@82.5532%), 1.96821% (@70.7801%), 1.15297% (@60.5674%)
			mva = (0.0967392240618) +  frac*(0.219458184349) + lon*(-0.217044189681) + sL*(4.64730542003e-07) + lat*(-0.089098059736) + sR*(2.80447996118e-05) + cL*(-0.000190343977829);
		break;

		default:
			return 0;
			break;
		}


	//Calculate Fisher MVA and compare with cuts
	if(mva <= looseCut) return 0;
	else if(mva <= mediumCut) return 1;
	else if(mva <= tightCut) return 3;
	else return 7;

}


