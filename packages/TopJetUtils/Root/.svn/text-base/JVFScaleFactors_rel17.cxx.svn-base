#include "TopJetUtils/JVFScaleFactors_rel17.h"

#include <iostream>
#include <cmath>
#include "TMath.h"

using std::fabs;

double topjetutils::JVF_SF(const double pT, const double JVF, const bool matched, const std::string year, const int ps)
{
	double SF = -1.;
	double JVFcut = -1.;
	
	//Defining JVF working point
	if (year=="2011") JVFcut = 0.75;
	else if(year=="2012") JVFcut = 0.50;
	else JVFcut = 0.50;
	
	if (fabs(JVF)>=JVFcut){
		if (matched)
			SF = topjetutils::JVFSasS_SF(pT, year, ps);
		else
			SF = 1.; 
	} else {
		if (matched)
			SF = topjetutils::JVFSasB_SF(pT, year, ps);
		else
			SF = 1.; 
	}
	return SF;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// JVF scale factors functions:
//SasS = hard-scatter jet selection efficiency
//SasB = hard-scatter jet selection inefficiency
//BasB = pile-up rejection selection efficiency
//BasS = pile-up rejection selection inefficiency
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double topjetutils::JVFSasS_PythiaSF(const double /*jet_pT*/){//constant scale factor
	double a = 0.9984;
	return a;
}

double topjetutils::JVFSasS_HerwigSF(const double jet_pT){//exponential scale factor
	double a = 0.03769;
	double b = -0.00957;
       	double c = 0.9913;
	return a*TMath::Exp(b*jet_pT)+c;
}

double topjetutils::JVFSasS_SF(const double jet_pT, const std::string year, const int ps)
{
	if (year=="2011"){
		double a = 0.0298319;
		double b = -0.0119208;
       		double c = 1.00479;
			
		return a*TMath::Exp(b*jet_pT)+c;
	
	}
	else if(year=="2012"){
		double a = 0.0382635;
		double b = -0.0121108;
       		double c = 0.998046;
			
		return a*TMath::Exp(b*jet_pT)+c;
	}
	else{
	
		if(ps==0){//Pythia parton showering
			return(JVFSasS_PythiaSF(jet_pT));
		}
		else if(ps==1){//Herwig parton showering
			return(JVFSasS_HerwigSF(jet_pT));
		}
		else {//Other parton showering - SF is the averaged one
			double SF = (JVFSasS_HerwigSF(jet_pT)+JVFSasS_PythiaSF(jet_pT))/2;
			return(SF);
		}
		
	}
}

double topjetutils::JVFSasB_PythiaSF(const double /*jet_pT*/){
	double a = 1.022;
	return a;
}

double topjetutils::JVFSasB_HerwigSF(const double jet_pT){
	double a = 0.5789;
	double b = -0.03239;
       	double c = 0.5486;
	return a*TMath::Exp(b*jet_pT)+c;
}

double topjetutils::JVFSasB_SF(const double jet_pT, const std::string year, const int ps)
{

	if (year=="2011"){
		double a = 0.436658;
		double b = -0.0167345;
		double c = 0.5376;
		
		return a*TMath::Exp(b*jet_pT)+c;
	
	}
	else if(year=="2012"){
	
		double a = 0.961611;
		double b = -0.0389405;
		double c = 0.00281162;
		double d = 0.362232;
		
		if (jet_pT > 250.)
			return 1.;
		else
			return a*TMath::Exp(b*jet_pT)+c*jet_pT+d;
	
	}
	else{
		if (jet_pT > 250.){
			return 1.;
		}
		else{
			if(ps==0){//Pythia parton showering
				return(JVFSasB_PythiaSF(jet_pT));
			}
			else if(ps==1){//Herwig parton showering
				return(JVFSasB_HerwigSF(jet_pT));
			}
			else {//Other parton showering - SF is the averaged one
				double SF = (JVFSasB_PythiaSF(jet_pT)+JVFSasB_HerwigSF(jet_pT))/2;
				return(SF);
			}
		}
	}
}

double topjetutils::JVFBasS_SF( const double /*jet_pT*/, const std::string /*year*/, const int /*ps*/ )
{
	return 1.;
}


double topjetutils::JVFBasB_SF( const double /*jet_pT*/, const std::string /*year*/, const int /*ps*/ )
{
	return 1.;

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// General tools
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double topjetutils::exp_er(double x, double * par, double *matrix)
{   
  const int npar = 3; 
  
  double derivate[npar] = 
    {
      exp(x*par[1]),
      x*par[0]*exp(x*par[1]),
      1,
    };
  
  float er=0.;
  for(int i=0; i<npar;i++)
    { 
      for(int j=0; j<npar;j++)
	{  
	  er =er+derivate[i]*( *(matrix+npar*i+j) )*derivate[j];
	}
    }
  return sqrt(er);
}

double topjetutils::explinear_er(double x, double * par, double *matrix)
{   
  const int npar = 4; 
  
  double derivate[npar] = 
    {
      exp(x*par[1]),
      x*par[0]*exp(x*par[1]),
      x,
      1,
    };
  
  float er=0.;
  for(int i=0; i<npar;i++)
    { 
      for(int j=0; j<npar;j++)
	{  
	  er =er+derivate[i]*( *(matrix+npar*i+j) )*derivate[j];
	}
    }
  return sqrt(er);
}

double topjetutils::erfPlus_er(double x, double * par, double *matrix)
{   
  const int npar = 3; 
  
  double derivate[npar] = 
    {
      par[1]+TMath::Erf((par[2]-x)/20.),
      par[0],
      (2*par[0]*exp( -pow( (par[2]-x)/20. ,2) ))/(20.*sqrt(TMath::Pi()))
      
    };
  
  float er=0.;
  for(int i=0; i<npar;i++)
    { 
      for(int j=0; j<npar;j++)
	{  
	  er =er+derivate[i]*( *(matrix+npar*i+j) )*derivate[j];
	}
    }
  return sqrt(er);
}

double topjetutils::erfMinus_er(double x, double * par, double *matrix)
{   
  const int npar = 3; 
  
  double derivate[npar] = 
    {
      par[1]-TMath::Erf((par[2]-x)/20.),
      par[0],
      (-2*par[0]*exp( -pow( (par[2]-x)/20. ,2) ))/(20.*sqrt(TMath::Pi()))
      
    };
  
  float er=0.;
  for(int i=0; i<npar;i++)
    { 
      for(int j=0; j<npar;j++)
	{  
	  er =er+derivate[i]*( *(matrix+npar*i+j) )*derivate[j];
	}
    }
  return sqrt(er);
}

double topjetutils::exp_func(double x, double *par)
{
  double arg =0;
  if(x!=0 )   arg=par[0]*exp(x*par[1])+par[2];
  else arg=999999;
  return arg;
}

double topjetutils::explinear_func(double x, double *par)
{
  double arg =0;
  if(x!=0 )   arg=par[0]*exp(x*par[1])+par[2]*x+par[3];
  else arg=999999;
  return arg;
}

double topjetutils::erfPlus_func(double x, double *par)
{
  double arg =0;
  if(x!=0 )   arg=par[0]*(par[1]+TMath::Erf((par[2]-x)/20.));

  else arg=999999;
  return arg;
}

double topjetutils::erfMinus_func(double x, double *par)
{
  double arg =0;
  if(x!=0 )   arg=par[0]*(par[1]-TMath::Erf((par[2]-x)/20.));

  else arg=999999;
  return arg;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Total systematic uncertainty functions
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


double topjetutils::JVFSasS_Averaged_Sys_SF( const int pT, const std::string /*year*/, const int up){
	
	double sys_SF = 0;
	
	if(up>0){
	
		double table_up_pT[66] = {	20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,160,170,180,190,200,
						210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,
						440,450,460,470,480,490,500,550,600,650,700,750,800,850,900,950	};
		double table_up_SF[66] = {	1.02294,1.0214,1.01985,1.01848,1.01711,1.0159,1.01469,1.01363,1.01257,1.01164,1.0107,1.00989,
						1.00908,1.00837,1.00766,1.00706,1.00645,1.0054,1.00451,1.00374,1.00307,1.00248,1.00195,1.00147,
						1.00104,1.00065,1.00032,1.00005,0.999831,0.999682,0.999597,0.999577,0.999616,0.99971,0.999853,
						1.00004,1.00027,1.00054,1.00083,1.00116,1.00151,1.00189,1.00229,1.00271,1.00315,1.0036,1.00407,
						1.00456,1.00506,1.00557,1.00609,1.00663,1.00717,1.00772,1.00828,1.00885,1.00942,1.01239,1.01548,
						1.01865,1.02187,1.02513,1.0284,1.03168,1.03496,1.03824	};
		
		for( unsigned int i = 0; i<65; ++i){
			if( pT>=table_up_pT[i] && pT<table_up_pT[i+1]){
				sys_SF = table_up_SF[i];
			}
		}
		if( pT >= table_up_pT[65] ){ return table_up_SF[65]; }	
	}	
	else{
	
		double table_dw_pT[66] = { 	20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,160,170,180,190,200,
						210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,
						440,450,460,470,480,490,500,550,600,650,700,750,800,850,900,950};
		double table_dw_SF[66] = {	0.996849,0.997027,0.997204,0.997342,0.99748,0.997583,0.997686,0.997759,0.997832,0.997878,0.997923,
						0.997946,0.997969,0.997974,0.997978,0.997968,0.997958,0.997916,0.997856,0.99778,0.997687,0.997575,
						0.997433,0.997249,0.997008,0.996694,0.996308,0.99586,0.995369,0.994854,0.994326,0.993796,0.993267,
						0.992742,0.992222,0.991708,0.991201,0.9907,0.990204,0.989713,0.989228,0.988746,0.988269,0.987796,0.987326,
						0.98686,0.986396,0.985935,0.985477,0.985021,0.984567,0.984116,0.983667,0.983219,0.982774,0.98233,0.981889,
						0.979706,0.977564,0.975464,0.973405,0.971388,0.969413,0.967455,0.965497,0.963539};		
		for( unsigned int i = 0; i<65; ++i){
			if( pT>=table_dw_pT[i] && pT<table_dw_pT[i+1]){
				sys_SF = table_dw_SF[i];
			}
		}
		if( pT >= table_dw_pT[65] ){ return table_dw_SF[65]; }	
	}
	
	return sys_SF;
}

double topjetutils::JVFSasS_Herwig_Sys_SF( const int pT, const std::string /*year*/, const int up){
	
	double sys_SF = 0;
	
	if(up>0){
		double table_up_pT[66] = {	20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,160,170,180,190,200,
						210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,
						440,450,460,470,480,490,500,550,600,650,700,750,800,850,900,950	};
		double table_up_SF[66] = {	1.02906,1.02661,1.02416,1.02211,1.02005,1.01839,1.01674,1.01549,1.01425,1.01333,1.01241,1.01166,
						1.0109,1.01022,1.00954,1.00889,1.00825,1.00701,1.00582,1.00469,1.0036,1.00258,1.00163,1.00075,0.999945,
						0.999233,0.998617,0.998099,0.997678,0.997346,0.997091,0.996901,0.996763,0.996668,0.996606,0.996571,
						0.996557,0.996558,0.996573,0.996597,0.996628,0.996665,0.996706,0.99675,0.996795,0.996842,0.996888,
						0.996935,0.99698,0.997025,0.997068,0.99711,0.997151,0.99719,0.997227,0.997262,0.997295,0.997437,0.997541,
						0.997614,0.997665,0.9977,0.997724,0.997743,0.997763,0.997783};
		
		for( unsigned int i = 0; i<65; ++i){
			if( pT>=table_up_pT[i] && pT<table_up_pT[i+1]){
				sys_SF = table_up_SF[i];
			}
		}
		if( pT >= table_up_pT[65] ){ return table_up_SF[65]; }	
	}	
	else{
	
		double table_dw_pT[66] = {	20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,160,170,180,190,200,210,
						220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,
						460,470,480,490,500,550,600,650,700,750,800,850,900,950 };
		double table_dw_SF[66] = {	1.01101,1.01085,1.01069,1.01044,1.01019,1.00985,1.00952,1.00909,1.00865,1.00811,1.00756,1.00689,
						1.00621,1.00545,1.00469,1.00391,1.00313,1.00165,1.00032,0.999146,0.998137,0.997279,0.996561,0.995966,
						0.99548,0.995088,0.994774,0.994519,0.994303,0.994096,0.993859,0.993547,0.993138,0.992649,0.992111,
						0.991552,0.990987,0.990425,0.989873,0.989333,0.988806,0.988293,0.987796,0.987313,0.986846,0.986393,
						0.985955,0.985531,0.98512,0.984724,0.98434,0.98397,0.983611,0.983265,0.98293,0.982606,0.982293,0.980879,
						0.979685,0.978675,0.977821,0.977097,0.976483,0.975908,0.975333,0.974759};

		for( unsigned int i = 0; i<65; ++i){
			if( pT>=table_dw_pT[i] && pT<table_dw_pT[i+1]){
				sys_SF = table_dw_SF[i];
			}
		}
		if( pT >= table_dw_pT[65] ){ return table_dw_SF[65]; }	
	}
	
	return sys_SF;
}


double topjetutils::JVFSasS_Pythia_Sys_SF( const int /*pT*/, const std::string /*year*/, const int up){
	
	double sys_SF = 0;
	
	if(up>0){
		return 0.999219;
	}	
	else{
		return 0.995913;
	}
	return sys_SF;
}




double topjetutils::JVFSasB_Averaged_Sys_SF( const int pT, const std::string /*year*/, const int up){
	
	double sys_SF = 0;
	
	if(up>0){
	
		double table_up_pT[66] = {	20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,
						150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,
						320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,
						490,500,550,600,650,700,750,800,850,900,950 	};
		double table_up_SF[66] = {	1.0491,1.04178,1.03474,1.0378,1.04086,1.04367,1.04642,1.04844,1.05041,
						1.05201,1.05358,1.05496,1.05632,1.05754,1.05873,1.05978,1.0608,1.06248,
						1.06379,1.06475,1.06538,1.06576,1.06591,1.06588,1.06572,1.06546,1.06512,
						1.06472,1.06427,1.0638,1.06331,1.06281,1.0623,1.06178,1.06127,1.06075,
						1.06024,1.05974,1.05924,1.05874,1.05825,1.05777,1.05729,1.05682,1.05636,
						1.05591,1.05546,1.05502,1.05459,1.05416,1.05375,1.05334,1.05293,1.05254,
						1.05215,1.05177,1.0514,1.04965,1.04808,1.04671,1.04551,1.04451,1.04368,
						1.04293,1.04218,1.04143 	};
		for( unsigned int i = 0; i<65; ++i){
			if( pT>=table_up_pT[i] && pT<table_up_pT[i+1]){
				sys_SF = table_up_SF[i];
			}
		}
		if( pT >= table_up_pT[65] ){ return table_up_SF[65]; }
			
	}	
	else{
	
		double table_dw_pT[66] = { 	20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,
						150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,
						320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,
						490,500,550,600,650,700,750,800,850,900,950 };
		double table_dw_SF[66] = {	0.830183,0.797842,0.76556,0.735557,0.705945,0.683683,0.661868,0.645498,
						0.629544,0.617296,0.605422,0.596122,0.587154,0.580038,0.573217,0.567763,
						0.562569,0.554451,0.548287,0.543629,0.540126,0.537503,0.535547,0.534094,
						0.533018,0.532224,0.531638,0.531208,0.530892,0.530661,0.530491,0.530368,
						0.530277,0.530211,0.530163,0.530128,0.530103,0.530084,0.530071,0.530061,
						0.530054,0.530049,0.530045,0.530042,0.53004,0.530039,0.530038,0.530037,
						0.530037,0.530036,0.530036,0.530036,0.530036,0.530036,0.530035,0.530035,
						0.530035,0.530035,0.530035,0.530035,0.530035,0.530035,0.530035,0.530035,
						0.530035,0.530035 };

		for( unsigned int i = 0; i<65; ++i){
			if( pT>=table_dw_pT[i] && pT<table_dw_pT[i+1]){
				sys_SF = table_dw_SF[i];
			}
		}
		if( pT >= table_dw_pT[65] ){ return table_dw_SF[65]; }	
	}
	
	return sys_SF;
}

double topjetutils::JVFSasB_Herwig_Sys_SF( const int pT, const std::string /*year*/, const int up){
	
	double sys_SF = 0;
	
	if(up>0){
	
		double table_up_pT[66] = {	20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,
						150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,
						320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,
						490,500,550,600,650,700,750,800,850,900,950};
		double table_up_SF[66] = {	0.928098,0.873497,0.819773,0.799832,0.780095,0.764192,0.748421,0.734253,
						0.72033,0.708594,0.697217,0.688315,0.679832,0.673684,0.66795,0.664135,
						0.660678,0.65686,0.655405,0.655446,0.656363,0.657741,0.659319,0.660939,
						0.66251,0.663983,0.665337,0.666564,0.667666,0.66865,0.669526,0.670303,
						0.670992,0.671601,0.67214,0.672616,0.673037,0.673409,0.673738,0.674028,
						0.674284,0.67451,0.67471,0.674887,0.675043,0.675181,0.675302,0.675409,
						0.675504,0.675588,0.675662,0.675727,0.675784,0.675835,0.67588,0.675919,
						0.675954,0.676075,0.67614,0.676174,0.676193,0.676202,0.676208,0.676211,
						0.676215,0.676219};
						
		for( unsigned int i = 0; i<65; ++i){
			if( pT>=table_up_pT[i] && pT<table_up_pT[i+1]){
				sys_SF = table_up_SF[i];
			}
		}
		if( pT >= table_up_pT[65] ){ return table_up_SF[65]; }	
	}	
	else{
		double table_dw_pT[66] = {	20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,
						160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,
						340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490,500,550,
						600,650,700,750,800,850,900,950};
		double table_dw_SF[66] = {	0.759259,0.73971,0.720001,0.694132,0.668465,0.646586,0.625123,0.608726,
						0.592627,0.579168,0.565885,0.553813,0.541928,0.531087,0.520499,0.511084,
						0.501971,0.486432,0.473686,0.463396,0.455186,0.448697,0.443606,0.439635,
						0.436555,0.434176,0.432344,0.43094,0.429866,0.429046,0.428423,0.427949,
						0.42759,0.427319,0.427113,0.426958,0.426842,0.426754,0.426688,0.426638,
						0.426601,0.426573,0.426552,0.426537,0.426525,0.426516,0.42651,0.426505,
						0.426501,0.426499,0.426497,0.426495,0.426494,0.426493,0.426492,0.426492,
						0.426492,0.426491,0.426491,0.426491,0.426491,0.426491,0.426491,0.426491,
						0.426491,0.426491};

		for( unsigned int i = 0; i<65; ++i){
			if( pT>=table_dw_pT[i] && pT<table_dw_pT[i+1]){
				sys_SF = table_dw_SF[i];
			}
		}
		if( pT >= table_dw_pT[65] ){ return table_dw_SF[65]; }	
	}
	
	return sys_SF;
}


double topjetutils::JVFSasB_Pythia_Sys_SF( const int /*pT*/, const std::string /*year*/, const int up){

	if(up>0){
		return 1.09865;
	}	
	else{
		return 0.986779;
	}

}


double topjetutils::JVFSasS_SF_err(double pT, int up, const std::string year, const int ps){

	if(year=="2011" || year=="2012"){
		double err = TMath::Sqrt(TMath::Power(TMath::Abs( JVFSasS_SF_errSelection(pT, up, year) ),2)+TMath::Power(TMath::Abs( JVFSasS_SF_errFit(pT, year) ),2));	

		if (up>0)
			return JVFSasS_SF(pT, year, ps) + err;
		else
			return JVFSasS_SF(pT, year, ps) - err;
	}
	else{
		if(ps==0){//Pythia specific systematics
			return JVFSasS_Pythia_Sys_SF(pT,year,up);
		}
		else if(ps==1){//Herwig specific systematics
			return JVFSasS_Herwig_Sys_SF(pT,year,up);
		}
		else{
			return JVFSasS_Averaged_Sys_SF(pT,year,up);
		}
	}
				
}

double topjetutils::JVFSasB_SF_err(double pT, int up, const std::string year, const int ps){


	if( year=="2011" || year=="2012"){
		double err = TMath::Sqrt(TMath::Power(TMath::Abs( JVFSasB_SF_errSelection(pT, up, year) ),2)+TMath::Power(TMath::Abs( JVFSasB_SF_errFit(pT, year) ),2));;

		if (up>0)
			return JVFSasB_SF(pT, year, ps) + err;
		else
			return JVFSasB_SF(pT, year, ps) - err;
	}
	else{
		if(ps==0){//Pythia specific systematics
			return JVFSasB_Pythia_Sys_SF(pT,year,up);
		}
		else if(ps==1){//Herwig specific systematics
			return JVFSasB_Herwig_Sys_SF(pT,year,up);
		}
		else{
			return JVFSasB_Averaged_Sys_SF(pT,year,up);	
		}
	}
}

double topjetutils::JVFBasB_SF_err(double pT, int up, const std::string year, const int ps){

	if(year=="2011" || year=="2012"){
	
		double err = 2*(TMath::Sqrt(TMath::Power(TMath::Abs(JVFSasS_SF_errSelection(pT, up, year) ),2)+TMath::Power(TMath::Abs( JVFSasS_SF_errFit(pT, year) ),2)));

		if (up>0)
			return 1. + err;
		else
			return 1. - err;
	}
	else{
		//double error = 0;
		if(ps==0){//Pythia specific systematics
			if(up>0){
				return (1.05);			
			}
			else{
				return (0.95);				
			}
		}
		else if(ps==1){//Herwig specific systematics
			if(up>0){
				return (1.05);			
			}
			else{
				return (0.95);				
			}
		}
		else{
			if(up>0){
				return (1.05);			
			}
			else{
				return (0.95);				
			}
		}		
	}
}	

double topjetutils::JVFBasS_SF_err(double pT, int up, const std::string year, const int ps ){


	if(year=="2011" || year=="2012"){

		double err = 2*(TMath::Sqrt(TMath::Power(TMath::Abs( JVFSasB_SF_errSelection(pT, up, year) ),2)+TMath::Power(TMath::Abs( JVFSasB_SF_errFit(pT, year) ),2)));

		if (up>0)
			return 1. + err;
		else
			return 1. - err;
	
	}
	else{
		double error = 0;
		if(ps==0){//Pythia specific systematics
			error = TMath::Abs(JVFSasB_Pythia_Sys_SF(pT,year,up)-JVFSasB_SF(pT, year, ps));
			if(up>0){
				return (1+2*error);			
			}
			else{
				return (1-2*error);				
			}
		}
		else if(ps==1){//Herwig specific systematics
			error = TMath::Abs(JVFSasB_Herwig_Sys_SF(pT,year,up)-JVFSasB_SF(pT, year, ps));
			if(up>0){
				return (1+2*error);			
			}
			else{
				return (1-2*error);				
			}
		}
		else{	
			error = TMath::Abs(JVFSasB_Averaged_Sys_SF(pT,year,up)-JVFSasB_SF(pT, year, ps));
			if(up>0){
				return (1+1*error);			
			}
			else{
				return (1-1*error);				
			}
		}	
	}
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Fit systematic uncertainty functions
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double topjetutils::JVFSasS_SF_errFit(double pT, std::string year){

	double err=-1.;
	
	if (year=="2011"){
		double par_fit[3] = {0.0298319, -0.0119208, 1.00479};
	
		double chi2 = 0.965373;

		double matrix_fit[3][3] = 	{{1.6416e-05,-2.20442e-07,-6.91319e-06},
                        	   		{-2.20442e-07,3.93417e-05,-3.10102e-05},
                        	   		{-6.91319e-06,-3.10102e-05,2.80743e-05}
				  		};
	
		if(chi2>1.) err = TMath::Sqrt(chi2)*exp_er(pT, par_fit, &matrix_fit[0][0]);
		else err = exp_er(pT, par_fit, &matrix_fit[0][0]);
	}else{
		double par_fit[3] = {0.0382635, -0.0121108, 0.998046};
	
		double chi2 = 2.18524;

		double matrix_fit[3][3] = 	{{1.08076e-05,5.59519e-07,-5.42049e-06},
                        			{5.59519e-07,8.02464e-06,-8.56388e-06},
                        			{-5.42049e-06,-8.56388e-06,1.17015e-05}
						};

	
		if(chi2>1.) err = TMath::Sqrt(chi2)*exp_er(pT, par_fit, &matrix_fit[0][0]);
		else err = exp_er(pT, par_fit, &matrix_fit[0][0]);
	}	
							
	return err;
}

double topjetutils::JVFSasB_SF_errFit(double pT, std::string year){

	double err=1.;

	if (year=="2011"){
		double par_fit[3] = {0.436658, -0.0167345, 0.5376};
	
		double chi2 = 0.983234;

		double matrix_fit[3][3] = 	{{0.0039227,0.000170431,-0.00335489},
                        			{0.000170431,8.71661e-05,-0.000815901},
                        			{-0.00335489,-0.000815901,0.00858169}
						};

	
		if(chi2>1.) err = TMath::Sqrt(chi2)*exp_er(pT, par_fit, &matrix_fit[0][0]);
		else err = exp_er(pT, par_fit, &matrix_fit[0][0]);
	}else{
		double par_fit[4] = {0.961611, -0.0389405, 0.00281162, 0.362232};

		double chi2 = 0.961223;
		
		double matrix_fit[4][4] = 	{{0.0378161,-0.00148858,1.13215e-05,0.00314679},
                        			{-0.00148858,0.000245578,1.72816e-05,-0.0023322},
                        			{1.13215e-05,1.72816e-05,1.87664e-06,-0.000217394},
                        			{0.00314679,-0.0023322,-0.000217394,0.0268588}
						};
						
		if(chi2>1.) err = TMath::Sqrt(chi2)*explinear_er(pT, par_fit, &matrix_fit[0][0]);
		else err = explinear_er(pT, par_fit, &matrix_fit[0][0]);				
	
	}
									
	return err;
}

double topjetutils::JVFBasB_SF_errFit(double pT, std::string /*year*/){

	double par_fit[3] = {0.188125, 4.24933, 48.7055};
	
	double chi2 = 3.84489;

	double matrix_fit[3][3] = 	{{0.000968296,-0.0274005,0.0674636},
                        		{-0.0274005,0.776175,-1.94882},
                        		{0.0674636,-1.94882,7.46003}
					};
	
	double err;
	if(chi2>1.) err = TMath::Sqrt(chi2)*erfPlus_er(pT, par_fit, &matrix_fit[0][0]);
	else err = erfPlus_er(pT, par_fit, &matrix_fit[0][0]);
								
	return err;
}

double topjetutils::JVFBasS_SF_errFit(double pT, std::string /*year*/){

	double par_fit[3] = {0.228856, 5.41818, 41.871};
	
	double chi2 = 2.66638;

	double matrix_fit[3][3] = 	{{0.00172011,-0.0324136,0.100606},
                        		{-0.0324136,0.616389,-1.71014},
                        		{0.100606,-1.71014,12.9077}
					};

	double err;
	if(chi2>1.) err = TMath::Sqrt(chi2)*erfMinus_er(pT, par_fit, &matrix_fit[0][0]);
	else err = erfMinus_er(pT, par_fit, &matrix_fit[0][0]);
							
	return err;
	
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Selection systematic uncertainty functions
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double topjetutils::JVFSasS_SF_errSelection(double pT, int up, std::string year){

	double err = -1.;
		
	if (year=="2011"){
		double par_nominal[3] = {0.0298319, -0.0119208, 1.00479};
		double par_deltaphiUp[3] = {0.030911, -0.0107551, 1.00357};
		double par_deltaphiDw[3] = {0.0372879, -0.0218913, 1.00655};
		double par_pTUp[3] = {0.0349472, -0.0168269, 1.00688};
		double par_pTDw[3] = {0.0264818, -0.00972147, 1.00381};
	
						
		if (up>0)									
			err = TMath::Sqrt(TMath::Power(TMath::Abs(exp_func(pT,par_deltaphiUp)-exp_func(pT,par_nominal)),2)+TMath::Power(TMath::Abs(exp_func(pT,par_pTUp)-exp_func(pT,par_nominal)),2));
		else
			err = TMath::Sqrt(TMath::Power(TMath::Abs(exp_func(pT,par_deltaphiDw)-exp_func(pT,par_nominal)),2)+TMath::Power(TMath::Abs(exp_func(pT,par_pTDw)-exp_func(pT,par_nominal)),2));
	}else{
		double par_nominal[3] = { 0.0376064, -0.00983575, 0.995614,  };
		double par_deltaphiUp[3] = { 0.0387806, -0.012209, 0.997475,  };
		double par_deltaphiDw[3] = { 0.0456764, -0.00826356, 0.987249,  };
		double par_pTUp[3] = { 0.0470831, -0.0147344, 0.998215,  };
		double par_pTDw[3] = { 0.0366557, -0.0104705, 0.996454,  };
		
		if (up>0)									
			err = TMath::Sqrt(TMath::Power(TMath::Abs(exp_func(pT,par_deltaphiUp)-exp_func(pT,par_nominal)),2)+TMath::Power(TMath::Abs(exp_func(pT,par_pTUp)-exp_func(pT,par_nominal)),2));
		else
			err = TMath::Sqrt(TMath::Power(TMath::Abs(exp_func(pT,par_deltaphiDw)-exp_func(pT,par_nominal)),2)+TMath::Power(TMath::Abs(exp_func(pT,par_pTDw)-exp_func(pT,par_nominal)),2));
	
	}
		
	return	err;	
}

double topjetutils::JVFSasB_SF_errSelection(double pT, int up, std::string year){

	double err = -1.;	
	
	if (year=="2011"){
		double par_nominal[3] = {0.436658, -0.0167345, 0.5376};
		double par_deltaphiUp[3] = {0.303512, -0.0151723, 0.6};
		double par_deltaphiDw[3] = {0.454962, -0.0201396, 0.561454};
		double par_pTUp[3] = {0.282541, -0.0188132, 0.6};
		double par_pTDw[3] = {0.522894, -0.0258739, 0.581227};
	
						
		if (up>0)									
			err = TMath::Sqrt(TMath::Power(TMath::Abs(exp_func(pT,par_deltaphiUp)-exp_func(pT,par_nominal)),2)+TMath::Power(TMath::Abs(exp_func(pT,par_pTUp)-exp_func(pT,par_nominal)),2));
		else
			err = TMath::Sqrt(TMath::Power(TMath::Abs(exp_func(pT,par_deltaphiDw)-exp_func(pT,par_nominal)),2)+TMath::Power(TMath::Abs(exp_func(pT,par_pTDw)-exp_func(pT,par_nominal)),2));
	}else{
		double par_nominal[4] = { 1.05915, -0.0407948, 0.00288741, 0.357726,  };
		double par_deltaphiUp[4] = { 1.03204, -0.0421817, 0.00344174, 0.3488,  };
		double par_deltaphiDw[4] = { 1.01214, -0.048076, 0.00228067, 0.451086,  };
		double par_pTUp[4] = { 1.12715, -0.0543705, 0.00259316, 0.408258,  };
		double par_pTDw[4] = { 0.972341, -0.0329967, 0.00316192, 0.306292,  };
		
		if (up>0)									
			err = TMath::Sqrt(TMath::Power(TMath::Abs(explinear_func(pT,par_deltaphiUp)-explinear_func(pT,par_nominal)),2)+TMath::Power(TMath::Abs(explinear_func(pT,par_pTUp)-explinear_func(pT,par_nominal)),2));
		else
			err = TMath::Sqrt(TMath::Power(TMath::Abs(explinear_func(pT,par_deltaphiDw)-explinear_func(pT,par_nominal)),2)+TMath::Power(TMath::Abs(explinear_func(pT,par_pTDw)-explinear_func(pT,par_nominal)),2));
	
	}
		
	return	err;	
}

double topjetutils::JVFBasB_SF_errSelection(double pT, int up, std::string /*year*/){

	double par_nominal[3] = {0.188125, 4.24933, 48.7055};
	double par_pTUp[3] = {0.117526, 7.64243, 52.2137};
	double par_pTDw[3] = {0.205827, 3.71907, 45.6386}; 
	
	double err = -1.;						
	if (up>0)									
		err = TMath::Abs(erfPlus_func(pT,par_pTUp)-erfPlus_func(pT,par_nominal));
	else
		err = TMath::Abs(erfPlus_func(pT,par_pTDw)-erfPlus_func(pT,par_nominal));
		
	return	err;	
}

double topjetutils::JVFBasS_SF_errSelection(double pT, int up, std::string /*year*/){

	double par_nominal[3] = {0.228856, 5.41818, 41.871};
	double par_pTUp[3] = {0.199202, 6.22488, 39.3176};
	double par_pTDw[3] = {0.0966727, 10.9839, 44.0509};
	
	double err = -1.;						
	if (up>0)									
		err = TMath::Abs(erfMinus_func(pT,par_pTUp)-erfMinus_func(pT,par_nominal));
	else
		err = TMath::Abs(erfMinus_func(pT,par_pTDw)-erfMinus_func(pT,par_nominal));
		
	return	err;	
}

