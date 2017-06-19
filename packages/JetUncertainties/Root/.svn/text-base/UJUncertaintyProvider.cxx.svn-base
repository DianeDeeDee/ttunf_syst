#include "JetUncertainties/UJUncertaintyProvider.h"
#include "TEnv.h"
#include <stdlib.h>

/////////////////////
// Constructor/destructor
/////////////////////

// Constructor
//compatibility with MultijetUJ: LC needs to be both here and in the config file
//the function parameter is ignored later on in the UJUncertaintyProvider
UJUncertaintyProvider::UJUncertaintyProvider(TString inputFile, TString jetAlgo, TString MCtype):
  m_file(0), m_GeV(1000.0), m_isInit(false), m_NPVRef(11.8), m_muRef(20.7)
{
  initUJ(inputFile, jetAlgo, MCtype);
}

//for persistification, doesn't call initUJ() automatically
UJUncertaintyProvider::UJUncertaintyProvider():
  m_file(0), m_GeV(1000.0), m_isInit(false), m_NPVRef(11.8), m_muRef(20.7) { }

// Destructor
UJUncertaintyProvider::~UJUncertaintyProvider()
{
  if (m_file) {
    m_file->Close();
    delete m_file;
  }
}

int UJUncertaintyProvider::getComponentIndex(TString compName) {
  if (m_compIndex.find(compName)==m_compIndex.end())
    Fatal("UJUncertaintyProvider::getComponentIndex()",
	  "ERROR: No such uncertainty component: %s",compName.Data());
  return m_compIndex[compName];
}

void UJUncertaintyProvider::warning(TString m1, TString m2){
  static int nWarn_b = 0;
  if (m_permissive){
    if (++nWarn_b<m_warning) 
      Warning(m1,m2+TString::Format("\n (only first %i messages will be printed)",m_warning));
  } else {
    //Fatal(m1,m2+TString("\n (use UJUncertaintyProvider.usePermissive(true); to ignore this error)"));
    if (++nWarn_b<m_warning) 
      Warning(m1,m2
	      +TString::Format("\n (%f will be returned, cf. TWiki on how to change this)",m_unphysical)
	      +TString::Format("\n (only first %i messages will be printed)",m_warning));
    m_return_unphysical=true;
  }
}

double UJUncertaintyProvider::boundaries(TAxis *axis, double value, TString variable, bool checkRegion){
  //ranges are defined inside the histogram by the left edge of the first bin and the right edge of the last bin
  //if not inside the range, a warning is prompt and the value at the edge is returned
  //if between the edge and the center of the bin, the interpolation function is not able to compute something
  //so, in this case, the value at bin center is returned

  static int nWarn_b_cr = 0;

  double new_value = value;
  double epsilon = 0.001;

  if ((m_validarea == 0) or (checkRegion)){
    if ((axis->GetBinLowEdge(1) > value) or (value  > axis->GetBinUpEdge(axis->GetNbins()))){
      warning("UJUncertaintyProvider::boundaries()",TString::Format("jet %s=%.3f outside range.",(const char*)variable,value));
      if (checkRegion){
	if (++nWarn_b_cr<m_warning) 
	  Warning("UJUncertaintyProvider::boundaries()",
		  "component outside range, you should not use this component in this region"
		  +TString::Format("\n (only first %i messages will be printed)",m_warning));
      }
    }
  }
  
  if (axis->GetBinCenter(1) >= value)
    new_value = (axis->GetBinCenter(1))+(epsilon*0.5*axis->GetBinWidth(1));
  else if (value >= axis->GetBinCenter(axis->GetNbins()) )
    new_value = (axis->GetBinCenter(axis->GetNbins()))-(epsilon*0.5*axis->GetBinWidth(axis->GetNbins()));

  return new_value; 
}

double UJUncertaintyProvider::readPtEtaHisto(TH2 *h2d, double pt, double eta, bool checkRegion) {
  //ranges are defined inside the TH2 by the left edge of the first bin and the right edge of the last bin
  //if not inside the range, a warning is prompt and the value at the edge is returned

  m_return_unphysical = false;
  if (h2d == 0){
    Warning("UJUncertaintyProvider::readPtEtaHisto","histogram does not exist");
    return m_unphysical;
  }

  pt = boundaries(h2d->GetXaxis(), pt/m_GeV, "pt", checkRegion)*m_GeV;
  eta = boundaries(h2d->GetYaxis(), fabs(eta), "abs(eta)", checkRegion);
  if (m_return_unphysical)
    return m_unphysical;
  else
    return h2d->Interpolate(pt/m_GeV, fabs(eta));
}

double UJUncertaintyProvider::readPtEtaHisto(TH3 *h3d, double pt, double mass, double eta, bool checkRegion) {
  double mopt = mass / pt;

  m_return_unphysical = false;
  if (h3d == 0){
    Warning("UJUncertaintyProvider::readPtEtaHisto","histogram does not exist");
    return m_unphysical;
  }

  //use the validarea histo to check if pt,mass,eta are in
  //the region where the method is valid
  if (m_validarea){
    int ix = m_validarea->GetXaxis()->FindBin(pt/m_GeV);
    int iy = m_validarea->GetYaxis()->FindBin(mopt);
    int iz = m_validarea->GetZaxis()->FindBin(fabs(eta));
    if (m_validarea->GetBinContent(ix,iy,iz)==0)
      warning("UJUncertaintyProvider::readPtEtaHisto()",TString::Format("jet outside range: (pT=%.1f GeV, eta=%.2f, mass=%.1f).",pt/m_GeV,eta,mass/m_GeV));
  }

  //retro compatibility
  //in order to have the same behavior 
  if (m_old){
    if ( pt/m_GeV<100.0 || pt/m_GeV>=1000.0 ) {
      warning("UJUncertaintyProvider::readPtEtaHisto()",TString::Format("jet pT outside 100-1000 range: pT=%.1f GeV.",pt/m_GeV));
      if (pt/m_GeV<100.0) pt=100.0*m_GeV;
      if (pt/m_GeV>=1000.0) pt=1000.0*m_GeV;
    }
    if (fabs(eta)>=1.2) {
      warning("UJUncertaintyProvider::readPtEtaHisto()",TString::Format("jet eta=%.2f outside |eta|<1.2.",eta));
      eta=1.2;
    }
    if (mass/m_GeV<0.0 || mass/m_GeV>=400.0) {
      warning("UJUncertaintyProvider::readPtEtaHisto()",TString::Format("jet mass=%.2f outside 0-400 range.",mass/m_GeV));
      if (mass/m_GeV<0.0) { mass = 0.0*m_GeV; }
      if (mass/m_GeV>400.0) { mass = 400.0*m_GeV; }
    }
    mopt = mass / pt;
  }
  else{
    pt = boundaries(h3d->GetXaxis(), pt/m_GeV, "pt", checkRegion)*m_GeV;
    mopt = boundaries(h3d->GetYaxis(), mopt, "mass/pt", checkRegion);
    eta = boundaries(h3d->GetZaxis(), fabs(eta), "abs(eta)", checkRegion);
  }

  double value = h3d->Interpolate(pt/m_GeV, mopt,fabs(eta));
  if (m_return_unphysical)
    return m_unphysical;
  else
    return value;
}

/////////////////////
// Initialisation
/////////////////////
void UJUncertaintyProvider::initUJ(TString inputFile, TString jetAlgo, TString MCtype) 
{

  //prevent multiple initializations
  if (m_isInit == false) {

    m_warning = 100;
    m_permissive = false;
    m_return_unphysical = false;//used if m_permissive = false
    //m_unphysical = std::numeric_limits<double>::quiet_NaN();
    m_unphysical = -1E33;//quiet_NaN can crash other scripts

    m_file = new TFile(inputFile, "READ");
    if (!m_file->IsOpen()) {
      Fatal("UJUncertaintyProvider::initUJ()", "Cannot open input file %s!", (const char*)inputFile);
    }

    //retro compatibility
    //the previous behavior with the previous rootfile is not changed
    m_old = false;
    if (inputFile.Contains("UJInput2012.root")) m_old = true;

    m_jetAlgo = jetAlgo;
    m_MCtype = MCtype;

    // Load the histogram defining the valid area
    // If this histogram does not exist (m_validarea==0),
    //  the borders of the TH3F are used
    // The values of "m_compCheckRegion" supersede Valid_area for components
    m_validarea = (TH3C*)m_file->Get("Valid_area");

    // Load the uncertainty shift
    // Shifts: hard-coded
    const int _SHIFT_SIZE = 2;
    static const char* shiftnames[_SHIFT_SIZE] = {"Pileup_NPV","Pileup_Mu"};
    for (int i = 0; i < _SHIFT_SIZE; ++i){
      TString hkey;
      hkey = TString::Format("Shift_%s_%s", (const char*)jetAlgo, shiftnames[i]);
      TH1F* unch_shift = (TH1F*)m_file->Get(hkey);
      m_shift_unc.push_back(unch_shift);
    }


    // Always check this matches e_prop in the headers if you change it
    static const char* propstr[_PROP_SIZE] = {"pT", "m", "tau21", "tau32", "d12", "d23", "E"};

    // Components: hard-coded
    const int _COMP_SIZE = 18;

    static const char* compnames[_COMP_SIZE] = {"Total",
						"Component_dataMC",
						"Component_pt2Cut",
						"Component_dPhiCut",
						"Component_photonPurity",
						"Component_PES",
						"Component_generator",
						"Component_kterm",
						"Component_JER",
						"Component_akt4insideOutsideLargeR",
						"Component_more1smallJetInsideLargeR",
						"Component_stats",
						"Component_MoverPt",
						"Component_topology",// used in compCheckRegion, don't move it
						"Component_DoubleRatioInterpolation",
						"Component_DoubleRatioInterpolation700to900",
						"Component_DoubleRatioInterpolation800to1000",
						"Component_DoubleRatioInterpolation900to1100"
    };
    static const char* compdescr[_COMP_SIZE] = {"Total Uncertainty",
						"dataMC (Gamma+jet method)",
						"pt2Cut (Gamma+jet method)",
						"dPhiCut (Gamma+jet method)",
						"photonPurity (Gamma+jet method)",
						"PES (Gamma+jet method)",
						"generator (Gamma+jet method)",
						"kterm (Gamma+jet method)",
						"JER (Gamma+jet method)",
						"akt4insideOutsideLargeR (Gamma+jet method)",
						"more1smallJetInsideLargeR (Gamma+jet method)",
						"stats (Gamma+jet method)",
						"M/pt (<0.15 vs. >0.15)  (Gamma+jet method)",
						"topology (top content) (MC method)",
						"Interpolation to double-ratio method (for Gamma+jet method)",
						"Interpolation to double-ratio method (for Gamma+jet method, 700 -> 900)",
						"Interpolation to double-ratio method (for Gamma+jet method, 800 -> 1000)",
						"Interpolation to double-ratio method (for Gamma+jet method, 900 -> 1100)"};
    //m_compIndex["Total Uncertainty"] = 0;

    for (int j = 0; j < _COMP_SIZE; j++) {

      std::vector<TH3*> temp_m_unc;
      bool temp_found=false;

      for (int i = 0; i < _PROP_SIZE; i++) {
	TString hkey;
	if (j == 0){
	  hkey = TString::Format("%s_%s", (const char*)jetAlgo, propstr[i]);
	}
	else{
	  hkey = TString::Format("%s_%s_%s", compnames[j], (const char*)jetAlgo, propstr[i]);
	}
	TH3F* unch = (TH3F*)m_file->Get(hkey);
	//if (!unch) {
	//  Warning("UJUncertaintyProvider::initUJ()", "Algorithm %s and variable %s is not present in input file!", (const char*)jetAlgo, propstr[i]);
	//  Warning("UJUncertaintyProvider::initUJ()", "It will crash if someone try to use it!");
	//}
	temp_m_unc.push_back(unch);
	if (unch) temp_found=true;	
      }

      if (temp_found){
	for (unsigned int i=0; i< temp_m_unc.size();++i){
	  m_unc.push_back(temp_m_unc.at(i));
	  if (temp_m_unc.at(i))
	    Info("UJUncertaintyProvider::init()","  Store histogram (%s) in position %i",(const char*)temp_m_unc.at(i)->GetName(),(int) m_unc.size()-1);
	  else
	    Info("UJUncertaintyProvider::init()","  Store histogram (empty) in position %i",(int) m_unc.size()-1);
	}
	m_compIndex[compnames[j]]=m_compNames.size();
	m_compNames.push_back(compnames[j]);
	m_compDesc.push_back(compdescr[j]);
	//fill m_compCheckRegion
	//if you look at this component, you have to limit yourself to this component region
	if (j == 13){//13 is: Component_topology
	  m_compCheckRegion.push_back(true);
	} else {
	  m_compCheckRegion.push_back(false);
	}
      }
    }

    Info("UJUncertaintyProvider::init()","================================================");
    Info("UJUncertaintyProvider::init()","  Initialized the UJUncertaintyProvider tool");
    //Info("UJUncertaintyProvider::init()","  Using track-jet uncertainties" );
    Info("UJUncertaintyProvider::init()","  For %s",(const char*)jetAlgo);
    
  }
  else {
    Warning( "UJUncertaintyProvider::init()", "WARNING: UJUncertaintyProvider already initialized, skipping re-initialization");
  }
}

//bool UJUncertaintyProvider::setInputCollection(TString jetAlgo) {} 

TString UJUncertaintyProvider::getMCtype() {
  if (m_MCtype.Contains("MC11b", TString::kIgnoreCase)) return "MC11b";
  if (m_MCtype.Contains("MC11c", TString::kIgnoreCase)) return "MC11c";
  if (m_MCtype.Contains("FrozenShowers", TString::kIgnoreCase)) return "FrozenShowers";
  if (m_MCtype.Contains("AFII", TString::kIgnoreCase)) return "AFII";
  return "Error";
}

double UJUncertaintyProvider::getRelUncertComponent(TString compName, double pT, double mass, double eta, e_prop p) {
  return getRelUncertComponent(getComponentIndex(compName),pT,mass,eta,p);
}

double UJUncertaintyProvider::getRelUncertComponent(int iComp, double pT, double mass, double eta, e_prop p) {
  if (iComp>getNUncertaintyComponents())
    Fatal("UJUncertaintyProvider::getRelUncertComponent","UJ uncertainty component index out-of-range.");

  TH3 *h3d = m_unc.at(p+(_PROP_SIZE*iComp));

  if (m_compCheckRegion.at(iComp))
    return readPtEtaHisto(h3d,pT,mass,eta,true);
  else
    return readPtEtaHisto(h3d,pT,mass,eta);
}

// This method should be removed at some point
double UJUncertaintyProvider::getRelUncertComponent(TString compName, double pT, double mass, double eta, 
						     double NPV, double mu, e_prop p) {
  return getRelUncertComponent(getComponentIndex(compName),pT,mass,eta,NPV,mu, p);
}

// This method should be removed at some point
double UJUncertaintyProvider::getRelUncertComponent(int iComp, double pT, double mass, double eta, 
						     double /*NPV*/, double /*mu*/, e_prop p) {

  return getRelUncertComponent(iComp,pT,mass,eta,p);
}


/////////////////////
// Uncertainty 
/////////////////////

// Absolute Uncertainty
double UJUncertaintyProvider::getAbsUncert(double pT, double mass, double eta, double NPV, double mu, double val, e_prop p) {
  return getRelUncert(pT, mass, eta, NPV, mu, p)*val;
}

// Shortcut for GammaJet method
double UJUncertaintyProvider::getRelUncert_GammaJet_Apr042014(double pT, double mass, double eta){

  static const char* components[] = {"Component_dataMC",
				     "Component_pt2Cut",
				     "Component_dPhiCut",
				     "Component_photonPurity",
				     "Component_PES",
				     "Component_generator",
				     "Component_kterm",
				     "Component_JER",
				     "Component_akt4insideOutsideLargeR",
				     "Component_more1smallJetInsideLargeR",
				     "Component_stats",
				     "Component_MoverPt"};
  double total_square = 0;
  for (int i=0; i<12; ++i){
    total_square += pow(getRelUncertComponent(components[i], pT, mass, eta, UJUncertaintyProvider::JPTS),2);
  }

  return sqrt(total_square);
}


// Shortcut for GammaJet method                                                                                                                                 
double UJUncertaintyProvider::getRelUncert_GammaJet_ExcludingPileup(double pT, double mass, double eta){
  e_prop p = UJUncertaintyProvider::JPTS;
  return fabs(getRelUncertComponent("Total", pT, mass, eta, p)) + fabs(getRelUncertBias(pT,mass,eta, p));
}

// Shortcut for GammaJet method                                                                                                                                 
double UJUncertaintyProvider::getRelUncert_GammaJet_ShiftNPV(double pT, double NPV){
  double value=0;
  for (unsigned int i=0;i < m_shift_unc.size(); ++i){
    if (m_shift_unc[i]!=0){
      if (TString(m_shift_unc[i]->GetName()).Contains("Pileup_NPV")){
	value = m_shift_unc[i]->GetBinContent((m_shift_unc[i]->FindBin(pT/m_GeV)));
	break;
      }
    }
  }
  return value * (NPV - m_NPVRef);
}
double UJUncertaintyProvider::getRelUncert_GammaJet_ShiftMu(double pT, double mu){
  double value=0;
  for (unsigned int i=0;i < m_shift_unc.size(); ++i){
    if (m_shift_unc[i]!=0){
      if (TString(m_shift_unc[i]->GetName()).Contains("Pileup_Mu")){
	value = m_shift_unc[i]->GetBinContent((m_shift_unc[i]->FindBin(pT/m_GeV)));
	break;
      }
    }
  }
  return value * (mu - m_muRef);
}

// Shortcut for DoubleRatio method                                                                                                                              
double UJUncertaintyProvider::getRelUncert_DoubleRatio(double pT, double mass, double eta, e_prop p){
  return fabs(getRelUncertComponent("Total", pT, mass, eta, p)) + fabs(getRelUncertBias(pT,mass,eta, p));
}


double UJUncertaintyProvider::getRelUncert(double pT, double mass, double eta, double NPV, double mu, e_prop p) {
  return sqrt( pow(getRelUncert(pT,mass,eta,p),2) + pow(getRelOffsetUncert(pT,mass,eta,NPV,mu,p),2) );
}

double UJUncertaintyProvider::getRelUncert(double pT, double mass, double eta, e_prop p) {
  return fabs(getRelUncertUncorr(pT,mass,eta, p)) + fabs(getRelUncertBias(pT,mass,eta, p));
}

double UJUncertaintyProvider::getRelUncertUncorr(double pT, double mass, double eta, e_prop p) {
  return getRelUncertComponent(0, pT, mass, eta, p); 
}

double UJUncertaintyProvider::getRelUncertBias(double pT, double mass, double eta, e_prop p) {
  // Cast parameters to void to tell the compiler that we know the variables aren't used right now
  (void)pT;
  (void)eta;
  (void)p;
  (void)mass;
  return 0.0;
}

double UJUncertaintyProvider::getRelNPVOffsetTerm(double pT, double mass, double eta, double NPV, e_prop p) {
  // Cast parameters to void to tell the compiler that we know the variables aren't used right now
  (void)pT;
  (void)mass;
  (void)eta;
  (void)p;
  (void)NPV;
  return 0.0;
}

double UJUncertaintyProvider::getRelMuOffsetTerm(double pT, double mass, double eta, double mu, e_prop p) {
  // Cast parameters to void to tell the compiler that we know the variables aren't used right now
  (void)pT;
  (void)mass;
  (void)eta;
  (void)mu;
  (void)p;
  return 0.0;
}

double UJUncertaintyProvider::getRelOffsetUncert(double pT, double mass, double eta, double NPV, double mu, e_prop p) {
  double muUncert = getRelMuOffsetTerm(pT,mass,eta,mu,p);
  double npvUncert = getRelNPVOffsetTerm(pT,mass,eta,NPV,p);
  return sqrt( muUncert*muUncert + npvUncert*npvUncert );
}

//function for easier looping on nuisance parameters
//these are the names that should be used as components when calling getRelUncert
std::vector<TString> UJUncertaintyProvider::getComponentNames() {  
  return m_compNames;
}

TString UJUncertaintyProvider::getComponentName(int iComp) {
  if (iComp>=int(m_compNames.size()))
    Fatal("UJUncertaintyProvider::getComponentName()",
	  "You are asking for comp. %d. Only %d components available.",iComp+1,int(m_compNames.size()));
  return m_compNames.at(iComp);
}



//these are the titles of the nuisance parameters (index-parallel with the vector of names)
std::vector<TString> UJUncertaintyProvider::getComponentDescriptions() {
  return m_compDesc;
}

bool UJUncertaintyProvider::isIntercalibrationComponent (const TString & component) {
  return component.Contains("etaIntercalibration", TString::kIgnoreCase);
}

bool UJUncertaintyProvider::isInSituComponent (const TString & component) {
  return ( component.Contains("InSitu", TString::kIgnoreCase) || 
	   component.Contains("Zjet", TString::kIgnoreCase) || 
	   component.Contains("MPF", TString::kIgnoreCase) ||
	   component.Contains("MJB", TString::kIgnoreCase) );
}

bool UJUncertaintyProvider::isNonClosureComponent (const TString & component) {
  return ( component.Contains("NonClosure", TString::kIgnoreCase) );
}

bool UJUncertaintyProvider::isPileupComponent (const TString & component) {
  return ( component.Contains("Pileup", TString::kIgnoreCase) );
}

StrV UJUncertaintyProvider::Vectorize(TString str, TString sep) {
  StrV result; TObjArray *strings = str.Tokenize(sep.Data());
  if (strings->GetEntries()==0) return result;
  TIter istr(strings);
  while (TObjString* os=(TObjString*)istr())
    if (os->GetString()[0]!='#') result.push_back(os->GetString());
    else break;
  delete strings;
  return result;
}

