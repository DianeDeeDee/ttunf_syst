/**
 * This class provides the interface to obtain an event weight corresponding to the correction of the Wjet light and heavy
 * flavor components with a normalization base from Charge Asymmetry.
 *
 * Two input files are needed to be defined in the constructor of the tool :
 *
 *   (1) a *.root file containing Fbb,Fcc,Fc,Fll and ca_norm factors for nominal and each systematic for electron and muon channel
 *   (2) a *.txt file containing the PRETAGED Wjet event yields (no NOT use TAGGED yields)
 *
 * The user has to provide the paths to these two files in the constructor, afterwards only the function GetWjetSFWeight() is to be called.
 *
 * Per default the Fxx factors are taken from the 2exclusive jet bin for all jet multiplicities. Only the charge asymmetry normalization is 
 * taken per jet bin. The default setup herefore uses the 1ex, 2ex, 3ex, 4ex and 5in (ex=exclusive, in=inclusive) jet bins. One can trigger the
 * usage of the 3in or 4in bins with the UseInclusive3JetBin() or UseInclusive4JetBin() function respectively (do not use both in parallel!) 
 */

#include "ttResoSingleLepton/HFsys.h"

//==============================================================================================================================================
// Class constructor calling Init_* functions
//
HFsys::HFsys(TString yields_file, TString sf_file, HFSYS::SELECTION sel) : wjyieldfile(yields_file),
									   sffile(sf_file),
									   selection(sel),
									   wjy_init(false), 
									   sf_init(false), 
									   sys_init(false),
									   use3incl(false),
									   use4incl(false) {
  this->Init();
}


//==============================================================================================================================================
// Trigger usage of 3 inclusive jet bin (3in) instead of exclusive jet bins (3ex+4ex+5in) separately
//
void HFsys::UseInclusive3JetBin() {
  use3incl = true;
  if(use4incl) {
    cout << " WARNING in ttResoSingleLepton/HFsys::UseInclusive3JetBin() : 4jet inclusive already enabled, must not use both in parallel. Disabling usage of inclusive 4jet bin" << endl;
    use4incl = false;
  }
  cout << " INFO in ttResoSingleLepton/HFsys : Using inclusive 3jet bin (3in) instead of exclusive jet bins (3ex+4ex+5in) separately" << endl;
}

//==============================================================================================================================================
// Trigger usage of 4 inclusive jet bin (4in) instead of exclusive jet bins (4ex+5in) separately
//
void HFsys::UseInclusive4JetBin() {
  use4incl = true;
  if(use3incl) {
    cout << " WARNING in ttResoSingleLepton/HFsys::UseInclusive4JetBin() : 3jet inclusive already enabled, must not use both in parallel. Disabling usage of inclusive 3jet bin" << endl;
    use3incl = false;
  }
  cout << " INFO in ttResoSingleLepton/HFsys : Using inclusive 4jet bin (4in) instead of exclusive jet bins (4ex+5in) separately" << endl;
}

//==============================================================================================================================================
// Gets the event weight for a given systematic, lepton channel, hfor flag and jetbin
//
double HFsys::GetWjetSFWeight(TString sys, TString chan, int hfor, int njets, bool applyCA) {

  double weight = 1.0;

  HFSYS::FLAVOR flav   = this->Flavor(hfor);
  int           jetbin = this->JetBin(njets);

  if(!wjy_init) {
    if(nwarnings < 5) {
      cout << " WARNING in ttResoSingleLepton/HFsys::GetWjetSFWeight() : Wjet Yields not properly initialized, returning default weight = 1.0" << endl;
      nwarnings++;
    } else if(nwarnings == 5) {
      cout << " WARNING in ttResoSingleLepton/HFsys::GetWjetSFWeight() : Wjet Yields not properly initialized, returning default weight = 1.0. Last warning!" << endl;
      nwarnings++;
    }
  } else if(!sf_init) {
    if(nwarnings < 5) {
      cout << " WARNING in ttResoSingleLepton/HFsys::GetWjetSFWeight() : HFSys SFs not properly initialized, returning default weight = 1.0" << endl;
      nwarnings++;
    } else if(nwarnings == 5) {
      cout << " WARNING in ttResoSingleLepton/HFsys::GetWjetSFWeight() : HFSys SFs not properly initialized, returning default weight = 1.0. Last warning!" << endl;
      nwarnings++;
    }
  } else if(flav == HFSYS::INDEF) {
    if(nwarnings < 5) {
      cout << " WARNING in ttResoSingleLepton/HFsys::GetWjetSFWeight() : Wjet flavor not properly initialized, returning default weight = 1.0" << endl;
      nwarnings++;
    } else if(nwarnings == 5) {
      cout << " WARNING in ttResoSingleLepton/HFsys::GetWjetSFWeight() : Wjet flavor not properly initialized, returning default weight = 1.0. Last warning!" << endl;
      nwarnings++;
    }
  } else if(find(syslist.begin(), syslist.end(), sys) == syslist.end()) {
    if(nwarnings < 5) {
      cout << " WARNING in ttResoSingleLepton/HFsys::GetWjetSFWeight() : Unkown systematic = " << sys << " provided, returning default weight = 1.0" << endl;
      nwarnings++;
    } else if(nwarnings == 5) {
      cout << " WARNING in ttResoSingleLepton/HFsys::GetWjetSFWeight() : Unkown systematic = " << sys << " provided, returning default weight = 1.0. Last warning!" << endl;
      nwarnings++;
    }
  } else if(find(lepchanlist.begin(), lepchanlist.end(), chan) == lepchanlist.end()) {
    if(nwarnings < 5) {
      cout << " WARNING in ttResoSingleLepton/HFsys::GetWjetSFWeight() : Unkown lepton channel = " << chan << " provided, returning default weight = 1.0" << endl;
      nwarnings++;
    } else if(nwarnings == 5) {
      cout << " WARNING in ttResoSingleLepton/HFsys::GetWjetSFWeight() : Unkown lepton channel = " << chan << " provided, returning default weight = 1.0. Last warning!" << endl;
      nwarnings++;
    }
  } else if(jetbin == 0) {
    //not fatal, but also returns default value weight=1.0
  } else if(flav == HFSYS::NOW || flav == HFSYS::KILL) {
    //not fatal, but also returns default value weight=1.0
  } else { //all checks passed, proceed

    //Read out scalings for given setup
    double Fbb, Fcc, Fc, Fll, cafac;
    this->Scalings(sys, chan, jetbin, Fbb, Fcc, Fc, Fll, cafac);

    //Read out wjet pretag yields for given setup
    map<HFSYS::FLAVOR, double> yields_in;
    Yields(chan, jetbin, yields_in);
    
    map<HFSYS::FLAVOR, double> yields_out;
    yields_out[HFSYS::BB] = Fbb * yields_in[HFSYS::BB];
    yields_out[HFSYS::CC] = Fcc * yields_in[HFSYS::CC];
    yields_out[HFSYS::C]  = Fc  * yields_in[HFSYS::C];
    yields_out[HFSYS::LL] = Fll * yields_in[HFSYS::LL];

    //Normalize
    double sum_in  = yields_in [HFSYS::BB] + yields_in [HFSYS::CC] + yields_in [HFSYS::C] + yields_in [HFSYS::LL];
    double sum_out = yields_out[HFSYS::BB] + yields_out[HFSYS::CC] + yields_out[HFSYS::C] + yields_out[HFSYS::LL];
    double norm = 1.0;
    if(sum_out > 0) {
      norm = sum_in/sum_out;
      yields_out[HFSYS::BB] *= norm;
      yields_out[HFSYS::CC] *= norm;
      yields_out[HFSYS::C]  *= norm;
      yields_out[HFSYS::LL] *= norm;
    } else {
      cout << " WARNING in ttResoSingleLepton/HFsys::GetWjetSFWeight() : Sum of scaled WjetYields is ZERO! Skipping normalization" << endl;
    }
    
    //calculate actual weight
    weight = yields_out[flav] / yields_in[flav];

    //apply charge asymmetry normalization only if flag enabled
    if(applyCA) weight *= cafac;
  }

  return weight;
}

//##############################################################################################################################################
// PRIVATE
//##############################################################################################################################################


//==============================================================================================================================================
// Determine W+jets flavor from hfor flag
//
HFSYS::FLAVOR HFsys::Flavor(int hfor) {

  HFSYS::FLAVOR flav = HFSYS::INDEF;

  switch(hfor) {
    case -1 : flav = HFSYS::NOW;  break;
    case  0 : flav = HFSYS::BB;   break;
    case  1 : flav = HFSYS::CC;   break;
    case  2 : flav = HFSYS::C;    break;
    case  3 : flav = HFSYS::LL;   break;
    case  4 : flav = HFSYS::KILL; break;
    default : flav = HFSYS::INDEF; 
              cout << " WARNING in ttResoSingleLepton/HFsys :: Unkown hfor = " << hfor << " provided (expecting E [-1, 0, 1, 2, 3, 4])" << endl;    
  }

  return flav;
}


//==============================================================================================================================================
// Determine jet bin from number of jets
//
int HFsys::JetBin(int njets) {

  int jetbin = 0;

  //exclusive jetbin usage
  if(njets < 5) {
    jetbin = njets; //0ex/1ex/2ex/3ex/4ex
  } else {
    jetbin = 5; //5in
  }

  //inclusive jetbin usage
  if(use3incl && njets >= 3) {
    jetbin = 6; //3in
  }
  if(use4incl && njets >= 4) {
    jetbin = 7; //4in
  }

  return jetbin;
}


//==============================================================================================================================================
// Set W+Jet input yields for given lepton channel and jet bin
//
void HFsys::Yields(TString chan, int jbin, map<HFSYS::FLAVOR, double> &yields) {

  if(chan == "el") {  
    yields[HFSYS::BB] = WJYields_el.at(jbin-1)[HFSYS::BB];
    yields[HFSYS::CC] = WJYields_el.at(jbin-1)[HFSYS::CC];
    yields[HFSYS::C]  = WJYields_el.at(jbin-1)[HFSYS::C];
    yields[HFSYS::LL] = WJYields_el.at(jbin-1)[HFSYS::LL];
  } else if(chan == "mu") {
    yields[HFSYS::BB] = WJYields_mu.at(jbin-1)[HFSYS::BB];
    yields[HFSYS::CC] = WJYields_mu.at(jbin-1)[HFSYS::CC];
    yields[HFSYS::C]  = WJYields_mu.at(jbin-1)[HFSYS::C];
    yields[HFSYS::LL] = WJYields_mu.at(jbin-1)[HFSYS::LL];
  } else {
    cout << " WARNING in ttResoSingleLepton/HFsys :: Unkown lepton channel = " << chan << " provided (expecting E [el, mu])" << endl;    
  }
}

//==============================================================================================================================================
// Set scalings for given systematic and jet bin
//
void HFsys::Scalings(TString sys, TString chan, int jbin, double &Fbb, double &Fcc, double &Fc, double &Fll, double &cafac) {

  if(chan == "el") {
    Fbb = Fxx_el[sys][HFSYS::BB];
    Fcc = Fxx_el[sys][HFSYS::CC];
    Fc  = Fxx_el[sys][HFSYS::C];
    Fll = Fxx_el[sys][HFSYS::LL];
    cafac = ca_el[sys][jbin-1];
  }

  if(chan == "mu") {
    Fbb = Fxx_mu[sys][HFSYS::BB];
    Fcc = Fxx_mu[sys][HFSYS::CC];
    Fc  = Fxx_mu[sys][HFSYS::C];
    Fll = Fxx_mu[sys][HFSYS::LL];
    cafac = ca_mu[sys][jbin-1];
  }
}

//==============================================================================================================================================
// Initialize tool
//
void HFsys::Init() {

  nwarnings = 0;

  lepchanlist.push_back("el");
  lepchanlist.push_back("mu");

  Init_SystematicList();
  Init_Scalings();
  Init_Yields();
}

//==============================================================================================================================================
// Initializes the list of available systematics
//
void HFsys::Init_SystematicList() {

  syslist.push_back("nominal");                                          // nominal without any variation

  syslist.push_back("nlo");                                              // no NNLO corrections applied

  syslist.push_back("qcdup");        syslist.push_back("qcddw");         // QCD variation

  if(selection == HFSYS::RESOLVED || selection == HFSYS::BOOSTED) {
    //WShape
    syslist.push_back("wshapeup");     syslist.push_back("wshapedw");      // combined WShape variation
    syslist.push_back("wshape1up");    syslist.push_back("wshape1dw");     // iqopt3 WShape variation
    syslist.push_back("wshape2up");    syslist.push_back("wshape2dw");     // ptjmin10 WShape variation

    //Btag
    syslist.push_back("btagup");       syslist.push_back("btagdw");        // combined Btag variation
    syslist.push_back("btagcup");      syslist.push_back("btagcdw");       // combined Btag c-flavor variation
    syslist.push_back("btaglup");      syslist.push_back("btagldw");       // combined Btag light-flavor variation
    syslist.push_back("btagbup");      syslist.push_back("btagbdw");       // combined Btag b-flavor variation
    syslist.push_back("btagb6up");     syslist.push_back("btagb6dw");      // Btag NumVariation6
    syslist.push_back("btagb7up");     syslist.push_back("btagb7dw");      // Btag NumVariation7
    syslist.push_back("btagb8up");     syslist.push_back("btagb8dw");      // Btag NumVariation8
    syslist.push_back("btagb9up");     syslist.push_back("btagb9dw");      // Btag NumVariation9
    syslist.push_back("btagb10up");    syslist.push_back("btagb10dw");     // Btag NumVariation10
  
    //Anti-kt 0.4 jet energy scale (JES)
    syslist.push_back("jesup");        syslist.push_back("jesdw");         // combined JES variation
    syslist.push_back("jes4up");       syslist.push_back("jes4dw");        // JES component # 4 (EffectiveNP_Modelling1)
    syslist.push_back("jes8up");       syslist.push_back("jes8dw");        // JES component # 8 (EffectiveNP_Detector1)
    syslist.push_back("jes13up");      syslist.push_back("jes13dw");       // JES component #13 (EtaIntercalibration_Modelling)
    syslist.push_back("jes20up");      syslist.push_back("jes20dw");       // JES component #20 (PileupRhoTopology)
    syslist.push_back("jes22up");      syslist.push_back("jes22dw");       // JES component #22 (FlavorComp)
    syslist.push_back("jes23up");      syslist.push_back("jes23dw");       // JES component #23 (FlavorResponse)
    syslist.push_back("jes24up");      syslist.push_back("jes24dw");       // JES component #24 (BJES)
    syslist.push_back("jes99up");      syslist.push_back("jes99dw");       // combination of all non-explicit JES components

    //Anti-kt 0.4 jet energy resolution (JER)
    syslist.push_back("jerup");        syslist.push_back("jerdw");         // JER variation

    //Anti-kt 0.4 jet vertex fraction (JVF)
    syslist.push_back("jvfup");        syslist.push_back("jvfdw");         // JVF variation
  }
 
  //##### BOOSTED specific systematics ##############################################################################################
  //##                                                                                                                             ##
  if(selection == HFSYS::BOOSTED) {

    //Anti-kt 1.0 jet energy scale (JES)
    syslist.push_back("jesboostup");   syslist.push_back("jesboostdw");    // combined boosted JES variation
    syslist.push_back("jesboost1up");  syslist.push_back("jesboost1dw");   // boosted JES component # 1 (dataMC)
    syslist.push_back("jesboost13up"); syslist.push_back("jesboost13dw");  // boosted JES component #13 (topology)
    syslist.push_back("jesboost99up"); syslist.push_back("jesboost99dw");  // combination of all non-explicit boosted JES components

    //Anti-kt 1.0 jet energy resolution (JER)
    syslist.push_back("jerboostup");   syslist.push_back("jerboostdw");    // boosted JER variation
  }
  //##                                                                                                                             ##
  //##### BOOSTED specific systematics ##############################################################################################

  sys_init = true;
}


//==============================================================================================================================================
// Initializes scaling factors, reading factors for each systematic from *.root file
//
void HFsys::Init_Scalings() {

  TFile *f = TFile::Open(sffile, "READ");

  if(!f) {
    cout << " WARNING in ttResoSingleLepton/HFsys :: Can not open SF file " << sffile << endl;
    return;
  }

  //for each lepton channel
  for(vector<TString>::iterator it = lepchanlist.begin() ; it != lepchanlist.end(); ++it) {
    TString channel = *it;

    //for each systematic
    for(vector<TString>::iterator jt = syslist.begin() ; jt != syslist.end(); ++jt) {
      TString sys = *jt;

      TString hname = sys+"_"+channel;
      TH2F *h2 = (TH2F*)f->Get(hname);
      
      bool sys_aval;
      if(h2 == 0) {
	cout << " WARNING :: ttResoSingleLepton::HFsys::Init_Scalings() : Systematic " << sys << " not (yet) available in SF file " << sffile << endl;
	sys_aval = false;
      } else {
	sys_aval = true;
      }

      //Set factor if sys available, else use default 1.0
      //y_idx = 6 for 2ex jet bin
      if(channel == "el") {
	Fxx_el[sys][HFSYS::BB] = sys_aval ? h2->GetBinContent(1, 6) : 1.0;
	Fxx_el[sys][HFSYS::CC] = sys_aval ? h2->GetBinContent(2, 6) : 1.0;
	Fxx_el[sys][HFSYS::C]  = sys_aval ? h2->GetBinContent(3, 6) : 1.0;
	Fxx_el[sys][HFSYS::LL] = sys_aval ? h2->GetBinContent(4, 6) : 1.0;
      } 
      if(channel == "mu") {
	Fxx_mu[sys][HFSYS::BB] = sys_aval ? h2->GetBinContent(1, 6) : 1.0;
	Fxx_mu[sys][HFSYS::CC] = sys_aval ? h2->GetBinContent(2, 6) : 1.0;
	Fxx_mu[sys][HFSYS::C]  = sys_aval ? h2->GetBinContent(3, 6) : 1.0;
	Fxx_mu[sys][HFSYS::LL] = sys_aval ? h2->GetBinContent(4, 6) : 1.0;
      }      

      int nybins = 7;
      for(int ibin = 0; ibin < nybins; ibin++) {
	//x_idx = 5 for CAnorm bin, y_idx reversed
	if(channel == "el") ca_el[sys].push_back(sys_aval ? h2->GetBinContent(5, nybins-ibin) : 1.0);
	if(channel == "mu") ca_mu[sys].push_back(sys_aval ? h2->GetBinContent(5, nybins-ibin) : 1.0);
      }
    }//end sys loop
  }//end lepchan loop
  
  f->Close();

  sf_init = true;
}

//==============================================================================================================================================
// Initialize pretag wjet event yields provided in yield file
//
void HFsys::Init_Yields() {

  ifstream ifs(wjyieldfile);

 if(!ifs.good()) {
    cout << " WARNING :: Can not open file pretag wjet yields file " << wjyieldfile << endl;
    return;
  }
  
  string line;
  TString chan = "none";
  while(getline(ifs, line)) {

    vector<string> tokens;
    Tokenize(line, tokens, " \t");

    //skip empty lines
    if(tokens.size() == 0) continue;
    
    //trigger electorn/muon channel factors
    if(tokens.at(0) == "el" || tokens.at(0) == "mu") { 
      chan = tokens.at(0);
    }

    //channel needs to be triggered before reading numbers
    if(chan == "none") {
      cout << " WARNING :: Unkown format in Yields file" << endl;
      return;
    }
    
    int jbin = -1;
    if(tokens.at(0) == "1ex") jbin = 1;
    if(tokens.at(0) == "2ex") jbin = 2;
    if(tokens.at(0) == "3ex") jbin = 3;
    if(tokens.at(0) == "4ex") jbin = 4;
    if(tokens.at(0) == "5in") jbin = 5;
    if(tokens.at(0) == "3in") jbin = 6;
    if(tokens.at(0) == "4in") jbin = 7;

    if(jbin < 0) continue;

    if(tokens.size() != 5) {
      cout << " WARNING :: Unkown format in Yields file" << endl;
      return;
    }

    if(chan == "el") {
      map<HFSYS::FLAVOR, double> yields_el;
      yields_el[HFSYS::BB] = (double)atof(tokens.at(1).c_str());
      yields_el[HFSYS::CC] = (double)atof(tokens.at(2).c_str());
      yields_el[HFSYS::C]  = (double)atof(tokens.at(3).c_str());
      yields_el[HFSYS::LL] = (double)atof(tokens.at(4).c_str());
      WJYields_el.push_back(yields_el);
    }
    if(chan == "mu") {
      map<HFSYS::FLAVOR, double> yields_mu;
      yields_mu[HFSYS::BB] = (double)atof(tokens.at(1).c_str());
      yields_mu[HFSYS::CC] = (double)atof(tokens.at(2).c_str());
      yields_mu[HFSYS::C]  = (double)atof(tokens.at(3).c_str());
      yields_mu[HFSYS::LL] = (double)atof(tokens.at(4).c_str());
      WJYields_mu.push_back(yields_mu);
    }    
  }//end of getline
  ifs.close();

  wjy_init = true;
}

//==============================================================================================================================================
// Helper function to split TString
//
void HFsys::Tokenize(const string& str, vector<string>& tokens, const string& delimiters) {
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type pos = str.find_first_of(delimiters, lastPos);
  while(string::npos != pos || string::npos != lastPos) {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}
