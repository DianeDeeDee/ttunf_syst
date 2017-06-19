#ifndef HFSYS_H
#define HFSYS_H

#include "TString.h"
#include "TFile.h"
#include "TH2F.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <stdlib.h>
#include <string>
#include <algorithm>

using namespace std;

namespace HFSYS {
  enum FLAVOR {BB, CC, C, LL, KILL, NOW, INDEF};
  enum SELECTION {RESOLVED, BOOSTED};
};

class HFsys {

  //-----------------------------------------------------------------------------------------------------------------------
  // PUBLIC 
  //-----------------------------------------------------------------------------------------------------------------------

 public:

  /**
   * Class constructor calling Init_* functions
   *
   * @param yields_file  full path to the *.txt file containing pretag wjet yields
   * @param sf_file      full path to the *.root file containing scaling factors
   */
  HFsys(TString yields_file, TString sf_file, HFSYS::SELECTION sel);

  /**
   * Gets the event weight for a given systematic, lepton channel, hfor flag and jetbin
   *
   * @param sys   systematic
   * @param chan  lepton channel
   * @param hfor  heavy flavor overlap removal flag
   * @param jbin  jet bin
   *
   * @return eventweight
   */
  double GetWjetSFWeight(TString sys, TString chan, int hfor, int njets, bool applyCA = true);

  /**
   * Trigger usage of 3 inclusive jet bin (3in) instead of exclusive jet bins (3ex+4ex+5in) separately
   */
  void UseInclusive3JetBin();

  /**
   * Trigger usage of 4 inclusive jet bin (4in) instead of exclusive jet bins (4ex+5in) separately
   */
  void UseInclusive4JetBin();


  //-----------------------------------------------------------------------------------------------------------------------
  // PRIVATE
  //-----------------------------------------------------------------------------------------------------------------------

 private:

  TString wjyieldfile;                               ///< full path to wjet pretag *.txt yield file
  TString sffile;                                    ///< full path to *.root sf file
  HFSYS::SELECTION selection;                        ///< Selection enum
  bool wjy_init;                                     ///< flag for initialization of wjet pretag yields
  bool sf_init;                                      ///< flag for initialization of scaling factors
  bool sys_init;                                     ///< flag for initialization of systematic list
  bool use3incl;                                     ///< flag for usage of 3 inclusive bin instead of 3ex+4ex+5in separately
  bool use4incl;                                     ///< flag for usage of 4 inclusive bin instead of 4ex+5in separately

  int nwarnings;                                     ///< counting number of warnings

  vector<TString> syslist;                           ///< list of valid systematics
  vector<TString> lepchanlist;                       ///< list of valid lepton channels

  map<TString, map<HFSYS::FLAVOR, double> > Fxx_el;  ///< Fbb, Fcc, Fc and Fll factors for 2ex jetbin for each systematic and flavor for electron channel
  map<TString, map<HFSYS::FLAVOR, double> > Fxx_mu;  ///< Fbb, Fcc, Fc and Fll factors for 2ex jetbin for each systematic and flavor for muon channel

  map<TString, vector<double> > ca_el;               ///< charge asymmetry factors for each systematic and jetbin for electron channel
  map<TString, vector<double> > ca_mu;               ///< charge asymmetry factors for each systematic and jetbin for muon channel

  vector<map<HFSYS::FLAVOR, double> > WJYields_el;   ///< pretag Wjet yields for each jetbin and heavy flavor for electron channel
  vector<map<HFSYS::FLAVOR, double> > WJYields_mu;   ///< pretag Wjet yields for each jetbin and heavy flavor for muon channel


  /**
   * Determine W+jets flavor from hfor flag
   *
   * @param hfor  heavy flavor overlap removal flag
   * 
   * @return flavor enum
   */
  HFSYS::FLAVOR Flavor(int hfor);

  /**
   * Determine jet bin from number of jets
   *
   * @param  number of jets
   *
   * @return jetbin
   */
  int JetBin(int njets);

  /**
   * Gets the wjet pretag yields for the given lepton channel and jetbin
   *
   * @param chan    lepton channel enum
   * @param jbin    jet bin
   * @param yields  will be filled with wjet pretag yields for given lepton channel and jetbin
   */
  void Yields(TString chan, int jbin, map<HFSYS::FLAVOR, double> &yields);  

  /**
   * Reads out Fxx factors and charge asymmetry factors for a given systematic
   *
   * @param sys   systematic enum
   * @param chan  lepton channel enum
   * @param jbin  jet bin
   * @param Fbb   will be filled with Fbb factor for 2ex jbin for given systematic
   * @param Fcc   will be filled with Fcc factor for 2ex jbin for given systematic
   * @param Fc    will be filled with Fc  factor for 2ex jbin for given systematic
   * @param Fll   will be filled with Fll factor for 2ex jbin for given systematic
   * @param ca    will be filled with charge asymmetry factors for all jbins for given systematic  
   */
  void Scalings(TString sys, TString chan, int jbin, double &Fbb, double &Fcc, double &Fc, double &Fll, double &cafac);

  /**
   * Initialize tool
   */
  void Init();

  /**
   * Initializes the list of available systematics
   */
  void Init_SystematicList();

  /**
   * Initializes scaling factors, reading factors for each systematic provided in *.root sffile
   */
  void Init_Scalings();

  /**
   * Initializes pretag wjet event yields provided in *.txt yieldfile
   */
  void Init_Yields();

  /**
   * Helper function to split TString
   *
   * @param str        string to be tokenized
   * @param tokens     vector to be filled with tokens
   * @param delimiters list of delimiters used to split
   */
  void Tokenize(const string& str, vector<string>& tokens, const string& delimiters);  
};

#endif

