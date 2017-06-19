#ifndef LARGERJETRESOTOOL_H
#define LARGERJETRESOTOOL_H

#include <utility>
#include "TFile.h"
#include "TH1F.h"
#include "TRandom3.h"

using namespace std;

class LargeRJetResoTool {
 public:	
  enum resType { JER = 0, JMR };
  LargeRJetResoTool(const string &fname="LargeRJetReso.root");
  ~LargeRJetResoTool();	
  
  // only two method at this time. but plan to orgainze a bit more. 
  float GetResolution(float pTGeV, float eta, resType iType);

  float GetSmearFactor(float pTGeV, float eta, resType iType, UInt_t seed=0 );

  void SetSeed(UInt_t seed) { rdn.SetSeed(seed); }

 private:

  TFile* inputFile;
  TH1F* JEREta1;
  TH1F* JEREta2;
  TH1F* JEREta3;
  TH1F* JMREta1;
  TH1F* JMREta2;
  TH1F* JMREta3;

  TRandom3 rdn;
  
};

#endif
