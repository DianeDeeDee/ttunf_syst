#ifndef FAKEQCDMMBOOSTED_H
#define FAKEQCDMMBOOSTED_H

#include <string>
#include "TFile.h"
#include "FakeQCDMMWeight.h"

class FakeQCDMMBoosted : public FakeQCDMMWeight {

  public:
    FakeQCDMMBoosted(const std::string &filename, const std::string &filename_pre = "");
    virtual ~FakeQCDMMBoosted();

    // not really necessary, but useful for debugging
    double getEff(bool electron, bool btagging, float lepPt, float lepJetdR, bool nom = true, bool up = true); // x1 = lepPt
    double getFake(bool electron, bool btagging, float lepPt, float lepJetdR, float closeJetPt, float sd0, bool nom = true, bool up = true); // y1 = lepPt

    // useful prototype format: this is the function to be used
    double getWeight(bool electron, bool passTight, bool btagging,
                     float lepPt, float lepJetdR, float closeJetPt, float sd0,
                     bool nom = true, bool eff_error = true, bool up = true);

  protected:
    // must be overridden
    TH1 *getEffHist(bool electron, float *x, int nx);
    TH1 *getFakeHist(bool electron, float *x, int nx);

    void getEffBin(bool electron, float *x, int nx, float &a, float &b, float &c, int &n);
    void getFakeBin(bool electron, float *x, int nx, float &a, float &b, float &c, int &n);

};

#endif

