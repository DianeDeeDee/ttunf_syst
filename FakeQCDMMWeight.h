#ifndef FAKEQCDMMWEIGHT_H
#define FAKEQCDMMWEIGHT_H

#include <string>
#include "TFile.h"
#include "TH1.h"

#include "TAxis.h"

class FakeQCDMMWeight {

  public:

    FakeQCDMMWeight(const std::string &filename, const std::string &filename_pre = "");

    virtual ~FakeQCDMMWeight();

  protected:

    // this is quite general, but it could be easier to override it to a simpler function prototype that calls it
    double weight(bool electron, bool passTight, float *x, float *y, int nx, int ny, bool nom = true, bool eff_error = true, bool up = true);

    // these are quite general given the "general" binning of the efficiency/fake rate histogram: x and y
    // no need to override those
    double eff(bool electron, float *x, int nx = 1, bool nom = true, bool up = true);
    double fake(bool electron, float *y, int ny = 1, bool nom = true, bool up = true);

    // The following six functions must be coded by a chield of FakeQCDMMWeight
    // use parameters in x to get the proper histogram used for the efficiency/fake rate estimate
    // this exists because we might select one histogram or another for different parameter sets
    virtual TH1 *getEffHist(bool electron, float *x, int nx) = 0;
    virtual TH1 *getFakeHist(bool electron, float *x, int nx) = 0;

    // some of the parameters in x might be used to select the histogram for the efficiency/fake rate estimate
    // this returns the parameters actually used in the histogram (at most three: a, b, c)
    virtual void getEffBin(bool electron, float *x, int nx, float &a, float &b, float &c, int &n) = 0;
    virtual void getFakeBin(bool electron, float *x, int nx, float &a, float &b, float &c, int &n) = 0;


    TFile *m_file;
    TFile *m_file_pre;

  private:

    double restrictRange(double a, TAxis *axis);

};

#endif

