#ifndef PLOTSEMILEP_H
#define PLOTSEMILEP_H

#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "Event.h"
#include "Plot.h"

class PlotSemilep : public Plot {
  public:
    PlotSemilep(const std::string &filename, bool electron, const std::vector<std::string> &systs);
    virtual ~PlotSemilep();

    void run(const Event &e, double weight, double pweight, const std::string &s);
    bool getWFromLeptonicDecay(TLorentzVector &momLepton, float NuX, float NuY, TLorentzVector &momNu, TLorentzVector &momW, float lepMass);
  protected:
    bool m_electron;
};

#endif

