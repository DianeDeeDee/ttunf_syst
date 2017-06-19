#ifndef PLOT_H
#define PLOT_H

#include <string>
#include "TH1F.h"
#include "TFile.h"
#include "Event.h"
#include "HistogramService.h"

class Plot {
  public:
    Plot(const std::string &filename, const std::vector<std::string> &systs);
    virtual ~Plot();

    virtual void run(const Event &e, double weight, double pweight, const std::string &s) = 0;

  protected:
    std::string m_filename;
    HistogramService m_hSvc;

    float topPtWeight(const Event &e);

};

#endif

