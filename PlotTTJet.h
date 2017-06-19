#ifndef PLOTTTJET_H
#define PLOTTTJET_H

#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "Event.h"
#include "Plot.h"

#include "SD-tosvn/TopGluonModel.h"
#include "SD-tosvn/BackgroundModel.h"
#include "SD-tosvn/Deconstruct.h"
#include "SD-tosvn/ISRModel.h"
#include "SD-tosvn/AnalysisParameters.h"

#include "LargeJet.h"

class PlotTTJet : public Plot {
  public:
    PlotTTJet(const std::string &filename, const std::vector<std::string> &systs);
    virtual ~PlotTTJet();

    void run(const Event &e, double weight, double pweight, const std::string &s);

  protected:

    //std::vector<TLorentzVector> applySjUnc(const std::vector<TLorentzVector> &sj, const std::string &s);
    std::vector<LargeJet::Subjet> applySjUnc(const std::vector<LargeJet::Subjet> &sj, const std::vector<LargeJet::Subjet> &sju, const std::string &s);

    Deconstruction::TopGluonModel *signal;
    Deconstruction::BackgroundModel *background;
    Deconstruction::ISRModel *isr;
    Deconstruction::Deconstruct *deconstruct;
    AnalysisParameters param;

    Deconstruction::TopGluonModel *signalb;
    Deconstruction::BackgroundModel *backgroundb;
    Deconstruction::ISRModel *isrb;
    Deconstruction::Deconstruct *deconstructb;
    AnalysisParameters paramb;

    int    m_tree_syst;
    int    m_tree_nsub;
    double m_tree_chi;
    double m_tree_chib;
    int    m_tree_bmatch;
    double m_tree_pt;
    double m_tree_eta;
    double m_tree_truept;
    double m_tree_trueeta;
    double m_tree_m;
    double m_tree_d12;

    int    m_tree_cansub;
    double m_tree_cachi;
    double m_tree_cachib;
    int    m_tree_cabmatch;
    double m_tree_capt;
    double m_tree_caeta;
    double m_tree_catruept;
    double m_tree_catrueeta;
    double m_tree_cahtt;

    double m_tree_mu;
    int    m_tree_npv;
    double m_tree_w;
};

#endif

