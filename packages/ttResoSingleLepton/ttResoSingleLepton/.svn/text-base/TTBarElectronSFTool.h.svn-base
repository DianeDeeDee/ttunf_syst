#ifndef __TTBARELECTRONSFTOOL_H__
#define __TTBARELECTRONSFTOOL_H__

#include <utility>

#include "ElectronEfficiencyCorrection/TElectronEfficiencyCorrectionTool.h"

#include "TFile.h"

// use pairs as value +/- uncertainty
typedef std::pair<double, double> ValUncert;

class TTBarElectronSFTool {
    protected:
        // files and histograms derived specifically for the l+jets ttbar
        // resonances analysis
        TFile *fIDSFFile;
        TFile *fTrigIsoSFFile;
        TFile *fFarJetTrigIsoSFFile;

        TH1D *fPtIDSFHist;
        TH1D *fEtaIDSFHist;

        TH1D *fPtTrigIsoSFHist;
        TH1D *fEtaTrigIsoSFHist;

        TH1D *fFarJetPtTrigIsoSFHist;
        TH1D *fFarJetEtaTrigIsoSFHist;

        // e/gamma tools
        Root::TElectronEfficiencyCorrectionTool *fRecoSFTool;
        Root::TElectronEfficiencyCorrectionTool *fIDSFTool;

        ValUncert ContentAndError(
                TH1D *h,
                double x,
                bool allow_overflow);

        bool HaveSF(double pt,
                double eta,
                double jet_dr,
                double jet_pt);

    public:
        TTBarElectronSFTool();

        ~TTBarElectronSFTool();

        ValUncert GetSF(
                int RunNumber,
                double pt,
                double eta,
                double jet_dr,
                double jet_pt,
                bool afii);

        ValUncert GetRecoSF(
                int RunNumber,
                double pt,
                double eta,
                bool afii);

        ValUncert GetIDSF(
                int RunNumber,
                double pt,
                double eta,
                double jet_dr,
                double jet_pt,
                bool afii);

        ValUncert GetTrigIsoSF(
                int RunNumber,
                double pt,
                double eta,
                double jet_dr,
                double jet_pt,
                bool afii);

        void SetIDSFFile(const std::string &fname);
        void SetTrigIsoSFFile(const std::string &fname);
        void SetFarJetTrigIsoSFFile(const std::string &fname);

        void SetRecoSFTool(
                Root::TElectronEfficiencyCorrectionTool *sftool);

        void SetIDSFTool(
                Root::TElectronEfficiencyCorrectionTool *sftool);
};


// simple numeric operators on uncertain values.
ValUncert VUAdd(const ValUncert &u, const ValUncert &v);
ValUncert VUNeg(const ValUncert &u);
ValUncert VUSub(const ValUncert &u, const ValUncert &v);
ValUncert VUMult(const ValUncert &u, const ValUncert &v);
ValUncert VUDiv(const ValUncert &u, const ValUncert &v);

// cast double as a value with no uncertainty.
ValUncert VUConst(const double c);

#endif
