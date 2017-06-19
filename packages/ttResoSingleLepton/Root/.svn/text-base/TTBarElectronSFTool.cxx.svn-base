#include <string>
#include <iostream>

#include "ttResoSingleLepton/TTBarElectronSFTool.h"

#include "TMath.h"

using std::string;

using std::cout;
using std::cerr;
using std::endl;

TTBarElectronSFTool::TTBarElectronSFTool() {

    fIDSFFile = NULL;
    fTrigIsoSFFile = NULL;
    fFarJetTrigIsoSFFile = NULL;

    return;
}

TTBarElectronSFTool::~TTBarElectronSFTool() {

    if (fIDSFFile)
        fIDSFFile->Close();

    if (fTrigIsoSFFile)
        fTrigIsoSFFile->Close();

    if (fFarJetTrigIsoSFFile)
        fFarJetTrigIsoSFFile->Close();

    return;
}

void TTBarElectronSFTool::SetIDSFFile(const string &fname) {
    // clean up old file.
    if (fIDSFFile)
        fIDSFFile->Close();

    fIDSFFile = TFile::Open(fname.c_str());

    if (!fIDSFFile) {
        cerr << "TTBarElectronSFTool: Could not open file " <<
            fname << " for ID scale factor histograms. Exiting." <<
            endl;

        exit(-1);
    }

    TH2D *id_pt_dr_hist = (TH2D *) fIDSFFile->Get("pt_drjet");
    fPtIDSFHist = id_pt_dr_hist->ProjectionX("id_nearjet_pt", 1, 1, "e");

    delete id_pt_dr_hist;

    fEtaIDSFHist = (TH1D *) fIDSFFile->Get("eta");

    return;
}


void TTBarElectronSFTool::SetTrigIsoSFFile(const string &fname) {
    // clean up old file.
    if (fTrigIsoSFFile)
        fTrigIsoSFFile->Close();

    fTrigIsoSFFile = TFile::Open(fname.c_str());

    if (!fTrigIsoSFFile) {
        cerr << "TTBarElectronSFTool: Could not open file " <<
            fname << " for Trig+Iso scale factor histograms. " <<
            "Exiting." << endl;

        exit(-1);
    }

    TH2D *trigiso_pt_dr_hist = (TH2D *) fTrigIsoSFFile->Get("pt_drjet");
    fPtTrigIsoSFHist = trigiso_pt_dr_hist->ProjectionX("trigiso_nearjet_pt", 1, 1, "e");

    delete trigiso_pt_dr_hist;

    fEtaTrigIsoSFHist = (TH1D *) fTrigIsoSFFile->Get("eta");

    return;
}


void TTBarElectronSFTool::SetFarJetTrigIsoSFFile(const string &fname) {
    // clean up old file.
    if (fFarJetTrigIsoSFFile)
        fFarJetTrigIsoSFFile->Close();

    fFarJetTrigIsoSFFile = TFile::Open(fname.c_str());

    if (!fFarJetTrigIsoSFFile) {
        cerr << "TTBarElectronSFTool: Could not open file " <<
            fname << " for Trig+Iso scale factor histograms far " <<
            "from jets. Exiting." << endl;

        exit(-1);
    }

    fFarJetPtTrigIsoSFHist = (TH1D *) fFarJetTrigIsoSFFile->Get("pt");
    fFarJetEtaTrigIsoSFHist = (TH1D *) fFarJetTrigIsoSFFile->Get("eta");

    return;
}

void TTBarElectronSFTool::SetRecoSFTool(
        Root::TElectronEfficiencyCorrectionTool *sftool) {

    fRecoSFTool = sftool;
    return;
}

void TTBarElectronSFTool::SetIDSFTool(
        Root::TElectronEfficiencyCorrectionTool *sftool) {

    fIDSFTool = sftool;
    return;
}

bool TTBarElectronSFTool::HaveSF(double pt, double eta, double jet_dr,
        double jet_pt) {

    double abseta = TMath::Abs(eta);

    if (abseta > 2.47 || (1.37 < abseta && abseta < 1.52)) {

        cerr << "TTBarElectronSFTool: no electron scale factors "
            "available for electrons with |eta| > 2.47 or "
            "1.37 < |eta| 1.52." << endl;

        return false;
    } else if (pt < 25000) {

        cerr << "TTBarElectronSFTool: no electron scale factors "
            "available for electrons with pt < 25 GeV." << endl;

        return false;
    } else if (jet_dr < 0.2) {

        cerr << "TTBarElectronSFTool: no electron scale factors "
            "available for electrons within 0.2 of the nearest jet."
            << endl
            << "The scale factor for \\Delta R = 0.2 will be used."
            << endl;
    }

    return true;
}

ValUncert TTBarElectronSFTool::ContentAndError(
        TH1D *h, double x, bool allow_overflow) {

    int bin = h->FindBin(x);

    if (allow_overflow)
        bin = bin <= h->GetNbinsX() ? bin : h->GetNbinsX();

    double cont = h->GetBinContent(bin);
    double err = h->GetBinError(bin);

    return ValUncert(cont, err);
}

ValUncert TTBarElectronSFTool::GetRecoSF(
        int RunNumber, double pt, double eta, bool afii) {

    if (!HaveSF(pt, eta, 1.0, 100))
        return ValUncert(-999, 999);

    if (!fRecoSFTool) {
        cerr << "TTBarElectronSFTool: please initialize an " <<
            "electron reconstruction SF tool and set it by " <<
            "calling SetRecoSFTool()" << endl;
        return ValUncert(-999, 999);
    }

    PATCore::ParticleDataType::DataType dt;
    if (afii)
        dt = PATCore::ParticleDataType::Fast;
    else
        dt = PATCore::ParticleDataType::Full;

    Root::TResult r = fRecoSFTool->calculate(
            dt, RunNumber, eta, pt);

    return ValUncert(r.getScaleFactor(), r.getTotalUncertainty());
}

ValUncert TTBarElectronSFTool::GetIDSF(
        int RunNumber, double pt, double eta, double jet_dr,
        double jet_pt, bool afii) {

    if (!HaveSF(pt, eta, jet_dr, jet_pt))
        return ValUncert(-999, 999);

    if (jet_dr < 0.2)
        jet_dr = 0.2;

    if (!fIDSFTool) {
        cerr << "TTBarElectronSFTool: please initialize an " <<
            "electron identification SF tool and set it by " <<
            "calling SetIDSFTool()" << endl;
        return ValUncert(-999, 999);
    }


    ValUncert sf;

    // ID efficiency SF parameterized in dr(jet) bins out to 0.4.
    if (jet_dr < 0.4) {
        sf = ContentAndError(fPtIDSFHist, pt, true);

        // no eta overflow
        sf = VUMult(sf,
                ContentAndError(fEtaIDSFHist, eta, false));

    } else {
        PATCore::ParticleDataType::DataType dt;
        if (afii)
            dt = PATCore::ParticleDataType::Fast;
        else
            dt = PATCore::ParticleDataType::Full;

        Root::TResult r = fRecoSFTool->calculate(
                dt, RunNumber, eta, pt);

        sf = ValUncert(r.getScaleFactor(), r.getTotalUncertainty());
    }

    return sf;
}

ValUncert TTBarElectronSFTool::GetTrigIsoSF(
        int RunNumber, double pt, double eta, double jet_dr,
        double jet_pt, bool afii) {

    // Trig+Iso efficiency parameterized in dr(jet) bins out to 0.6.
    ValUncert sf;

    if (!HaveSF(pt, eta, jet_dr, jet_pt))
        return ValUncert(-999, 999);

    if (jet_dr < 0.2)
        jet_dr = 0.2;

    if (jet_dr < 0.4) {
        sf = ContentAndError(fPtTrigIsoSFHist, pt, true);

        // no eta overflow
        sf = VUMult(sf,
                ContentAndError(fEtaTrigIsoSFHist, eta, false));

    } else {
        sf = ContentAndError(fFarJetPtTrigIsoSFHist, pt, true);

        sf = VUMult(sf,
                ContentAndError(fFarJetEtaTrigIsoSFHist, eta, false));
    }

    return sf;
}

ValUncert TTBarElectronSFTool::GetSF(int RunNumber,
        double pt, double eta, double jet_dr, double jet_pt, bool afii) {

    ValUncert sf;

    // have to check, otherwise might get two -1s which multiply to 1.
    if (!HaveSF(pt, eta, jet_dr, jet_pt))
        return ValUncert(-999, 999);

    sf = GetRecoSF(RunNumber, pt, eta, afii);

    sf = VUMult(sf, GetIDSF(RunNumber, pt, eta, jet_dr, jet_pt, afii));

    sf = VUMult(sf, GetTrigIsoSF(RunNumber, pt, eta, jet_dr, jet_pt, afii));

    return sf;
}


// operations on ValUncerts

double QuadSum(double x, double y) {
    return sqrt(x*x + y*y);
}

ValUncert VUAdd(const ValUncert &u, const ValUncert &v) {
    return ValUncert(u.first+v.first, QuadSum(u.second, v.second));
}

ValUncert VUNeg(const ValUncert &u) {
    return ValUncert(-u.first, u.second);
}

ValUncert VUSub(const ValUncert &u, const ValUncert &v) {
    return VUAdd(u, VUNeg(v));
}

ValUncert VUMult(const ValUncert &u, const ValUncert &v) {
    return ValUncert(u.first*v.first,
            QuadSum(u.first*v.second, v.first*u.second));
}

ValUncert VUDiv(const ValUncert &u, const ValUncert &v) {
    return ValUncert(u.first/v.first,
            QuadSum(u.first*v.second/v.first/v.first, u.second/v.first));
}

ValUncert VUConst(const double c) {
    return ValUncert(c, 0.0);
}
