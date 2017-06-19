#ifndef __TTBARELECTRONJETOVERLAP_H__
#define __TTBARELECTRONJETOVERLAP_H__

#include <vector>
#include <set>
#include <utility>
#include "TLorentzVector.h"

class TTBarElectronJetOverlap {
    private:
        std::vector<TLorentzVector> fJetTLVs;
        std::vector<float> fJetJVFs;
        std::vector<float> fJetD3PDTrkPtPVSums;
        std::vector<float> fJetD3PDTrkPtSums;
        std::vector<std::set<int> > fJetTrkAssocIndices;

        std::vector<bool> fGoodJets;
        std::vector<std::set<int> > fJetAssocElCls;

        std::vector<TLorentzVector> fElTLVs;
        std::vector<TLorentzVector> fElClTLVs;

        std::vector<bool> fGoodEls;
        std::vector<int> fElClAssocJet;
        std::vector<int> fElTrkIdx;
        std::vector<int> fElGSFTrkIndices;

        std::vector<int> fGSFTrkIndices;
        std::vector<TLorentzVector> fTrkTLVs;
        std::vector<float> fTrkThetas;
        std::vector<float> fTrkZ0s;
        std::vector<float> fTrkD0s;

        bool fAnalyzed;
        bool fDebug;

        void FindElTrks();
        void FindAssocEls();
        void SubtractEls();
        void FindGoodObjects();
        void RecalcJVF();

        void AnalyzeEvent();

    public:
        TTBarElectronJetOverlap() {
            fDebug = false;
            fAnalyzed = false;
            return;
        }

        ~TTBarElectronJetOverlap() { }

        // sets the debug level (info will print if set true).
        void SetDebug(bool db) {
            fDebug = db;
            return;
        }

        bool GetDebug() {
            return fDebug;
        }

        // load all anti-kt 0.4 LCTopo jets in the event.
        // these variables should correspond to the *corrected* jet
        // quantities. Every jet in the D3PD should be passed.
        void LoadJets(
                const std::vector<float> &jet_pt,
                const std::vector<float> &jet_eta,
                const std::vector<float> &jet_phi,
                const std::vector<float> &jet_E,
                const std::vector<float> &jet_jvtxf,
                const std::vector<float> &jet_sumPtTrk_pv0_500MeV,
                const std::vector<std::vector<int> > &jet_TrackAssoc_index);

        // load selected electrons.
        // these variables should only be filled for /selected/
        // electrons (i.e. pass ID and isolation cuts).
        void LoadSelectedElectrons(
                const std::vector<float> &el_cl_E,
                const std::vector<float> &el_cl_eta,
                const std::vector<float> &el_cl_phi,
                const std::vector<float> &el_trk_eta,
                const std::vector<float> &el_trk_phi,
                const std::vector<int> &el_GSF_trk_index);

        // load all tracks from the event.
        // these variables should come directly from the D3PD.
        void LoadTracks(
                const std::vector<int> &GSF_trk_trk_index,
                const std::vector<float> &trk_pt,
                const std::vector<float> &trk_eta,
                const std::vector<float> &trk_phi_wrtPV,
                const std::vector<float> &trk_theta_wrtPV,
                const std::vector<float> &trk_z0_wrtPV,
                const std::vector<float> &trk_d0_wrtPV);

        // return the recalculated JVFs for all jets in the event.
        std::vector<float> NewJetJVFs();

        // return the recalculated TLorentzVectors for all jets in the event.
        std::vector<TLorentzVector> NewJetTLVs();

        // return the indices for electrons whose 4-vectors were
        // subtraced from the given jet.
        std::vector<std::set<int> > JetAssocEls();

        // return a boolean for each electron which represents whether or
        // not it should be considered a selected electron (e.g. it
        // wasn't removed by the overlap procedure).
        std::vector<bool> GoodEls();
};

#endif
