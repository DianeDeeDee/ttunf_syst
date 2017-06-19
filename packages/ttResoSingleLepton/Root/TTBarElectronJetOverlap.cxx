#include "ttResoSingleLepton/TTBarElectronJetOverlap.h"

#include <iostream>

#include "TMath.h"

using namespace std;

void TTBarElectronJetOverlap::LoadJets(
        const vector<float> &jet_pt,
        const vector<float> &jet_eta,
        const vector<float> &jet_phi,
        const vector<float> &jet_E,
        const vector<float> &jet_jvf,
        const vector<float> &jet_sumPtTrk_pv0_500MeV,
        const vector<vector<int> > &jet_TrackAssoc_index) {

    fAnalyzed = false;

    size_t s = jet_pt.size();
    if (s != jet_eta.size() ||
            s != jet_phi.size() ||
            s != jet_E.size() ||
            s != jet_jvf.size() ||
            s != jet_sumPtTrk_pv0_500MeV.size() ||
            s != jet_TrackAssoc_index.size()) {

        cout << "ElectronJetOverlap Error: not all vectors in " <<
            "LoadJets() are the same size." << endl;

        return;
    }

    fJetJVFs = jet_jvf;
    fJetD3PDTrkPtPVSums = jet_sumPtTrk_pv0_500MeV;

    fJetTLVs.resize(s);
    fJetD3PDTrkPtSums.resize(s);
    fJetTrkAssocIndices.resize(s);
    for (size_t i = 0; i < s; i++) {
        fJetTLVs[i].SetPtEtaPhiE(jet_pt[i], jet_eta[i], jet_phi[i], jet_E[i]);
        fJetD3PDTrkPtSums[i] = fJetJVFs[i] ? fJetD3PDTrkPtPVSums[i]/fJetJVFs[i] : -1;
        fJetTrkAssocIndices[i] = set<int>();
        for (size_t j = 0; j < jet_TrackAssoc_index[i].size(); j++)
            fJetTrkAssocIndices[i].insert(jet_TrackAssoc_index[i][j]);
    }

    return;
}

void TTBarElectronJetOverlap::LoadSelectedElectrons(
        const vector<float> &el_cl_E,
        const vector<float> &el_cl_eta,
        const vector<float> &el_cl_phi,
        const vector<float> &el_trk_eta,
        const vector<float> &el_trk_phi,
        const vector<int> &el_GSF_trk_index) {

    fAnalyzed = false;

    size_t s = el_cl_E.size();
    if (s != el_cl_eta.size() ||
            s != el_cl_phi.size() ||
            s != el_trk_eta.size() ||
            s != el_trk_phi.size() ||
            s != el_GSF_trk_index.size()) {

        cout << "ElectronJetOverlap Error: not all vectors in " <<
            "LoadSelectedElectrons() are the same size." << endl;

        return;
    }

    fElGSFTrkIndices = el_GSF_trk_index;

    fElClTLVs.resize(s);
    fElTLVs.resize(s);

    for (size_t i = 0; i < s; i++) {
        fElClTLVs[i].SetPtEtaPhiM(el_cl_E[i]/TMath::CosH(el_cl_eta[i]),
                el_cl_eta[i], el_cl_phi[i], 0.511);
        fElTLVs[i].SetPtEtaPhiM(el_cl_E[i]/TMath::CosH(el_trk_eta[i]),
                    el_trk_eta[i], el_trk_phi[i], 0.511);
    }

    return;
}

void TTBarElectronJetOverlap::LoadTracks(
        const vector<int> &GSF_trk_trk_index,
        const vector<float> &trk_pt,
        const vector<float> &trk_eta,
        const vector<float> &trk_phi_wrtPV,
        const vector<float> &trk_theta_wrtPV,
        const vector<float> &trk_z0_wrtPV,
        const vector<float> &trk_d0_wrtPV) {

    fAnalyzed = false;

    size_t s = trk_pt.size();
    if (s != trk_eta.size() ||
            s != trk_phi_wrtPV.size() ||
            s != trk_theta_wrtPV.size() ||
            s != trk_z0_wrtPV.size() ||
            s != trk_d0_wrtPV.size()) {

        cout << "ElectronJetOverlap Error: not all vectors in " <<
            "LoadTracks() are the same size." << endl;

        return;
    }

    fTrkTLVs.resize(s);
    for (size_t i = 0; i < s; i++)
        fTrkTLVs[i].SetPtEtaPhiM(trk_pt[i], trk_eta[i], trk_phi_wrtPV[i], 0);


    fTrkThetas = trk_theta_wrtPV;
    fTrkZ0s = trk_z0_wrtPV;
    fTrkD0s = trk_d0_wrtPV;
    fGSFTrkIndices = GSF_trk_trk_index;

    return;
}

void TTBarElectronJetOverlap::FindAssocEls() {
    size_t nJets = fJetTLVs.size();
    fJetAssocElCls = vector<set<int> >(nJets);

    size_t nEls = fElTLVs.size();
    fElClAssocJet = vector<int>(nEls, -1);


    // find the associated electron clusters in each jet
    double drmin, dr;
    int drmin_idx;
    TLorentzVector elcl, jet;
    for (size_t iEl = 0; iEl < nEls; iEl++) {
        elcl = fElClTLVs[iEl];

        drmin = 0.4;
        drmin_idx = -1;
        for (size_t iJet = 0; iJet < nJets; iJet++) {
            jet = fJetTLVs[iJet];
            dr = elcl.DeltaR(jet);

            // attempt to match the jet to the electron.
            if (dr < drmin) {
                drmin = dr;
                drmin_idx = iJet;
            }
        }

        if (drmin_idx < 0)
            continue;

        fJetAssocElCls[drmin_idx].insert(iEl);
        fElClAssocJet[iEl] = drmin_idx;
    }

    return;
}

void TTBarElectronJetOverlap::SubtractEls() {
    size_t nJets = fJetTLVs.size();
    int ElIdx;
    for (size_t iJet = 0; iJet < nJets; iJet++) {
        for (set<int>::iterator iEl = fJetAssocElCls[iJet].begin();
                iEl != fJetAssocElCls[iJet].end(); iEl++) {
            ElIdx = *iEl;

            fJetTLVs[iJet] -= fElTLVs[ElIdx];
        }
    }

    return;
}

void TTBarElectronJetOverlap::FindGoodObjects() {
    size_t nEls = fElTLVs.size();
    fGoodEls = vector<bool>(nEls, true);

    size_t nJets = fJetTLVs.size();
    vector<TLorentzVector> TmpJetTLVs = fJetTLVs;

    TLorentzVector el, jet;
    for (size_t iEl = 0; iEl < nEls; iEl++) {
        el = fElTLVs[iEl];

        if (fDebug) {
            cout << "El # " << iEl <<
                " Pt Eta Phi: " <<
                el.Pt() << " " <<
                el.Eta() << " " <<
                el.Phi() << endl;
        }

        for (size_t iJet = 0; iJet < nJets; iJet++) {
            jet = fJetTLVs[iJet];

            if (fDebug) {
                cout << "  Jet # " << iJet <<
                    " Pt Eta Phi M: " <<
                    jet.Pt() << " " <<
                    jet.Eta() << " " <<
                    jet.Phi() << " " <<
                    jet.M() << endl;
            }

            if (jet.Pt() < 25000)
                continue;

            // if the electron is too close to a jet...
            if (jet.DeltaR(el) < 0.2) {

                if (fDebug)
                    cout << "    El too close to jet. Removing." << endl;

                // remove from good electrons list.
                fGoodEls[iEl] = false;

                if (fElClAssocJet[iEl] >= 0) {
                    // add electron 4-vector back to jet.
                    TmpJetTLVs[fElClAssocJet[iEl]] += el; 

                    // remove this electron from the association.
                    fJetAssocElCls[iJet].erase(iEl);
                }
            }
        }
    }

    fJetTLVs = TmpJetTLVs;

    fGoodJets.resize(nJets);
    for (size_t iJet = 0; iJet < nJets; iJet++)
        fGoodJets[iJet] = fJetTLVs[iJet].Pt() > 25000;

    return;
}

void TTBarElectronJetOverlap::RecalcJVF() {

    size_t nJets = fJetTLVs.size();
    int ElIdx;
    int TrkIdx;
    TLorentzVector trk;

    for (size_t iJet = 0; iJet < nJets; iJet++) {

        // if jvf == 0 or -1, no need to recalulate...
        if (!fJetJVFs[iJet] ||
                fJetJVFs[iJet] < 0)
            continue;

        for (set<int>::iterator iEl = fJetAssocElCls[iJet].begin();
                iEl != fJetAssocElCls[iJet].end(); iEl++) {

            ElIdx = *iEl;
            TrkIdx = fGSFTrkIndices[fElGSFTrkIndices[ElIdx]];

            // not in the associated tracks.
            if (fJetTrkAssocIndices[iJet].find(TrkIdx) ==
                    fJetTrkAssocIndices[iJet].end())
                continue;

            // recompute JVF.
            trk = fTrkTLVs[TrkIdx];
            fJetD3PDTrkPtSums[iJet] -= trk.Pt();

            if (TMath::Abs(TMath::Sin(fTrkThetas[TrkIdx])*fTrkZ0s[TrkIdx]) < 1.0 &&
                    TMath::Abs(fTrkD0s[TrkIdx]) < 1.0)
                fJetD3PDTrkPtPVSums[iJet] -= trk.Pt();
        }

        // we subtracted too much for some reason?
        if (fJetD3PDTrkPtSums[iJet] < 0 ||
                fJetD3PDTrkPtPVSums[iJet] < 0)
            fJetJVFs[iJet] = 0;
        // no tracks associated with the jet.
        else if (fJetD3PDTrkPtSums[iJet] == 0)
            fJetJVFs[iJet] = -1;
        // all good. recalculate JVF.
        else
            fJetJVFs[iJet] = fJetD3PDTrkPtPVSums[iJet]/fJetD3PDTrkPtSums[iJet];
    }

    return;
}

void TTBarElectronJetOverlap::AnalyzeEvent() {
    FindAssocEls();
    SubtractEls();
    FindGoodObjects();
    RecalcJVF();

    fAnalyzed = true;

    return;
}

vector<float> TTBarElectronJetOverlap::NewJetJVFs() {
    if (!fAnalyzed)
        AnalyzeEvent();

    return fJetJVFs;
}

vector<TLorentzVector> TTBarElectronJetOverlap::NewJetTLVs() {
    if (!fAnalyzed)
        AnalyzeEvent();

    return fJetTLVs;
}

vector<set<int> > TTBarElectronJetOverlap::JetAssocEls() {
    if (!fAnalyzed)
        AnalyzeEvent();

    return fJetAssocElCls;
}

vector<bool> TTBarElectronJetOverlap::GoodEls() {
    if (!fAnalyzed)
        AnalyzeEvent();

    return fGoodEls;
}
