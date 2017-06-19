// weak corrections for powheg ttbar MC
// using parametrisations based on the HATHOR
// implementation of calculations from J. Kuehn, P. Uwer
// Sebastian.Fleischmann -at- cern.ch

#include <iostream>

#include "ttResoSingleLepton/WeakCorrScaleFactorParam.h"

#include "TLorentzVector.h"
#include "TFile.h"

WeakCorr::WeakCorrScaleFactorParam::WeakCorrScaleFactorParam(TString mapfile) :
    m_fuu(0),
    m_fdd(0),
    m_fgg(0),
    m_mt(172.5)
{
    if (!init(mapfile)) {
        std::cout << "ERROR: Could not initialise histograms!" << std::endl;
        delete m_fuu; m_fuu = 0;
        delete m_fdd; m_fdd = 0;
        delete m_fgg; m_fgg = 0;
    }
}

WeakCorr::WeakCorrScaleFactorParam::~WeakCorrScaleFactorParam() {
    delete m_fuu;
    delete m_fdd;
    delete m_fgg;
}

bool WeakCorr::WeakCorrScaleFactorParam::init(TString mapfile) {
    TFile infile(mapfile.Data(), "READ");
    TH2F* hist = dynamic_cast<TH2F*>(infile.Get("Fuu"));
    if (!hist) return false;
    m_fuu = (TH2F*) hist->Clone("myfuu");
    m_fuu->SetDirectory(0);
    hist = dynamic_cast<TH2F*>(infile.Get("Fdd"));
    if (!hist) return false;
    m_fdd = (TH2F*) hist->Clone("myfdd");
    m_fdd->SetDirectory(0);
    hist = dynamic_cast<TH2F*>(infile.Get("Fgg"));
    if (!hist) return false;
    m_fgg = (TH2F*) hist->Clone("myfgg");
    m_fgg->SetDirectory(0);
    return true;
}

WeakCorr::ScaleFactor WeakCorr::WeakCorrScaleFactorParam::getScaleFactor(const double& shat, const double& z, const InitialStateType& type) {
    ScaleFactor sf;
    sf.nominal = getWeight(shat, z, type);
    // 10% uncertainty on correction (F being the relative correction):
    // nominal = 1 + F
    // up   = 1 + F + 0.1*F = nominal + 0.1*(nominal - 1) = 1.1*nominal - 0.1
    // down = 1 + F - 0.1*F = nominal - 0.1*(nominal - 1) = 0.9*nominal + 0.1
    // up/down mean stronger/less correction, i.e. down > nominal in case of negative correction
    sf.up = 1.1 * sf.nominal - 0.1;
    sf.down = 0.9 * sf.nominal + 0.1;
    return sf;
}

WeakCorr::ScaleFactor WeakCorr::WeakCorrScaleFactorParam::getScaleFactor(
        const std::vector<int>& mc_status,
        const std::vector<int>& mc_pdgId,
        const std::vector<float>& mc_pt,
        const std::vector<float>& mc_eta,
        const std::vector<float>& mc_phi,
        const std::vector<float>& mc_m,
        const std::vector<std::vector<int> >& mc_parent_index,
        const std::vector<std::vector<int> >& mc_child_index
    ) {
    ScaleFactor sf;
    sf.nominal = getWeight(mc_status, mc_pdgId, mc_pt, mc_eta, mc_phi, mc_m, mc_parent_index, mc_child_index);
    // 10% uncertainty on correction (F being the relative correction):
    // nominal = 1 + F
    // up   = 1 + F + 0.1*F = nominal + 0.1*(nominal - 1) = 1.1*nominal - 0.1
    // down = 1 + F - 0.1*F = nominal - 0.1*(nominal - 1) = 0.9*nominal + 0.1
    // up/down mean stronger/less correction, i.e. down > nominal in case of negative correction
    sf.up = 1.1 * sf.nominal - 0.1;
    sf.down = 0.9 * sf.nominal + 0.1;
    return sf;    
}

double WeakCorr::WeakCorrScaleFactorParam::getWeight(const double& shat, const double& z, const InitialStateType& type) {
    if (!m_fuu || !m_fdd  || !m_fgg) {
        return -1.;
    }
    double weight = 1.;
    const double mtt = sqrt(shat);
    if (mtt < m_mt*2.) {
        return weight;
    }
    if (type == GG) {
//         weight = m_fgg->Interpolate(shat, z);
        weight = m_fgg->Interpolate(mtt, z);
    } else if (type == UU) {
        weight = m_fuu->Interpolate(mtt, z);
//         weight = m_fuu->Interpolate(shat, z);
    } else if (type == DD) {
        weight = m_fdd->Interpolate(mtt, z);
//         weight = m_fdd->Interpolate(shat, z);
    } else {
        std::cout << "ERROR: Wrong initial state type given" << std::endl;
    }
    
    return weight;
}

double WeakCorr::WeakCorrScaleFactorParam::getWeight(
        const std::vector<int>& mc_status,
        const std::vector<int>& mc_pdgId,
        const std::vector<float>& mc_pt,
        const std::vector<float>& mc_eta,
        const std::vector<float>& mc_phi,
        const std::vector<float>& mc_m,
        const std::vector<std::vector<int> >& mc_parent_index,
        const std::vector<std::vector<int> >& mc_child_index
    ) {
    double weight = 1.;
    
    int topIndex = -1;
    int antitopIndex = -1;
    
    for (unsigned int partIndex = 0;  partIndex < mc_pt.size(); partIndex++) {
        // take particles from original Les Houches Event record of powheg
        if (mc_status.at(partIndex) != 3) {
            continue;
        }
        if (mc_pdgId.at(partIndex) == 6) {
            topIndex = partIndex;
            if (antitopIndex > -1) {
                break;
            }
        }
        if (mc_pdgId.at(partIndex) == -6) {
            antitopIndex = partIndex;
            if (topIndex > -1) {
                break;
            }
        }
    }
    if ((topIndex < 0) || (antitopIndex < 0)) {
        std::cout << "ERROR: Could not find MC tops!" << std::endl;
        return weight;
    }
    ////////
    // use direct parent of top and anti-top:
    if (mc_parent_index.at(topIndex).size() < 2) {
        std::cout << "ERROR: Could not get top parents!" << std::endl;
        return weight;        
    }
    const int initialPartons0 = mc_parent_index.at(topIndex).at(0);
    const int initialPartons1 = mc_parent_index.at(topIndex).at(1);
    if ((initialPartons0 < 0) || (initialPartons1 < 0)) {
        std::cout << "ERROR: Could not get top parents!" << std::endl;
        return weight;        
    }
    const int pdg1 = abs(mc_pdgId.at(initialPartons0));
    const int pdg2 = abs(mc_pdgId.at(initialPartons1));
    int incomingQuark = -1;
    
    // search outgoing parton (ttg/ttq):
    int outgoingParton = -1;
    std::vector<int>::const_iterator child = mc_child_index.at(initialPartons0).begin();
    for ( ; child != mc_child_index.at(initialPartons0).end(); child++) {
        const int pdgID = abs(mc_pdgId.at(*child));
        if (pdgID != 6) {
            outgoingParton = (*child);
            if (!(pdgID == 1 || pdgID == 3 || pdgID == 5 || pdgID == 2 || pdgID == 4 || pdgID == 21)) {
                std::cout << "ERROR: Could not determine outgoing parton!" << std::endl;
                return weight;
            }
            break;
        }
    }
    if (outgoingParton < 0) {
        std::cout << "ERROR: Could not determine outgoing parton!" << std::endl;
        return weight;
    }

    InitialStateType eventType = Undefined;
    InitialStateType factorToUse = Undefined;
    if (pdg1 == pdg2) {
        switch (pdg1) {
            case 21:
                eventType = GG;
                break;
            case 1:
            case 3:
            case 5:
                eventType = DD;
                break;
            case 2:
            case 4:
                eventType = UU;
                break;
        }
        factorToUse = eventType;
    } else {
        int pdgQ = pdg1;
        if (pdg1 == 21) {
            pdgQ = pdg2;
            incomingQuark = initialPartons1;
        } else if (pdg2 == 21) {
            incomingQuark = initialPartons0;
        } else {
            std::cout << "ERROR: Could not determine type of initial state!" << std::endl;
            return weight;
        }
        switch (pdgQ) {
            case 1:
            case 3:
            case 5:
                eventType = GD;
                break;
            case 2:
            case 4:
                eventType = GU;
                break;
            default:
                std::cout << "ERROR: Could not determine type of initial state!" << std::endl;
                return weight;
        }
        TLorentzVector incomingQuarkMomentum;
        incomingQuarkMomentum.SetPtEtaPhiM( mc_pt.at(incomingQuark)*0.001,
                                            mc_eta.at(incomingQuark),
                                            mc_phi.at(incomingQuark),
                                            mc_m.at(incomingQuark)*0.001  );
        TLorentzVector hadronMomentum(0., 0., -8000., 8000.);
        if (incomingQuarkMomentum.Z() > 0.) {
            hadronMomentum.SetZ(8000.);
        }
        TLorentzVector outgoingPartonMomentum;
        outgoingPartonMomentum.SetPtEtaPhiM( mc_pt.at(outgoingParton)*0.001,
                                             mc_eta.at(outgoingParton),
                                             mc_phi.at(outgoingParton),
                                             mc_m.at(outgoingParton)*0.001  );
        const double thetaQ = hadronMomentum.Angle(outgoingPartonMomentum.Vect());
        if (thetaQ < (TMath::Pi()*0.5)) {
            factorToUse = GG;
        } else {
            if (eventType == GD) {
                factorToUse = DD;
            } else if (eventType == GU) {
                factorToUse = UU;
            }
        }
    }
    TLorentzVector topMomentum;
    topMomentum.SetPtEtaPhiM( mc_pt.at(topIndex)*0.001,
                              mc_eta.at(topIndex),
                              mc_phi.at(topIndex),
                              mc_m.at(topIndex)*0.001  );
    TLorentzVector antitopMomentum;
    antitopMomentum.SetPtEtaPhiM( mc_pt.at(antitopIndex)*0.001,
                                  mc_eta.at(antitopIndex),
                                  mc_phi.at(antitopIndex),
                                  mc_m.at(antitopIndex)*0.001   );

    TLorentzVector ttbarSystem = topMomentum + antitopMomentum;
    const double shat = ttbarSystem.M2();
    
    TLorentzVector beam1(0., 0.,  8000., 8000.);
    TLorentzVector beam2(0., 0., -8000., 8000.);
    TVector3 boostVec = ttbarSystem.BoostVector();
    boostVec *= -1.;
    
    topMomentum.Boost(boostVec);
    beam1.Boost(boostVec);
    beam2.Boost(boostVec);
    
    TVector3 scatteringAxis = (beam1 - beam2).Vect().Unit();
    
    const double z = topMomentum.Vect().Unit() * scatteringAxis;
    
    weight = getWeight(shat, z, factorToUse);
    
    return weight;
}
