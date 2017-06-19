#ifndef WEAKCORRSCALEFACTORPARAM_H
#define WEAKCORRSCALEFACTORPARAM_H
// weak corrections for powheg ttbar MC
// using parametrisations based on the HATHOR
// implementation of calculations from J. Kuehn, P. Uwer
// Sebastian.Fleischmann -at- cern.ch

// if using this code, please cite:
// * J. H. Kuehn, A. Scharf and P. Uwer,
//    "Electroweak corrections to top-quark pair production in quark-antiquark annihilation",
//     Eur.Phys.J. C45, 139 (2006), hep-ph/0508092.
// * J. H. Kuehn, A. Scharf and P. Uwer, 
//     "Electroweak effects in top-quark pair production at hadron colliders”, 
//     Eur.Phys.J. C51, 37 (2007), hep-ph/0610335.
// * J. Kuehn, A. Scharf and P. Uwer, 
//     "Weak Interactions in Top-Quark Pair Production at Hadron Colliders: An Update”,
//     arxiv:1305.5773.
// * G. van Oldenborgh, "FF: A Package to evaluate one loop Feynman diagrams",
//     Comput.Phys.Commun. 66, 1 (1991).
// * M. Luscher,
//     "A Portable high quality random number generator for lattice field theory simulations",
//     Comput.Phys.Commun. 79, 100 (1994), hep-lat/9309020.
// * G. P. Lepage, "A New Algorithm for Adaptive Multidimensional Integration",
//     J.Comput.Phys. 27, 192 (1978).
// * G. P. Lepage, "VEGAS: AN ADAPTIVE MULTIDIMENSIONAL INTEGRATION PROGRAM".


#include "TH2F.h"
#include "ttResoSingleLepton/ScaleFactor.h"

namespace WeakCorr {

class WeakCorrScaleFactorParam {
  
public:
    enum InitialStateType {
        Undefined,
        UU,
        DD,
        GG,
        GD,
        GU
    };
    
    WeakCorrScaleFactorParam(TString mapfile="ewParam-mtt.root");
    ~WeakCorrScaleFactorParam();
    ScaleFactor getScaleFactor(const double& shat, const double& z, const InitialStateType& type);
    ScaleFactor getScaleFactor(
        const std::vector<int>& mc_status,
        const std::vector<int>& mc_pdgId,
        const std::vector<float>& mc_pt,
        const std::vector<float>& mc_eta,
        const std::vector<float>& mc_phi,
        const std::vector<float>& mc_m,
        const std::vector<std::vector<int> >& mc_parent_index,
        const std::vector<std::vector<int> >& mc_child_index
    );
    double getWeight(const double& shat, const double& z, const InitialStateType& type);
    double getWeight(
        const std::vector<int>& mc_status,
        const std::vector<int>& mc_pdgId,
        const std::vector<float>& mc_pt,
        const std::vector<float>& mc_eta,
        const std::vector<float>& mc_phi,
        const std::vector<float>& mc_m,
        const std::vector<std::vector<int> >& mc_parent_index,
        const std::vector<std::vector<int> >& mc_child_index
    );
private:
    bool init(TString mapfile);
    TH2F* m_fuu;
    TH2F* m_fdd;
    TH2F* m_fgg;
    double m_mt;
};
} // namespace WeakCorr

#endif //WEAKCORRSCALEFACTORPARAM_H
