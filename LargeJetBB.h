#ifndef LARGEJETBB_H
#define LARGEJETBB_H

#include "MObject.h"
#include "TLorentzVector.h"
#include <fastjet/PseudoJet.hh> 
#include <fastjet/ClusterSequence.hh> 

class LargeJetBB : public MObject {
  public:
    LargeJetBB();
    LargeJetBB(const TLorentzVector &v);
    LargeJetBB(const LargeJetBB &l);
    virtual ~LargeJetBB();

    const std::vector<TLorentzVector> &subjet() const;
    std::vector<TLorentzVector> &subjet();

    const std::vector<TLorentzVector> &subjet_area() const;
    std::vector<TLorentzVector> &subjet_area();

    bool pass() const;

    double &detE();
    double &detEta();
    double &detPhi();
    double &detM();
    double &Ax();
    double &Ay();
    double &Az();
    double &Ae();
    double &ug_tau1();
    double &ug_tau2();
    double &bdrs_m();
    double &bdrs_pt();
    double &bdrs_tau1();
    double &bdrs_tau2();


  protected:
    std::vector<TLorentzVector> m_subjet;
    std::vector<TLorentzVector> m_subjet_area;

    double m_detE;
    double m_detEta;
    double m_detPhi;
    double m_detM;
    double m_Ax;
    double m_Ay;
    double m_Az;
    double m_Ae;
    double m_split12;
    double m_ug_tau1;
    double m_ug_tau2;
    double m_bdrs_m;
    double m_bdrs_pt;
    double m_bdrs_tau1;
    double m_bdrs_tau2;
};

#endif
