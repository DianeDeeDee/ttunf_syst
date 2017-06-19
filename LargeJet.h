#ifndef LARGEJET_H
#define LARGEJET_H

#include "MObject.h"
#include "TLorentzVector.h"
#include <fastjet/PseudoJet.hh> 
#include <fastjet/ClusterSequence.hh> 

class LargeJet : public MObject {
  public:

    class Subjet {
      public:
      TLorentzVector v;
      float lcpt;
    };

    LargeJet();
    LargeJet(const TLorentzVector &v);
    LargeJet(const LargeJet &l);
    virtual ~LargeJet();

    int &trueFlavour();
    const int trueFlavour() const;

    const std::vector<Subjet> &subjet() const;
    std::vector<Subjet> &subjet();

    const std::vector<Subjet> &subjetJER() const;
    std::vector<Subjet> &subjetJER();

    const std::vector<TLorentzVector> &subjet_area() const;
    std::vector<TLorentzVector> &subjet_area();

    bool pass() const;
    bool passLoose() const;

    double &detE();
    double &detEta();
    double &detPhi();
    double &detM();
    double &Ax();
    double &Ay();
    double &Az();
    double &Ae();

    TLorentzVector &bdrs();
    const TLorentzVector bdrs() const;

    double &trimmed_tau1();
    const double trimmed_tau1() const;
    double &trimmed_tau2();
    const double trimmed_tau2() const;
    double &trimmed_am_tau1();
    const double trimmed_am_tau1() const;
    double &trimmed_am_tau2();
    const double trimmed_am_tau2() const;


    double &ug_tau1();
    const double ug_tau1() const;
    double &ug_tau2();
    const double ug_tau2() const;
    double &ug_am_tau1();
    const double ug_am_tau1() const;
    double &ug_am_tau2();
    const double ug_am_tau2() const;



    const double bdrs_m() const;
    const double bdrs_pt() const;
    const double bdrs_eta() const;
    const double bdrs_phi() const;

    double &ug_bdrs_tau1();
    const double ug_bdrs_tau1() const;
    double &ug_bdrs_tau2();
    const double ug_bdrs_tau2() const;
    /*
    double &trimmed_bdrs_tau1();
    const double trimmed_bdrs_tau1() const;
    double &trimmed_bdrs_tau2();
    const double trimmed_bdrs_tau2() const;
    */


    double &split12();
    const double split12() const;

    double &htt();
    const double htt() const;

  protected:
    std::vector<Subjet> m_subjet;
    std::vector<Subjet> m_subjetJER;
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
    double m_trimmed_tau1;
    double m_trimmed_tau2;
    double m_ug_tau1;
    double m_ug_tau2;
    TLorentzVector m_bdrs;
    //   double m_trimmed_bdrs_tau1;
    //double m_trimmed_bdrs_tau2;
    double m_ug_bdrs_tau1;
    double m_ug_bdrs_tau2;
    double m_trimmed_am_tau1;
    double m_trimmed_am_tau2;
    double m_ug_am_tau1;
    double m_ug_am_tau2;

    int m_trueFlavour;

    double m_htt;
};

#endif
