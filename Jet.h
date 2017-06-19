#ifndef JET_H
#define JET_H

#include "MObject.h"
#include "TLorentzVector.h"

class Jet : public MObject {
  public:
    Jet();
    Jet(const TLorentzVector &v);
    virtual ~Jet();

    int &trueFlavour();
    const int trueFlavour() const;

    const float mv1() const;
    float &mv1();

    const float jvf() const;
    float &jvf();

    const float sumPtTrk() const;
    float &sumPtTrk();

    const bool isBadLoose() const;
    bool &isBadLoose();
   
    const bool isBadLooseMinus() const;
    bool &isBadLooseMinus();
    
    const std::vector<int> &TrackAssoc_index() const;
    std::vector<int> &TrackAssoc_index();
const float &jfit_deltaEta() const;
	    float &jfit_deltaEta();
	
	    const float &jfit_deltaPhi() const;
	    float &jfit_deltaPhi();

    bool pass() const;
    bool btag() const;

    bool passBadLoose() const;
    bool passBadLooseMinus() const;
    double &detE();
    double &detEta();
    double &detPhi();
    double &detM();
    double &Ax();
    double &Ay();
    double &Az();
    double &Ae();

  protected:
    int m_trueflavour;
    float m_mv1;

    float m_jvf;

    bool m_isBadLoose;
    bool m_isBadLooseMinus;
float m_jfit_deltaEta;
	    float m_jfit_deltaPhi;
    double m_detE;
    double m_detEta;
    double m_detPhi;
    double m_detM;
    double m_Ax;
    double m_Ay;
    double m_Az;
    double m_Ae;

    std::vector<int> m_TrackAssoc_index;
    float m_sumPtTrk;
};

#endif
