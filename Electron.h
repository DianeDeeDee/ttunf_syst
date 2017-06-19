#ifndef ELECTRON_H
#define ELECTRON_H

#include "MObject.h"
#include "TLorentzVector.h"

class Electron : public MObject {
  public:
    Electron();
    Electron(const TLorentzVector &v);
    virtual ~Electron();

    void setMI(float iso);
    float mi() const;
    void setTightPP(bool isTightPP);
    bool isTightPP() const;

    void setMediumPP(bool isMediumPP);
    bool isMediumPP() const;

    void setLoosePP(bool isLoosePP);
    bool isLoosePP() const;

    const TLorentzVector &caloMom() const;
    TLorentzVector &caloMom();
    const TLorentzVector &trkMom() const;
    TLorentzVector &trkMom();
    const float z0() const;
    float &z0();
    const int author() const;
    int &author();

    int &nSiHits();
    const int nSiHits() const;

    int &oq();
    const int oq() const;

    int &isEM();
    const int isEM() const;

    int &GSF_trk_index();
    const int GSF_trk_index() const;

    bool pass() const;
    bool passLoose(bool looser = false) const;

    bool &passOR();
  protected:
    float m_mi;
    bool m_isTightPP;
    bool m_isMediumPP;
    bool m_isLoosePP;
    TLorentzVector m_mom_calo;
    TLorentzVector m_mom_trk;
    float m_z0;
    int m_author;

    int m_nSiHits;
    int m_oq;
    int m_isEM;
    int m_GSF_trk_index;

    bool m_passOR;
};

#endif
