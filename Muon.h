#ifndef MUON_H
#define MUON_H

#include "MObject.h"
#include "TLorentzVector.h"

class Muon : public MObject {
  public:
    Muon();
    Muon(const TLorentzVector &v);
    virtual ~Muon();

    void setMI(float iso);
    float mi() const;

    bool isTight() const;
    void setTight(bool t);

    const float z0() const;
    float &z0();

    const float z0_exPV() const;
    float &z0_exPV();

    const float d0() const;
    float &d0();

    const float sd0() const;
    float &sd0();

    const int author() const;
    int &author();

    const bool passTrkCuts() const;
    bool &passTrkCuts();

    bool pass() const;
    bool passLoose() const;

    TLorentzVector &momME();
    const TLorentzVector &momME() const;

    TLorentzVector &momMECorr();
    const TLorentzVector &momMECorr() const;

    TLorentzVector &momMS();
    const TLorentzVector &momMS() const;

    TLorentzVector &momID();
    const TLorentzVector &momID() const;

    TLorentzVector &momTrk();
    const TLorentzVector &momTrk() const;

    int charge() const;
    int &charge();

    bool setST(bool b);
    bool st() const;

    bool setSA(bool b);
    bool sa() const;

    bool setCB(bool b);
    bool cb() const;

  protected:
    float m_mi;
    bool m_tight;
    float m_z0;
    float m_z0_exPV;
    float m_d0;
    float m_sd0;
    int m_author;
    bool m_passTrkCuts;

    TLorentzVector m_momME;
    TLorentzVector m_momMECorr;
    TLorentzVector m_momMS;
    TLorentzVector m_momID;
    TLorentzVector m_momTrk;

    int m_charge;
    bool m_st;
    bool m_sa;
    bool m_cb;

};

#endif
