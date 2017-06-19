#ifndef EVENT_H
#define EVENT_H

#include <vector>
#include "Electron.h"
#include "Muon.h"
#include "Jet.h"
#include "LargeJet.h"


class Reader;
//class RawReader;
#include "Particle.h"

class Correction;
class SkimReader;

class Event {
  public:
    Event();
    virtual ~Event();

    void clear();
    //void read(RawReader &rr);
    void read(Reader &rr);
    void correct(Correction &c);


    std::vector<Electron> &electron();
    std::vector<Muon> &muon();
    std::vector<Jet> &jet();
    std::vector<LargeJet> &largeJet();
    std::vector<LargeJet> &largeJetBB();


    const std::vector<Electron> &electron() const;
    const std::vector<Muon> &muon() const;
    const std::vector<Jet> &jet() const;
    const std::vector<LargeJet> &largeJet() const;
    const std::vector<LargeJet> &largeJetBB() const;


    void clearWeights();
    std::vector<std::string> weightNames() const;
    const std::vector<double> &weight(const std::string &name) const;
    std::vector<double> &weight(const std::string &name, bool create = false);

    void met(float met_x, float met_y);
    TLorentzVector met() const;

    const TLorentzVector &metCorr(const std::string &s) const;
    TLorentzVector &metCorr(const std::string &s, bool create = false);

    const bool triggerElectron() const;
    const bool triggerMuon() const;
    const bool triggerLargeJet() const; 

    const TLorentzVector &triggerJet4mom() const;
    TLorentzVector &triggerJet4mom();

    bool &triggerElectron();
    bool &triggerMuon();
    bool &triggerLargeJet(); 

    bool &passReco();
    bool &passPart();

    bool passReco() const;
    bool passPart() const;

    std::vector<float> &cutFlow();

    std::vector<Particle> &partMom();
    const std::vector<Particle> &partMom() const;

    std::vector<Particle> &partElectron();
    std::vector<Particle> &partMuon();
    std::vector<Jet> &partJet();
    std::vector<LargeJet> &partLargeJet();
    std::vector<LargeJet> &partLargeJetBB();

    const std::vector<Particle> &partElectron() const;
    const std::vector<Particle> &partMuon() const;
    const std::vector<Jet> &partJet() const;
    const std::vector<LargeJet> &partLargeJet() const;
    const std::vector<LargeJet> &partLargeJetBB() const;

    void partMet(float met_x, float met_y);
    TLorentzVector partMet() const;

    int &runNumber();
    const int runNumber() const;

    int &eventNumber();
    const int eventNumber() const;

    bool &isData();
    const bool isData() const;

    bool &isBB();
    const bool isBB() const;

    int &channelNumber();
    const int channelNumber() const;

    float &mu();
    const float mu() const;

    int &npv_good();
    const int npv_good() const;
    int &npv();
    const int npv() const;
    int &npv_met();
    const int npv_met() const;
    float &rho();
    const float rho() const;
    int &lerr();
    const int lerr() const;
    int &terr();
    const int terr() const;
    int &cfl();
    const int cfl() const;

    bool &atlFastII();
    const bool atlFastII() const;

    float &vxZ();
    const float vxZ() const;

    float &mcWeight();
    const float mcWeight() const;
    unsigned int &lbn();
    const unsigned int lbn() const;

    enum DataPeriod {
      perUnknown = 0,
      perAtoB3,
      perB4toB14,
      perC1toC5,
      perC6toD3,
      perD4toX,
      perXon
    };
    const DataPeriod period() const;
    DataPeriod &period();

    const SkimReader *sr() const;

    void setupSystVariation(const std::string &s);

    const int hfor() const;
    int &hfor();

    const bool isTight() const;
    bool &isTight();

  protected:

    bool m_isTight;

    int m_hfor;

    std::vector<float> m_cutFlow;

    std::vector<Electron> m_electron;
    std::vector<Muon> m_muon;
    std::vector<Jet> m_jet;
    std::vector<LargeJet> m_largeJet;
    std::vector<LargeJet> m_largeJetBB;

    TLorentzVector m_triggerJet4mom;

    std::vector<Particle> m_partMom;

    TLorentzVector m_met;
    std::map<std::string, TLorentzVector> m_metCorr;

    bool m_triggerElectron;
    bool m_triggerMuon;
    bool m_triggerLargeJet; //*

    // particle level
    std::vector<Particle> m_partElectron;
    std::vector<Particle> m_partMuon;
    std::vector<Jet> m_partJet;
    std::vector<LargeJet> m_partLargeJet;
    std::vector<LargeJet> m_partLargeJetBB;

    TLorentzVector m_partMet;

    int m_runNumber;
    int m_eventNumber;
    int m_channelNumber;
    bool m_isData;
    bool m_isBB;
    float m_mu;

    int m_npv_good;
    int m_npv;
    int m_npv_met;
    float m_rho;
    int m_lerr;
    int m_terr;
    int m_cfl;
    bool m_atlFastII;

    float m_vxZ;
    DataPeriod m_period;

    std::vector<std::pair<std::string, std::vector<double> > > m_weight;

    float m_mcWeight;
    unsigned int m_lbn;

    bool m_passReco;
    bool m_passPart;
    SkimReader *m_sr;
};

#endif

