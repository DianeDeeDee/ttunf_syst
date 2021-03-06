#ifndef RAWREADER_H
#define RAWREADER_H

#include <vector>
#include "Electron.h"
#include "Muon.h"
#include "Jet.h"
#include "LargeJet.h"

#include "SkimReader.h"
#include "Particle.h"

#include "Event.h"

#include "Reader.h"

class RawReader : public Reader {
  public:
    RawReader(SkimReader &sr);
    virtual ~RawReader();

    void electron(std::vector<Electron> &v);
    void muon(std::vector<Muon> &v);
    void jet(std::vector<Jet> &v, std::vector<LargeJet> &u);//, std::vector<LargeJet> &w);
    void met(float &met_x, float &met_y);
    void trigger(bool &triggerElectron, bool &triggerMuon, bool &triggerLargeJet);//, bool &triggerElectron24, bool &triggerElectron60);

    void partElectron(std::vector<Particle> &v);
    void partMuon(std::vector<Particle> &v);
    void partMom(std::vector<Particle> &v);
    void partJet(std::vector<Jet> &v, std::vector<LargeJet> &u);//, std::vector<LargeJet> &w);
    void partMet(float &met_x, float &met_y);

    void readBB(const bool readBB);
    void doParticleLevelSelection(const bool p);

    void general(int &channelNumber, bool &isData, int &runNumber, int &eventNumber, float &mu, int &npv, int &npv_met, float &rho, int &lerr, int &terr, int &cfl,  bool &atlFastII, Event::DataPeriod &p, float &vxZ, float &mcWeight, unsigned int &lbn,  int &npv_good);

  private:
    Event::DataPeriod getPeriodFromDataRun(int run);
    bool m_readBB;
    bool m_doParticleLevelSelection;
};

#endif

