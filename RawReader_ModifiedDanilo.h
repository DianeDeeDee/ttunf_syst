#ifndef RAWREADER_H
#define RAWREADER_H

#include <vector>
#include "Electron.h"
#include "Muon.h"
#include "Jet.h"
#include "LargeJet.h"

#include "SkimReader.h"
#include "Particle.h"
#include "Reader.h"
#include "Event.h"

class RawReader {
  public:
    RawReader(SkimReader &sr);
    virtual ~RawReader();

    void electron(std::vector<Electron> &v);
    void muon(std::vector<Muon> &v);
    void jet(std::vector<Jet> &v, std::vector<LargeJet> &u);
    //void jet(std::vector<Jet> &v, std::vector<LargeJet> &u, std::vector<LargeJet> &w);
    void met(float &met_x, float &met_y);
    void trigger(bool &triggerElectron, bool &triggerMuon, bool &triggerLargeJet);
    void readBB(const bool readBB);
    void partElectron(std::vector<Particle> &v);
    void partMuon(std::vector<Particle> &v);
    void partJet(std::vector<Jet> &v, std::vector<LargeJet> &u, std::vector<LargeJet> &w);
    void partMet(float &met_x, float &met_y);
    void partMom(std::vector<Particle> &v);
    void general(int &channelNumber, bool &isData, int &runNumber, int &eventNumber, float &mu, int &npv, int &npv_met, float &rho, int &lerr, int &terr, int &cfl,  bool &atlFastII, Event::DataPeriod &p, float &vxZ, float &mcWeight, unsigned int &lbn, int &npv_good);
   void doParticleLevelSelection(const bool p);
    //void general(int &channelNumber, bool &isData, int &runNumber, int &eventNumber, float &mu, int &npv, int &npv_met, float &rho, int &lerr, int &terr, int &cfl,  bool &atlFastII, Event::DataPeriod &p, float &vxZ, float &mcWeight, unsigned int &lbn, int &hfor, int &npv_good);
    //void sr(SkimReader *&_sr);
     void sr(SkimReader *_sr);
      SkimReader *sr();
  protected:
    SkimReader *m_sr;

  private:
    Event::DataPeriod getPeriodFromDataRun(int run);
 bool m_readBB;
    bool m_doParticleLevelSelection;
};

#endif

