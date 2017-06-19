#ifndef SEEDTOOL_H
#define SEEDTOOL_H

#include "TRandom3.h"

class SeedTool {
  public:

    enum SeedType {
      forWhichPeriod = 0,
      forEnergyRescaler,
      forMuonSmear,
      forJetEnergyResolution,
      forPileUpRandomRun,
      lastSeedType
    };

    SeedTool();
    virtual ~SeedTool();

    void setEventNumber(int eventNumber);
    int  getSeed(SeedType seedType);

    void setDebug(bool debug = true);
    void setMaximumSeed(int maximum = 1000000);

  protected:
    int m_eventNumber;
    int m_maximumValue;
    int m_seeds[lastSeedType];

    bool m_debug;

    TRandom3 m_randomGenerator;
};

#endif

