#ifndef ALLWEIGHTCALCS_H
#define ALLWEIGHTCALCS_H

#include "WeightCalc.h"
#include <vector>

class Event;

namespace WeightCalculator {

  class MC : public WeightCalc {
    public:
      MC() { }
      virtual ~MC() { }
      std::vector<double> calc(Event &e);
  };

  class Pileup : public WeightCalc {
    public:
      Pileup() { }
      virtual ~Pileup() { }
      std::vector<double> calc(Event &e);
  };

  class ElectronSF : public WeightCalc {
    public:
      ElectronSF() { }
      virtual ~ElectronSF() { }
      std::vector<double> calc(Event &e);
  };

  class MuonSF : public WeightCalc {
    public:
      MuonSF() { }
      virtual ~MuonSF() { }
      std::vector<double> calc(Event &e);
  };

  class JVFSF : public WeightCalc {
    public:
     JVFSF() { }
     virtual ~JVFSF() { }
     std::vector<double> calc(Event &e);
  };

  class BtagSF : public WeightCalc {
    public:
      BtagSF() { }
      virtual ~BtagSF() { }
      std::vector<double> calc(Event &e);
  };

  class ZVertexWeight : public WeightCalc {
    public:
      ZVertexWeight() { }
      virtual ~ZVertexWeight() { }
      std::vector<double> calc(Event &e);
  };
}

#endif

