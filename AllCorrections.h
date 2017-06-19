#ifndef ALLCORRECTIONS_H
#define ALLCORRECTIONS_H

#include "Correction.h"

namespace Corrections {

  class ElectronCorr : public Correction {
    public:
      ElectronCorr() { }
      void run(Event &e);
  };
 class ElJetORLoose : public Correction {
    public:
      ElJetORLoose() { }
      void run(Event &e);
  };
  // for the el-jet OR with electron subtraction
  class ElJetOR : public Correction {
    public:
      ElJetOR() { }
      void run(Event &e);
  };

  class MuonCorr : public Correction {
    public:
      MuonCorr() { }
      void run(Event &e);
  };

  class JetCalib : public Correction {
    public:
      JetCalib() { }
      void run(Event &e);
  };

  class BoostedJES : public Correction {
    public:
      BoostedJES() { }
      void run(Event &e);
  };

  class JES : public Correction {
    public:
      JES() { }
      void run(Event &e);
  };

  class JER : public Correction {
    public:
      JER() { }
      void run(Event &e);
  };

  class JEE : public Correction {
    public:
      JEE() { }
      void run(Event &e);
  };

  class METRecalc : public Correction {
    public:
      METRecalc() { }
      void run(Event &e);
  };

}

#endif

