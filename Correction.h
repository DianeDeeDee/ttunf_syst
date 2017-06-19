#ifndef CORRECTION_H
#define CORRECTION_H

#include <vector>
#include "Electron.h"
#include "Muon.h"
#include "Jet.h"
#include "LargeJet.h"
#include "Tools.h"

class Event;

class Correction {
  public:
    static Tools *globalTools;
    Correction();
    virtual ~Correction();

    virtual void run(Event &e) = 0;
  protected:
    Tools *tools;
};

#endif

