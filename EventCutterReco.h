#ifndef EVENTCUTTERRECO_H
#define EVENTCUTTERRECO_H

#include "Event.h"
#include "EventCutter.h"

class EventCutterReco : public EventCutter {
  public:
    //EventCutterReco(); //if not QCD
    EventCutterReco(bool loose = false); //for QCD
    virtual ~EventCutterReco();
    bool select(const Event &e, Event &sel);
    bool ef_mu36Hypo(float eta, float pt);
    bool ef_mu24Hypo(float eta, float pt);
      //  int muonTriggerMatch(const Muon* muon,const Event &e, Event &sel); 
   private:
    bool m_loose; //for QCD
};

#endif

