#ifndef EVENTCUTTERRECO_H
#define EVENTCUTTERRECO_H

#include "Event.h"
#include "EventCutter.h"

class EventCutterReco : public EventCutter {
  public:
    EventCutterReco();
    virtual ~EventCutterReco();
    bool select(const Event &e, Event &sel);
};

#endif

