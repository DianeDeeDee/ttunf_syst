#ifndef EVENTCUTTERRECO14_H
#define EVENTCUTTERRECO14_H

#include "Event.h"
#include "EventCutter.h"

class EventCutterReco14 : public EventCutter {
  public:
    EventCutterReco14();
    virtual ~EventCutterReco14();
    bool select(const Event &e, Event &sel);
};

#endif

