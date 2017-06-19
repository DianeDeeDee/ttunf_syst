#ifndef EVENTCUTTERPARTTT_H
#define EVENTCUTTERPARTTT_H

#include "Event.h"
#include "EventCutter.h"

class EventCutterPartTT : public EventCutter {
  public:
    EventCutterPartTT();
    virtual ~EventCutterPartTT();
    bool select(const Event &e, Event &sel);
};

#endif

