#ifndef EVENTCUTTERRECODJ_H
#define EVENTCUTTERRECODJ_H

#include "Event.h"
#include "EventCutter.h"

class EventCutterRecoDJ : public EventCutter {
  public:
    EventCutterRecoDJ();
    virtual ~EventCutterRecoDJ();
    bool select(const Event &e, Event &sel);
};

#endif

