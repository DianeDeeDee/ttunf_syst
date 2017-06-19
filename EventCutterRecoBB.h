#ifndef EVENTCUTTERRECOBB_H
#define EVENTCUTTERRECOBB_H

#include "Event.h"
#include "EventCutter.h"

class EventCutterRecoBB : public EventCutter {
  public:
    EventCutterRecoBB();
    virtual ~EventCutterRecoBB();
    bool select(const Event &e, Event &sel);
};

#endif

