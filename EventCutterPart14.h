#ifndef EVENTCUTTERPART14_H
#define EVENTCUTTERPART14_H

#include "Event.h"
#include "EventCutter.h"

class EventCutterPart14 : public EventCutter {
  public:
    EventCutterPart14();
    virtual ~EventCutterPart14();
    bool select(const Event &e, Event &sel);
};

#endif

