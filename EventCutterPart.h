#ifndef EVENTCUTTERPART_H
#define EVENTCUTTERPART_H

#include "Event.h"
#include "EventCutter.h"

class EventCutterPart : public EventCutter {
  public:
    EventCutterPart();
    virtual ~EventCutterPart();
    bool select(const Event &e, Event &sel);
};

#endif

