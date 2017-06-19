#ifndef EVENTCUTTERPARTDJ_H
#define EVENTCUTTERPARTDJ_H

#include "Event.h"
#include "EventCutter.h"

class EventCutterPartDJ : public EventCutter {
  public:
    EventCutterPartDJ();
    virtual ~EventCutterPartDJ();
    bool select(const Event &e, Event &sel);
};

#endif

