#ifndef EVENTCUTTERPARTBB_H
#define EVENTCUTTERPARTBB_H

#include "Event.h"
#include "EventCutter.h"

class EventCutterPartBB : public EventCutter {
  public:
    EventCutterPartBB();
    virtual ~EventCutterPartBB();
    bool select(const Event &e, Event &sel);
};

#endif

