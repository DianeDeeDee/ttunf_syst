#ifndef EVENTCUTTER_H
#define EVENTCUTTER_H

#include "Event.h"

class EventCutter {
  public:
    EventCutter();
    virtual ~EventCutter();
    virtual bool select(const Event &e, Event &sel) = 0;
};

#endif

