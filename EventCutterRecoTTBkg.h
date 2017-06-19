#ifndef EVENTCUTTERRECOTTBKG_H
#define EVENTCUTTERRECOTTBKG_H

#include "Event.h"
#include "EventCutter.h"

class EventCutterRecoTTBkg : public EventCutter {
  public:
    EventCutterRecoTTBkg();
    virtual ~EventCutterRecoTTBkg();
    bool select(const Event &e, Event &sel);
};

#endif

