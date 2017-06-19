#ifndef EVENTCUTTERRECOTTBKGDATA_H
#define EVENTCUTTERRECOTTBKGDATA_H

#include "Event.h"
#include "EventCutter.h"

class EventCutterRecoTTBkgData : public EventCutter {
  public:
    EventCutterRecoTTBkgData();
    virtual ~EventCutterRecoTTBkgData();
    bool select(const Event &e, Event &sel);
};

#endif

