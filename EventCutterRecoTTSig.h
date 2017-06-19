#ifndef EVENTCUTTERRECOTTSIG_H
#define EVENTCUTTERRECOTTSIG_H

#include "Event.h"
#include "EventCutter.h"

class EventCutterRecoTTSig : public EventCutter {
  public:
    EventCutterRecoTTSig();
    virtual ~EventCutterRecoTTSig();
    bool select(const Event &e, Event &sel);
};

#endif

