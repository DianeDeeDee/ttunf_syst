#ifndef EVENTCUTTERRECOTT_H
#define EVENTCUTTERRECOTT_H

#include "Event.h"
#include "EventCutter.h"

class EventCutterRecoTT : public EventCutter {
  public:
    EventCutterRecoTT(bool loose = false);
    virtual ~EventCutterRecoTT();
    bool select(const Event &e, Event &sel);

  private:
    bool m_loose;
};

#endif

