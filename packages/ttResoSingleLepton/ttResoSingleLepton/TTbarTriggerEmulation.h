// Dear emacs, this is -*- c++ -*-
#ifndef TTBARTRIGGEREMULATION_H
#define TTBARTRIGGEREMULATION_H

#include <string>
#include <vector>
#include "TObject.h"

class TTbarTriggerEmulation : public TObject{
  
 public:
  TTbarTriggerEmulation();
  TTbarTriggerEmulation(const TTbarTriggerEmulation &ref);
  virtual ~TTbarTriggerEmulation();
  
  bool GetTrigger(const int jet_n, 
                  const std::vector<float> &trig_EF_jet_emscale_E, 
		  const std::vector<float> &trig_EF_jet_emscale_eta, 
		  const std::vector<std::string> &trig_EF_jet_calibtags,
                  const float treshold, const std::string calibtag, const bool L1, const bool L2, 
		  std::vector<int> *trig_EF_jet_decision=0); // trig_EF_jet_decision is filled with decision of each trigger object

  inline void SetDebug(const bool value){m_debug = value;};
  
 private:
  bool m_debug;
};

#endif
