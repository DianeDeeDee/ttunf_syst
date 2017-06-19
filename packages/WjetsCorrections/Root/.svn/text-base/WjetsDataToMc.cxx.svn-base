#include "WjetsCorrections/WjetsDataToMc.h"
#include <iostream>

// Data block
const float WjetsDataToMc::m_elecPretagNominal[NJET_BINS] =     {0.948,0.907,0.881,0.839,0.906,1.098};
const float WjetsDataToMc::m_elecPretagUncertainty[NJET_BINS] = {0.080,0.058,0.123,0.166,0.163,0.331};
const float WjetsDataToMc::m_elecTaggedNominal[NJET_BINS] =     {0.948,0.907,0.881,0.839,0.906,1.098};
const float WjetsDataToMc::m_elecTaggedUncertainty[NJET_BINS] = {0.209,0.184,0.243,0.403,0.271,0.583};
const float WjetsDataToMc::m_muonPretagNominal[NJET_BINS] =     {0.983,0.942,0.870,0.849,0.814,0.687};
const float WjetsDataToMc::m_muonPretagUncertainty[NJET_BINS] = {0.034,0.076,0.097,0.142,0.121,0.180};
const float WjetsDataToMc::m_muonTaggedNominal[NJET_BINS] =     {0.983,0.942,0.870,0.849,0.814,0.687};
const float WjetsDataToMc::m_muonTaggedUncertainty[NJET_BINS] = {0.197,0.194,0.229,0.406,0.230,0.353};


float WjetsDataToMc::correction(unsigned int jetMultiplicity, bool isElec, bool isTagged, int shift) {
  if(jetMultiplicity >= NJET_BINS) {
    std::cerr << "Error: jetMultiplicity in WjetsCorrections::correction out of range " << NJET_1 << " to " << NJET_5_INC << std::endl;
    return -1;
  }
  if(shift < -1 || shift > 1) {
    std::cerr << "Error: shift in WjetsCorrections::correction out of range -1 to +1" << std::endl;
    return -1;
  }

  float scaleFactor = -1.0;
  if(isElec) {
    if(isTagged) {
      scaleFactor = m_elecTaggedNominal[jetMultiplicity];
      if(shift != 0) {
        scaleFactor += shift*m_elecTaggedUncertainty[jetMultiplicity];
      }
    }
    else {
      scaleFactor = m_elecPretagNominal[jetMultiplicity];
      if(shift != 0) {
        scaleFactor += shift*m_elecPretagUncertainty[jetMultiplicity];
      }
    }
  }
  else {
    if(isTagged) {
      scaleFactor = m_muonTaggedNominal[jetMultiplicity];
      if(shift != 0) {
        scaleFactor += shift*m_muonTaggedUncertainty[jetMultiplicity];
      }
    }
    else{
      scaleFactor = m_muonPretagNominal[jetMultiplicity];
      if(shift != 0) {
        scaleFactor += shift*m_muonPretagUncertainty[jetMultiplicity];
      }
    }
  }
  return scaleFactor;
}
