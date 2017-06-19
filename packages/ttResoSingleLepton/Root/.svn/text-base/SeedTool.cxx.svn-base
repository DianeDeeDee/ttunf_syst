#include "ttResoSingleLepton/SeedTool.h"

#include <iostream>
#include "TRandom3.h"
#include <cmath>

SeedTool::SeedTool()
  : m_eventNumber(0), m_maximumValue(1000000), m_debug(false) {
  for (int k = 0; k < lastSeedType; ++k)
    m_seeds[k] = 0;
}

SeedTool::~SeedTool() {
}

void SeedTool::setEventNumber(int eventNumber) {
  m_eventNumber = eventNumber;
  m_randomGenerator.SetSeed(m_eventNumber);

  if (m_debug) std::cout << "Seeds for event " << eventNumber << ":" << std::endl;
  for (int k = 0; k < lastSeedType; ++k) {
    m_seeds[k] = (int) std::floor(m_randomGenerator.Uniform(m_maximumValue));
    if (m_debug) std::cout << "Seed " << k << ": " << m_seeds[k] << std::endl;
  }

}

int SeedTool::getSeed(SeedType seedType) {
  if (seedType >= lastSeedType) {
    std::cout << "SeedTool::getSeed: The seed type cannot be bigger or equal to " << lastSeedType << ". Returning 0." << std::endl;
    return 0;
  }
  return m_seeds[seedType];
}

void SeedTool::setDebug(bool debug) {
  m_debug = debug;
}

void SeedTool::setMaximumSeed(int maximum) {
  m_maximumValue = maximum;
}

