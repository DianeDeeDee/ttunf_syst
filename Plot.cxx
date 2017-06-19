#include "Plot.h"

Plot::Plot(const std::string &filename, const std::vector<std::string> &systs )
  : m_filename(filename), m_hSvc(filename) {
  for (int i = 0; i < systs.size(); ++i) {
    m_hSvc.addSystematics(systs[i]);
  }
  m_hSvc.addTrigger("");
}

Plot::~Plot() {
}

