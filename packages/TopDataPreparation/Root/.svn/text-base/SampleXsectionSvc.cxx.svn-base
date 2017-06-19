#include "TopDataPreparation/SampleXsectionSvc.h"
#include <iostream>
#include <cstdlib>

SampleXsectionSvc *SampleXsectionSvc::s_instance = 0;

SampleXsectionSvc::SampleXsectionSvc(std::string inputFile):
  m_inputFile(inputFile),
  m_sampleXsection(0) {
}

//-------------------------------------------------------------------------------

SampleXsectionSvc::~SampleXsectionSvc(){
  if(m_sampleXsection) delete m_sampleXsection;
}

//-------------------------------------------------------------------------------

int SampleXsectionSvc::initialize() {
  m_sampleXsection = new SampleXsection();

  char *rootcoreDir = getenv("ROOTCOREDIR");
  if(rootcoreDir && m_inputFile == "") {
    //m_inputFile = std::string(rootcoreDir) + "/data/TopDataPreparation/AMIXSection-MC11.data";
    m_inputFile = std::string(rootcoreDir) + "/data/TopDataPreparation/XSection-MC12-8TeV.data";
    if (!m_sampleXsection->readFromFile(m_inputFile.c_str())) {
      std::cerr << "SampleXsection::unable to read input file " << m_inputFile << std::endl;
      return 1;
    }
    m_inputFile = std::string(rootcoreDir) + "/data/TopDataPreparation/XSection-MC12-8TeV-4gt.data";
    if (!m_sampleXsection->readFromFile(m_inputFile.c_str())) {
      std::cerr << "SampleXsection::unable to read input file " << m_inputFile << std::endl;
      return 1;
    }
    m_inputFile = std::string(rootcoreDir) + "/data/TopDataPreparation/XSection-MC12-8TeV-Higgs.data";
    if (!m_sampleXsection->readFromFile(m_inputFile.c_str())) {
      std::cerr << "SampleXsection::unable to read input file " << m_inputFile << std::endl;
      return 1;
    }
    return 0;
  }

  if (!m_sampleXsection->readFromFile(m_inputFile.c_str())) {
    std::cerr << "SampleXsection::unable to read input file " << m_inputFile << std::endl;
    return 1;
  }

  return 0;
}

//----------------------------------------------------------------

SampleXsectionSvc *SampleXsectionSvc::svc(std::string inputFile) {
  if(s_instance == 0) {
    s_instance = new SampleXsectionSvc(inputFile);
    
    if(s_instance->initialize() != 0) {
      std::cerr << "Error: initialize failed.  Could not create SampleXsectionSvc." << std::endl;
      delete s_instance;
      s_instance = 0;
    }
  }
  
  return s_instance;
}
