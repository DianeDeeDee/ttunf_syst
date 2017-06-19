#include "TopGoodRunsList/GrlSvc.h"
#include "GoodRunsLists/TGoodRunsListReader.h"
#include <cstdlib>

GrlSvc *GrlSvc::s_instance = 0;

GrlSvc::GrlSvc(std::string grlFileName):
  m_grlFileName(grlFileName),
  m_grl() {
}

//-------------------------------------------------------------------------------

GrlSvc::~GrlSvc(){}

//-------------------------------------------------------------------------------

int GrlSvc::initialize() {
  Root::TGoodRunsListReader reader = Root::TGoodRunsListReader();
  std::string rootcoredir = getenv("ROOTCOREBIN");
  rootcoredir += "/";
  std::cout << "Info: GrlSvs using GRL " << (rootcoredir + m_grlFileName).c_str() << std::endl;
  reader.SetXMLFile(rootcoredir + m_grlFileName);
  if(!reader.Interpret()) return 1;
  m_grl = reader.GetMergedGoodRunsList();

  return 0;
}

//----------------------------------------------------------------

GrlSvc *GrlSvc::svc(std::string grlFileName) {
  if(s_instance == 0) {
    s_instance = new GrlSvc(grlFileName);
    
    if(s_instance->initialize() != 0) {
      std::cerr << "Error: initialize failed.  Could not create GrlSvc." << std::endl;
      delete s_instance;
      s_instance = 0;
    }   
  }
  
  return s_instance;
}
