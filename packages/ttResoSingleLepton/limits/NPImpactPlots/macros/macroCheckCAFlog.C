#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdio.h>

#include "TString.h"

void macroCheckCAFlog(std::string logFilename, std::string systvar="Nominal") {

  std::cout << "Reading in log file: " << logFilename.c_str() << std::endl;
  std::ifstream inFile(logFilename.c_str());
  if (inFile.fail())  {
    std::cout << "ERROR::Couldn't open file: " << logFilename.c_str() << std::endl;
    exit(1);
  }


  //Reading the log (starting really after Wjets and QCD)
  std::string curline="";
  std::string prevline="";
  bool startChecks=false;
  std::vector<std::string> vecIDs;
  while(inFile) {
    prevline=curline;
    getline(inFile,curline);
    if (!startChecks && curline.find("/diboson/NonWW") != std::string::npos)
	startChecks=true;
    if (startChecks && curline.find("failed") != std::string::npos) {
      TString strID = prevline;
      strID=strID.ReplaceAll("+ ","");
      strID=strID.ReplaceAll("[  [1m....[0m  ]","");
      strID=strID.ReplaceAll(" ","");
      strID=strID.ReplaceAll("[1A","");
      //std::cout<< "\n found a failed corresponding to " << strID.Data() << std::endl;

      vecIDs.push_back(strID.Data());
      //std::cout << " Will dump " << vecIDs[vecIDs.size()-1] << std::endl;
      //if (vecIDs[vecIDs.size()-1]==0) std::cout << prevline.c_str() << std::endl;
    }

  }	    
  inFile.close();


  //Dumping the output
  if (vecIDs.size()) {
    //Removing duplicates
    std::sort( vecIDs.begin(), vecIDs.end() );
    vecIDs.erase( unique( vecIDs.begin(), vecIDs.end() ), vecIDs.end() );
    for (unsigned int iID=0; iID<vecIDs.size(); iID++) {
      //std::cout << "Missing #" << iID <<" : " << vecIDs[iID].c_str() << "\n" << std::endl;
      TString command = TString::Format("echo 'Missing sample %s for %s' >> missingSamples.txt", vecIDs[iID].c_str(), systvar.c_str());
      system(command.Data());
    }
  }

}
