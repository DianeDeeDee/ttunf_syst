#include "Event.h"
#include "RawReader14.h"
#include "TChain.h"
#include "SkimReader.h"
#include "EventCutterReco14.h"
#include "EventCutterPart14.h"
#include "MiniTree.h"

#include <iostream>

#include "Tools.h"

#include "Correction.h"
#include "WeightCalc.h"
#include "AllWeightCalcs.h"

#include "Cintex/Cintex.h"

#include "ParseUtils.h"

int main(int argc, char **argv) {

  ROOT::Cintex::Cintex::Enable();

  // input parameters
  std::string files = "skim.root";
  std::string output = "tree.root";
  int maxevents = -1;
  int help = 0;

  static struct extendedOption extOpt[] = {
        {"help",          no_argument,       &help,   1, "Display help", &help, extendedOption::eOTInt},
        {"files",         required_argument,     0, 'f', "Input list of comma-separated D3PD files to apply the selection on", &files, extendedOption::eOTString},
        {"output",       required_argument,     0, 'o', "Name of the ROOT file to use to save the TTree.", &output, extendedOption::eOTString},
        {"maxevents",       required_argument,     0, 'z', "Maximum number of events to presents (-1 means all events )", &maxevents, extendedOption::eOTInt},

        {0, 0, 0, 0, 0, 0, extendedOption::eOTInt}
      };

  if (!parseArguments(argc, argv, extOpt) || help) {
    dumpHelp("preselect14", extOpt, "preselect\nOnly select events and dump them in trees for future post-processing.\n");
    return 0;
  } else {
    std::cout << "Dumping options:" << std::endl;
    dumpOptions(extOpt);
  }

  // parse list of input files
  std::vector<std::string> fileList;
  if ( ((files.find("input") != std::string::npos) && (files.find(".txt") != std::string::npos)) ) {
    std::cout << "Using file given as text list." << std::endl;
    std::vector<std::string> inputList;
    // split by ','
    std::string argStr = files;
    for (size_t i = 0,n; i <= argStr.length(); i=n+1) {
      n = argStr.find_first_of(',',i);
      if (n == std::string::npos)
        n = argStr.length();
      std::string tmp = argStr.substr(i,n-i);
      if (tmp != "")
        inputList.push_back(tmp);
    }

    for (int k = 0; k < inputList.size(); ++k) {

      ifstream f(inputList[k].c_str());
      std::string thePathStr;
      while (f.good()) {
        std::stringstream ss;
        f.get(*ss.rdbuf(), '\n');
        if (f.get() == EOF)
          break;

       thePathStr = ss.str();
         if (thePathStr != "") {
          // bug in Ganga/Grid nodes: \\n is substituted instead of a newline
          size_t idx = std::string::npos;
          do {
            idx = thePathStr.find("\\n");
            std::string aFile = thePathStr.substr(0, idx);
            if (aFile.find(".root") != std::string::npos)
              fileList.push_back(aFile);
            if (idx != std::string::npos) {
              thePathStr = thePathStr.substr(idx+2);
            }
          } while (idx != std::string::npos);
        }
      }
    }
    for (std::vector<std::string>::const_iterator it = fileList.begin(); it != fileList.end(); ++it) {
      std::cout << "Input file \""<<*it<<"\""<< std::endl;
    }
  } else {
    // split by ','
    std::string argStr = files;
    for (size_t i = 0,n; i <= argStr.length(); i=n+1) {
      n = argStr.find_first_of(',',i);
      if (n == std::string::npos)
        n = argStr.length();
      std::string tmp = argStr.substr(i,n-i);
      fileList.push_back(tmp);
    }
  }

  bool isData = false;
  bool isAtlFastII = false;
  Correction::globalTools = new Tools(isData, isAtlFastII);

  //TChain *c = new TChain("physics");
  TChain *c = new TChain("truth");
  for (unsigned int iFile=0; iFile<fileList.size(); ++iFile) {
    std::cout << "preselect: Open " << fileList[iFile].c_str() << std::endl;
    if (c->Add(fileList[iFile].c_str(), -1) == 0) {
      std::cout << "ERROR: Failed to open file \""<<fileList[iFile]<<"\". Abort." << std::endl;
      delete c;
      return -1;
    }
  }

  SkimReader sr(c);
  RawReader14 rr(sr);

  Event e;

  MiniTree mt(true, output.c_str());
  Event sel; // selected objects
  EventCutter *se = 0;
  EventCutter *sep = 0;
  //se = new EventCutterReco14;
  sep = new EventCutterPart14;

  std::map<std::string, WeightCalc *> mapWeightCalcPart;
  mapWeightCalcPart.insert(std::make_pair("mc", new WeightCalculator::MC));

  std::vector<std::string> systNames;
  systNames.push_back("");

  for (std::map<std::string, WeightCalc *>::iterator it = mapWeightCalcPart.begin(); it != mapWeightCalcPart.end(); ++it) {
    mt.weightNames()->push_back(it->first);
  }

  for (std::vector<std::string>::iterator it = systNames.begin(); it != systNames.end(); ++it) {
    mt.systematicsNames()->push_back(*it);
  }

  if (maxevents > c->GetEntries())
    maxevents = c->GetEntries();

  if (maxevents == -1)
    maxevents = c->GetEntries();

  for (int k = 0; k < maxevents; ++k) {
    if (k % 10 == 0)
      std::cout << "Entry " << k << "/" << maxevents << std::endl;
    c->GetEntry(k);

    e.clear();
    e.read(rr);

    mt.sumWeights() += (mapWeightCalcPart["mc"]->calc(e)[0]);

    bool selectedPart = false;
    sel.clear();

    sel.lerr() = e.lerr();
    sel.terr() = e.terr();
    sel.cfl() = e.cfl();
    sel.channelNumber() = e.channelNumber();
    sel.isData() = e.isData();
    sel.runNumber() = e.runNumber();
    sel.eventNumber() = e.eventNumber();
    sel.mu() = e.mu();
    sel.npv() = e.npv();
    sel.period() = e.period();
    sel.vxZ() = e.vxZ();
    sel.mcWeight() = e.mcWeight();
    sel.lbn() = e.lbn();
    sel.atlFastII() = e.atlFastII();
    sel.rho() = e.rho();

    // substitute nominal 4-vectors with syst ones
    //e.setupSystVariation(systNames[0]);

    bool selectedReco = false; //se->select(e, sel);
    selectedPart = sep->select(e, sel);

    // add the SFs and weights into the Event::weight for the
    // selected event
    for (std::map<std::string, WeightCalc *>::iterator it = mapWeightCalcPart.begin(); it != mapWeightCalcPart.end(); ++it) {
      if (selectedPart||selectedReco) {
        sel.weight(it->first, true) = it->second->calc(sel);
      } else {
        sel.weight(it->first, true) = std::vector<double>(1,1);
      }
    }

    sel.passReco() = selectedReco;
    sel.passPart() = selectedPart;
    if (selectedReco || selectedPart) {
      mt.write(sel, 0);
    }
  }

  delete se;
  delete sep;

  return 0;
}

