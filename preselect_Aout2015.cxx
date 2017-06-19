#include "Event.h"
#include "RawReader.h"
#include "TChain.h"
#include "SkimReader.h"
#include "EventCutterReco.h"
#include "EventCutterPart.h"
#include "EventCutterRecoDJ.h"
#include "EventCutterPartDJ.h"
#include "EventCutterRecoBB.h"
#include "EventCutterPartBB.h"
#include "MiniTree.h"

#include <iostream>

#include "AllCorrections.h"
#include "Tools.h"

#include "WeightCalc.h"
#include "AllWeightCalcs.h"

#include "Cintex/Cintex.h"

#include "ParseUtils.h"

int main(int argc, char **argv) {

  ROOT::Cintex::Cintex::Enable();

  // input parameters
  int isData = 0;
  int isAtlFastII = 0;
  std::string files = "skim.root";
  std::string output = "tree.root";
  int syst = 0;
  int mode = 0;
  int maxevents = -1;
  // int maxevents = 1000;
  int help = 0;
  int doParticle = 1;
  int doLoose = 0; //for QCD
  
  static struct extendedOption extOpt[] = {
        {"help",          no_argument,       &help,   1, "Display help", &help, extendedOption::eOTInt},
        {"data",         required_argument,     0, 'd', "Is this data?", &isData, extendedOption::eOTInt},
        {"doParticle",         required_argument,     0, 'P', "Do particle-level selection?", &doParticle, extendedOption::eOTInt},
        {"atlFastII",         required_argument,     0, 'a', "Is this AtlFastII? (0/1)", &isAtlFastII, extendedOption::eOTInt},
        {"files",         required_argument,     0, 'f', "Input list of comma-separated D3PD files to apply the selection on", &files, extendedOption::eOTString},
        {"syst",   required_argument,     0, 's', "Do all systematic variations? (0/1)", &syst, extendedOption::eOTInt},
        {"mode",   required_argument,     0, 'D', "0 = ttbar selection, 1 = dijets selection, 2 = boson tagging", &mode, extendedOption::eOTInt},
        {"output",       required_argument,     0, 'o', "Name of the ROOT file to use to save the TTree.", &output, extendedOption::eOTString},
        {"maxevents",       required_argument,     0, 'z', "Maximum number of events to presents (-1 means all events )", &maxevents, extendedOption::eOTInt},
        {"doLoose",       required_argument,     0, 'L', "Set to 1 to run the lepton loose selection too.", &doLoose, extendedOption::eOTInt}, //for QCD
        
        {0, 0, 0, 0, 0, 0, extendedOption::eOTInt}
      };

  if (!parseArguments(argc, argv, extOpt) || help) {
    dumpHelp("preselect", extOpt, "preselect\nOnly select events and dump them in trees for future post-processing.\n");
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
DQ::SetXMLFile("extra/data12_8TeV.periodAllYear_DetStatus-v61-pro14-02_DQDefects-00-01-00_PHYS_StandardGRL_All_Good.xml");
  //Correction::globalTools = new Tools(isData, isAtlFastII);
  Correction::globalTools = new Tools(isData, isAtlFastII, syst);
  TChain *c = new TChain("physics");
  for (unsigned int iFile=0; iFile<fileList.size(); ++iFile) {
    std::cout << "preselect: Open " << fileList[iFile].c_str() << std::endl;
    if (c->Add(fileList[iFile].c_str(), -1) == 0) {
      std::cout << "ERROR: Failed to open file \""<<fileList[iFile]<<"\". Abort." << std::endl;
      delete c;
      return -1;
    }
  }

  SkimReader sr(c);
  RawReader rr(sr);

  Event e;

  MiniTree mt(true, output.c_str());
  Event sel; // selected objects
  EventCutter *se_loose = 0; //for QCD
  //if (doLoose) se_loose = new EventCutterReco(doLoose);//for QCD
  
  EventCutter *se = 0;
  EventCutter *sep = 0;
//  if (mode == 1) {
//    se = new EventCutterRecoDJ;
//    sep = new EventCutterPartDJ;
//    int ncuts=11;
//    for (int iC=0; iC<ncuts; iC++){
//      mt.passCuts().push_back(0.0);
//    }
//  } else if (mode == 2) {
//    se = new EventCutterRecoBB;
//    sep = new EventCutterPartBB;
//    int ncuts=11;
//    for (int iC=0; iC<ncuts; iC++){
//      mt.passCuts().push_back(0.0);
//    }
//  } else if (mode == 0) {
    se = new EventCutterReco;
    sep = new EventCutterPart;
    /*int ncuts=40;
    for (int iC=0; iC<ncuts; iC++){
      mt.passCuts().push_back(0.0);
    }*/
    if (doLoose) se_loose = new EventCutterReco(doLoose);//for QCD
 // }
  int ncuts=40;
  for (int iC=0; iC<ncuts; iC++){
      mt.passCuts().push_back(0.0);
    }
 //rr.doParticleLevelSelection((bool) doParticle);

  // corrections
  std::vector<Correction *> vecCorr;
  Correction *metRecalc = 0;
  Correction *ejetOR = 0;
  Correction *ejetORLoose = 0;
  vecCorr.push_back(new Corrections::ElectronCorr);
  vecCorr.push_back(new Corrections::MuonCorr);
  vecCorr.push_back(new Corrections::JetCalib);
  vecCorr.push_back(new Corrections::ElJetOR); //New
  vecCorr.push_back(new Corrections::BoostedJES);
  vecCorr.push_back(new Corrections::JES);
  vecCorr.push_back(new Corrections::JER);
  vecCorr.push_back(new Corrections::JEE);
  
  ejetOR = new Corrections::ElJetOR;
  ejetORLoose = new Corrections::ElJetORLoose;
  metRecalc = new Corrections::METRecalc;
    
 // Correction *metRecalc = new Corrections::METRecalc;
 // vecCorr.push_back(metRecalc);

  std::map<std::string, WeightCalc *> mapWeightCalcReco;
  mapWeightCalcReco.insert(std::make_pair("esf", new WeightCalculator::ElectronSF));
  mapWeightCalcReco.insert(std::make_pair("musf", new WeightCalculator::MuonSF));
  mapWeightCalcReco.insert(std::make_pair("bsf", new WeightCalculator::BtagSF));
  std::map<std::string, WeightCalc *> mapWeightCalcPart;
  mapWeightCalcPart.insert(std::make_pair("mc", new WeightCalculator::MC));
  mapWeightCalcPart.insert(std::make_pair("pileup", new WeightCalculator::Pileup));
  mapWeightCalcPart.insert(std::make_pair("zsf", new WeightCalculator::ZVertexWeight));

  std::vector<std::string> systNames;
  systNames.push_back("");
  if (!isData && syst) {
    systNames.push_back("eSmearUp");
    systNames.push_back("eSmearDown");
    systNames.push_back("eRescaleUp");
    systNames.push_back("eRescaleDown");
    systNames.push_back("muSmearIDUP");
    systNames.push_back("muSmearIDLOW");
    systNames.push_back("muSmearMSUP");
    systNames.push_back("muSmearMSLOW");
    systNames.push_back("muSmearSCALEUP");
    systNames.push_back("boostedJESUp");
    systNames.push_back("boostedJESDown");
    systNames.push_back("jesUp");
    systNames.push_back("jesDown");
    systNames.push_back("jer");
    systNames.push_back("jee");
    systNames.push_back("metResoSoftTermsUp");
    systNames.push_back("metResoSoftTermsDown");
    systNames.push_back("metScaleSoftTermsUp");
    systNames.push_back("metScaleSoftTermsDown");
  }

  for (std::map<std::string, WeightCalc *>::iterator it = mapWeightCalcPart.begin(); it != mapWeightCalcPart.end(); ++it) {
    mt.weightNames()->push_back(it->first);
  }
  for (std::map<std::string, WeightCalc *>::iterator it = mapWeightCalcReco.begin(); it != mapWeightCalcReco.end(); ++it) {
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
    if (k % 1000 == 0)
      std::cout << "Entry " << k << "/" << maxevents << std::endl;
    c->GetEntry(k);
     e.clear();
    e.read(rr);
//sel.clear();
	
//    sel.clear();
//    sel.read(rr);



    // Apply corrections using e.correct(), if there are any
    for (int c = 0; c < vecCorr.size(); ++c) {
      e.correct(*vecCorr[c]);
    }
    //e.clear();
    mt.sumWeights() += (mapWeightCalcPart["mc"]->calc(e)[0]) * (mapWeightCalcPart["pileup"]->calc(e)[0]) * (mapWeightCalcPart["zsf"]->calc(e)[0]);

    for (int s = 0; s < systNames.size(); ++s) { // for each syst.
      // s == 0 represents nominal always!
      //std::cout << "preselect: syst loop " << std::endl;
      bool selectedPart = false;
      sel.clear();
     // sel.read(rr); //Test

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
      //sel.rho() = e.rho();
      sel.isTight() = true; //for QCD

      // substitute nominal 4-vectors with syst ones: commented for preselect over Data
      e.setupSystVariation(systNames[s]); //comment??

      // (re-)apply MET
      // if the el-jet OR subtraction scheme is used, this
      // would also need to be re-done
       if (ejetOR) {
        e.correct(*ejetOR);
      }
      if (metRecalc) {
        e.correct(*metRecalc);
        // this is necessary for the MET-related systs:
        std::string mets = "";
        if (systNames[s].find("met") != std::string::npos)
        mets = systNames[s];
        TLorentzVector m = e.metCorr(mets);
        e.met(m.Px(), m.Py());
       }

                                                         
      /*e.correct(*metRecalc);
      // this is necessary for the MET-related systs:
      TLorentzVector m = e.metCorr(systNames[s]);
      e.met(m.Px(), m.Py());
*/
      sel.cutFlow().clear();
      sel.cutFlow().resize(40, 0);
      bool selectedReco = se->select(e, sel);
         //for QCD
      if (selectedReco) {
//      	sel.clear();
//	    sel.read(rr);
        // mark tight
        sel.isTight() = true;
      } else if (doLoose && se_loose) {
      	e.clear();
	    e.read(rr);  
	    
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
      //sel.rho() = e.rho();
      sel.isTight() = true; //for QCD

      // substitute nominal 4-vectors with syst ones: commented for preselect over Data
      e.setupSystVariation(systNames[s]); //comment??

      // (re-)apply MET
      // if the el-jet OR subtraction scheme is used, this
      // would also need to be re-done
       if (ejetORLoose) {
        e.correct(*ejetORLoose);
      }
      if (metRecalc) {
        e.correct(*metRecalc);
        // this is necessary for the MET-related systs:
        std::string mets = "";
        if (systNames[s].find("met") != std::string::npos)
        mets = systNames[s];
        TLorentzVector m = e.metCorr(mets);
        e.met(m.Px(), m.Py());
       }
       
        selectedReco = se_loose->select(e, sel);
        // mark loose
        sel.isTight() = false;
      }
//       e.clear();
//	   e.read(rr);
	  // sel.clear();
	
//	      sel.lerr() = e.lerr();
//	      sel.terr() = e.terr();
//	      sel.cfl() = e.cfl();
//	      sel.channelNumber() = e.channelNumber();
//	      sel.isData() = e.isData();
//	      sel.runNumber() = e.runNumber();
//	      sel.eventNumber() = e.eventNumber();
//	      sel.mu() = e.mu();
//	      sel.npv() = e.npv();
//	      sel.period() = e.period();
//	      sel.vxZ() = e.vxZ();
//	      sel.mcWeight() = e.mcWeight();
//	      sel.lbn() = e.lbn();
//	      sel.atlFastII() = e.atlFastII();
//	      //sel.rho() = e.rho();
//	      sel.isTight() = true;
//	
//	      // substitute nominal 4-vectors with syst ones
//	      e.setupSystVariation(systNames[s]);
//	
//	      // (re-)apply MET
//	      // if the el-jet OR subtraction scheme is used, this419	      // would also need to be re-done
//	      if (ejetOR) {
//	        e.correct(*ejetOR);
//	      }
//	      if (metRecalc) {
//	        e.correct(*metRecalc);
//	        // this is necessary for the MET-related systs:
//	        std::string mets = "";
//	        if (systNames[s].find("met") != std::string::npos)
//	          mets = systNames[s];
//	        TLorentzVector m = e.metCorr(mets);
//	        e.met(m.Px(), m.Py());
//         }
      //End for QCD

      for (int l=0; l<sel.cutFlow().size(); l++){
	    mt.passCuts().at(l)+= sel.cutFlow().at(l);
      }
      // no need to apply particle-level selection
      // more than once per systematic variation
      //if (s == 0) {
      //  selectedPart = sep->select(e, sel);
      //}
      if (s == 0 && doParticle) {
        selectedPart = sep->select(e, sel);
      }
      // add the SFs and weights into the Event::weight for the
      // selected event
      for (std::map<std::string, WeightCalc *>::iterator it = mapWeightCalcPart.begin(); it != mapWeightCalcPart.end(); ++it) {
        if (selectedPart||selectedReco) {
          sel.weight(it->first, true) = it->second->calc(sel);
        } else {
          sel.weight(it->first, true) = std::vector<double>(1,1);
        }
      }
      for (std::map<std::string, WeightCalc *>::iterator it = mapWeightCalcReco.begin(); it != mapWeightCalcReco.end(); ++it) {
        if (selectedReco) {
          sel.weight(it->first, true) = it->second->calc(sel);
        } else {
          sel.weight(it->first, true) = std::vector<double>(1,1);
        }
      }

      sel.passReco() = selectedReco;
      sel.passPart() = selectedPart;
      if (selectedReco || selectedPart) {
        mt.write(sel, s);
      }
      e.setupSystVariation(systNames[0]);
    
    } //systNames.size loop
    
  }//event loop

  for (int c = 0; c < vecCorr.size(); ++c) {
    delete vecCorr[c];
  }
  //e.clear();
  vecCorr.clear();
  metRecalc = 0;
  ejetOR = 0;

  delete se;
  if (sep)  delete sep;
  if (se_loose) delete se_loose;

  return 0;
}

