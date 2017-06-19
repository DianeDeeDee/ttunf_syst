#include "Event.h"
#include "RawReader.h"
#include "TChain.h"
#include "SkimReader.h"
#include "EventCutterReco.h"
#include "EventCutterPart.h"
/*#include "EventCutterRecoDJ.h"
#include "EventCutterPartDJ.h"
#include "EventCutterRecoBB.h"
#include "EventCutterPartBB.h"
#include "EventCutterRecoTT.h"
#include "EventCutterPartTT.h"
#include "EventCutterRecoTTSig.h"
#include "EventCutterRecoTTBkg.h"
#include "EventCutterRecoTTBkgData.h"
*/#include "MiniTree.h"

#include <iostream>

#include "AllCorrections.h"
#include "Tools.h"

#include "WeightCalc.h"
#include "AllWeightCalcs.h"

#include "Cintex/Cintex.h"

#include "ParseUtils.h"

#include "TInterpreter.h"

#include "PDFReweighter.h"

/*void pdfInit(const std::string &pdfName, const int nPdf) {
  PDFReweighter::Instance()->SetOutputPDFSet(pdfName.c_str());
  std::vector<int> sets;
  for (int z = 0; z < nPdf; ++z) sets.push_back(z);
  PDFReweighter::Instance()->SetErrorSets(sets);
  PDFReweighter::Instance()->PrintInfo();
}

double getPdfWeight(SkimReader *sr) {
  if ( (!sr->mcevt_pdf_x1)    || (!sr->mcevt_pdf_x2) ||
       (!sr->mcevt_pdf_scale) || (!sr->mcevt_pdf_id1) ||
       (!sr->mcevt_pdf_id2)   || (!sr->mcevt_pdf1) ||
       (!sr->mcevt_pdf2) )
    return 0;

  return PDFReweighter::Instance()->GetReweightingCentral(sr->mcevt_pdf_x1->at(0),    sr->mcevt_pdf_x2->at(0),
                                                          sr->mcevt_pdf_scale->at(0), sr->mcevt_pdf_id1->at(0),
                                                          sr->mcevt_pdf_id2->at(0),   sr->mcevt_pdf1->at(0),
                                                          sr->mcevt_pdf2->at(0));
}

std::vector<double> getPdfWeightError(SkimReader *sr) {
  if ( (!sr->mcevt_pdf_x1)    || (!sr->mcevt_pdf_x2) ||
       (!sr->mcevt_pdf_scale) || (!sr->mcevt_pdf_id1) ||
       (!sr->mcevt_pdf_id2)   || (!sr->mcevt_pdf1) ||
       (!sr->mcevt_pdf2) )
    return std::vector<double>();

  return PDFReweighter::Instance()->GetReweightingErrorSets(sr->mcevt_pdf_x1->at(0),    sr->mcevt_pdf_x2->at(0),
                                                            sr->mcevt_pdf_scale->at(0), sr->mcevt_pdf_id1->at(0),
                                                            sr->mcevt_pdf_id2->at(0),   sr->mcevt_pdf1->at(0),
                                                            sr->mcevt_pdf2->at(0));
}
*/
int main(int argc, char **argv) {

  ROOT::Cintex::Cintex::Enable();

  gROOT->ProcessLine("#include <vector>");
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");


  // input parameters
  int isData = 0;
  int isAtlFastII = 0;
  int doParticle = 1;
  std::string files = "skim.root";
  std::string output = "tree.root";
  int syst = 0;
  int mode = 3;
  int maxevents = -1;
  int help = 0;
  int doCorr = 1;
  int useR03 = 0;
  int useCA12 = 0;

  int doLoose = 0;

  std::string pdfName = "";
  int nPdf = -1;

  static struct extendedOption extOpt[] = {
        {"help",          no_argument,       &help,   1, "Display help", &help, extendedOption::eOTInt},
        {"data",         required_argument,     0, 'd', "Is this data?", &isData, extendedOption::eOTInt},
        {"doParticle",         required_argument,     0, 'P', "Do particle-level selection?", &doParticle, extendedOption::eOTInt},
        {"atlFastII",         required_argument,     0, 'a', "Is this AtlFastII? (0/1)", &isAtlFastII, extendedOption::eOTInt},
        {"files",         required_argument,     0, 'f', "Input list of comma-separated D3PD files to apply the selection on", &files, extendedOption::eOTString},
        {"syst",   required_argument,     0, 's', "Do all systematic variations? (0/1)", &syst, extendedOption::eOTInt},
        {"mode",   required_argument,     0, 'D', "0 = ttbar selection, 1 = dijets selection, 2 = boson tagging, 3 = top tagging signal data, 4 = top tagging signal Z' MC, 5 = top tagging bkg MC, 6 = top tagging bkg in data", &mode, extendedOption::eOTInt},
        {"output",       required_argument,     0, 'o', "Name of the ROOT file to use to save the TTree.", &output, extendedOption::eOTString},
        {"maxevents",       required_argument,     0, 'z', "Maximum number of events to presents (-1 means all events )", &maxevents, extendedOption::eOTInt},
        {"corrections",       required_argument,     0, 'C', "Apply corrections in MC?", &doCorr, extendedOption::eOTInt},
        //{"useR03",       required_argument,     0, 'R', "Set to 1 to use R=0.3 for sub-jets and to 0 to use R=0.2.", &useR03, extendedOption::eOTInt},
        //{"useCA12",       required_argument,     0, 'c', "Set to 1 to use CA R=1.2 instead of CA R=1.0 for the largeJetBB.", &useCA12, extendedOption::eOTInt},
        {"doLoose",       required_argument,     0, 'L', "Set to 1 to run the lepton loose selection too.", &doLoose, extendedOption::eOTInt},
        //{"pdfName",       required_argument,     0, 'p', "PDF name for reweighting.", &pdfName, extendedOption::eOTString},
        //{"nPdf",       required_argument,     0, 'n', "Number of PDF error sets.", &nPdf, extendedOption::eOTInt},

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
  EventCutter *se = 0;
  EventCutter *sep = 0;
  EventCutter *se_loose = 0;

  //bool readBB = false;
  /*if (mode == 1) {
    se = new EventCutterRecoDJ;
    sep = new EventCutterPartDJ;
  } else if (mode == 2) {
    se = new EventCutterRecoBB;
    sep = new EventCutterPartBB;
    readBB = true;
  } else if (mode == 0) {*/
    se = new EventCutterReco;
   sep = new EventCutterPart;
   if (doLoose) se_loose = new EventCutterReco(doLoose);//for QCD
   /* sep = new EventCutterPart;
  } else if (mode == 3) {
    se = new EventCutterRecoTT;
    sep = new EventCutterPartTT;
    if (doLoose) se_loose = new EventCutterRecoTT(doLoose);
  } else if (mode == 4) {
    se = new EventCutterRecoTTSig;
    sep = new EventCutterPartTT;
  } else if (mode == 5) {
    se = new EventCutterRecoTTBkg;
    sep = new EventCutterPartTT;
  } else if (mode == 6) {
    se = new EventCutterRecoTTBkgData;
    sep = new EventCutterPartTT;
  }*/
  int ncuts = 30;
  for (int iC = 0; iC < ncuts; iC++) {
    mt.passCuts().push_back(0);
  }
 // rr.readBB(readBB);
 // rr.doParticleLevelSelection((bool) doParticle);



  // corrections
  std::vector<Correction *> vecCorr;
  Correction *metRecalc = 0;
  Correction *ejetOR = 0;
  if (doCorr) {
    vecCorr.push_back(new Corrections::ElectronCorr);
    vecCorr.push_back(new Corrections::MuonCorr);
    vecCorr.push_back(new Corrections::JetCalib);
    vecCorr.push_back(new Corrections::BoostedJES);
    vecCorr.push_back(new Corrections::JES);
    vecCorr.push_back(new Corrections::JER);
    vecCorr.push_back(new Corrections::JEE);

    if (mode <= 3) ejetOR = new Corrections::ElJetOR;

    metRecalc = new Corrections::METRecalc;

    //vecCorr.push_back(ejetOR);
    //vecCorr.push_back(metRecalc);
  }

  std::map<std::string, WeightCalc *> mapWeightCalcReco;
  if (mode <= 3) {
    mapWeightCalcReco.insert(std::make_pair("esf", new WeightCalculator::ElectronSF));
    mapWeightCalcReco.insert(std::make_pair("musf", new WeightCalculator::MuonSF));
    mapWeightCalcReco.insert(std::make_pair("bsf", new WeightCalculator::BtagSF));
  }
  std::map<std::string, WeightCalc *> mapWeightCalcPart;
  mapWeightCalcPart.insert(std::make_pair("mc", new WeightCalculator::MC));
  mapWeightCalcPart.insert(std::make_pair("pileup", new WeightCalculator::Pileup));
  mapWeightCalcPart.insert(std::make_pair("zsf", new WeightCalculator::ZVertexWeight));

  std::vector<std::string> systNames;
  systNames.push_back("");
  if (!isData && syst && mode <= 3) {
    systNames.push_back("jee");
    systNames.push_back("eSmearUp");
    systNames.push_back("eSmearDown");
    systNames.push_back("eRescaleUp");
    systNames.push_back("eRescaleDown");
    systNames.push_back("muSmearIDUP");
    systNames.push_back("muSmearIDLOW");
    systNames.push_back("muSmearMSUP");
    systNames.push_back("muSmearMSLOW");
    systNames.push_back("muSmearSCALEUP");
    systNames.push_back("boostedJER");
    systNames.push_back("boostedJESUp");
    systNames.push_back("boostedJESDown");
    systNames.push_back("boostedJMSUp");
    systNames.push_back("boostedJMSDown");
    systNames.push_back("jesUp");
    systNames.push_back("jesDown");
    systNames.push_back("jer");
    systNames.push_back("metResoSoftTermsUp");
    systNames.push_back("metResoSoftTermsDown");
    systNames.push_back("metScaleSoftTermsUp");
    systNames.push_back("metScaleSoftTermsDown");
  }
  if (!isData && syst && (mode == 4 || mode == 5 || mode == 6)) {
    systNames.push_back("jee");
    systNames.push_back("boostedJER");
    systNames.push_back("boostedJESUp");
    systNames.push_back("boostedJESDown");
    systNames.push_back("boostedJMSUp");
    systNames.push_back("boostedJMSDown");
    systNames.push_back("jesUp");
    systNames.push_back("jesDown");
    systNames.push_back("jer");
  }

  for (std::map<std::string, WeightCalc *>::iterator it = mapWeightCalcPart.begin(); it != mapWeightCalcPart.end(); ++it) {
    mt.weightNames()->push_back(it->first);
  }
  for (std::map<std::string, WeightCalc *>::iterator it = mapWeightCalcReco.begin(); it != mapWeightCalcReco.end(); ++it) {
    mt.weightNames()->push_back(it->first);
  }
  mt.weightNames()->push_back("zzs");
  /*if (pdfName != "") {
    mt.weightNames()->push_back("pdf");
    mt.weightNames()->push_back("pdf_err");
  }*/

  for (std::vector<std::string>::iterator it = systNames.begin(); it != systNames.end(); ++it) {
    mt.systematicsNames()->push_back(*it);
  }

  if (maxevents > c->GetEntries())
    maxevents = c->GetEntries();

  if (maxevents == -1)
    maxevents = c->GetEntries();

  // scale variations
  mt.sumWeights_var().push_back(0);
  mt.sumWeights_var().push_back(0);
  /*if (pdfName != "") {
    pdfInit(pdfName, nPdf);
    for (size_t z = 0; z < nPdf+1; ++z) mt.sumWeights_var().push_back(0);
  }*/

  for (int k = 0; k < maxevents; ++k) {
    if (k % 10 == 0)
      std::cout << "Entry " << k << "/" << maxevents << std::endl;
    c->GetEntry(k);

    e.clear();
    e.read(rr);

    // Apply corrections using e.correct(), if there are any
    for (int c = 0; c < vecCorr.size(); ++c) {
      e.correct(*vecCorr[c]);
    }

    mt.sumWeights() += (mapWeightCalcPart["mc"]->calc(e)[0]) * (mapWeightCalcPart["pileup"]->calc(e)[0]) * (mapWeightCalcPart["zsf"]->calc(e)[0]);
    double su = 1;
    double sd = 1;
    if (e.channelNumber() == 117050) {
      int idxTop = -1;
      int idxAntiTop = -1;
      for (int k = 0; k < e.partMom().size(); ++k) {
        if (e.partMom()[k].id() == 6) {
          idxTop = k;
        }
        if (e.partMom()[k].id() == -6) {
          idxAntiTop = k;
        }
      }
      if (idxTop != -1 && idxAntiTop != -1) {
        float pt1 = e.partMom()[idxTop].mom().Perp();
        float pt2 = e.partMom()[idxAntiTop].mom().Perp();
        if (pt1 < pt2) {
          float tmp = pt2;
          pt2 = pt1;
          pt1 = tmp;
        }
        su = Correction::globalTools->topScaleWeight(pt1*1e-3, pt2*1e-3, true);
        sd = Correction::globalTools->topScaleWeight(pt1*1e-3, pt2*1e-3, false);
      }
    }
    mt.sumWeights_var()[0] += su*(mapWeightCalcPart["mc"]->calc(e)[0]) * (mapWeightCalcPart["pileup"]->calc(e)[0]) * (mapWeightCalcPart["zsf"]->calc(e)[0]);
    mt.sumWeights_var()[1] += sd*(mapWeightCalcPart["mc"]->calc(e)[0]) * (mapWeightCalcPart["pileup"]->calc(e)[0]) * (mapWeightCalcPart["zsf"]->calc(e)[0]);

    /*std::vector<double> pdfWeights;
    double pdfCentral = 1;
    if (pdfName != "") {
      pdfCentral = getPdfWeight(&sr);
      pdfWeights = getPdfWeightError(&sr);
      mt.sumWeights_var()[2] += pdfCentral*(mapWeightCalcPart["mc"]->calc(e)[0]) * (mapWeightCalcPart["pileup"]->calc(e)[0]) * (mapWeightCalcPart["zsf"]->calc(e)[0]);
      for (size_t z = 0; z < nPdf; ++z) {
        mt.sumWeights_var()[3+z] += pdfWeights[z]*(mapWeightCalcPart["mc"]->calc(e)[0]) * (mapWeightCalcPart["pileup"]->calc(e)[0]) * (mapWeightCalcPart["zsf"]->calc(e)[0]);
      }
    }*/

    for (int s = 0; s < systNames.size(); ++s) { // for each syst.
      // s == 0 represents nominal always!

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
      sel.isTight() = true;

      // substitute nominal 4-vectors with syst ones
      e.setupSystVariation(systNames[s]);

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

      bool selectedReco = se->select(e, sel);
      for (int l = 0; l < sel.cutFlow().size(); l++) {
        mt.passCuts().at(l) += sel.cutFlow().at(l);
      }
      if (selectedReco) {
        // mark tight
        sel.isTight() = true;
      } else if (doLoose && se_loose) {
        selectedReco = se_loose->select(e, sel);
        // mark loose
        sel.isTight() = false;
      }

      if (doParticle && sep) {
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

      sel.weight("zzs", true) = std::vector<double>(2, 1);
      sel.weight("zzs")[0] = su;
      sel.weight("zzs")[1] = sd;
      /*if (pdfName != "") {
        sel.weight("pdf", true) = std::vector<double>(1, pdfCentral);
        sel.weight("pdf_err", true) = pdfWeights;
      }*/

      sel.passReco() = selectedReco;
      sel.passPart() = selectedPart;
      if (selectedReco || selectedPart) {
        mt.write(sel, s);
      }
      e.setupSystVariation(systNames[0]);
    }
  }

  for (int c = 0; c < vecCorr.size(); ++c) {
    delete vecCorr[c];
  }
  vecCorr.clear();
  metRecalc = 0;
  ejetOR = 0;

  delete se;
  if (sep) delete sep;
  if (se_loose) delete se_loose;

  return 0;
}

