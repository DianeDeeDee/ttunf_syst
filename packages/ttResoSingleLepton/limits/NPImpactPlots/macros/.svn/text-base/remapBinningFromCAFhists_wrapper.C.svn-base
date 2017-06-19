void remapBinningFromCAFhists_wrapper(std::string rebinningListFile, std::string version, std::string anaPath = "HWWAnalysisCode/analysis/HWWlvlv_2012") {
  gSystem->Load("HWWAnalysisCode/lib/libQFramework.so");
  gROOT->ProcessLine(".L macros/optimizeBinning.C+");
  gROOT->ProcessLine(".L macros/remapBinningFromCAFhists.C+");
  remapBinningFromCAFhists(rebinningListFile, version, anaPath);
}

