{
  gROOT->ProcessLine(".L macros/runPulls.C+");
  gROOT->ProcessLine(".L macros/runBreakdown.C+");
  gROOT->ProcessLine(".L macros/drawPlot_pulls.C+");
  
  runPulls("Z500_Boosted_data.root","mu_SIG","combined","ModelConfig","obsData","Z500_Boosted_obsData");
  runBreakdown("Z500_Boosted_data.root","combined","ModelConfig","obsData","mu_SIG","config_lepjets","add","total",0.005,0.0,"Z500_Boosted_obsData");
  drawPlot_pulls("125","Z500_Boosted_obsData",1);
  
}