//choose asimovData or obsData
runPulls("ttH_combined_AllSyst_lepjets_8TeV_PP_valerio_noCut_125_model.root","SigXsecOverSM","combined","ModelConfig","asimovData","Jan22_AllSyst_lepjets_8TeV_PP_valerio_noCut")
// the parameter config_lepjets_8TeV_PP_valerio_noCut is suposed to be an
// optional xml file but it doesn't matter if you skip it
runBreakdown("ttH_combined_AllSyst_lepjets_8TeV_PP_valerio_noCut_125_model.root","combined","ModelConfig","asimovData","SigXsecOverSM","config_lepjets_8TeV_PP_valerio_noCut","add","total",0.005,0.0,"Jan22_AllSyst_lepjets_8TeV_PP_valerio_noCut")
drawPlot_pulls("125","Jan22_AllSyst_lepjets_8TeV_PP_valerio_noCut",1)
