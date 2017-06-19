#!/bin/bash


btagConfigs=(  "default" )
excludedSysts=("Nominal")
#excludedSysts=("MC_Gen__EWS__IFSR_BoostedJES0")
#excludedSysts=("BoostedJES0")
#excludedSysts=("smallJES")
#excludedSysts=("ALLSYST")
#excludedBackgrounds="qcd_ttV"
#excludedBackgrounds="qcd_ttV_vv_zjets_stop"
excludedBackgrounds="_"

#jobType="Limit"
jobType="Fit"


sigSamples=( "" )

channels=("Resolved" "Boosted" "Combined")
#channels=("Resolved" "Combined")

#channels=("Resolved" )
#channels=("Combined")

#extension="Corr"

set=Nominal

