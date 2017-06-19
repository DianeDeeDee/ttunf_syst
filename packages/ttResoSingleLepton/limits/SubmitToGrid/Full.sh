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

jobType="Limit"
#jobType="Fit"

scanIntervalMod="Data"


Zsamples=( "Z400" "Z500" "Z750" "Z1000" "Z1250" "Z1500" "Z1750" "Z2000" "Z2250" "Z2500" "Z3000" )
KKgSamples=( "KKg400" "KKg500" "KKg600" "KKg700" "KKg800" "KKg900" "KKg1000" "KKg1150" "KKg1300" "KKg1600" "KKg1800" "KKg2000" "KKg2250" "KKg2500" "KKg2750" "KKg3000" )
KKgWidthSamples1=( "KKg1000_width10pc" "KKg1000_width15pc" "KKg1000_width20pc" "KKg1000_width25pc" "KKg1000_width30pc" "KKg1000_width35pc" "KKg1000_width40pc" )
KKgWidthSamples2=( "KKg2000_width10pc" "KKg2000_width15pc" "KKg2000_width20pc" "KKg2000_width25pc" "KKg2000_width30pc" "KKg2000_width35pc" "KKg2000_width40pc" )
KKgWidthSamples3=( "KKg3000_width10pc" "KKg3000_width15pc" "KKg3000_width20pc" "KKg3000_width25pc" "KKg3000_width30pc" "KKg3000_width35pc" "KKg3000_width40pc" )
RSGSamples=("RSG400" "RSG500" "RSG600" "RSG700" "RSG800" "RSG900" "RSG1000" "RSG1200" "RSG1400" "RSG1600" "RSG1800" "RSG2000" "RSG2500" )
HHSamples=( "HH400" "HH500" "HH750" "HH1000" "HH1250" "HH1500" "HH1750" "HH2000" "HH2250" "HH2500" "HH2750" "HH3000")

sigSamples=("${Zsamples[@]}" "${KKgSamples[@]}" "${KKgWidthSamples1[@]}" "${KKgWidthSamples2[@]}" "${KKgWidthSamples3[@]}" "${RSGSamples[@]}" "${HHSamples[@]}")
# sigSamples=("${KKgSamples[@]}" "${KKgWidthSamples1[@]}" "${KKgWidthSamples2[@]}" "${KKgWidthSamples3[@]}" "${RSGSamples[@]}" "${HHSamples[@]}")
# sigSamples=("${HHSamples[@]}")

# echo ${sigSamples[@]}


channels=("Combined" "Resolved" "Boosted" )
#channels=("Resolved" "Combined")
# channels=("Boosted" )
#channels=("Resolved" )
#channels=("Combined")


extension=""

set=Nominal

