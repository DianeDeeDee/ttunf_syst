#!/bin/bash

#Code to extract limits from log files
#By A. Altheimer
#Modified for grid output by K. Behr

spectrumID=$1

currentdir=$(pwd)

#path to directory that contains the files downloaded from the grid
outputdir="/data/atlas/atlasdata2/behr/TTbarResonanceSearch/LimitSetting/results_${spectrumID}"
cd ${outputdir}


Zsamples=( "Z400" "Z500" "Z750" "Z1000" "Z1250" "Z1500" "Z1750" "Z2000" "Z2250" "Z2500" "Z3000" )
KKgSamples=( "KKg400" "KKg500" "KKg600" "KKg700" "KKg800" "KKg900" "KKg1000" "KKg1150" "KKg1300" "KKg1600" "KKg1800" "KKg2000" "KKg2250" "KKg2500" "KKg2750" "KKg3000" )
KKgWidthSamples1=( "KKg1000_width10pc" "KKg1000_width15pc" "KKg1000_width20pc" "KKg1000_width25pc" "KKg1000_width30pc" "KKg1000_width35pc" "KKg1000_width40pc" )
KKgWidthSamples2=( "KKg2000_width10pc" "KKg2000_width15pc" "KKg2000_width20pc" "KKg2000_width25pc" "KKg2000_width30pc" "KKg2000_width35pc" "KKg2000_width40pc" )
KKgWidthSamples3=( "KKg3000_width10pc" "KKg3000_width15pc" "KKg3000_width20pc" "KKg3000_width25pc" "KKg3000_width30pc" "KKg3000_width35pc" "KKg3000_width40pc" )
RSGSamples=("RSG400" "RSG500" "RSG600" "RSG700" "RSG800" "RSG900" "RSG1000" "RSG1200" "RSG1400" "RSG1600" "RSG1800" "RSG2000" "RSG2500")
HHSamples=( "HH400" "HH500" "HH750" "HH1000" "HH1250" "HH1500" "HH1750" "HH2000" "HH2250" "HH2500" "HH2750" "HH3000")

sigSamples=("${Zsamples[@]}" "${KKgSamples[@]}" "${KKgWidthSamples1[@]}" "${KKgWidthSamples2[@]}" "${KKgWidthSamples3[@]}" "${RSGSamples[@]}" "${HHSamples[@]}")
#sigSamples=("${HHSamples[@]}")


getLimits='grep -A 5 " The computed upper limit is" */lim*'


mkdir Results/Limits/LimitNumbers

for chan in Resolved Boosted Combined
  do

  mkdir Results/Limits/LimitNumbers/${chan}
  
  for systType in stat syst
    do
#   systType="syst"

    echo $chan $systType
    
    for signal in "${sigSamples[@]}"
      do
	  echo Signal: ${signal}
	  
	  if [[ $signal == Z* ]]
	      then
	      file="Zprime-lim_pretag_"${systType}"_exp.dat"
	      fileObs="Zprime-lim_pretag_"${systType}"_obs.dat"
	      mass=${signal#"Z"}
	      
	  elif [[ "$signal" == KKg1000*width* ]]; then
	      file="KKg-Width-1TeV-lim_pretag_"${systType}"_exp.dat"
	      fileObs="KKg-Width-1TeV-lim_pretag_"${systType}"_obs.dat"
	      mass=${signal#"KKg1000_width"}
	      mass=${mass%"pc"}
# 	      mass=`echo ${signal[14,15]} |  sed -r 's/width10pc/_/g' `

	  elif [[ "$signal" == KKg2000*width* ]]; then
	      file="KKg-Width-2TeV-lim_pretag_"${systType}"_exp.dat"
	      fileObs="KKg-Width-2TeV-lim_pretag_"${systType}"_obs.dat"
	      mass=${signal#"KKg2000_width"}
	      mass=${mass%"pc"}
# 	      mass=`echo ${signal[14,15]} |  sed -r 's/width10pc/_/g' `

	  elif [[ "$signal" == KKg3000*width* ]]; then
	      file="KKg-Width-3TeV-lim_pretag_"${systType}"_exp.dat"
	      fileObs="KKg-Width-3TeV-lim_pretag_"${systType}"_obs.dat"
	      mass=${signal#"KKg3000_width"}
	      mass=${mass%"pc"}
# 	      mass=`echo ${signal[14,15]} |  sed -r 's/width10pc/_/g' `

	  elif [[ "$signal" == KKg* ]]; then
	      file="KKg-lim_pretag_"${systType}"_exp.dat"
	      fileObs="KKg-lim_pretag_"${systType}"_obs.dat"
	      mass=${signal#"KKg"}
	      echo Current signal is ${signal}!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  elif [[ "$signal" == RSG* ]]; then
	      file="KKGrav-lim_pretag_"${systType}"_exp.dat"
	      fileObs="KKGrav-lim_pretag_"${systType}"_obs.dat"
	      mass=${signal#"RSG"}
	      
	  elif [[ "$signal" == HH* ]]; then
	      file="HH-lim_pretag_"${systType}"_exp.dat"
	      fileObs="HH-lim_pretag_"${systType}"_obs.dat"
	      mass=${signal#"HH"}
	      
	  else
	      file=${signal[1,3]}"-lim_pretag_"${systType}"_exp.dat"
	      fileObs=${signal[1,3]}"-lim_pretag_"${systType}"_obs.dat"
	      mass="Unknown_Signal"
	  fi
	  
	  systConfig=""
	  if [[ "$systType" == "stat" ]]
	      then
	      systConfig="ALLSYST"
	  elif  [[ "$systType" == "syst" ]]
	      then
	      systConfig="Nominal"
# 	      systConfig="JVF"
	  fi
	  
	  
	  
	  line0=(`grep 'The computed upper limit is' ${outputdir}/user*Limit*${systConfig}*${chan}*${signal}.*/myResults_*/limit*.log | tail -1`)
	  line1=(`grep 'expected limit (median)' ${outputdir}/user*Limit*${systConfig}*${chan}*${signal}.*/myResults_*/limit*.log | tail -1`)
	  line2=(`grep 'expected limit (-1 sig)' ${outputdir}/user*Limit*${systConfig}*${chan}*${signal}.*/myResults_*/limit*.log | tail -1`)
	  line3=(`grep 'expected limit (+1 sig)' ${outputdir}/user*Limit*${systConfig}*${chan}*${signal}.*/myResults_*/limit*.log | tail -1`)
	  line4=(`grep 'expected limit (-2 sig)' ${outputdir}/user*Limit*${systConfig}*${chan}*${signal}.*/myResults_*/limit*.log | tail -1`)
	  line5=(`grep 'expected limit (+2 sig)' ${outputdir}/user*Limit*${systConfig}*${chan}*${signal}.*/myResults_*/limit*.log | tail -1`)
	  
	  echo What directories to expect: ${outputdir}/user*Limit*${systConfig}*${chan}*${signal}.*
	  echo
	  echo Line ${line1[@]} ${line2[@]} ${line3[@]} ${line4[@]} ${line5[@]}
	  echo 
	  echo Compare line1: $mass ${line0[7]} ${line1[5]} ${line2[6]}  ${line3[6]}  ${line4[6]} ${line5[6]}
	  echo Compare line2: $mass ${line1[5]} ${line2[6]}  ${line3[6]}  ${line4[6]}  ${line5[6]} 
	  echo $mass ${line1[5]} ${line2[6]}  ${line3[6]}  ${line4[6]}   ${line5[6]} >> Results/Limits/LimitNumbers/${chan}/$file
	  echo $mass ${line0[7]} >> Results/Limits/LimitNumbers/${chan}/$fileObs
	  
    done
  done
done

cd ${currentdir}
