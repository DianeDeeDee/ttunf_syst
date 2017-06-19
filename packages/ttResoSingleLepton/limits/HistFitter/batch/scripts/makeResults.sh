
Zsamples=( "Z400" "Z500" "Z750" "Z1000" "Z1250" "Z1500" "Z1750" "Z2000" "Z2250" "Z2500" "Z3000" )
KKgSamples=( "KKg400" "KKg500" "KKg600" "KKg700" "KKg800" "KKg900" "KKg1000" "KKg1150" "KKg1300" "KKg1600" "KKg1800" "KKg2000" "KKg2250" "KKg2500" "KKg2750" "KKg3000" )
KKgWidthSamples1=( "KKg1000_width10pc" "KKg1000_width15pc" "KKg1000_width20pc" "KKg1000_width25pc" "KKg1000_width30pc" "KKg1000_width35pc" "KKg1000_width40pc" )
KKgWidthSamples2=( "KKg2000_width10pc" "KKg2000_width15pc" "KKg2000_width20pc" "KKg2000_width25pc" "KKg2000_width30pc" "KKg2000_width35pc" "KKg2000_width40pc" )
KKgWidthSamples3=( "KKg3000_width10pc" "KKg3000_width15pc" "KKg3000_width20pc" "KKg3000_width25pc" "KKg3000_width30pc" "KKg3000_width35pc" "KKg3000_width40pc" )
RSGSamples=("RSG400" "RSG500" "RSG600" "RSG700" "RSG800" "RSG900" "RSG1000" "RSG1200" "RSG1400" "RSG1600" "RSG1800" "RSG2000" "RSG2500")

sigSamples=(${Zsamples} ${KKgSamples} ${KKgWidthSamples1} ${KKgWidthSamples2} ${KKgWidthSamples3} ${RSGSamples})



getLimits='grep -A 5 " The computed upper limit is" */lim*'

mkdir Results
mkdir Results/limits
mkdir Results/CLSPlots
mkdir Results/CLSPlots_stat
for chan in Resolved Boosted Combined
  do
  


  mkdir Results/limits/${chan}
  mkdir Results/CLSPlots/${chan}
  mkdir Results/CLSPlots_stat/${chan}

  source ../scripts/makeCLSpage.sh > Results/CLSPlots/${chan}/showAll.html
  source ../scripts/makeCLSpage.sh > Results/CLSPlots_stat/${chan}/showAll.html
#  cp ../scripts/showAllCLS.html Results/CLSPlots/${chan}
  
  for systType in stat syst
    do
    #echo ${chan} ${systType}


    echo $chan $systType
    for signal in "${sigSamples[@]}"
      do
      #echo ${chan} ${systType} ${signal}
    #echo ${signal[1,1]}
	  if [[ "$signal[1,1]" == "Z" ]]
	      then
	      file="Zprime-lim_pretag_"${systType}"_exp.dat"
	      fileObs="Zprime-lim_pretag_"${systType}"_obs.dat"
	      mass=${signal[2,100]}
	      
	  elif [[ "$signal" == KKg1000*width* ]]; then
	      file="KKg-Width-1TeV-lim_pretag_"${systType}"_exp.dat"
	      fileObs="KKg-Width-1TeV-lim_pretag_"${systType}"_obs.dat"
	      mass=`echo ${signal[14,15]} |  sed -r 's/width10pc/_/g' `
	  elif [[ "$signal" == KKg2000*width* ]]; then
	      file="KKg-Width-2TeV-lim_pretag_"${systType}"_exp.dat"
	      fileObs="KKg-Width-2TeV-lim_pretag_"${systType}"_obs.dat"
	      mass=`echo ${signal[14,15]} |  sed -r 's/width10pc/_/g' `
	      
	  elif [[ "$signal" == KKg3000*width* ]]; then
	      file="KKg-Width-3TeV-lim_pretag_"${systType}"_exp.dat"
	      fileObs="KKg-Width-3TeV-lim_pretag_"${systType}"_obs.dat"
	      mass=`echo ${signal[14,15]} |  sed -r 's/width10pc/_/g' `
	  elif [[ "$signal" == RSG* ]]; then
	      file="KKGrav-lim_pretag_"${systType}"_exp.dat"
	      fileObs="KKGrav-lim_pretag_"${systType}"_obs.dat"
              mass=${signal[4,100]}
	      

	  else
	      file=${signal[1,3]}"-lim_pretag_"${systType}"_exp.dat"
	      fileObs=${signal[1,3]}"-lim_pretag_"${systType}"_obs.dat"
	      mass=${signal[4,100]}
	  fi
	#  echo $file $mass 
	  
	  systConfig=""
	  if [[ "$systType" == "stat" ]]
	      then
	      systConfig="ALLSYST"
	      cp ${chan}/*${systConfig}*${signal}/upperlimit_cls_poi_*_*_CLs_grid_ts3.root.png Results/CLSPlots_stat/${chan}
	  elif  [[ "$systType" == "syst" ]]
	      then
	      systConfig="Nominal"
	      #systConfig="muonSF_electronSF"
	      cp ${chan}/*${systConfig}*${signal}/upperlimit_cls_poi_*_*_CLs_grid_ts3.root.png Results/CLSPlots/${chan}
	  fi
	  
	  #echo test $systConfig
	  #ls ${chan}/*${systConfig}*${signal}/limit*.log 
	  line0=(`grep 'The computed upper limit is' ${chan}/*${systConfig}*${signal}/limit*.log | tail -1`)
	  line1=(`grep 'expected limit (median)' ${chan}/*${systConfig}*${signal}/limit*.log | tail -1`)
	  line2=(`grep 'expected limit (-1 sig)' ${chan}/*${systConfig}*${signal}/limit*.log | tail -1`)
	  line3=(`grep 'expected limit (+1 sig)' ${chan}/*${systConfig}*${signal}/limit*.log | tail -1`)
	  line4=(`grep 'expected limit (-2 sig)' ${chan}/*${systConfig}*${signal}/limit*.log | tail -1`)
	  line5=(`grep 'expected limit (+2 sig)' ${chan}/*${systConfig}*${signal}/limit*.log | tail -1`)
	  
	  #echo $line1 $line2 $line3 $line4
	  
	  #echo ${signal}
	  echo $mass $line0[8] $line1[6] $line2[7]  $line3[7]  $line4[7]   $line5[7]
	  echo $mass $line1[6] $line2[7]  $line3[7]  $line4[7]   $line5[7] >> Results/limits/${chan}/$file
	  echo $mass $line0[8] >> Results/limits/${chan}/$fileObs
	  
    done
  done
done
