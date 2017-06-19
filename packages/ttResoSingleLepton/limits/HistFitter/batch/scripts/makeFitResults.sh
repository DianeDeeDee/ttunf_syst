



getLimits='grep -A 5 " The computed upper limit is" */lim*'


mkdir -p Results/fits

mainDir=`pwd`

for chan in Resolved Boosted Combined
  do
  echo ${chan}
  mkdir Results/fits/${chan}

  for systType in Nominal
    do
    for signal in Bgr #Z400 Z500  Z750 Z1000 Z1250 Z1500 Z1750 Z2000 Z2250 Z2500 Z3000  KKg500 KKg600 KKg700 KKg800 KKg900 KKg1000 KKg1150 KKg1300 KKg1600 KKg1800 KKg2000 KKg2250 KKg2500 KKg2750 KKg3000
      do
      mkdir  Results/fits/${chan}/${signal}
      

      resultDir=${mainDir}/Results/fits/${chan}/${signal}
      echo ${resultDir}
      cd $chan/default_${systType}_${signal}
      cp fitresult* ${resultDir}
      
      source ${mainDir}/../scripts/makeFitPage.sh > ${resultDir}/showAll.html
      cp results.tgz ${resultDir}/.
      tar -xzf results.tgz
      cp results/*/*.pdf ${resultDir}
      cp results/*/*.eps ${resultDir}
      cp results/*/*.png ${resultDir}
      cp results/*/*.table ${resultDir}
	    
      cd $mainDir

    #echo ${signal[1,1]}
#    if [[ "$signal[1,1]" == "Z" ]]
#	then
#     	file=Zprime-lim_pretag_syst_exp.dat
#	mass=${signal[2,100]}
#    else 
#	file=${signal[1,3]}"-lim_pretag_"${systType}"_exp.dat"
#        mass=${signal[4,100]}
#    fi
#    echo $file $mass 

#    systConfig=""
#    if [[ "$systType" == "stat" ]]
#	then
#	systConfig="ALLSYST"
#    elif  [[ "$systType" == "syst" ]]
#	then
#	systConfig="Nominal"
#    fi
#    cp ${chan}/*${signal}/upperlimit_cls_poi_*_Asym_CLs_grid_ts3.root.png Results/CLSPlots/${chan}

#  line1=(`grep 'expected limit (median)' ${chan}/*${signal}/limit*.log | tail -1`)
#  line2=(`grep 'expected limit (-1 sig)' ${chan}/*${signal}/limit*.log | tail -1`)
#  line3=(`grep 'expected limit (+1 sig)' ${chan}/*${signal}/limit*.log | tail -1`)
#  line4=(`grep 'expected limit (-2 sig)' ${chan}/*${signal}/limit*.log | tail -1`)
#  line5=(`grep 'expected limit (+2 sig)' ${chan}/*${signal}/limit*.log | tail -1`)
  
  # echo $line1 $line2 $line3 $line4

  #echo $mass GeV:  x-section = $line1[6] 16% = $line2[7]   84% = $line3[7]  2.5% = $line4[7]  97.5% = $line5[7]
 # echo $mass $line1[6] $line2[7]  $line3[7]  $line4[7]   $line5[7] >> Results/limits/${chan}/$file
      
    done
  done
done
