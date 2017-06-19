#!/bin/bash

source configs/Full.sh
#source configs/Stat.sh
#source configs/Fit.sh
#source configs/SigSys.sh
#source configs/Freq.sh


#channels=("Resolved")
#channels=("Boosted")
#channels=( "Combined" "Boosted" )
#channels=( "Combined"  )

#set=FullBJES
#set=Minimal
#set=noSigSys
echo $set

#extension="data"
extension="BgrNom"


#inputFile="MassSpectra_Paper2014_may9_noQCD.root"
#inputFile="MassSpectra_Paper2014.root"
inputFile="MassSpectra_Paper2014_june11.root"
#inputFile="MassSpectra_Paper2014_bjes_may27.root"

#sigSamples=(${Zsamples} ${KKgSamples} ${KKgWidthSamples1} ${KKgWidthSamples2} ${KKgWidthSamples3} ${RSGSamples})
#sigSamples=( ${Zsamples} ${KKgSamples} ) 
#sigSamples=( ${KKgSamples} ) 
#sigSamples=( ${Zsamples}  ) 


rdir=results_${set}_${jobType}_${extension}


ln -sf paperMay27/CombinedTTbar_paper_${set}.py ../HistFitter-00-00-39/python/myConfig.py
tar -czvf Histfitter.tgz ../HistFitter-00-00-39 > tar.log


echo -e "\n \n submitting jobs: "
for channel in "${channels[@]}"
  do
#jobDir=${channel}_JER_bcmtag_IFSR_mcgen_ps_topmass_jes3_12_20
#jobDir=${channel}_JER_bcmtag_IFSR_mcgen_ps_EWS_topmass_jes3_12_20
  jobDir=${channel}
  home=`pwd`

  



  if [ ! -d $rdir/$jobDir  ] ; then
      mkdir -p $rdir/$jobDir
  fi
  
  
  
  for sig in "${sigSamples[@]}"
    do
    for excludedsyst in "${excludedSysts[@]}"
      do
      for btagConf in "${btagConfigs[@]}"
        do

	echo configs/${channel}Settings${scanIntervalMod}.txt
	line=(`grep "${sig}:" configs/${channel}Settings${scanIntervalMod}.txt`)
	if [ -z $line ]; then
	    line=(${sig}: -1 -1 1)
	    echo "using defaults"
	fi

	echo $line
	echo $line[1] : $line[2] : $line[3] : $line[4]
	export ScanUpperLimit=$line[2]
	export NumPoints=$line[3]
	
	export  SignalScaleFactor=1
	if [[  $line[4] != "" ]]; then
	    export  SignalScaleFactor=$line[4]
	    echo $SignalScaleFactor 
	fi
	export SIG=$sig
	export EXCLUDEDSYST=$excludedsyst
	export EXCLUDEDBGR=$excludedBackgrounds
	export BTAGCONF=$btagConf
	export CHANNEL=$channel
	export INPUTFILE=$inputFile
	export JOBTYPE=$jobType
	echo "Limit: signal $sig config $btagConf excluded $excludedsyst"
	
	
	folder=$rdir/${jobDir}/${btagConf}_${excludedsyst}_${sig}
	echo $folder
	if [[ ${sig} == "" ]] ;then
	    folder=${folder}Bgr
	fi
	if [ -d $folder  ] ; then
	    
	    echo output directory already exists, change conditions!
	    #echo $folder
	      #rm -rf $folder
	    continue
	fi
	cp -r template $folder
	cp Histfitter.tgz $folder
	cp -L ../HistFitter-00-00-39/python/myConfig.py $folder
	cd $folder
	  #condor_submit LimitSetting.cmd
	#condor_submit runFit.cmd
	condor_submit runLimit.cmd
	
	cd $home
      done
    done
  done
done
  