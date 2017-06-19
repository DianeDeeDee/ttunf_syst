#!/bin/bash

# source configs/Full.sh
source configs/Stat.sh
# source configs/Fit.sh
#source configs/SigSys.sh
#source configs/Freq.sh
# source configs/QuickTest.sh
# source configs/Retry.sh


#======================
# Unblind here
#======================
blindingSwitch="dataNom"  #unblinded analysis
# blindingSwitch="BgrNom" #blinded analysis


# inputFile="MassSpectra_Paper2014_June19.root"
# inputFile="MassSpectra_Paper2014_2014Jul07.root"
# inputFile="MassSpectra_Paper2014_Jul20_modErrors.root"
# inputFile="MassSpectra_Paper2014_Jul29Alldata.root"
# inputFile="MassSpectra_Paper2014_Aug06Alldata.root"
# inputFile="MassSpectra_Paper2014_Aug11.root"
# inputFile="MassSpectra_Paper2014_Aug15_unsmoothed.root"
# inputFile="MassSpectra_Paper2014_Aug15_Skimmed.root"
# inputFile="MassSpectra_Paper2014_Aug15_unsmoothed_Skimmed.root"
inputFile="MassSpectra_Paper2014_Nov07_smoothed_Skimmed.root"
# inputFile="MassSpectra_Paper2014_Oct06_smoothed_Skimmed_AsimovSigInjZ2000.root"
# inputFile="MassSpectra_Paper2014_Aug26_Unsmoothed_Skimmed.root"

HistFitterConfigFile="paperOct06/CombinedTTbar_paper_Nominal.py"
# HistFitterConfigFile="paperAug15/CombinedTTbar_paper_Nominal.py"
# HistFitterConfigFile="paperAug07/CombinedTTbar_paper_Nominal.py"
# HistFitterConfigFile="paperAug11_JESBreakdown/CombinedTTbar_paper_Nominal.py"

settingsdir="settings"

#Make settingsdir if it doesn't already exist
if [ -d ${settingsdir} ]; then
  echo "Settings directory ${settingsdir} exists. Proceeding..."
else
  mkdir ${settingsdir}
fi

rdir=results_${set}_${jobType}_${blindingSwitch}

echo -e "Submitting jobs: "
for channel in "${channels[@]}"
  do
  for sig in "${sigSamples[@]}"
    do
    for excludedsyst in "${excludedSysts[@]}"
      do
      for btagConf in "${btagConfigs[@]}"
        do

        settingsID="${set}_${jobType}_${blindingSwitch}_${btagConf}_${excludedsyst}_${channel}_${sig}"

	echo configs/${channel}Settings${scanDataOrAsimov}.txt
	line=(`grep "${sig}:" configs/${channel}Settings${scanDataOrAsimov}.txt`)
	if [ -z $line ]; then
	    line=(${sig}: -1 -1 1)
	    echo "using defaults"
	fi

	echo $line
 	echo ${line[0]} : ${line[1]} : ${line[2]} : ${line[3]}
	
	#=================================================================================
	#Write current settings to settings file which will be sent to the grid
	settingsfile="settings_${settingsID}.txt"
	
	if [ -f "${settingsdir}/${settingsfile}" ] ; then
	    #delete file
	    rm "${settingsdir}/${settingsfile}"
	    #recreate the file
	    touch "${settingsdir}/${settingsfile}"
	else
	    #create the file
	    touch "${settingsdir}/${settingsfile}"
	fi
	
	export  SignalScaleFactor=1
	if [[  ${line[3]} != "" ]]; then
	    export  SignalScaleFactor=$line[3]
	fi
	
# 	echo "export SignalScaleFactor=$SignalScaleFactor" >> "${settingsdir}/${settingsfile}"
# 	echo "export SIG=$sig" >> "${settingsdir}/${settingsfile}"
# 	echo "export EXCLUDEDSYST=$excludedsyst" >> "${settingsdir}/${settingsfile}"
# 	echo "export EXCLUDEDBGR=$excludedBackgrounds" >> "${settingsdir}/${settingsfile}"
# 	echo "export BTAGCONF=$btagConf" >> "${settingsdir}/${settingsfile}"
# 	echo "export CHANNEL=$channel" >> "${settingsdir}/${settingsfile}"
# 	echo "export INPUTFILE=$inputFile" >> "${settingsdir}/${settingsfile}"
# 	echo "export CONFIGFILE=$HistFitterConfigFile" >> "${settingsdir}/${settingsfile}"
# 	echo "export JOBTYPE=$jobType" >> "${settingsdir}/${settingsfile}"

	echo "#!/bin/bash" >> "${settingsdir}/${settingsfile}"
	echo export ScanUpperLimit=\"${line[1]}\" >> "${settingsdir}/${settingsfile}"
	echo export NumPoints=\"${line[2]}\" >> "${settingsdir}/${settingsfile}"
	echo export SignalScaleFactor=\"${SignalScaleFactor}\" >> "${settingsdir}/${settingsfile}"
	echo export SIG=\"$sig\" >> "${settingsdir}/${settingsfile}"
	echo export EXCLUDEDSYST=\"$excludedsyst\" >> "${settingsdir}/${settingsfile}"
	echo export EXCLUDEDBGR=\"$excludedBackgrounds\" >> "${settingsdir}/${settingsfile}"
	echo export BTAGCONF=\"$btagConf\" >> "${settingsdir}/${settingsfile}"
	echo export CHANNEL=\"$channel\" >> "${settingsdir}/${settingsfile}"
	echo export INPUTFILE=\"$inputFile\" >> "${settingsdir}/${settingsfile}"
	echo export CONFIGFILE=\"$HistFitterConfigFile\" >> "${settingsdir}/${settingsfile}"
	echo export JOBTYPE=\"$jobType\" >> "${settingsdir}/${settingsfile}"
	echo export BLINDING=\"$blindingSwitch\" >> "${settingsdir}/${settingsfile}"
	
	echo "Done writing config file named ${settingsdir}/${settingsfile}"
	#=================================================================================
	
	echo "Limit: signal $sig config $btagConf excluded $excludedsyst"
	
	eval source submitHistFitter.sh ${settingsID}
#  	./runHistFitter.sh ${settingsID}

      done
    done
  done
done
  
