#!/bin/bash

echo "Current directory: "
pwd

ls -ld *
ls -ld */*
ls -ld */*/*

#Optional
2>&1

#==================================================
#Reading in the settings file
#Settings made available as environment variables
#==================================================
settingsID=$1
source settings/"settings_${settingsID}.txt"

#Check the environment setup
linuxVersion=`lsb_release -d`
echo $linuxVersion
echo `lsb_release -d`

echo ScanUpperLimit $ScanUpperLimit
echo NumPoints $NumPoints
echo SignalScaleFactor $SignalScaleFactor
echo SIG $SIG
echo EXCLUDEDSYST $EXCLUDEDSYST
echo EXCLUDEDBGR $EXCLUDEDBGR
echo BTAGCONF $BTAGCONF
echo CHANNEL $CHANNEL
echo INPUTFILE $INPUTFILE
echo CONFIGFILE $CONFIGFILE
echo JOBTYPE $JOBTYPE
echo BLINDING $BLINDING

sig=$SIG
echo "Signal: $sig"

btagsyst=$BTAGCONF
echo "btag systematics type: $btagsyst"

inclBtag="False"

echo "run inclusive btagging [True/False]: $inclBtag"

excludesyst=$EXCLUDEDSYST
excludebgr=$EXCLUDEDBGR

echo "excluded systematic: $excludesyst"
echo "excluded backgrounds: $excludebgr"

infile=$INPUTFILE
configfile=$CONFIGFILE

btagBins="'btag_cat1','btag_cat2','btag_cat3'"

if [ "$inclBtag" = "True" ]; then
    inclBtagString="Incl"
else
    inclBtagString="Excl"
fi

channel=$CHANNEL

catString=${btagBins//[-._,\']/}
echo catString: $catString


#Make sure the fit mode is consistent with the settings
#e.g. discovery fit only makes sense if signal model available
FitMode="disc"
if [[ $sig == "" ]]; then
    FitMode="bkg"
fi
echo FitMode: $FitMode

#=====================================================================
#Setup and compile HistFitter
tar -xzvf HistFitter.tar.gz > HFunzip.log
cd HistFitter
ls
ls -a results
#clean up results directory 
#(should have already been done - just a safety measure)
if [ "$(ls -A results/)" ] ; then
     echo "Non-empty results/ directory! Cleaning now..."
     rm -rf results/*
fi

#Setup the HistFitter environment (SL6 environment by default)
source setup.sh

#Compile HistFitter code
cd src/
make
cd ..
#=====================================================================
#For log files
resultdir=`pwd`
echo ${JOBTYPE}

if [ "${JOBTYPE}" == "Limit" ] ; then
    echo "Running limit scan" #Add -f -D to make pull plots and correlation matrix
    echo "HistFitter.py -w -d -l -m ALL -c \"dataSampleName='$BLINDING';sigSampleName='$sig';btagSystematicsType='$btagsyst';btag_bins=[$btagBins];doInclusive='$inclBtag';excludeSyst='$excludesyst';excludeBackgrounds='$excludebgr';channel='$channel';inputFile='$infile'\" python/${configfile} > $resultdir/limitBoosted-$inclBtagString-$catString-$excludesyst-$btagsyst-$sig$logsuffix.log "
    
    echo ScanUpperLimit $ScanUpperLimit
    echo NumPoints $NumPoints

    date
#     HistFitter.py -w -d -l -m ALL -c "sigSampleName='$sig';btagSystematicsType='$btagsyst';btag_bins=[$btagBins];doInclusive=$inclBtag;dataSampleName='PseudoData';excludeSyst='$excludesyst'"  python/BoostedTTbar_CONFnote14fb.py  > $resultdir/limitBoosted-$inclBtagString-$catString-$excludesyst-$btagsyst-$sig$logsuffix.log 2>&1
    HistFitter.py -w -d -l -m ALL -c "dataSampleName='$BLINDING';sigSampleName='$sig';btagSystematicsType='$btagsyst';btag_bins=[$btagBins];doInclusive=$inclBtag;excludeSyst='$excludesyst';excludeBackgrounds='$excludebgr';channel='$channel';inputFile='$infile'" python/${configfile}  > $resultdir/limitBoosted-$inclBtagString-$catString-$excludesyst-$btagsyst-$sig$logsuffix.log 2>&1
    date
    if [ -f upperlimit_cls_poi_${sig}_Asym_CLs_grid_ts3.root.eps ] ; then
	mv results/upperlimit_cls_poi_${sig}_Asym_CLs_grid_ts3.root.eps results/upperlimit_cls_poi_${inclBtagString}_${catString}_${excludesyst}_${btagsyst}_${sig}${logsuffix}_Asym_CLs_grid_ts3.root.eps
    fi
fi

if [ "${JOBTYPE}" == "Fit" ] ; then
    echo "Running fit"
    echo "HistFitter.py -F $FitMode -w -f -D allPlots -c \"dataSampleName='$BLINDING';sigSampleName='$sig';btagSystematicsType='$btagsyst';btag_bins=[$btagBins];doInclusive='$inclBtag';excludeSyst='$excludesyst';excludeBackgrounds='$excludebgr';channel='$channel';inputFile='$infile'\" python/${configfile} > $resultdir/fitBoosted-$inclBtagString-$catString-$excludesyst-$btagsyst-$sig$logsuffix.log  2>&1"
    date
    scripts/HistFitter.py -F $FitMode -w -f -d  -D allPlots -c "dataSampleName='$BLINDING';sigSampleName='$sig';btagSystematicsType='$btagsyst';btag_bins=[$btagBins];doInclusive=$inclBtag;excludeSyst='$excludesyst';excludeBackgrounds='$excludebgr';channel='$channel';inputFile='$infile'" python/${configfile}  > $resultdir/fitBoosted-$inclBtagString-$catString-$excludesyst-$btagsyst-$sig$logsuffix.log  2>&1
    date

    echo "Grepping line from log file"
    LINE=`grep "have written workspace to file" $resultdir/fitBoosted-$inclBtagString-$catString-$excludesyst-$btagsyst-$sig$logsuffix.log`
    echo LINE=$LINE
    outfilename=""
    
    for L in `echo $LINE`;
     do
      outfilename=$L;
    done
    echo $outfilename
    
    echo -e "\n \n"
    echo "scripts/PrintFitResult.py -w $outfilename -c Boosted_MassSpectraFit_${inclBtagString}_${catString}_${excludesyst}_${btagsyst}_${sig} > $resultdir/limitPlots-$inclBtagString-$catString-$excludesyst-$btagsyst-$sig$logsuffix.log 2>&1"
    scripts/PrintFitResult.py -w $outfilename -c Boosted_MassSpectraFit_${inclBtagString}_${catString}_${excludesyst}_${btagsyst}_${sig} > $resultdir/limitPlots-$inclBtagString-$catString-$excludesyst-$btagsyst-$sig$logsuffix.log 2>&1
    
    echo test
    echo "\\input{fitresult_Boosted_MassSpectraFit_${inclBtagString}_${catString}_${excludesyst}_${btagsyst}_${sig}} \\clearpage \n"
    echo "\\input{fitresult_Boosted_MassSpectraFit_${inclBtagString}_${catString}_${excludesyst}_${btagsyst}_${sig}} \\clearpage" >> $resultdir/fittable-all.tex
fi

echo "Done running HistFitter. Current directory:"
pwd
echo "Now collecting and zipping results"

#Collect other files
if [ -d myResults ] ; then
  echo "myResults/ exists. Cleaning now..."
  rm -rf myResults/*
else
  mkdir myResults
fi

if [ "$(ls -A *.pdf)" ] ; then
  mv *.pdf myResults 
fi
if [ "$(ls -A *.png)" ] ; then
mv *.png myResults 
fi
if [ "$(ls -A *.eps)" ] ; then
mv *.eps myResults 
fi
if [ "$(ls -A *.tex)" ] ; then
mv *.tex myResults 
fi
if [ "$(ls -A *.log)" ] ; then
mv *.log myResults 
fi
if [ "$(ls -A *.txt)" ] ; then
mv *.txt myResults 
fi

#Rename directories to reflect the run ID (for better bookkeeping after download from grid)
mv results results_${settingsID}
mv myResults myResults_${settingsID}

tar -czvf results.tar.gz results_${settingsID} myResults_${settingsID} > ResultsTar.log

#Copy to head directory so that it can be collected by the grid job
cp results.tar.gz ../
