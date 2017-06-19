#!/bin/bash

2>&1

linuxVersion=`lsb_release -d`
echo $linuxVersion
echo `lsb_release -d`



echo SIG $SIG
echo EXCLUDEDSYST $EXCLUDEDSYST
echo EXCLUDEDBGR $EXCLUDEDBGR
echo BTAGCONF $BTAGCONF
echo CHANNEL $CHANNEL
echo INPUTFILE $INPUTFILE
echo JOBTYPE $JOBTYPE

FitMode="disc"
if [[ SIG == "" ]] 
    then
    FitMode="bkg"
fi

echo -e "\n \n tar -xzf Histfitter.tgz"
tar -xzvf Histfitter.tgz > tarlog.txt
ls



echo -e "\n \n"

sig=$SIG
echo "Signal: $sig"

btagsyst=$BTAGCONF
echo "btag systematics type: $btagsyst"

inclBtag="False"

echo "run inclusive: $inclBtag"

excludesyst=$EXCLUDEDSYST
excludebgr=$EXCLUDEDBGR

echo "excluded systematic: $excludesyst"
echo "excluded backgrounds: $excludebgr"

infile=$INPUTFILE

btagBins="'btag_cat1','btag_cat2','btag_cat3'"

if [ "$inclBtag" = "True" ]; then
    inclBtagString="Incl"
else
    inclBtagString="Excl"
fi
channel=$CHANNEL


catString=${btagBins//[-._,\']/}
echo catString: $catString

histfitdir=`pwd`
histfitterversion=HistFitter-00-00-39

resultdir=${histfitdir} #/result
echo home directory: $histfitdir
#echo making $resultdir
#mkdir $resultdir

cd $histfitdir/$histfitterversion
rm -rf results/*
ls 
if [[ $linuxVersion == *Carbon* ]] ; then
    echo -e "\n mySetup-slc6.sh"
    ls mySetup-slc6.sh
    ls setup.sh
    source mySetup-slc6.sh
    source setup.sh
    #cd src 
    #make clean
    #make
    #cd ..
else
    echo -e "\n mySetup-slc5.sh"
    source mySetup-slc5.sh
fi

if [ $JOBTYPE == "Limit" ] ; then
    echo "running Limit Scan"
    echo "HistFitter.py -w -d -l -m ALL -c \"sigSampleName='$sig';btagSystematicsType='$btagsyst';btag_bins=[$btagBins];doInclusive='$inclBtag';excludeSyst='$excludesyst';excludeBackgrounds='$excludebgr';channel='$channel';inputFile='$infile'\" python/myConfig.py > $resultdir/limitBoosted-$inclBtagString-$catString-$excludesyst-$btagsyst-$sig$logsuffix.log "

    date
    #HistFitter.py -w -d -l -m ALL -c "sigSampleName='$sig';btagSystematicsType='$btagsyst';btag_bins=[$btagBins];doInclusive=$inclBtag;dataSampleName='PseudoData';excludeSyst='$excludesyst'"  python/BoostedTTbar_CONFnote14fb.py  > $resultdir/limitBoosted-$inclBtagString-$catString-$excludesyst-$btagsyst-$sig$logsuffix.log 2>&1
    HistFitter.py -w -d -l -m ALL -c "sigSampleName='$sig';btagSystematicsType='$btagsyst';btag_bins=[$btagBins];doInclusive=$inclBtag;excludeSyst='$excludesyst';excludeBackgrounds='$excludebgr';channel='$channel';inputFile='$infile'" python/myConfig.py  > $resultdir/limitBoosted-$inclBtagString-$catString-$excludesyst-$btagsyst-$sig$logsuffix.log 2>&1
    date
    mv results/upperlimit_cls_poi_${sig}_Asym_CLs_grid_ts3.root.eps results/upperlimit_cls_poi_${inclBtagString}_${catString}_${excludesyst}_${btagsyst}_${sig}${logsuffix}_Asym_CLs_grid_ts3.root.eps
fi

if [ $JOBTYPE == "Fit" ] ; then
    echo "running Fit"
    echo "HistFitter.py -F $FitMode -w -f -d -c \"sigSampleName='$sig';btagSystematicsType='$btagsyst';btag_bins=[$btagBins];doInclusive='$inclBtag';excludeSyst='$excludesyst';excludeBackgrounds='$excludebgr';channel='$channel';inputFile='$infile'\" python/myConfig.py > $resultdir/fitBoosted-$inclBtagString-$catString-$excludesyst-$btagsyst-$sig$logsuffix.log  2>&1"
    date
    scripts/HistFitter.py -F $FitMode -w -f -d  -D allPlots -c "sigSampleName='$sig';btagSystematicsType='$btagsyst';btag_bins=[$btagBins];doInclusive=$inclBtag;excludeSyst='$excludesyst';excludeBackgrounds='$excludebgr';channel='$channel';inputFile='$infile'" python/myConfig.py > $resultdir/fitBoosted-$inclBtagString-$catString-$excludesyst-$btagsyst-$sig$logsuffix.log  2>&1
    date
    
    

    echo grepping
    LINE=`grep "have written workspace to file" $resultdir/fitBoosted-$inclBtagString-$catString-$excludesyst-$btagsyst-$sig$logsuffix.log`
    echo LINE=$LINE
    outfilename=""
    
    for L in `echo $LINE`;
     do
      outfilename=$L;
    done
    echo $outfilename
    
    echo -e "\n \n"
    echo "scripts/PrintFitResult.py -w $outfilename -c Boosted_MassSpectraFit_${inclBtagString}_${catString}_${excludesyst}_${btagsyst}_${sig}"
    scripts/PrintFitResult.py -w $outfilename -c Boosted_MassSpectraFit_${inclBtagString}_${catString}_${excludesyst}_${btagsyst}_${sig}
    
     echo test
     echo "\\input{fitresult_Boosted_MassSpectraFit_${inclBtagString}_${catString}_${excludesyst}_${btagsyst}_${sig}} \\clearpage \n"
     echo "\\input{fitresult_Boosted_MassSpectraFit_${inclBtagString}_${catString}_${excludesyst}_${btagsyst}_${sig}} \\clearpage" >> $resultdir/fittable-all.tex
fi


cp *.pdf ..
cp *.png ..
cp *.eps ..
cp *.tex ..
cp  results/* ..
tar -czvf ../results.tgz results > tar.log
echo -e "\n job terminating normally"

#sleep 50