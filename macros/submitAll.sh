#!/bin/bash

TAG1=Electrons #exclue GOEGRID_SCRATCHDISK
TAG2=Muons
CERN_USER=dshoaleh
TTSUF=" --mode 0 --doParticle 1 --useR03 0 "
SUF=" --mode 0 --doParticle 0 --useR03 0 "
#DATASUF=" --data 1 --mode 0 --doParticle 0 --useR03 0 "
DATASUF="--data 1 --mode 0 --doParticle 0 --syst 0 --doLoose 1" #for data
python macros/maketar.py


python macros/submit.py - --tag $TAG1 --input datasets/dataElGoodPeriod.txt --data True --user $CERN_USER  --nFilesPerJob 2 --suffix " $DATASUF "> FollowdataEl.txt
python macros/submit.py - --tag $TAG2 --input datasets/dataMuonGoodPeriod.txt --data True --user $CERN_USER --nFilesPerJob 2 --suffix " $DATASUF "> FollowdataMuon.txt

#python macros/submit.py - --tag $TAG --input datasets/dataElectronAB3.txt --user $CERN_USER --nFilesPerJob -1 --suffix " $DATASUF "

#CA-MCGILL-CLUMEQ-T2_SCRATCHDISK

#--destSE 'CA-MCGILL-CLUMEQ-T2_LOCALGROUPDISK'
