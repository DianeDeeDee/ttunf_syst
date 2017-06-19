#!/bin/bash

#=================================================
# File to set prun options and launch grid jobs
#=================================================

#tar -czvf HistFitter.tar.gz ../HistFitter

export settingsID=$1

if [[ ${settingsID} == "" ]] ;then
  echo "Empty outputFileName given --> using default: HistFitter_Dummy"
  settingsID="dummy_settings"
fi

export tag="v088_07Nov_smoothed"

user="kbehr"

# site="--excludedSite ANALY_RALPP_SLC6,ANALY_NIKHEF-ELPROD_SHORT,ANALY_BNL_SHORT,ANALY_ROMANIA07,ANALY_FZK_SHORT,ANALY_LUNARC,ANALY_GOEGRID,ANALY_GRIF-IRFU,ANALY_GRIF-LPNHE,ANALY_HEPHY-UIBK,ANALY_RHUL_SL6,ANALY_RRC-KI-T1,ANALY_RRC-KI,ANALY-GRIF-LAL,ANALY_INFN-LECCE,ANALY_wuppertalprod,ANALY_NCG-INGRID-PT_SL6,ANALY_INFN-BOLOGNA-T3,ANALY_IL-TAU-HEP-CREAM,ANALY_WEIZMANN-CREAM"
site="--site ANALY_BNL_LONG"

jobs="--nJobs 1"

export cmdLT="prun --exec \"./runHistFitter.sh ${settingsID}\" --outputs results.tar.gz "

export cmdOpt="--rootVer=5.34/14 --cmtConfig=x86_64-slc6-gcc47-opt --match \"*root*\" "${jobs}" "${site}" "

export cmdIn="--inTarBall ../${tag}.tar.gz "
export cmdOut="--extFile HistFitter.tar.gz --outTarBall ../${tag}.tar.gz --noSubmit"
# eval $cmdLT$cmdOpt$cmdOut --outDS user.${user}.HistFitter.Test.${tag}

# export cmdOut="--extFile HistFitter.tar.gz "

if [[ ! -f "../${tag}.tar.gz" ]] ; then
  echo "Input tarball doesn't exist. Creating one..."
  echo $cmdLT$cmdOpt$cmdOut --outDS user.${user}.${settingsID}.${tag} ${cmdOut}
  eval $cmdLT$cmdOpt$cmdOut --outDS user.${user}.${settingsID}.${tag} ${cmdOut}
fi

echo $cmdLT$cmdOpt$cmdIn --outDS user.${user}.${settingsID}.${tag} --memory=10240
eval $cmdLT$cmdOpt$cmdIn --outDS user.${user}.${settingsID}.${tag} --memory=10240
