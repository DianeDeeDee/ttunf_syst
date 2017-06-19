#!/bin/bash

spectrumID=$1

#path to directory that contains the files downloaded from the grid
outputdir="/data/atlas/atlasdata2/behr/TTbarResonanceSearch/LimitSetting/results_${spectrumID}"

currentdir=$(pwd)

#===============================================
Unpack the results tarball in each directory
for dir in ${outputdir}/user*
do
    cd ${dir}
    tar -xzvf user*results*.tar.gz*
    ls
done

cd ${outputdir}

#===============================================
#Copy Cls plots into a common Results/Limits/ClsPlots/ directory
mkdir Results
mkdir Results/Limits
mkdir Results/Limits/ClsPlots

mkdir Results/Limits/ClsPlots/Stats
mkdir Results/Limits/ClsPlots/Stats/Boosted
mkdir Results/Limits/ClsPlots/Stats/Combined
mkdir Results/Limits/ClsPlots/Stats/Resolved

cd Results/Limits/ClsPlots/Stats/Boosted
cp ${outputdir}/user*Limit*ALLSYST*Boosted*/results_*/*.png .
cd ${outputdir}

cd Results/Limits/ClsPlots/Stats/Combined
cp ${outputdir}/user*Limit*ALLSYST*Combined*/results_*/*.png .
cd ${outputdir}

cd Results/Limits/ClsPlots/Stats/Resolved
cp ${outputdir}/user*Limit*ALLSYST*Resolved*/results_*/*.png .
cd ${outputdir}

mkdir Results/Limits/ClsPlots/Syst
mkdir Results/Limits/ClsPlots/Syst/Boosted
mkdir Results/Limits/ClsPlots/Syst/Combined
mkdir Results/Limits/ClsPlots/Syst/Resolved

cd Results/Limits/ClsPlots/Syst/Boosted
cp ${outputdir}/user*Limit*Nominal*Boosted*/results_*/*.png .
cd ${outputdir}

cd Results/Limits/ClsPlots/Syst/Combined
cp ${outputdir}/user*Limit*Nominal*Combined*/results_*/*.png .
cd ${outputdir}

cd Results/Limits/ClsPlots/Syst/Resolved
cp ${outputdir}/user*Limit*Nominal*Resolved*/results_*/*.png .
cd ${outputdir}

# ===============================================
# #Copy Fit plots into a common Results/Fit directory
# mkdir Results/Fit
# mkdir Results/Fit/Boosted
# mkdir Results/Fit/Combined
# mkdir Results/Fit/Resolved
# 
# cd Results/Fit/Boosted
# # cp ${outputdir}/user*Fit*Boosted*/results_*/MassSpectra_Paper_*/c_corrMatrix* .
# cp ${outputdir}/user*Fit*Boosted*/results_*/MassSpectra_Paper_*/* .
# cp ${outputdir}/user*Fit*Boosted*/myResults_*/fitresult_Boosted_MassSpectraFit_Excl_btagcat1btagcat2btagcat3_Nominal_default_.* .
# cd ${outputdir}
# 
# cd Results/Fit/Combined
# # cp ${outputdir}/user*Fit*Combined*/results_*/MassSpectra_Paper_*/c_corrMatrix* .
# cp ${outputdir}/user*Fit*Combined*/results_*/MassSpectra_Paper_*/* .
# cp ${outputdir}/user*Fit*Combined*/myResults_*/fitresult_Boosted_MassSpectraFit_Excl_btagcat1btagcat2btagcat3_Nominal_default_.* .
# cd ${outputdir}
# 
# cd Results/Fit/Resolved
# # cp ${outputdir}/user*Fit*Resolved*/results_*/MassSpectra_Paper_*/c_corrMatrix* .
# cp ${outputdir}/user*Fit*Resolved*/results_*/MassSpectra_Paper_*/* .
# cp ${outputdir}/user*Fit*Resolved*/myResults_*/fitresult_Boosted_MassSpectraFit_Excl_btagcat1btagcat2btagcat3_Nominal_default_.* .
# cd ${outputdir}
# ===============================================
cd ${currentdir}
