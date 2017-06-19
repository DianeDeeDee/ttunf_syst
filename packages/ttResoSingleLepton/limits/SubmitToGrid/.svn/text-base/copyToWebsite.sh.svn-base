#!/bin/bash

spectrumID=$1

projectName="Spectrum_${spectrumID}_allIn"
sourceDir="/data/atlas/atlasdata2/behr/TTbarResonanceSearch/LimitSetting/results_${spectrumID}/Results/"
currentdir=$(pwd)

cd ~/public_html/LimitSetting

mkdir ${projectName}
cd ${projectName}

cp -r ${sourceDir}/* .

rm Limits/LimitPlots/*.log

cd ${currentdir}
