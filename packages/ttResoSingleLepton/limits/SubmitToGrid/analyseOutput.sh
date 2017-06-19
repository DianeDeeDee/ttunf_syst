#!/bin/bash

spectrumID=$1

source sortGridOutput.sh ${spectrumID}
source extractLimitsFromLogs.sh ${spectrumID}
cd ..
cd limits/plotting
source draw_all_8TeV.sh ${spectrumID}
mkdir ../../results_${spectrumID}/Results/Limits/LimitPlots
mv limits* *log ../../results_${spectrumID}/Results/Limits/LimitPlots
cd ..
cd ..
cd SubmitToGrid
source copyToWebsite.sh ${spectrumID}
