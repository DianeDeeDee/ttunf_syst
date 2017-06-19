#!/bin/bash

folder=${1}
mass=125

root -b -q macros/runBreakdown.C+\(\"workspaces/${folder}/${mass}.root\",\"combined\",\"ModelConfig\",\"asimovData_1\",\"SigXsecOverSM_HWW\",\"config/HSG3_muhat_breakdown.xml\",\"add\",\"total\",0.005,0.0,\"Breakdown_${folder}\"\)
root -b -q macros/runBreakdown.C+\(\"workspaces/${folder}/${mass}.root\",\"combined\",\"ModelConfig\",\"asimovData_1\",\"SigXsecOverSM_HWW\",\"config/HSG3_muhat_breakdown.xml\",\"add\",\"statistical\",0.005,0.0,\"Breakdown_${folder}\"\)
root -b -q macros/runBreakdown.C+\(\"workspaces/${folder}/${mass}.root\",\"combined\",\"ModelConfig\",\"asimovData_1\",\"SigXsecOverSM_HWW\",\"config/HSG3_muhat_breakdown.xml\",\"add\",\"theoSigInclPlusAccept\",0.005,0.0,\"Breakdown_${folder}\"\)

root -b -q macros/getFiducialCrossSectionMuhatUncertainty.C+\(\"${folder}\"\)
