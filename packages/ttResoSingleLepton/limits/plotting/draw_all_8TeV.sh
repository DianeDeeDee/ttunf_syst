#!/bin/bash

#To run, do
# bash draw_all_boosted.sh

echo "Plotting all results for the note"

spectrumID=$1

# doObserved=true
doObserved=false

Nmass_Zp=11
Nmass_KKg=16
Nmass_KKgWidth=7
Nmass_KKGrav=13
Nmass_HH=12

atlas_label[1]=false 
atlas_label[2]=true

t_label[1]="nolabel";
t_label[2]="labelInternal";
#t_label[1]="labelPreliminary";
logscale=true #false
shape[1]=false;
shape[2]=true;
t_shape[1]="quadratic";
t_shape[2]="rectangular";


dtype[1]="Resolved"
dtype[2]="Boosted"
dtype[3]="Combined"


sel_z[1]="Zprime-lim_pretag_syst"
sel_z[2]="Zprime-lim_pretag_stat"

sel_kkg[1]="KKg-lim_pretag_syst"
sel_kkg[2]="KKg-lim_pretag_stat"

sel_kkgWidth[1]="KKg-Width-1TeV-lim_pretag_syst"
sel_kkgWidth[2]="KKg-Width-1TeV-lim_pretag_stat"
sel_kkgWidth[3]="KKg-Width-2TeV-lim_pretag_syst"
sel_kkgWidth[4]="KKg-Width-2TeV-lim_pretag_stat"
sel_kkgWidth[5]="KKg-Width-3TeV-lim_pretag_syst"
sel_kkgWidth[6]="KKg-Width-3TeV-lim_pretag_stat"

sel_kkgrav[1]="KKGrav-lim_pretag_syst"
sel_kkgrav[2]="KKGrav-lim_pretag_stat"

sel_hh[1]="HH-lim_pretag_syst"
sel_hh[2]="HH-lim_pretag_stat"

sel[1]="syst"
sel[2]="stat"

echo ${dtype[1]} ${dtype[2]} ${dtype[3]}
echo ${sel_z[1]} ${sel_z[2]}
echo ${sel[1]} ${sel[2]}

# currentdir=$(pwd)
# datadir="/data/atlas/atlasdata2/behr/TTbarResonanceSearch/LimitSetting/results_v019/Results/Limits/"
# cd ${datadir}
# mkdir LimitPlots
# cd LimitPlots

for l in `seq 2 2`; #1 2 loop over label # nolabel, label
do
    for rq in `seq 2 2`; #loop over plot shape
    do
	for d in `seq 1 3`; #1 3 loop over dir types (reso, boost, comb)
	do
	    dir_w="/data/atlas/atlasdata2/behr/TTbarResonanceSearch/LimitSetting/results_${spectrumID}/Results/Limits/LimitNumbers/"${dtype[d]}
# 	    dir_w="Results/19June2014/LimitNumbers/"${dtype[d]}
# 	    dir_w="Results/May27/limits/"${dtype[d]}

	    for s in `seq 1 2`; #1 2 loop over Z' selections (syst/stat)
	    do
		echo "Zprime"
		echo directory: $dir_w ${sel_z[s]}
		echo  root -l -q -b plotLimits_r17.C+\(\"${dir_w}\"\,\"${sel_z[s]}\"\,\"multi\"\,$Nmass_Zp\,${doObserved}\,${atlas_label[l]}\,$logscale\,${shape[rq]}\) 
		root -l -q -b plotLimits_r17.C+\(\"${dir_w}\"\,\"${sel_z[s]}\"\,\"multi\"\,$Nmass_Zp\,${doObserved}\,${atlas_label[l]}\,$logscale\,${shape[rq]}\) &> log_zprime_${dtype[d]}_${t_label[l]}_${t_shape[rq]}_${sel[s]}.log
	    done
	    
	    for s in `seq 1 2`; #1 2 loop over HH selections (syst/stat)
	    do
		echo "HH"
		echo directory: $dir_w ${sel_hh[s]}
		echo  root -l -q -b plotLimits_r17.C+\(\"${dir_w}\"\,\"${sel_hh[s]}\"\,\"multi\"\,$Nmass_HH\,${doObserved}\,${atlas_label[l]}\,$logscale\,${shape[rq]}\) 
		root -l -q -b plotLimits_r17.C+\(\"${dir_w}\"\,\"${sel_hh[s]}\"\,\"multi\"\,$Nmass_HH\,${doObserved}\,${atlas_label[l]}\,$logscale\,${shape[rq]}\) #&> log_hh_${dtype[d]}_${t_label[l]}_${t_shape[rq]}_${sel[s]}.log
	    done
	    
	    for s in `seq 1 2`; #loop over KKg selections(syst/stat)
	    do
	      echo "KKg"
  		root -l -q -b plotLimits_r17.C+\(\"${dir_w}\"\,\"${sel_kkg[s]}\"\,\"multi\"\,$Nmass_KKg\,${doObserved}\,${atlas_label[l]}\,$logscale\,${shape[rq]}\) &> log_kkg_${dtype[d]}_${t_label[l]}_${t_shape[rq]}_${sel[s]}.log 
	    done
	    
	     for s in `seq 1 2`; #loop over KKGrav selections(syst/stat)
            do
              echo "KKGrav"
                root -l -q -b plotLimits_r17.C+\(\"${dir_w}\"\,\"${sel_kkgrav[s]}\"\,\"multi\"\,$Nmass_KKGrav\,${doObserved}\,${atlas_label[l]}\,$logscale\,${shape[rq]}\) &> log_kkgrav_${dtype[d]}_${t_label[l]}_${t_shape[rq]}_${sel[s]}.log
            done

	     for s in 1 2 3 4 5 6; #loop over KKGrav selections(syst/stat)
	       do
	       echo "KKg width"
	       root -l -q -b plotLimits_r17.C+\(\"${dir_w}\"\,\"${sel_kkgWidth[s]}\"\,\"multi\"\,$Nmass_KKgWidth\,${doObserved}\,${atlas_label[l]}\,$logscale\,${shape[rq]}\) &> log_${sel_kkgWidth[s]}_${dtype[d]}_${t_label[l]}_${t_shape[rq]}.log
	     done

	done
	zip limits_${dir_w_end}_${t_label[l]}_${t_shape[rq]}.zip  limits_*eps log*log
    done
done

# mv ${currentdir}/limits_* .
# mv ${currentdir}/*log .
# cd ${currentdir}






