#echo "Step1"
#root -l -b -q prepareSpectraFullBinkkg.C++
#cp PreSpectra_2014Jul29Alldata.root PreSpectra_2014Jul29Alldata.root.bkp

cp inputs/MttSpectra_Oct06.root MttSpectra_Oct06.root
echo "Step1: mergeBkg"
root -l -b -q mergeBackgrounds_full.C++
echo "Step2: Smooth"
#root -l -b -q smoothSpectra.C++
echo "Step3"
#root -l -b -q mergeBackgrounds_full.C\(2\)++
###root -l -b -q mergeCategories.C++

#echo "Step4"
python ReNameTopStatInput.py -m 1 > mylog1 2>&1 &
python ReNameTopStatInput.py -m 2 > mylog1 2>&1 &
python ReNameTopStatInput.py -m 3 > mylog1 2>&1 &
python ReNameTopStatInput.py -m 4 > mylog1 2>&1 &

#echo "Done"
#hadd MassSpectra_Paper2014_Oct06.root MassSpectra_Paper2014_Boosted_e.root  MassSpectra_Paper2014_Resolved_e.root MassSpectra_Paper2014_Boosted_mu.root MassSpectra_Paper2014_Resolved_mu.root
