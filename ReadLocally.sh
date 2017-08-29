rm *histograms_*.root
rm Read_*
#rm *Elhistograms_e.root
#rm 1DataMuhistograms_mu.root
#rm 1qcdDataMuhistograms_mu.root
#rm Read_*

./read --data 1 --syst 0 --qcd 1 --files input_Electrons.txt --output qcdDataElhistograms_e.root,qcdDataElhistograms_mu.root >&Read_qcdEl.txt&
./read --data 1 --syst 0 --qcd 0 --files input_Electrons.txt --output DataElhistograms_e.root,DataElhistograms_mu.root >&Read_El.txt&
./read --data 1 --syst 0 --qcd 1 --files input_Muons.txt --output qcdDataMuhistograms_e.root,qcdDataMuhistograms_mu.root >&Read_qcdMu.txt&
./read --data 1 --syst 0 --qcd 0 --files input_Muons.txt --output DataMuhistograms_e.root,DataMuhistograms_mu.root >&Read_Mu.txt&

./read --data 0 --qcd 0 --syst 1 --files input_ttbar.txt --output ttbarhistograms_e.root,ttbarhistograms_mu.root >&Read_ttbar.txt&
