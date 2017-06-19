rm *histograms_*.root
rm Read_*txt

./runsd --data 0 --qcd 0 --syst 1 --files input_GstartT_MG2000MT1400_TtZ.txt --output GstartT_MG2000MT1400_TtZhistograms_e.root,GstartT_MG2000MT1400_TtZhistograms_mu.root >&Read_GstartT_MG2000MT1400_TtZ.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GstartT_MG1000MT600_TtZ.txt --output GstartT_MG1000MT600_TtZhistograms_e.root,GstartT_MG1000MT600_TtZhistograms_mu.root >&Read_GstartT_MG1000MT600_TtZ.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GstartT_MG2000MT1400_TtH.txt --output GstartT_MG2000MT1400_TtHhistograms_e.root,GstartT_MG2000MT1400_TtHhistograms_mu.root >&Read_GstartT_MG2000MT1400_TtH.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GstartT_MG1000MT600_TtH.txt --output GstartT_MG1000MT600_TtHhistograms_e.root,GstartT_MG1000MT600_TtHhistograms_mu.root >&Read_GstartT_MG1000MT600_TtH.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GstartT_MG2000MT1400_TWb.txt --output GstartT_MG2000MT1400_TWbhistograms_e.root,GstartT_MG2000MT1400_TWbhistograms_mu.root >&Read_GstartT_MG2000MT1400_TWb.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GstartT_MG1000MT600_TWb.txt --output GstartT_MG1000MT600_TWbhistograms_e.root,GstartT_MG1000MT600_TWbhistograms_mu.root >&Read_GstartT_MG1000MT600_TWb.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_ttbarBoson.txt --output ttBosonhistograms_e.root,ttBosonhistograms_mu.root >&Read_ttbarBoson.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_SingleTop.txt --output SingleTophistograms_e.root,SingleTophistograms_mu.root >&Read_SingleTop.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_Dibosons.txt --output Dibosonhistograms_e.root,Dibosonhistograms_mu.root > Read_Diboson.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_Zj.txt --output Zjhistograms_e.root,Zjhistograms_mu.root >&Read_Zj.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_Wj.txt --output Wjhistograms_e.root,Wjhistograms_mu.root >&Read_Wj.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_ttbar.txt --output ttbarhistograms_e.root,ttbarhistograms_mu.root >&Read_ttbar.txt&
./runsd --data 1 --syst 0 --qcd 1 --files input_Electrons.txt --output qcdDataElhistograms_e.root,qcdDataElhistograms_mu.root >&Read_qcdEl.txt&
./runsd --data 1 --syst 0 --qcd 0 --files input_Electrons.txt --output DataElhistograms_e.root,DataElhistograms_mu.root >&Read_El.txt&
./runsd --data 1 --syst 0 --qcd 1 --files input_Muons.txt --output qcdDataMuhistograms_e.root,qcdDataMuhistograms_mu.root >&Read_qcdMu.txt&
./runsd --data 1 --syst 0 --qcd 0 --files input_Muons.txt --output DataMuhistograms_e.root,DataMuhistograms_mu.root >&Read_Mu.txt&

./runsd --data 0 --qcd 0 --syst 1 --files input_Gtt_MG500.txt --output GKKtt_115550histograms_e.root,GKKtt_115550histograms_mu.root >&Read_GKKtt_115550.txt& 
./runsd --data 0 --qcd 0 --syst 1 --files input_Gtt_MG600.txt --output GKKtt_115551histograms_e.root,GKKtt_115551histograms_mu.root >&Read_GKKtt_115551.txt& 
./runsd --data 0 --qcd 0 --syst 1 --files input_Gtt_MG700.txt --output GKKtt_115552histograms_e.root,GKKtt_115552histograms_mu.root >&Read_GKKtt_115552.txt& 
./runsd --data 0 --qcd 0 --syst 1 --files input_Gtt_MG800.txt --output GKKtt_115553histograms_e.root,GKKtt_115553histograms_mu.root >&Read_GKKtt_115553.txt& 
./runsd --data 0 --qcd 0 --syst 1 --files input_Gtt_MG900.txt --output GKKtt_119318histograms_e.root,GKKtt_119318histograms_mu.root >&Read_GKKtt_119318.txt& 
./runsd --data 0 --qcd 0 --syst 1 --files input_Gtt_M1000.txt --output GKKtt_115554histograms_e.root,GKKtt_115554histograms_mu.root >&Read_GKKtt_115554.txt& 
./runsd --data 0 --qcd 0 --syst 1 --files input_Gtt_MG1150.txt --output GKKtt_119319histograms_e.root,GKKtt_119319histograms_mu.root >&Read_GKKtt_119319.txt& 
./runsd --data 0 --qcd 0 --syst 1 --files input_Gtt_MG1300.txt --output GKKtt_115555histograms_e.root,GKKtt_115555histograms_mu.root >&Read_GKKtt_115555.txt& 
./runsd --data 0 --qcd 0 --syst 1 --files input_Gtt_MG1600.txt --output GKKtt_115556histograms_e.root,GKKtt_115556histograms_mu.root >&Read_GKKtt_115556.txt& 
./runsd --data 0 --qcd 0 --syst 1 --files input_Gtt_MG1800.txt --output GKKtt_115799histograms_e.root,GKKtt_115799histograms_mu.root >&Read_GKKtt_115799.txt& 
./runsd --data 0 --qcd 0 --syst 1 --files input_Gtt_MG2000.txt --output GKKtt_119582histograms_e.root,GKKtt_119582histograms_mu.root >&Read_GKKtt_119582.txt& 
./runsd --data 0 --qcd 0 --syst 1 --files input_Gtt_MG2550.txt --output GKKtt_158768histograms_e.root,GKKtt_158768histograms_mu.root >&Read_GKKtt_158768.txt& 
./runsd --data 0 --qcd 0 --syst 1 --files input_Gtt_MG2500.txt --output GKKtt_158769histograms_e.root,GKKtt_158769histograms_mu.root >&Read_GKKtt_158769.txt& 
./runsd --data 0 --qcd 0 --syst 1 --files input_Gtt_MG2750.txt --output GKKtt_180575histograms_e.root,GKKtt_180575histograms_mu.root >&Read_GKKtt_180575.txt& 
./runsd --data 0 --qcd 0 --syst 1 --files input_Gtt_MG3000.txt --output GKKtt_180576histograms_e.root,GKKtt_180576histograms_mu.root >&Read_GKKtt_180576.txt& 

./runsd --data 0 --qcd 0 --syst 1 --files input1.txt --output 2ttbarhistograms_e.root,2ttbarhistograms_mu.root
./runsd --data 1 --syst 0 --qcd 1 --files input.txt --output 1qcdDataMuhistograms_e.root,1qcdDataMuhistograms_mu.root >&Read_qcdMu.txt&
./runsd --data 1 --syst 0 --qcd 0 --files input.txt --output 1DataMuhistograms_e.root,1DataMuhistograms_mu.root >&Read_Mu.txt



./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_115550.txt --output GKKtt_115550histograms_e.root,GKKtt_115550histograms_mu.root >&Read_GKKtt_115550.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_115551.txt --output GKKtt_115551histograms_e.root,GKKtt_115551histograms_mu.root >&Read_GKKtt_115551.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_115552.txt --output GKKtt_115552histograms_e.root,GKKtt_115552histograms_mu.root >&Read_GKKtt_115552.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_115553.txt --output GKKtt_115553histograms_e.root,GKKtt_115553histograms_mu.root >&Read_GKKtt_115553.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_115554.txt --output GKKtt_115554histograms_e.root,GKKtt_115554histograms_mu.root >&Read_GKKtt_115554.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_115555.txt --output GKKtt_115555histograms_e.root,GKKtt_115555histograms_mu.root >&Read_GKKtt_115555.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_115556.txt --output GKKtt_115556histograms_e.root,GKKtt_115556histograms_mu.root >&Read_GKKtt_115556.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_115799.txt --output GKKtt_115799histograms_e.root,GKKtt_115799histograms_mu.root >&Read_GKKtt_115799.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_119318.txt --output GKKtt_119318histograms_e.root,GKKtt_119318histograms_mu.root >&Read_GKKtt_119318.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_119319.txt --output GKKtt_119319histograms_e.root,GKKtt_119319histograms_mu.root >&Read_GKKtt_119319.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_119582.txt --output GKKtt_119582histograms_e.root,GKKtt_119582histograms_mu.root >&Read_GKKtt_119582.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_145583.txt --output GKKtt_145583histograms_e.root,GKKtt_145583histograms_mu.root >&Read_GKKtt_145583.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_145585.txt --output GKKtt_145585histograms_e.root,GKKtt_145585histograms_mu.root >&Read_GKKtt_145585.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_145586.txt --output GKKtt_145586histograms_e.root,GKKtt_145586histograms_mu.root >&Read_GKKtt_145586.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_158765.txt --output GKKtt_158765histograms_e.root,GKKtt_158765histograms_mu.root >&Read_GKKtt_158765.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_158766.txt --output GKKtt_158766histograms_e.root,GKKtt_158766histograms_mu.root >&Read_GKKtt_158766.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_158767.txt --output GKKtt_158767histograms_e.root,GKKtt_158767histograms_mu.root >&Read_GKKtt_158767.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_158768.txt --output GKKtt_158768histograms_e.root,GKKtt_158768histograms_mu.root >&Read_GKKtt_158768.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_158769.txt --output GKKtt_158769histograms_e.root,GKKtt_158769histograms_mu.root >&Read_GKKtt_158769.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_158770.txt --output GKKtt_158770histograms_e.root,GKKtt_158770histograms_mu.root >&Read_GKKtt_158770.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_158771.txt --output GKKtt_158771histograms_e.root,GKKtt_158771histograms_mu.root >&Read_GKKtt_158771.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_158772.txt --output GKKtt_158772histograms_e.root,GKKtt_158772histograms_mu.root >&Read_GKKtt_158772.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_158773.txt --output GKKtt_158773histograms_e.root,GKKtt_158773histograms_mu.root >&Read_GKKtt_158773.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_158774.txt --output GKKtt_158774histograms_e.root,GKKtt_158774histograms_mu.root >&Read_GKKtt_158774.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_158775.txt --output GKKtt_158775histograms_e.root,GKKtt_158775histograms_mu.root >&Read_GKKtt_158775.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_158776.txt --output GKKtt_158776histograms_e.root,GKKtt_158776histograms_mu.root >&Read_GKKtt_158776.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_145584.txt --output GKKtt_145584histograms_e.root,GKKtt_145584histograms_mu.root &>Read_GKKtt_145584.txt &

./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_180384.txt --output GKKtt_180384histograms_e.root,GKKtt_180384histograms_mu.root >&Read_GKKtt_180384.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_180385.txt --output GKKtt_180385histograms_e.root,GKKtt_180385histograms_mu.root >&Read_GKKtt_180385.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_180386.txt --output GKKtt_180386histograms_e.root,GKKtt_180386histograms_mu.root >&Read_GKKtt_180386.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_180387.txt --output GKKtt_180387histograms_e.root,GKKtt_180387histograms_mu.root >&Read_GKKtt_180387.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_180389.txt --output GKKtt_180389histograms_e.root,GKKtt_180389histograms_mu.root >&Read_GKKtt_180389.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_180390.txt --output GKKtt_180390histograms_e.root,GKKtt_180390histograms_mu.root >&Read_GKKtt_180390.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_180575.txt --output GKKtt_180575histograms_e.root,GKKtt_180575histograms_mu.root >&Read_GKKtt_180575.txt&
./runsd --data 0 --qcd 0 --syst 1 --files input_GKKtt_180576.txt --output GKKtt_180576histograms_e.root,GKKtt_180576histograms_mu.root >&Read_GKKtt_180575.txt&
