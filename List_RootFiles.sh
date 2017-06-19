ls user.dshoaleh.ttunf8TeV_periodA.Electrons_*/* > input_Electrons.txt
ls user.dshoaleh.ttunf8TeV_periodB.Electrons_*/* >> input_Electrons.txt
ls user.dshoaleh.ttunf8TeV_periodC.Electrons_*/* >> input_Electrons.txt
ls user.dshoaleh.ttunf8TeV_periodD.Electrons_*/*>> input_Electrons.txt
ls user.dshoaleh.ttunf8TeV_periodE.Electrons_*/* >> input_Electrons.txt
ls user.dshoaleh.ttunf8TeV_periodG.Electrons_*/* >> input_Electrons.txt
ls user.dshoaleh.ttunf8TeV_periodH.Electrons_*/* >> input_Electrons.txt
ls user.dshoaleh.ttunf8TeV_periodI.Electrons_*/* >> input_Electrons.txt
ls user.dshoaleh.ttunf8TeV_periodJ.Electrons_*/* >> input_Electrons.txt
ls user.dshoaleh.ttunf8TeV_periodL.Electrons_*/* >> input_Electrons.txt
ls user.dshoaleh.ttunf8TeV_periodA.Muons_*/* > input_Muons.txt
ls user.dshoaleh.ttunf8TeV_periodB.Muons_*/* >> input_Muons.txt
ls user.dshoaleh.ttunf8TeV_periodC.Muons_*/* >> input_Muons.txt
ls user.dshoaleh.ttunf8TeV_periodD.Muons_*/* >> input_Muons.txt
ls user.dshoaleh.ttunf8TeV_periodE.Muons_*/* >> input_Muons.txt
ls user.dshoaleh.ttunf8TeV_periodG.Muons_*/* >> input_Muons.txt
ls user.dshoaleh.ttunf8TeV_periodH.Muons_*/* >> input_Muons.txt
ls user.dshoaleh.ttunf8TeV_periodI.Muons_*/* >> input_Muons.txt
ls user.dshoaleh.ttunf8TeV_periodJ.Muons_*/* >> input_Muons.txt
ls user.dshoaleh.ttunf8TeV_periodL.Muons_*/* >> input_Muons.txt
find  user.dshoaleh.ttunf8TeV_*.BosonJet_*/ -iname "*roo*" > input_BosonJet.txt 
find  user.dshoaleh.ttunf8TeV_*.ttbarBoson_*/ -iname "*roo*" > input_ttbarBoson.txt
find  user.dshoaleh.ttunf8TeV_*.SingleTop_*/ -iname "*roo*" > input_SingleTop.txt
find user.dshoaleh.ttunf8TeV_117050.ttbar_*/ -iname "*roo*" > input_ttbar.txt
ls user.dshoaleh.ttunf8TeV_203054.GStartT_*/* > input_GstartT_MG1000MT600_TWb.txt
ls user.dshoaleh.ttunf8TeV_203057.GStartT_*/* > input_GstartT_MG2000MT1400_TWb.txt
ls user.dshoaleh.ttunf8TeV_203055.GStartT_*/* > input_GstartT_MG1000MT600_TtZ.txt
ls user.dshoaleh.ttunf8TeV_203058.GStartT_*/* > input_GstartT_MG2000MT1400_TtZ.txt
ls user.dshoaleh.ttunf8TeV_203056.GStartT_*/* > input_GstartT_MG1000MT600_TtH.txt
ls user.dshoaleh.ttunf8TeV_203059.GStartT_*/* > input_GstartT_MG2000MT1400_TtH.txt
ls user.dshoaleh.ttunf8TeV*DiBosons_*/*>input_DiBosons.txt
find  user.dshoaleh.ttunf8TeV_*.BosonJet_*/ -iname "*root" > input_Bkg.txt
find  user.dshoaleh.ttunf8TeV_*.ttbarBoson_*/ -iname "*root" >> input_Bkg.txt
find  user.dshoaleh.ttunf8TeV_*.SingleTop_*/ -iname "*root" >> input_Bkg.txt
find  user.dshoaleh.ttunf8TeV_*.BosonJet_*/ -iname "*root.2" > input_Bkg.txt
find  user.dshoaleh.ttunf8TeV_*.ttbarBoson_*/ -iname "*root.2" >> input_Bkg.txt
find  user.dshoaleh.ttunf8TeV_*.SingleTop_*/ -iname "*root.2" >> input_Bkg.txt
find  user.dshoaleh.ttunf8TeV_*.BosonJet_*/ -iname "*root".1 > input_Bkg.txt
find  user.dshoaleh.ttunf8TeV_*.ttbarBoson_*/ -iname "*root.1" >> input_Bkg.txt
find  user.dshoaleh.ttunf8TeV_*.SingleTop_*/ -iname "*root.1" >> input_Bkg.txt
ls user.dshoaleh.ttunf8TeV*DiBosons_*/*>>input_Bkg.txt
