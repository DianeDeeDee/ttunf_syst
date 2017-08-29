#!/bin/bash
#rm MC*.txt
#touch MCZPrime1000.txt
#rm MCZPrime2000.txt
#touch MCZPrime2000.txt


#TAG=Truth_M60D1220_Sup1LjLooseInf3Sj_vtest002
#TAG=Truth_M60D1220_Sup1LjLooseInf3Sj_vtest5p
TAGGtt=Gstartt
TAGMG1000MT600_TWb=GstartT_MG1000MT600_TWb
TAGMG1000MT600_TtH=GstartT_MG1000MT600_TtH
TAGMG1000MT600_TtZ=GstartT_MG1000MT600_TtZ
TAGMG2000MT1400_TWb=GstartT_MG2000MT1400_TWb
TAGMG2000MT1400_TtH=GstartT_MG2000MT1400_TtH
TAGMG2000MT1400_TtZ=GstartT_MG2000MT1400_TtZ
TAGDiBoson=DiBoson
TAGtt117050=ttbarDan
TAGtt117050_0=ttbar_0
TAGtt110404=ttbarSys
TAGttW=ttW
TAGttZ=ttZ
TAGSingleTop=SingleTop
TAGWj=Wj
TAGZj=Zj0
CERN_USER=dshoaleh
SUF="--doParticle 0  --syst 0 --mode 0 --doLoose 0 --data 0"
#DATASUF="--data 1 --mode 0 " #for data
TTSUF="--data 0 --doParticle 0 --mode 0 --syst 1 --doLoose 0  "
TTSUFMSTW="--data 0 --doParticle 0 --mode 0 --syst 1 --doLoose 0 --nPdf=41 --pdfName MSTW2008nlo68cl"
TTSUFCT="--data 0 --doParticle 0 --mode 0 --syst 1 --doLoose 0 --nPdf=53 --pdfName CT10"
TTSUFNNPDF="--data 0 --doParticle 0 --mode 0 --syst 1 --doLoose 0 --nPdf=101 --pdfName NNPDF30_nlo_as_0118"
DATASUF="--data 1 --doParticle 0 --mode 0 --syst 0 --doLoose 1" #for data

TF="--data 0 --doParticle 0 --mode 0 --syst 1 --doLoose 0 --nPdf=53 --pdfName CT10"
python macros/maketar.py

python macros/submit.py - --tag $TAGtt110404 --input datasets/ttbarPowhegPythiaSys.txt --user $CERN_USER --nFilesPerJob 2 --suffix " --data 0 --doParticle 0 --mode 0 --syst 1 --doLoose 0 "  > MCttbar.txt
python macros/submit.py - --tag $TAGGtt --input datasets/KKGluon_ttbar.txt --user $CERN_USER --nFilesPerJob 2 --suffix " $TTSUF "  > MCT_Gtt.txt --destSE 'CA-MCGILL-CLUMEQ-T2_SCRATCHDISK'
python macros/submit.py - --tag $TAGMG1000MT600_TWb --input datasets/GstartT_MG1000MT600_TWb.txt --user $CERN_USER --nFilesPerJob 2 --suffix " $TTSUF "  > MCT_MG1000MT600_TWb.txt
python macros/submit.py - --tag $TAGMG1000MT600_TtH --input datasets/GstartT_MG1000MT600_TtH.txt  --user $CERN_USER --nFilesPerJob 2 --suffix " $TTSUF "  > MCT_MG1000MT600_TtH.txt
python macros/submit.py - --tag $TAGMG1000MT600_TtZ --input datasets/GstartT_MG1000MT600_TtZ.txt  --user $CERN_USER --nFilesPerJob 2 --suffix " $TTSUF " > MCT_MG1000MT600_TtZ.txt
python macros/submit.py - --tag $TAGMG2000MT1400_TWb --input datasets/GstartT_MG2000MT1400_TWb.txt  --user $CERN_USER --nFilesPerJob 2 --suffix " $TTSUF " > MCT_MG2000MT1400_TWb.txt
python macros/submit.py - --tag $TAGMG2000MT1400_TtH --input datasets/GstartT_MG2000MT1400_TtH.txt  --user $CERN_USER --nFilesPerJob 2 --suffix " $TTSUF "> MCT_MG2000MT1400_TtH.txt
python macros/submit.py - --tag $TAGMG2000MT1400_TtZ --input datasets/GstartT_MG2000MT1400_TtZ.txt  --user $CERN_USER --nFilesPerJob 2 --suffix " $TTSUF " > MCT_MG2000MT1400_TtZ.txt:
python macros/submit.py - --tag $TAGtt117050_0 --input datasets/ttbarPowhegPythia.txt --data False --user $CERN_USER --destSE --nFilesPerJob 2 --suffix " $TTSUF " > MCttbar117050_0.txt
python macros/submit.py - --tag $TAGtt110404 --input datasets/ttbarPowhegPythia2.txt --data False --user $CERN_USER --nFilesPerJob 2 --suffix  " $TTSUF " > MCttbar110404.txt
python macros/submit.py - --tag $TAGtt117050 --input datasets/ttbarDanilo.txt --data False --user $CERN_USER  --nFilesPerJob 2 --suffix " $TTSUF " > MCttbar117050.txt
python macros/submit.py - --tag $TAGttZ --input datasets/ttbarZ.txt --data False  --user $CERN_USER --nFilesPerJob 2 --suffix " $TTSUF " > MCttbarZ.txt
python macros/submit.py - --tag $TAGttW --input datasets/ttbarW.txt --data False  --user $CERN_USER --nFilesPerJob 2 --suffix " $TTSUF " > MCttbarW.txt
python macros/submit.py - --tag $TAGZj --input datasets/Zj1.txt --data False --user $CERN_USER  --nFilesPerJob 2 --suffix " $TTSUF " > MCZj.txt
python macros/submit.py - --tag $TAGWj --input datasets/Wj.txt --data False  --user $CERN_USER --nFilesPerJob 2 --suffix " $TTSUF " > MCWj.txt
python macros/submit.py - --tag $TAGDiBoson --input datasets/DiBoson.txt --data False  --user $CERN_USER --nFilesPerJob 2 --suffix " $TTSUF " > MC2Bosons.txt
python macros/submit.py - --tag $TAGSingleTop --input datasets/SingleTop.txt --data False  --user $CERN_USER --nFilesPerJob 2 --suffix " $TTSUF " > MCSingleTop.txt
#--nGBPerJob 5 \
#--destSE CA-MCGILL-CLUMEQ-T2_SCRATCHDISK >> RDO_to_ESD.217999.txt
#--nSkipFiles 16
#--mergeOutput
