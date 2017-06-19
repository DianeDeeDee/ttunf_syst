
################################################################
## In principle all you have to setup is defined in this file ##
################################################################


from configManager import configMgr
from ROOT import kBlack,kWhite,kGray,kRed,kPink,kMagenta,kViolet,kBlue,kAzure,kCyan,kTeal,kGreen,kSpring,kYellow,kOrange
from configWriter import Measurement,Sample
from systematic import Systematic
from math import sqrt
import os

import logger
from logger import Logger

from ROOT import gROOT
gROOT.LoadMacro("./macros/AtlasStyle.C")
import ROOT
ROOT.SetAtlasStyle()

print "testing os get"
print os.getenv('SIG', -1)
print os.getenv('EXCLUDEDSYST', "default")

print os.getenv('ScanUpperLimit', "-1"), "  ", os.getenv('NumPoints', "40")
print os.getenv('SignalScaleFactor', "1")

try:
    theScanInterval=float(os.getenv('ScanUpperLimit', "-1"))
    NumPoints=int(os.getenv('NumPoints', "80"))
except:
    theScanInterval=-1
    NumPoints = 80
        
if NumPoints == -1 :
    NumPoints = 80

if NumPoints < 80 :
    NumPoints = 80

try:
    theSignalScaleFactor=float(os.getenv('SignalScaleFactor', "-1"))
except:
    theSignalScaleFactor=-1


    

try:
    print sigSampleName
except:
    sigSampleName=""
                

try:
    print "channel:", channel
except:
    channel = "Combined"
    #channel = "Boosted"
    #channel = "Resolved"
    print "channel:", channel

try:
    #excludeSyst = os.environ["EXCLUDE_SYST"]
    print "exclude Syst:", excludeSyst
except:
    excludeSyst = "Nominal"
    print "exclude Syst:", excludeSyst

try:
    print "exclude backgrounds:", excludeBackgrounds
except:
    excludeBackgrounds = ""
    print "exclude backgrounds:", excludeBackgrounds

try:
    print "Options:", Options
except:
    if sigSampleName in ["Z400","Z500","Z750","KKg400","KKg500","KKg600","KKg700"]:
        #Options = "Short"
        Options="Clipped"
    else :
        #Options = "None"
        Options= "Clipped"
    print "Default Options:", Options 
                        

if myFitType==FitType.Exclusion:
    print "FitType: Exclusion"

if myFitType==FitType.Discovery:
    print "FitType: Discovery"

if myFitType==FitType.Background:
    print "FitType: Background"
                        
#doSignalScale=True
doSignalScale=False

doSummedBgr=True
#doSummedBgr=False

noSyst = (excludeSyst == "ALLSYST")
doMCStat = True

JESTYPE=4

#doCombinedBtag = True
doCombinedBtag = False


    
                    
#theScanInterval=-1



#if not noSyst : 
#try:
#    theScanInterval=ScanUpperLimits[sigSampleName_temp]
#except:
#    theScanInterval=-1

#if theScanInterval == -1 :
#    print "scan range for signal sample not found, checking for similar signals"
#    sigSampleName_temp=sigSampleName.split("_")[0]
#    sigSampleName_temp=sigSampleName_temp.replace("RSG","KKg")
#    print "temporary signSampleName: ", sigSampleName_temp
#    try:
#        theScanInterval=ScanUpperLimits[sigSampleName_temp]
#    except:
#        theScanInterval=-1
                    


#this is a terrible idea, don't have a better solution atm... I cant use the same scan ranges on syst and stat only limits on these points
#if noSyst and channel == "Boosted" :
#    if sigSampleName == "Z750" : theScanInterval=5
#    #if sigSampleName == "Z1250" : theScanInterval=.4
#    if sigSampleName == "KKg1000" : theScanInterval=1.5
#    if sigSampleName == "KKg1300" : theScanInterval=0.7
    
    


#try:
#    if doSignalScale : theSignalScaleFactor=SignalScaleFactor[sigSampleName]
#    else : theSignalScaleFactor=-1
#except:
#    theSignalScaleFactor=-1
                

if doSignalScale and not (theSignalScaleFactor == -1): #and not (theScanInterval == -1) :
    theScanInterval=theScanInterval/theSignalScaleFactor
    #theScanInterval=theScanInterval/theSignalScaleFactor

print "the scan interval is : ", theScanInterval
print "the number of scanning points is : ", NumPoints
print "the signal scale factor is : ", theSignalScaleFactor
                            

print "Initial Scan Maximum: ", theScanInterval



log = Logger("MyShapeFitExample")
#log.setLevel(logger.WARNING)
#log.setLevel(logger.INFO) #should have no effect if -L is used
#log.warning("example warning from python")
#log.error("example error from python")


#-------------------------------
# Parameters for hypothesis test
#-------------------------------
#configMgr.doHypoTest=False
configMgr.nTOYs=5000
configMgr.calculatorType=2 # 2=asymptotic calculator, 0=frequentist calculator
configMgr.testStatType=3   # 3=one-sided profile likelihood test statistic (LHC default)
#configMgr.nPoints=200       # number of values scanned of signal-strength for upper-limit determination of signal strength.
configMgr.nPoints=NumPoints  # number of values scanned of signal-strength for upper-limit determination of signal strength.
configMgr.InitialScanInterval=theScanInterval

configMgr.useCacheToTreeFallback = True # enable the fallback to trees
configMgr.useHistBackupCacheFile = True # enable the use of an alternate data file
configMgr.histBackupCacheFile =  "foo.root" # the data file of your background fit (= the backup cache file)


try :
    print "fix mu: ", fixMu
except :
    fixMu = False

try :
    print "no mu: ", noMu
except :
    noMu = False

#sfx = ""


# try if signal name set from command line (with -c "sigSampleName='zprime1500_tt';" ):
try :
    print "signal: ", sigSampleName
except :
    print "no signal set, do bkg only"
    sigSampleName = ""
sfx=sigSampleName

# try if splitting set from command line:
try : 
    print "type of b-tag systematics: ", btagSystematicsType
except :
    btagSystematicsType = "default"

if btagSystematicsType == "splitCategory" :
    sfx = sfx + "_splitBtag"
elif btagSystematicsType == "splitPt" :
    sfx = sfx + "_splitBtagPt"
elif btagSystematicsType == "splitMtt" :
    sfx = sfx + "_splitBtagMtt"
elif btagSystematicsType == "eigenvector" :
    sfx = sfx + "_btagEigenvec"
else :
    print "use default b-tagging systematics"


try :
    print "data sample:", dataSampleName
    sfx = sfx + "_" + dataSampleName
except :
    #print "no data sample set, use",
    #dataSampleName = "DataNom"
    dataSampleName = "BgrNom"
    #dataSampleName = "PseudoDatabtagHigh"
    print "no data sample set, use",dataSampleName


doInclusive = False

srPrefix = "SRexc"
if doInclusive :
    srPrefix = "SRinc"
    sfx = sfx + "_Inclusive"
            
try : 
    print "btag bins to use: ", btag_bins
except :
    btag_bins=["btag_cat1","btag_cat2","btag_cat3"]
if doInclusive :
    btag_bins=[""]
else :
    print "btag bins to use: ", btag_bins
    sfx = sfx + "_ExclusiveCat"
    for btagBin in btag_bins :
        if btagBin.find("btag_cat") > -1 :
            sfx = sfx + "c"
        sfx = sfx + btagBin[-1]

if fixMu :
    sfx = sfx + "_fixedMu"
elif noMu :
    sfx = sfx + "_noMu"
    

if "Short" in Options:
    inputFile="short_"+inputFile
#if "Clipped" in Options:
#    inputFile="clipped_"+inputFile
        
print "input file: ", inputFile

                


#https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TopSystematicUncertainties
#wjetsNormUncert = 0.2 # 10-20% when using data driven SF, depending on jet bin; MC only (4-jets bin): 48%
#ttNormUncert = 0.1
#stopNormUncert = 0.05 # 4-7% depending on process
#othersNormUncert = 0.4 # z+jets: MC only (4-jets bin): 48%

configMgr.analysisName = "%s_%s_%s" % ("MassSpectra_Paper_Boosted",sfx,excludeSyst)
configMgr.histCacheFile = "data/"+inputFile
configMgr.outputFileName = "results/"+configMgr.analysisName+"_Output.root"

# Scaling calculated by outputLumi / inputLumi
configMgr.inputLumi = 1.0 # Luminosity of input TTree after weighting
configMgr.outputLumi = 1.0 # Luminosity required for output histograms
configMgr.setLumiUnits("fb-1")

print "A:  configMgr.cutsDict[...]"
for btagBin in btag_bins :
    if (channel == "Combined" or channel == "Boosted"):
        configMgr.cutsDict[srPrefix+btagBin+"ElectronsBoosted"] = "1.0"
        configMgr.cutsDict[srPrefix+btagBin+"MuonsBoosted"] = "1.0"
    if (channel == "Combined" or channel == "Resolved"):
        configMgr.cutsDict[srPrefix+btagBin+"ElectronsResolved"] = "1.0"
        configMgr.cutsDict[srPrefix+btagBin+"MuonsResolved"] = "1.0"
            
ttSample = Sample("tt",kWhite)
if not (noSyst or excludeSyst == "Norm_tt"):
    ##ttSample.addSystematic(Systematic("Norm",None,1.098,0.894,"user","overallHistoSys"))
    ttSample.addSystematic(Systematic("Norm_ttbar","Nom","Norm_ttbarHigh","Norm_ttbarLow","tree","overallSys"))
ttSample.setStatConfig(doMCStat)
ttSample.setNormByTheory(True)

WjetsSample = Sample("W",kOrange)
##if not (noSyst or excludeSyst == "NormW"):
##    WjetsSample.addSystematic(Systematic("NormW","Nom","NormWHigh","NormWLow","tree","overallSys"))
WjetsSample.setStatConfig(doMCStat)
WjetsSample.setNormByTheory(False)

stopSample = Sample("single-top",kBlue)
#if not (noSyst or excludeSyst == "Norm"):
#    stopSample.addSystematic(Systematic("singleTop_xsec",None,1.077,0.923,"user","overallSys"))
stopSample.setStatConfig(doMCStat)
stopSample.setNormByTheory(True)

dibosonSample = Sample("Diboson",kYellow)
#if not (noSyst or excludeSyst == "Norm"):
#   dibosonSample.addSystematic(Systematic("Norm_Diboson_","Nom","Norm_Diboson_High","Norm_Diboson_Low","tree","overallSys"))
dibosonSample.setStatConfig(doMCStat)
dibosonSample.setNormByTheory(True)

zjetsSample = Sample("Z",kYellow)
#if not (noSyst or excludeSyst == "Norm"):
    #zjetsSample.addSystematic(Systematic("Norm_Zjets_","Nom","Norm_Zjets_High","Norm_Zjets_Low","tree","overallSys"))
zjetsSample.setStatConfig(doMCStat)
zjetsSample.setNormByTheory(True)

#othersSample = Sample("others",kYellow)
#if not (noSyst or excludeSyst == "Norm"):
    #othersSample.addSystematic(Systematic("Norm_others_","Nom","Norm_others_High","Norm_others_Low","tree","overallSys"))
#othersSample.setStatConfig(doMCStat)
#othersSample.setNormByTheory(True)


ttVSample =  Sample("ttV",kYellow)
ttVSample.setStatConfig(doMCStat)
ttVSample.setNormByTheory(True)

smallBgrSample =  Sample("smallBgr",kBlue)
smallBgrSample.setStatConfig(doMCStat)
smallBgrSample.setNormByTheory(True)


dataSample = Sample(dataSampleName,kBlack)
#dataSample.setData()
#dataSample = Sample("BgrNom",kBlack)
#dataSample = Sample("datahistoNom",kBlack)
dataSample.setData()

commonSamples = []
 
qcdSample = Sample("QCD", kViolet)
##if not (noSyst or excludeSyst == "Norm_QCDe"):
##    qcdSample.addSystematic(Systematic("Norm_QCDe","Nom","Norm_QCDeHigh","Norm_QCDeLow","tree","overallSys"))
##if not (noSyst or excludeSyst == "Norm_QCDmu"):
##    qcdSample.addSystematic(Systematic("Norm_QCDmu","Nom","Norm_QCDmuHigh","Norm_QCDmuLow","tree","overallSys"))
qcdSample.setStatConfig(doMCStat)
qcdSample.setNormByTheory(False)


commonSamples += [ttSample]
if "wjets" not in excludeBackgrounds:
    commonSamples += [WjetsSample]
    
if doSummedBgr :
    commonSamples += [smallBgrSample]
            

else :
    if "ttV" not in excludeBackgrounds:
        commonSamples += [ttVSample]
    if "qcd" not in excludeBackgrounds:
        commonSamples += [qcdSample]
    if "vv" not in excludeBackgrounds:
        commonSamples += [dibosonSample]
    if "zjets" not in excludeBackgrounds:
        commonSamples += [zjetsSample]
    if "stop" not in excludeBackgrounds:
        commonSamples += [stopSample]

                    


commonSamples.append(dataSample)

#**************
# Exclusion fit
#**************

#Fit config instance
print "A :configMgr.addFitConfig(MyFitConfig)"
exclusionFitConfig = configMgr.addFitConfig("MyFitConfig")
exclusionFitConfig.statErrThreshold=0.05
meas=exclusionFitConfig.addMeasurement(name="NormalMeasurement",lumi=1.0,lumiErr=0.028)

#if not myFitType==FitType.Discovery:
#    meas.addPOI("mu_SIG")
meas.addPOI("mu_SIG")  #what does this mean if there is no signal?

meas.addParamSetting("Lumi",True,1)  #true= constant, false= floating
#meas.addParamSetting("mu_tt",True,1.0246e+00)
#meas.addParamSetting("mu_Wjets",True,8.7260e-01)

#Samples
print "A :configMgr.addSamples(commonSamples)"
exclusionFitConfig.addSamples(commonSamples)

sigSampleNameList = []
if not sigSampleName=="" :
    sigSample = Sample(sigSampleName,kPink)
    #sigSample.setNormByTheory()
    sigSample.setNormByTheory(True) # SF: test!
    sigSample.setStatConfig(doMCStat)# SF: test! -> no influence on p-val here!
#    sigSample.setNormFactor("mu_SIG",1.,0.,5.)
    sigSample.setNormFactor("mu_SIG",0 ,0.,5.)                    
    if doSignalScale and not (theSignalScaleFactor == -1) :
        print "applying signal scale factor"
        sigSample.addNormFactor("Scale",theSignalScaleFactor,theSignalScaleFactor,theSignalScaleFactor,True)
    exclusionFitConfig.addSamples(sigSample)
    exclusionFitConfig.setSignalSample(sigSample)
    sigSampleNameList = [sigSampleName]




nBins_boosted=15
nBins_resolved=17
mttMaxBoosted=1.0
mttMaxResolved=1.0

if "Short" in Options:
    print "using Short Config"
    nBins_boosted=7
    nBins_resolved=9
    mttMaxBoosted=0.45
    mttMaxResolved=0.52

if "Clipped" in Options:
    print "using clipped Config"
    nBins_boosted=15
    nBins_resolved=15
    mttMaxBoosted=1.0
    mttMaxResolved=0.88
                        

print "nbins_boosted:",nBins_boosted
print "nbins_resolved:",nBins_resolved
print "mttMaxBoosted:",mttMaxBoosted
print "mttMaxResolved:",mttMaxResolved



#Channels
sRegs = []
sRegsBoosted = []
sRegsResolved = []
sRegsDict = {}
print "A :  exclusionFitConfig.addChannel(...)"
for btagBin in btag_bins :
    if (channel == "Combined" or channel == "Boosted"):
        myChannel = exclusionFitConfig.addChannel("cuts",[srPrefix+btagBin+"ElectronsBoosted"],nBins_boosted,0.,mttMaxBoosted)
        myChannel.useOverflowBin=False
        myChannel.useUnderflowBin=False
        sRegsDict[srPrefix+btagBin+"ElectronsBoosted"] = myChannel
        sRegs.append(myChannel)
        sRegsBoosted.append(myChannel)
        myChannel = exclusionFitConfig.addChannel("cuts",[srPrefix+btagBin+"MuonsBoosted"],nBins_boosted,0.,mttMaxBoosted)
        myChannel.useOverflowBin=False
        myChannel.useUnderflowBin=False
        sRegsDict[srPrefix+btagBin+"MuonsBoosted"] = myChannel
        sRegs.append(myChannel)
        sRegsBoosted.append(myChannel)
    if (channel == "Combined" or channel == "Resolved"):
        myChannel = exclusionFitConfig.addChannel("cuts",[srPrefix+btagBin+"ElectronsResolved"],nBins_resolved,0.,mttMaxResolved)
        #myChannel = exclusionFitConfig.addChannel("cuts",[srPrefix+btagBin+"ElectronsResolved"],nBins_resolved,0.,0.75)
        myChannel.useOverflowBin=False
        myChannel.useUnderflowBin=False
        sRegsDict[srPrefix+btagBin+"ElectronsResolved"] = myChannel
        sRegs.append(myChannel)
        sRegsResolved.append(myChannel)
        myChannel = exclusionFitConfig.addChannel("cuts",[srPrefix+btagBin+"MuonsResolved"],nBins_resolved,0.,mttMaxResolved)
        #myChannel = exclusionFitConfig.addChannel("cuts",[srPrefix+btagBin+"MuonsResolved"],nBins_resolved,0.,0.75)
        myChannel.useOverflowBin=False
        myChannel.useUnderflowBin=False
        sRegsDict[srPrefix+btagBin+"MuonsResolved"] = myChannel
        sRegs.append(myChannel)
        sRegsResolved.append(myChannel)
        
exclusionFitConfig.setSignalChannels(sRegs)
#allBgr=["tt","Diboson","Z","single-top", "W"]
allBgr=["tt"]

if "wjets" not in excludeBackgrounds:
    allBgr += ["W"]

if doSummedBgr :
    allBgr += ["smallBgr"]

else :
    if "ttV" not in excludeBackgrounds:
        allBgr += ["ttV"]

    if "qcd" not in excludeBackgrounds:
        allBgr += ["QCD"]
    
    if "vv" not in excludeBackgrounds:
        allBgr += ["Diboson"]
    if "zjets" not in excludeBackgrounds:
        allBgr += ["Z"]
    if "stop" not in excludeBackgrounds:
        allBgr += ["single-top"]

                                
        

#print "A: addSystematics "
#if not (noSyst or excludeSyst == "JER"):
    #for mySampleName in (allBgr + sigSampleNameList ) :
    ##for mySampleName in (allBgr  ) :
        #for myChannel in sRegs :
            #print "A: JER  ", mySampleName
            #myChannel.getSample(mySampleName).addSystematic(Systematic("JER","Nom","JERHigh","JERLow","tree","overallHistoSys"))
            
            
print "A: addSystematics JVF"
if not (noSyst or excludeSyst == "JVF"):
    for mySampleName in (allBgr + sigSampleNameList ) :
    #for mySampleName in (allBgr  ) :
        for myChannel in sRegs :
            print "A: JVF  ", mySampleName
            myChannel.getSample(mySampleName).addSystematic(Systematic("JVF","Nom","JVFHigh","JVFLow","tree","overallHistoSys"))
            
            
print "A: addSystematics PDF"
if not (noSyst or excludeSyst == "PDF"):
    for mySampleName in (["tt"] + sigSampleNameList) :
    #for mySampleName in (allBgr + sigSampleNameList ) :
    #for mySampleName in (allBgr  ) :
        for myChannel in sRegs :
            print "A: PDF  ", mySampleName
            myChannel.getSample(mySampleName).addSystematic(Systematic("PDF","Nom","PDFHigh","PDFLow","tree","overallHistoSys"))



print "A: addSystematics BoostedJES0"
if not (noSyst or "BoostedJES0" in excludeSyst):
    for mySampleName in (allBgr + sigSampleNameList) :
    #for mySampleName in (allBgr ) :
        #for myChannel in sRegs :
        for myChannel in sRegsBoosted :
            myChannel.getSample(mySampleName).addSystematic(Systematic("BoostedJES0","Nom","BoostedJES0High","BoostedJES0Low","tree","overallHistoSys"))

print "A: addSystematics BoostedJES13"
if not (noSyst or "BoostedJES13" in excludeSyst):
    for mySampleName in (allBgr + sigSampleNameList) :
        #for myChannel in sRegs :
        for myChannel in sRegsBoosted :
            myChannel.getSample(mySampleName).addSystematic(Systematic("BoostedJES13","Nom","BoostedJES13High","BoostedJES13Low","tree","overallHistoSys"))

print "A: addSystematics BoostedJES14"
if not (noSyst or "BoostedJES14" in excludeSyst):
    for mySampleName in (allBgr + sigSampleNameList) :
        #for myChannel in sRegs :
        for myChannel in sRegsBoosted :
            myChannel.getSample(mySampleName).addSystematic(Systematic("BoostedJES14","Nom","BoostedJES14High","BoostedJES14Low","tree","overallHistoSys"))

print "A: addSystematics BoostedJES15"
if not (noSyst or "BoostedJES15" in excludeSyst):
    for mySampleName in (allBgr + sigSampleNameList) :
        #for myChannel in sRegs :
        for myChannel in sRegsBoosted :
            myChannel.getSample(mySampleName).addSystematic(Systematic("BoostedJES15","Nom","BoostedJES15High","BoostedJES15Low","tree","overallHistoSys"))

print "A: addSystematics BoostedJES16"
if not (noSyst or "BoostedJES16" in excludeSyst):
    for mySampleName in (allBgr + sigSampleNameList) :
        #for myChannel in sRegs :
        for myChannel in sRegsBoosted :
            myChannel.getSample(mySampleName).addSystematic(Systematic("BoostedJES16","Nom","BoostedJES16High","BoostedJES16Low","tree","overallHistoSys"))
                                
print "A: addSystematics BoostedJER"
if not (noSyst or "BoostedJER" in excludeSyst):
    for mySampleName in (allBgr + sigSampleNameList) :
        #for myChannel in sRegs :
        for myChannel in sRegsBoosted :
            myChannel.getSample(mySampleName).addSystematic(Systematic("BoostedJER","Nom","BoostedJERHigh","BoostedJERLow","tree","overallHistoSys"))
            
print "A: addSystematics BoostedJMR"
if not (noSyst or "BoostedJMR" in excludeSyst):
    for mySampleName in (allBgr + sigSampleNameList) :
        #for myChannel in sRegs :
        for myChannel in sRegsBoosted :
            myChannel.getSample(mySampleName).addSystematic(Systematic("BoostedJMR","Nom","BoostedJMRHigh","BoostedJMRLow","tree","overallHistoSys"))
            
print "A: addSystematics BoostedJMS"
if not (noSyst or "BoostedJMS" in excludeSyst):
    for mySampleName in (allBgr + sigSampleNameList) :
        #for myChannel in sRegs :
        for myChannel in sRegsBoosted :
            myChannel.getSample(mySampleName).addSystematic(Systematic("BoostedJMS","Nom","BoostedJMSHigh","BoostedJMSLow","tree","overallHistoSys"))
                                


if doCombinedBtag :
    print "A: addSystematics btag"
    if not (noSyst or excludeSyst == "btag"):
        #for mySampleName in (allBgr ) :
        for mySampleName in (allBgr + sigSampleNameList) :
            for myChannel in sRegs :
                myChannel.getSample(mySampleName).addSystematic(Systematic("btag","Nom","btagHigh","btagLow","tree","overallHistoSys"))

else :
    print "A: addSystematics btag componenets"
    # for N in range(0, 22):
    #for N in [0,1,2,3,4,5,6,7,8,9,10] :
    for N in [7,8,9,10] :
        btagN="btag"+str(N)
        if not (noSyst or "btag" in excludeSyst ):
            #for mySampleName in (allBgr ) :
            for mySampleName in (allBgr  + sigSampleNameList):
                for myChannel in sRegs :
                    myChannel.getSample(mySampleName).addSystematic(Systematic(btagN,"Nom",btagN+"High",btagN+"Low","tree","overallHistoSys"))

    #if not (noSyst or "btag" in excludeSyst ):
    #    #for mySampleName in (allBgr ) :
    #    for mySampleName in (allBgr  + sigSampleNameList):
    #        for myChannel in sRegs :
    #            myChannel.getSample(mySampleName).addSystematic(Systematic("smallBtag","Nom","smallBtagHigh","smallBtagLow","tree","overallHistoSys"))
                                                
    

print "A: addSystematics ctag"
if not (noSyst or excludeSyst == "ctag"):
    for mySampleName in  (allBgr ):
    #for mySampleName in  (allBgr + sigSampleNameList):
        for myChannel in sRegs :
            myChannel.getSample(mySampleName).addSystematic(Systematic("ctag","Nom","ctagHigh","ctagLow","tree","overallHistoSys"))
#syst2.5

print "A: addSystematics mistag"
if not (noSyst or excludeSyst == "mistag"):
    for mySampleName in  (allBgr ):
    #for mySampleName in  (allBgr  + sigSampleNameList):
        for myChannel in sRegs :
            myChannel.getSample(mySampleName).addSystematic(Systematic("mistag","Nom","mistagHigh","mistagLow","tree","overallHistoSys"))

print "A: addSystematics IFSR"
if not (noSyst or excludeSyst == "IFSR"):
  ### tt only:
  for myChannel in sRegs :
    myChannel.getSample("tt").addSystematic(Systematic("IFSR","Nom","IFSRHigh","IFSRLow","tree","overallHistoSys"))
#syst 3

print "A: addSystematics MC_GEN"
if not (noSyst or "MC_Gen" in excludeSyst):
##if False:
  ##tt only:
  for myChannel in sRegs :
    myChannel.getSample("tt").addSystematic(Systematic("MC_Gen","Nom","MC_GenHigh","MC_GenLow","tree","overallHistoSys"))

print "A: addSystematics PS"
if not (noSyst or excludeSyst == "PS"):
    for myChannel in sRegs :
        myChannel.getSample("tt").addSystematic(Systematic("PS","Nom","PSHigh","PSLow","tree","overallHistoSys"))



#print "A: addSystematics EWS"
#if not (noSyst or "EWS" in excludeSyst):
## if False:
    ## tt only
#    for mySampleName in ["tt"] :
#        for myChannel in sRegs :
        #for myChannel in sRegsBoosted :
#            myChannel.getSample(mySampleName).addSystematic(Systematic("EWS","Nom","EWSHigh","EWSLow","tree","overallHistoSys"))


#resolved channel Z750- pickes excessively large range for first scan, fails
print "A: addSystematics top_mass"
if not (noSyst or excludeSyst == "top_mass"):
    for myChannel in sRegs :
    #for myChannel in sRegsBoosted :
        myChannel.getSample("tt").addSystematic(Systematic("top_mass","Nom","top_massHigh","top_massLow","tree","overallHistoSys"))
#syst 4

#print "A: addSystematics MuonSF"
#if not (noSyst or  "MuonSF" in excludeSyst):
#    for mySampleName in  (allBgr + sigSampleNameList):
#        for myChannel in sRegs :
#            myChannel.getSample(mySampleName).addSystematic(Systematic("MuonSF","Nom","MuonSFHigh","MuonSFLow","tree","overallHistoSys"))

#print "A: addSystematics ElectronSF"
#if not (noSyst or  "ElectronSF" in excludeSyst):
#    for mySampleName in  (allBgr  + sigSampleNameList):
#        for myChannel in sRegs :
#            myChannel.getSample(mySampleName).addSystematic(Systematic("ElectronSF","Nom","ElectronSFHigh","ElectronSFLow","tree","overallHistoSys"))
                        
print "A: addSystematics luminosity"
if not (noSyst or excludeSyst == "luminosity"):
    for mySampleName in  (allBgr + sigSampleNameList):
    #for mySampleName in  (allBgr):
        for myChannel in sRegs :
            myChannel.getSample(mySampleName).addSystematic(Systematic("luminosity","Nom","luminosityHigh","luminosityLow","tree","overallSys"))
                        
if JESTYPE == 1 :
    print "A: addSystematics JESALL"
    if not (noSyst or excludeSyst == "JESALL"):
        for mySampleName in  (allBgr + sigSampleNameList):
        #for mySampleName in  (allBgr ):
            for myChannel in sRegs :
                myChannel.getSample(mySampleName).addSystematic(Systematic("JESALL","Nom","JESALLHigh","JESALLLow","tree","overallHistoSys"))



elif JESTYPE == 4 :
    print "A: addSystematics small JES "
    if not (noSyst or  "JES" in excludeSyst or  "JESsmall" in excludeSyst):
        #for mySampleName in  (allBgr ):
        for mySampleName in  (allBgr + sigSampleNameList):
            for myChannel in sRegs :
                myChannel.getSample(mySampleName).addSystematic(Systematic("JESsmall","Nom","JESsmallHigh","JESsmallLow","tree","overallHistoSys"))
                
    print "A: addSystematics JES"
    # for N in range(0, 22):
    #for N in [3,7,12,18,20,21,22]:
    #for N in [3,7,12,20,22]:
    for N in [3,7,12,20]: 
        jesN="JES"+str(N)
        if not (noSyst or excludeSyst == jesN or excludeSyst == "JES"):
            #    ##for mySampleName in (["tt", "Wjets", "others"] + sigSampleNameList) :
            #for mySampleName in  (allBgr ):
            for mySampleName in (allBgr + sigSampleNameList) :
                #for myChannel in sRegsResolved :
                for myChannel in sRegs :
                    myChannel.getSample(mySampleName).addSystematic(Systematic(jesN,"Nom",jesN+"High",jesN+"Low","tree","overallHistoSys"))

elif JESTYPE == 22 :
    print "A: addSystematics JES"
    # for N in range(0, 22):
    for N in [0,1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,20,21,22]:
        jesN="JES"+str(N)
        if not (noSyst or excludeSyst == jesN or excludeSyst == "JES"):
            #    ##for mySampleName in (["tt", "Wjets", "others"] + sigSampleNameList) :
            for mySampleName in  (allBgr ):
            #for mySampleName in (allBgr  + sigSampleNameList):
                for myChannel in sRegs :
                    myChannel.getSample(mySampleName).addSystematic(Systematic(jesN,"Nom",jesN+"High",jesN+"Low","tree","overallHistoSys"))
                                                                     


print "A: finished creating systematic list"


############below this line is garbage

#if not (noSyst or excludeSyst == "PDF"):
#    for mySampleName in (["tt","Diboson","Z","single-top", "W"] + sigSampleNameList) :
#        for myChannel in sRegs :
#            myChannel.getSample(mySampleName).addSystematic(Systematic("PDF","Nom","PDFHigh","PDFLow","tree","overallHistoSys"))


 
            
            


#if not (noSyst or excludeSyst == "JES1"):
#    ##for mySampleName in (["tt", "Wjets", "others"] + sigSampleNameList) :
#    for mySampleName in (["tt","Diboson","Z","single-top", "W"] + sigSampleNameList) :
#        for myChannel in sRegs :
#            myChannel.getSample(mySampleName).addSystematic(Systematic("JES1","Nom","JES1High","JES1Low","tree","overallHistoSys"))




#if not (noSyst or excludeSyst == "Ptjmin10"):
  #for myChannel in sRegs :
    #myChannel.getSample("W").addSystematic(Systematic("Ptjmin10","Nom","Ptjmin10High","Ptjmin10Low","tree","overallHistoSys"))

#if not (noSyst or excludeSyst == "Iqopt3"):
  #for myChannel in sRegs :
    #myChannel.getSample("W").addSystematic(Systematic("Iqopt3","Nom","Iqopt3High","Iqopt3Low","tree","overallHistoSys"))

#if not (noSyst or excludeSyst == "WHFC3"):
  #for myChannel in sRegs :
    #myChannel.getSample("W").addSystematic(Systematic("WHFC3","Nom","WHFC3High","WHFC3Low","tree","overallHistoSys"))

#if not (noSyst or excludeSyst == "WHFC4"):
  #for myChannel in sRegs :
    #myChannel.getSample("W").addSystematic(Systematic("WHFC4","Nom","WHFC4High","WHFC4Low","tree","overallHistoSys"))

#if not (noSyst or excludeSyst == "WHFC0"):
  #for myChannel in sRegs :
    #myChannel.getSample("W").addSystematic(Systematic("WHFC0","Nom","WHFC0High","WHFC0Low","tree","overallHistoSys"))

######################################################################################################
#if not (noSyst or excludeSyst == "JEff"):
    #for mySampleName in [] :
        #for myChannel in sRegs :
            #myChannel.getSample(mySampleName).addSystematic(Systematic("JEff","Nom","JEffHigh","JEffLow","tree","overallHistoSys"))


##if not (noSyst or excludeSyst == "MUSC"):
### if False:
    ###for mySampleName in ["tt", "Diboson", "stop", "Wjets", "Zjets"] :
    ### exlude others:
    ##for mySampleName in (["tt", "stop", "Wjets"] + sigSampleNameList) :
        ##for myChannel in sRegs :
            ##myChannel.getSample(mySampleName).addSystematic(Systematic("MUSC","Nom","MUSCHigh","MUSCLow","tree","overallHistoSys"))

#if not (noSyst or excludeSyst == "muonSF"):
## if False:
    ##for mySampleName in ["tt", "Diboson", "stop", "Wjets", "Zjets"] :
    ## exlude others:
    ##for mySampleName in (["tt", "others", "Wjets"] + sigSampleNameList) :
    #for mySampleName in (["tt","Diboson","Zjets","stop", "Wjets"] + sigSampleNameList) :
        #sRegs[1].getSample(mySampleName).addSystematic(Systematic("muonSF","Nom","muonSFHigh","muonSFLow","tree","overallHistoSys"))

#if not (noSyst or excludeSyst == "electronSF"):
## if False:
    ##for mySampleName in ["tt", "Diboson", "stop", "Wjets", "Zjets"] :
    ## exlude others:
    ##for mySampleName in (["tt", "others", "Wjets"] + sigSampleNameList) :
    #for mySampleName in (["tt","Diboson","Zjets","stop", "Wjets"] + sigSampleNameList) :
        #sRegs[0].getSample(mySampleName).addSystematic(Systematic("electronSF","Nom","electronSFHigh","electronSFLow","tree","overallHistoSys"))

##if not (noSyst or excludeSyst == "EER"):
### if False:
    ###for mySampleName in ["tt", "Diboson", "stop", "Wjets", "Zjets"] :
    ### exlude others:
    ##for mySampleName in (["tt", "stop", "Wjets"] + sigSampleNameList) :
        ##for myChannel in sRegs :
            ##myChannel.getSample(mySampleName).addSystematic(Systematic("EER","Nom","EERHigh","EERLow","tree","overallHistoSys"))

##if not (noSyst or excludeSyst == "EES"):
### if False:
    ###for mySampleName in ["tt", "Diboson", "stop", "Wjets", "Zjets"] :
    ### exlude others:
    ##for mySampleName in (["tt", "stop", "Wjets"] + sigSampleNameList) :
        ##for myChannel in sRegs :
            ##myChannel.getSample(mySampleName).addSystematic(Systematic("EES","Nom","EESHigh","EESLow","tree","overallHistoSys"))

##if not (noSyst or excludeSyst == "MET_SJS"):
### if False:
    ###for mySampleName in ["tt", "Diboson", "stop", "Wjets", "Zjets"] :
    ### exlude others:
    ##for mySampleName in (["tt", "stop", "Wjets"] + sigSampleNameList) :
        ##for myChannel in sRegs :
            ##myChannel.getSample(mySampleName).addSystematic(Systematic("MET_SJS","Nom","MET_SJSHigh","MET_SJSLow","tree","overallHistoSys"))

##if not (noSyst or excludeSyst == "MET_RSJ"):
### if False:
    ###for mySampleName in ["tt", "Diboson", "stop", "Wjets", "Zjets"] :
    ### exlude others:
    ##for mySampleName in (["tt", "stop", "Wjets"] + sigSampleNameList) :
        ##for myChannel in sRegs :
            ##myChannel.getSample(mySampleName).addSystematic(Systematic("MET_RSJ","Nom","MET_RSJHigh","MET_RSJLow","tree","overallHistoSys"))

##if not (noSyst or excludeSyst == "muonSmear"):
### if False:
    ###for mySampleName in ["tt", "Diboson", "stop", "Wjets", "Zjets"] :
    ### exlude others:
    ##for mySampleName in (["tt", "stop", "Wjets"] + sigSampleNameList) :
        ##for myChannel in sRegs :
            ##myChannel.getSample(mySampleName).addSystematic(Systematic("muonSmear","Nom","muonSmearHigh","muonSmearLow","tree","overallHistoSys"))

## Define measurement
#meas = ana.addMeasurement(name="NormalMeasurement",lumi=1.0,lumiErr=lumiError)
#meas.addPOI("mu_Sig")

