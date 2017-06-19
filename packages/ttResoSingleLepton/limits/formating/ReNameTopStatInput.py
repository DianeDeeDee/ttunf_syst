from ROOT import *
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-m", "--mode", action="store", dest="mode", type="int", default=0, 
                                    help="Run Mode: 1=Boosted El; 2=Resolved El; 3=Boosted Mu;  4=Resolved Mu; 5=Boosted Both; 6=Resolved Both; 0=All", metavar="MODE")
(options, args) = parser.parse_args()

######## author: Dominik Duda (dominik.duda@cern.ch)
###############################################################################################
#############                      Only change this lines !!!!

filename = "MttSpectra_Oct06.root";
lepton_channel="Muons"
#lepton_channel="Electrons"
isBoosted=true
ReScaleXAxis=True
###############################################################################################
print filename

for lepton_channel in ["Muons", "Electrons", "Leptons"]:
    if lepton_channel=="Electrons" and ( options.mode!=1 and options.mode!=2 and options.mode!=0):  continue;
    if lepton_channel=="Muons"     and ( options.mode!=3 and options.mode!=4 and options.mode!=0):  continue;
    if lepton_channel=="Leptons"   and ( options.mode!=5 and options.mode!=6 and options.mode!=0):  continue;
    for isBoosted in [true, false]:
        if isBoosted==True  and ( options.mode==2 or options.mode==4 or options.mode==6):  continue;
        if isBoosted==False and ( options.mode==1 or options.mode==3 or options.mode==5):  continue;
        print lepton_channel
        if isBoosted==True:
            BoostedOrResolved="Boosted"
            mass_ttbar="masstT"
        elif isBoosted==False:
            BoostedOrResolved="Resolved"
            mass_ttbar="massTTbarChi2LPC"
        if lepton_channel=="Muons":
            lep_channel="mu"
        elif lepton_channel=="Electrons":
            lep_channel="e"
        elif lepton_channel=="Leptons":
            lep_channel="lep"

        def CopyHisto(histo):
            h2 = histo.Clone()
            hname=histo.GetName()
            h2.SetName(hname.replace("_obs_cuts","_obs_cutsNorm"))
            h2.SetTitle(hname.replace("_obs_cuts","_obs_cutsNorm"))
            hlist.Add(h2)

        def SetRangeofXaxis(histo):
            if ReScaleXAxis==False :
                return histo
            hname=histo.GetName()
            histo.SetName("temp_name")
            if isBoosted :
                n_unfilled_bins=5
                reduced_bins=0
            else :
                n_unfilled_bins=3
                reduced_bins=2
            bins=histo.GetNbinsX()-reduced_bins
            xbins=histo.GetXaxis().GetXbins().GetArray()
            for i in range(0,bins+1) :
                xbins[i]=float(i)/float(bins-n_unfilled_bins)
               # print xbins[i], "  ", i
            h = TH1F(hname,"", bins-n_unfilled_bins,xbins)
            h.Sumw2()
            for i in range(n_unfilled_bins+1, bins+1):
                _bc = max(histo.GetBinContent(i), 0.)
                _be = histo.GetBinError(i)
                if _bc and _be == 0: 
                    _be = sqrt(histo.GetBinContent(i))
                #print "setting bin ", i 
                h.SetBinContent(i-n_unfilled_bins, _bc)
                h.SetBinError(i-n_unfilled_bins, _be)
            return h 


        Syst_list={
            "EleSFHigh":"ElectronSFHigh",
            "EleSFLow":"ElectronSFLow",
            "MuSFHigh":"MuonSFHigh",
            "MuSFLow":"MuonSFLow",

            "smallBtagHigh"  : "smallBtagHigh",
            "smallBtagLow"   : "smallBtagLow",

            "BtagHigh":"btagHigh",
            "BtagLow":"btagLow",
            "BtagCHigh":"ctagHigh",
            "BtagCLow":"ctagLow",
            "BtagLHigh":"mistagHigh",
            "BtagLLow": "mistagLow",
            "Btag0High":"btag0High",
            "Btag0Low":"btag0Low",
            "Btag1High":"btag1High",
            "Btag1Low":"btag1Low",
            "Btag2High":"btag2High",
            "Btag2Low":"btag2Low",
            "Btag3High":"btag3High",
            "Btag3Low":"btag3Low",
            "Btag4High":"btag4High",
            "Btag4Low":"btag4Low",
            "Btag5High":"btag5High",
            "Btag5Low":"btag5Low",
            "Btag6High":"btag6High",
            "Btag6Low":"btag6Low",
            "Btag7High":"btag7High",
            "Btag7Low":"btag7Low",
            "Btag8High":"btag8High",
            "Btag8Low":"btag8Low",
            "Btag9High":"btag9High",
            "Btag9Low":"btag9Low",
            "Btag10High":"btag10High",
            "Btag10Low":"btag10Low",
            "BtagC0Low":"ctag0Low",
            "BtagC0High":"ctag0High",
            "BtagC1High":"ctag1High",
            "BtagC1Low":"ctag1Low",
            "BtagC2High":"ctag2High",
            "BtagC2Low":"ctag2Low",
            "BtagC3High":"ctag3High",
            "BtagC3Low":"ctag3Low",
            "BtagC4High":"ctag4High",
            "BtagC4Low":"ctag4Low",
            "BtagC5High":"ctag5High",
            "BtagC5Low":"ctag5Low",
            "BtagC6High":"ctag6High",
            "BtagC6Low":"ctag6Low",
            "BtagL0High":"mistag0High",
            "BtagL0Low":"mistag0Low",
            "BtagL1High":"mistag1High",
            "BtagL1Low":"mistag1Low",
            "BtagL2High":"mistag2High",
            "BtagL2Low":"mistag2Low",
            "BtagL3High":"mistag3High",
            "BtagL3Low":"mistag3Low",
            "BtagL4High":"mistag4High",
            "BtagL4Low":"mistag4Low",
            "BtagL5High":"mistag5High",
            "BtagL5Low":"mistag5Low",
            "BtagL6High":"mistag6High",
            "BtagL6Low":"mistag6Low",
            "BtagL7High":"mistag7High",
            "BtagL7Low":"mistag7Low",
            "BtagL8High":"mistag8High",
            "BtagL8Low":"mistag8Low",
            "BtagL9High":"mistag9High",
            "BtagL9Low":"mistag9Low",
            "BtagL10High":"mistag10High",
            "BtagL10Low":"mistag10Low",
            "BtagL11High":"mistag11High",
            "BtagL11Low":"mistag11Low",
            "BtagL12High":"mistag12High",
            "BtagL12Low":"mistag12Low",
            "MCGenLow":"MC_GenLow",
            "MCGenHigh":"MC_GenHigh",
            "PartonShowerHigh":"PSHigh",
            "PartonShowerLow":"PSLow",
            "EWSHigh":"EWSHigh",
            "EWSLow":"EWSLow",
            "PtttHigh":"PtttHigh",
            "PtttLow":"PtttLow",
            "hdampHigh":"hdampHigh",
            "hdampLow":"hdampLow",
            "IFSRLow":"IFSRHigh",
            "IFSRHigh":"IFSRLow",
            "topmassHigh":"top_massHigh",
            "topmassLow":"top_massLow",
            "JetEnerResHigh" :  "JERHigh",
            "JetEnerResLow"  :  "JERLow",

            #"JESHigh"  : "JESHigh",
            #"JESLow"   : "JESLow",
            "JES_ALLHigh"  : "JESALLHigh",
            "JES_ALLLow"   : "JESALLLow",
            "smallJESHigh"  : "smallJESHigh",
            "smallJESLow"   : "smallJESLow",
            "JESsmallHigh"  : "JESsmallHigh",
            "JESsmallLow"   : "JESsmallLow",
            "JES0High" : "JES0High",
            "JES0Low"  : "JES0Low",
            "JES1High" : "JES1High",
            "JES1Low"  : "JES1Low",
            "JES2High" : "JES2High",
            "JES2Low"  : "JES2Low",
            "JES3High" : "JES3High",
            "JES3Low"  : "JES3Low",
            "JES4High" : "JES4High",
            "JES4Low"  : "JES4Low",
            "JES5High" : "JES5High",
            "JES5Low"  : "JES5Low",
            "JES6High" : "JES6High",
            "JES6Low"  : "JES6Low",
            "JES7High" : "JES7High",
            "JES7Low"  : "JES7Low",
            "JES8High" : "JES8High",
            "JES8Low"  : "JES8Low",
            "JES9High" : "JES9High",
            "JES9Low"  : "JES9Low",
            "JES10High" : "JES10High",
            "JES10Low"  : "JES10Low",
            "JES11High" : "JES11High",
            "JES11Low"  : "JES11Low",
            "JES12High" : "JES12High",
            "JES12Low"  : "JES12Low",
            "JES13High" : "JES13High",
            "JES13Low"  : "JES13Low",
            "JES14High" : "JES14High",
            "JES14Low"  : "JES14Low",
            "JES15High" : "JES15High",
            "JES15Low"  : "JES15Low",
            "JES16High" : "JES16High",
            "JES16Low"  : "JES16Low",
            "JES17High" : "JES17High",
            "JES17Low"  : "JES17Low",
            "JES18High" : "JES18High",
            "JES18Low"  : "JES18Low",
            "JES19High" : "JES19High",
            "JES19Low"  : "JES19Low",
            "JES20High" : "JES20High",
            "JES20Low"  : "JES20Low",
            "JES21High" : "JES21High",
            "JES21Low"  : "JES21Low",
            "JES22High" : "JES22High",
            "JES22Low"  : "JES22Low",
            "BoostedJESHigh" :  "BoostedJESHigh",
            "BoostedJESLow"  :  "BoostedJESLow",
            "BoostedJES0High" : "BoostedJES0High",
            "BoostedJES0Low"  : "BoostedJES0Low",
            "BoostedJES1High" : "BoostedJES1High",
            "BoostedJES1Low"  : "BoostedJES1Low",
            "BoostedJES2High" : "BoostedJES2High",
            "BoostedJES2Low"  : "BoostedJES2Low",
            "BoostedJES3High" : "BoostedJES3High",
            "BoostedJES3Low"  : "BoostedJES3Low",
            "BoostedJES4High" : "BoostedJES4High",
            "BoostedJES4Low"  : "BoostedJES4Low",
            "BoostedJES5High" : "BoostedJES5High",
            "BoostedJES5Low"  : "BoostedJES5Low",
            "BoostedJES6High" : "BoostedJES6High",
            "BoostedJES6Low"  : "BoostedJES6Low",
            "BoostedJES7High" : "BoostedJES7High",
            "BoostedJES7Low"  : "BoostedJES7Low",
            "BoostedJES8High" : "BoostedJES8High",
            "BoostedJES8Low"  : "BoostedJES8Low",
            "BoostedJES9High" : "BoostedJES9High",
            "BoostedJES9Low"  : "BoostedJES9Low",
            "BoostedJES10High" : "BoostedJES10High",
            "BoostedJES10Low"  : "BoostedJES10Low",
            "BoostedJES11High" : "BoostedJES11High",
            "BoostedJES11Low"  : "BoostedJES11Low",
            "BoostedJES12High" : "BoostedJES12High",
            "BoostedJES12Low"  : "BoostedJES12Low",
            "BoostedJES13High" : "BoostedJES13High",
            "BoostedJES13Low"  : "BoostedJES13Low",
            "BoostedJES14High" : "BoostedJES14High",
            "BoostedJES14Low"  : "BoostedJES14Low",
            "BoostedJES15High" : "BoostedJES15High",
            "BoostedJES15Low"  : "BoostedJES15Low",
            "BoostedJES16High" : "BoostedJES16High",
            "BoostedJES16Low"  : "BoostedJES16Low",
            "BoostedJESothersHigh" : "BoostedJESothersHigh",
            "BoostedJESothersLow"  : "BoostedJESothersLow",
            "BoostedJMSHigh" : "BoostedJMSHigh",
            "BoostedJMSLow"  : "BoostedJMSLow",
            "BoostedJERHigh" : "BoostedJERHigh",
            "BoostedJERLow"  : "BoostedJERLow",
            "BoostedJMRHigh" : "BoostedJMRHigh",
            "BoostedJMRLow"  : "BoostedJMRLow",

            "BoostedJES_ALLHigh" : "BoostedJES_ALLHigh",
            "BoostedJES_ALLLow"  : "BoostedJES_ALLLow",
            

            "PDFLow":"PDFLow",
            "PDFHigh":"PDFHigh",


            #"muonScaleLow" : "MUSCLow",
            #"muonScaleHigh" : "MUSCHigh",
            "WHFC0High":"WHFC0High",
            "WHFC3High":"WHFC3High",
            "WHFC0Low":"WHFC0Low",
            "WHFC3Low":"WHFC3Low",
            "WHFC4High":"WHFC4High",
            "NormWHigh":"NormWHigh",
            "WHFC4Low":"WHFC4Low",
            "NormWLow":"NormWLow",




            "JVFCutHigh": "JVFHigh",
            "JVFCutLow" : "JVFLow",
            #"electronEnergyRescaleHigh" : "EERHigh",
            #"electronEnergyRescaleLow"  :  "EERLow",
            #"electronSmearHigh" : "EESHigh",
            #"electronSmearLow"  : "EESLow",
            #"BoostedJES0High" : "FATJPTSHigh",
            #"BoostedJES0Low"  : "FATJPTSLow",
            #"BoostedJES1High" : "FATJMSHigh",
            #"BoostedJES1Low"  : "FATJMSLow",
            #"BoostedJES2High" : "FATD12High",
            #"BoostedJES2Low"  : "FATD12Low",

            #"JEELow"  :  "JEffLow",
            #"JEEHigh"  : "JEffHigh",

            "JER0High" :    "JER0High",
            "JER0Low"  :    "JER0Low",
            "JER1High" :    "JER1High",
            "JER1Low"  :    "JER1Low",
            "JER2High" :    "JER2High",
            "JER2Low"  :    "JER2Low",
            #"Btag_r1High": "Btag_r1High",
            #"Btag_r2High": "Btag_r2High",
            #"Btag_r3High": "Btag_r3High",
            #"Btag_b1High": "Btag_b1High",
            #"Btag_b2High": "Btag_b2High",
            #"Btag_b3High": "Btag_b3High",
            #"BtagC_r1High": "BtagC_r1High",
            #"BtagC_r2High": "BtagC_r2High",
            #"BtagC_r3High": "BtagC_r3High",
            #"BtagC_b1High": "BtagC_b1High",
            #"BtagC_b2High": "BtagC_b2High",
            #"BtagC_b3High": "BtagC_b3High",
            #"BtagL_r1High": "BtagL_r1High",
            #"BtagL_r2High": "BtagL_r2High",
            #"BtagL_r3High": "BtagL_r3High",
            #"BtagL_b1High": "BtagL_b1High",
            #"BtagL_b2High": "BtagL_b2High",
            #"BtagL_b3High": "BtagL_b3High",
            #"Btag_r1Low": "Btag_r1Low",
            #"Btag_r2Low": "Btag_r2Low",
            #"Btag_r3Low": "Btag_r3Low",
            #"Btag_b1Low": "Btag_b1Low",
            #"Btag_b2Low": "Btag_b2Low",
            #"Btag_b3Low": "Btag_b3Low",
            #"BtagC_r1Low": "BtagC_r1Low",
            #"BtagC_r2Low": "BtagC_r2Low",
            #"BtagC_r3Low": "BtagC_r3Low",
            #"BtagC_b1Low": "BtagC_b1Low",
            #"BtagC_b2Low": "BtagC_b2Low",
            #"BtagC_b3Low": "BtagC_b3Low",
            #"BtagL_r1Low": "BtagL_r1Low",
            #"BtagL_r2Low": "BtagL_r2Low",
            #"BtagL_r3Low": "BtagL_r3Low",
            #"BtagL_b1Low": "BtagL_b1Low",
            #"BtagL_b2Low": "BtagL_b2Low",
            #"BtagL_b3Low": "BtagL_b3Low",
            #"METScaleSoftJetHigh": "MET_SJSHigh",
            #"METScaleSoftJetLow" : "MET_SJSLow",
            #"METResoSoftJetHigh" : "MET_RSJHigh",
            #"METResoSoftJetLow"  : "MET_RSJLow",
            #"muonSmearSystMSLow"   :  "MUMSDOWN",
            #"muonSmearSystMSHigh"  :  "MUMSUP",
            #"muonSmearSystIDLow"   :  "MUIDDOWN",
            #"muonSmearSystIDHigh"  :  "MUIDUP",
            #"muonSmearHigh" : "muonSmearHigh",
            #"muonSmearLow"  : "muonSmearLow",
            "norm_ttHigh" : "Norm_ttbarHigh",
            "norm_ttLow"  : "Norm_ttbarLow",
            #"dibosonLow"    : "Norm_dibosonLow",
            #"dibosonHigh"   : "Norm_dibosonHigh",
            #"singletopLow"  : "Norm_stopLow",
            #"singletopHigh" : "Norm_stopHigh",
            #"ZjetsMCNormHigh" : "Norm_zjetsHigh",
            #"ZjetsMCNormLow"  : "Norm_zjetsLow",
            "iqopt3Low"		: "Iqopt3Low",
            "iqopt3High"	: "Iqopt3High",
            "ptjmin10Low"	:"Ptjmin10Low",
            "ptjmin10High"	:"Ptjmin10High",
            "WhfsfLow"		:"WhfsfLow",
            "WhfsfHigh"		:"WhfsfHigh",
            "luminosityHigh"	:"luminosityHigh",
            "luminosityLow"	:"luminosityLow",
            "norm_QCDeHigh":"Norm_QCDeHigh",
            "norm_QCDeLow" :"Norm_QCDeLow",
            "norm_QCDmuHigh":"Norm_QCDmuHigh",
            "norm_QCDmuLow" :"Norm_QCDmuLow",
            "norm_QCDHigh":"Norm_QCDHigh",
            "norm_QCDLow" :"Norm_QCDLow",
        }

        def CleanTheList(name):
            if name.find("DRMin")!= -1 :
                return True
            #if name.find("QCD")!= -1 :
            #    return True
            if name.find(mass_ttbar)==-1 :
                return True
            if lepton_channel=="Electrons":
                if name.find("_e") == -1 :
                    return True
                if name.find("_mu_") != -1 :
                    return True
                if name.find("_lep_") != -1 :
                    return True
            elif lepton_channel=="Muons" :
                if name.find("_mu") == -1 :
                    return True
                if name.find("_e_") != -1 :
                    return True
                if name.find("_lep_") != -1 :
                    return True
            elif lepton_channel=="Leptons" :
                if name.find("_lep") == -1 :
                    return True
                if name.find("_e_") != -1 :
                    return True
                if name.find("_mu_") != -1 :
                    return True
                
            #if name.find("ttbar_") == -1 and name.find("_ttbarnorm") != -1:
                #return True
            #if name.find("VV_") == -1 and name.find("_diboson") != -1:
                #return True
            #if name.find("Zjets_") == -1 and name.find("_ZjetsMCNorm") != -1:
                #return True
            #if name.find("Wjets_") == -1 and name.find("_Wjets") != -1:
                #return True
            #if name.find("stop_") == -1 and name.find("_singletop") != -1:
                #return True
            return False

        def ReNameHisto(name):
            btagC=""
            if name.find("_cat1") != -1:
                name=name.replace("_cat1","")
                btagC="SRexcbtag_cat1"
            if name.find("_cat2") != -1:
                name=name.replace("_cat2","")
                btagC="SRexcbtag_cat2"
            if name.find("_cat3") != -1:
                name=name.replace("_cat3","")
                btagC="SRexcbtag_cat3"
            if name.find("binkkg_") != -1:
                name=name.replace("binkkg_","")
            if lepton_channel=="Electrons":
                if name.find(mass_ttbar+"_e") != -1 :
                    name=name.replace("_"+mass_ttbar+"_e","Nom_"+btagC+lepton_channel+BoostedOrResolved)
            elif lepton_channel=="Muons":
                if name.find(mass_ttbar+"_mu") != -1 :
                    name=name.replace("_"+mass_ttbar+"_mu","Nom_"+btagC+lepton_channel+BoostedOrResolved)
            elif lepton_channel=="Leptons":
                if name.find(mass_ttbar+"_lep") != -1 :
                    name=name.replace("_"+mass_ttbar+"_lep","Nom_"+btagC+lepton_channel+BoostedOrResolved)
            if name.find("_up") != -1:
                name=name.replace("_up","High")
            elif name.find("_dw") != -1:
                name=name.replace("_dw","Low")
            name="h"+name+"_obs_cuts"
            #name="h"+name+"_obs_cutsNorm"
            if name.find("_up") != -1:
                name=name.replace("_up","High")
            if name.find("_smeared") !=-1:
                name=name.replace("_smeared","")
                name=name.replace("Nom","smearedNom")
            for syst in Syst_list:
                if name.find("_"+syst) !=-1 :
                    name=name.replace(syst+"_","")
                    name=name.replace("Nom",Syst_list[syst])
                    #name=name.replace("obs_cuts", "obs_cutsNorm")
                    
                    break
            if name.find("hQCDe") != -1:
                name=name.replace("hQCDe","hQCD")
            elif name.find("hQCDmu") != -1:
                name=name.replace("hQCDmu","hQCD")
            return name



        print "open file"
        myfile = TFile.Open(filename,"READ")
        next = myfile.GetListOfKeys().MakeIterator()
        n_syst_hist=0
        key = next()
        hlist = TObjArray(0)
        n_histos=0

        print "loop over keys"
        while (key and n_histos<1e10) :
            obj = myfile.Get(key.GetName())
            if (obj.InheritsFrom("TH1")) :
                histo_name=obj.GetName()
                #print "key ",histo_name
                if histo_name.find("smeared")!=-1:
                    key = next()
                    continue
                #print "t1"
                #if histo_name.find("binkkg")==-1:
                #    key = next()
                #    continue
                #print "t2"
                if histo_name.find("QCDe")!=-1 and ( histo_name.find("_mu_")!=-1 or histo_name.endswith("_mu")) :
                    key = next()
                    continue
                #print "t3"
                if histo_name.find("QCDmu")!=-1 and ( histo_name.find("_e_")!=-1 or histo_name.endswith("_e")) :
                    key = next()
                    continue
                #print "t4"
                if CleanTheList(obj.GetName())==True:
                    key = next()
                    continue
                #print "t5"
                #if histo_name.find("Diboson")==-1:
                    #key = next()
                    #continue
                newHisto_name=obj.GetName()
                newHisto_name=newHisto_name.replace("binkkg_","")
                #print "1 ", newHisto_name
                newHisto_name=ReNameHisto(newHisto_name)
                #print "2 ", newHisto_name
                obj.SetName(newHisto_name)
                obj=SetRangeofXaxis(obj)
                hlist.Add(obj)
                n_histos+=1
            key = next()

        print "loop over histograms"
        for h in hlist:
            hname=h.GetName()
            if hname.find("hQCDe")!= -1 :
                hname=hname.replace("hQCDe","hQCD")
                h.SetName(hname)
            if hname.find("hQCDmu")!= -1 :
                hname=hname.replace("hQCDmu","hQCD")
                h.SetName(hname)
            for sample in ["ttbar_POWHEG","VV","stop","Zjets","Wjets","others"] :
                if h.GetName().find(sample+"NormHigh")!= -1 or h.GetName().find(sample+"NormLow")!= -1 :
                    name=h.GetName()
                    h.SetName(name.replace(sample+"Norm",sample+"Norm_"+sample+"_"))
                    h.SetTitle(name.replace(sample+"Norm",sample+"Norm_"+sample+"_"))
            # no longer want to have 2 copies of every histo, but no fix for now...
            if not "_obs_cutsNorm" in h.GetName():
                if h.GetName().find("_obs_cuts")== -1 :
                    continue
                CopyHisto(h)


        #if lepton_channel=="Electrons":
            #histName ="binkkg_QCDe_"+mass_ttbar+"_"+lep_channel
            #histNom = myfile.Get(histName)
            #histNom.SetName(ReNameHisto(histName))
            #hlist.Add(histNom)
            #histName ="binkkg_QCDe_"+mass_ttbar+"_"+lep_channel+"_norm_QCDe_up"
            #histVarUp = myfile.Get(histName)
            #histVarUp.SetName(ReNameHisto(histName))
            #hlist.Add(histVarUp)
            #histName ="binkkg_QCDe_"+mass_ttbar+"_"+lep_channel+"_norm_QCDe_dw"
            #histVarDown = myfile.Get(histName)
            #histVarDown.SetName(ReNameHisto(histName))
            #hlist.Add(histVarDown)
            #hist_QCDmuNormHigh_e=histNom.Clone()
            #hist_QCDmuNormLow_e=histNom.Clone()
            #hist_QCDmuNormHigh_e.SetName("hQCDNorm_QCDmuHigh_SRincElectronsBoosted_obs_cuts")
            #hist_QCDmuNormLow_e.SetName("hQCDNorm_QCDmuLow_SRincElectronsBoosted_obs_cuts")
            #hlist.Add(hist_QCDmuNormHigh_e)
            #hlist.Add(hist_QCDmuNormLow_e)

        #if lepton_channel=="Muons":
            #histName ="binkkg_QCDmu_"+mass_ttbar+"_"+lep_channel
            #histNom = myfile.Get(histName)
            #histNom.SetName(ReNameHisto(histName))
            #hlist.Add(histNom)
            #histName ="binkkg_QCDmu_"+mass_ttbar+"_"+lep_channel+"_norm_QCDmu_up"
            #histVarUp = myfile.Get(histName)
            #histVarUp.SetName(ReNameHisto(histName))
            #hlist.Add(histVarUp)
            #histName ="binkkg_QCDmu_"+mass_ttbar+"_"+lep_channel+"_norm_QCDmu_dw"
            #histVarDown = myfile.Get(histName)
            #histVarDown.SetName(ReNameHisto(histName))
            #hlist.Add(histVarDown)
            #hist_QCDeNormHigh_mu=histNom.Clone()
            #hist_QCDeNormLow_mu=histNom.Clone()
            #hist_QCDeNormHigh_mu.SetName("hQCDNorm_QCDeHigh_SRincMuonsBoosted_obs_cuts")
            #hist_QCDeNormLow_mu.SetName("hQCDNorm_QCDeLow_SRincMuonsBoosted_obs_cuts")
            #hlist.Add(hist_QCDeNormHigh_mu)
            #hlist.Add(hist_QCDeNormLow_mu)
    
        print "Opening outfile"
        outf = TFile.Open("MassSpectra_Paper2014_"+BoostedOrResolved+"_"+lep_channel+".root","RECREATE")

        print "writing"
        outf.cd()
        hlist.Write()
        print "closing outfile"
        outf.Close()

        print "closing infile"
        myfile.Close()
        print "next..."

print "finished!"
