#!/usr/bin/env python
from ROOT import gROOT,gSystem,gDirectory
gSystem.Load("libSusyFitter.so")
from ROOT import ConfigMgr,FitConfig #this module comes from gSystem.Load("libSusyFitter.so")
gROOT.Reset()

from ROOT import TFile, RooWorkspace, TObject, TString, RooAbsReal, RooRealVar, RooFitResult, RooDataSet, RooAddition, RooArgSet,RooAbsData,RooRandom 
from ROOT import Util, TMath
from ROOT import RooFit
from ROOT import RooExpandedFitResult
    
import os
import sys
from sys import exit

from PrintFitResultTex import *
from PrintFitResultPlot import *
import pickle


def getnamemap():

  namemap = {}
  namemap['alpha_BoostedJES0'] 		= 'Akt10 JES0 (Gamma-jet)                 '
  namemap['alpha_BoostedJES1'] 		= 'Akt10 JES1 (dataMC)                    '
  namemap['alpha_BoostedJES2'] 		= 'Akt10 JES2 (pt2Cut)                    '
  namemap['alpha_BoostedJES3'] 		= 'Akt10 JES3 (dPhiCut)                   '
  namemap['alpha_BoostedJES4'] 		= 'Akt10 JES4 (photonPurity)              '
  namemap['alpha_BoostedJES5'] 		= 'Akt10 JES5 (PES)                       '
  namemap['alpha_BoostedJES6'] 		= 'Akt10 JES6 (photonPurity)              '
  namemap['alpha_BoostedJES7'] 		= 'Akt10 JES7 (kterm)                     '
  namemap['alpha_BoostedJES8'] 		= 'Akt10 JES8 (JER)                       '
  namemap['alpha_BoostedJES9'] 		= 'Akt10 JES9 (akt4insideOutsideLargeR)   '
  namemap['alpha_BoostedJES10'] 	= 'Akt10 JES10 (more1smallJetInsideLargeR)'
  namemap['alpha_BoostedJES11'] 	= 'Akt10 JES11 (stats)                    '
  namemap['alpha_BoostedJES12'] 	= 'Akt10 JES12 (MoverPt)                  '
  namemap['alpha_BoostedJES13'] 	= 'Akt10 JES13 (Topology)                 '
  namemap['alpha_BoostedJES14'] 	= 'Akt10 JES14 (Extrapl.)                 '
  namemap['alpha_BoostedJES15'] 	= 'Akt10 JES15 (NPV)                      '
  namemap['alpha_BoostedJES16'] 	= 'Akt10 JES16 ($\mu$)                    '
  namemap['alpha_BoostedJMS'] 		= 'Akt10 jet m/d12 scale                  '
  namemap['alpha_BoostedJER'] 		= 'Akt10 jet energy res.                  '
  namemap['alpha_BoostedJMR'] 		= 'Akt10 jet mass res.                    '
  namemap['alpha_BoostedJESothers'] 	= 'Akt10 JES (Others)                     '
  namemap['alpha_JES3'] 		= 'Akt4 JES3 (Modelling1)                 '
  namemap['alpha_JES7'] 		= 'Akt4 JES7 (Detector1)                  '
  namemap['alpha_JES12'] 		= 'Akt4 JES12 (EtaInter)                  '
  namemap['alpha_JES18'] 		= 'Akt4 JES18                             '
  namemap['alpha_JES20'] 		= 'Akt4 JES20 (PU Rho)                    '
  namemap['alpha_JES21'] 		= 'Akt4 JES21                             '
  namemap['alpha_JES22'] 		= 'Akt4 JES22                             '
  namemap['JESsmall'] 			= 'Akt4 JES (Others)                      '
  namemap['alpha_JER'] 			= 'Akt4 jet energy res.                   '
  namemap['alpha_JVF'] 			= 'Jet vertex fraction                    '
  namemap['alpha_EWS'] 			= '\\ttbar\\ EW Sudakov                   '
  namemap['alpha_ElectronSF'] 		= 'Electron scale factor                  '
  namemap['alpha_MuonSF'] 		= 'Muon scale factor                      '
  namemap['alpha_IFSR'] 		= '\\ttbar{} ISR, FSR                     '
  namemap['alpha_Iqopt3'] 		= '$W$ shape, ``iqopt3''                  '
  namemap['alpha_Ptjmin10'] 		= '$W$ shape, ``ptjmin10''                '
  namemap['alpha_MC_Gen'] 		= '\\ttbar{}gGenerator dependence         '
  namemap['alpha_Norm_ttbar'] 		= '\\ttbar{} normalisation                '
  namemap['alpha_luminosity'] 		= 'Luminosity                             '
  namemap['alpha_PS'] 			= '\\ttbar\\ Parton shower                '
  namemap['alpha_top_mass'] 		= 'Top mass                               '
  namemap['alpha_PDF'] 			= 'PDF uncertainty                        '
  namemap['alpha_Whfsf'] 		= '$W$ normalisation                      '
  namemap['alpha_btag6'] 		= '$b$-tag EV6                            '
  namemap['alpha_btag7'] 		= '$b$-tag EV7                            '
  namemap['alpha_btag8'] 		= '$b$-tag EV8                            '
  namemap['alpha_btag9'] 		= '$b$-tag EV9                            '
  namemap['alpha_btag10'] 		= '$b$-tag high \pt                       '
  namemap['alpha_ctag'] 		= '$c$-tag                                '
  namemap['alpha_mistag'] 		= 'Mistag                                 '
  namemap['alpha_hdamp']		= 'hdamp parameter                        '
  namemap['alpha_Norm_QCDe'] 		= 'Multi-jets norm, $e$+jets              '
  namemap['alpha_Norm_QCDmu'] 		= 'Multi-jets norm, $\\mu$+jets           '

  return namemap

  

def latexfitresults( filename, resultName="RooExpandedFitResult_afterFit", outName="test.tex" ):

  namemap = {}
  namemap = getnamemap()

  ############################################
  workspacename = 'w'
  w = Util.GetWorkspaceFromFile(filename,workspacename)

  if w==None:
    print "ERROR : Cannot open workspace : ", workspacename
    sys.exit(1) 

  result = w.obj(resultName)
  if result==None:
    print "ERROR : Cannot open fit result ", resultName
    sys.exit(1)

  #####################################################

  regSys = {}

  # calculate error per parameter on  fitresult
  fpf = result.floatParsFinal() 
  fpi = result.floatParsInit()

  '''
  // P r i n t   l a t ex   t a b l e   o f   p a r a m e t e r s   o f   p d f 
  // --------------------------------------------------------------------------


  // Print parameter list in LaTeX for (one column with names, one column with values)
  params->printLatex() ;

  // Print parameter list in LaTeX for (names values|names values)
  params->printLatex(Columns(2)) ;

  // Print two parameter lists side by side (name values initvalues)
  params->printLatex(Sibling(*initParams)) ;

  // Print two parameter lists side by side (name values initvalues|name values initvalues)
  params->printLatex(Sibling(*initParams),Columns(2)) ;

  // Write LaTex table to file
  params->printLatex(Sibling(*initParams),OutputFile("rf407_latextables.tex")) ;
  '''

  ####fpf.printLatex(RooFit.Format("NE",RooFit.AutoPrecision(2),RooFit.VerbatimName()),RooFit.Sibling(fpi),RooFit.OutputFile(outName)) 

  # set all floating parameters constant
  for idx in range(fpf.getSize()):
    parname = fpf[idx].GetName()
    if parname.find("gamma_stat_") > -1 :
        continue
    ip = fpi[idx]
    ipv  = ip.getVal()
    ipe  = ip.getError()
    ipel = ip.getErrorLo()
    ipeh = ip.getErrorHi()

    fp = fpf[idx]
    fpv  = fp.getVal()
    fpe  = fp.getError()
    fpel = fp.getErrorLo()
    fpeh = fp.getErrorHi()

    name = parname
    if namemap.has_key(name): name = namemap[name]

    regSys[name] = (ipv,ipe,ipel,ipeh,fpv,fpe,fpel,fpeh)

  return regSys




##################################

# MAIN

if __name__ == "__main__":
  
  import os, sys
  import getopt
  def usage():
    print "Usage:"
    print "PrintFitResult.py [-c channel] [-w workspace_afterFit] [-o outputFileName]\n"
    print "Minimal set of inputs [-c channels] [-w workspace_afterFit]"
    print "*** Options are: "
    print "-c <analysis name>: single name accepted only (OBLIGATORY) "
    print "-w <workspaceFileName>: single name accepted only (OBLIGATORY) ;   if multiple channels/regions given in -c, assumes the workspace file contains all channels/regions"
    sys.exit(0)        

  wsFileName='/results/MyOneLeptonKtScaleFit_HardLepR17_BkgOnlyKt_combined_NormalMeasurement_model_afterFit.root'
  try:
    opts, args = getopt.getopt(sys.argv[1:], "o:c:w:m:f:s:%b")
  except:
    usage()
  if len(opts)<1:
    usage()

  analysisName = ''
  outputFileName="default"
  method="1"
  showAfterFitError=True
  showPercent=False
  for opt,arg in opts:
    if opt == '-c':
      analysisName=arg
    if opt == '-w':
      wsFileName=arg

  resultName = 'RooExpandedFitResult_afterFit'
  if not showAfterFitError:
    resultName =  'RooExpandedFitResult_beforeFit'

  regSys = latexfitresults(wsFileName,resultName,outputFileName)

  line_chanSysTight = tablefragment(regSys,analysisName)
  plot = plotfragment(regSys,analysisName)
  #plot.Print("fitresult_" + analysisName + ".pdf")

  outputFileName = "fitresult_" + analysisName + ".tex"
  printname = analysisName
  #printname = printname.replace('_','\_')
  printname = printname.replace('_',' ')

  f = open(outputFileName, 'w')
  f.write( line_chanSysTight )
  f.write( "\n\\begin{figure*}[htb!]\n\\includegraphics[width=\\textwidth]{fitresult_" + analysisName + ".pdf}\n\\caption[]{Fit parameters for " + printname + ".}\n\\end{figure*}\n" )
  f.close()
  print "\nwrote results in file: %s"%(outputFileName)

