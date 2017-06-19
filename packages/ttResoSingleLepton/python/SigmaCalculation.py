#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python

# script to calculate LO cross section for Z'->TC 
# J. Ferrando, University of Glasgow 29/07/2013

import TopColor as tc2
from ROOT import TFile, gROOT, TCanvas, gStyle, TLatex, TLegend, Double

execfile("./root_set.py")

gROOT.SetStyle("Plain")
gROOT.SetStyle("zeus_pub")
gROOT.ForceStyle()
gStyle.SetOptTitle(0)

tcpar=tc2.TopColorParams()
widthcalc=tc2.ZpWidthCalculator(tcpar)
brcalc=tc2.ZpBrCalculator(tcpar)
sigcalc=tc2.ZpSigmaCalculator(widthcalc)
sigcalc.initPDF()

#IMPORTANT INPUT PARAMETERS
# follow arXiv 1112:4928 and only consider two flavours
sigcalc.flavours=[1,2]
# increase this number for more precision
nsteps=40000
# width (in percent)
width=1.2
# centre-of-mass-energy
ecom=8000.0
#ecom=13000.0
#pdf choice
pdf="cteq6ll.LHpdf"

print "LEADING ORDER RESULTS FOR TOPCOLOR MODEL IV SCENARIO"


print "mass (GeV)", "width (%)", "cos^2(theta_H)", "cos(theta_H)", "BR(Z'->ttbar)","sigma (pb)"
for i in range(8,81):
    mzp=i*50.0
    topwidth=widthcalc.Width(mzp)
    ratio= width/(topwidth)
    newc2t= (tcpar.cot2thetah)*ratio
    tcpar.cot2thetah=newc2t
    topwidth=widthcalc.Width(mzp)
    sigcalc.nsteps=nsteps
    sigcalc.sqrts=ecom
    sigcalc.pdf=pdf
    res=sigcalc.Sigma(mzp)
    # Gev/C to picobarns
    hbar=1.054571726e-34
    c=299792458
    gevj=1.602176565e-10
    #        sifac=1.0e6*1.0e34*((hbar*c)**2)/(gevj**2)
    #       1e6 already applied
    sifac=1.0e34*((hbar*c)**2)/(gevj**2)
#    print mzp, topwidth, newc2t, newc2t**0.5, brcalc.Br(6,mzp),"sigma:" , sifac*res[0], "+-", sifac*res[1], "pb"
    print mzp, topwidth, newc2t, newc2t**0.5, brcalc.Br(6,mzp),"sigma:" , sifac*res, "pb"



