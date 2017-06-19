#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python

import TopColor as tc2
import XsecTools as xst
from ROOT import TFile, gROOT, TCanvas, gStyle, TLatex, TLegend, Double

execfile("./root_set.py")

gROOT.SetStyle("Plain")
gROOT.SetStyle("zeus_pub")
gROOT.ForceStyle()
gStyle.SetOptTitle(0)


tcpar=tc2.TopColorParams()
brcalc=tc2.ZpBrCalculator(tcpar)


print "LEADING ORDER BRANCHING RATIOS FOR TOPCOLOR MODEL IV SCENARIO"

print "Mass\tBR(tt)\tBR(uu)\tBR(bb)\tBR(dd)\tSum"

for i in range(10,80):
    mzp=i*50.0
    topbr=brcalc.Br(6,mzp)
    upbr=brcalc.Br(2,mzp)
    bottombr=brcalc.Br(5,mzp)
    downbr=brcalc.Br(1,mzp)
    sum =topbr+upbr+bottombr+downbr
    print str(mzp)+"\t"+str(topbr)[0:9]+"\t"+str(upbr)[0:9]+"\t"+str(bottombr)[0:9]+"\t"+str(downbr)[0:9]+"\t"+str(sum)


ups=[]
downs=[]
bottoms=[]
tops=[]

for mzp in [350.0,360.0,370.0,380.0,390.0,400.0,420.0,440.0,460.0,480.0,500.0,550.0,600.0,650.0,700.0,750.0,800.0,900.0,1000.0,1200.0,1400.0,1800.0,2000.0]:
    m_up=xst.CrossSectionPDFErrors()
    m_up.description=str(mzp/1000.0)
    m_up.SetValueAndStatError(brcalc.Br(2,mzp),0.0)
    m_down=xst.CrossSectionPDFErrors()
    m_down.description=str(mzp/1000.0)
    m_down.SetValueAndStatError(brcalc.Br(1,mzp),0.0)
    m_b=xst.CrossSectionPDFErrors()
    m_b.description=str(mzp/1000.0)
    m_b.SetValueAndStatError(brcalc.Br(5,mzp),0.0)
    m_top=xst.CrossSectionPDFErrors()
    m_top.description=str(mzp/1000.0)
    m_top.SetValueAndStatError(brcalc.Br(6,mzp),0.0)
    ups.append(m_up)
    downs.append(m_down)
    bottoms.append(m_b)
    tops.append(m_top)

#    topbr=brcalc.Br(6,mzp)
#    upbr=brcalc.Br(2,mzp)
#    bottombr=brcalc.Br(5,mzp)
#    downbr=brcalc.Br(1,mzp)
#    sum =topbr+upbr+bottombr+downbr
#    print str(mzp)+"\t"+str(topbr)[0:6]+"\t"+str(upbr)[0:6]+"\t"+str(bottombr)[0:6]+"\t"+str(downbr)[0:6]+"\t"+str(sum)



graph_up=xst.MakeGraphs(ups,"uBR")[0]
graph_down=xst.MakeGraphs(downs,"dBR")[0]
graph_top=xst.MakeGraphs(tops,"bBR")[0]
graph_bottom=xst.MakeGraphs(bottoms,"tBR")[0]


c1=c1=TCanvas("SelEff","SelEff",800,600)


ghist=graph_up.GetHistogram()
ghist.SetAxisRange(0.355,2,"X")
ghist.SetMaximum(0.5)
ghist.SetTitle("")
ghist.Draw("X")

graph_up.SetLineWidth(3)
graph_up.SetLineColor(11)
graph_up.Draw("l")

graph_down.SetLineWidth(3)
graph_down.SetLineStyle(2)
graph_down.SetLineColor(2)
graph_down.Draw("l")


graph_bottom.SetLineWidth(3)
graph_bottom.SetLineStyle(3)
graph_bottom.SetLineColor(4)
graph_bottom.Draw("l")

graph_top.SetLineWidth(3)
graph_top.SetLineStyle(4)
graph_top.SetLineColor(1)
graph_top.Draw("l")

xTit = TLatex()
xTit.SetTextSize(0.06)
xTit.DrawLatex(1.4,ghist.GetMinimum()-(1.0-ghist.GetMinimum())*0.1,"m(Z'_{TC2}) [TeV]")

yTit = TLatex()
yTit.SetTextSize(0.06)
yTit.SetTextAngle(90)
#yTit.DrawLatex(-25,my_ghist.GetMinimum()+(1.0-my_ghist.GetMinimum())*0.2,"#sigma(tW')/(#sigma(tW')+#sigma(#bar{t}W'))")
yTit.DrawLatex(0.1,ghist.GetMinimum()+(ghist.GetMaximum()-ghist.GetMinimum())*0.25,"Branching Ratio")
#yTit.DrawLatex(0,my_ghist.GetMaximum()*0.05,"#sigma(tZ)/(#sigma(tZ)+#sigma(#bar{t}Z))")

xTit.SetTextSize(0.07)
xTit.DrawLatex(0.5,ghist.GetMaximum()+(ghist.GetMaximum()-ghist.GetMinimum())*0.05,"LO Z'_{TC2} Branching Ratios")

xTit.SetTextSize(0.03)
xTit.DrawLatex(0.5,ghist.GetMinimum()+(ghist.GetMaximum()-ghist.GetMinimum())*0.20,"hep-ph/9911288 - Harris et. al, Model IV, f_{1}=1.0, f_{2}=0.0")
xTit.DrawLatex(0.5,ghist.GetMinimum()+(ghist.GetMaximum()-ghist.GetMinimum())*0.14,"Including correction from Ferrando and Frandsen used in")
xTit.DrawLatex(0.5,ghist.GetMinimum()+(ghist.GetMaximum()-ghist.GetMinimum())*0.08,"Eur. Phys. J. C (2012) 72, 2072 - Harris and Jain")

leg=TLegend(0.55,0.75,0.85,0.85)

leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextSize(0.05)
leg.SetNColumns(2)

leg.AddEntry(graph_up,"u#bar{u}","l")
leg.AddEntry(graph_top,"t#bar{t}","l")
leg.AddEntry(graph_down,"d#bar{d}","l")
leg.AddEntry(graph_bottom,"b#bar{b}","l")


leg.Draw()

c1.RedrawAxis()

c1.Print("ZpBRs.eps")
c1.Print("ZpBRs.pdf")
c1.Print("ZpBRs.png")
c1.Print("ZpBRs.C")

