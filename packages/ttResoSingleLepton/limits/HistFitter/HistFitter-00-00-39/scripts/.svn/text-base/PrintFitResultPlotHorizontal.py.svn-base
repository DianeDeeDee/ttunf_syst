from ROOT import TH1F
from ROOT import TCanvas
from ROOT import TLine
from ROOT import gStyle

def plotfragment(m,aname):
  m_listofkeys = m.keys()
  m_listofkeys.sort()
  
  print "PrintFitResultPlot.py!!!!!!!!!!!!!!!!!!!!!!!"
  
  #graph = TGraphErrors(len(m_listofkeys))
  graph = TH1F("fitres"+aname, "", len(m_listofkeys), -0.5, len(m_listofkeys)-0.5)
  count = 0
  for name in m_listofkeys:
      count += 1
      #graph.SetPoint(count, count, m[name][4])
      #graph.SetPointError(count, count, m[name][5])
      graph.SetBinContent(count, m[name][4])
      graph.SetBinError(count, m[name][5])
      graph.GetXaxis().SetBinLabel(count, name)
      #printname = name
      #printname = printname.replace('syserr_','')
      #printname = printname.replace('_','\_')

      #phantom = ""
      #if m[name][4]>=0: phantom = "\phantom{-}"
      #if str(("%.2f" %m[name][4]))=="0.00": phantom = "\phantom{-}"

      #tableline += "\n" + printname + " & $" + str(("%.2f" %m[name][0])) + "\\pm " + str(("%.2f" %m[name][1])) + "$  & $ " + phantom + str(("%.2f" %m[name][4])) + "\\pm " + str(("%.2f" %m[name][5])) + "$  \\\\"
      #tableline += '''
#%%'''
      #pass
  canvas = TCanvas("canFitres"+aname)
  graph.Draw()
  graph.GetYaxis().SetRangeUser(-2., 2.)
  graph.SetMarkerStyle(21)
  xlow= -0.5
  xhigh=len(m_listofkeys)-0.5
  line1=TLine(xlow, -1., xhigh, -1.)
  line1.SetLineStyle(2)
  line1.Draw()
  line1.DrawLine(xlow, 1., xhigh, 1.)
  line1.DrawLine(xlow, 0., xhigh, 0.)
  gStyle.SetOptStat(0)
  canvas.Print("fitresult_" + aname + ".pdf")
  canvas.Print("fitresult_" + aname + ".eps")
  canvas.Print("fitresult_" + aname + ".png")
  
  return (graph, canvas)

