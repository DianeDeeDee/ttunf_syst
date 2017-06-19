from ROOT import TH1F
from ROOT import TH2F
from ROOT import TGraphErrors
from ROOT import TCanvas
from ROOT import TPad
from ROOT import TLine
from ROOT import TText
from ROOT import gStyle

def plotfragment(m,aname):
  m_listofkeys = m.keys()
  m_listofkeys.sort()
  
  print "PrintFitResultPlot.py!!!!!!!!!!!!!!!!!!!!!!!"
  
  #Write pois (pre- and post-fit) with errors and variations to file
  out_file = open('poiNumbers.txt','w')
  
  #Define canvas and pads
  canvas = TCanvas("canvas","gerrors2",1024,1448)
  canvas.GetFrame().SetBorderSize(12)

  pad1 = TPad("pad1", "pad1", 0.0  , 0.0  , 1.0 , 1.0  , 0)
  pad1.SetLeftMargin(0.2)
  pad1.SetRightMargin(0.05)
  pad1.SetBottomMargin(0.09)
  pad1.Draw()
  pad1.cd()

  #Set histogram to get the range and nuisance parameter labels correct
  h = TH2F("h", "", 1, -2.5, 2.5, len(m_listofkeys), 0., len(m_listofkeys))

  count = 0
  for name in m_listofkeys:
    count += 1
    h.GetYaxis().SetBinLabel(count, name)
    print name
    
  h.LabelsOption("h")
  labelSize = 1./len(m_listofkeys)
  h.SetLabelSize(labelSize,"Y")
  h.GetXaxis().SetLabelColor(1)
  h.GetXaxis().SetAxisColor(1)
  h.GetXaxis().SetTickLength(0.)
  h.GetYaxis().SetLabelColor(1)
  h.GetYaxis().SetAxisColor(1)
  h.GetYaxis().SetTickLength(0.)
  h.SetStats(0)
  h.Draw("h")
  
  #Define histograms for colour bands
  h1sig = TH1F("h1sig", "", 1, -2., 2.)
  h1sig.SetBinContent(1,len(m_listofkeys))
  h1sig.SetFillColor(3)
  h1sig.SetMaximum(len(m_listofkeys))
  h1sig.SetMinimum(0.)
  h1sig.Draw("same")
  
  h2sig = TH1F("h2sig", "", 1, -1., 1.)
  h2sig.SetBinContent(1,len(m_listofkeys))
  h2sig.SetFillColor(5)
  h2sig.SetMaximum(len(m_listofkeys))
  h2sig.SetMinimum(0.)
  h2sig.Draw("same")
  

  #Define and fill TGraph
  graph = TGraphErrors(len(m_listofkeys))
   
  count = 0
  for name in m_listofkeys:
      graph.SetPoint(count, m[name][4], count+0.5)
      graph.SetPointError(count, m[name][5], 0.)
      count += 1
      
      out_file.write( name+"\t"+str(m[name][0])+"\t"+str(m[name][1])+"\t"+str(m[name][2])+"\t"+str(m[name][3])+"\t"+str(m[name][4])+"\t"+str(m[name][5])+"\t"+str(m[name][6])+"\t"+str(m[name][7])+"\n" )

  graph.SetMarkerStyle(21)
  graph.Draw("P same")

  #Draw horizontal 1- and 2-sigma lines
  ylow= 0.
  yhigh=len(m_listofkeys)
  line1=TLine(-1., ylow, -1., yhigh)
  line2=TLine(-2., ylow, -2., yhigh)
  line1.SetLineStyle(2)
  line2.SetLineStyle(2)
  line1.Draw()
  line2.Draw()
  line1.DrawLine(1., ylow, 1., yhigh)
  line1.DrawLine(0., ylow, 0., yhigh)
  line2.DrawLine(2., ylow, 2., yhigh)
  
  gStyle.SetOptStat(0)
  
  #Save results
  canvas.Print("fitresult_" + aname + ".pdf")
  canvas.Print("fitresult_" + aname + ".eps")
  canvas.Print("fitresult_" + aname + ".png")
  
  
  #Close outputfile
  out_file.close()
  
  return (graph, canvas)

