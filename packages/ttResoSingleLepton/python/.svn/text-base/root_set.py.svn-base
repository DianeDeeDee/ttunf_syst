from ROOT import TStyle
zeus_pub = TStyle("zeus_pub", "Style for ZEUS Publications");
# stuff from "plain" style
zeus_pub.SetPadColor(0)
zeus_pub.SetCanvasColor(0)
zeus_pub.SetTitleColor(0)
zeus_pub.SetStatColor(0)

zeus_pub.SetFrameBorderMode(0)
zeus_pub.SetCanvasBorderMode(0)
zeus_pub.SetPadBorderMode(0)
# Your stuff
zeus_pub.SetOptStat(000000)
zeus_pub.SetOptTitle(000000)
zeus_pub.SetLabelSize(0.06,"x")
zeus_pub.SetLabelSize(0.06,"y")
zeus_pub.SetLabelOffset(0.02,"x")
zeus_pub.SetLabelOffset(0.02,"y")
zeus_pub.SetTitleOffset(0.04,"x")
zeus_pub.SetTitleOffset(0.04,"y")
zeus_pub.SetLabelFont(22,"x")
zeus_pub.SetLabelFont(22,"y")
zeus_pub.SetErrorX(0.0000)
zeus_pub.SetTickLength(0.05,"x")
zeus_pub.SetTickLength(0.05,"y")
#zeus_pub.SetLineWidth(0.8)
zeus_pub.SetLineWidth(1)
zeus_pub.SetPadTickX(1)
zeus_pub.SetPadTickY(1)
zeus_pub.SetPadLeftMargin(0.15)
zeus_pub.SetPadBottomMargin(0.2)



# based on a style file from BaBar
#

#..BABAR style from RooLogon.C in workdir
atlasStyle= TStyle("ATLAS","Atlas style");




# use plain black on white colors
icol=0;
atlasStyle.SetFrameBorderMode(icol);
atlasStyle.SetCanvasBorderMode(icol);
atlasStyle.SetPadBorderMode(icol);
atlasStyle.SetPadColor(icol);
atlasStyle.SetCanvasColor(icol);
atlasStyle.SetStatColor(icol);
#atlasStyle.SetFillColor(icol);

# set the paper & margin sizes
atlasStyle.SetPaperSize(20,26);
atlasStyle.SetPadTopMargin(0.05);
atlasStyle.SetPadRightMargin(0.05);
atlasStyle.SetPadBottomMargin(0.16);
atlasStyle.SetPadLeftMargin(0.18);

# use large fonts
#Int_t font=72;
font=42;
tsize=0.05;
atlasStyle.SetTextFont(font);


atlasStyle.SetTextSize(tsize);
atlasStyle.SetLabelFont(font,"x");
atlasStyle.SetTitleFont(font,"x");
atlasStyle.SetLabelFont(font,"y");
atlasStyle.SetTitleFont(font,"y");
atlasStyle.SetLabelFont(font,"z");
atlasStyle.SetTitleFont(font,"z");

atlasStyle.SetLabelSize(tsize,"x");
atlasStyle.SetTitleSize(tsize,"x");
atlasStyle.SetLabelSize(tsize,"y");
atlasStyle.SetTitleSize(tsize,"y");
atlasStyle.SetLabelSize(tsize,"z");
atlasStyle.SetTitleSize(tsize,"z");


#use bold lines and markers
atlasStyle.SetMarkerStyle(20);
atlasStyle.SetMarkerSize(1.2);
atlasStyle.SetHistLineWidth(2);
atlasStyle.SetLineStyleString(2,"[12 12]"); # postscript dashes

#get rid of X error bars and y error bar caps
atlasStyle.SetErrorX(0.001);

atlasStyle.SetTitleOffset(1.6,"y")

#do not display any of the standard histogram decorations
atlasStyle.SetOptTitle(0);
#atlasStyle.SetOptStat(1111);
atlasStyle.SetOptStat(0);
#atlasStyle.SetOptFit(1111);
atlasStyle.SetOptFit(0);

# put tick marks on top and RHS of plots
atlasStyle.SetPadTickX(1);
atlasStyle.SetPadTickY(1);
 
# gROOT.SetStyle("Plain");


 

#gROOT.SetStyle("zeus_pub")




