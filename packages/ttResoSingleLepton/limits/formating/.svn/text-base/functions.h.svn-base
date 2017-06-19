#include <iostream> 
#include <fstream> 
#include <sstream> 
#include <vector> 
#include <iomanip> 
#include "TString.h" 
#include "TROOT.h" 
#include "TStyle.h" 
#include "TGraph.h" 
#include "TGraphErrors.h" 
#include "TGraphAsymmErrors.h" 
#include "TCanvas.h" 
#include "TLegend.h" 
#include "TAxis.h" 
#include "TLatex.h" 
#include "TPaveText.h" 
#include <assert.h> 

using namespace std; 


void getMass(TString selection, double *mass, const int Nmasses) { 

  double dummy[9]={500,600,700,800,1000,1300,1600,2000,3000}; 
  if(selection=="syst4j" || selection=="systdr" || selection=="newlum4j" || selection=="newlumdr") { 
    dummy[0]=500; 
    dummy[1]=600; 
    dummy[2]=700; 
    dummy[3]=800; 
    dummy[4]=1000; 
  } 
  else if(selection.Contains("r164j") ||selection.Contains("r16dr") ) { 
    dummy[0]=500; 
    dummy[1]=700; 
    dummy[2]=1000; 
  } 
  else if(selection.Contains("Zprime-lims")){ 
    dummy[0]=500; 
    dummy[1]=600; 
    dummy[2]=700; 
    dummy[3]=800; 
    dummy[4]=1000; 
  } 
  else if(selection.Contains("BH-lim")){ 
    //dummy[0]=750; 
    //dummy[1]=1500; 
    //dummy[2]=2250; 
    dummy[0]=750; 
    dummy[1]=1500; 
    dummy[2]=2000; 
    dummy[3]=2250; 
    dummy[4]=2500; 
    //dummy[5]=2750; 
  } 
  else if(selection.Contains("qbh")){ 
    //dummy[0]=750; 
    //dummy[1]=1500; 
    //dummy[2]=2250; 
    dummy[0]=750; 
    dummy[1]=1500; 
    //dummy[2]=2000; 
    dummy[2]=2250; 
    dummy[3]=2500; 
    dummy[4]=2750; 
  } 
  else if(selection.Contains("_sc")){ 
    dummy[0]=500; 
    dummy[1]=600; 
    dummy[2]=700; 
    dummy[3]=800; 
    dummy[4]=1000; 
    dummy[5]=1300; 
    dummy[6]=1600; 
    dummy[7]=2000; 
    dummy[8]=3000; 
  } 
  else if(selection.Contains("_th")){ 
    dummy[0]=700; 
    dummy[1]=1000; 
    dummy[2]=1600; 
  } 
  else if(selection.Contains("_Zp")){ 
    dummy[0]=500; 
    dummy[1]=600; 
    dummy[2]=700; 
    dummy[3]=800; 
    dummy[4]=1000; 
    dummy[5]=1300; 
    dummy[6]=1600; 
    dummy[7]=2000; 
  } 
  else if(selection.Contains("_kkgsc")){ 
    dummy[0]=500; 
    dummy[1]=600; 
    dummy[2]=700; 
    dummy[3]=800; 
    dummy[4]=1000; 
    dummy[5]=1300; 
    dummy[6]=1600; 
    dummy[7]=1800; 
  } 
  else if(selection.Contains("_KKg")){ 
    dummy[0]=500; 
    dummy[1]=600; 
    dummy[2]=700; 
    dummy[3]=800; 
    dummy[4]=1000; 
    dummy[5]=1300; 
    dummy[6]=1600; 
  } 
  else if(selection.Contains("kkgox")){ 
    dummy[0] =700; 
    dummy[1] =800; 
    dummy[2] =1000; 
    dummy[3] =1150; 
    dummy[4] =1300; 
    dummy[5] =1600; 
    dummy[6] =1800; 
    dummy[7] =2000; 

  } 
  else if(selection.Contains("kkgazox")){ 
    dummy[0]=700; 
    dummy[1]=800; 
    dummy[2]=900; 
    dummy[3]=1000; 
    dummy[4]=1150; 
    dummy[5]=1300; 
    dummy[6]=1600; 
    dummy[7]=1800; 
    dummy[8]=2000; 
  } 
  else if(selection.Contains("_kkgaz")){ 
    dummy[0]=700; 
    dummy[1]=800; 
    dummy[2]=900; 
    dummy[3]=1000; 
    dummy[4]=1150; 
    dummy[5]=1300; 
    dummy[6]=1600; 
    dummy[7]=1800; 
    dummy[8]=2000; 
  } 
  else if(selection.Contains("azox")){ 
    dummy[0]=500; 
    dummy[1]=600; 
    dummy[2]=700; 
    dummy[3]=800; 
    dummy[4]=1000; 
    dummy[5]=1300; 
    dummy[6]=1600; 
    dummy[7]=2000; 
    dummy[8]=3000; 
  } 

  for(int i=0;i<Nmasses;++i) { 
    mass[i]=dummy[i]; 
  } 

} 


void getXsec(TString selection, double *xsec, const int Nmasses); 
void getXsec(TString selection, double *xsec, const int Nmasses) { 

  //double dummy[5]={500,700,1000,1600,2000}; 
  //double dummy[5]={3.552, 1.373, 0.310, 0.027, 0.0063}; 
  double dummy[5]; 
  dummy[0] = 3.552; dummy[1] = 2.246; 
  dummy[2] = 1.373; dummy[3] = 0.826; 
  dummy[4] = 0.310; 

  if(selection=="10tev") { 
    dummy[0]=1; 
    dummy[1]=1; 
    dummy[2]=1; 
  } 
  if(selection=="venkat") { 
    dummy[0]=3.552; 
    dummy[1]=1.373; 
    dummy[2]=0.310; 
  } 
  if(selection=="clermond") { 
    dummy[0] = 35.52;  
    dummy[1] = 137.3;  
    dummy[2] = 31.0;  
  } 
  if(selection=="spano") { 
    dummy[0] = 1.373;  
    dummy[1] = 0.310;  
    dummy[2] = 0.027;  
  } 
  if(selection=="lucia4jet" || selection=="luciadrcut") { 
    dummy[0] = 3.552;  
    dummy[1] = 1.373;  
    dummy[2] = 0.310;  
  } 
  if(selection=="wuppertal") { 
    dummy[0] = 3.552;  
    dummy[1] = 1.373;  
    dummy[2] = 0.310;  
    dummy[3] = 0.027;  
    dummy[4] = 0.0063;  
  } 
  if(selection.Contains("syst4j") || selection.Contains("systdr") || selection.Contains("newlum4j") || selection.Contains("newlumdr")) { 
    dummy[0] = 3.552;  
    dummy[1] = 2.246; 
    dummy[2] = 1.373; 
    dummy[3] = 0.826; 
    dummy[4] = 0.310; 
  } 
  if(selection.Contains("r164j") || selection.Contains("r16dr")) { 
    dummy[0] = 3.552;  
    dummy[1] = 1.373; 
    dummy[2] = 0.310; 
  } 

  for(int i=0;i<Nmasses;++i) { 
    xsec[i]=dummy[i]; 
  } 


} 

TString getName(TString selection); 
TString getName(TString selection) { 

  TString Name=""; 
  if(selection=="10tev") Name="10 TeV MC";   
  else if(selection.Contains("-lims")) Name="CLs, dRmin"; 
  else if(selection.Contains("-lim_dr")) Name="dRmin"; 
  else if(selection.Contains("-lim_4j")) Name="4 jets"; 
  else if(selection.Contains("-lim_nc")) Name="#chi^{2} no cal"; 
  else if(selection.Contains("_sc4j")) Name="GCBT, 4 jets"; 
  else if(selection.Contains("_scdr")) Name="GCBT, dRmin"; 
  else if(selection.Contains("_scnc")) Name="GCBT, #chi^{2} no cal"; 
  else if(selection.Contains("_scx2")) Name="GCBT, #chi^{2} rescale"; 
  else if(selection.Contains("_kkgsc4j")) Name="GCBT, 4 jets"; 
  else if(selection.Contains("_kkgscdr")) Name="GCBT, dRmin"; 
  else if(selection.Contains("_kkgscnc")) Name="GCBT, #chi^{2} no cal"; 
  else if(selection.Contains("_kkgscx2")) Name="GCBT, #chi^{2} rescale"; 
  else if(selection.Contains("_th0")) Name="Cut 0"; 
  else if(selection.Contains("_th1")) Name="Cut 1"; 
  else if(selection.Contains("_th2")) Name="Cut 2"; 
  else if(selection.Contains("_th3")) Name="Cut 3"; 
  else if(selection.Contains("_th4")) Name="Cut 4"; 
  else if(selection.Contains("_th5")) Name="Cut 5"; 
  else if(selection.Contains("ox")) Name="Oxford"; 
  else if(selection.Contains("az")) Name="Arizona"; 


  return Name; 
} 

TGraph * getErrorBand(const int N, double* x, double* yUp, double* yDown); 
TGraph * getErrorBand(const int N, double* x, double* yUp, double* yDown) { 
  // Returns a TGraph with the error band. 
  //const int N=x.size;//(x); 
  //const int Ne = (N>1 ? 2*N+1 : 2*N); 
  assert(N>0); 

  const int Ne = 2*N+1; 

  double xe[Ne], ye[Ne]; 
  if(N>1) { 
    for(int i=0;i<N;i++) { 
      xe[i] = x[i]; 
      xe[2*N-1-i] = x[i]; 
      ye[i] = yUp[i]; 
      ye[2*N-1-i] = yDown[i]; 
    } 
    xe[Ne-1] = xe[0]; 
    ye[Ne-1] = yUp[0]; 
  } 
  else { 
    xe[0] = x[0];      
    xe[1] = x[0]; 
    xe[2] = x[0]; 
    ye[0] = yUp[0]; 
    ye[1] = yDown[0]; 
    ye[2] = yUp[0]; 
  } 

  TGraph *errorBand = new TGraph(Ne,xe,ye); 

  //cout << "Point  x   y" << endl; 
  //for(int i=0;i<Ne;++i) { 
  //  cout << i << "  " << xe[i] << "  " << ye[i] << endl; 
  //} 

  return errorBand; 

} 

TGraph * getErrorBand(const int N, double* x, double* yCentral, double* yPlus, double* yMinus); 
TGraph * getErrorBand(const int N, double* x, double* yCentral, double* yPlus, double* yMinus) { 
  // Returns a TGraph with the error band. 
  //yCental is the central value, with plus and minus limits 
  //const int N=x.size;//(x); 
  //const int Ne = (N>1 ? 2*N+1 : 2*N); 
  assert(N>0); 

  const int Ne = 2*N+1; 

  double xe[Ne], ye[Ne]; 
  if(N>1) { 
    for(int i=0;i<N;i++) { 
      xe[i] = x[i]; 
      xe[2*N-1-i] = x[i]; 
      ye[i] = yCentral[i] + yPlus[i]; 
      ye[2*N-1-i] = yCentral[i] - fabs(yMinus[i]); 
    } 
    xe[Ne-1] = xe[0]; 
    ye[Ne-1] = yCentral[0] + yPlus[0]; 
  } 
  else { 
    xe[0] = x[0];      
    xe[1] = x[0]; 
    xe[2] = x[0]; 
    ye[0] = yCentral[0] + yPlus[0]; 
    ye[1] = yCentral[0] - fabs(yMinus[0]); 
    ye[2] = yCentral[0] + yPlus[0]; 
  } 

  TGraph *errorBand = new TGraph(Ne,xe,ye); 

  //cout << "Point  x   y" << endl; 
  //for(int i=0;i<Ne;++i) { 
  //  cout << i << "  " << xe[i] << "  " << ye[i] << endl; 
  //} 

  return errorBand; 

} 

TString getSystName(TString systematics); 
TString getSystName(TString systematics) { 

  TString Name=systematics; 
  if(systematics=="") Name ="Nominal";  
  if(systematics=="jeru") Name ="JER up";  
  if(systematics=="jerd") Name ="JER 'down'";  
  if(systematics=="jesu") Name ="JES+JMS+d12s up";  
  if(systematics=="jesd") Name ="JES+JMS+d12s down";  
  if(systematics=="rjesu") Name ="JES up"; //resolved 
  if(systematics=="rjesd") Name ="JES down";  
  if(systematics=="a10jesu") Name ="a10 JES/R+JMS/R up";  
  if(systematics=="a10jesd") Name ="a10 JES/R+JMS/R down";  
  if(systematics=="bjesu") Name ="BJES up"; 
  if(systematics=="bjesd") Name ="BJES down"; 
  if(systematics=="jrec") Name ="Jet reco"; 
  if(systematics=="jeffu") Name ="Jet eff up"; 
  if(systematics=="jeffd") Name ="Jet eff 'down'"; 
  if(systematics=="btagd") Name ="Btag down";  
  if(systematics=="btagu") Name ="Btag up";  
  if(systematics=="midptd") Name ="MuID pt down";  
  if(systematics=="midptu") Name ="MuID pt up";  
  if(systematics=="mmsptd") Name ="MuMS pt down";  
  if(systematics=="mmsptu") Name ="MuMS pt up"; 
  if(systematics=="musfd") Name ="MuSF down"; 
  if(systematics=="musfu") Name ="MuSF up"; 
  if(systematics=="muidu") Name ="Mu ID SF up"; 
  if(systematics=="muidd") Name ="Mu ID SF down"; 
  if(systematics=="muridu") Name ="Mu res ID up"; 
  if(systematics=="muridd") Name ="Mu res ID down"; 
  if(systematics=="murmsu") Name ="Mu res MS up"; 
  if(systematics=="murmsd") Name ="Mu res MS down"; 
  if(systematics=="mutrgu") Name ="Mu trig up"; 
  if(systematics=="mutrgd") Name ="Mu trig down"; 
  if(systematics=="muptu") Name ="Mu mom scale up"; 
  if(systematics=="muptd") Name ="Mu mom scale down"; 
  if(systematics=="mutrigu") Name ="Mu trig up"; 
  if(systematics=="mutrigd") Name ="Mu trig down"; 
  if(systematics=="murecu") Name ="Mu reco up"; 
  if(systematics=="murecd") Name ="Mu reco down"; 
  if(systematics=="muresmsu") Name ="Mu res MS up"; 
  if(systematics=="muresmsd") Name ="Mu res MS down"; 
  if(systematics=="muresidu") Name ="Mu res ID up"; 
  if(systematics=="muresidd") Name ="Mu res ID down"; 
  if(systematics=="vvxsd") Name ="VV xs down"; 
  if(systematics=="vvxsu") Name ="VV xs up"; 
  if(systematics=="stxsd") Name ="Sing top xs down"; 
  if(systematics=="stxsu") Name ="Sing top xs up"; 
  if(systematics=="ttxsd") Name ="ttbar xs down"; 
  if(systematics=="ttxsu") Name ="ttbar xs up"; 
  if(systematics=="lumd") Name ="Lumi down"; 
  if(systematics=="lumu") Name ="Lumi up"; 
  if(systematics=="elsfd") Name ="ElSF down"; 
  if(systematics=="elsfu") Name ="ElSF up"; 
  if(systematics=="elidd") Name ="ElID down"; 
  if(systematics=="elidu") Name ="ElID up"; 
  if(systematics=="elesd") Name ="El ES down"; 
  if(systematics=="elesu") Name ="El ES up"; 
  if(systematics=="elerd") Name ="El ER down"; 
  if(systematics=="eleru") Name ="El ER up"; 
  if(systematics=="elptd") Name ="El pt down"; 
  if(systematics=="elptu") Name ="El pt up"; 
  if(systematics=="elresu") Name ="El res up"; 
  if(systematics=="eleresu") Name ="El ER up"; 
  if(systematics=="eleresd") Name ="El ER down"; 
  if(systematics=="eltrgu") Name ="El trig up"; 
  if(systematics=="eltrgd") Name ="El trig down"; 
  if(systematics=="elescu") Name ="El ES up"; 
  if(systematics=="elescd") Name ="El ES down"; 
  if(systematics=="fsrd") Name ="FSR down"; 
  if(systematics=="fsru") Name ="FSR up"; 
  if(systematics=="isrd") Name ="ISR down"; 
  if(systematics=="isru") Name ="ISR up"; 
  if(systematics=="ifsrd") Name ="I/FSR down"; 
  if(systematics=="ifsru") Name ="I/FSR up"; 
  if(systematics=="wjsfd") Name ="Wjet SF down"; 
  if(systematics=="wjsfu") Name ="Wjet SF up"; 
  if(systematics=="wjhf") Name ="Wjet HF"; 
  if(systematics=="whfcu") Name ="W HF c up"; 
  if(systematics=="whfcd") Name ="W HF c down"; 
  if(systematics=="whfbbccu") Name ="W HF bbcc up"; 
  if(systematics=="whfbbccd") Name ="W HF bbcc down"; 
  if(systematics=="wjallu") Name ="Wjet all up"; 
  if(systematics=="wjalld") Name ="Wjet all down"; 
  if(systematics=="qcde") Name ="QCD anti-ele"; 
  if(systematics=="crosscheck") Name ="Nom. (Mainz)"; 
  if(systematics=="ttnoru") Name ="tt norm up"; 
  if(systematics=="ttnord") Name ="tt norm down"; 
  if(systematics=="scaleu") Name ="Scale e up"; 
  if(systematics=="scaled") Name ="Scale e down"; 
  if(systematics=="scalmu") Name ="Scale mu up"; 
  if(systematics=="scalmd") Name ="Scale mu down"; 
  if(systematics=="wnoru") Name ="W norm up"; 
  if(systematics=="wnord") Name ="W norm down"; 
  if(systematics=="qcdnoru") Name ="QCD norm up"; 
  if(systematics=="qcdnord") Name ="QCD norm down"; 
  if(systematics=="znormu") Name ="Z norm up"; 
  if(systematics=="znormd") Name ="Z norm down"; 
  if(systematics=="stnoru") Name ="s-top norm up"; 
  if(systematics=="stnord") Name ="s-top norm down"; 
  if(systematics=="dibnoru") Name ="di-bos norm up"; 
  if(systematics=="dibnord") Name ="di-bos norm down"; 
  if(systematics=="iqoptu") Name ="Wiqopt3 up"; 
  if(systematics=="iqoptd") Name ="Wiqopt3 'down'"; 
  if(systematics=="wptjminu") Name ="Wptjmin10 up"; 
  if(systematics=="wptjmind") Name ="Wptjmin10 'down'"; 
  if(systematics=="wqktfacu") Name ="WqfacKtfac up"; 
  if(systematics=="wqktfacd") Name ="WqfacKtfac 'down'"; 
  if(systematics=="ptjminu") Name ="W ptjmin10 up"; 
  if(systematics=="ptjmind") Name ="W ptjmin10 'down'"; 
  if(systematics=="mcgenu") Name ="tt MC gen up"; 
  if(systematics=="mcgend") Name ="tt MC gen 'down'"; 
  if(systematics=="fgpsu") Name ="Frag+PS up"; 
  if(systematics=="fgpsd") Name ="Frag+PS 'down'"; 
  if(systematics=="fragu") Name ="Frag up"; 
  if(systematics=="fragd") Name ="Frag 'down'"; 
  if(systematics=="psu") Name ="PS up"; 
  if(systematics=="psd") Name ="PS 'down'"; 
  if(systematics=="pdfu") Name ="PDF up"; 
  if(systematics=="pdfd") Name ="PDF down"; 
  if(systematics=="qcdshu") Name ="QCD shape up"; 
  if(systematics=="qcdshd") Name ="QCD shape 'down'"; 
  if(systematics=="qcdpuu") Name ="QCD PU up"; 
  if(systematics=="qcdpud") Name ="QCD PU 'down'"; 
  if(systematics=="qcdhtu") Name ="QCD HT up"; 
  if(systematics=="qcdhtd") Name ="QCD HT 'down'"; 
  if(systematics=="puu") Name ="PU up"; 
  if(systematics=="pud") Name ="PU 'down'"; 
  if(systematics=="lholeu") Name ="LAr hole up"; 
  if(systematics=="lholed") Name ="LAr hole 'down'"; 
  if(systematics=="larhu") Name ="LAr hole up"; 
  if(systematics=="larhd") Name ="LAr hole down"; 
  if(systematics=="metcou") Name ="MET cell out up"; 
  if(systematics=="metcod") Name ="MET cell out 'down'"; 
  if(systematics=="") Name =""; 

  //if(systematics=="lumd1") Name ="Lumi down 1%"; 
  //if(systematics=="lumd2") Name ="Lumi down 2%"; 
  //if(systematics=="lumd3") Name ="Lumi down 3%"; 
  //if(systematics=="lumd4") Name ="Lumi down 4%"; 
  //if(systematics=="lumd5") Name ="Lumi down 5%"; 
  //if(systematics=="lumd6") Name ="Lumi down 6%"; 
  //if(systematics=="lumd7") Name ="Lumi down 7%"; 
  //if(systematics=="lumd8") Name ="Lumi down 8%"; 
  //if(systematics=="lumd9") Name ="Lumi down 9%"; 
  //if(systematics=="lumd10") Name ="Lumi down 10%"; 
  //if(systematics=="lumd11") Name ="Lumi down 11%"; 
  // 
  //if(systematics=="lumu2") Name ="Lumi up 2.2%"; 
  //if(systematics=="lumd2") Name ="Lumi down 2.2%"; 
  //if(systematics=="lumu6") Name ="Lumi up 6.6%"; 
  //if(systematics=="lumd6") Name ="Lumi down 6.6%"; 
  //if(systematics=="lumu10") Name ="Lumi up 11%"; 
  //if(systematics=="lumd10") Name ="Lumi down 11%"; 
  //if(systematics=="lumu12") Name ="Lumi up 13.2%"; 
  //if(systematics=="lumd12") Name ="Lumi down 13.2%"; 
  //if(systematics=="lumu16") Name ="Lumi up 17.6%"; 
  //if(systematics=="lumd16") Name ="Lumi down 17.6%"; 
  //if(systematics=="lumu20") Name ="Lumi up 22%"; 
  //if(systematics=="lumd20") Name ="Lumi down 22%"; 
  //if(systematics=="lumu22") Name ="Lumi up 24.2%"; 
  //if(systematics=="lumd22") Name ="Lumi down 24.2%"; 
  //if(systematics=="lumu26") Name ="Lumi up 28.6%"; 
  //if(systematics=="lumd26") Name ="Lumi down 28.6%"; 
  //if(systematics=="lumu30") Name ="Lumi up 33%"; 
  //if(systematics=="lumd30") Name ="Lumi down 33%"; 
 


  return Name; 

} 

void ATLAS_LABEL(Double_t x,Double_t y,Color_t color,TCanvas* cv)  
{ 
  if (cv==NULL) return; 
  cv->cd(); 
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize);  
  l.SetNDC(); 
  l.SetTextFont(72); 
  l.SetTextColor(color); 
  l.DrawLatex(x,y,"ATLAS"); 
  return; 
} 

void ATLASPRELIM_LABEL(Double_t x,Double_t y,Color_t color,TCanvas* cv)  
{ 
  if (cv==NULL) return; 
  cv->cd(); 
  TLatex* l = new TLatex(x,y,""); //l.SetTextAlign(12); l.SetTextSize(tsize);  
  l->SetNDC(); 
  l->SetTextFont(72); 
  l->SetTextColor(color); 
  l->DrawLatex(x,y,"ATLAS"); 
  TLatex* l2 = new TLatex(x,y,""); //l.SetTextAlign(12); l.SetTextSize(tsize);  
  l2->SetNDC(); 
  l2->SetTextColor(color); 
  l2->DrawLatex(x+0.175,y,"Preliminary"); 
  return; 
} 

void ATLASLabel(Double_t x,Double_t y,TCanvas* cv,bool Preliminary=false,Color_t color=1) { 
  if (cv==NULL) return; 
  cv->cd(); 
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize);  
  l.SetNDC(); 
  l.SetTextFont(72); 
  l.SetTextColor(color); 

  double delx = 0.115*696*gPad->GetWh()/(472*gPad->GetWw()); 
  double dely = 0.029*696*gPad->GetWh()/(472*gPad->GetWw()); 
  TPaveText *bk = new TPaveText(x,y-0.01,x+delx,y+dely,"NDC"); 
  //cout << "gPad->GetWh()=" << gPad->GetWh() << endl; 
  //cout << "gPad->GetWw()=" << gPad->GetWw() << endl; 
  //cout << 696.0*gPad->GetWh()/(472.0*gPad->GetWw()) << endl; 
  //cout << "Coordinates: " << x << ", " << y << endl; 
  //cout << "           : " << x+delx << ", " << y+dely << endl; 
  //bk->AddText("ATLAS"); 
  bk->SetTextColor(kWhite); 
  bk->SetFillColor(kWhite); 
  bk->SetBorderSize(0); 
  bk->Draw(); 
  l.DrawLatex(x,y,"ATLAS"); 
  //l.SetTextFont(42);  l.SetTextSize(0.04);  l.DrawLatex(x,y-0.08,"#sqrt{#font[12]{s}} = 7 TeV");  l.SetTextFont(72);  l.SetTextSize(0.05); 
  if (Preliminary) { 
    TLatex p;  
    p.SetNDC(); 
    p.SetTextFont(42); 
    p.SetTextColor(color); 
    double delx2 = 0.115*696*gPad->GetWh()/(472*gPad->GetWw()); 
    TPaveText *bk2 = new TPaveText(x+delx,y-0.01,x+delx+delx2,y+dely,"NDC"); 
    bk2->SetTextColor(kWhite); 
    bk2->SetFillColor(kWhite); 
    bk2->SetBorderSize(0); 
    bk2->Draw(); 
    p.DrawLatex(x+delx,y,"Preliminary"); 
    //    p.DrawLatex(x,y,"#sqrt{s}=900GeV"); 
  } 
  return; 
} 

void ATLASCMEIntL(Double_t x,Double_t y,Double_t IntL,TCanvas* cv,Color_t color=1) { 
  if (cv==NULL) return; 
  cv->cd(); 
  TLatex l; 
  l.SetNDC(); 
  l.SetTextFont(42); 
  l.SetTextColor(color); 
  l.SetTextSize(0.04); 
  char charIntL[200]; 
  sprintf(charIntL,"%.2f",IntL); 
  //sprintf(charIntL,"%.1f",IntL); 
  //l.DrawLatex(x,y,TString("#scale[0.65]{#int} Ldt = "+string(charIntL)+" nb^{-1}      #sqrt{s} = 7 TeV")); 
  //l.DrawLatex(x,y,TString("#scale[0.6]{#int}_{  }#font[12]{L_{ }dt} = "+string(charIntL)+" pb^{-1}")); 
  l.DrawLatex(x,y,TString("#scale[0.6]{#int}_{  }#font[12]{L_{ }dt} = "+string(charIntL)+" fb^{-1}")); 
  //l.DrawLatex(x,y+0.2,"#sqrt{s} = 7 TeV"); 
  l.DrawLatex(x,y+0.08,"#sqrt{#font[12]{s}} = 7 TeV"); 
  return; 
} 

void ATLASCME(Double_t x,Double_t y,TCanvas* cv,Color_t color=1) { 
  if (cv==NULL) return; 
  cv->cd(); 
  TLatex l; 
  l.SetNDC(); 
  l.SetTextFont(42); 
  l.SetTextColor(color); 
  l.SetTextSize(0.04); 
  l.DrawLatex(x,y,"#sqrt{s} = 7 TeV"); 
  return; 
} 

TString getSystLongName(TString sshort) { 

  TString longlabel=""; 

  if(sshort=="lum")          longlabel="Luminosity                   "; 
  else if(sshort=="larh")    longlabel="LAr ``hole''                 "; 
  else if(sshort=="ttnor")   longlabel="\\ttbar{} normalisation      ";    
  else if(sshort=="ifsr")    longlabel="\\ttbar{} ISR, FSR           "; 
  else if(sshort=="fgps")    longlabel="\\ttbar{} parton shower      "; 
  else if(sshort=="mcgen")   longlabel="\\ttbar{} generator dependence"; 
  else if(sshort=="pdf")     longlabel="PDF uncertainty              "; 
  else if(sshort=="wnor")    longlabel="$W+\\mathrm{jets}$ normalisation"; 
  else if(sshort=="qcdnor")  longlabel="QCD normalisation            "; 
  else if(sshort=="qcdsh")   longlabel="QCD shape                    "; 
  else if(sshort=="znorm")   longlabel="$Z+\\mathrm{jets}$ normalisation"; 
  else if(sshort=="stnor")   longlabel="single top normalisation     "; 
  else if(sshort=="dibnor")  longlabel="Di-boson normalisation       "; 
  else if(sshort=="jes")     longlabel="Jet energy and mass scale    ";  
  else if(sshort=="jer")     longlabel="Jet energy and mass resolution"; 
  else if(sshort=="jeff")    longlabel="Jet reconstruction efficiency"; 
  else if(sshort=="eltrg")   longlabel="Electron trigger SF          "; 
  else if(sshort=="elid")    longlabel="Electron ID and reco eff SF  "; 
  else if(sshort=="elesc")   longlabel="Electron energy scale        "; 
  else if(sshort=="eleres")  longlabel="Electron energy resolution   "; 
  else if(sshort=="mutrg")   longlabel="Muon trigger                 "; 
  else if(sshort=="murec")   longlabel="Muon reco efficiency SF      "; 
  else if(sshort=="mupt")    longlabel="Muon momentum scale          "; 
  else if(sshort=="muid")    longlabel="Muon ID SF                   "; 
  else if(sshort=="murid")   longlabel="Muon momentum resolution (ID)"; 
  else if(sshort=="murms")   longlabel="Muon momentum resolution (MS)"; 
  else if(sshort=="metco")   longlabel="\\MET{} uncertainty           "; 
  else if(sshort=="JES")          longlabel="Jet energy scale     "; 
  else if(sshort=="BoostedJES")   longlabel="JES, \\akt1.0 jets    "; 
  else if(sshort=="JEE")	  longlabel="Jet efficiency       "; 
  else if(sshort=="JetEnerRes")	  longlabel="Akt4 jet energy res. "; 
  else if(sshort=="NormW")	  longlabel="$W$+jets norm        "; 
  else if(sshort=="norm_W")	  longlabel="$W$+jets norm        "; 
  else if(sshort=="WjetsAsym")	  longlabel="W asymmetry          "; 
  else if(sshort=="iqopt3")       longlabel="$W$ shape, ``iqopt3''"; 
  else if(sshort=="ptjmin10")     longlabel="$W$ shape, ``ptjmin10''"; 
  else if(sshort=="Whfsf")        longlabel="W data-driven SFs    ";
  else if(sshort=="Btag")	  longlabel="$b$-tag              "; 
  else if(sshort=="BtagC")	  longlabel="$c$-tag              "; 
  else if(sshort=="BtagL")	  longlabel="Mistag               "; 
  else if(sshort=="Btag6")	  longlabel="$b$-tag EV6          "; 
  else if(sshort=="Btag7")	  longlabel="$b$-tag EV7          "; 
  else if(sshort=="Btag8")	  longlabel="$b$-tag EV8          "; 
  else if(sshort=="Btag9")	  longlabel="$b$-tag EV9          "; 
  else if(sshort=="Btag10")	  longlabel="$b$-tag high \\pt    "; 
  else if(sshort=="EleReco")	  longlabel="Electron reco        "; 
  else if(sshort=="EleTrig")	  longlabel="Electron trigger     "; 
  else if(sshort=="EleID")	  longlabel="Electron ID          "; 
  else if(sshort=="EleSF")	  longlabel="Electron scale factor"; 
  else if(sshort=="EnerScale")	  longlabel="Electron energy scale"; 
  else if(sshort=="EleEnerRes")   longlabel="Electron energy res. "; 
  else if(sshort=="scale_mu")	  longlabel="Muon \\pt\\ scale      "; 
  else if(sshort=="MuScale")	  longlabel="Muon \\pt\\ scale      "; 
  else if(sshort=="MuTrig")	  longlabel="Muon trigger SF      "; 
  else if(sshort=="MuTrigMatch")  longlabel="Muon trigger match   "; 
  else if(sshort=="MuID")	  longlabel="Muon ID              "; 
  else if(sshort=="MuSF")	  longlabel="Muon scale factor    "; 
  else if(sshort=="Mureco")	  longlabel="Muon reco            "; 
  else if(sshort=="MuResMs")	  longlabel="Muon res. MS         "; 
  else if(sshort=="MuResId")	  longlabel="Muon res. ID         "; 
  else if(sshort=="MetCellOut")   longlabel="\\MET\\ cell out       "; 
  else if(sshort=="METPileUp")	  longlabel="\\MET\\ pile-up        "; 
  else if(sshort=="MetPileUp")	  longlabel="\\MET\\ pile-up        "; 
  else if(sshort=="EWS")	  longlabel="\\ttbar\\ EW Sudakov   "; 
  else if(sshort=="Pttt")	  longlabel="\\ttbar\\ pt RW   "; 
  else if(sshort=="MI")	          longlabel="Mini-iso SF          "; 
  else if(sshort=="PDF")	  longlabel="PDF                  "; 
  else if(sshort=="norm_Diboson")    longlabel="Di-boson norm        "; 
  else if(sshort=="norm_single-top") longlabel="Single-top norm      "; 
  else if(sshort=="norm_QCD")     longlabel="QCD norm             "; 
  else if(sshort=="norm_QCDe")    longlabel="Multi-jets norm, $e$+jets  "; 
  else if(sshort=="norm_QCDmu")   longlabel="Multi-jets norm, $\\mu$+jets"; 
  else if(sshort=="norm_Z")	  longlabel="$Z$+jets norm        "; 
  else if(sshort=="norm_tt")	  longlabel="\\ttbar\\ norm         "; 
  else if(sshort=="NormZ")	  longlabel="$Z$+jets norm        "; 
  else if(sshort=="JVFCut")	  longlabel="Jet vertex fraction  "; 
  else if(sshort=="PartonShower") longlabel="\\ttbar\\ Parton shower"; 
  else if(sshort=="luminosity")   longlabel="Luminosity           "; 
  else if(sshort=="IFSR")         longlabel="\\ttbar\\ ISR/FSR      "; 
  else if(sshort=="MCGen")        longlabel="\\ttbar\\ generator    "; 
  else if(sshort=="MCGen")        longlabel="\\ttbar\\ generator    "; 
  else if(sshort=="ScaleTtbar")   longlabel="\\ttbar\\ gen. (scale) "; 

 
  else if(sshort=="BoostedJES0" ) longlabel ="Akt10 JES (Gamma-jet)"; 
  else if(sshort=="BoostedJES1" ) longlabel ="Akt10 JES (GJ:dataMC)   "; 
  else if(sshort=="BoostedJESothers" ) longlabel ="Akt10 JES (GJ:syst)   "; 
  else if(sshort=="BoostedJES2" ) longlabel ="Akt10 JES2 (GJ:pt2Cut)   "; 
  else if(sshort=="BoostedJES3" ) longlabel ="Akt10 JES3 (GJ:dPhiCut)  "; 
  else if(sshort=="BoostedJES4" ) longlabel ="Akt10 JES4 (GJ:photonPurity)"; 
  else if(sshort=="BoostedJES5" ) longlabel ="Akt10 JES5 (GJ:PES)      "; 
  else if(sshort=="BoostedJES6" ) longlabel ="Akt10 JES6 (GJ:photonPurity)"; 
  else if(sshort=="BoostedJES7" ) longlabel ="Akt10 JES7 (GJ:kterm)    "; 
  else if(sshort=="BoostedJES8" ) longlabel ="Akt10 JES8 (GJ:JER)      "; 
  else if(sshort=="BoostedJES9" ) longlabel ="Akt10 JES9 (GJ:akt4in/out)"; 
  else if(sshort=="BoostedJES10") longlabel ="Akt10 JES10 (GJ:more1smallJet)"; 
  else if(sshort=="BoostedJES11") longlabel ="Akt10 JES11 (GJ:stats)   "; 
  else if(sshort=="BoostedJES12") longlabel ="Akt10 JES12 (GJ:MoverPt) "; 
  else if(sshort=="BoostedJES13") longlabel ="Akt10 JES (Topology)"; 
  else if(sshort=="BoostedJES14") longlabel ="Akt10 JES (Extrapl.)"; 
  else if(sshort=="BoostedJES15") longlabel ="Akt10 JES (NPV)     "; 
  else if(sshort=="BoostedJES16") longlabel ="Akt10 JES ($\\mu$)   "; 
  else if(sshort=="BoostedJMS"  ) longlabel ="Akt10 jet m/d12 scale "; 
  else if(sshort=="BoostedJER"  ) longlabel ="Akt10 jet energy res. "; 
  else if(sshort=="BoostedJMR"  ) longlabel ="Akt10 jet mass res.   "; 


  else if(sshort=="JES3"    )     longlabel ="Akt4 JES (Modelling1)    ";              
  else if(sshort=="JES7"    )     longlabel ="Akt4 JES (Detector1)     ";
  else if(sshort=="JES12"   )     longlabel ="Akt4 JES (EtaInter\_Model)";
  else if(sshort=="JES18"   )     longlabel ="Akt4 JES (PileupRhoTop)  ";
  else if(sshort=="JES20"   )     longlabel ="Akt4 JES (FlavComp)      ";
  else if(sshort=="JES21"   )     longlabel ="Akt4 JES (FlavResp)      ";
  else if(sshort=="JES22"   )     longlabel ="Akt4 JES (BJES)          ";
  else if(sshort=="JESsmall")     longlabel ="Akt4 JES (Others)        ";

  else longlabel=sshort;

  //cout<<sshort<<"\t=> "<<longlabel<<endl;
  return longlabel;

}

