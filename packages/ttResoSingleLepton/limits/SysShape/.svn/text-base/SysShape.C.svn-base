#include <fstream>
#include <map>
using namespace std;
void SysShape(TString Chan1, TString cond,TString sys)
{
  gStyle->SetOptStat(0);
  gStyle->SetLineWidth(2);                     // width of ticks
  gStyle->SetPadTickX(1); //0//1:ticks on upper,2: ticks+labels on upper xaxis
  gStyle->SetPadTickY(1);
  double lsize=0.01;
  double tsize=0.03;
  gStyle->SetTextSize(tsize);
  gStyle->SetLabelSize(lsize,"x");
  gStyle->SetLabelSize(lsize,"y");


  //andrew
  gStyle->SetOptTitle(0);
  gStyle->SetStatBorderSize(0); 
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  
  
  gStyle->SetLineWidth(3); // width of ticks
  gStyle->SetPadTickX(1); //0//1:ticks on upper,2: ticks+labels on upper xaxis
  gStyle->SetPadTickY(1); //0

  gStyle->SetPadLeftMargin(0.16); // 0.21// 0.18
  gStyle->SetPadRightMargin(0.05); //0.08
  gStyle->SetPadTopMargin(0.05); //0.07
  gStyle->SetPadBottomMargin(0.17);
  gStyle->SetFrameLineWidth(3);
  gStyle->SetFrameFillColor(0);
  gStyle->SetPaintTextFormat(".2f");

  
  //TFile* input = TFile::Open("07_04_2013_glgw_Systematics.root");//andrew
  TFile* input = TFile::Open("glgw_Systematics.root");//andrew
  TString HistName1("massTTbarChi2LPC");
  //TString HistName1("masstT");
  //   TString Chan1("e");
  // TString Chan1("mu");
  //TString cond("resolvedChi2");
  //TString cond("boosted");
  if (cond=="boosted") HistName1="masstT";
  if (cond=="resolvedChi2") HistName1="massTTbarChi2LPC";

  if (cond=="boosted_cat1") HistName1="masstT_cat1";
  if (cond=="boosted_cat2") HistName1="masstT_cat2";
  if (cond=="boosted_cat3") HistName1="masstT_cat3";
  if (cond=="resolved_cat1") HistName1="massTTbarChi2LPC_cat1";
  if (cond=="resolved_cat2") HistName1="massTTbarChi2LPC_cat2";
  if (cond=="resolved_cat3") HistName1="massTTbarChi2LPC_cat3";

  //TString mass("Z1750");
  //TString mass("KKg2000");
  TString mass("may27");
  //convert nan bin error in this input
  TIter next_key(input->GetListOfKeys());
  TKey* key(0);
  while ((key = dynamic_cast<TKey*>( next_key() )) !=0 ) {
    TObject* obj = key->ReadObj();
    TH1* h = dynamic_cast<TH1*>( obj );
    string hname=h->GetName();
    if (hname.find("IFSR")==string::npos && hname.find("PartonShower")==string::npos) continue;
    for (int i=0; i<h->GetNbinsX()+2;i++)
      if (TMath::IsNaN(h->GetBinError(i))) h->SetBinError(i,1e-6);
  }
  map<string, string> sysNameMap;
  sysNameMap["JER"]="JetEnerRes";
  sysNameMap["MC_Gen"]="MCGen";
  sysNameMap["Norm_ttbar"]="norm_tt";
  sysNameMap["PS"]="PartonShower";
  sysNameMap["ctag"]="BtagC";
  sysNameMap["mistag"]="BtagL";
  sysNameMap["top_mass"]="topmass";

  sysNameMap["btag7"]="Btag7";
  sysNameMap["btag8"]="Btag8";
  sysNameMap["btag9"]="Btag9";
  sysNameMap["btag10"]="Btag10";


  //
  map<string,double> sysshift;
  ifstream fin;
  fin.open(("sysposterior."+mass+".log").Data());
  string line;
  while(getline(fin, line)) {    
    istringstream iss(line);
    string name,sysvalue;
    iss>>name>>sysvalue;
    if(sysNameMap.count(name)!=0) name=sysNameMap[name];
    if (name==string(sys) || sys=="all") {sysshift[name]=atof(sysvalue.c_str());
      cout<<name<<" "<<sysshift[name]<<endl;}
  }

  cout<<"getting histograms \n";
  map<string,double>::iterator iter; 
  TH1F* hsys;  
  cout<<"binkkg_Bgr_"+HistName1+"_"+Chan1<<endl;
  TH1F* hcentral=input->Get(("binkkg_Bgr_"+HistName1+"_"+Chan1).Data());
  cout<<"binkkg_datahisto_"+HistName1+"_"+Chan1<<endl;
  cout<<"data disabled!! \n";
  //TH1F* hdata=input->Get(("binkkg_datahisto_"+HistName1+"_"+Chan1).Data());
  hdata=hcentral;
  
  for (iter=sysshift.begin(), int nsys=0; iter!=sysshift.end(); iter++, nsys++)
    {

      TH1F* htmp;
      TString sysname(iter->first);    
      double shift=iter->second;
      cout<<nsys<<" "<<HistName1<<" "<<Chan1<<" "<<sysname<<" "<<shift<<endl;
      cout<<"binkkg_Bgr_"+HistName1+"_"+Chan1+"_"+sysname+"_up"<<endl;
      if( shift>0) {
	if (nsys==0) { hsys=(TH1F*)input->Get(("binkkg_Bgr_"+HistName1+"_"+Chan1+"_"+sysname+"_up").Data()); hsys->Add(hcentral,-1); hsys->Scale(fabs(shift)); hsys->Add(hcentral);} 
	else { htmp=(TH1F*)input->Get(("binkkg_Bgr_"+HistName1+"_"+Chan1+"_"+sysname+"_up").Data()); htmp->Add(hcentral,-1); hsys->Add(htmp,fabs(shift));}
      }
      else {
	if (nsys==0) { hsys=(TH1F*)input->Get(("binkkg_Bgr_"+HistName1+"_"+Chan1+"_"+sysname+"_dw").Data()); hsys->Add(hcentral,-1); hsys->Scale(fabs(shift));  hsys->Add(hcentral);} 
	else { htmp=(TH1F*)input->Get(("binkkg_Bgr_"+HistName1+"_"+Chan1+"_"+sysname+"_dw").Data()); 
	  htmp->Add(hcentral,-1);
	  hsys->Add(htmp,fabs(shift));}
      }
      delete htmp;
      //  if (nsys==13)  { hsys->Draw();break;} 
    }
  cout<<"get "<<nsys<<" systematics"<<endl;
  for (int i=0;i<hsys->GetNbinsX();i++) hsys->SetBinError(i+1,hcentral->GetBinError(i+1));
  hsys->SetLineColor(2);
  hsys->SetLineWidth(2);
  TH1F* ratio= hdata->Clone();
  TH1F* ratio_2= hdata->Clone();
  // ratio->Divide(hsys,hcentral,1.,1.);
  ratio->Divide(hdata,hsys,1.,1.);
  ratio_2->Divide(hdata,hcentral,1.,1);
  for (int i=1; i<=ratio->GetNbinsX();i++) if(hdata->GetBinContent(i)!=0) {
    ratio->SetBinError(i,hdata->GetBinError(i)/hsys->GetBinContent(i));
    ratio_2->SetBinError(i,hdata->GetBinError(i)/hcentral->GetBinContent(i));
  }
  ratio->SetMaximum(1.4);
  ratio->SetMinimum(0.6);
  // ratio->SetMaximum(1.05);
  // ratio->SetMinimum(0.95);
  TCanvas* c1=new TCanvas("mycanvas","",700,700);
  TPad *pad1 = new TPad("pad1","pad1",0,0.23,1.,1); //upper pad 
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.23);
  //  pad1->SetRightMargin(0.05);
  pad1->SetTopMargin(0.05);
  pad1->SetBottomMargin(0.05); 
  pad1->SetLogy(1);
  pad2->SetTopMargin(0.05);
  //  pad2->SetRightMargin(0.05);
  pad2->SetBottomMargin(0.3);
  pad2->SetFillColor(0);
  c1->Clear();
  pad1->Draw();
  pad1->cd();
  // if (HistName1=="massTTbarChi2LPC") hsys->SetTitle(("resolved "+Chan1).Data());
  // else hsys->SetTitle(("boosted "+Chan1).Data());
  hdata->SetTitle("");
  //hdata->Draw(); //andrew
  //hsys->Draw("same"); //andrew
  hsys->GetYaxis()->SetTitle("Events");
  hsys->Draw();
  hcentral->SetLineWidth(2);
  hcentral->SetLineColor(4);
  hcentral->SetLineStyle(2);
  hcentral->SetMarkerStyle(5);
  hcentral->SetMarkerColor(4);
  hdata->SetMarkerStyle(20);
  hdata->SetLineWidth(2);
  hdata->SetLineColor(1);  
  hcentral->GetYaxis()->SetTitle("Events");
  hcentral->Draw("same");
  //  for (int i = 0; i<=hsys->GetNbinsX(); i++) cout<<hsys->GetBinContent(i)<<endl;
  double lxlow=0.5, lylow= 0.7, lxup=0.92, lyup=0.92;
  TLegend* leg = new TLegend(lxlow, lylow, lxup, lyup);
  leg->SetNColumns(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(tsize);
  leg->SetHeader((cond+"_"+Chan1+" for "+mass).Data());
  leg->AddEntry(hdata, "data");
  if (sys=="all") leg->AddEntry(hsys,"bkg after systematic variations");
  else leg->AddEntry(hsys,("bkg after "+sys+" variation").Data());
  leg->AddEntry(hcentral,"nominal bkg");
  //leg->AddEntry(hcentral, "Pseudo-data (Bgr central)");
  
  leg->Draw();
  c1->cd();
  pad2->SetTicks(1,1);
  pad2->SetGrid(1,1);
  pad2->Draw();
  pad2->cd();
  ratio->SetTitle();
  // ratio->GetYaxis()->SetTitle("sys/central");
  ratio->GetYaxis()->SetTitle("data/bkg");
  ratio->SetNdivisions(505,"y");
  ratio->GetYaxis()->SetLabelSize(0.1);
  ratio->GetXaxis()->SetLabelSize(0.1);
  ratio->GetYaxis()->SetTitleOffset(0.35);
  ratio->GetYaxis()->SetTitleSize(0.1);
  ratio->GetXaxis()->SetTitleSize(0.1);
  ratio->SetLineColor(2);
  ratio->SetLineWidth(2);
  ratio_2->SetLineWidth(2);
  ratio_2->SetLineStyle(2);
  ratio_2->SetLineColor(4);
  //  ratio_2->SetMarkerStyle(5);
  //ratio_2->SetMarkerColor(4);
  ratio->Draw();
  ratio_2->Draw("same");
  
  c1->SaveAs(("sysShape/obs_sys_"+mass+"_"+cond+"_"+Chan1+"_"+sys+".pdf").Data());
  c1->SaveAs(("sysShape/obs_sys_"+mass+"_"+cond+"_"+Chan1+"_"+sys+".eps").Data());
  c1->SaveAs(("sysShape/obs_sys_"+mass+"_"+cond+"_"+Chan1+"_"+sys+".png").Data());
  
}
