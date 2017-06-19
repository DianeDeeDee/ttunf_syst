#include<vector>
#include<TString.h>

using namespace std;

const string InputSpectraName="MttSpectra_Oct06.root";

void buildSystList(vector<string> * syst, int mode=1) { 
  // mode==1: full list
  // mode==2: no lumi, norm_tt/QCD
  // mode==3: no SFs (for smoothing)

  syst->clear();
  if (mode<2) {
    syst->push_back("luminosity");
  }

  if (mode<3) {
    syst->push_back("EleSF");
    syst->push_back("MuSF");
  }

//  syst->push_back("JES0");  
//  syst->push_back("JES1");
//  syst->push_back("JES2");  
//  syst->push_back("JES3");  
//  syst->push_back("JES4");  
//  syst->push_back("JES5");  
//  syst->push_back("JES6");  
//  syst->push_back("JES7");  
//  syst->push_back("JES8");  
//  syst->push_back("JES9");
//  syst->push_back("JES10");
//  syst->push_back("JES11");
//  syst->push_back("JES12");
//  syst->push_back("JES13");
//  syst->push_back("JES14"); 
//  syst->push_back("JES15");
//  syst->push_back("JES16");
//  syst->push_back("JES17");
//  syst->push_back("JES18");
////syst->push_back("JES19"); //zero!
//  syst->push_back("JES20");
//  syst->push_back("JES21");
//  syst->push_back("JES22");

//  syst->push_back("JES3");
//  syst->push_back("JES7");
//  syst->push_back("JES12");
//  syst->push_back("JES20");

  syst->push_back("JES3");
  syst->push_back("JES7");
  syst->push_back("JES12");
  syst->push_back("JES18");
  syst->push_back("JES20");
  syst->push_back("JES21");
  syst->push_back("JES22");

  syst->push_back("JESsmall");
  syst->push_back("JetEnerRes");
  syst->push_back("JVFCut");
  
  //  syst->push_back("BoostedJES0");
  syst->push_back("BoostedJES1");
  syst->push_back("BoostedJESothers");
//  syst->push_back("BoostedJES2");
//  syst->push_back("BoostedJES3");
//  syst->push_back("BoostedJES4");
//  syst->push_back("BoostedJES5");
//  syst->push_back("BoostedJES6");
//  syst->push_back("BoostedJES7");
//  syst->push_back("BoostedJES8");
//  syst->push_back("BoostedJES9");
//  syst->push_back("BoostedJES10");
//  syst->push_back("BoostedJES11");
//  syst->push_back("BoostedJES12");
  syst->push_back("BoostedJES13");
//  syst->push_back("BoostedJES14");
//  syst->push_back("BoostedJES15");
//  syst->push_back("BoostedJES16");
  syst->push_back("BoostedJMS");
  syst->push_back("BoostedJER");
  syst->push_back("BoostedJMR");

  if (mode<3) {
    //syst->push_back("Btag5");
    syst->push_back("Btag6");
    syst->push_back("Btag7");
    syst->push_back("Btag8");
    syst->push_back("Btag9");
    syst->push_back("Btag10");
    syst->push_back("BtagC");
    syst->push_back("BtagL");
  }

  if (mode<2) {
    syst->push_back("norm_tt");
  }
  syst->push_back("MCGen");
  syst->push_back("PartonShower");
  syst->push_back("topmass");
  syst->push_back("IFSR");
  syst->push_back("EWS");
  //syst->push_back("Pttt");
  //syst->push_back("hdamp");

  if (mode<3) {
    syst->push_back("iqopt3");
    syst->push_back("ptjmin10");
    //syst->push_back("Whfsf");
    syst->push_back("PDF");
  }

  if (mode<2) {
    syst->push_back("norm_QCDe");
    syst->push_back("norm_QCDmu");
  }

}


void buildSigList(vector<string> * signals) {
  signals->push_back("Z400");
  signals->push_back("Z500");
  signals->push_back("Z750");
  signals->push_back("Z1000");
  signals->push_back("Z1250");
  signals->push_back("Z1500");
  signals->push_back("Z1750");
  signals->push_back("Z2000");
  signals->push_back("Z2250");
  signals->push_back("Z2500");
  signals->push_back("Z3000");
  
  
  signals->push_back("KKg400");
  signals->push_back("KKg500");
  signals->push_back("KKg600");
  signals->push_back("KKg700");
  signals->push_back("KKg800");
  signals->push_back("KKg900");
  signals->push_back("KKg1000");
  signals->push_back("KKg1150");
  signals->push_back("KKg1300");
  signals->push_back("KKg1600");
  signals->push_back("KKg1800");
  signals->push_back("KKg2000");
  signals->push_back("KKg2250");
  signals->push_back("KKg2500");
  signals->push_back("KKg2750");
  signals->push_back("KKg3000");

  signals->push_back("KKg1000_width10pc");
  signals->push_back("KKg1000_width15pc");
  signals->push_back("KKg1000_width20pc");
  signals->push_back("KKg1000_width25pc");
  signals->push_back("KKg1000_width30pc");
  signals->push_back("KKg1000_width35pc");
  signals->push_back("KKg1000_width40pc");

  signals->push_back("KKg2000_width10pc");
  signals->push_back("KKg2000_width15pc");
  signals->push_back("KKg2000_width20pc");
  signals->push_back("KKg2000_width25pc");
  signals->push_back("KKg2000_width30pc");
  signals->push_back("KKg2000_width35pc");
  signals->push_back("KKg2000_width40pc");

  signals->push_back("KKg3000_width10pc");
  signals->push_back("KKg3000_width15pc");
  signals->push_back("KKg3000_width20pc");
  signals->push_back("KKg3000_width25pc");
  signals->push_back("KKg3000_width30pc");
  signals->push_back("KKg3000_width35pc");
  signals->push_back("KKg3000_width40pc");

  signals->push_back("RSG400");
  signals->push_back("RSG500");
  signals->push_back("RSG600");
  signals->push_back("RSG700");
  signals->push_back("RSG800");
  signals->push_back("RSG900");
  signals->push_back("RSG1000");
  signals->push_back("RSG1200");
  signals->push_back("RSG1400");
  signals->push_back("RSG1600");
  signals->push_back("RSG1800");
  signals->push_back("RSG2000");
  signals->push_back("RSG2500");
  
  signals->push_back("HH400");
  signals->push_back("HH500");
  signals->push_back("HH750");
  signals->push_back("HH1000");
  signals->push_back("HH1250");
  signals->push_back("HH1500");
  signals->push_back("HH1750");
  signals->push_back("HH2000");
  signals->push_back("HH2250");
  signals->push_back("HH2500");
  signals->push_back("HH2750");
  signals->push_back("HH3000");
}

void buildChList(vector<string> * channels) {
  channels->push_back("massTTbarChi2LPC_cat1_e");
  channels->push_back("massTTbarChi2LPC_cat2_e");
  channels->push_back("massTTbarChi2LPC_cat3_e");
  channels->push_back("masstT_cat1_e");
  channels->push_back("masstT_cat2_e");
  channels->push_back("masstT_cat3_e");
  channels->push_back("massTTbarChi2LPC_cat1_mu");
  channels->push_back("massTTbarChi2LPC_cat2_mu");
  channels->push_back("massTTbarChi2LPC_cat3_mu");
  channels->push_back("masstT_cat1_mu");
  channels->push_back("masstT_cat2_mu");
  channels->push_back("masstT_cat3_mu");
  channels->push_back("massTTbarChi2LPC_cat1_lep");
  channels->push_back("massTTbarChi2LPC_cat2_lep");
  channels->push_back("massTTbarChi2LPC_cat3_lep");
  channels->push_back("masstT_cat1_lep");
  channels->push_back("masstT_cat2_lep");
  channels->push_back("masstT_cat3_lep");
}
