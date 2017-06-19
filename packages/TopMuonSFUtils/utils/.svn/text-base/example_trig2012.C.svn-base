
/**
 * Small example script to access 2012 trigger scale factors.
 * Needs TrigMuonEfficiency-00-02-08
 */
void example_trig2012() {
  TString dir = getenv("ROOTCOREBIN");
  gSystem->Load("libHist.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load( dir + "/lib/libegammaAnalysisUtils.so");
  gSystem->Load( dir + "/lib/libTrigMuonEfficiency.so");
  
  gInterpreter->GenerateDictionary("vector<TLorentzVector>","TLorentzVector.h;vector");
  
  dir += "/data/TopMuonSFUtils/mc12a/periodsAE/";
  LeptonTriggerSF* sftool = new LeptonTriggerSF(2012, dir, "muon_trigger_sf_muid_2012_AtoE.root");
  
  int runnumber_A = 200804; //this gives 2012 period A
  //runnumber_A = 202660; //this gives 2012 period B
  TRandom3 rand(3218);

  TH2D* h_sf_eta_phi = new TH2D("h_sf_eta_phi","h_sf_eta_phi",20,-2.5,2.5,10,-3.14159,3.14159);
  TH2D* h_sf_err_eta_phi = new TH2D("h_sf_err_eta_phi","h_sf_err_eta_phi",20,-2.5,2.5,10,-3.14159,3.14159);

  for(int i=0; i<100000; ++i) {

    TLorentzVector testmu1;    
    double pt = -1.0;
    while(pt < 25.0) {
      pt = rand.Gaus(60.0,40.0);
    }
    double eta = -3.0;
    while(eta < -2.5 || eta > 2.5) {
      eta = rand.Gaus(0.0,1.0);
    }
    double phi = rand.Uniform(TMath::Pi());
    if(rand.Rndm() > 0.5) phi *= -1.0;

    if(i==0) {
      pt = 25.0;
      eta = -0.9;
      phi = -2.9;
    }

    testmu1.SetPtEtaPhiM( pt, eta, phi, 0.106);
    


    vector<TLorentzVector> muons;
    muons.push_back(testmu1);

    pair<double, double> res = sftool->GetTriggerSF(runnumber_A, true, muons, 1);
    
    const int xbin = h_sf_eta_phi->GetXaxis()->FindFixBin( eta );
    const int ybin = h_sf_eta_phi->GetYaxis()->FindFixBin( phi );
    
    h_sf_eta_phi->SetBinContent(xbin, ybin, res.first);
    h_sf_err_eta_phi->SetBinContent(xbin, ybin, res.second);
    
    if(i<10) {
      cout << "pt, eta, phi = " << pt << ", " << eta << ", " << phi << endl;
      cout << "SF period A = " << res.first << " +- " << res.second << endl;
    }

    if(res.second < 0.0000001) {
      cout << "WARNING error close to 0:" << endl;
      cout << "pt, eta, phi = " << pt << ", " << eta << ", " << phi << endl;
      cout << "SF period A = " << res.first << " +- " << res.second << endl;
    }

  }//event loop
  gStyle->SetPaintTextFormat("1.2f");
  TCanvas* c = new TCanvas("c","c",1500,500);
  c->Divide(2,1);
  c->cd(1);
  h_sf_eta_phi->Draw("TEXT");
  c->cd(2);
  h_sf_err_eta_phi->Draw("TEXT");
}
