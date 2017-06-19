
{
  gSystem->Load("libHist.so");
  gSystem->Load("../StandAlone/libTopMuonSFUtils.so");

  TRandom3 rand(3521);

  MuonSF muonsf("../data/mc11b/");

  cout << "period2011BtoD = " << TopMuon::period2011BtoD << endl;
  cout << "period2011EtoH = " << TopMuon::period2011EtoH << endl;
  cout << "period2011I = " << TopMuon::period2011I << endl;
  cout << "period2011JtoK = " << TopMuon::period2011JtoK << endl;


  cout << "period2011BtoI = " << TopMuon::period2011BtoI << endl;
  cout << "period2011LtoM = " << TopMuon::period2011LtoM << endl;

  TH2D* h_sf_eta_phi = new TH2D("h_sf_eta_phi","h_sf_eta_phi",20,-2.5,2.5,10,-3.14159,3.14159);
  TH2D* h_sf_err_eta_phi = new TH2D("h_sf_err_eta_phi","h_sf_err_eta_phi",20,-2.5,2.5,10,-3.14159,3.14159);

  for(int i=0; i<10000; ++i) {
    TLorentzVector testmu1;    
    double pt = -1.0;
    while(pt < 15.0) {
      pt = rand.Gaus(60.0,40.0);
    }
    double eta = -3.0;
    while(eta < -2.5 || eta > 2.5) {
      eta = rand.Gaus(0.0,1.0);
    }
    double phi = rand.Uniform(TMath::Pi());
    if(rand.Rndm() > 0.5) phi *= -1.0;
    testmu1.SetPtEtaPhiM( pt, eta, phi, 0.106);



    const int xbin = h_sf_eta_phi->GetXaxis()->FindFixBin( eta );
    const int ybin = h_sf_eta_phi->GetYaxis()->FindFixBin( phi );

    h_sf_eta_phi->SetBinContent(xbin, ybin, muonsf.mu_trigger_SF(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011BtoI));
    h_sf_err_eta_phi->SetBinContent(xbin, ybin, muonsf.mu_trigger_SF_err(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011BtoI));
				

    if(i < 10) {
    cout << "pt, eta, phi = " << pt << ", " << eta << ", " << phi << endl;
    cout << "\tSF to D  = " << muonsf.mu_trigger_SF(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011BtoD) << endl;
    cout << "\tSF to I = " << muonsf.mu_trigger_SF(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011BtoI ) << endl;
    cout << "\tSF in I  = " << muonsf.mu_trigger_SF(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011I ) << endl;
    cout << "\tSF Lbadrpc   = " << muonsf.mu_trigger_SF(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011L_badrpc) << endl;
    cout << "\tSF Lgoodrpc   = " << muonsf.mu_trigger_SF(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011L_goodrpc) << endl;
    cout << "\tSF J-K   = " << muonsf.mu_trigger_SF(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011JtoK) << endl;
    cout << "\tSF L-M   = " << muonsf.mu_trigger_SF(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011LtoM) << endl;
    cout << endl;
    cout << "\tSFerr to D  = " << muonsf.mu_trigger_SF_err(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011BtoD) << endl;
    cout << "\tSFerr to I = " << muonsf.mu_trigger_SF_err(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011BtoI ) << endl;
    cout << "\tSFerr in I  = " << muonsf.mu_trigger_SF_err(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011I ) << endl;
    cout << "\tSFerr Lbadrpc   = " << muonsf.mu_trigger_SF_err(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011L_badrpc) << endl;
    cout << "\tSFerr Lgoodrpc   = " << muonsf.mu_trigger_SF_err(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011L_goodrpc) << endl;
    cout << "\tSFerr J-K   = " << muonsf.mu_trigger_SF_err(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011JtoK) << endl;
    cout << "\tSFerr L-M   = " << muonsf.mu_trigger_SF_err(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011LtoM) << endl;
    cout << endl;
    cout << "\teff data to D  = " << muonsf.mu_trigger_eff_data(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011BtoD) << endl;
    cout << "\teff data to I = " << muonsf.mu_trigger_eff_data(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011BtoI ) << endl;
    cout << "\teff data in I  = " << muonsf.mu_trigger_eff_data(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011I ) << endl;
    cout << "\teff data in Lbadrpc  = " << muonsf.mu_trigger_eff_data(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011L_badrpc ) << endl;
    cout << "\teff data in Lgoodrpc  = " << muonsf.mu_trigger_eff_data(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011L_goodrpc ) << endl;
    cout << "\teff data J-K   = " << muonsf.mu_trigger_eff_data(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011JtoK) << endl;
    cout << "\teff data L-M   = " << muonsf.mu_trigger_eff_data(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011LtoM) << endl;
    cout << endl;
    cout << "\teff mc to D  = " << muonsf.mu_trigger_eff_mc(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011BtoD) << endl;
    cout << "\teff mc to I = " << muonsf.mu_trigger_eff_mc(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011BtoI ) << endl;
    cout << "\teff mc in I  = " << muonsf.mu_trigger_eff_mc(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011I ) << endl;
    cout << "\teff mc in Lbadrpc  = " << muonsf.mu_trigger_eff_mc(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011L_badrpc ) << endl;
    cout << "\teff mc in Lgoodrpc  = " << muonsf.mu_trigger_eff_mc(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011L_goodrpc ) << endl;
    cout << "\teff mc J-K   = " << muonsf.mu_trigger_eff_mc(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011JtoK) << endl;
    cout << "\teff mc L-M   = " << muonsf.mu_trigger_eff_mc(testmu1.Eta(), testmu1.Phi(),  TopMuon::period2011LtoM) << endl;
    }//print first 10
  }//muon loop
  
  cout << "Muon ID SF to I = " << muonsf.mu_ID_SF(TopMuon::period2011BtoI) << endl;
  cout << "Muon ID SF J-K  = " << muonsf.mu_ID_SF(TopMuon::period2011JtoK) << endl;
  cout << "Muon ID SF L-M  = " << muonsf.mu_ID_SF(TopMuon::period2011LtoM) << endl;

  cout << "Muon ID SF err to I = " << muonsf.mu_ID_SF_err(TopMuon::period2011BtoI) << endl;
  cout << "Muon ID SF err J-K  = " << muonsf.mu_ID_SF_err(TopMuon::period2011JtoK) << endl;
  cout << "Muon ID SF err L-M  = " << muonsf.mu_ID_SF_err(TopMuon::period2011LtoM) << endl;

  gStyle->SetPaintTextFormat("1.2f");
  TCanvas* c = new TCanvas("c","c",1500,500);
  c->Divide(2,1);
  c->cd(1);
  h_sf_eta_phi->Draw("TEXT");
  c->cd(2);
  h_sf_err_eta_phi->Draw("TEXT");

}

