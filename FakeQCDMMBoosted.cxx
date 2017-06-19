#include "FakeQCDMMBoosted.h"
#include "FakeQCDMMWeight.h"
#include <cmath>

FakeQCDMMBoosted::FakeQCDMMBoosted(const std::string &filename, const std::string &filename_pre)
  : FakeQCDMMWeight(filename, filename_pre) {
}

FakeQCDMMBoosted::~FakeQCDMMBoosted() {
}

TH1 *FakeQCDMMBoosted::getEffHist(bool electron, float *x, int nx) {
  // parametrised as a function of (lepton p_T [GeV])
  if ( electron && (nx != 3)) return 0;
  if (!electron && (nx != 3)) return 0;
  bool btagging = (x[0] > 0);
  if (electron) {
    float dR = x[2];
    if (btagging)
      return (TH1 *) m_file->Get("eff_boo_lepPtdREl");
    // else
    return (TH1 *) m_file_pre->Get("eff_boo_lepPtdREl");
  } else {
    float dR = x[2];
    if (btagging)
      return (TH1 *) m_file->Get("eff_boo_lepPtdRMu");
    // else
    return (TH1 *) m_file_pre->Get("eff_boo_lepPtdRMu");
  }

  return 0; // should not happen
}

TH1 *FakeQCDMMBoosted::getFakeHist(bool electron, float *y, int ny) {
  // parametrised as a function of (lepton p_T [GeV])
  if ( electron && (ny != 5)) return 0;
  if (!electron && (ny != 5)) return 0;
  float sd0 = std::fabs(y[4]);
  bool btagging = (y[0] > 0);
  float dR = y[2];
  if (electron) {
    if (btagging) {
      if (dR <= 0.4) {
        return (TH1 *) m_file->Get("fkr_boo_lepPtJetPtdRl04El");
        //if (sd0 < 3) {
        //  return (TH1 *) m_file->Get("fkr_boo_lepPtJetPtdRl04lsd0El");
        //} else {
        //  return (TH1 *) m_file->Get("fkr_boo_lepPtJetPtdRl04hsd0El");
        //}
      } else {
        return (TH1 *) m_file->Get("fkr_boo_lepPtdR04El");
        //if (sd0 < 3) {
        //  return (TH1 *) m_file->Get("fkr_boo_lepPtdR04lsd0El");
        //  //return (TH1 *) m_file->Get("fkr_boo_lepPtJetPtdRg04lsd0El");
        //} else {
        //  return (TH1 *) m_file->Get("fkr_boo_lepPtdR04hsd0El");
        //  //return (TH1 *) m_file->Get("fkr_boo_lepPtJetPtdRg04hsd0El");
        //}
      }
    } else {
      if (dR <= 0.4) {
        return (TH1 *) m_file_pre->Get("fkr_boo_lepPtJetPtdRl04El");
      } else {
        return (TH1 *) m_file_pre->Get("fkr_boo_lepPtdR04El");
      }
    }
  } else {
    if (btagging) {
      if (dR <= 0.4) {
        return (TH1 *) m_file->Get("fkr_boo_lepPtJetPtdRl04Mu");
        //if (sd0 < 3) {
        //  return (TH1 *) m_file->Get("fkr_boo_lepPtJetPtdRl04lsd0Mu");
        //} else {
        //  return (TH1 *) m_file->Get("fkr_boo_lepPtJetPtdRl04hsd0Mu");
        //}
      } else {
        return (TH1 *) m_file->Get("fkr_boo_lepPtdR04Mu");
        //return (TH1 *) m_file->Get("fkr_boo_lepPtJetPtdRg04Mu");
        //if (sd0 < 3) {
        //  return (TH1 *) m_file->Get("fkr_boo_lepPtJetPtdRg04lsd0Mu");
        //} else {
        //  return (TH1 *) m_file->Get("fkr_boo_lepPtJetPtdRg04hsd0Mu");
        //}
      }
    } else {
      if (dR <= 0.4) {
        return (TH1 *) m_file_pre->Get("fkr_boo_lepPtJetPtdRl04Mu");
      } else {
        return (TH1 *) m_file_pre->Get("fkr_boo_lepPtdR04Mu");
      }
    }
  }
  return 0; // should not happen
}

// expect x = {btagging, lepton p_T [GeV]}
void FakeQCDMMBoosted::getEffBin(bool electron, float *x, int nx, float &a, float &b, float &c, int &n) {
  if ( electron && (nx != 3)) return;
  if (!electron && (nx != 3)) return;
  
  if (electron) {
    float dR = x[2];
    n = 2;
    a = x[1];
    b = x[2];
  } else {
    float dR = x[2];
    n = 2;
    a = x[1];
    b = x[2];
  }
}

// expect y = {btagging, lepton p_T [GeV]}
void FakeQCDMMBoosted::getFakeBin(bool electron, float *y, int ny, float &a, float &b, float &c, int &n) {
  if ( electron && (ny != 5)) return;
  if (!electron && (ny != 5)) return;
  
  if (electron) {
    float dR = y[2];
    if (dR > 0.4) {
      //n = 2;
      //a = y[3]; // close jet pT
      //b = y[1]; // lepton pT
      n = 1;
      a = y[1];
    } else {
      n = 2;
      a = y[3]; // close jet pT
      b = y[1]; // lepton pT
    }
  } else {
    float dR = y[2];
    if (dR > 0.4) {
      //n = 2;
      //a = y[3]; // close jet pT
      //b = y[1]; // close jet pT
      n = 1;
      a = y[1];
    } else {
      n = 2;
      a = y[3]; // close jet pT
      b = y[1]; // lepton pT
    }
  }
}

double FakeQCDMMBoosted::getWeight(bool electron, bool passTight, bool btagging,
                                   float lepPt, float lepJetdR, float closeJetPt, float sd0,
                                   bool nom, bool eff_error, bool bup) {
  float f_btagging = 1.0;
  if (!btagging) f_btagging = -1.0;

  double w = 1;

  if (electron) {
    float x[3] = {f_btagging, lepPt, lepJetdR};
    float y[5] = {f_btagging, lepPt, lepJetdR, closeJetPt, sd0};

    float *mx = x;
    float *my = y;

    w = weight(electron, passTight, mx, my, 3, 5, nom, eff_error, bup);
  } else {
    float x[3] = {f_btagging, lepPt, lepJetdR};
    float y[5] = {f_btagging, lepPt, lepJetdR, closeJetPt, sd0};

    float *mx = x;
    float *my = y;

    w = weight(electron, passTight, mx, my, 3, 5, nom, eff_error, bup);
  }

  return w;
}

double FakeQCDMMBoosted::getEff(bool electron, bool btagging, float lepPt, float lepJetdR, bool nom, bool bup) {
  float f_btagging = 1.0;
  if (!btagging) f_btagging = -1.0;

  if (electron) {
    float x[3] = {f_btagging, lepPt, lepJetdR};
    float *mx = x;
    return eff(electron, mx, 3, nom, bup);
  } else {
    float x[3] = {f_btagging, lepPt, lepJetdR};
    float *mx = x;
    return eff(electron, mx, 3, nom, bup);
  }
  return 1;
}

double FakeQCDMMBoosted::getFake(bool electron, bool btagging, float lepPt, float lepJetdR, float closeJetPt, float sd0, bool nom, bool bup) {
  float f_btagging = 1.0;
  if (!btagging) f_btagging = -1.0;

  if (electron) {
    float y[5] = {f_btagging, lepPt, lepJetdR, closeJetPt, sd0};
    float *my = y;
    return fake(electron, my, 5, nom, bup);
  } else {
    float y[5] = {f_btagging, lepPt, lepJetdR, closeJetPt, sd0};
    float *my = y;
    return fake(electron, my, 5, nom, bup);
  }
}

