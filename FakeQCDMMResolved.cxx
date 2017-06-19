#include "FakeQCDMMResolved.h"
#include "FakeQCDMMWeight.h"
#include <cmath>
#include <iostream>

FakeQCDMMResolved::FakeQCDMMResolved(const std::string &filename, const std::string &filename_pre)
  : FakeQCDMMWeight(filename, filename_pre) {
}

FakeQCDMMResolved::~FakeQCDMMResolved() {
}

TH1 *FakeQCDMMResolved::getEffHist(bool electron, float *x, int nx) {
  // parametrised as a function of (lepton p_T [GeV], delta R)
  if (nx != 3) return 0;
  bool btagging = (x[0] > 0);
  if (btagging) {
    if (electron) {
      return (TH1 *) m_file->Get("eff_res_lepPtdREl");
    } else {
      return (TH1 *) m_file->Get("eff_res_lepPtdRMu");
    }
  } else {
    if (electron) {
      return (TH1 *) m_file_pre->Get("eff_res_lepPtdREl");
    } else {
      return (TH1 *) m_file_pre->Get("eff_res_lepPtdRMu");
    }
  }
  return 0; // should not happen
}

TH1 *FakeQCDMMResolved::getFakeHist(bool electron, float *y, int ny) {
  // parametrised as a function of (lepton p_T [GeV]) for electrons
  // parametrised as a function of (lepton p_T [GeV], delta R) for muon if delta R <= 0.4 or (lepton p_T [GeV] if delta R > 0.4
  if (electron && (ny != 6)) return 0;
  if (!electron && (ny != 6)) return 0;
  bool btagging = (y[0] > 0);
  float sd0 = std::fabs(y[4]);
  float dR = y[2];
  if (electron) {
    if (btagging) {
      if (dR <= 0.4) {
        return (TH1 *) m_file->Get("fkr_res_lepPtJetPtdRl04El");
        /*
        if (std::fabs(sd0) < 1.5) {
          return (TH1 *) m_file->Get("fkr_res_lepPtJetPtdRl04asd0El");
        } else if (std::fabs(sd0) < 3) {
          return (TH1 *) m_file->Get("fkr_res_lepPtJetPtdRl04bsd0El");
        } else {
          return (TH1 *) m_file->Get("fkr_res_lepPtJetPtdRl04hsd0El");
        }*/
      } else {
        //return (TH1 *) m_file->Get("fkr_res_lepPtJetPtdRg04El");
        return (TH1 *) m_file->Get("fkr_res_lepPtdR04El");
        /*
        if (std::fabs(sd0) < 1.5) {
          return (TH1 *) m_file->Get("fkr_res_lepPtdR04asd0El");
        } else if (std::fabs(sd0) < 3) {
          //return (TH1 *) m_file->Get("fkr_res_lepPtJetPtdRg04lsd0El");
          return (TH1 *) m_file->Get("fkr_res_lepPtdR04bsd0El");
        } else {
          //return (TH1 *) m_file->Get("fkr_res_lepPtJetPtdRg04hsd0El");
          return (TH1 *) m_file->Get("fkr_res_lepPtdR04hsd0El");
        }*/
      }
    } else {
      if (dR <= 0.4) {
        return (TH1 *) m_file_pre->Get("fkr_res_lepPtJetPtdRl04El");
      } else {
        return (TH1 *) m_file_pre->Get("fkr_res_lepPtdR04El");
      }
    }
  } else { // muon
    if (btagging) {
      if (dR <= 0.4) {
        return (TH1 *) m_file->Get("fkr_res_lepPtJetPtdRl04Mu");
      } else  {
        return (TH1 *) m_file->Get("fkr_res_lepPtdR04Mu");
        //return (TH1 *) m_file->Get("fkr_res_lepPtJetPtdRg04Mu");
      }
        //return (TH1 *) m_file->Get("fkr_res_closeJetPtdR04Mu");
        //return (TH1 *) m_file->Get("fkr_res_lepPtdR04Mu");
        //return (TH1 *) m_file->Get("fkr_res_lepPtJetPtdRg04Mu");
        /*
        if (std::fabs(sd0) < 1.5) {
          return (TH1 *) m_file->Get("fkr_res_lepPtdR04asd0Mu");
        } else if (std::fabs(sd0) < 3) {
          return (TH1 *) m_file->Get("fkr_res_lepPtdR04bsd0Mu");
        } else {
          return (TH1 *) m_file->Get("fkr_res_lepPtdR04hsd0Mu");
        }
        */
        //return (TH1 *) m_file->Get("fkr_res_lepPtJetPtdRl04Mu");
        //return (TH1 *) m_file->Get("fkr_res_lepPtdRMu");
        /*
        if (std::fabs(sd0) < 1.5) {
          return (TH1 *) m_file->Get("fkr_res_lepPtJetPtdRl04asd0Mu");
        } else if (std::fabs(sd0) < 3) {
          return (TH1 *) m_file->Get("fkr_res_lepPtJetPtdRl04bsd0Mu");
        } else {
          return (TH1 *) m_file->Get("fkr_res_lepPtJetPtdRl04hsd0Mu");
        }
        */
    } else {
      if (dR <= 0.4) {
        return (TH1 *) m_file_pre->Get("fkr_res_lepPtJetPtdRl04Mu");
      } else {
        return (TH1 *) m_file_pre->Get("fkr_res_lepPtdR04Mu");
        //return (TH1 *) m_file_pre->Get("fkr_res_lepPtdR04Mu");
      }
    }
  }
  return 0; // should not happen
}
// expect x = {btagging, lepton p_T [GeV], delta R between lepton and closest jet}
void FakeQCDMMResolved::getEffBin(bool electron, float *x, int nx, float &a, float &b, float &c, int &n) {
  if (nx != 3) // ERROR!
    return;
  n = 2;
  a = x[1];
  b = x[2];
}


// expect y = {btagging, lepton p_T [GeV]} for electrons or
//        y = {btagging, lepton p_T [GeV], delta R} for muons
void FakeQCDMMResolved::getFakeBin(bool electron, float *y, int ny, float &a, float &b, float &c, int &n) {
  if (electron && (ny != 6)) return; // ERROR
  if (!electron && (ny != 6)) return; // ERROR
  float leadJetPt = y[5];
  if (electron) {
    float dR = y[2];
    if (dR > 0.4) {
      //n = 2;
      //a = y[3];
      //b = y[1];
      n = 1;
      a = y[1];
    } else {
      n = 2;
      a = y[3];
      b = y[1];
    }
  } else {
    float dR = y[2];
    if (dR > 0.4) {
      //n = 2;
      //a = y[3];
      //b = y[1];
      n = 1;
      a = y[1];
    } else {
      n = 2;
      a = y[3];
      b = y[1];
      //n = 2;
      //a = y[1];
      //b = y[2];
      //n = 1;
      //a = y[5];
    }
  }
}

double FakeQCDMMResolved::getWeight(bool electron, bool passTight, bool btagging,
                                    float lepPt, float lepJetdR, float closeJetPt, float sd0, float leadJetPt,
                                    bool nom, bool eff_error, bool bup) {
  float f_btagging = 1.0;
  if (!btagging) f_btagging = -1.0;

  double w = 1;

  if (electron) {
    float x[3] = {f_btagging, lepPt, lepJetdR};
    float y[6] = {f_btagging, lepPt, lepJetdR, closeJetPt, sd0, leadJetPt};

    float *mx = x;
    float *my = y;

    w = weight(electron, passTight, mx, my, 3, 6, nom, eff_error, bup);
  } else {
    float x[3] = {f_btagging, lepPt, lepJetdR};
    float y[6] = {f_btagging, lepPt, lepJetdR, closeJetPt, sd0, leadJetPt};

    float *mx = x;
    float *my = y;

    w = weight(electron, passTight, mx, my, 3, 6, nom, eff_error, bup);
  }
  return w;
}

double FakeQCDMMResolved::getEff(bool electron, bool btagging, float lepPt, float lepJetdR, bool nom, bool bup) {
  float f_btagging = 1.0;
  if (!btagging) f_btagging = -1.0;

  float x[3] = {f_btagging, lepPt, lepJetdR};
  float *mx = x;
  return eff(electron, mx, 3, nom, bup);
}

double FakeQCDMMResolved::getFake(bool electron, bool btagging, float lepPt, float lepJetdR, float closeJetPt, float sd0, float leadJetPt, bool nom, bool bup) {
  float f_btagging = 1.0;
  if (!btagging) f_btagging = -1.0;

  if (electron) {
    float y[6] = {f_btagging, lepPt, lepJetdR, closeJetPt, sd0, leadJetPt};
    float *my = y;
    return fake(electron, my, 6, nom, bup);
  } else {
    float y[6] = {f_btagging, lepPt, lepJetdR, closeJetPt, sd0, leadJetPt};
    float *my = y;
    return fake(electron, my, 6, nom, bup);
  }
}

