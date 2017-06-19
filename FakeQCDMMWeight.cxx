#include "FakeQCDMMWeight.h"

#include "TH1.h"
#include "TFile.h"
#include <string>
#include <algorithm>

#include <iostream>

FakeQCDMMWeight::FakeQCDMMWeight(const std::string &filename, const std::string &filename_pre) {
  m_file = new TFile(filename.c_str(), "READ");
  if (filename_pre == "")
    m_file_pre = 0;
  else
    m_file_pre = new TFile(filename_pre.c_str(), "READ");
}

double FakeQCDMMWeight::restrictRange(double a, TAxis *axis) {
  if (!axis) return a;
  double amax = axis->GetXmax()-0.1;
  double amin = axis->GetXmin()+0.1;
  if (a > amax) a = amax;
  if (a < amin) a = amin;
  return a;
}

double FakeQCDMMWeight::eff(bool electron, float *x, int nx, bool nom, bool bup) {

  // many histograms might be used, select the right one
  TH1 *h_eff = getEffHist(electron, x, nx);

  // some parameters might be used to select the histogram above, so select the parameters used in the histogram itself
  float a = 0;
  float b = 0;
  float c = 0;
  int n = 0; // number of bins
  getEffBin(electron, x, nx, a, b, c, n);

  int ibin = 0;
  if (n == 1) {
    a = restrictRange(a, h_eff->GetXaxis());
    ibin = h_eff->FindFixBin(a);
  } else if (n == 2) {
    a = restrictRange(a, h_eff->GetXaxis());
    b = restrictRange(b, h_eff->GetYaxis());
    ibin = h_eff->FindFixBin(a, b);
  } else if (n == 3) {
    a = restrictRange(a, h_eff->GetXaxis());
    b = restrictRange(b, h_eff->GetYaxis());
    c = restrictRange(c, h_eff->GetZaxis());
    ibin = h_eff->FindFixBin(a, b, c);
  }
  double eff = h_eff->GetBinContent(ibin);
  if (!nom) {
    if (bup)
      eff += h_eff->GetBinError(ibin);
    else
      eff -= h_eff->GetBinError(ibin);
  }

  return eff;
}

double FakeQCDMMWeight::fake(bool electron, float *x, int nx, bool nom, bool bup) {

  // many histograms might be used, select the right one
  TH1 *h_eff = getFakeHist(electron, x, nx);

  // some parameters might be used to select the histogram above, so select the parameters used in the histogram itself
  float a = 0;
  float b = 0;
  float c = 0;
  int n = 0; // number of bins
  getFakeBin(electron, x, nx, a, b, c, n);


  int ibin = 0;
  if (n == 1) {
    a = restrictRange(a, h_eff->GetXaxis());
    ibin = h_eff->FindFixBin(a);
  } else if (n == 2) {
    a = restrictRange(a, h_eff->GetXaxis());
    b = restrictRange(b, h_eff->GetYaxis());
    ibin = h_eff->FindFixBin(a, b);
  } else if (n == 3) {
    a = restrictRange(a, h_eff->GetXaxis());
    b = restrictRange(b, h_eff->GetYaxis());
    c = restrictRange(c, h_eff->GetZaxis());
    ibin = h_eff->FindFixBin(a, b, c);
  }

  double eff = h_eff->GetBinContent(ibin);
  if (!nom) {
    if (bup)
      eff += h_eff->GetBinError(ibin);
    else
      eff -= h_eff->GetBinError(ibin);
  }

  return eff;
}

FakeQCDMMWeight::~FakeQCDMMWeight() {
  delete m_file;
  if (m_file_pre != 0) delete m_file_pre;
}

double FakeQCDMMWeight::weight(bool electron, bool passTight, float *x, float *y, int nx, int ny, bool nom, bool eff_error, bool bup) {
  double effic = 1;
  double fkr = 0;

  if (nom) {
    effic = eff(electron, x, nx, true);
    fkr = fake(electron, y, ny, true);
  } else {
    if (eff_error) {
      effic = eff(electron, x, nx, false, bup);
      fkr = fake(electron, y, ny, true);
    } else {
      effic = eff(electron, x, nx, true);
      fkr = fake(electron, y, ny, false, bup);
    }
  }

  double wA = effic*fkr/(effic - fkr);
  double wT = (effic - 1)*fkr/(effic - fkr);

  if (wA != wA) wA = 0;
  if (wT != wT) wT = 0;

  if (passTight)
    return wT;

  return wA;
}

