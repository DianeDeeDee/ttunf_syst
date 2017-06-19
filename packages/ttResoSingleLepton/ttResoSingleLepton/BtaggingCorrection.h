#include <utility>

class BtaggingCorrection{
  
  enum display{NONE=0, INFO, DEBUG, VERBOSE};
  BtaggingCorrection::display m_display;
  
 public:
  BtaggingCorrection(BtaggingCorrection::display t_display=NONE);
  ~BtaggingCorrection();
  std::pair<double, double> GetBtaggingCorrection_Eigenvector_GeV(double jet_pT, double jet_eta, double jet_jfit_deltaR, int flavor, bool isBtagged, double BtagDataEffiency, double BtagMCEffiency, double SFin_err_up /*1+delta*/, double SFin_err_dw /*1-delta*/);
  std::pair<double, double> GetBtaggingCorrection_Eigenvector_MeV(double jet_pT, double jet_eta, double jet_jfit_deltaR, int flavor, bool isBtagged, double BtagDataEffiency, double BtagMCEffiency, double SFin_err_up /*1+delta*/, double SFin_err_dw /*1-delta*/);

  std::pair<double, double> GetBtaggingCorrection_GeV(double jet_pT, double jet_eta, double jet_jfit_deltaR, int flavor, bool isBtagged, double SFeff, double SFeff_err /*delta*/, double BtagDataEffiency, double BtagMCEffiency);
  std::pair<double, double> GetBtaggingCorrection_MeV(double jet_pT, double jet_eta, double jet_jfit_deltaR, int flavor, bool isBtagged, double SFeff, double SFeff_err /*delta*/, double BtagDataEffiency, double BtagMCEffiency);
  
  inline void SetDisplay(int i){m_display = (BtaggingCorrection::display)i;} // a bit dangerous, but allow to quickly change the level of message, while debugging
};

