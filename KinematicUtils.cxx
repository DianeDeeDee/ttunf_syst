#include "KinematicUtils.h"
#include <cmath>

namespace KinematicUtils {

  double deltaR(double eta1, double eta2, double phi1, double phi2) {
    double dPhi=deltaPhi(phi1,phi2);
    double dEta=std::fabs(eta1-eta2);
    double dR=std::sqrt(std::pow(dEta,2)+std::pow(dPhi,2));
    return dR;
  }
//----------------------------------------------------------------------------------------------------------------------


  double deltaPhi(double phi1, double phi2) {
    double dPhi=std::fabs(phi1-phi2);
    if (dPhi>M_PI) dPhi=2*M_PI-dPhi;
    return dPhi;
  }

//----------------------------------------------------------------------------------------------------------------------
  double transMass(double ptLep, double phiLep, double met, double phiMet) {
//     return std::sqrt(2.0*ptLep*met*( 1 - std::cos( phiLep-phiMet ) ) );
//
    return std::sqrt(2.0*ptLep*met*( 1 - std::cos( KinematicUtils::deltaPhi(phiLep,phiMet) ) ) );
  }

}










