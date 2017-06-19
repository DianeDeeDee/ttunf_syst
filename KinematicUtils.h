#ifndef KINEMATICUTILS_H
 #define KINEMATICUTILS_H
 
 namespace KinematicUtils {
 
  double deltaR(double eta1, double eta2, double phi1, double phi2);
 
  double deltaPhi(double phi1, double phi2);
 
  double transMass(double ptLep, double phiLep, double met, double phiMet);
 
  float mcEnergy(float pt, float eta, float m);
 }
 
 #endif 
