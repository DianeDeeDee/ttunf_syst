#ifndef JETEFFICIENCYESTIMATOR_H
#define JETEFFICIENCYESTIMATOR_H

#include <vector>
#include <map>
#include <TROOT.h>
#include <TRandom3.h>

class JetEfficiencyEstimator{
public:
  JetEfficiencyEstimator(unsigned int seed = 8364323);
  ~JetEfficiencyEstimator(){}
  inline void useICHEPResult(){m_useICHEPResult = true;}
  bool isGoodJet(double, double);
  bool checkEvent(std::vector<double>, std::vector<double>, int&);
  bool checkEvent(std::vector<double>, std::vector<double>, int&, std::vector<bool>&);
  void setSeed(unsigned int seed) { m_rndmGen.SetSeed(seed); }

 private:
  TRandom3 m_rndmGen;
  std::map<double, double> m_effiMap; // efficiency
  std::map<double, double> m_effiMapLow; //  low
  std::map<double, double> m_effiMapHigh; // high
  bool m_useICHEPResult;
};

#endif
