#include "AllWeightCalcs.h"
#include "Event.h"

#include "ttResoSingleLepton/BtaggingCorrection.h"

#include "egammaAnalysisUtils/VertexPositionReweightingTool.h"
#include "TrigMuonEfficiency/LeptonTriggerSF.h"
#include "MuonEfficiencyCorrections/AnalysisMuonConfigurableScaleFactors.h"
#include "CalibrationDataInterface/CalibrationDataInterfaceROOT.h"
#include "TopElectronSFUtils/electron_MI_SF_R17.h"

#include <utility>
#include "ElectronEfficiencyCorrection/TElectronEfficiencyCorrectionTool.h"
#include "ttResoSingleLepton/TTBarElectronSFTool.h"

std::vector<double> WeightCalculator::JVFSF::calc(Event &e) {
  std::vector<double> x;
  x.push_back(1);
  return x;
}

std::vector<double> WeightCalculator::MC::calc(Event &e) {
  std::vector<double> x;
  x.push_back(e.mcWeight());
  return x;
}

std::vector<double> WeightCalculator::Pileup::calc(Event &e) {
  std::vector<double> x;
  x.push_back(1);
  if (tools->m_isData) return x;
  float m = (!tools->m_isData && e.lbn() == 1 && int(e.mu()+0.5)==1) ? 0. : e.mu();
  int rn = e.runNumber();
  if (rn == 212399)
    rn = 195848;
  x[0] = tools->m_pileupTool.GetCombinedWeight(rn, e.channelNumber(), m);
  return x;
}

std::vector<double> WeightCalculator::ElectronSF::calc(Event &e) {
  std::vector<double> x;
  x.push_back(1);
  if (tools->m_isData) return x;
  x.push_back(1);
  x.push_back(1);


  for (int k = 0; k < e.electron().size(); ++k) {
    double sf = 1;

    double minDr = 99;
    double minZ = -1;
    for (int z = 0; z < e.jet().size(); ++z) {
      if (e.jet()[k].pass()) {
        double dr = e.jet()[z].mom().DeltaR(e.electron()[k].mom());
        if (dr < minDr) {
          minZ = z;
          minDr = dr;
        }
      }
    }
    std:pair<double, double> sf_err = tools->m_ettbarsfTool.GetSF(e.runNumber(), e.electron()[k].mom().Perp(), e.electron()[k].caloMom().Eta(),
        minDr, e.jet()[minZ].mom().Perp(), false);
    x[0] *= sf_err.first;
    x[1] *= sf_err.first*(1 + sf_err.second);
    x[2] *= sf_err.first*(1 - sf_err.second);
  }
  /*
  for (int k = 0; k < e.electron().size(); ++k) {
    double sf = 1;
    double sferr = 1;
    sf *= ele_ID_SF(e.electron()[k].mom().Eta(), e.electron()[k].mom().Perp());
    sferr = std::pow(ele_ID_SF_err(e.electron()[k].mom().Eta(), e.electron()[k].mom().Perp()), 2);

    sf *= ele_reco_SF(e.electron()[k].mom().Eta());
    sferr += std::pow(ele_reco_SF_err(e.electron()[k].mom().Eta()), 2);

    int set = 0;
    if (e.period() == Event::perAtoB3) set = 0;
    else if ( (e.period() == Event::perB4toB14) || (e.period() == Event::perC6toD3) ) set = 1;
    else if (e.period() == Event::perC1toC5) set = 2;
    else if ( (e.period() == Event::perD4toX) || (e.period() == Event::perXon) ) set = 3;
    sf *= ele_trigger_SF(e.electron()[k].mom().Eta(), e.electron()[k].mom().Perp(), set);
    sferr += std::pow(ele_trigger_SF_err(e.electron()[k].mom().Eta(), e.electron()[k].mom().Perp(), set), 2);

    sferr = std::sqrt(sferr);
    x[0] *= sf;
    x[1] *= sf*(1 + sferr);
    x[2] *= sf*(1 - sferr);
  }
  */
  return x;
}

std::vector<double> WeightCalculator::MuonSF::calc(Event &e) {
  std::vector<double> x;
  x.push_back(1);
  if (tools->m_isData) return x;
  x.push_back(1);
  x.push_back(1);
  for (int k = 0; k < e.muon().size(); ++k) {
    double sf = 1;
    double sferr = 0;

    std::vector<TLorentzVector> muonVec;
    muonVec.push_back(e.muon()[k].mom());
    std::vector<TLorentzVector> eleVec;
    std::pair<double, double> tSF = tools->m_leptonTriggerSF.GetTriggerSF(tools->m_pileupTool.GetRandomRunNumber((int) e.period()), false, muonVec, combined, eleVec, tightpp, 0);

    sf *= tSF.first;
    sferr += std::pow(tSF.second/tSF.first, 2);

    sf *= 1; // ID SF
    sferr += std::pow(0.02, 2); // 2% error

    if (std::fabs(e.muon()[k].mom().Perp()) > 1e10) {
      std::cout << "DANILO: Crazy muon smearing bug again! Event "<<e.runNumber()<<" muon "<<k<<": pt " << e.muon()[k].mom().Perp() << std::endl;
      std::exit(1);
    }

    double rSF = tools->m_MCPsf.scaleFactor(e.muon()[k].charge(), e.muon()[k].mom());
    double rSF_err = 0;
    if (rSF != 0) {
      rSF_err = (tools->m_MCPsf.scaleFactorUncertainty(e.muon()[k].charge(), e.muon()[k].mom()) + tools->m_MCPsf.scaleFactorSystematicUncertainty(e.muon()[k].charge(), e.muon()[k].mom()))/rSF;
    }
    sf *= rSF;
    sferr += std::pow(rSF_err, 2);

    sferr = std::sqrt(sferr);
    sferr *= sf;

    x[0] *= sf;
    x[1] *= (sf+sferr);
    x[2] *= (sf-sferr);
  }
  return x;
}

std::vector<double> WeightCalculator::BtagSF::calc(Event &e) {
  std::vector<double> x;
  x.push_back(1);
  if (tools->m_isData) return x;

  std::string cutString = "0_7892";
  double btagSF = 1;
  double btagSFWithPlusUncertaintyForbjets = 1;
  double btagSFWithMinusUncertaintyForbjets = 1;
  double btagSFWithPlusUncertaintyForcjets = 1;
  double btagSFWithMinusUncertaintyForcjets = 1;
  double btagSFWithPlusUncertaintyForLightJets = 1;
  double btagSFWithMinusUncertaintyForLightJets = 1;

  Analysis::CalibrationDataVariables ajet;
  ajet.jetAuthor = "AntiKt4TopoLCJVF0_5";
  Analysis::Uncertainty uncertainty = Analysis::Total;
  for (int k = 0; k < e.jet().size(); ++k) {
    ajet.jetPt = e.jet()[k].mom().Perp();
    ajet.jetEta = e.jet()[k].mom().Eta();

    std::string str_flavor = "";
    if(e.jet()[k].trueFlavour()<4){      str_flavor = "Light";  }
    else if(e.jet()[k].trueFlavour()==4){   str_flavor = "C";  }
    else if(e.jet()[k].trueFlavour()==5){   str_flavor = "B";  }
    else if(e.jet()[k].trueFlavour()==15){   str_flavor = "T";  }

    Analysis::CalibResult weight;
    weight.first = 1;
    weight.second = 0;
    Analysis::CalibResult weightMC;
    weightMC.first = 1;
    weightMC.second = 0;

    bool isBtagged = e.jet()[k].btag();

    if (isBtagged) weight = tools->m_btagSFTool->getScaleFactor(ajet, str_flavor.c_str(), cutString.c_str(), uncertainty); 
    else weight = tools->m_btagSFTool->getInefficiencyScaleFactor(ajet, str_flavor.c_str(), cutString.c_str(), uncertainty); 

    double jfit_deltaR = std::sqrt(std::pow(e.jet()[k].jfit_deltaEta(),2)+std::pow(e.jet()[k].jfit_deltaPhi(),2));

    weightMC = tools->m_btagSFTool->getMCEfficiency(ajet, str_flavor.c_str(), cutString.c_str(), uncertainty);
  
    double btagMCEfficiency = weightMC.first;
    double btagDataEfficiency = btagMCEfficiency*weight.first;

    BtaggingCorrection btaggingCorrection;
    weight = btaggingCorrection.GetBtaggingCorrection_MeV(ajet.jetPt, ajet.jetEta, jfit_deltaR, e.jet()[k].trueFlavour(), isBtagged, weight.first, weight.second, btagDataEfficiency, btagMCEfficiency);

    if (isBtagged) {
      if (std::fabs(e.jet()[k].trueFlavour()) == 5) { // true bottom
        btagSFWithPlusUncertaintyForbjets *= weight.first+weight.second;
        btagSFWithMinusUncertaintyForbjets *= weight.first-weight.second;
        btagSFWithPlusUncertaintyForcjets *= weight.first;
        btagSFWithMinusUncertaintyForcjets *= weight.first;
        btagSFWithPlusUncertaintyForLightJets *= weight.first;
        btagSFWithMinusUncertaintyForLightJets *= weight.first;
      } else if (std::fabs(e.jet()[k].trueFlavour()) == 4) { // true charm
        btagSFWithPlusUncertaintyForbjets *= weight.first;
        btagSFWithMinusUncertaintyForbjets *= weight.first;
        btagSFWithPlusUncertaintyForcjets *= weight.first+weight.second;
        btagSFWithMinusUncertaintyForcjets *= weight.first-weight.second;
        btagSFWithPlusUncertaintyForLightJets *= weight.first;
        btagSFWithMinusUncertaintyForLightJets *= weight.first;
      } else if (std::fabs(e.jet()[k].trueFlavour()) == 15) { // true tau
        btagSFWithPlusUncertaintyForbjets *= weight.first;
        btagSFWithMinusUncertaintyForbjets *= weight.first;
        btagSFWithPlusUncertaintyForcjets *= weight.first+weight.second;
        btagSFWithMinusUncertaintyForcjets *= weight.first-weight.second;
        btagSFWithPlusUncertaintyForLightJets *= weight.first;
        btagSFWithMinusUncertaintyForLightJets *= weight.first;
      } else { // light jet
        btagSFWithPlusUncertaintyForbjets *= weight.first;
        btagSFWithMinusUncertaintyForbjets *= weight.first;
        btagSFWithPlusUncertaintyForcjets *= weight.first;
        btagSFWithMinusUncertaintyForcjets *= weight.first;
        btagSFWithPlusUncertaintyForLightJets *= weight.first+weight.second;
        btagSFWithMinusUncertaintyForLightJets *= weight.first-weight.second;
      }
    } else {
      if (std::fabs(e.jet()[k].trueFlavour()) == 5) { // true bottom
        btagSFWithPlusUncertaintyForbjets *= weight.first-weight.second;
        btagSFWithMinusUncertaintyForbjets *= weight.first+weight.second;
        btagSFWithPlusUncertaintyForcjets *= weight.first;
        btagSFWithMinusUncertaintyForcjets *= weight.first;
        btagSFWithPlusUncertaintyForLightJets *= weight.first;
        btagSFWithMinusUncertaintyForLightJets *= weight.first;
      } else if (std::fabs(e.jet()[k].trueFlavour()) == 4) { // true charm
        btagSFWithPlusUncertaintyForbjets *= weight.first;
        btagSFWithMinusUncertaintyForbjets *= weight.first;
        btagSFWithPlusUncertaintyForcjets *= weight.first-weight.second;
        btagSFWithMinusUncertaintyForcjets *= weight.first+weight.second;
        btagSFWithPlusUncertaintyForLightJets *= weight.first;
        btagSFWithMinusUncertaintyForLightJets *= weight.first;
      } else if (std::fabs(e.jet()[k].trueFlavour()) == 15) { // true tau
        btagSFWithPlusUncertaintyForbjets *= weight.first;
        btagSFWithMinusUncertaintyForbjets *= weight.first;
        btagSFWithPlusUncertaintyForcjets *= weight.first-weight.second;
        btagSFWithMinusUncertaintyForcjets *= weight.first+weight.second;
        btagSFWithPlusUncertaintyForLightJets *= weight.first;
        btagSFWithMinusUncertaintyForLightJets *= weight.first;
      } else { // light jet
        btagSFWithPlusUncertaintyForbjets *= weight.first;
        btagSFWithMinusUncertaintyForbjets *= weight.first;
        btagSFWithPlusUncertaintyForcjets *= weight.first;
        btagSFWithMinusUncertaintyForcjets *= weight.first;
        btagSFWithPlusUncertaintyForLightJets *= weight.first-weight.second;
        btagSFWithMinusUncertaintyForLightJets *= weight.first+weight.second;
      }
    }

    btagSF *= weight.first;
  }

  x[0] = btagSF;
  x.push_back(btagSFWithPlusUncertaintyForbjets);
  x.push_back(btagSFWithMinusUncertaintyForbjets);
  x.push_back(btagSFWithPlusUncertaintyForcjets);
  x.push_back(btagSFWithMinusUncertaintyForcjets);
  x.push_back(btagSFWithPlusUncertaintyForLightJets);
  x.push_back(btagSFWithMinusUncertaintyForLightJets);
  return x;
}

std::vector<double> WeightCalculator::ZVertexWeight::calc(Event &e) {
  std::vector<double> x;
  x.push_back(1);
  if (tools->m_isData) return x;
  x[0] = tools->m_vertexPositionReweightingTool.GetWeight(e.vxZ());
  return x;
}

