#ifndef ROOTCORE_PACKAGES_H
#define ROOTCORE_PACKAGES_H

// This file contains one define statement for each package detected by
// RootCore.  It is meant to allow to write package A in a form that
// it uses the services of package B when available, but otherwise
// skips it without failing.  For this to work properly you need to list
// package B in the PACKAGE_TRYDEP of package A.

#define ROOTCORE_PACKAGE_RootCore
#define ROOTCORE_PACKAGE_ApplyJetCalibration
#define ROOTCORE_PACKAGE_JetResolution
#define ROOTCORE_PACKAGE_ApplyJetResolutionSmearing
#define ROOTCORE_PACKAGE_CalibrationDataInterface
#define ROOTCORE_PACKAGE_PATCore
#define ROOTCORE_PACKAGE_egammaAnalysisUtils
#define ROOTCORE_PACKAGE_ElectronEfficiencyCorrection
#define ROOTCORE_PACKAGE_GoodRunsLists
#define ROOTCORE_PACKAGE_JVFUncertaintyTool
#define ROOTCORE_PACKAGE_JetEffiProvider
#define ROOTCORE_PACKAGE_JetUncertainties
#define ROOTCORE_PACKAGE_MuonMomentumCorrections
#define ROOTCORE_PACKAGE_MissingETUtility
#define ROOTCORE_PACKAGE_MultiJesInputFiles
#define ROOTCORE_PACKAGE_MuonEfficiencyCorrections
#define ROOTCORE_PACKAGE_PileupReweighting
#define ROOTCORE_PACKAGE_TopDataPreparation
#define ROOTCORE_PACKAGE_TopElectronSFUtils
#define ROOTCORE_PACKAGE_TopGoodRunsList
#define ROOTCORE_PACKAGE_TopJetUtils
#define ROOTCORE_PACKAGE_TopMuonSFUtils
#define ROOTCORE_PACKAGE_TrigMuonEfficiency
#define ROOTCORE_PACKAGE_WjetsCorrections
#define ROOTCORE_PACKAGE_egammaEvent
#define ROOTCORE_PACKAGE_ttResoSingleLepton

#endif