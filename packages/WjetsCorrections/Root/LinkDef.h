#ifdef __CINT__
#include "WjetsCorrections/ScaleWjets.h"
#include "WjetsCorrections/HFsys_factor_ttbar_emu.h"
#include "WjetsCorrections/Wasymmetry_rel17.h"
#include "WjetsCorrections/WjetsDataToMc.h"
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;

// Class from ScaleWjets
#pragma link C++ class ScaleWjets;

// Functions from HFsys_factor_ttbar_emu.h
#pragma link C++ function GetFFactors;
#pragma link C++ function SetWflavors;

// Functions from Wasymmetry_rel17.h
#pragma link C++ function GetWscaleFactorPretag;
#pragma link C++ function GetWscaleFactorTag;
#pragma link C++ function GetWstatRelPretag;
#pragma link C++ function GetWstatRelf2j;
#pragma link C++ function GetWstatRelf2toN;
#pragma link C++ function GetWweight;
#pragma link C++ function GetWweightSqErr;
#pragma link C++ function GetRvalue;
#pragma link C++ function GetWvalue;
#pragma link C++ function GetSFvalue;
#pragma link C++ function GetSFuncRel;
#pragma link C++ function GetSFuncRelSingle;

// Class from WjetsDataToMc.h
#pragma link C++ class WjetsDataToMc;

#endif
