// Author      : Stefan Gadatsch
// Email       : stefan.gadatsch@nikhef.nl
// Date        : 2013-10-25
// Description : Perform a likelihood scan

#include "TFile.h"
#include "TTime.h"
#include "TSystem.h"
#include "TTree.h"

#include "RooWorkspace.h"
#include "RooRealSumPdf.h"
#include "RooDataSet.h"
#include "RooSimultaneous.h"
#include "RooRealVar.h"
#include "RooMinimizer.h"
#include "RooGaussModel.h"

#include "RooStats/ModelConfig.h"

#include "macros/minimize.C"
#include "macros/log.C"
#include "macros/parseString.C"

using namespace std;
using namespace RooFit;
using namespace RooStats;

// ____________________________________________________________________________|__________
void PrintResourcesUsed(const TTime& progStart);

// ____________________________________________________________________________|__________
void runLLHScan() {
	// for compiling only
}

// ____________________________________________________________________________|__________
void runLLHScan(const char* inFileName,
	const char* poiName = "SigXsecOverSM_HWW",
	const char* poiVal = "1.0",
	const char* wsName = "combined",
	const char* modelConfigName = "ModelConfig",
	const char* dataName = "asimovData_1",
	const char* snapshot = "nominalNuis",
	const char* folder = "test",
	int numCPU = 1,
	bool binnedevaluation = kFALSE,
	string loglevel = "DEBUG",
	const char* profileName = "")
{
	TTime thistime = gSystem->Now();

	// DEBUG OUTPUT
	// - ERROR
	// - WARNING
	// - INFO
	// - DEBUG
	LOG::ReportingLevel() = LOG::FromString(loglevel);

	// some settings
	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
	ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
	if (loglevel == "DEBUG") {
		ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1);
	} else {
		ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
		RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
	}

	// loading the workspace etc.
	LOG(logINFO) << "Running over workspace: " << inFileName;
	system(("mkdir -vp root-files/" + string(folder)).c_str());
	vector<string> parsedPoiName = parseString(poiName, ",");
	vector<string> parsedPoiVal = parseString(poiVal, ",");
	if (parsedPoiName.size() != parsedPoiVal.size()) {
		LOG(logERROR) << "Number of POIs does not match number of given values!";
		exit(1);
	}
	vector<string> parsedProfileName = parseString(profileName, ",");

	TFile* file = new TFile(inFileName);

	RooWorkspace* ws = (RooWorkspace*)file->Get(wsName);
	if (!ws) {
		LOG(logERROR) << "Workspace: " << wsName << " doesn't exist!";
		exit(1);
	}

	// Activate binned likelihood calculation for binned models
	if (binnedevaluation) {
		LOG(logINFO) << "Activating binned likelihood evaluation";
		RooFIter iter = ws->components().fwdIterator() ;
		RooAbsArg* arg ;
		while((arg = iter.next())) {
			if (arg->IsA() == RooRealSumPdf::Class()) {
				arg->setAttribute("BinnedLikelihood");
			}
		}
	}

	ModelConfig* mc = (ModelConfig*)ws->obj(modelConfigName);
	if (!mc) {
		LOG(logERROR) << "ModelConfig: " << modelConfigName << " doesn't exist!";
		exit(1);
	}

	RooDataSet* data = (RooDataSet*)ws->data(dataName);
	if (!data) {
		LOG(logERROR) << "Dataset: " << dataName << " doesn't exist!";
		exit(1);
	}

	const RooArgSet* globs = mc->GetGlobalObservables();
	if (!globs) {
		LOG(logERROR) << "Global observables don't exist!";
		exit(1);
	}

	bool loadedsnapshot = ws->loadSnapshot(snapshot);
	if (string(dataName).find("asimovData_muhat") != string::npos) ws->loadSnapshot("conditionalGlobs_muhat");
	
	TIterator* nitr = mc->GetParametersOfInterest()->createIterator();
	RooRealVar* var;
	nitr->Reset();
	while ((var = (RooRealVar*)nitr->Next())) {
    	cout << "setting " << var->GetName() << " constant" << endl;
		var->setRange(-5.0,5.0);
		var->setVal(1.0);
		var->setConstant(1);
	}

	const RooArgSet* nuis = mc->GetNuisanceParameters();    

	vector<RooRealVar*> pois;
	vector<double> poiVals;
	for (int i = 0; i < parsedPoiName.size(); i++) {
		RooRealVar* thisPoi = (RooRealVar*)ws->var(parsedPoiName[i].c_str());
		double thisPoiVal = atof(parsedPoiVal[i].c_str());
		poiVals.push_back(thisPoiVal);
		LOG(logINFO) << "Getting POI " << thisPoi->GetName() << " and set value to " << thisPoiVal;
		if (!thisPoi) {
			LOG(logERROR) << "POI: " << thisPoi->GetName() << " doesn't exist!";
			exit(1);
		}
		thisPoi->setRange(-1000.,1000.);
		thisPoi->setError(0.2);
		thisPoi->setVal(thisPoiVal);
		thisPoi->setConstant(1);
		pois.push_back(thisPoi);
	}


	for (int i = 0; i < parsedProfileName.size(); i++) {
		TString thisName = parsedProfileName[i];
		int sign = 0;
		if (thisName.Contains("+")) {
			sign = +1;
			thisName.ReplaceAll("+","");
		} else if (thisName.Contains("-")) {
			sign = -1;
			thisName.ReplaceAll("-","");
		}
		RooRealVar* thisPoi = (RooRealVar*)ws->var(thisName);
		LOG(logINFO) << "Getting POI to profile " << thisPoi->GetName() << " and adjusting range " << ((sign > 0)?"positive":"negative");
		if (!thisPoi) {
			LOG(logERROR) << "POI: " << thisPoi->GetName() << " doesn't exist!";
			exit(1);
		}
		thisPoi->setRange(-5.0, 5.0);
		if (sign > 0) {
			thisPoi->setRange(0.0, thisPoi->getMax());
			thisPoi->setVal(1.0);
		} else if (sign < 0) {
			thisPoi->setRange(thisPoi->getMin(), 0.0);
			thisPoi->setVal(-1.0);
		}
		thisPoi->setError(0.2);
		thisPoi->setConstant(0);
	}

/*
	ws->var("PZ")->setVal(1.0);
	ws->var("PW")->setVal(1.0);
	ws->var("Ptop")->setVal(1.0);
	ws->var("Pb")->setVal(1.0);
	ws->var("Ptau")->setVal(1.0);
	ws->var("Pmu")->setVal(1.0);
	ws->var("PZga")->setVal(1.0);
	
	ws->var("PZ")->setConstant(1);
	ws->var("PW")->setConstant(1);
	ws->var("Ptop")->setConstant(1);
	ws->var("Pb")->setConstant(1);
	ws->var("Ptau")->setConstant(1);
	ws->var("Pmu")->setConstant(1);
	ws->var("PZga")->setConstant(1);
*/

	RooNLLVar* nll = (RooNLLVar*)mc->GetPdf()->createNLL(*data, Constrain(*nuis), GlobalObservables(*globs), Offset(1), NumCPU(numCPU, RooFit::Hybrid));

	cout.precision(15);
	cout << "NLL before minimisation: " << nll->getVal() << endl;
	RooFitResult* fitresult = minimize(nll);
	
	nll->enableOffsetting(0);

	double minnll = nll->getVal();
	cout.precision(15);
	cout << "NLL after minimisation: " << minnll << endl;

	stringstream filename;
	filename << "root-files/" << folder << "/scan";
	for (int i = 0; i < pois.size(); ++i) {
		filename << "_" << pois[i]->GetName() << pois[i]->getVal();
	}
	filename << ".root";
	TFile f(filename.str().c_str(), "recreate");

	TTree* resultTree = new TTree("result", "result");
	resultTree->SetDirectory(0);
	// resultTree->Branch("fitresult", &fitresult);
	resultTree->Branch("NLL", &minnll);
	resultTree->Branch("POIname", &parsedPoiName);
	resultTree->Branch("POIval", &poiVals);
	resultTree->Fill();
	resultTree->ResetBranchAddresses();
	resultTree->Write("",TObject::kOverwrite);
  	f.Close();

  	// fitresult->Print("v");

	PrintResourcesUsed(thistime);
}

// ____________________________________________________________________________|__________
// Print used resources 
void PrintResourcesUsed(const TTime& progStart)
{
	ProcInfo_t info;
	if (gSystem->GetProcInfo(&info)<0) return;
	Long_t cput= TMath::CeilNint(info.fCpuUser);
	Long_t wall= Long64_t(gSystem->Now()-progStart+TTime(500))/Long64_t(1000);
	LOG(logINFO) << Form("resources used: cput=%02ld:%02ld:%02ld, mem=%ldkb, vmem=%ldkb, walltime=%02ld:%02ld:%02ld",
		cput/3600, (cput/60)%60, cput%60,
		info.fMemResident, info.fMemVirtual,
		wall/3600, (wall/60)%60, wall%60);
}
