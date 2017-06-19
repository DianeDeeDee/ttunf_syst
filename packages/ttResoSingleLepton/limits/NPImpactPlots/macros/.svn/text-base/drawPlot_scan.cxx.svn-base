// Author      : Stefan Gadatsch
// Email       : stefan.gadatsch@nikhef.nl
// Date        : 2013-10-25
// Description : Plot a likelihood scan with some information for debugging

#include "TFile.h"
#include "TTime.h"
#include "TSystem.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TColor.h"

#include "RooWorkspace.h"
#include "RooRealSumPdf.h"
#include "RooDataSet.h"
#include "RooSimultaneous.h"
#include "RooRealVar.h"
#include "RooMinimizer.h"

#include "RooStats/ModelConfig.h"

#include "inc/include_all.h"
#include "GlobIter.h"
#include "macros/drawPlot.C"

using namespace std;
using namespace RooFit;
using namespace RooStats;

// ____________________________________________________________________________|__________
const char* Scan(const char* a, std::vector<std::string>& v, char sep = ',', bool append = kTRUE);
void ShiftToZero(std::vector<double>& numbers);
double* getAry(std::vector<double> numbers);
TGraph* findIntersection(TGraph &a, TGraph &b);
TGraph* findShadedArea(TGraph &g, TGraph &interg, double maxpoint, double firstx);
void drawPlot_scan2_1D(std::string folder, std::string parameternames, TCanvas* c1, TPad* pad1);
void drawPlot_scan2_2D(std::string folder, std::string parameternames, TCanvas* c1, TPad* pad1);

// ____________________________________________________________________________|__________
void drawPlot_scan(std::string folder = "test", std::string parameternames = "SigXsecOverSM_HWW", bool remakeAscii = 0, std::string filenames = "")
{
	std::vector<string> parameterName;
	Scan(parameternames.c_str(), parameterName);

	if (remakeAscii) {
		system("mkdir -vp ascii/");
		stringstream outFileName;
		outFileName << "ascii/scan_" << folder << "_" << parameternames << ".txt";
		ofstream outFile(outFileName.str().c_str());

		std::vector<string> fileName;
		Scan(filenames.c_str(), fileName);

		std::vector< std::vector<double> > poival;
		std::vector<double> nllval;
		std::vector<int> status;

		for (size_t i = 0; i < fileName.size(); i++) {
			for (FileGlobIter files(fileName[i].c_str()); const char* fname = files();) {
				TFile* f = TFile::Open(fname);
				if (!f) {
					cout << "could not open file " << fname << endl;
					exit(-1);
				} else {
					cout << "open file " << fname << endl;

					TTree* resultTree = NULL;
					f->GetObject("result", resultTree);
					if (!resultTree) exit(-1);

					float NLL;
					RooFitResult*  fitresult = NULL;
					std::vector<std::string>* POIname = NULL;
					std::vector<double>* POIval = NULL;

					TBranch *b_NLL = NULL;
					TBranch *b_fitresult = NULL;
					TBranch *b_POIname = NULL;
					TBranch *b_POIval = NULL;

					resultTree->SetBranchAddress("NLL", &NLL, &b_NLL);
					resultTree->SetBranchAddress("fitresult", &fitresult, &b_fitresult);
					resultTree->SetBranchAddress("POIname", &POIname, &b_POIname);
					resultTree->SetBranchAddress("POIval", &POIval, &b_POIval);

					Long64_t tentry = resultTree->LoadTree(0);

					b_NLL->GetEntry(tentry);
					b_fitresult->GetEntry(tentry);
					b_POIname->GetEntry(tentry);
					b_POIval->GetEntry(tentry);

					int thisstatus = fitresult->status();
					float tmp = 2.0*NLL;

					std::vector<double> thisVals;

					cout << thisstatus;
					for (size_t j = 0; j < POIname->size(); j++) {
						thisVals.push_back(POIval->at(j));
						cout << " " << POIname->at(j) << "=" << thisVals[j];
					}
					cout << " NLL=" << NLL << endl;

					poival.push_back(thisVals);
					nllval.push_back(tmp);
					status.push_back(thisstatus);

					f->Close();
					delete f;
				}
			}
		}

		size_t n = poival.size();
		ShiftToZero(nllval);

		// sort the vector
		for (size_t i = 0; i < n-1; i++) {
			for (size_t j = 0; j < n-1-i; j++) {
				if (poival[j][0] > poival[j+1][0]) {
						swap(poival[j], poival[j+1]);
						swap(nllval[j], nllval[j+1]);
				}
			}
		}

		for (size_t i = 0; i < n; ++i) {
			outFile << i << " ";
			for (size_t j = 0; j < parameterName.size(); ++j) {
				outFile << poival[i][j] << " ";
			}
			outFile << nllval[i] << " " << status[i] << "\n";
		}
		outFile.close();
	}

	TCanvas* c1 = new TCanvas("c1","c1",1024,768);
	TPad *pad1 = new TPad("pad1", "pad1", 0.0  , 0.0  , 1.0 , 1.0  , 0);
	pad1->Draw();

	minMass = -1000.0;
	maxMass = 1000.0;
	if (parameterName.size() == 1) drawPlot_scan2_1D(folder, parameternames, c1, pad1);
	else if (parameterName.size() == 2) drawPlot_scan2_2D(folder, parameternames, c1, pad1);

	// pad1->cd();

	stringstream saveName;
	saveName << "scan_" << folder << "_" << parameternames;
	save(saveName.str(), "eps", c1);
	save(saveName.str(), "pdf", c1);
	save(saveName.str(), "C", c1);
}

// ____________________________________________________________________________|__________
void drawPlot_scan2_1D(std::string folder, std::string parameternames, TCanvas* c1, TPad* pad1)
{
	// read numbers to plot from text file
	ifstream testFile(("ascii/scan_"+folder+"_"+parameternames+".txt").c_str());
	if (testFile.fail()) {
		cout << "ERROR::file " << ("ascii/scan_"+folder+"_"+parameternames+".txt").c_str() << "does not exist.";
		exit(1);
	}
	fileHolder scan;
	drawPlot("ascii/scan_"+folder+"_"+parameternames+".txt", 3, scan);

	size_t nrPoints = scan.massPoints.size();
	vector<double> index = scan.massPoints;
	vector<double> poival = scan.getCol(0);
	vector<double> nllval = scan.getCol(1);
	vector<double> status = scan.getCol(2);

	// find maximum nll value
	double maxpoint = 0;
	for (size_t i = 0; i < nrPoints; ++i) {
		if (nllval[i] > maxpoint) maxpoint = nllval[i];
	}

	// multigraph that holds all objects to plot	
	TMultiGraph *mg = new TMultiGraph();
	mg->SetTitle((";"+parameternames+"; -2 ln #Lambda").c_str());

	// -2 ln Lambda graph
	TGraph* g = new TGraph(nrPoints, getAry(poival), getAry(nllval));
	g->SetMarkerSize(0.8);
	g->SetMarkerStyle(20);

	// bad fits use different marker
	std::vector<double> poival_stat1;
	std::vector<double> nllval_stat1;
	std::vector<double> poival_stat2;
	std::vector<double> nllval_stat2;

	for (size_t i = 0; i < nrPoints; ++i) {
		if (status[i] == 1) {
			poival_stat1.push_back(poival[i]);
			nllval_stat1.push_back(nllval[i]);
		}
		if (status[i] != 0 && status[i] != 1) {
			poival_stat2.push_back(poival[i]);
			nllval_stat2.push_back(nllval[i]);
		}
	}

	TGraph* g_stat1 = new TGraph(poival_stat1.size(), getAry(poival_stat1), getAry(nllval_stat1));
	g_stat1->SetMarkerSize(0.8);
	g_stat1->SetMarkerStyle(20);
	g_stat1->SetMarkerColor(kBlue-4);

	TGraph* g_stat2 = new TGraph(poival_stat2.size(), getAry(poival_stat2), getAry(nllval_stat2));
	g_stat2->SetMarkerSize(0.8);
	g_stat2->SetMarkerStyle(20);
	g_stat2->SetMarkerColor(kMagenta-3);

	// find best fit value as intersection with straight line at 0
	std::vector<double> vec0s;
	for (size_t i = 0; i < nrPoints; ++i) {
		vec0s.push_back(0.0);
	}

	TGraph* sig0 = new TGraph(nrPoints, getAry(poival), getAry(vec0s));
	TGraph *interg0 = findIntersection(*g,*sig0);
	interg0->SetMarkerStyle(21);
	interg0->SetMarkerColor(kRed+1);
	double firstx, firsty = 0;
	interg0->GetPoint(0, firstx, firsty);

	// find 1 sigma boundaries as intersection with straight line at 1
	std::vector<double> vec1s;
	for (size_t i = 0; i < nrPoints; ++i) {
		vec1s.push_back(1.0);
	}
	TGraph* sig1 = new TGraph(nrPoints, getAry(poival), getAry(vec1s));
	TGraph *interg1 = findIntersection(*g,*sig1);
	interg1->SetMarkerStyle(21);
	interg1->SetMarkerColor(kRed+1);

	// find 2 sigma boundaries as intersection with straight line at 4
	std::vector<double> vec2s;
	for (size_t i = 0; i < nrPoints; ++i) {
		vec2s.push_back(4.0);
	}
	TGraph* sig2 = new TGraph(nrPoints, getAry(poival), getAry(vec2s));
	TGraph *interg2 = findIntersection(*g,*sig2);
	interg2->SetMarkerStyle(21);
	interg2->SetMarkerColor(kRed+1);

	// get shaded 1 sigma area
	TGraph *grshade_1s = findShadedArea(*g,*interg1,maxpoint,firstx);
	grshade_1s->SetFillColor(kGreen-10);

	// get shaded 2 sigma area
	TGraph *grshade_2s = findShadedArea(*g,*interg2,maxpoint,firstx);
	grshade_2s->SetFillColor(kYellow-10);

	TLine l;
	l.SetLineWidth(2);
	l.SetLineColor(kRed+1);
	l.SetLineStyle(2);

	mg->Add(grshade_2s,"f");
	mg->Add(grshade_1s,"f");
	mg->Add(g,"LP");
	mg->Add(g_stat1,"P");
	mg->Add(g_stat2,"P");
	mg->Add(interg0,"P");
	mg->Add(interg1,"P");

	mg->Draw("A");
	gPad->Modified();
	mg->SetMinimum(0.);
	mg->SetMaximum(maxpoint);
	mg->GetXaxis()->SetRangeUser(poival[0]-(poival[poival.size()-1] - poival[0])/20, poival[poival.size()-1]+(poival[poival.size()-1] - poival[0])/20);

	for(int i = 0; i < interg1->GetN(); ++i) {
		double x, y = 0;
		interg1->GetPoint(i, x, y);
		l.DrawLine(x, 0.0, x, maxpoint);
	}

	for(int i = 0; i < interg2->GetN(); ++i) {
		double x, y = 0;
		interg2->GetPoint(i, x, y);
		l.DrawLine(x, 0.0, x, maxpoint);
	}

	TLatex t;
	double nsig = 1.0;
	while (nsig*nsig < maxpoint) {
		l.DrawLine(mg->GetXaxis()->GetXmin(), nsig*nsig, mg->GetXaxis()->GetXmax(), nsig*nsig);
		stringstream str;
		str << "#color[2]{" << static_cast<int>(nsig) << " #sigma}";
		t.DrawLatex(mg->GetXaxis()->GetXmax()+0.01, nsig*nsig, str.str().c_str());
		nsig+=1.0;
	}

	TLatex ll;
	ll.SetTextSize(0.015);
	ll.SetTextColor(kBlack);

	double high = poival[0]+(poival[poival.size()-1] - poival[0])/20;
	double low = poival[0]-(poival[poival.size()-1] - poival[0])/20;

	for(int i = 0; i < interg0->GetN(); ++i) {
		double x, y = 0;
		interg0->GetPoint(i, x, y);
		ll.DrawLatex(x - 0.0275*(high-low), y + 0.025*maxpoint, TString::Format("(%0.2f,%0.2f)", x, y));
	}

	bool evenodd = false;
	for(int i = 0; i < interg1->GetN(); ++i) {
		double x, y = 0;
		interg1->GetPoint(i, x, y);
		ll.DrawLatex((evenodd?(x + 0.02*(high-low)):(x - 0.08*(high-low))), y + 0.0125*maxpoint, TString::Format("(%0.2f,%0.2f)", x, y));
		evenodd = !evenodd;
	}

	evenodd = false;
	for(int i = 0; i < interg2->GetN(); ++i) {
		double x, y = 0;
		interg2->GetPoint(i, x, y);
		ll.DrawLatex((evenodd?(x + 0.02*(high-low)):(x - 0.08*(high-low))), y + 0.0125*maxpoint, TString::Format("(%0.2f,%0.2f)", x, y));
		evenodd = !evenodd;
	}

	c1->Update();
}

// ____________________________________________________________________________|__________
void drawPlot_scan2_2D(std::string folder, std::string parameternames, TCanvas* c1, TPad* pad1)
{
	gPad->SetRightMargin(0.175);

	std::vector<string> parameterName;
	Scan(parameternames.c_str(), parameterName);

	gStyle->SetOptStat(0);

	// read numbers to plot from text file
	ifstream testFile(("ascii/scan_"+folder+"_"+parameternames+".txt").c_str());
	if (testFile.fail()) {
		cout << "ERROR::file " << ("ascii/scan_"+folder+"_"+parameternames+".txt").c_str() << "does not exist.";
		exit(1);
	}
	fileHolder scan;
	drawPlot("ascii/scan_"+folder+"_"+parameternames+".txt", 4, scan);

	size_t nrPoints = scan.massPoints.size();
	vector<double> index = scan.massPoints;
	vector<double> poiAval = scan.getCol(0);
	vector<double> poiBval = scan.getCol(1);
	vector<double> nllval = scan.getCol(2);
	vector<double> status = scan.getCol(3);

	std::vector<double> minA;
	std::vector<double> minB;
	for (size_t i = 0; i < nrPoints; ++i) {
		if (nllval[i] == 0) {
			minA.push_back(poiAval[i]);
			minB.push_back(poiBval[i]);
			break;
		}
	}

	// finding bins
	std::vector<double> xbins;
	std::vector<double> ybins;
	for (size_t i = 0; i < nrPoints; ++i) {
		if (std::find(xbins.begin(), xbins.end(), poiAval[i]) == xbins.end()) {
			xbins.push_back(poiAval[i]);
		}
		if (std::find(ybins.begin(), ybins.end(), poiBval[i]) == ybins.end()) {
			ybins.push_back(poiBval[i]);
		}
	}

	std::sort (xbins.begin(), xbins.end());
	std::sort (ybins.begin(), ybins.end());

	size_t nrBinsX = xbins.size()-1;
	size_t nrBinsY = ybins.size()-1;

	double minX = 1000;
	double maxX = -1000;
	for (size_t i = 0; i < nrBinsX+1; ++i)	{
		if (xbins[i] > maxX) maxX = xbins[i];
		if (xbins[i] < minX) minX = xbins[i];
	}

	double minY = 1000;
	double maxY = -1000;
	for (size_t i = 0; i < nrBinsY+1; ++i)	{
		if (ybins[i] > maxY) maxY = ybins[i];
		if (ybins[i] < minY) minY = ybins[i];
	}

	double xspacing = fabs(xbins[1] - xbins[0]);
	for (size_t i = 0; i < nrBinsX+1; ++i)	{
		xbins[i] = xbins[i] - xspacing / 2;
	}
	
	double yspacing = fabs(ybins[1] - ybins[0]);
	for (size_t i = 0; i < nrBinsY+1; ++i)	{
		ybins[i] = ybins[i] - yspacing / 2;
	}

	// Create the 2D histogram for plotting
	TH2F *h2 = new TH2F("", "", xbins.size()-1, getAry(xbins), ybins.size()-1, getAry(ybins));
	h2->SetTitle((";"+parameterName[0]+";"+parameterName[1]+"; -2 ln #Lambda").c_str());

	for (size_t i = 0; i < nrPoints; ++i) {
		int bin = h2->FindBin(poiAval[i], poiBval[i]);
		h2->SetBinContent(bin, nllval[i]);
	}

	// palette settings
	const Int_t NRGBs = 4;
	const Int_t NCont = 999;

	const double max = 15.0;
	const double min = -max/NCont;
	
	double def1s = ROOT::Math::chisquared_quantile( 0.68, 2 );
	double def2s = ROOT::Math::chisquared_quantile( 0.95, 2 );
	double def3s = ROOT::Math::chisquared_quantile( 0.997, 2 );

	double stops[NRGBs] = { def1s/max, def2s/max, def3s/max, 1.00 };
	double red[NRGBs]   = { 1.00, 0.00, 0.40, 1.00 };
	double green[NRGBs] = { 1.00, 1.00, 0.80, 1.00 };
	double blue[NRGBs]  = { 0.40, 0.50, 1.00, 1.00 };

	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	gStyle->SetNumberContours(NCont);
  
	const Int_t nLevels = NCont;
	double levels[nLevels];

	for(int i = 1; i < nLevels; i++) {
		levels[i] = min + (max - min) / (nLevels - 1) * (i);
	}
	levels[0] = min;

	h2->SetContour((sizeof(levels)/sizeof(double)), levels);
	h2->DrawClone("COL");
	h2->GetZaxis()->SetRangeUser(min, max);
	h2->Draw("COL Z SAME");

	// find and draw n sigma contours on top
	double contours[3];
	contours[0] = def1s;
	contours[1] = def2s;
	contours[2] = def3s;

	TCanvas* tmpc = new TCanvas("c","Contour List",0,0,600,600);
	TH2F *h3 = new TH2F(*h2);
	h3->SetContour(3, contours);
	
	tmpc->cd();
	h3->Draw("CONT COL Z LIST");
	tmpc->Update();
	
	// Get Contours
	TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	TList* contLevel = NULL;
	TGraph* curv     = NULL;
	TGraph* gc       = NULL;

	Int_t nGraphs    = 0;
	Int_t TotalConts = 0;

	if (conts == NULL) {
		cout << "No Contours Were Extracted!" << endl;
		TotalConts = 0;
		return;
	} else {
		TotalConts = conts->GetSize();
	}

	for(int i = 0; i < TotalConts; i++){
		contLevel = (TList*)conts->At(i);
		nGraphs += contLevel->GetSize();
	}
	
	nGraphs = 0;

	// back to the original canvas
	c1->cd();

	double x0, y0, z0;
	TLatex l;
	l.SetTextSize(0.03);
	string val;

	for(int i = 0; i < TotalConts; i++){
		contLevel = (TList*)conts->At(i);
		if (i<3) z0 = contours[2-i];
		else     z0 = contours[i];

		// Get first graph from list on curves on this level
		curv = (TGraph*)contLevel->First();
		for(int j = 0; j < contLevel->GetSize(); j++){
			curv->GetPoint(0, x0, y0);
			curv->SetLineColor(kBlack);
			if (z0 == def3s) curv->SetLineStyle(1);
			else if (z0 == def2s) curv->SetLineStyle(2);
			else if (z0 == def1s) curv->SetLineStyle(3);
			nGraphs ++;

			// Draw clones of the graphs to avoid deletions in case the 1st pad is redrawn.
			gc = (TGraph*)curv->Clone();
			gc->Draw("C SAME");
			
			if (z0 == def3s) val = "1 #sigma";
			else if (z0 == def2s) val = "2 #sigma";
			else if (z0 == def1s) val = "3 #sigma";

			l.DrawLatex(x0,y0,val.c_str());
			curv = (TGraph*)contLevel->After(curv); // Get Next graph
		}
	}

	TGraph* bestfit = new TGraph(1, getAry(minA), getAry(minB));
	bestfit->SetMarkerStyle(34);
	bestfit->SetMarkerColor(kBlack);
	bestfit->Draw("P SAME");
	l.DrawLatex(minA[0],minB[0],TString::Format("(%0.2f,%0.2f)", minA[0], minB[0]));

	c1->Update();
}

// ____________________________________________________________________________|__________
const char* Scan(const char* a, std::vector<std::string>& v, char sep, bool append)
{
	if (!append) v.clear();
	if (!*a) return a;

	for (;; a++) {
		const char* b = strchr (a, sep);
		if (!b) b = a+strlen(a);
		v.push_back (std::string(a,b));
		if (!*b) return b;
		a = b;
	}
}

// ____________________________________________________________________________|__________
void ShiftToZero(std::vector<double>& numbers)
{
	// find minimum
	double min = numbers[0];
	for (size_t i = 0; i < numbers.size(); ++i) {
		if (numbers[i] < min) {
			min = numbers[i];
		}
	}

	for (size_t i = 0; i < numbers.size(); ++i) {
		numbers[i] -= min;
	}
}

// ____________________________________________________________________________|__________
// Return a TGraph with the points of intersection
TGraph* findIntersection(TGraph &a, TGraph &b)
{
	TGraph *interPoint = new TGraph();
	Int_t i = 0;
   
	// Loop over all points in this TGraph
	for(int a_i = 0; a_i < a.GetN()-1; ++a_i) {
		// Loop over all points in the other TGraph
		for(int b_i = 0; b_i < b.GetN()-1; ++b_i) {
		 
			// Get the current point, and the next point for each of the objects
			double x1, y1, x2, y2 = 0;
			double ax1, ay1, ax2, ay2 = 0;
			a.GetPoint(a_i, x1, y1);
			a.GetPoint(a_i+1, x2, y2);
			b.GetPoint(b_i, ax1, ay1);
			b.GetPoint(b_i+1, ax2, ay2);

			// Calculate the intersection between two straight lines, x axis
			double x = (ax1 *(ay2 *(x1-x2)+x2 * y1 - x1 * y2 )+ ax2 * (ay1 * (-x1+x2)- x2 * y1+x1 * y2)) / (-(ay1-ay2) * (x1-x2)+(ax1-ax2)* (y1-y2));

			// Calculate the intersection between two straight lines, y axis
			double y = (ax1 * ay2 * (y1-y2)+ax2 * ay1 * (-y1+y2)+(ay1-ay2) * (x2 * y1-x1 * y2))/(-(ay1-ay2) * (x1-x2)+(ax1-ax2) * (y1-y2));

			// Find the tightest interval along the x-axis defined by the four points
			double xrange_min = max(min(x1, x2), min(ax1, ax2));
			double xrange_max = min(max(x1, x2), max(ax1, ax2));

			if ((x1 == ax1 and y1 == ay1) or (x2 == ax2 and y2 == ay2)) {
				// If points from the two lines overlap, they are trivially intersecting
				interPoint->SetPoint(i, (x1 == ax1 and y1 == ay1) ? x1 : x2, (x1 == ax1 and y1 == ay1) ? y1 : y2);
				i++;
			} else if(x > xrange_min && x < xrange_max) {
				// If the intersection between the two lines is within the tight range, add it to the list of intersections.
				interPoint->SetPoint(i,x, y);
				i++;
			}
		}
	}

	return interPoint;
}

// ____________________________________________________________________________|__________
TGraph* findShadedArea(TGraph &g, TGraph &interg, double maxpoint, double firstx)
{
	std::vector<double> poival;
	std::vector<double> nllval;

	for (int i = 0; i < g.GetN(); ++i) {
		double x, y = 0;
		g.GetPoint(i, x, y);
		poival.push_back(x);
		nllval.push_back(y);
	}

	// add n sigma intersections not nominal llh scan values
	std::vector<double> poival_lo = poival;
	std::vector<double> nllval_lo = nllval;

	for (int i = 0; i < interg.GetN(); ++i) {
		double x, y = 0;
		interg.GetPoint(i, x, y);
		poival_lo.push_back(x);
		nllval_lo.push_back(y);
	}

	// sort the vector
	for (size_t i = 0; i < poival_lo.size()-1; i++) {
		for (size_t j = 0; j < poival_lo.size()-1-i; j++) {
			if (poival_lo[j] > poival_lo[j+1]) {
					swap(poival_lo[j], poival_lo[j+1]);
					swap(nllval_lo[j], nllval_lo[j+1]);
			}
		}
	}

	// clean points outside n sigma range
	int id_poival_lo = 0;
	int id_nllval_lo = 0;

	double interx, intery = 0;
	interg.GetPoint(0, interx, intery);
	if (intery == 0) intery = maxpoint;
	while (nllval_lo[id_poival_lo] > intery) {
		id_poival_lo++;
	}
	interg.GetPoint(interg.GetN()-1, interx, intery);
	if (intery == 0) intery = maxpoint;
	while (nllval_lo[nllval_lo.size()-id_nllval_lo-1] > intery) {
		id_nllval_lo++;
	}

	if (id_poival_lo > 0) {
		poival_lo.erase(poival_lo.begin(), poival_lo.begin()+id_poival_lo-1);
		nllval_lo.erase(nllval_lo.begin(), nllval_lo.begin()+id_poival_lo-1);
	}
	if (id_nllval_lo > 0) {
		poival_lo.erase(poival_lo.end()-id_nllval_lo+1, poival_lo.end());
		nllval_lo.erase(nllval_lo.end()-id_nllval_lo+1, nllval_lo.end());
	}

	// make n sigma bottom delimiter
	TGraph* grmin = new TGraph(poival_lo.size(), getAry(poival_lo), getAry(nllval_lo));
	
	// make n sigma top delimiter
	std::vector<double> poival_hi;
	std::vector<double> nllval_hi;
	for (int i = 0; i < interg.GetN(); ++i) {
		double x, y = 0;
		interg.GetPoint(i, x, y);
		poival_hi.push_back(x);
		nllval_hi.push_back(maxpoint);
	}
	
	if (poival_hi.size()==0) {
		poival_hi.push_back(poival[0]);
		nllval_hi.push_back(maxpoint);
	}
	if (poival_hi[0] > firstx) {
		poival_hi.push_back(poival[0]);
		nllval_hi.push_back(maxpoint);
	}

	// sort the vector
	for (size_t i = 0; i < poival_hi.size()-1; i++) {
		for (size_t j = 0; j < poival_hi.size()-1-i; j++) {
			if (poival_hi[j] > poival_hi[j+1]) {
				swap(poival_hi[j], poival_hi[j+1]);
				swap(nllval_hi[j], nllval_hi[j+1]);
			}
		}
	}

	if (poival_hi[poival_hi.size()-1] < firstx) {
		poival_hi.push_back(poival[poival.size()-1]);
		nllval_hi.push_back(maxpoint);
	}

	TGraph* grmax = new TGraph(poival_hi.size(), getAry(poival_hi), getAry(nllval_hi));
	
	// construct shaded n sigma area
	TGraph *grshade = new TGraph(grmin->GetN() + grmax->GetN());
	for (int i = 0; i < grmax->GetN(); ++i) {
		grshade->SetPoint(i,poival_hi[i],nllval_hi[i]);
		cout << i << " " << poival_hi[i] << " " << nllval_hi[i] << endl;
	}
	for (int i = 0; i < grmin->GetN(); ++i) {
		grshade->SetPoint(grmax->GetN()+i,poival_lo[grmin->GetN()-i-1],nllval_lo[grmin->GetN()-i-1]);
		cout << grmax->GetN()+i << " " << poival_lo[grmin->GetN()-i-1] << " " << nllval_lo[grmin->GetN()-i-1] << endl;
	}
	
	return grshade;
}
