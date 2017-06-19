#ifndef WASYMMETRY_REL17_H
#define WASYMMETRY_REL17_H

#include <cmath>
#include <string>

using namespace std;

double GetWscaleFactorPretag(int DataPlus, int DataMinus, double SumBkgPlus, double SumBkgMinus, double WmcPlus, double WmcMinus);

double GetWscaleFactorTag(int DataPlus, int DataMinus, double SumBkgPlus, double SumBkgMinus, double WmcPlus, double WmcMinus,
                          int DataPre2j, int DataTag2j, double SumBkgPre2j, double SumBkgTag2j, double ttPre2j, double ttTag2j,
                          double WPre2j, double WTag2j, double WPreNj, double WTagNj);

double GetWstatRelPretag(int DataPlus, int DataMinus, double SumBkgPlus, double SumBkgMinus);

double GetWstatRelf2j(int DataPre2j, int DataTag2j, double SumBkgPre2j, double SumBkgTag2j, double ErrSumBkgPre2j, 
                      double ErrSumBkgTag2j, double ttPre2j, double ttTag2j, double ErrttPre2j, double ErrttTag2j);

double GetWstatRelf2toN(double WPre2j, double WTag2j, double WPreNj, double WTagNj, double ErrWPre2j, double ErrWTag2j, 
                        double ErrWPreNj, double ErrWTagNj);

double GetWweight(int charge, double R, bool isData=true, bool isBkg=false);

double GetWweightSqErr(int charge, double R, bool isData=true);

double GetRvalue(int channel, float njet, int btag);     // R = W+/W- value from MC
double GetWvalue(int channel, float njet, int btag);    // W(data-driven) value, from 2fb-1 2011 data
double GetSFvalue(int channel, float njet, int btag);    // SF = W(data-driven)/W(MC) value, from 2fb-1 2011 data

double GetSFuncRel(string type="all", int channel=0, float njet=4.1, int btag=0);         // uncertainties are summed in quadrature
double GetSFuncRelSingle(string type="stat1", int channel=0, float njet=4.1, int btag=0); // single uncertainty, with sign

#endif
