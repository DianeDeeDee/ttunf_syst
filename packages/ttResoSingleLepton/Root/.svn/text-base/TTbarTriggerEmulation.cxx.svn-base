// Dear emacs, this is -*- c++ -*-
#include "ttResoSingleLepton/TTbarTriggerEmulation.h"

#include "TMath.h"
#include <iostream>
#include <string>
#include "TString.h"

using namespace std;

//___________________________________________________
TTbarTriggerEmulation::TTbarTriggerEmulation():
  TObject(),
  m_debug(false)
{
}
//___________________________________________________
TTbarTriggerEmulation::TTbarTriggerEmulation(const TTbarTriggerEmulation &ref):
  TObject(ref),
  m_debug(ref.m_debug)
{
}
//___________________________________________________
TTbarTriggerEmulation::~TTbarTriggerEmulation()
{
}
//___________________________________________________
bool TTbarTriggerEmulation::GetTrigger(const int jet_n, const vector<float> &trig_EF_jet_emscale_E, const vector<float> &trig_EF_jet_emscale_eta, const vector<string> &trig_EF_jet_calibtags, const float treshold, const string calibtag, const bool L1, const bool L2, std::vector<int> *trig_EF_jet_decision){


  if (trig_EF_jet_decision) {
    trig_EF_jet_decision->resize(jet_n);
    for(int i=0; i<jet_n; i++) trig_EF_jet_decision->at(i)=0;
    if (m_debug) std::cout<< "TTbarTriggerEmulation::GetTrigger : size of trig_EF_jet_decision="<<trig_EF_jet_decision->size()<<endl;
  }

  if (!L1 || !L2) return false;
  if(m_debug) std::cout << "TTbarTriggerEmulation::GetTrigger : pass L1 & L2" << std::endl;

  bool passEF=false;

  if(calibtag.compare("AntiKt10_topo")==0){

    for(int i=0; i<jet_n; i++)
      {
	if (trig_EF_jet_decision) trig_EF_jet_decision->at(i)=false;
	
	float eT = trig_EF_jet_emscale_E.at(i) / TMath::CosH( trig_EF_jet_emscale_eta.at(i) ) ;
	if(m_debug) std::cout << "TTbarTriggerEmulation::GetTrigger : eT=" << eT << "\t treshold *1000=" << treshold *1000 << std::endl;
	if( eT > treshold *1000 /* GeV->MeV */ 
	    && TMath::Abs( trig_EF_jet_emscale_eta.at(i) ) < 3.2) { 
	  if (trig_EF_jet_decision) trig_EF_jet_decision->at(i)=true;
	  if (calibtag.compare(trig_EF_jet_calibtags.at(i))==0) passEF=true;
	}
	
      }
  }
  else{
    TString dummy(calibtag);
    printf("  This calibtag %s has not been configurured",dummy.Data());
    return false;
  }

  return passEF;
}
