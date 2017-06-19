#include "MObject.h"
#include "LargeJet.h"
#include <cmath>
#include <fastjet/PseudoJet.hh>
#include "TLorentzVector.h"

LargeJet::LargeJet()
  : MObject() {
  m_type = MObject::largejet;
}

LargeJet::LargeJet(const TLorentzVector &v)
  : MObject(v, MObject::largejet), m_bdrs(v) {
}

double &LargeJet::htt() {
  return m_htt;
}

const double LargeJet::htt() const {
  return m_htt;
}

LargeJet::LargeJet(const LargeJet &l)
  : MObject(l.mom(), MObject::largejet) {
  m_subjet = l.m_subjet;
  m_subjetJER = l.m_subjetJER;
  m_subjet_area = l.m_subjet_area;
  m_detE = l.m_detE;
  m_detEta = l.m_detEta;
  m_detPhi = l.m_detPhi;
  m_detM = l.m_detM;
  m_Ax = l.m_Ax;
  m_Ay = l.m_Ay;
  m_Az = l.m_Az;
  m_Ae = l.m_Ae;
  m_split12 = l.m_split12;
  m_trimmed_tau1 = l.m_trimmed_tau1;
  m_trimmed_tau2 = l.m_trimmed_tau2;
  m_ug_tau1 = l.m_ug_tau1;
  m_ug_tau2 = l.m_ug_tau2;
  m_bdrs = l.m_bdrs;
  //  m_trimmed_bdrs_tau1 = l.m_trimmed_bdrs_tau1;
  // m_trimmed_bdrs_tau2 = l.m_trimmed_bdrs_tau2;
  m_ug_bdrs_tau1 = l.m_ug_bdrs_tau1;
  m_ug_bdrs_tau2 = l.m_ug_bdrs_tau2;
  m_trimmed_am_tau1 = l.m_trimmed_am_tau1;
  m_trimmed_am_tau2 = l.m_trimmed_am_tau2;
  m_ug_am_tau1 = l.m_ug_am_tau1;
  m_ug_am_tau2 = l.m_ug_am_tau2;

  m_trueFlavour = l.m_trueFlavour;
  m_htt = l.m_htt;
}

LargeJet::~LargeJet() {
}

const int LargeJet::trueFlavour() const {
  return m_trueFlavour;
}

int &LargeJet::trueFlavour() {
  return m_trueFlavour;
}

const std::vector<LargeJet::Subjet> &LargeJet::subjet() const {
  return m_subjet;
}

std::vector<LargeJet::Subjet> &LargeJet::subjet() {
  return m_subjet;
}

const std::vector<LargeJet::Subjet> &LargeJet::subjetJER() const {
  return m_subjetJER;
}

std::vector<LargeJet::Subjet> &LargeJet::subjetJER() {
  return m_subjetJER;
}

const std::vector<TLorentzVector> &LargeJet::subjet_area() const {
  return m_subjet_area;
}

std::vector<TLorentzVector> &LargeJet::subjet_area() {
  return m_subjet_area;
}

bool LargeJet::pass() const {
  if (std::fabs(mom().Eta()) > 1.2) return false;
  if (mom().Perp() < 200e3) return false;
  return true;
}

bool LargeJet::passLoose() const {
  if (std::fabs(mom().Eta()) > 2.0) return false;
  //if (mom().Perp() < 200e3) return false;
  if (mom().Perp() < 160e3) return false;
  return true;
}

double &LargeJet::detE() {
  return m_detE;
}

double &LargeJet::detEta() {
  return m_detEta;
}
double &LargeJet::detPhi() {
  return m_detPhi;
}
double &LargeJet::detM() {
  return m_detM;
}
double &LargeJet::Ax() {
  return m_Ax;
}
double &LargeJet::Ay() {
  return m_Ay;
}
double &LargeJet::Az() {
  return m_Az;
}

double &LargeJet::Ae() {
  return m_Ae;
}

double &LargeJet::split12() {
  return m_split12;
}
const double LargeJet::split12() const {
  return m_split12;
}


double &LargeJet::trimmed_tau1() {
  return m_trimmed_tau1;
}
const double LargeJet::trimmed_tau1() const {
  return m_trimmed_tau1;
}

double &LargeJet::trimmed_tau2() {
  return m_trimmed_tau2;
}
const double LargeJet::trimmed_tau2() const {
  return m_trimmed_tau2;
}



double &LargeJet::trimmed_am_tau1() {
  return m_trimmed_am_tau1;
}
const double LargeJet::trimmed_am_tau1() const {
  return m_trimmed_am_tau1;
}

double &LargeJet::trimmed_am_tau2() {
  return m_trimmed_am_tau2;
}
const double LargeJet::trimmed_am_tau2() const {
  return m_trimmed_am_tau2;
}


/*

double &LargeJet::trimmed_bdrs_tau1() {
  return m_trimmed_bdrs_tau1;
}
const double LargeJet::trimmed_bdrs_tau1() const {
  return m_trimmed_bdrs_tau1;
}

double &LargeJet::trimmed_bdrs_tau2() {
  return m_trimmed_bdrs_tau2;
}
const double LargeJet::trimmed_bdrs_tau2() const {
  return m_trimmed_bdrs_tau2;
}

*/


double &LargeJet::ug_tau1() {
  return m_ug_tau1;
}
const double LargeJet::ug_tau1() const {
  return m_ug_tau1;
}

double &LargeJet::ug_tau2() {
  return m_ug_tau2;
}
const double LargeJet::ug_tau2() const {
  return m_ug_tau2;
}



double &LargeJet::ug_am_tau1() {
  return m_ug_am_tau1;
}
const double LargeJet::ug_am_tau1() const {
  return m_ug_am_tau1;
}

double &LargeJet::ug_am_tau2() {
  return m_ug_am_tau2;
}
const double LargeJet::ug_am_tau2() const {
  return m_ug_am_tau2;
}


double &LargeJet::ug_bdrs_tau1() {
  return m_ug_bdrs_tau1;
}
const double LargeJet::ug_bdrs_tau1() const {
  return m_ug_bdrs_tau1;
}

double &LargeJet::ug_bdrs_tau2() {
  return m_ug_bdrs_tau2;
}
const double LargeJet::ug_bdrs_tau2() const {
  return m_ug_bdrs_tau2;
}





TLorentzVector &LargeJet::bdrs() {
  return m_bdrs;
}
const TLorentzVector LargeJet::bdrs() const {
  return m_bdrs;
}

const double LargeJet::bdrs_m() const {
  return m_bdrs.M();
}


const double LargeJet::bdrs_pt() const {
  return m_bdrs.Perp();
}


const double LargeJet::bdrs_phi() const {
  return m_bdrs.Phi();
}


const double LargeJet::bdrs_eta() const {
  return m_bdrs.Eta();
}



