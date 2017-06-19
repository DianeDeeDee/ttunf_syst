#include "MObject.h"
#include "LargeJetBB.h"
#include <cmath>
#include <fastjet/PseudoJet.hh>
#include "TLorentzVector.h"

LargeJetBB::LargeJetBB()
  : MObject() {
  m_type = MObject::largejetBB;
}

LargeJetBB::LargeJetBB(const TLorentzVector &v)
  : MObject(v, MObject::largejetBB) {
}

LargeJetBB::LargeJetBB(const LargeJetBB &l)
  : MObject(l.mom(), MObject::largejet) {
  m_subjet = l.m_subjet;
  m_subjet_area = l.m_subjet_area;
}

LargeJetBB::~LargeJetBB() {
}

const std::vector<TLorentzVector> &LargeJetBB::subjet() const {
  return m_subjet;
}

std::vector<TLorentzVector> &LargeJetBB::subjet() {
  return m_subjet;
}

const std::vector<TLorentzVector> &LargeJetBB::subjet_area() const {
  return m_subjet_area;
}

std::vector<TLorentzVector> &LargeJetBB::subjet_area() {
  return m_subjet_area;
}

bool LargeJetBB::pass() const {
  if (std::fabs(mom().Eta()) > 1.2) return false;
  if (mom().Perp() < 200e3) return false;
  return true;
}




double &LargeJetBB::detE() {
  return m_detE;
}

double &LargeJetBB::detEta() {
  return m_detEta;
}
double &LargeJetBB::detPhi() {
  return m_detPhi;
}
double &LargeJetBB::detM() {
  return m_detM;
}
double &LargeJetBB::Ax() {
  return m_Ax;
}
double &LargeJetBB::Ay() {
  return m_Ay;
}
double &LargeJetBB::Az() {
  return m_Az;
}

double &LargeJetBB::Ae() {
  return m_Ae;
}



double &LargeJetBB::ug_tau1() {
  return m_ug_tau1;
}


double &LargeJetBB::ug_tau2() {
  return m_ug_tau2;
}


double &LargeJetBB::bdrs_m() {
  return m_bdrs_m;
}


double &LargeJetBB::bdrs_pt() {
  return m_bdrs_pt;
}


double &LargeJetBB::bdrs_tau1() {
  return m_bdrs_tau1;
}


double &LargeJetBB::bdrs_tau2() {
  return m_bdrs_tau2;
}
