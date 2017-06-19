#include "MObject.h"
#include "Electron.h"
#include <cmath>
#include "Jet.h"
#include <iostream>

Electron::Electron()
  : MObject(), m_mi(-1), m_isTightPP(false), m_mom_calo(0,0,0,0), m_mom_trk(0,0,0,0) {
  m_type = MObject::e;
  m_passOR = true;
}

Electron::Electron(const TLorentzVector &v)
  : MObject(v, MObject::e), m_mi(-1), m_isTightPP(false), m_mom_calo(v), m_mom_trk(v) {
}

Electron::~Electron() {
}

bool &Electron::passOR() {
  return m_passOR;
}

void Electron::setMI(float iso) {
  m_mi = iso;
}

float Electron::mi() const {
  return m_mi;
}

void Electron::setTightPP(bool t) {
  m_isTightPP = t;
}

bool Electron::isTightPP() const {
  return m_isTightPP;
}

void Electron::setMediumPP(bool t) {
  m_isMediumPP = t;
}

bool Electron::isMediumPP() const {
  return m_isMediumPP;
}

void Electron::setLoosePP(bool t) {
  m_isLoosePP = t;
}

bool Electron::isLoosePP() const {
  return m_isLoosePP;
}

const TLorentzVector &Electron::caloMom() const {
  return m_mom_calo;
}

TLorentzVector &Electron::caloMom() {
  return m_mom_calo;
}

const TLorentzVector &Electron::trkMom() const {
  return m_mom_trk;
}

TLorentzVector &Electron::trkMom() {
  return m_mom_trk;
}

const float Electron::z0() const {
  return m_z0;
}

float &Electron::z0() {
  return m_z0;
}

const int Electron::author() const {
  return m_author;
}

int &Electron::author() {
  return m_author;
}

const int Electron::nSiHits() const {
  return m_nSiHits;
}

int &Electron::nSiHits() {
  return m_nSiHits;
}

const int Electron::oq() const {
  return m_oq;
}

int &Electron::oq() {
  return m_oq;
}

const int Electron::isEM() const {
  return m_isEM;
}

int &Electron::isEM() {
  return m_isEM;
}

const int Electron::GSF_trk_index() const {
  return m_GSF_trk_index;
}

int &Electron::GSF_trk_index() {
  return m_GSF_trk_index;
}

bool Electron::pass() const {
  if (!m_passOR) return false;
  if (mom().Perp() < 25e3) return false;
  if ( (author() != 1) && (author() != 3) ) return false;
  float eta = std::fabs(caloMom().Eta());
  //if (eta == 0) return false;
  if (eta > 1.37 && eta < 1.52) return false;
  if (eta > 2.47) return false;
  if (!isTightPP()) return false;
  if (std::fabs(z0()) > 2) return false;
  if (mi()/mom().Perp() > 0.05) return false;
  if ( (oq() & 1446) != 0 ) return false;
  return true;
}

bool Electron::passLoose(bool looser) const {
  if (!m_passOR) return false;
  if (mom().Perp() < 25e3) return false;
  if ( (author() != 1) && (author() != 3) ) return false;
  float eta = std::fabs(caloMom().Eta());
  if (eta > 1.37 && eta < 1.52) return false;
  if (eta > 2.47) return false;
  if (std::fabs(z0()) > 2) return false;
  if ( (oq() & 1446) != 0 ) return false;

  //if (mi()/mom().Perp() > 0.05) return false;
  //if (!isTightPP()) return false;
  if (looser)
    if (!isLoosePP()) return false;
  else
    if (!isMediumPP()) return false;
  if (isEM() & (1<<0x1) != 0) return false; 
  return true;
}

