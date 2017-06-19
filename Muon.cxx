#include "MObject.h"
#include "Muon.h"
#include <cmath>
//If no QCD
//Muon::Muon()
//  : MObject(), m_mi(-1), m_tight(false), m_z0(0), m_d0(0), m_sd0(1), m_author(0), m_passTrkCuts(false) {
//  m_type = MObject::mu;
//}
//
//Muon::Muon(const TLorentzVector &v) //If no QCD
//  : MObject(v, MObject::mu), m_mi(-1), m_tight(false), m_z0(0), m_d0(0), m_sd0(1), m_author(0), m_passTrkCuts(false) {
//}
//For QCD
Muon::Muon()
  : MObject(), m_mi(-1), m_tight(false), m_z0(0), m_d0(0), m_sd0(1),m_z0_exPV(0), m_author(0), m_passTrkCuts(false) {
  m_type = MObject::mu;
}

Muon::Muon(const TLorentzVector &v)
  : MObject(v, MObject::mu), m_mi(-1), m_tight(false), m_z0(0), m_d0(0), m_sd0(1),m_z0_exPV(0), m_author(0), m_passTrkCuts(false) {
}
//End For QCD
Muon::~Muon() {
}

void Muon::setMI(float iso) {
  m_mi = iso;
}

float Muon::mi() const {
  return m_mi;
}

void Muon::setTight(bool t) {
  m_tight = t;
}

bool Muon::isTight() const {
  return m_tight;
}

const float Muon::z0() const {
  return m_z0;
}

float &Muon::z0() {
  return m_z0;
}

const float Muon::d0() const {
  return m_d0;
}

float &Muon::d0() {
  return m_d0;
}

const float Muon::sd0() const {
  return m_sd0;
}
float &Muon::sd0() {
  return m_sd0;
}
//For QCD
const float Muon::z0_exPV() const {
  return m_z0_exPV;
}
float &Muon::z0_exPV() {
  return m_z0_exPV;
}
//End For QCD
const int Muon::author() const {
  return m_author;
}
int &Muon::author() {
  return m_author;
}

const bool Muon::passTrkCuts() const {
  return m_passTrkCuts;
}

bool &Muon::passTrkCuts() {
  return m_passTrkCuts;
}


bool Muon::pass() const { //if no QCD
  if (mom().Perp() < 25e3) return false;
  float eta = std::fabs(mom().Eta());
  if (eta > 2.5) return false;
  if (!isTight()) return false;
  if (author() != 12) return false;// _author means which atlas algorithm is used to reconstruct them
  if (!passTrkCuts()) return false;
  if (mi()/mom().Perp() > 0.05) return false; //for mu isolation
  if (std::fabs(d0()/sd0()) > 3) return false; //to reduce multijet bkg
  if (std::fabs(z0_exPV()) > 2) return false;
 // if (std::fabs(z0()) > 2) return false;

  return true;
}

//For QCD
bool Muon::passLoose() const {
  if (mom().Perp() < 25e3) return false;
  float eta = std::fabs(mom().Eta());
  if (eta > 2.5) return false;
  if (!isTight()) return false; 
  if (author() != 12) return false;// _author means which atlas algorithm is used to reconstruct them. [2015-04-01 16:51:46] Zhong Jia-hang: 12 is this "Muid algorithm using combined information of ID/MS"[2015-04-01 16:52:02] Zhong Jia-hang: the other family is Staco
  if (!passTrkCuts()) return false;
  if (std::fabs(d0()/sd0()) > 3) return false; //to reduce multijet bkg
  if (std::fabs(z0_exPV()) > 2) return false;
  //if (mi()/mom().Perp() > 0.05) return false;
  return true;
}

//End for QCD

TLorentzVector &Muon::momME() {
  return m_momME;
}

const TLorentzVector &Muon::momME() const {
  return m_momME;
}

TLorentzVector &Muon::momMECorr() {
  return m_momMECorr;
}
const TLorentzVector &Muon::momMECorr() const {
  return m_momMECorr;
}

TLorentzVector &Muon::momMS() {
  return m_momMS;
}
const TLorentzVector &Muon::momMS() const {
  return m_momMS;
}

TLorentzVector &Muon::momID() {
  return m_momID;
}
const TLorentzVector &Muon::momID() const {
  return m_momID;
}

TLorentzVector &Muon::momTrk() {
  return m_momTrk;
}
const TLorentzVector &Muon::momTrk() const {
  return m_momTrk;
}

int Muon::charge() const {
  return m_charge;
}
int &Muon::charge() {
  return m_charge;
}

bool Muon::setST(bool b) {
  m_st = b;
}
bool Muon::st() const {
  return m_st;
}

bool Muon::setSA(bool b) {
  m_sa = b;
}
bool Muon::sa() const {
  return m_sa;
}

bool Muon::setCB(bool b) {
  m_cb = b;
}

bool Muon::cb() const {
  return m_cb;
}
