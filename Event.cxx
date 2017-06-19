#include "Event.h"
#include <vector>
#include "RawReader.h"
#include "Correction.h"
#include "Particle.h"
#include <cmath>
#include "SkimReader.h"

Event::Event() {
  m_sr = 0;
}

Event::~Event() {
}


const TLorentzVector &Event::triggerJet4mom() const {
  return m_triggerJet4mom;
}

TLorentzVector &Event::triggerJet4mom() {
  return m_triggerJet4mom;
}

const bool Event::isTight() const {
  return m_isTight;
}

bool &Event::isTight() {
  return m_isTight;
}
/*
const int Event::hfor() const {
  return m_hfor;
}

int &Event::hfor() {
  return m_hfor;
}
*/
const TLorentzVector &Event::metCorr(const std::string &s) const {
  std::map<std::string, TLorentzVector>::const_iterator it = m_metCorr.find(s);
  if (it == m_metCorr.end())
    return m_met;
  return it->second;
}

TLorentzVector &Event::metCorr(const std::string &s, bool create) {
  if (create)
    return m_metCorr[s];

  if (m_metCorr.find(s) == m_metCorr.end())
    return m_met;
  return m_metCorr[s];
}

void Event::clear() {
  m_cutFlow.clear();
  m_electron.clear();
  m_muon.clear();
  m_jet.clear();
  m_largeJet.clear();
  m_largeJetBB.clear();
  m_met.SetPxPyPzE(0,0,0,0);
  m_metCorr.clear();

  m_partMom.clear();
  m_partElectron.clear();
  m_partMuon.clear();
  m_partJet.clear();
  m_partLargeJet.clear();
  m_partLargeJetBB.clear();
  m_partMet.SetPxPyPzE(0,0,0,0);
  clearWeights();
}

std::vector<float> &Event::cutFlow() {
  return m_cutFlow;
}

bool &Event::passReco() {
  return m_passReco;
}

bool &Event::passPart() {
  return m_passPart;
}

bool Event::passReco() const {
  return m_passReco;
}
bool Event::passPart() const {
  return m_passPart;
}

void Event::read(Reader &rr) {
  rr.electron(m_electron);
  rr.muon(m_muon);
  rr.jet(m_jet, m_largeJet);//, m_largeJetBB);
  float mx, my;
  rr.met(mx, my);
  m_met.SetPxPyPzE(mx, my, 0, std::sqrt(mx*mx+my*my));

  rr.trigger(m_triggerElectron, m_triggerMuon, m_triggerLargeJet);//, m_triggerElectron24, m_triggerElectron60); 

  rr.general(m_channelNumber, m_isData, m_runNumber, m_eventNumber, m_mu, m_npv, m_npv_met, m_rho, m_lerr, m_terr, m_cfl, m_atlFastII, m_period, m_vxZ, m_mcWeight, m_lbn, m_npv_good);
  m_sr = rr.sr();

  if (!m_isData) {
    rr.partMom(m_partMom);
    rr.partElectron(m_partElectron);
    rr.partMuon(m_partMuon);
    rr.partJet(m_partJet, m_partLargeJet);//, m_partLargeJetBB);
    float pmx, pmy;
    rr.partMet(pmx, pmy);
    m_partMet.SetPxPyPzE(pmx, pmy, 0, std::sqrt(pmx*pmx+pmy*pmy));
  }

  //rr.sr(m_sr);

}

const SkimReader *Event::sr() const {
  return m_sr;
}

void Event::setupSystVariation(const std::string &s) {
  for (int k = 0; k < electron().size(); ++k) {
    electron()[k].mom() = electron()[k].corr(s);
  }
  for (int k = 0; k < muon().size(); ++k) {
    muon()[k].mom() = muon()[k].corr(s);
  }
  for (int k = 0; k < jet().size(); ++k) {
    jet()[k].mom() = jet()[k].corr(s);
    if (s == "jee" || s == "jer" || s == "jesUp" || s == "jesDown")
      jet()[k].corr("original") = jet()[k].corr(Form("original%s", s.c_str()));
    else
      jet()[k].corr("original") = jet()[k].corr(Form("original%s", "nom"));
  }
  for (int k = 0; k < largeJet().size(); ++k) {
    largeJet()[k].mom() = largeJet()[k].corr(s);
  }
  TLorentzVector m = metCorr(s);
  met(m.Px(), m.Py());
}

void Event::clearWeights() {
  m_weight.clear();
}

std::vector<std::string> Event::weightNames() const {
  std::vector<std::string> n;
  for (std::vector< std::pair<std::string, std::vector<double> > >::const_iterator it = m_weight.begin(); it != m_weight.end(); ++it) {
    n.push_back(it->first);
  }
  return n;
}

float &Event::mcWeight() {
  return m_mcWeight;
}

const float Event::mcWeight() const {
  return m_mcWeight;
}

void Event::correct(Correction &c) {
  c.run(*this);
}

std::vector<Electron> &Event::electron() {
  return m_electron;
}

std::vector<Muon> &Event::muon() {
  return m_muon;
}

std::vector<Jet> &Event::jet() {
  return m_jet;
}

std::vector<LargeJet> &Event::largeJet() {
  return m_largeJet;
}

std::vector<LargeJet> &Event::largeJetBB() {
  return m_largeJetBB;
}


void Event::met(float met_x, float met_y) {
  m_met.SetPxPyPzE(met_x, met_y, 0, std::sqrt(std::pow(met_x, 2) + std::pow(met_y, 2)));
}

TLorentzVector Event::met() const {
  return m_met;
}
/*
const bool Event::triggerElectron24() const {
  return m_triggerElectron24;
}

const bool Event::triggerElectron60() const {
  return m_triggerElectron60;
}
*/
const bool Event::triggerElectron() const {
  return m_triggerElectron;
}

const bool Event::triggerMuon() const {
  return m_triggerMuon;
}

const bool Event::triggerLargeJet() const {
  return m_triggerLargeJet;
}
/*
bool &Event::triggerElectron60() {
  return m_triggerElectron60;
}

bool &Event::triggerElectron24() {
  return m_triggerElectron24;
}
*/
bool &Event::triggerElectron() {
  return m_triggerElectron;
}
bool &Event::triggerMuon() {
  return m_triggerMuon;
}
bool &Event::triggerLargeJet() {
  return m_triggerLargeJet;
}



const std::vector<Electron> &Event::electron() const {
  return m_electron;
}

const std::vector<Muon> &Event::muon() const {
  return m_muon;
}

const std::vector<Jet> &Event::jet() const {
  return m_jet;
}

const std::vector<LargeJet> &Event::largeJet() const {
  return m_largeJet;
}

const std::vector<LargeJet> &Event::largeJetBB() const {
  return m_largeJetBB;
}


std::vector<Particle> &Event::partElectron() {
  return m_partElectron;
}

std::vector<Particle> &Event::partMom() {
  return m_partMom;
}

const std::vector<Particle> &Event::partMom() const {
  return m_partMom;
}

std::vector<Particle> &Event::partMuon() {
  return m_partMuon;
}

std::vector<Jet> &Event::partJet() {
  return m_partJet;
}

std::vector<LargeJet> &Event::partLargeJet() {
  return m_partLargeJet;
}

std::vector<LargeJet> &Event::partLargeJetBB() {
  return m_partLargeJetBB;
}



void Event::partMet(float met_x, float met_y) {
  m_partMet.SetPxPyPzE(met_x, met_y, 0, std::sqrt(std::pow(met_x, 2) + std::pow(met_y, 2)));
}

TLorentzVector Event::partMet() const {
  return m_partMet;
}

const std::vector<Particle> &Event::partElectron() const {
  return m_partElectron;
}

const std::vector<Particle> &Event::partMuon() const {
  return m_partMuon;
}

const std::vector<Jet> &Event::partJet() const {
  return m_partJet;
}

const std::vector<LargeJet> &Event::partLargeJet() const {
  return m_partLargeJet;
}

const std::vector<LargeJet> &Event::partLargeJetBB() const {
  return m_partLargeJetBB;
}


int &Event::runNumber() {
  return m_runNumber;
}

const int Event::runNumber() const {
  return m_runNumber;
}

int &Event::eventNumber() {
  return m_eventNumber;
}

const int Event::eventNumber() const {
  return m_eventNumber;
}

bool &Event::isData() {
  return m_isData;
}
const bool Event::isData() const {
  return m_isData;
}

int &Event::channelNumber() {
  return m_channelNumber;
}
const int Event::channelNumber() const {
  return m_channelNumber;
}

float &Event::mu() {
  return m_mu;
}
const float Event::mu() const {
  return m_mu;
}

float &Event::rho() {
  return m_rho;
}


int &Event::lerr() {
  return m_lerr;
}

int &Event::terr() {
  return m_terr;
}


int &Event::cfl() {
  return m_cfl;
}



bool &Event::atlFastII() {
  return m_atlFastII;
}


int &Event::npv() {
  return m_npv;
}
const int Event::npv() const {
  return m_npv;
}

int &Event::npv_good() {
  return m_npv_good;
}
const int Event::npv_good() const {
  return m_npv_good;
}


int &Event::npv_met() {
  return m_npv_met;
}
const int Event::npv_met() const {
  return m_npv_met;
}

const float Event::rho() const {
  return m_rho;
}

const int Event::lerr() const {
  return m_lerr;
}

const int Event::terr() const {
  return m_terr;
}

const int Event::cfl() const {
  return m_cfl;
}



const bool Event::atlFastII() const {
  return m_atlFastII;
}

const std::vector<double> &Event::weight(const std::string &name) const {
  static const std::vector<double> empty_const;
  std::vector< std::pair<std::string, std::vector<double> > >::const_iterator it = m_weight.begin();
  for (; it != m_weight.end(); ++it) {
    if (it->first == name) break;
  }
  if (it != m_weight.end())
    return it->second;
  return empty_const;
}

std::vector<double> &Event::weight(const std::string &name, bool create) {
  static std::vector<double> empty;
  std::vector< std::pair<std::string, std::vector<double> > >::iterator it = m_weight.begin();
  for (; it != m_weight.end(); ++it) {
    if (it->first == name) break;
  }
  if (it != m_weight.end())
    return it->second;
  if (create) {
    m_weight.push_back(std::make_pair(name, std::vector<double>()));
    return m_weight[m_weight.size()-1].second;
  }
  return empty;
}

float &Event::vxZ() {
  return m_vxZ;
}

const float Event::vxZ() const {
  return m_vxZ;
}

const Event::DataPeriod Event::period() const {
  return m_period;
}

Event::DataPeriod &Event::period() {
  return m_period;
}

unsigned int &Event::lbn() {
  return m_lbn;
}

const unsigned int Event::lbn() const {
  return m_lbn;
}

