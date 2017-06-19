#include "MObject.h"
#include "Particle.h"
#include <cmath>

Particle::Particle()
  : MObject(), m_isol(0), m_pdgId(11)  {
  m_type = MObject::part;
}

Particle::Particle(const TLorentzVector &v)
  : MObject(v, MObject::part), m_isol(0), m_pdgId(11) {
}

Particle::~Particle() {
}

void Particle::setIsol(float iso) {
  m_isol = iso;
}

float Particle::isol() const {
  return m_isol;
}

void Particle::setId(int pdgId) {
  m_pdgId = pdgId;
}

int Particle::id() const {
  return m_pdgId;
}


bool Particle::pass() const {
  if (mom().Perp() < 25e3) return false;
  float eta = std::fabs(mom().Eta());
  if (std::fabs(id()) == 11) {
    if (eta > 1.37 && eta < 1.52) return false;
    if (eta > 2.47) return false;
  } else if (std::fabs(id()) == 13) {
    if (eta > 2.5) return false;
  }
  if (isol()/mom().Perp() > 0.05) return false;
  return true;
}

