#include "MObject.h"
#include "Jet.h"
#include <cmath>

Jet::Jet()
  : MObject() {
  m_type = MObject::jet;
  m_trueflavour = 1;
  m_mv1 = 0;
}

Jet::Jet(const TLorentzVector &v)
  : MObject(v, MObject::jet) {
  m_trueflavour = 1;
  m_mv1 = 0;
}

Jet::~Jet() {
}

const float &Jet::jfit_deltaEta() const {
	  return m_jfit_deltaEta;
	}
	
	float &Jet::jfit_deltaEta() {
  return m_jfit_deltaEta;
	}

const float &Jet::jfit_deltaPhi() const {
	  return m_jfit_deltaPhi;
	}
	float &Jet::jfit_deltaPhi() {
	  return m_jfit_deltaPhi;
	}

int &Jet::trueFlavour() {
  return m_trueflavour;
}
const int Jet::trueFlavour() const {
  return m_trueflavour;
}

bool Jet::pass() const {
  if (std::fabs(mom().Eta()) > 2.5) return false;
  if (mom().Perp() < 25e3) return false;
  if ( (mom().Perp() <= 50e3) && (std::fabs(mom().Eta()) < 2.4) && \
       (std::fabs(jvf()) <= 0.5) ) //jvf=jet vertex fraction: By combining tracks and their primary vertices with calorimeter jets we define a discriminant, the Jet Vertex Fraction (or JVF) which measures the probability that a jet originated from a particular vertex = it's Pt weighted fraction of tracks in jet that are compatible with primary vertex
    return false;
  return true;
}

bool Jet::btag() const {
  if (mv1() < 0.7892) return false;
  return true;
}

const float Jet::mv1() const {
  return m_mv1;
}
float &Jet::mv1() {
  return m_mv1;
}

const float Jet::jvf() const {
  return m_jvf;
}
float &Jet::jvf() {
  return m_jvf;
}


double &Jet::detE() {
  return m_detE;
}

double &Jet::detEta() {
  return m_detEta;
}
double &Jet::detPhi() {
  return m_detPhi;
}
double &Jet::detM() {
  return m_detM;
}
double &Jet::Ax() {
  return m_Ax;
}
double &Jet::Ay() {
  return m_Ay;
}
double &Jet::Az() {
  return m_Az;
}

double &Jet::Ae() {
  return m_Ae;
}


const float Jet::sumPtTrk() const {
  return m_sumPtTrk;
}
float &Jet::sumPtTrk() {
  return m_sumPtTrk;
}

const std::vector<int> &Jet::TrackAssoc_index() const {
  return m_TrackAssoc_index;
}
std::vector<int> &Jet::TrackAssoc_index() {
  return m_TrackAssoc_index;
}

const bool Jet::isBadLoose() const {
  return m_isBadLoose;
}
bool &Jet::isBadLoose() {
  return m_isBadLoose;
}

const bool Jet::isBadLooseMinus() const {
	  return m_isBadLooseMinus;
	}
	bool &Jet::isBadLooseMinus() {
	  return m_isBadLooseMinus;
	}

bool Jet::passBadLoose() const {
  bool x = true;
  if (mom().Perp() > 20e3 && std::fabs(mom().Eta()) < 2.5) { //(mom().Perp() > 25e3 for G*
    x = !isBadLoose();
  }
  return x;
}
bool Jet::passBadLooseMinus() const {
	  bool x = true;
	  if (mom().E() > 0 && mom().Perp() > 20e3 && std::fabs(mom().Eta()) < 2.5) {
	    x = !isBadLooseMinus();
	  }
	  return x;
	}
