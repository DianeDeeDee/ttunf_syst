#include "MObject.h"

MObject::MObject()
  : m_mom(0,0,0,0), m_mom_corr(0,0,0,0), m_type(MObject::jet) {
}

MObject::MObject(const TLorentzVector &v, const Type t)
  :m_mom(v), m_hasCorr(false), m_mom_corr(v), m_type(t) {
}

MObject::~MObject() {
}

MObject::Type &MObject::type() {
  return m_type;
}

const MObject::Type &MObject::type() const {
  return m_type;
}

void MObject::setCorrection(const TLorentzVector &v) {
  m_hasCorr = true;
  m_mom_corr = v;
}

const TLorentzVector &MObject::mom(bool original) const {
  if (original || !m_hasCorr) return m_mom;
  return m_mom_corr;
}

TLorentzVector &MObject::mom(bool original) {
  if (original || !m_hasCorr) return m_mom;
  return m_mom_corr;
}

// must be in the header file
/*
template <class C>
float MObject::minDeltaR(const std::vector<C> &o) const {
  float dR = 99;
  for (int k = 0; k < o.size(); ++k) {
    float tdR = mom().DeltaR(o[k].mom());
    if (tdR < dR) {
      dR = tdR;
    }
  }
  return dR;
}
*/
TLorentzVector &MObject::corr(const std::string &syst, bool create) {
  if (create)
    return m_corrMap[syst];

  if (m_corrMap.find(syst) == m_corrMap.end())
  //  return m_mom_corr;
    return mom();
  return m_corrMap[syst];
}

TLorentzVector MObject::corr(const std::string &syst) const {
  std::map<std::string, TLorentzVector>::const_iterator it = m_corrMap.find(syst);
  if (it == m_corrMap.end())
//    return m_mom_corr;
return mom();
  return it->second;
}

