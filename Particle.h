#ifndef PARTICLE_H
#define PARTICLE_H

#include "MObject.h"
#include "TLorentzVector.h"

class Particle : public MObject {
  public:
    Particle();
    Particle(const TLorentzVector &v);
    virtual ~Particle();

    void setIsol(float iso);
    float isol() const;

    void setId(int pdgId);
    int id() const;

    bool pass() const;
  protected:
    float m_isol;
    int m_pdgId;
};

#endif
