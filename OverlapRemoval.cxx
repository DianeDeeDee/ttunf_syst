#include "TopD3PDSelection/OverlapRemoval.h"
#include "TopD3PDSelection/KinematicUtils.h"
#include <algorithm>
#include <set>

// just for debugging.
#include <iostream>

OverlapRemoval::OverlapRemoval():
  CachedSelection(),
  m_jets(),
  m_muons(),
  m_melectrons(),
  m_jelectrons()
{ }

int OverlapRemoval::classify(const EventInfo* /*eventInfo*/,
                             const ElectronContainer *electrons, 
                             const MuonContainer *muons, 
                             const JetContainer *jets,
                             ElectronSelection *electronSelection, 
                             MuonSelection *muonSelection, 
                             JetSelection *jetSelection,
                             const PhotonContainer *photons,
                             PhotonSelection *photonSelection)
{
  //if(CachedSelection::useCache(eventInfo)) return 0;

  int ret_val;

  m_jets.clear();
  m_muons.clear();
  m_melectrons.clear();
  m_jelectrons.clear();

  if((ret_val = classifyMuons(jets, muons, jetSelection, muonSelection)) != 0) return ret_val;
  if((ret_val = classifyJets(electrons, jets, electronSelection, jetSelection, photons, photonSelection)) != 0) return ret_val;
  if((ret_val = classifyElectrons(muons, jets, electrons, muonSelection, jetSelection, electronSelection)) != 0) return ret_val; 

  return 0;
}

int OverlapRemoval::classify(const EventInfo* /*eventInfo*/,
                             const ElectronContainer *electrons,
                             const MuonContainer *muons,
                             const JetContainer *jets,
                             const TauContainer *taus,
                             ElectronSelection *electronSelection,
                             MuonSelection *muonSelection,
                             JetSelection *jetSelection,
                             TauSelection *tauSelection,
                             const PhotonContainer* photons,
                             PhotonSelection* photonSelection)
{
   //if(CachedSelection::useCache(eventInfo)) return 0;

   int ret_val;

   m_jets.clear();
   m_muons.clear();
   m_melectrons.clear();
   m_jelectrons.clear();
   m_taus.clear();
   m_electrons.clear();

   if((ret_val = classifyMuons(jets, muons, jetSelection, muonSelection)) != 0) return ret_val;

   if((ret_val = classifyJets(electrons, jets, electronSelection, jetSelection, photons, photonSelection)) != 0) return ret_val;
   //remove tau in el
   if((ret_val = classifyTauElectron(taus,electrons,tauSelection,electronSelection)) != 0) return ret_val;
   //remove tau in mu
   if((ret_val = classifyTauMuon(taus,muons,tauSelection,muonSelection)) != 0) return ret_val;

   if((ret_val = classifyElectrons(muons, jets, electrons, muonSelection, jetSelection, electronSelection)) != 0) return ret_val; 

   //remove jet in tau 
   if((ret_val = classifyTauJet(taus,jets,tauSelection,jetSelection)) != 0) return ret_val;
   
   return 0;
}

int OverlapRemoval::classifyMuons(const JetContainer *jets,
                                  const MuonContainer *muons,
                                  JetSelection*,
                                  MuonSelection *muonSelection)
{
  // Loop over the muon indices and remove any muons which fail
  // the overlap removal.  This class is a friend of the MuonSelection
  // class to update the list of indices.
  std::set<unsigned int> idxSet;
  // push the indices onto the set (so we only do each index once).
  //idxSet.insert(muonSelection->m_looseMuons.begin(), muonSelection->m_looseMuons.end());
  idxSet.insert(muonSelection->m_goodMuons.begin(), muonSelection->m_goodMuons.end());
  std::set<unsigned int>::iterator itrMuonIndex = idxSet.begin();
  std::set<unsigned int>::iterator itrMuonIndexEnd = idxSet.end();
  for(;itrMuonIndex!=itrMuonIndexEnd;++itrMuonIndex) {
    double muonEta = muons->at(*itrMuonIndex)->eta();
    double muonPhi = muons->at(*itrMuonIndex)->phi();
    
    // Loop over all jets above the pT cut until a delta R value of
    // less than the overlap cut is found or there are no more jets.
    double dR = 1.0;
    JetContainer::const_iterator itrJet = jets->begin();
    JetContainer::const_iterator itrJetEnd = jets->end();
    while(dR >= 0.4 && itrJet != itrJetEnd) {
      
      // Select jets that pass the overlap pT cut (for 17 increased from 20-25 GeV) and the JVF cut
//      if(!(((*itrJet)->correctedPt()/1.0e3) > 25.0 && (fabs((*itrJet)->jvtxf())) > 0.50)) { itrJet++; continue; }      
       
       // new JVF approach 
        if (!((*itrJet)->correctedPt() > 25000 && // to pass jet is above 25 GeV and pass one of the following (or of any of the 3)
              (
                 ((*itrJet)->correctedPt() > 50000 ) ||  // if pt > 50 GeV no JVF cut applied
                 (std::fabs((*itrJet)->correctedEta()) > 2.4) ||  // if eta > 2.4 no JVF cut applied
                 (std::fabs((*itrJet)->jvtxf()) > (*itrJet)->jvfcut()) // JVF cut is passed
                 )))
        { itrJet++; continue; }   
      
      // Select jets that pass  E>=0 
      if((*itrJet)->m_orig.E() < 0.) { itrJet++; continue; }
      
      // Calculate delta R between the jet and the muon
      dR = KinematicUtils::deltaR(muonEta, (*itrJet)->correctedEta(), muonPhi, (*itrJet)->correctedPhi());

      itrJet++;
    }

    // If the muon fails the delta R requirement remove it from the
    // list of good muons.
    if(dR < 0.4) { 
      m_muons.push_back(*itrMuonIndex);
    }
  }

  // Remove the overlapping indices.
//   removeOverlap(&m_muons, &(muonSelection->m_looseMuons));
  removeOverlap(&m_muons, &(muonSelection->m_goodMuons));

  return 0;
}

int OverlapRemoval::classifyJets(const ElectronContainer *electrons,
                                 const JetContainer *jets,
                                 ElectronSelection *electronSelection,
                                 JetSelection *jetSelection,
                                 const PhotonContainer *photons,
                                 PhotonSelection *photonSelection)
{
  // (1) Loop over all good electrons
  if (electronSelection) {
    std::vector<unsigned int>::iterator itrElec = electronSelection->m_goodElectrons.begin();
    std::vector<unsigned int>::iterator itrElecEnd = electronSelection->m_goodElectrons.end();
    for(; itrElec != itrElecEnd; ++itrElec) {
      double eleEta = electrons->at(*itrElec)->tracketa();
      double elePhi = electrons->at(*itrElec)->trackphi();
    
      int dRminIndex = -1;
      unsigned int ijet = 0;
      double dR = 0., dRmin = 0.;
      JetContainer::const_iterator itrJet = jets->begin();
      JetContainer::const_iterator itrJetEnd = jets->end();
      while(itrJet != itrJetEnd) {
        dR = KinematicUtils::deltaR(eleEta, (*itrJet)->correctedEta(),
                                    elePhi, (*itrJet)->correctedPhi());
      
        if(dRminIndex == -1 || dR<dRmin) { 
          dRmin = dR;
          dRminIndex = ijet;
        }
      
        itrJet++;
        ijet++;
      }

      // If the jet with the smallest delta R with respect to the
      // electron is within the delta R cut assume this is the electron
      // and store the index of this jet.
      if(dRmin < 0.2 && dRminIndex >= 0) { 
        m_jets.push_back(dRminIndex);
      }
    }
  }

  // (2) Loop over all good photons (in case there are)
  if (photonSelection!=0) {
    std::vector<unsigned int>::iterator itrPhot    = photonSelection->m_goodPhotons.begin();
    std::vector<unsigned int>::iterator itrPhotEnd = photonSelection->m_goodPhotons.end();
    for(; itrPhot!=itrPhotEnd; ++itrPhot) {
      double phEta = photons->at(*itrPhot)->cl_eta();
      double phPhi = photons->at(*itrPhot)->cl_phi();
      
      int dRminIndex = -1;
      unsigned int ijet = 0;
      double dRmin = 0.;      
      JetContainer::const_iterator itrJet    = jets->begin();
      JetContainer::const_iterator itrJetEnd = jets->end();
      while(itrJet != itrJetEnd) {
        double dR = KinematicUtils::deltaR(phEta, (*itrJet)->correctedEta(),
                                           phPhi, (*itrJet)->correctedPhi());	
        if(dRminIndex == -1 || dR<dRmin) { 
          dRmin = dR;
          dRminIndex = ijet;
        }	
        itrJet++;
        ijet++;
      }      
      if(dRmin<0.1 && dRminIndex>=0)
        m_jets.push_back(dRminIndex);
    }
    std::sort(m_jets.begin(), m_jets.end()); // sort vector
    std::vector<unsigned int>::iterator it = std::unique(m_jets.begin(),m_jets.end()); // remove duplicates
    m_jets.resize(it-m_jets.begin()); // resize it
  }
    
  // If at least one jet is found that matches an electron, check if
  // this is a good jet or a tagged jet.
  if(m_jets.size() > 0) {
    removeOverlap(&m_jets, &(jetSelection->m_goodJets));
    removeOverlap(&m_jets, &(jetSelection->m_noJvfJets));
    removeOverlap(&m_jets, &(jetSelection->m_taggedJets));
    removeOverlap(&m_jets, &(jetSelection->m_badJets)); // Jet cleaning after overlap removal
  }

  return 0;
}

// overlap removal between electrons and muons
// remove the event if the 'good' electron and 'good' muon share a track
// where 'good' is passing all final lepton cuts except the jet-mu dR cut
// here pick a matching in phi and theta of 0.005 of the ID track
// Also required to look for jets (after overlap removal) within .4 of
// the electrons
int OverlapRemoval::classifyElectrons(const MuonContainer *muons,
				      const JetContainer *jets,
                                      const ElectronContainer *electrons,
                                      MuonSelection* muonSelection,
				      JetSelection *jetSelection,
                                      ElectronSelection *electronSelection)
{
  if (!electronSelection) return 0;

  // Loop over all good electrons
  std::vector<unsigned int>::iterator itrElecIndex = electronSelection->m_goodElectrons.begin();
  std::vector<unsigned int>::iterator itrElecIndexEnd = electronSelection->m_goodElectrons.end();
  for(;itrElecIndex!=itrElecIndexEnd; ++itrElecIndex) {
    MuonContainer::const_iterator itrMuon = muons->begin();
    MuonContainer::const_iterator itrMuonEnd = muons->end();
    bool found = false;
    // loop over all good muons
    while(!found && itrMuon != itrMuonEnd) {
      const Muon* muon = (*itrMuon); 
      //if ( !muonSelection->classifyGoodMuon(muon) ) { itrMuon++; continue; }
      if ( !muonSelection->classifyGoodMuon((*itrMuon)) ) { itrMuon++; continue; }
      // do they share the same track?
      if(std::fabs(electrons->at(*itrElecIndex)->trackphi() - muon->id_phi())<0.005 
         && std::fabs(electrons->at(*itrElecIndex)->tracktheta() - muon->id_theta())<0.005) {
        found = true;
        continue;
      }

      itrMuon++;
    }
    if (found) m_melectrons.push_back(*itrElecIndex); // Save this electron index

    std::vector<unsigned int>::iterator itrJetIndex = jetSelection->m_goodJets.begin();
    std::vector<unsigned int>::iterator itrJetIndexEnd = jetSelection->m_goodJets.end();
    for(;itrJetIndex!=itrJetIndexEnd; ++itrJetIndex) {
      if(KinematicUtils::deltaR(electrons->at(*itrElecIndex)->tracketa(),jets->at(*itrJetIndex)->correctedEta(),
				electrons->at(*itrElecIndex)->trackphi(),jets->at(*itrJetIndex)->correctedPhi()) < 0.4) {
	m_jelectrons.push_back(*itrElecIndex);
      }
    }
  }

  // remove only the jet overlaps (muon-ele events get rejected during
  // the cut flow)
  removeOverlap(&m_jelectrons, &(electronSelection->m_goodElectrons));
  return 0;
}

void OverlapRemoval::removeOverlap(std::vector<unsigned int>* overlapIndices, 
                                   std::vector<unsigned int>* selectionIndices)
{
  std::vector<unsigned int>::iterator itrSelectionIndex = selectionIndices->begin();
  std::vector<unsigned int>::iterator itrSelectionIndexEnd = selectionIndices->end();
  std::vector<unsigned int> remainingIndices;
  for(;itrSelectionIndex!=itrSelectionIndexEnd;++itrSelectionIndex) {
    if(std::find(overlapIndices->begin(), overlapIndices->end(), (*itrSelectionIndex)) == overlapIndices->end()) {
      remainingIndices.push_back(*itrSelectionIndex);
    }
  }
  *selectionIndices = remainingIndices; // Overwrite the original vector.
}

int OverlapRemoval::classifyElectronMuon(const MuonContainer *muons, const ElectronContainer *electrons,
                                         MuonSelection* /*muonSelection*/, ElectronSelection *electronSelection){
   // do not use the full muon definition here

   // Loop over all good electrons
   std::vector<unsigned int>::iterator itrElecIndex = electronSelection->m_goodElectrons.begin();
   std::vector<unsigned int>::iterator itrElecIndexEnd = electronSelection->m_goodElectrons.end();
   for(;itrElecIndex!=itrElecIndexEnd; ++itrElecIndex) {

      // The overlap removal does not specify that the muons should pass
      // the track quality cuts.  Therefore, loop over all muons in this event.
      MuonContainer::const_iterator itrMuon = muons->begin();
      MuonContainer::const_iterator itrMuonEnd = muons->end();
      bool found = false;
      while(!found && itrMuon != itrMuonEnd) {
         const Muon* muon = (*itrMuon);
         if(!muon->tight()) { itrMuon++; continue; } // tight
         if(muon->author() != 12) { itrMuon++; continue; } // MuonParameters::MuidCo
         if(muon->pt()/1.0e3 <= 15.0) { itrMuon++; continue; }
         if(fabs(muon->eta()) >= 2.5) { itrMuon++; continue; }
         // if((muon->ptcone30()/1.e3 >= 2.5) || (muon->etcone20()/1.e3 >= 4.0)) { itrMuon++; continue; }// isolation
         if(muon->miniIso10_4()/muon->correctedPt() < 0.05) { itrMuon++; continue; }// isolation


         if(std::fabs(electrons->at(*itrElecIndex)->trackphi() - muon->id_phi())<0.005
            && std::fabs(electrons->at(*itrElecIndex)->tracktheta() - muon->id_theta())<0.005) {
            found = true;
            continue;
         }

         itrMuon++;
      }

      if(found) m_electrons.push_back(*itrElecIndex); // Save this electron index
   }

   return 0;
}

int OverlapRemoval::classifyTauElectron(const TauContainer* taus , const ElectronContainer *electrons,
                                        TauSelection* tauSelection, ElectronSelection *electronSelection){

   // Loop over the good tau indices and remove any tau which fail
   // the overlap removal.  This class is a friend of the TauSelection
   // class to update the list of indices.

   std::vector<unsigned int> tausRemaining; // new vector of tau indices
   std::vector<unsigned int>::iterator itrTauIndex = tauSelection->m_goodTaus.begin();
   std::vector<unsigned int>::iterator itrTauIndexEnd = tauSelection->m_goodTaus.end();

   for(;itrTauIndex!=itrTauIndexEnd;++itrTauIndex) {
      double tauEta = taus->at(*itrTauIndex)->eta();
      double tauPhi = taus->at(*itrTauIndex)->phi();

      // Loop over all electrons above the pT cut until a delta R value of
      // less than the overlap cut is found or there are no more electron.
      double dR = 1.0;
      std::vector<unsigned int>::iterator itrElectronIndex = electronSelection->m_goodElectrons.begin();
      std::vector<unsigned int>::iterator itrElectronIndexEnd = electronSelection->m_goodElectrons.end();
      while(dR >= 0.2 && itrElectronIndex != itrElectronIndexEnd) {
         // Calculate delta R between the jet and the muon
         dR = KinematicUtils::deltaR(tauEta, electrons->at(*itrElectronIndex)->tracketa(), tauPhi, electrons->at(*itrElectronIndex)->trackphi());
         itrElectronIndex++;
      }

      // If the tau fails the delta R requirement remove it from the
      // list of good taus.
      if(dR < 0.2) {
         m_taus.push_back(*itrTauIndex);
      }
      else {
         tausRemaining.push_back(*itrTauIndex);
      }
   }

   tauSelection->m_goodTaus = tausRemaining;

   return 0;

}

int OverlapRemoval::classifyTauMuon(const TauContainer* taus, const MuonContainer *muons,
                                    TauSelection* tauSelection, MuonSelection *muonSelection){

   // Loop over the good tau indices and remove any tau which fail
   // the overlap removal.  This class is a friend of the TauSelection
   // class to update the list of indices.

   std::vector<unsigned int> tausRemaining; // new vector of tau indices
   std::vector<unsigned int>::iterator itrTauIndex = tauSelection->m_goodTaus.begin();
   std::vector<unsigned int>::iterator itrTauIndexEnd = tauSelection->m_goodTaus.end();

   for(;itrTauIndex!=itrTauIndexEnd;++itrTauIndex) {
      double tauEta = taus->at(*itrTauIndex)->eta();
      double tauPhi = taus->at(*itrTauIndex)->phi();

      // Loop over all muons above the pT cut until a delta R value of
      // less than the overlap cut is found or there are no more muons.
      double dR = 1.0;
      std::vector<unsigned int>::iterator itrMuonIndex = muonSelection->m_goodMuons.begin();
      std::vector<unsigned int>::iterator itrMuonIndexEnd = muonSelection->m_goodMuons.end();
      while(dR >= 0.2 && itrMuonIndex != itrMuonIndexEnd) {
         // Calculate delta R between the jet and the muon
         dR = KinematicUtils::deltaR(tauEta, muons->at(*itrMuonIndex)->eta(), tauPhi, muons->at(*itrMuonIndex)->phi());
         itrMuonIndex++;
      }

      // If the tau fails the delta R requirement remove it from the
      // list of good taus.
      if(dR < 0.2) {
         m_taus.push_back(*itrTauIndex);
      }
      else {
         tausRemaining.push_back(*itrTauIndex);
      }
   }

   tauSelection->m_goodTaus = tausRemaining;

   return 0;

}

int OverlapRemoval::classifyTauJet(const TauContainer* taus, const JetContainer *jets,
                                    TauSelection* tauSelection, JetSelection *jetSelection){
   // Loop over the good jet indices and remove any jet which fail
   // the overlap removal.

  std::vector<unsigned int> jetsRemaining; // new vector of jet indices 
  std::vector<unsigned int>::iterator itrJetIndex = jetSelection->m_goodJets.begin();
  std::vector<unsigned int>::iterator itrJetIndexEnd = jetSelection->m_goodJets.end();
  
   for(;itrJetIndex!=itrJetIndexEnd;++itrJetIndex) {

     double jetEta = jets->at(*itrJetIndex)->correctedEta();
     double jetPhi = jets->at(*itrJetIndex)->correctedPhi();

      // Loop over all taus until a delta R value of
      // less than the overlap cut is found or there are no more taus.
      double dR = 1.0;

      std::vector<unsigned int>::iterator itrTauIndex = tauSelection->m_goodTaus.begin();
      std::vector<unsigned int>::iterator itrTauIndexEnd = tauSelection->m_goodTaus.end();
      
      while(dR >= 0.2 && itrTauIndex != itrTauIndexEnd) {
         // Calculate delta R between the tau and the jet
         dR = KinematicUtils::deltaR(taus->at(*itrTauIndex)->eta(), jetEta, taus->at(*itrTauIndex)->phi(), jetPhi);
         itrTauIndex++;
      }

      // If the tau fails the delta R requirement remove it from the
      // list of good taus.
      if(dR < 0.2) {
         m_jets.push_back(*itrJetIndex);
      }
      else {
         jetsRemaining.push_back(*itrJetIndex);
      }
   }

   jetSelection->m_goodJets = jetsRemaining;


   return 0;
}

int OverlapRemoval::classifyTauBjet(const TauContainer* taus, const JetContainer *jets,
                                    TauSelection* tauSelection, JetSelection *jetSelection){
   // Loop over the good tau indices and remove any tau which fail
   // the overlap removal.  This class is a friend of the TauSelection
   // class to update the list of indices.

   std::vector<unsigned int> tausRemaining; // new vector of muon indices
   std::vector<unsigned int>::iterator itrTauIndex = tauSelection->m_goodTaus.begin();
   std::vector<unsigned int>::iterator itrTauIndexEnd = tauSelection->m_goodTaus.end();

   for(;itrTauIndex!=itrTauIndexEnd;++itrTauIndex) {
      double tauEta = taus->at(*itrTauIndex)->eta();
      double tauPhi = taus->at(*itrTauIndex)->phi();

      // Loop over all muons above the pT cut until a delta R value of
      // less than the overlap cut is found or there are no more muons.
      double dR = 1.0;
      std::vector<unsigned int>::iterator itrBjetIndex = jetSelection->m_taggedJets.begin();
      std::vector<unsigned int>::iterator itrBjetIndexEnd = jetSelection->m_taggedJets.end();
      while(dR >= 0.4 && itrBjetIndex != itrBjetIndexEnd) {
         // Calculate delta R between the jet and the muon
         dR = KinematicUtils::deltaR(tauEta, jets->at(*itrBjetIndex)->correctedEta(), tauPhi, jets->at(*itrBjetIndex)->correctedPhi());
         itrBjetIndex++;
      }

      // If the tau fails the delta R requirement remove it from the
      // list of good taus.
      if(dR < 0.4) {
         m_taus.push_back(*itrTauIndex);
      }
      else {
         tausRemaining.push_back(*itrTauIndex);
      }
   }

   tauSelection->m_goodTaus = tausRemaining;


   return 0;
}

