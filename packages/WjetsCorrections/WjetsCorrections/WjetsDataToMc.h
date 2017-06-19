#ifndef WJETSDATATOMC_H
#define WJETSDATATOMC_H

class WjetsDataToMc {
 public:
  /*
  * jetMultiplicity is defined by jetMult enum.
  * isElec is true for electrons and false for muons.
  * isTagged is truth for tagged and false for pretag.
  * shift is 0 for the nominal values, -1 for -ve shift and +1 for positive shift.
  * returns the scale factor documented at:
  *  https://twiki.cern.ch/twiki/bin/view/AtlasProtected/WplusJetsBackgroundsforTopAnalyses
  */
  float correction(unsigned int jetMultiplicity, bool isElec, bool isTagged, int shift);

  enum jetMult {
    NJET_1,
    NJET_2,
    NJET_3,
    NJET_4,
    NJET_4_INC,
    NJET_5_INC,
    NJET_BINS};

 private:
  const static float m_elecPretagNominal[NJET_BINS];
  const static float m_elecPretagUncertainty[NJET_BINS];
  const static float m_elecTaggedNominal[NJET_BINS];
  const static float m_elecTaggedUncertainty[NJET_BINS];
  const static float m_muonPretagNominal[NJET_BINS];
  const static float m_muonPretagUncertainty[NJET_BINS];
  const static float m_muonTaggedNominal[NJET_BINS];
  const static float m_muonTaggedUncertainty[NJET_BINS];
};

#endif
