#ifndef REGIONS_H
#define REGIONS_H

#include <vector>
#include <map>
#include <string>

class Regions{
    bool m_loadNames;
    /* number of regions */
    int m_nb_regions;
    /* eta boundaries of the regions */
    std::vector<float> m_eta_min;
    std::vector<float> m_eta_max;
    /* phi boundaries of the regions */
    std::vector<float> m_phi_min;
    std::vector<float> m_phi_max;
    std::vector<std::string> m_names;
    //In same cases I need to collect the simple regions in macroRegions according to specific criteria
    bool m_doMacroRegions;
    void m_collectMacroRegionsSL();//Small and large regions are collected together
    void m_collectMacroRegionsSL_UpDn();//Small,Large,Up,Down regions are collected together
    void m_collectMacroRegionsSL_SplitBAR();//Large,Small sectors split plus Feet(12+14) and 11+15 sector split in Barrel
    //This maps the standard regions indexes into the indexes of the macroRegions.
    std::map<int, int> m_macroRegionIdxMap; 
    // The macro regions themselves are stored in a vector with their name
    std::vector<std::string> m_macroRegionName;
    std::vector<double> m_macroRegionInnerEta;//I need this var in few occasions
public:
    void PrintRegions() const;
    unsigned int GetNRegions() const;
    bool verb;
    Regions(std::string inRegionFile, int doMacroRegionsFlag = 0);    
    int GetRegion(const double eta, const double phi) const;
    float GetRegionInnerEta(const int r_i) const;//Return Eta closer to the origin
    std::string GetRegionName(const int r_i) const;
    std::string GetRegionName(const double eta, const double phi) const;  
};

#endif
