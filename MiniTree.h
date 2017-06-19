#ifndef MINITREE_H
#define MINITREE_H

#include <string>
#include "TChain.h"
#include "TTree.h"
#include <vector>
#include "Event.h"

#include "TFile.h"

class MiniTree {
  public:

    MiniTree(bool toWrite = true, const std::string &file = "tree.root");
    virtual ~MiniTree();

    void addFileToRead(const std::string &fname);
    int GetEntries();
    double getSumWeights();
    std::vector<float> getPassCuts();
    void read(int event, Event &e, int &systIdx);
    void write(const Event &e, const int _systIdx);

    std::vector<std::string> *weightNames();
    std::vector<std::string> *systematicsNames();
    double &sumWeights();
    std::vector<double> &sumWeights_var();
    std::vector<float> &passCuts();


  private:
    double m_sumWeights;
    std::vector<double> *m_sumWeights_var;
    std::vector<float> *m_passCuts;

    void prepareBranches();

    TFile  *m_fileToWrite;
    TTree *m_chain;
    TTree *m_num;
    bool m_toWrite;

    std::vector<std::string> *m_sfs;
    std::vector<std::string> *m_systs;

    int mc_channel_number;

    int                   isTight;

    std::vector<std::vector<double> >   *weight;

    int                   systIdx;

    bool                  passReco;
    bool                  passPart;

    bool                  triggerElectron;
    bool                  triggerElectron24;
    bool                  triggerElectron60;
    bool                  triggerMuon;
    bool                  triggerLargeJet;

    int                   el_n;
    std::vector<float>   *el_pt;
    std::vector<float>   *el_eta;
    std::vector<float>   *el_phi;
    std::vector<float>   *el_E;

    int                   mu_n;
    std::vector<float>   *mu_pt;
    std::vector<float>   *mu_eta;
    std::vector<float>   *mu_phi;
    std::vector<float>   *mu_E;

    int                   jet_n;
    std::vector<float>   *jet_pt;
    std::vector<float>   *jet_eta;
    std::vector<float>   *jet_phi;
    std::vector<float>   *jet_E;
    std::vector<float>   *jet_mv1;
    std::vector<int>     *jet_trueflav;

    float   trigjet_pt;
    float   trigjet_eta;
    float   trigjet_phi;
    float   trigjet_E;

    int                   ljet_n;
    std::vector<float>   *ljet_pt;
    std::vector<float>   *ljet_eta;
    std::vector<float>   *ljet_phi;
    std::vector<float>   *ljet_E;
    std::vector<float>   *ljet_split12;
    std::vector<int>     *ljet_trueflav;
    std::vector<float>   *ljet_trimmed_tau1;
    std::vector<float>   *ljet_trimmed_tau2;
    std::vector<float>   *ljet_trimmed_am_tau1;
    std::vector<float>   *ljet_trimmed_am_tau2;
    std::vector<int>     *ljet_sub_n;
    std::vector<std::vector<float> > *ljet_sub_pt;
    std::vector<std::vector<float> > *ljet_sub_eta;
    std::vector<std::vector<float> > *ljet_sub_phi;
    std::vector<std::vector<float> > *ljet_sub_E;
    std::vector<std::vector<float> > *ljet_sub_lcpt;
    std::vector<int>     *ljet_subunc_n;
    std::vector<std::vector<float> > *ljet_subunc_pt;
    std::vector<std::vector<float> > *ljet_subunc_eta;
    std::vector<std::vector<float> > *ljet_subunc_phi;
    std::vector<std::vector<float> > *ljet_subunc_E;
    std::vector<std::vector<float> > *ljet_subunc_lcpt;
    std::vector<float>   *ljet_htt;
 
    int                   ljetBB_n;
    std::vector<float>   *ljetBB_pt;
    std::vector<float>   *ljetBB_eta;
    std::vector<float>   *ljetBB_phi;
    std::vector<float>   *ljetBB_E;
    std::vector<float>   *ljetBB_split12;
    std::vector<int>     *ljetBB_sub_n;
    std::vector<std::vector<float> > *ljetBB_sub_pt;
    std::vector<std::vector<float> > *ljetBB_sub_eta;
    std::vector<std::vector<float> > *ljetBB_sub_phi;
    std::vector<std::vector<float> > *ljetBB_sub_E;
    std::vector<std::vector<float> > *ljetBB_sub_lcpt;
    std::vector<int>     *ljetBB_subunc_n;
    std::vector<std::vector<float> > *ljetBB_subunc_pt;
    std::vector<std::vector<float> > *ljetBB_subunc_eta;
    std::vector<std::vector<float> > *ljetBB_subunc_phi;
    std::vector<std::vector<float> > *ljetBB_subunc_E;
    std::vector<std::vector<float> > *ljetBB_subunc_lcpt;
    std::vector<float>   *ljetBB_ug_tau1;
    std::vector<float>   *ljetBB_ug_tau2;
    std::vector<float>   *ljetBB_ug_am_tau1;
    std::vector<float>   *ljetBB_ug_am_tau2;
    std::vector<float>   *ljetBB_bdrs_m;
    std::vector<float>   *ljetBB_bdrs_pt;
    std::vector<float>   *ljetBB_bdrs_phi;
    std::vector<float>   *ljetBB_bdrs_eta;
    std::vector<float>   *ljetBB_ug_bdrs_tau1;
    std::vector<float>   *ljetBB_ug_bdrs_tau2;
    std::vector<float>   *ljetBB_htt;
   
    int hfor;


    float met_etx;
    float met_ety;

    int                   partel_n;
    std::vector<float>   *partel_pt;
    std::vector<float>   *partel_eta;
    std::vector<float>   *partel_phi;
    std::vector<float>   *partel_E;

    int                   partmu_n;
    std::vector<float>   *partmu_pt;
    std::vector<float>   *partmu_eta;
    std::vector<float>   *partmu_phi;
    std::vector<float>   *partmu_E;

    int                   partjet_n;
    std::vector<float>   *partjet_pt;
    std::vector<float>   *partjet_eta;
    std::vector<float>   *partjet_phi;
    std::vector<float>   *partjet_E;
    std::vector<int>     *partjet_trueflav;

    int                   partmom_n;
    std::vector<float>   *partmom_pt;
    std::vector<float>   *partmom_eta;
    std::vector<float>   *partmom_phi;
    std::vector<float>   *partmom_E;
    std::vector<int>     *partmom_trueflav;

    int                   partljet_n;
    std::vector<float>   *partljet_pt;
    std::vector<float>   *partljet_eta;
    std::vector<float>   *partljet_phi;
    std::vector<float>   *partljet_E;
    std::vector<float>   *partljet_split12;
    std::vector<int>     *partljet_trueflav;
    std::vector<int>     *partljet_sub_n;
    std::vector<std::vector<float> > *partljet_sub_pt;
    std::vector<std::vector<float> > *partljet_sub_eta;
    std::vector<std::vector<float> > *partljet_sub_phi;
    std::vector<std::vector<float> > *partljet_sub_E;

    int                   partljetBB_n;
    std::vector<float>   *partljetBB_pt;
    std::vector<float>   *partljetBB_eta;
    std::vector<float>   *partljetBB_phi;
    std::vector<float>   *partljetBB_E;
    std::vector<float>   *partljetBB_split12;
    std::vector<int>     *partljetBB_sub_n;
    std::vector<std::vector<float> > *partljetBB_sub_pt;
    std::vector<std::vector<float> > *partljetBB_sub_eta;
    std::vector<std::vector<float> > *partljetBB_sub_phi;
    std::vector<std::vector<float> > *partljetBB_sub_E;



    float partmet_etx;
    float partmet_ety;

    float avmu;
    float npv;
};

#endif

