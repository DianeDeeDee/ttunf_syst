#include "EventCount.h"

std::map<std::string, TFile *> f = std::map<std::string, TFile *>();
std::map<std::string, TTree *> t = std::map<std::string, TTree *>();

//float getEventCount(int mc_channel_number, const std::string &suff) {
float getEventCount(int mc_channel_number, const std::string &suff, int z) { 
  int type;
  int id;
  float value;
  t[suff]->SetBranchAddress("type", &type);
  t[suff]->SetBranchAddress("id", &id);
  t[suff]->SetBranchAddress("value", &value);
  for (int i = 0; i < t[suff]->GetEntries(); ++i) {
    t[suff]->GetEntry(i);
    //if (type == mc_channel_number && id == -1) {
      if (type == mc_channel_number && id == z) {
        return value;
    }
  }
  return 1;
}

void initEventCount(const std::string &file) {
  f.insert(std::pair<std::string, TFile *>("", new TFile(file.c_str())));
  for (std::map<std::string, TFile *>::const_iterator it = f.begin(); it != f.end(); ++it) {
    t.insert(std::pair<std::string, TTree *>(it->first, (TTree *) it->second->Get("count")));
  }
}
/*void initEventCount() {
  f.insert(std::pair<std::string, TFile *>("", new TFile("EventCount.root")));
  for (std::map<std::string, TFile *>::const_iterator it = f.begin(); it != f.end(); ++it) {
    t.insert(std::pair<std::string, TTree *>(it->first, (TTree *) it->second->Get("count")));
  }
}*/

