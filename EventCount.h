#ifndef EVENTCOUNT_H
#define EVENTCOUNT_H

#include "TFile.h"
#include "TTree.h"
#include <map>
#include <string>
float getEventCount(int mc_channel_number, const std::string &suff = "", int z = -1);
//float getEventCount(int mc_channel_number, const std::string &suff = "");
//void initEventCount();
void initEventCount(const std::string &file = "EventCount.root");
#endif
