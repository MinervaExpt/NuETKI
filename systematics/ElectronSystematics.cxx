#ifndef MUONSYSTEMATICS_CXX
#define MUONSYSTEMATICS_CXX

#include "universes/MuonSystematics.h"
#include <iostream>

// Helper functions -- get Weighters, containers of systematics universes
namespace PlotUtils {
  //=================================================================================
  // Minerva muon-momentum-shifted universe 
  //=================================================================================
  template <class T>
  std::vector<T*> GetMinervaMuonSystematics(typename T::config_t chain ) {
    std::vector<T*> ret;

    ret.push_back(new PlotUtils::MuonUniverseMinerva<T>(chain, -1.));
    ret.push_back(new PlotUtils::MuonUniverseMinerva<T>(chain, 1.));

    return ret;
  }
}
