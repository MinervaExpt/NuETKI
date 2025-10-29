#ifndef ELECTRONSYSTEMATICS_CPP
#define ELECTRONSYSTEMATICS_CPP

#include "universes/MuonSystematics.h"
#include <iostream>

// Helper functions -- get Weighters, containers of systematics universes
namespace PlotUtils {
  //=================================================================================
  // Minerva electron-momentum-shifted universe 
  //=================================================================================
  template <class T>
  std::vector<T*> GetMinervaElectronSystematics(typename T::config_t chain ) {
    std::vector<T*> ret;

    //ret.push_back(new PlotUtils::MuonUniverseMinerva<T>(chain, -1.));
    //ret.push_back(new PlotUtils::MuonUniverseMinerva<T>(chain, 1.));

    return ret;
  }
}
#endif // ELECTRONSYSTEMATICS_CPP
