#ifndef ELECTRONSYSTEMATICS_H
#define ELECTRONSYSTEMATICS_H

#include "PlotUtils/TreeWrapper.h"
#include "utilities/NSFDefaults.h"
#include "PlotUtils/FluxSystematics.cxx"

// Helper functions declared/defined in the .cxx
// GetMinervaMuonSystematicsMap( typename T::config_t chain );
// GetMinosMuonSystematicsMap( typename T::config_t chain );))

namespace PlotUtils{
  //=================================================================================
  // My electron specific energy shift universe (is this the leakage study one???
  //=================================================================================
  template<class T>
  class ElectronEnergyShiftUniverse: public T
  {
    public:
      ElectronEnergyShiftUniverse(typename T::config_t chw, double nsigma );

    //double GetMuonMomentumShiftMinerva() const;
    //virtual double GetPmuMinerva() const /*override*/;

    virtual std::string ShortName() const /*override*/;
    virtual std::string LatexName() const /*override*/;
    virtual bool IsVerticalOnly()  const  { return false; }/*override*/
  };
}
#endif // ELECTRONSYSTEMATICS_H
