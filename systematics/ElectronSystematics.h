#ifndef ELECTRONSYSTEMATICS_H
#define ELECTRONSYSTEMATICS_H

#include "event/CVUniverse.h"

#include <map>
#include <string>
#include <vector>

namespace PlotUtils{

  //=================================================================================
  // My electron specific energy shift universe (shifts electron energy up/down by the ECAL/HCAL energy * 0.015/0.05 respectively
  //=================================================================================
  template<class T>
  class ElectronEnergyShiftUniverse: public T
  {
    private: 
      std::string region;

    public:
    ElectronEnergyShiftUniverse(typename T::config_t chw, double nsigma, const std::string& region_) : T(chw, nsigma), region(region_) {}

    //Overrides CVUnivers::GetEMEnergyShift()
    double GetEMEnergyShift() const override {
      static const double EM_ENERGY_SCALE_SHIFT_ECAL = -0.058;    
      static const std::map<std::string, double > EM_ENERGY_SCALE_UNCERTAINTY = {
	{ "ECAL", 0.015 },
	{ "HCAL", 0.05 }
      };

      std::string branch_name = std::string("prong_") + region + "CalibE";
      //std::cout << "region = " << region << ", branch name = " << branch_name << std::endl;

      //This is DIFFERENT than what's in the nu_e python code, because I am almost 100% certain that that is a mistake.
      //In their code, the cv universe shifts electron energy raw down by prong_ECALCalibE*0.058
      //but the universe doesn't use that value at all, and instead varies around the unshifted value. The result is that the "central" value is always far than any of the universes.
      //So this first does the canonical downward 5.8% shift, and then shifts THAT value by the uncertainties above.
      double cv_shift = T::GetVecElem("prong_ECALCalibE", 0) * EM_ENERGY_SCALE_SHIFT_ECAL;
      return cv_shift + (T::m_nsigma * EM_ENERGY_SCALE_UNCERTAINTY.at(region) * T::GetVecElem(branch_name.c_str(), 0));
    }

    std::string ShortName() const override { return "Electron_Energy_Shift_" + region; }
    std::string LatexName() const override { return "Electron Energy Shift " + region; }
    bool IsVerticalOnly()  const override { return false; }  
  };

  //=================================================================================
  // Leakage correction universe, shifts the value of the leakage correction that gets applied to E_avail
  //=================================================================================
  template<class T>
  class LeakageUniverse: public T
  {
    public:
    
    //LeakageUniverse(typename T::config_t chw, double nsigma, string region ) : T(chw, nsigma) {}
    LeakageUniverse(typename T::config_t chw, double nsigma) : T(chw, nsigma) {}

    //These functions override functions from the base universe (my CVUniverse)
    double GetLeakageCorrection() const override {
      static const double LEAKAGE_UNCERTAINTY = 0.002; // +/- 2MeV to the leakage correction 

      double shift;
      //is this logic necessary? this is present in the python code but I'm not really sure why it's necessary...
      if ( abs(T::GetInt("mc_primaryLepton")) == 11 ){ shift = T::m_nsigma * LEAKAGE_UNCERTAINTY; }
      else { shift = 0; }
      
      //calls the cv universe to get the base leakage correction, then puts the systematic variation on top of that
      return CVUniverse::GetLeakageCorrection() + shift;
    }

    std::string ShortName() const override { return "Leakage_Uncertainty"; }
    std::string LatexName() const override { return "Leakage Uncertainty"; }
    bool IsVerticalOnly()  const override { return false; }  
  };



  //Initiation of the helper function to get electron energy shift universes
  template <class T>
  std::map<std::string, std::vector<T*> >
  GetElectronEnergyShiftSystematicsMap(typename T::config_t chain)
  {
    std::map<std::string, std::vector<T*> > ret;

    //Variations for both the ecal and the hcal
    ret["Electron_Energy_Shift_ECAL"].push_back( new ElectronEnergyShiftUniverse<T>(chain, -1., "ECAL"));
    ret["Electron_Energy_Shift_ECAL"].push_back( new ElectronEnergyShiftUniverse<T>(chain, 1., "ECAL"));
    ret["Electron_Energy_Shift_HCAL"].push_back( new ElectronEnergyShiftUniverse<T>(chain, -1., "HCAL"));
    ret["Electron_Energy_Shift_HCAL"].push_back( new ElectronEnergyShiftUniverse<T>(chain, 1., "HCAL"));
    
    return ret;
  }

  //Initiation of the helper function to get leakage universes
  template <class T>
  std::map<std::string, std::vector<T*> >
  GetLeakageSystematicsMap(typename T::config_t chain)
  {
    std::map<std::string, std::vector<T*> > ret;

    //Variations for both the ecal and the hcal
    ret["Leakage_Uncertainty"].push_back( new LeakageUniverse<T>(chain, -1.));
    ret["Leakage_Uncertainty"].push_back( new LeakageUniverse<T>(chain, 1.));
    
    return ret;
  }
  
} //end namespace PlotUtils

#endif // ELECTRONSYSTEMATICS_H
