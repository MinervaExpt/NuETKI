#ifndef ELECTRONSYSTEMATICS_H
#define ELECTRONSYSTEMATICS_H

//NuE TKI / PlotUtils includes
#include "event/CVUniverse.h"

//ROOT includes
#include "Math/Vector3D.h"
#include "Math/AxisAngle.h"
#include <TRandom3.h>

//C++ includes
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
    //Constructor
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
    //Constructor
    LeakageUniverse(typename T::config_t chw, double nsigma) : T(chw, nsigma) {}

    //These functions override functions from the base universe (my CVUniverse)
    double GetLeakageCorrection() const override {
      static const double LEAKAGE_UNCERTAINTY = 0.002; // +/- 2MeV to the leakage correction 

      double shift;
      //is this logic necessary? this is present in the python code but I'm not really sure why it's necessary...
      if ( abs(T::GetInt("mc_primaryLepton")) == 11 ){ shift = T::m_nsigma * LEAKAGE_UNCERTAINTY; }
      else { shift = 0; }
      
      //calls the cv universe to get the base leakage correction, then puts the systematic variation on top of that
      return T::GetLeakageCorrection() + shift;
    }

    std::string ShortName() const override { return "Leakage_Uncertainty"; }
    std::string LatexName() const override { return "Leakage Uncertainty"; }
    bool IsVerticalOnly()  const override { return false; }  
  };

  
  //=================================================================================
  // Electron Angle Shift Universe, shifts the angle (theta) of the electron track (aka cone)
  // 100 universes, with random axis angles around which to shift and random values by which to shift
  // axis angle is uniform from 0 to pi, shift_angle is gaussian centered at 0with a width of 1 sigma, i.e. 1e-3
  //=================================================================================
  template<class T>
  class ElectronAngleShiftUniverse: public T
  {
    private:
      double axis_angle;
      double shift_angle;

    public:
    //Constructor
    ElectronAngleShiftUniverse(typename T::config_t chw, double nsigma) : T(chw, nsigma) {
      static TRandom3 rng(0);  // seeded once      
      axis_angle  = rng.Uniform(0.0, M_PI);
      shift_angle = rng.Gaus( 0.0, this->m_nsigma * 0.001 ); //one mrad estimate for electron uncertainty?? it kinda seems REALLY small
    }
      
    ROOT::Math::XYZVector GetElectronP3D() const override {
      // Call base-class implementation (CV universe)
      ROOT::Math::XYZVector p = T::GetElectronP3D();

      // Unit vector along p
      ROOT::Math::XYZVector unit_3_vector = p / p.R();
      
      // First normal vector: p × x̂
      ROOT::Math::XYZVector normal_vector1 = p.Cross(ROOT::Math::XYZVector(1.0, 0.0, 0.0));

      // Normalize if non-zero
      if (normal_vector1.R() != 0.0) {
        normal_vector1 *= 1.0 / normal_vector1.R();
      }

      // First rotation: rotate normal_vector1 around unit_3_vector
      ROOT::Math::AxisAngle r1(unit_3_vector, axis_angle);
      ROOT::Math::XYZVector rotation_axis = r1(normal_vector1);
      
      // Second rotation: rotate p around the rotated axis
      ROOT::Math::AxisAngle r2(rotation_axis, shift_angle);

      return r2(p);
    }
    
    double GetElectronTheta() const override {
      //double ret = T::GetElectronTheta() + shift_angle;
      double ret = this->GetElectronP3D().Theta();
      //std::cout << "Electron Angle Shift Universe, shifted theta = " << ret << std::endl;
      return ret; 
      //return T::GetElectronTheta() + shift_angle; 
    }
    std::string ShortName() const override { return "Electron_Angle_Shift"; }
    std::string LatexName() const override { return "Electron Angle Shift"; }
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

  //Initiation of the helper function to get electron angle shift universes
  template <class T>
  std::map<std::string, std::vector<T*> >
  GetElectronAngleShiftSystematicsMap(typename T::config_t chain)
  {
    std::map<std::string, std::vector<T*> > ret;

    //maybe the number of electron angle universes shouldn't be hardcoded but I'll worry about that later
    for (int i=0; i<100; i++){
      ret["Electron_Angle_Shift"].push_back( new ElectronAngleShiftUniverse<T>(chain, 1.));
    }
    return ret;
  }

} //end namespace PlotUtils

#endif // ELECTRONSYSTEMATICS_H
