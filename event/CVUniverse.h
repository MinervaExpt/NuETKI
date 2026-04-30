// ============================================================================
// Base class for an un-systematically shifted (i.e. CV) universe. Implement
// "Get" functions for all the quantities that you need for your analysis.
//
// This class inherits from PU::MinervaUniverse, which in turn inherits from
// PU::BaseUniverse. PU::BU defines the interface with anatuples.
// 
// Within the class, "WeightFunctions" and "MuonFunctions" are included to gain
// access to standardized weight and muon variable getters. See:
// https://cdcvs.fnal.gov/redmine/projects/minerva-sw/wiki/MinervaUniverse_Structure_
// for a full list of standardized functions you can use. In general, if a
// standard version of a function is available, you should be using it.
// =============================================================================
#ifndef CVUNIVERSE_H
#define CVUNIVERSE_H

#include "PlotUtils/MinervaUniverse.h"
#include "PlotUtils/GeantHadronSystematics.h"
#include "PlotUtils/CaloCorrection.h"

#include <iostream>
#include <TMath.h>
#include <vector>
#include <TRandom.h>
#include <TMath.h>

//#include "Math/AxisAngle.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/Rotation3D.h"
#include "Math/AxisAngle.h"
#include "Math/RotationX.h"

#include <TVector3.h>
#include <TLorentzVector.h>

class CVUniverse : public PlotUtils::MinervaUniverse {
  private:
    mutable int m_cachedProtonIndex = -9999;
    mutable bool m_protonIndexCached = false; 
  public:
  #include "PlotUtils/MuonFunctions.h" // GetMinosEfficiencyWeight
  #include "PlotUtils/TruthFunctions.h" //Getq3True
  #include "PlotUtils/RecoilEnergyFunctions.h" //GetRecoilEnergy
  #include "PlotUtils/WeightFunctions.h"

  // ========================================================================
  // Constructor/Destructor
  // ========================================================================
  CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma = 0)
      : PlotUtils::MinervaUniverse(chw, nsigma) {}

  virtual ~CVUniverse() {}

  // ========================================================================
  // Quantities defined here as constants for the sake of below. Definition
  // matched to Dan's CCQENuInclusiveME variables from:
  // `/minerva/app/users/drut1186/cmtuser/Minerva_v22r1p1_OrigCCQENuInc/Ana/CCQENu/ana_common/include/CCQENuUtils.h`
  // ========================================================================
  static constexpr double M_n = 939.56536;  //MeV
  static constexpr double M_p = 938.272013; //MeV
  static constexpr double M_pi = 139.57039; //MeV
  //static constexpr double M_p = 938.27; //MeV
  //static constexpr double M_pi = 139.57; //MeV
  static constexpr double M_e = 0.51099895; //MeV
  static constexpr double M_nucleon = (1.5*M_n+M_p)/2.5;

  static constexpr int PDG_n = 2112;
  static constexpr int PDG_p = 2212;

  static constexpr double pi = 3.141592653589793;

  static constexpr double beam_tilt = 0.05887; //beam tilt downward with respect to z axis, in radians. 3.3 degrees, or 58 mrads.
  //According to the NuMI NIM paper, the exact angle is 3.34349, but at that precision you need to specify where along the target hall/minos hall you are
  //since the curvature of the earth becomes relevant at the ~1km scale lol

  //necessary for caching my true proton index on an event by event basis so I don't have to find it again for every variable
  //basically just says everytime we set a new entry, reset the "we've found the primary proton" flag
  void SetEntry(Long64_t entry) {
    PlotUtils::BaseUniverse::SetEntry(entry);  // call base class first
    m_protonIndexCached = false;  // invalidate cache
  }  
  // ========================================================================
  // Write a "Get" function for all quantities access by your analysis.
  // For composite quantities (e.g. Enu) use a calculator function.
  //
  // In order to properly calculate muon variables and systematics use the
  // various functions defined in MinervaUniverse.
  // E.g. GetPmu, GetEmu, etc.
  // ========================================================================

  // ========================================================================
  // Reco kinematics
  // ========================================================================  

  //total energy including rest mass in MeV
  virtual double GetElectronEnergy() const {
    return GetElectronEnergyRaw() + GetEMEnergyShift();
  }

  //uncorrected electron energy, in MeV
  virtual double GetElectronEnergyRaw() const {
    return GetVecElem("prong_part_E", 0, 3);
  }

  //the python nu_e people scale electron energy in the ECAL down 5.8% in the CV, this gets overridden in the ElectronEnergyShift universes though
  //Also in MeV
  virtual double GetEMEnergyShift() const {
    static const double EM_ENERGY_SCALE_SHIFT_ECAL = -0.058;
    return GetVecElem("prong_ECALCalibE", 0)*EM_ENERGY_SCALE_SHIFT_ECAL;    
  }

  //ElectronP3D returns a 3 momentum WITH RESPECT TO THE BEAM FRAME already, no need to rotate it
  //returns in rad
  virtual double GetElectronTheta() const { 
    ROOT::Math::XYZVector p = GetElectronP3D();
    //std::cout << "CVUniverse electron theta (rad) = " << p.Theta() << std::endl;
    return p.Theta();
  }

  //mostly for plotting
  virtual double GetElectronThetaDeg() const {
    double rad = GetElectronTheta();
    double deg = rad*180/pi;
    return deg;
  }

  //Electron kinetic energy (T) in MeV
  //Actually I think this is wrong, the E that comes from prong_part_E doesn't have electron mass added...
  //probably cause electron rest mass is basically negligible (0.5MeV, or 0.0005 GeV, or less than the average hit energy...)
  //But either way, makes no sense to subtract it back off if it was never added in the first place
  virtual double GetElectronT() const{
    double leptonT = GetElectronEnergy() - M_e; 
    if (leptonT > 0){
      return leptonT;
    }
    else return 0; 
    //return sqrt(std::max(0, (pow(GetElectronEnergy(), 2) - M_lep_sqr)));
  }
  
  //Electron momentum in MeV
  //And actually same story as above here, this doesn't really make sense. It should just be prong_part_E straight up
  //which again makes sense, for relativistic particles E = p (pretend its squiggly)
  //Doesn't really hurt to do it, but the value is gonna be the same basically.
  virtual double GetElectronP() const{
    double M_lep_sqr = pow(M_e, 2);
    double leptonP_sqr = pow(GetElectronEnergy(), 2)-M_lep_sqr;
    //sanity check, this should give the same result
    //double leptonP = GetElectronP3D().R();
    //std::cout << "lepton P from energy = " << sqrt(leptonP_sqr) << ", lepton P from vector = " << leptonP << std::endl;
    if (leptonP_sqr > 0){
      return sqrt(leptonP_sqr);
    }
    else return 0; 
    //return sqrt(std::max(0, (pow(GetElectronEnergy(), 2) - M_lep_sqr)));
  }

  //Returns electron 3 momentum in MeV, w.r.t the beam axis and scaled by appropriately scaled by GetEMEnergyShift
  virtual ROOT::Math::XYZVector GetElectronP3D() const{
    const std::vector<std::vector<double>>& electronp = GetVecOfVecDouble("prong_part_E");

    double scale = 0.0;
    if (electronp[0][3] > 0.0) { scale = GetEMEnergyShift() / electronp[0][3]; } //EMEnergyShift is a flat value in MeV, divide by electron energy to get scaling factor
    // create momentum vector scaled by (1 + scale)
    ROOT::Math::XYZVector p( electronp[0][0], electronp[0][1], electronp[0][2] );
    p *= (1.0 + scale);
    
    // rotate by beam angle about X
    ROOT::Math::RotationX r(-beam_tilt);
    ROOT::Math::XYZVector ret = r(p);

    return ret;
  } 

  //Transverse momentum of the lepton (electron for me), in MeV
  virtual double GetElectronPt() const{
    //return (GetElectronP() * std::sin(GetElectronTheta())); //older version, should give the same result either way
    return GetElectronP3D().Rho();
  }

  //Longitudinal momentum of the lepton (electron for me), in MeV
  virtual double GetElectronPParallel() const{
    //return (GetElectronP() * std::cos(GetElectronTheta())); //older version, should give the same result either way
    return GetElectronP3D().Z();
  }

  //proton total energy in MeV
  virtual double GetProtonEnergy() const{
    return GetDouble("MasterAnaDev_proton_E_fromdEdx");
  }
  //MasterAnaDev_proton_theta, in rad
  virtual double GetProtonTheta() const {
    return GetDouble("MasterAnaDev_proton_theta");
  }

  //same thing but in degrees, useful for plotting. Also if there's no proton candidate don't convert it
  virtual double GetProtonThetaDeg() const {
    double protonTheta = GetProtonTheta();
    //std::cout << "True (highest) Proton Momentum: " << protonP << std::end; 
    if (protonTheta > -9999){
      return protonTheta * (180/pi);
    }
    else {
      return protonTheta;
    }
  }

  //Kinetic energy (T) of ONLY the primary proton candidate (in MeV)
  virtual double GetProtonT() const {
    //double protonT = GetDouble("MasterAnaDev_proton_calib_energy"); //TO DO !! MAKE SURE THIS IS ACTUALLY T AND NOT TOTAL E
    double protonT = GetDouble("MasterAnaDev_proton_T_fromdEdx"); //what's the difference b/w these two I wonder... seems to be pretty small for most events (<1 MeV)
    return protonT; // = -9999 if there's no proton candidate
  }

  //Total Kinetic energy (sigma T) of all proton candidates (MeV). Dan uses this variable in his QELike analysis... but not in TKI?
  virtual double GetSumProtonT() const {
    double primary_proton_energy = GetDouble("MasterAnaDev_proton_calib_energy");
    if (primary_proton_energy < -1) { primary_proton_energy= 0; } //branch is set to -9999 by default, but tbh this shouldn't matter cause if there's no primary there shouldn't be any secondaries....
    
    //very few events have tracked secondary protons honestly... about 10% of my overall sample, fewer for signal

    //couple subtleties about this...
    // 1. make sure this is never -9999
    // 2. this is just kinetic energy, I want total right??? Because E_nu = E_e + E_avail, so it must include the rest mass right? same for above tbh, i think calib energy is just kinetic...
    double sum_secondary_proton_energy = 0;
    std::vector<double> secondary_protons = GetVecDouble("MasterAnaDev_sec_protons_T_fromCalo");
    for (int i = 0; i < secondary_protons.size(); i++){
      sum_secondary_proton_energy += secondary_protons[i];
    }
    //std::cout << "Sum Proton T RECO: " << primary_proton_energy + sum_secondary_proton_energy << std::endl;
    return primary_proton_energy + sum_secondary_proton_energy;
  }
  
  //MasterAnaDev_proton_P_fromdEdx, in MeV
  virtual double GetProtonP() const {
    double protonP = GetDouble("MasterAnaDev_proton_P_fromdEdx");
    return protonP; // = -9999 if no proton candidate
  }

  //proton transverse momentum, in MeV
  virtual double GetProtonPt() const{
    return (GetProtonP() * std::sin(GetProtonTheta()));
  }

  //proton longitudinal momentum, in MeV
  virtual double GetProtonPParallel() const{
    return (GetProtonP() * std::cos(GetProtonTheta()));
  }

  //Hadronic available energy (no neutrons) in MeV
  virtual double GetEavail() const {
    static const double AVAILABLE_E_CORRECTION = 1.17;   
    double E = GetDouble("blob_recoil_E_tracker") + GetDouble("blob_recoil_E_ecal");
    //double Eavail = (E * 1.17 - (0.008 * GetVecElem("prong_part_E", 0, 3)) + 5); 
    double Eavail = (E * AVAILABLE_E_CORRECTION) - GetLeakageCorrection();
    return Eavail;
  }

  //Flat correction to E_avail, based off electron shower energy. In MeV
  virtual double GetLeakageCorrection() const {
    static const double LEAKAGE_CORRECTION = 0.008; //0.8% on ElectronEnergyRaw
    static const double LEAKAGE_BIAS = 5; //flat -5 MeV to the e_avail
    // for some reason the python nu_e code applies this -5 shift ONLY to the CV and not the universes... but why???
    //according to david the systematics aren't based around the CV anyways, but instead the average of the universes so idk maybe it doesn't really matter
    
    return ( GetElectronEnergyRaw() * LEAKAGE_CORRECTION ) + LEAKAGE_BIAS;
  }

  //This is the same e-avail as above, but with all tracked proton calibrated kinetic energy subtracted off (
  //this is a completely reco variable by the way, in case I get lost lol
  double GetModifiedEavail() const {
    double Eavail = GetEavail();
    
    //Need to make sure there is a primary proton, otherwise we'll just add 9999 to Eavail by accident...
    double primary_proton_energy = GetDouble("MasterAnaDev_proton_calib_energy");
    if (primary_proton_energy < -1) { primary_proton_energy= 0; } //branch is set to -9999 by default

    //very few events have tracked secondary protons honestly... about 10% of my overall sample, fewer for signal
    double sum_secondary_proton_energy = 0;
    std::vector<double> secondary_protons = GetVecDouble("MasterAnaDev_sec_protons_T_fromCalo");
    for (int i = 0; i < secondary_protons.size(); i++){
      sum_secondary_proton_energy += secondary_protons[i];
    }

    double ModEavail = Eavail - primary_proton_energy - sum_secondary_proton_energy;
    //std::cout << "Reco Eavail:                                     " << Eavail << std::endl;
    //std::cout << "primary proton energy:                           " << primary_proton_energy << std::endl;
    //std::cout << "total secondary proton energy:                   " << sum_secondary_proton_energy << std::endl;
    //std::cout << "Modified (proton energy subtracted) Reco Eavail: " << ModEavail << " \n" << std::endl;
    return ModEavail;
  }

  //in MeV
  double GetEnu() const {
    double Enu = GetEavail() + GetElectronEnergy();
    //std::cout << "Reco Enu: " << Enu << std::endl;
    return Enu;
  }
  
  //Returns a root XYZVector object containing the lepton (electron) transverse 3 momentum (in MeV)
  ROOT::Math::XYZVector GetElectronPt3D() const{
    ROOT::Math::XYZVector electronP_vec = GetElectronP3D(); //GetElectronP3D ALREADY returns rotated to beam frame

    //Now set z to zero in the beam frame to keep only the transverse components
    electronP_vec.SetZ(0); 

    return electronP_vec;
  }
  
  //Returns a root XYZVector object containing the lepton (electron) transverse 3 momentum (in MeV), in BEAM FRAME
  ROOT::Math::XYZVector GetProtonPt3D() const{
    ROOT::Math::XYZVector protonP_vec(GetDouble("MasterAnaDev_proton_Px_fromdEdx"), GetDouble("MasterAnaDev_proton_Py_fromdEdx"), GetDouble("MasterAnaDev_proton_Pz_fromdEdx"));
    ROOT::Math::RotationX r(-beam_tilt);
    ROOT::Math::XYZVector protonPt_vec = r(protonP_vec);
    protonPt_vec.SetZ(0);

    return protonPt_vec;
  }
  
  //Returns a root XYZVector object containing the sum of the proton transverse 3 momentum and the lepton (electron) transverse 3 momentum in MeV
  //which is then used to calculate TKI variables
  //Remember: z direction != beam direction so transverse doesn't exactly mean z components are zero, although they should be small
  ROOT::Math::XYZVector GetDeltaPt3D() const{
    ROOT::Math::XYZVector electronPt_vec = GetElectronPt3D(); //returns in beam frame
    ROOT::Math::XYZVector protonPt_vec = GetProtonPt3D(); //returns in beam frame
    //ROOT::Math::XYZVector deltaP_vec = electronP_vec + protonP_vec; //sum of the full 3 momenta. Not sure if I need this so commenting it out for now
    ROOT::Math::XYZVector deltaPt_vec = electronPt_vec + protonPt_vec;
    //std::cout << "delta P total (kinda useless?): " << deltaP_vec.R() << std::endl;
    //std::cout << "delta Ptx: " << deltaPt_vec.X() << std::endl;
    //std::cout << "delta Pty: " << deltaPt_vec.Y() << std::endl;
    //std::cout << "delta Pt: " << deltaPt_vec.R() << std::endl;
    //std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n" << std::endl;
    return deltaPt_vec; //returns in beam frame
  }
  
  //delta pT (magnitude of the vector), in MeV
  double GetDeltaPt() const{
    ROOT::Math::XYZVector deltaPt_vec = GetDeltaPt3D();
    //std::cout << "delta Pt RECO: " << deltaPt_vec.R() << std::endl;
    return deltaPt_vec.R();
  }

  //Magnitude of the component of delta pT orthogonal to electron pT (MeV)
  double GetDeltaPtX() const{
    ROOT::Math::XYZVector electronPt_vec = GetElectronPt3D(); //returns in beam frame
    ROOT::Math::XYZVector deltaPt_vec = GetDeltaPt3D(); //returns in beam frame

    //expression for the scalar rejection, aka perpendicular dot product
    double num = deltaPt_vec.Y()*electronPt_vec.X() - deltaPt_vec.X()*electronPt_vec.Y();
    double denom = electronPt_vec.R();
    return num/denom;
  }

  //Magnitude of the component of  delta pT parallel to electron pT (MeV)
  double GetDeltaPtY() const{
    ROOT::Math::XYZVector electronPt_vec = GetElectronPt3D(); //returns in beam frame
    ROOT::Math::XYZVector deltaPt_vec = GetDeltaPt3D(); //returns in beam frame

    //std::cout << "\nelectron pT vec: ( " << electron_pt_in_beam_frame.X() << " , " << electron_pt_in_beam_frame.Y() << " , " << electron_pt_in_beam_frame.Z() << " ) " << std::endl;
    //std::cout << "delta pT vec: ( " << delta_pt_in_beam_frame.X() << " , " << delta_pt_in_beam_frame.Y() << " , " << delta_pt_in_beam_frame.Z() << " ) " << std::endl;
    //std::cout << "dot product of the two: " << delta_pt_in_beam_frame.Dot(electron_pt_in_beam_frame) << std::endl;
    //std::cout << "magnitude (norm) of e_pT vec: " << electron_pt_in_beam_frame.R() << std::endl;
    
    //expression for the scalar projection
    double num = deltaPt_vec.Dot(electronPt_vec);
    double denom = electronPt_vec.R();

    //std::cout << "final result: " << num/denom << "\n" << std::endl;
    return -num/denom; //negative cause I had it backwards before, this is the correct convention
  }
  
  //Alpha_t, the TKI boosting angle. it's the angle between inverted electron pT and delta pT
  double GetAlphaT() const{
    ROOT::Math::XYZVector electronPt_vec = GetElectronPt3D(); 
    ROOT::Math::XYZVector deltaPt_vec = GetDeltaPt3D();

    //angle between the two vectors
    double numerator = ((-1*electronPt_vec).Dot(deltaPt_vec));
    double denominator = ( electronPt_vec.R() * deltaPt_vec.R() );
    double alpha = std::acos(numerator/denominator);

    return alpha * 180/pi;  //return alpha in degrees cause that's how I've set up my bins for now
  }

  //Phi_t 
  double GetPhiT() const{
    ROOT::Math::XYZVector electronPt_vec = GetElectronPt3D();
    ROOT::Math::XYZVector protonPt_vec = GetProtonPt3D();

    //angle between the two vectors 
    double numerator = ((-1*electronPt_vec).Dot(protonPt_vec));
    double denominator = ( electronPt_vec.R() * protonPt_vec.R() );
    double phi = std::acos(numerator/denominator);

    //std::cout << "TKI phi: " << phi * 180/pi << std::endl;
    return phi * 180/pi;  //return in degrees cause that's how I've set up my bins for now
  }
  
  //This is the variable reported by Jeffrey & Xianguo, details in PhysRevLett.121.022504 (2018 Xianguo TKI paper)
  //
  //Jeffrey also had an alternate "CCQE" method of calculating this variable, found in:
  // /exp/minerva/app/users/kleykamp/cmtuser/Minerva_v22r1p1/Ana/NukeCCQE/scripts/code_templates/transverse_vars.cpp
  //Apparently that version assumes a CCQE interaction, and after some benchmarking it's NOT correct, so this is the right version
  //in MeV
  double GetDeltaPl() const{
    double DeltaPt = GetDeltaPt(); 
    double mass_nuke = 11188; //Mass of a carbon atom in MeV
    double binding_energy = 27.13; //this value comes from the Furmanski & Sobczyk paper: PhysRevC.95.065501
    
    double R = mass_nuke + GetElectronPParallel() + GetProtonPParallel() - ( GetElectronEnergy() + GetProtonEnergy() );
    //the first squared term is m_A', or the mass of the nuclear remnant after interaction
    //it's the carbon atom mass, minus neutron mass, plus binding energy (taken from Jeffrey's calculation)
    double numerator = pow(mass_nuke - M_p + binding_energy, 2) + pow(DeltaPt, 2);
    double DeltaPl = 0.5 * R - numerator/(2*R);
    return DeltaPl;
  }

  //This variable represents the initial struck nucleon momentum, determined by delta Pt and delta Pl.
  //in MeV
  double GetPn() const{
    double DeltaPt = GetDeltaPt();
    double DeltaPl = GetDeltaPl();
    return sqrt(pow(DeltaPt, 2) + pow(DeltaPl, 2));
  }
  
  virtual double GetQ2Reco() const{
    return GetDouble("qsquared_recoil");
  }

  //GetRecoilE is designed to match the NSF validation suite
  virtual double GetRecoilE() const {
    return GetVecElem("recoil_summed_energy", 0);
  }
  
  //in GeV? Definitely possible mistakes in here. ..
  virtual double Getq3() const{
    double eavail = GetEavail()/pow(10,3);
    double q2 = GetQ2Reco() / pow(10,6);
    double q3mec = sqrt(eavail*eavail + q2);
    return q3mec;
  }
  
  // ========================================================================
  // Truth kinematics
  // ========================================================================  

  //in MeV
  double GetElectronEnergyTrue() const {
    double out = GetVecElem("mc_primFSLepton", 3);
    //std::cout << "True Lepton Energy: " << out << std::endl;
    return out;
  }

  //in rad
  double GetElectronThetaTrue() const {
    std::vector<double> electronp = GetVecDouble("mc_primFSLepton");
    ROOT::Math::XYZVector p(electronp[0], electronp[1], electronp[2]);   
    ROOT::Math::RotationX r(-beam_tilt); //This is a slight rotation so that the theta we get is wrt to the beam direction and not the z axis, which points down to the ground 3.3 degrees.
    return (r(p)).Theta();    
  }

  //in deg
  double GetElectronThetaDegTrue() const {
    return GetElectronThetaTrue() * (180./pi);
  }
  
  //True lepton pT (in MeV)
  double GetElectronPtTrue() const {
    double leptonP = sqrt(pow(GetVecElem("mc_primFSLepton", 3), 2) - pow(M_e, 2));
    double out = leptonP * std::sin(GetElectronThetaTrue());
    //std::cout << "True lepton pT: " << out << std::endl;
    return out;
  }

  //True lepton p parallel (in MeV)
  double GetElectronPParallelTrue() const {
    double leptonP = sqrt(pow(GetVecElem("mc_primFSLepton", 3), 2) - pow(M_e, 2));
    double out = leptonP * std::cos(GetElectronThetaTrue());
    //std::cout << "True lepton pT: " << out << std::endl;
    return out;
  }

  //gets highest true proton kinetic energy in MeV, -9999 if no true protons
  double GetProtonTTrue() const {
    int highestProtonIndex = GetHighestEnergySignalProtonIndex();
    if (highestProtonIndex > -1){
      double out = GetVecElem("mc_FSPartE", highestProtonIndex) - M_p;
      //std::cout << "True (highest) Proton KE: " << out << std::endl;
      return out;
    }
    else {
      return -9999;
    }
  }

  //proton theta (spherical coords, in deg) wrt beam direction, -9999 if no true protons
  //in rad
  double GetProtonThetaTrue() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      ROOT::Math::XYZVector p(GetVecElem("mc_FSPartPx", i), GetVecElem("mc_FSPartPy", i), GetVecElem("mc_FSPartPz", i));
      ROOT::Math::RotationX r(-beam_tilt);
      double out = (r(p)).Theta();
      //std::cout << "True (highest energy) proton Theta: " << out << std::endl;
      return out;
    }
    else {
      return -9999;
    }
  }

  //to deg
  double GetProtonThetaDegTrue() const {
    double proton_theta = GetProtonThetaTrue();
    if (proton_theta == -9999) {
      return -9999;
    } else {
      return proton_theta * (180./pi);
    }
  }


  //Returns momentum of the highest energy true final state proton (in MeV)
  double GetProtonPTrue() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      //double protonP = sqrt(pow(GetVecElem("mc_FSPartE", i),2) - pow(M_p, 2));
      //return protonP;
      ROOT::Math::XYZVector p(GetVecElem("mc_FSPartPx", i), GetVecElem("mc_FSPartPy", i), GetVecElem("mc_FSPartPz", i));
      return p.R();
    }
    else {
      return -9999;
    }
  }

  //true proton pT (in MeV) of highest energy proton
  double GetProtonPtTrue() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      double protonP = GetProtonPTrue();
      double protonAngleRad = GetProtonThetaTrue();
      
      double out = protonP * std::sin(protonAngleRad);
      //std::cout << "True (highest energy) proton pT: " << out << std::endl;
      return out;
    }
    else {
      return -9999;
    }
  }

  //true proton pT (in MeV) of highest energy proton
  double GetProtonPParallelTrue() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      double protonP = GetProtonPTrue();
      double protonAngleRad = GetProtonThetaTrue();
      
      double out = protonP * std::cos(protonAngleRad);	
      //std::cout << "True (highest energy) proton pT: " << out << std::endl;
      return out;
    }
    else {
      return -9999;
    }
  }


  //Total Kinetic energy (sigma T) of all proton candidates
  double GetSumProtonTTrue() const {
    double T_p = 0; //sum of proton kinetic energies

    int nFSPart = GetInt("mc_nFSPart");
    std::vector<double> FSPartE = GetVecDouble("mc_FSPartE");
    std::vector<int> FSPartPDG = GetVecInt("mc_FSPartPDG");

    for (int i=0; i<nFSPart; i++){
      if (FSPartPDG[i] == 2212){ 
	T_p += FSPartE[i] - M_p; //add only kinetic (not total) energy for protons
      } 
    }
    return T_p;
  }

  //Angle (in deg) between true primary electron and true highest energy proton in rad, -9999 if no protons
  double GetOpeningAngleTrue() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      std::vector<double> electronp = GetVecDouble("mc_primFSLepton");

      ROOT::Math::XYZVector electron(electronp[0], electronp[1], electronp[2]);
      ROOT::Math::XYZVector proton(GetVecElem("mc_FSPartPx", i), GetVecElem("mc_FSPartPy", i), GetVecElem("mc_FSPartPz", i));
      
      double dotProd = electron.Dot(proton);

      //no shot this just works first try
      double out = std::acos(dotProd/( sqrt( electron.Mag2() * proton.Mag2() ))) * (180/pi);
      //std::cout << "True angle between electron and (highest energy) proton: " << out << std::endl;
      return out;
    }
    else {
      return -9999;
    }
  }

  //E_avail (MeV) as defined in truth, basically summing all sources of visible energy
  //This is copied from the python nu_e framework ( tools/KinematicsCalculator.py )
  double GetEavailTrue() const {
    double T_p = 0; //sum of proton kinetic energies
    double T_pi = 0; //sum of charged pion kinetic energies
    double E_pi0 = 0; //sume of neutral pion total energies
    double E_s = 0; //Sum of (strange baryon energy - proton mass)
    double E_sbar = 0; //Sum of (anti baryon energy + proton mass)
    double E_other = 0; //sum of total energy of any other particles, NOT INCLUDING NEUTRONS

    int nFSPart = GetInt("mc_nFSPart");
    std::vector<double> FSPartE = GetVecDouble("mc_FSPartE");
    std::vector<int> FSPartPDG = GetVecInt("mc_FSPartPDG");

    for (int i=0; i<nFSPart; i++){
      if (abs(FSPartPDG[i]) == 11 || abs(FSPartPDG[i]) == 13 || FSPartPDG[i] == 2112 || FSPartPDG[i] > 1000000000){
	continue; //don't want to count leptons, neutrons, or nuclear remnants
      } else if (FSPartPDG[i] == 2212){ 
	T_p += (FSPartE[i] - M_p); //add only kinetic (not total) energy for protons
      } else if (abs(FSPartPDG[i]) == 211){ 
	T_pi += (FSPartE[i] - M_pi); //add only kinetic (not total) energy for charged pions
      } else if (abs(FSPartPDG[i]) == 111){ 
	E_pi0 += FSPartE[i]; //add pi0 total energy
      } else if (FSPartPDG[i] > 2000){
	E_s += (FSPartE[i] - M_p); //add total energy - proton mass for strange baryons (why?)
      } else if (FSPartPDG[i] < -2000){
	E_sbar += (FSPartE[i] + M_p); //add total energy + proton mass for strange antibaryons (why???)
      } else {
	E_other += FSPartE[i]; //add total energy for anything else (mostly gammas, kaons, and eta mesons)
      }
    }
    double E_avail = T_p + T_pi + E_pi0 + E_s + E_sbar + E_other;
    return E_avail;
  }

  //MeV
  double GetDeltaPtTrue() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      ROOT::Math::XYZVector electronP_vec(GetVecElem("mc_primFSLepton", 0), GetVecElem("mc_primFSLepton", 1), GetVecElem("mc_primFSLepton", 2));
      ROOT::Math::XYZVector protonP_vec(GetVecElem("mc_FSPartPx",i), GetVecElem("mc_FSPartPy",i), GetVecElem("mc_FSPartPz",i));

      ROOT::Math::XYZVector neutrinoP_vec(GetVecElem("mc_incomingPartVec", 0), GetVecElem("mc_incomingPartVec", 1), GetVecElem("mc_incomingPartVec", 2));
      ROOT::Math::XYZVector unit_neutrino = neutrinoP_vec.Unit();      

      ROOT::Math::XYZVector electronPl_vec = unit_neutrino * electronP_vec.Dot(unit_neutrino);
      ROOT::Math::XYZVector protonPl_vec = unit_neutrino * protonP_vec.Dot(unit_neutrino);

      ROOT::Math::XYZVector electronPt_vec = electronP_vec - electronPl_vec;
      ROOT::Math::XYZVector protonPt_vec = protonP_vec - protonPl_vec;
      
      ROOT::Math::XYZVector deltaPt_vec = electronPt_vec + protonPt_vec;      
      return deltaPt_vec.R();
    }
    else {
      return -9999;
    }
  }

  //MeV
  double GetDeltaPtXTrue() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      ROOT::Math::XYZVector electronP_vec(GetVecElem("mc_primFSLepton", 0), GetVecElem("mc_primFSLepton", 1), GetVecElem("mc_primFSLepton", 2));
      ROOT::Math::XYZVector protonP_vec(GetVecElem("mc_FSPartPx",i), GetVecElem("mc_FSPartPy",i), GetVecElem("mc_FSPartPz",i));

      ROOT::Math::XYZVector neutrinoP_vec(GetVecElem("mc_incomingPartVec", 0), GetVecElem("mc_incomingPartVec", 1), GetVecElem("mc_incomingPartVec", 2));
      ROOT::Math::XYZVector unit_neutrino = neutrinoP_vec.Unit();

      ROOT::Math::XYZVector electronPl_vec = unit_neutrino * electronP_vec.Dot(unit_neutrino);
      ROOT::Math::XYZVector protonPl_vec = unit_neutrino * protonP_vec.Dot(unit_neutrino);

      ROOT::Math::XYZVector electronPt_vec = electronP_vec - electronPl_vec;
      ROOT::Math::XYZVector protonPt_vec = protonP_vec - protonPl_vec;

      ROOT::Math::XYZVector deltaPt_vec = electronPt_vec + protonPt_vec;      
      double num = deltaPt_vec.Y()*electronPt_vec.X() - deltaPt_vec.X()*electronPt_vec.Y();
      double denom = electronPt_vec.R();

      //double deltaPtX_1 = num/denom;
      //double deltaPtX_2 = unit_neutrino.Cross(electronPt_vec.Unit()).Dot(deltaPt_vec);
      //std::cout << "CARLOS TEST, deltaPtX_1 = " << deltaPtX_1 << ", deltaPtX_2 = " << deltaPtX_2 << std::endl;
      return num/denom;
      //return deltaPtX_2;
    }
    else {
      return -9999;
    }
  }

  //MeV
  double GetDeltaPtYTrue() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      ROOT::Math::XYZVector electronP_vec(GetVecElem("mc_primFSLepton", 0), GetVecElem("mc_primFSLepton", 1), GetVecElem("mc_primFSLepton", 2));
      ROOT::Math::XYZVector protonP_vec(GetVecElem("mc_FSPartPx",i), GetVecElem("mc_FSPartPy",i), GetVecElem("mc_FSPartPz",i));

      ROOT::Math::XYZVector neutrinoP_vec(GetVecElem("mc_incomingPartVec", 0), GetVecElem("mc_incomingPartVec", 1), GetVecElem("mc_incomingPartVec", 2));
      ROOT::Math::XYZVector unit_neutrino = neutrinoP_vec.Unit();

      ROOT::Math::XYZVector electronPl_vec = unit_neutrino * electronP_vec.Dot(unit_neutrino);
      ROOT::Math::XYZVector protonPl_vec = unit_neutrino * protonP_vec.Dot(unit_neutrino);

      ROOT::Math::XYZVector electronPt_vec = electronP_vec - electronPl_vec;
      ROOT::Math::XYZVector protonPt_vec = protonP_vec - protonPl_vec;
      ROOT::Math::XYZVector deltaPt_vec = electronPt_vec + protonPt_vec;      

      double num = deltaPt_vec.Dot(electronPt_vec);
      double denom = electronPt_vec.R();
      //std::cout << "TRUE delta PtY: " << num/denom. << std::endl;
      return -num/denom;
    }
    else {
      return -9999;
    }
  }

  //need true proton energy, true lepton energy, true lepton & proton longitudinal momentum, TRUE NUCLEUS MASS/BINDING, true delta pt
  double GetDeltaPlTrue() const{
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      ROOT::Math::XYZVector electronP_vec(GetVecElem("mc_primFSLepton", 0), GetVecElem("mc_primFSLepton", 1), GetVecElem("mc_primFSLepton", 2));
      ROOT::Math::XYZVector protonP_vec(GetVecElem("mc_FSPartPx",i), GetVecElem("mc_FSPartPy",i), GetVecElem("mc_FSPartPz",i));

      ROOT::Math::XYZVector neutrinoP_vec(GetVecElem("mc_incomingPartVec", 0), GetVecElem("mc_incomingPartVec", 1), GetVecElem("mc_incomingPartVec", 2));
      ROOT::Math::XYZVector unit_neutrino = neutrinoP_vec.Unit();

      ROOT::Math::XYZVector electronPl_vec = unit_neutrino * electronP_vec.Dot(unit_neutrino);
      ROOT::Math::XYZVector protonPl_vec = unit_neutrino * protonP_vec.Dot(unit_neutrino);

      double E_lep = GetVecElem("mc_primFSLepton", 3);
      double E_proton = GetVecElem("mc_FSPartE", i);
      double DeltaPt = GetDeltaPtTrue();

      double hit_nucleus = GetDouble("mc_targetZ");
      //double mass_nuke = 11.188*1000;
      double mass_nuke = 11.17486*1000;
      double binding_energy = 27.13;
      // Jeffrey's has an option to adjust TKI based on true nucleus... but it only checks for the nuclei in the target region
      // so it'll miss stuff like hydrogen, silicon, titanium
      // looking at some tuples, there's especially a decent amount of hydrogen interactions...
      // not sure how to handle that though? a normal hydrogen has no neutrons, no binding energy... it's just a proton basically
      // also I checked and it looks like all of my selected events that hit hydrogen are background, none are signal, which I guess makes sense...
      // ah well just dw bout it for now
      /*
      if (hit_nucleus == 8) { //oxygen
	mass_nuke = 14.903*1000;
	binding_energy = 24.1;
      } else if (hit_nucleus == 26) { //iron
	mass_nuke = 52.019*1000;
	binding_energy = 29.6;
      } else if (hit_nucleus == 82) { // lead
	mass_nuke = 193.000*1000;
	binding_energy = 22.8;
      } else {
	mass_nuke = 11.188*1000; // Return carbon as default
	binding_energy = 27.13;
	}*/
      
      double R = mass_nuke + electronPl_vec.R() + protonPl_vec.R() - ( E_lep + E_proton );
      double numerator = pow(mass_nuke - M_n + binding_energy, 2) + pow(DeltaPt, 2);
      
      double DeltaPl = 0.5 * R - numerator/(2*R);
      return DeltaPl;
    }
    else {
      return -9999;
    }
  }

  double GetPnTrue() const{
    double DeltaPt = GetDeltaPtTrue();
    if ( DeltaPt < -900 ) { //if this is below -900, or negative at all really, that means there were no true final state protons
      return -9999;
    }
    double DeltaPl = GetDeltaPlTrue();
    return sqrt(pow(DeltaPt, 2) + pow(DeltaPl, 2));
  }

  double GetAlphaTTrue() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      ROOT::Math::XYZVector electronP_vec(GetVecElem("mc_primFSLepton", 0), GetVecElem("mc_primFSLepton", 1), GetVecElem("mc_primFSLepton", 2));
      ROOT::Math::XYZVector protonP_vec(GetVecElem("mc_FSPartPx",i), GetVecElem("mc_FSPartPy",i), GetVecElem("mc_FSPartPz",i));

      ROOT::Math::XYZVector neutrinoP_vec(GetVecElem("mc_incomingPartVec", 0), GetVecElem("mc_incomingPartVec", 1), GetVecElem("mc_incomingPartVec", 2));
      ROOT::Math::XYZVector unit_neutrino = neutrinoP_vec.Unit();
      //ROOT::Math::XYZVector unit_neutrino(0, -0.05882, 0.99827); //assume neutrino came in exactly along beam direction, if I need to test for comparison

      ROOT::Math::XYZVector electronPl_vec = unit_neutrino * electronP_vec.Dot(unit_neutrino);
      ROOT::Math::XYZVector protonPl_vec = unit_neutrino * protonP_vec.Dot(unit_neutrino);

      ROOT::Math::XYZVector electronPt_vec = electronP_vec - electronPl_vec;
      ROOT::Math::XYZVector protonPt_vec = protonP_vec - protonPl_vec;

      ROOT::Math::XYZVector deltaPt_vec = electronPt_vec + protonPt_vec;

      double numerator = ((-1*electronPt_vec).Dot(deltaPt_vec));
      double denominator = ( electronPt_vec.R() * deltaPt_vec.R() );
      double alpha = std::acos(numerator/denominator);

      return alpha * 180/pi;
    }
    else {
      return -9999;
    }
  }

  double GetPhiTTrue() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      ROOT::Math::XYZVector electronP_vec(GetVecElem("mc_primFSLepton", 0), GetVecElem("mc_primFSLepton", 1), GetVecElem("mc_primFSLepton", 2));
      ROOT::Math::XYZVector protonP_vec(GetVecElem("mc_FSPartPx",i), GetVecElem("mc_FSPartPy",i), GetVecElem("mc_FSPartPz",i));

      ROOT::Math::XYZVector neutrinoP_vec(GetVecElem("mc_incomingPartVec", 0), GetVecElem("mc_incomingPartVec", 1), GetVecElem("mc_incomingPartVec", 2));
      ROOT::Math::XYZVector unit_neutrino = neutrinoP_vec.Unit();

      ROOT::Math::XYZVector electronPl_vec = unit_neutrino * electronP_vec.Dot(unit_neutrino);
      ROOT::Math::XYZVector protonPl_vec = unit_neutrino * protonP_vec.Dot(unit_neutrino);

      ROOT::Math::XYZVector electronPt_vec = electronP_vec - electronPl_vec;
      ROOT::Math::XYZVector protonPt_vec = protonP_vec - protonPl_vec;
      ROOT::Math::XYZVector deltaPt_vec = electronPt_vec + protonPt_vec;

      double numerator = ((-1*electronPt_vec).Dot(protonPt_vec));
      double denominator = ( electronPt_vec.R() * protonPt_vec.R() );
      double phi = std::acos(numerator/denominator);
      
      return phi * 180/pi;
    }
    else {
      return -9999;
    }
  }
  
  // ========================================================================
  // Other reco variables
  // Mostly used for evaluating cuts, mostly straight from the tuple (other than ESC)
  // ========================================================================  

  double GetImprovedNMichel() const {
    return GetInt("improved_nmichel");
  }
  
  //Calculates a chi^2 for proton end dE/dX if more than 5 nodes on proton track
  //Otherwise uses the original ESC method of checking node by node for last 5 nodes, and sets to an arbitrary, high chi^2 (75) if it fails any node cut
  double GetProtonESCNodeChi2() const {        
    double nodeEnergyVal = 0;
    double chi2=0;
    int n_nodes = GetInt("MasterAnaDev_proton_nodes_nodesNormE_sz");
    std::vector<double> nodesNormE = GetVecDouble("MasterAnaDev_proton_nodes_nodesNormE");

    if(n_nodes>5){ //Has more than 5 nodes, so use the chi^2 method
      double means[5] = {31.302, 11.418, 9.769, 8.675, 7.949};
      double sigmas[5] = {8.997, 3.075, 2.554, 2.484, 2.232};
      //These numbers represent a proton end track dE/dX fit for stopping protons, against which we calculate our chi^2 for the candidate proton track
      //If chi^2 is high, the fit doesn't accurately represent the dE/dX profile of the track in question
      //Which means it's probably not a stopping proton, rescattered inelastically or escaped or something
      //These values come from docdb 27170 from Dan

      for(int i=0;i<n_nodes;i++){
	if(i==6) break;
	if(i==0) nodeEnergyVal+=nodesNormE[0];
	else if(i==1) nodeEnergyVal+=nodesNormE[1];
	else nodeEnergyVal=nodesNormE[i];
	if(i>=1){
	  chi2+=(nodeEnergyVal-means[i-1])*(nodeEnergyVal-means[i-1])/(sigmas[i-1]*sigmas[i-1]);
	}
      }
    }
    else{ //Has 5 or fewer nodes, so use ESC node method
      bool pass = true;

      //Cut Values based on docdb 22302, where node 0 is the final node on the track and so on
      //These are horrible for my selection based on resolution plots, NEED TO ADJUST !
      //double cutval1 = 19; //Node 0+1
      //double cutval2 = 10; //Node 2
      //double cutval3 = 9;  //Node 3
      //double cutval4 = 8;  //Node 4
      //double cutval5 = 5;  //Node 5

      double cutval1 = 19; //Node 0+1
      double cutval2 = 1; //Node 2
      double cutval3 = 1;  //Node 3
      double cutval4 = 1;  //Node 4
      double cutval5 = 1;  //Node 5

      if(n_nodes==0) { pass = false; }///no nodes
      if(n_nodes>1) {
	if(nodesNormE[0]+nodesNormE[1] < cutval1) { pass = false; }
      }
      if(n_nodes>2) {
	if(nodesNormE[2] < cutval2) { pass = false; }
      }
      if(n_nodes>3) {
	if(nodesNormE[3] < cutval3) { pass = false; }
      }
      if(n_nodes>4) {
	if(nodesNormE[4] < cutval4) { pass = false; }
      }
      if(n_nodes>5) {
	if(nodesNormE[5] < cutval5){ pass = false; }
      }
      //survived loops?
      if(pass) { chi2=0; } //doesn't really matter, we just want to make sure it passes the cut
      else { chi2=75; } //Arbitrary, large value above the cut threshold
    }
    //std::cout << "chi2 = " << chi2 << std::endl;
    return chi2;
  }

  //Difference between furthest upstream blob start point z and vertex z (in mm)
  double GetBlobZDiffToVtxZ() const {
    std::vector<double> blobs_start_z = GetVecDouble("nonvtx_iso_blobs_start_position_z_in_prong");
    int n_blobs = GetInt("nonvtx_iso_blobs_start_position_z_in_prong_sz");
    std::vector<double> vtx = GetVecDouble("vtx");

    //find the furthest upstream blob start z, then calc the difference between that and vertex Z
    double minBlobStartZ = 100000; // some ridiculously large number
    if (n_blobs > 0){
      for (int i=0; i < n_blobs; i++){
	if (blobs_start_z[i] < minBlobStartZ){
	  minBlobStartZ = blobs_start_z[i];
	}
      }
      return minBlobStartZ - vtx[2];
    }
    else{ return 0; }
  }

  
  //This is E_lep * theta_lep^2, Ryan uses it for a cut but I don't at the moment
  double GetETheta() const {
    //double etheta = GetElectronEnergy() * pow(std::sin(GetElectronTheta()), 2);
    double etheta = GetElectronEnergy() * pow(GetElectronTheta(), 2);
    //std::cout << "E_lep*theta^2:      " << etheta << std::endl;
    return etheta;
  }
  
  //But he calls it (and it seems like it really should be) E_lep * sin^2(theta), so that's what im gonna calculate here...
  //not sure if this gets used anywhere
  double GetELepSin2Theta() const {
    //double etheta = GetElectronEnergy() * pow(std::sin(GetElectronTheta()), 2);
    double etheta = GetElectronEnergy() * pow(std::sin(GetElectronTheta()), 2);
    //std::cout << "E_lep*sin(theta)^2: " << etheta << std::endl;
    //std::cout << " " << std::endl;
    return etheta;
  }

  //Adaptation of Nimmy's side exiting muon cut that she's trying out
  bool GetSideExitingMuon() const {
    std::vector<int> isSideECAL = GetVecInt("MasterAnaDev_hadron_isSideECAL");
    std::vector<int> isODMatch = GetVecInt("MasterAnaDev_hadron_isODMatch");

    bool passCut = true;

    //Loop through hadron prongs, if any of them is in the side ecal OR has OD match, event fails
    //the idea being maybe one of them is actually the muon of a numu CC event that got misreco'd as a hadron, I guess?
    for (int i=0; i < isSideECAL.size(); i++){
      //std::cout << "CARLOS TEST SIDE EXITING MUON - prong " << i << ": isSideECAL=" << isSideECAL[i] << "; isODMatch=" << isODMatch[i] << std::endl;
      //if (isSideECAL[i] == 1 || isODMatch[i] == 1) { passCut = false; }
      if (isSideECAL[i] == 1 && isODMatch[i] == 1) { passCut = false; }
    }
    //std::cout << "pass cut = " << passCut << std::endl;
    return passCut;
  }

  //shower start Z in mm
  double GetShowerStartZ() const {
    return GetVecElem("prong_axis_vertex", 0, 2);
  }

  //proton track start Z in mm
  double GetProtonStartZ() const {
    return GetDouble("MasterAnaDev_proton_startPointZ");
  }

  //proton track end z in mm
  double GetProtonEndZ() const {
    return GetDouble("MasterAnaDev_proton_endPointZ");
  }

  //ratio of visible energy in DS_HCAL to DS_ECAL 
  double GetDSCalVisE() const {
    //std::cout << "prong_HCALVisE: " << GetVecElem("prong_HCALVisE", 0) << std::endl;
    //std::cout << "prong_ECALVisE: " << GetVecElem("prong_ECALVisE", 0) << std::endl;
    //if (GetVecElem("prong_HCALVisE", 0)==-999 && GetVecElem("prong_ECALVisE", 0)==-999) {	return 0; }
    if (GetVecElem("prong_HCALVisE", 0)<=0 || GetVecElem("prong_ECALVisE", 0)<=0) { return 0; }
    else { return GetVecElem("prong_HCALVisE", 0) / GetVecElem("prong_ECALVisE", 0); }
  }

  //ratio of visible energy in OD_HCAL to OD_ECAL 
  double GetODCalVisE() const {
    //std::cout << "prong_ODVisE: " << GetVecElem("prong_ODVisE", 0) << std::endl;
    //std::cout << "prong_SideECALVisE: " << GetVecElem("prong_SideECALVisE", 0) << std::endl;
    if (GetVecElem("prong_ODVisE", 0)<=0 || GetVecElem("prong_SideECALVisE", 0)<=0) { return 0; }
    else { return GetVecElem("prong_ODVisE", 0) / GetVecElem("prong_SideECALVisE", 0); }
  }

  //Jeremy's defined Psi variable, he used it for a cut but I don't because it biases what I get in terms of events with protons
  //something like ratio of energy outside the cone to energy inside the cone? Or something like that, maybe the other way around
  //can double check his thesis
  double GetPsi() const {
    return GetDouble("Psi");
  }

  int GetHasTracks() const {
    if (GetInt("n_prongs") > 0) return 1;
    else return 0;
  }

  //number of track start vertices, not including tracks that start from the end of another track
  double GetStartPointVertexMultiplicity() const {
    return GetInt("StartPointVertexMultiplicity");
  }

  double GetHasNoVertexMismatch() const {
    return GetInt("HasNoVertexMismatch");
  }

  //Number of tracks leaving the neutrino vertex
  double GetVertexTrackMultiplicity() const {
    return GetInt("VertexTrackMultiplicity");
  }

  //defined in Jeremy's thesis, I forget what it is right now though
  double GetNonMIPClusFrac() const {
    return GetVecElem("prong_NonMIPClusFrac", 0);
  }

  //For afterpulsing cut
  double GetFirstFireFraction() const {
    return GetVecElem("prong_FirstFireFraction", 0);
  }

  //still no 100% clear on how this is calculated in the Anatool...
  double GetEMLikeShowerScore() const {
    double score = GetVecElem("prong_part_score", 0);
    //std::cout << "emscore: " << score << std::endl;
    return score;
  }

  //defined in Jeremy's thesis, I forget what it is right now though
  double GetTransverseGapScore() const {
    return GetVecElem("prong_TransverseGapScore", 0);
  }

  //defined in Jeremy's thesis, something like average de/dx in the first 6 planes of the shower?
  double GetMeanFrontdEdx() const {
    return GetVecElem("prong_dEdXMeanFrontTracker", 0);
  }

  //
  double GetExitsBack() const {
    return GetInt("HasNoBackExitingTracks");
  }

  //Number of isolated blob prongs
  double GetNIsoBlobs() const {
    return GetInt("nonvtx_iso_blobs_energy_in_prong_sz");
  }

  //Summed energy of all isolated blob prongs
  double GetIsoBlobsEnergy() const {
    return GetDouble("nonvtx_iso_blobs_energy");
  }

  //in mm
  ROOT::Math::XYZTVector GetVertex() const {
    ROOT::Math::XYZTVector result;
    result.SetCoordinates(GetVec<double>("vtx").data());
    return result;
  }

  //This returns 0 or 1, 1 if its in it, 0 if its not
  double GetWithinFiducialApothem() const {
    const ROOT::Math::XYZTVector vertex = GetVertex();
    const double slope = -1./sqrt(3.);
    const double apothem = 850.;
    //So there are two checks here that if they're both true, means we have to be within the right apothem
    //take the abs value of both x and y, to move our point into the first quadrant just to not have to worry about sign
    //first we do the easy one-> is x less than the apothem (cause minerva's hexagon is flat on the left/right, pointy on top)
    //then, we check if y < (-1/sqrt(3))x + 2a/sqrt(3) which is the equation for the top right line segment (slope of -1/sqrt(3) and y offset of 2a/sqrt(3).
    //a being the apothem, and 2a/sqrt(3) being the distance from the center of the hexagon to the pointy top (or any other vertex)
    return (fabs(vertex.x()) < apothem)   &&   (fabs(vertex.y()) < slope*fabs(vertex.x()) + 2.*apothem/sqrt(3.));
  }

  //dead time branch/cut
  virtual int GetTDead() const {
    return GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj");
  }  

  //Needed for response/recoil systematics to work properly after Anezka's 2023 updates
  //Which correspond to reworked particle response branches in p4 and onwards
  //but like, this is target region specific...? also its just not working and my systematic universes are the same as my cv lol
  double ApplyCaloTuning(double calRecoilE) const{
    //Path to MParamFiles calibration files
    std::string pwd = "$MPARAMFILESROOT/data/Calibrations/energy_calib/CalorimetryTunings.txt";
    
    //THESE TWO ARE SUPERRRR DIFFERENT, how do I know which to use??? do I even need to worry about this? But how do I do response systematics without it??
    util::CaloCorrection Nu_Tracker(pwd.c_str(), "NukeCC_Nu_Tracker");
    //util::CaloCorrection Nu_Tracker(pwd.c_str(), "CCNuE_Nu_Tracker");

    return Nu_Tracker.eCorrection(calRecoilE/1000.)*1000.; //MeV
  }

  virtual double GetCalRecoilEnergy() const {
    double result = 0.0;
    for (const char* det : {"tracker","ecal","hcal","od","nucl"}) {
      std::string branchName = std::string("blob_recoil_E_") + det;
      result += GetDouble(branchName.c_str());
    }
    return result - GetLeakageCorrection();

    //this is what anezka's does...
    //return GetDouble("part_response_total_recoil_passive_allNonMuonClusters_id")+ GetDouble("part_response_total_recoil_passive_allNonMuonClusters_od"); // in MeV
  }
  virtual double GetCalRecoilEnergy2() const {
    //this is what anezka's does...
    return GetDouble("part_response_total_recoil_passive_allNonMuonClusters_id")+ GetDouble("part_response_total_recoil_passive_allNonMuonClusters_od"); // in MeV
  }

  double GetNonCalRecoilEnergy() const {
    return 0.0; //??? Anezka's is like this but is that correct?
  }

  
  // ========================================================================
  // Other truth variables
  // ========================================================================  
  int GetInteractionType() const {
    // it's mostly 1, 2, 3, lil bit of 8
    // and very lil bit of 4 
    // 1 = QE
    // 2 = Res
    // 3 = DIS
    // 4 = Coherent
    // 7 = nu + e elastic
    // 8 = MEC (meson exchange current?)
    return GetInt("mc_intType");
  }
  
  virtual bool IsNuEelastic() const {
    return GetInt("mc_intType") == 7;
  }
  
  //returns the index of the true primary proton (highest energy within kinematic signal constraints)
  //Only does the search once per event, after wards it caches it and that value gets used until CVUniverse->SetEntry is called again (aka we go to the next event)
  int GetHighestEnergySignalProtonIndex() const {
    if (m_protonIndexCached) { //this resets everytime cvUniv->SetEntry(i) is called
      return m_cachedProtonIndex; //if we've already found it for this event, don't need to do it again
    }    
    double highestEnergy = -9999;
    int index = -9999;
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    //fs particle energies in MeV
    std::vector<double> energies = GetVecDouble("mc_FSPartE");
    for (int i = 0; i < FSParticles.size(); i++){
      //So I can choose to count any protons, or only protons above our reco threshold. 
      //if (FSParticles[i] == 2212){
      if (FSParticles[i] == 2212 && energies[i] > highestEnergy){
	//require momentum between 450 and 1200 MeV/C , and angle under 70 degrees (same as other TKI analyses)
	
	double protonP = sqrt(pow(GetVecElem("mc_FSPartE", i),2) - pow(M_p, 2)); //in MeV/C
	ROOT::Math::XYZVector p(GetVecElem("mc_FSPartPx", i), GetVecElem("mc_FSPartPy", i), GetVecElem("mc_FSPartPz", i)); 
	ROOT::Math::RotationX r(-beam_tilt);
	double protonTheta = (r(p)).Theta()*(180/pi); //in degrees
	//if (protonP>450 && protonP<1200 && protonTheta<70){
	if (protonP>450 && protonP<1200 && (protonTheta<70 || protonTheta>110)){ //allow backwards protons
	  highestEnergy = energies[i];
	  index = i;
	}
      }
    }
    m_cachedProtonIndex = index;
    m_protonIndexCached = true;
    return index;
  }

  bool GetHasSignalFSProton() const{
    bool hasProton = false;
    int i = GetHighestEnergySignalProtonIndex();
    if (i>-1){
      hasProton = true;
    }
    return hasProton;
  }

  //For background categorization, how many events have protons that just don't land in my signal region requirements?
  bool GetHasAnyFSProton() const {
    bool hasProton = false;
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    for (int i = 0; i < FSParticles.size(); i++){
      if (FSParticles[i] == 2212){
	hasProton = true;
      }
    }
    return hasProton;
  }

  //in mm
  ROOT::Math::XYZTVector GetTrueVertex() const {
    ROOT::Math::XYZTVector result;
    result.SetCoordinates(GetVec<double>("mc_vtx").data());
    return result;
  }

  //which is which? double check
  virtual int GetCurrent() const { return GetInt("mc_current"); }

  // 12 = nu_e , 14 = nu_mu
  virtual int GetTruthNuPDG() const { return GetInt("mc_incoming"); }

  /*
  int GetHasFSProton() const {
    //mc_FSPartPDG
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    //fs particle energies in MeV
    std::vector<double> energies = GetVecDouble("mc_FSPartE");
    int hasProton = 0;
    for (int i = 0; i < FSParticles.size(); i++){
      //So I can choose to count any protons, or only protons above our reco threshold. 
      if (FSParticles[i] == 2212){
      //if (FSParticles[i] == 2212 && (energies[i]-938.272) > 80){

      //and then cause I might as well be an idiot, true proton 3 momentum. 
      //if (FSParticles[i] == 2212 && 
	hasProton = 1;
      }
    }
    //std::cout << "hasProton: " << hasProton<< std::endl;
    return hasProton;
  }
  */

  //Checks for pions & kaons & nothing else, maybe I'll change that but it's probably not a big deal
  bool GetHasFSMeson() const { 
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    bool hasMeson = false;
    for (int i = 0; i < FSParticles.size(); i++){
      //std::cout << "final state particle: " << i << ": " << FSParticles[i] << std::endl;
      if (abs(FSParticles[i]) == 211 || abs(FSParticles[i]) == 321 || abs(FSParticles[i]) == 311 || abs(FSParticles[i]) == 130 || abs(FSParticles[i]) == 111 || abs(FSParticles[i]) == 310){
	hasMeson = true;
      }
    }
    //std::cout << "hasMeson: " << hasMeson<< std::endl;
    return hasMeson;
  }

  //For background categorization, checks if the event has no mesons other than a single pi plus
  bool GetHasSinglePiPlus() const { 
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    bool hasSinglePiPlus = false;
    int nPiPlus = 0;
    int nPiMinus = 0;
    int nPi0 = 0;
    int nKaon = 0;
    for (int i = 0; i < FSParticles.size(); i++){
      if (FSParticles[i] == 211){ nPiPlus++; }
      else if (FSParticles[i] == -211){ nPiMinus++; }
      else if (FSParticles[i] == 111) { nPi0++; }
      else if (abs(FSParticles[i]) == 321 || abs(FSParticles[i]) == 311 || abs(FSParticles[i]) == 130 || abs(FSParticles[i]) == 310){ nKaon++; }
    }
    //if (nPiPlus==1 && nPiMinus==0 && nPi0==0 && nKaon==0) { hasSinglePiPlus = true; }
    if (nPiPlus==1 && nPiMinus==0 && nPi0==0) { hasSinglePiPlus = true; }
    return hasSinglePiPlus;
  }

  //For background categorization, checks if the event has no mesons other than a single pi minus
  bool GetHasSinglePiMinus() const { 
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    bool hasSinglePiMinus = false;
    int nPiPlus = 0;
    int nPiMinus = 0;
    int nPi0 = 0;
    int nKaon = 0;
    for (int i = 0; i < FSParticles.size(); i++){
      if (FSParticles[i] == 211){ nPiPlus++; }
      else if (FSParticles[i] == -211){ nPiMinus++; }
      else if (FSParticles[i] == 111) { nPi0++; }
      else if (abs(FSParticles[i]) == 321 || abs(FSParticles[i]) == 311 || abs(FSParticles[i]) == 130 || abs(FSParticles[i]) == 310){ nKaon++; }
    }
    //if (nPiPlus==0 && nPiMinus==1 && nPi0==0 && nKaon==0) { hasSinglePiMinus = true; }
    if (nPiPlus==0 && nPiMinus==1 && nPi0==0) { hasSinglePiMinus = true; }
    return hasSinglePiMinus;
  }

  //For background categorization, checks if the event has no mesons other than a single pi0 (lepton is still nu_e)
  bool GetHasSinglePiZero() const { 
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    bool hasSinglePiZero = false;
    int nPiPlus = 0;
    int nPiMinus = 0;
    int nPi0 = 0;
    int nKaon = 0;
    for (int i = 0; i < FSParticles.size(); i++){
      if (FSParticles[i] == 211){ nPiPlus++; }
      else if (FSParticles[i] == -211){ nPiMinus++; }
      else if (FSParticles[i] == 111) { nPi0++; }
      else if (abs(FSParticles[i]) == 321 || abs(FSParticles[i]) == 311 || abs(FSParticles[i]) == 130 || abs(FSParticles[i]) == 310){ nKaon++; }
    }
    //if (nPiPlus==0 && nPiMinus==0 && nPi0==1 && nKaon==0) { hasSinglePiZero = true; }
    if (nPiPlus==0 && nPiMinus==0 && nPi0==1) { hasSinglePiZero = true; }
    return hasSinglePiZero;
  }

  //For background categorization, checks if the event has multiple pions of any type
  bool GetHasMultiplePions() const { 
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    bool hasMultiPi = false;
    int nPiPlus = 0;
    int nPiMinus = 0;
    int nPi0 = 0;
    int nKaon = 0;
    for (int i = 0; i < FSParticles.size(); i++){
      if (FSParticles[i] == 211){ nPiPlus++; }
      else if (FSParticles[i] == -211){ nPiMinus++; }
      else if (FSParticles[i] == 111) { nPi0++; }
      else if (abs(FSParticles[i]) == 321 || abs(FSParticles[i]) == 311 || abs(FSParticles[i]) == 130 || abs(FSParticles[i]) == 310){ nKaon++; }
    }
    //if ( (nPiPlus + nPiMinus + nPi0)>1 && nKaon==0) { hasMultiPi = true; }
    if ( (nPiPlus + nPiMinus + nPi0)>1) { hasMultiPi = true; }
    return hasMultiPi;
  }

  bool GetHasKaons() const { 
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    bool hasKaons = false;
    for (int i = 0; i < FSParticles.size(); i++){
      if (abs(FSParticles[i]) == 321 || abs(FSParticles[i]) == 311 || abs(FSParticles[i]) == 130 || abs(FSParticles[i]) == 310){ hasKaons = true; }
    }
    return hasKaons;
  }
  
  //checks for photons above 10 MeV in final state
  bool GetHasFSPhoton() const {
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    std::vector<double> energies = GetVecDouble("mc_FSPartE");
    bool hasPhoton = false;
    for (int i = 0; i < FSParticles.size(); i++){
      //std::cout << "final state particle: " << i << ": " << FSParticles[i] << std::endl;
      if (abs(FSParticles[i]) == 22 && energies[i] > 10){
	hasPhoton = true;
      }
    }
    //std::cout << "hasPhoton: " << hasPhoton<< std::endl;
    return hasPhoton;
  }

  bool GetHasFSPi0() const { 
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    int hasPi0 = false;
    for (int i = 0; i < FSParticles.size(); i++){
      //std::cout << "final state particle: " << i << ": " << FSParticles[i] << std::endl;      
      if (abs(FSParticles[i]) == 111){ 
	hasPi0 = true;
      }
    }
    //std::cout << "hasPi0: " << hasPi0<< std::endl;
    return hasPi0;
  }

  // ========================================================================
  // Here are getters which return my analysis variables but in GeV for plotting
  // closure tests, etc. All cross sections are reported in terms of GeV it seems,
  // but I want to keep all of the output of my normal functions in MeV for
  // calculations. Before I was having trouble keeping track of MeV/GeV for every
  // single function...
  // ========================================================================
  double GetElectronEnergyGeV() const {
    return GetElectronEnergy()/1000.;
  }
  double GetElectronEnergyTrueGeV() const {
    return GetElectronEnergyTrue()/1000.;
  }

  double GetEavailGeV() const {
    return GetEavail()/1000.;
  }
  double GetEavailTrueGeV() const {
    return GetEavailTrue()/1000.;
  }

  double GetEnuGeV() const {
    return GetEnu()/1000.;
  }
  double GetEnuTrueGeV() const {
    return GetEnuTrue()/1000.;
  }

  double GetElectronPtGeV() const {
    return GetElectronPt()/1000.;
  }
  double GetElectronPtTrueGeV() const {
    return GetElectronPtTrue()/1000.;
  }

  double GetElectronPParallelGeV() const {
    return GetElectronPParallel()/1000.;
  }
  double GetElectronPParallelTrueGeV() const {
    return GetElectronPParallelTrue()/1000.;
  }

  double GetProtonPGeV() const {
    return GetProtonP()/1000.;
  }
  double GetProtonPTrueGeV() const {
    return GetProtonPTrue()/1000.;
  }

  double GetProtonPtGeV() const {
    return GetProtonPt()/1000.;
  }
  double GetProtonPtTrueGeV() const {
    return GetProtonPtTrue()/1000.;
  }

  double GetProtonTGeV() const {
    return GetProtonT()/1000.;
  }
  double GetProtonTTrueGeV() const {
    return GetProtonTTrue()/1000.;
  }

  double GetDeltaPtGeV() const {
    return GetDeltaPt()/1000.;
  }
  double GetDeltaPtTrueGeV() const {
    return GetDeltaPtTrue()/1000.;
  }

  double GetDeltaPtXGeV() const {
    return GetDeltaPtX()/1000.;
  }
  double GetDeltaPtXTrueGeV() const {
    return GetDeltaPtXTrue()/1000.;
  }

  double GetDeltaPtYGeV() const {
    return GetDeltaPtY()/1000.;
  }
  double GetDeltaPtYTrueGeV() const {
    return GetDeltaPtYTrue()/1000.;
  }

  double GetDeltaPlGeV() const {
    return GetDeltaPl()/1000.;
  }
  double GetDeltaPlTrueGeV() const {
    return GetDeltaPlTrue()/1000.;
  }

  double GetPnGeV() const {
    return GetPn()/1000.;
  }
  double GetPnTrueGeV() const {
    return GetPnTrue()/1000.;
  }



  //This is for if I want to have a truth variable, i.e. I want to plot true proton momentum vs true proton angle
  //I need to put something for the reco var, and I don't want to write the function for reco proton momentum and angle yet.
  //
  //Update jan 2025 - i actually don't think I need this, the VariableBase() constructor might do this automatically if you
  //don't provide a reco OR truth getter unction? need to double check tho
  double GetDummyVar() const {
    return -999; 
  }
  
  // ========================================================================
  // START OF PRE EXISTING GETTERS THAT I DON'T USE
  // mostly muon kinematic stuff, didn't feel like getting rid of it
  // also it's all in GeV while all of mine are in MeV (cause that's what the python code has), so that's annoying if I ever need them...
  // ========================================================================  

  // Quantities only needed for cuts
  // Although unlikely, in principle these quanties could be shifted by a
  // systematic. And when they are, they'll only be shifted correctly if we
  // write these accessor functions.

  //Muon kinematics
  double GetMuonPT() const //GeV/c
  {
    return GetPmu()/1000. * sin(GetThetamu());
  }

  double GetMuonPz() const //GeV/c
  {
    return GetPmu()/1000. * cos(GetThetamu());
  }

  double GetMuonPTTrue() const //GeV/c
  {
    return GetPlepTrue()/1000. * sin(GetThetalepTrue());
  }

  double GetMuonPzTrue() const //GeV/c
  {
    return GetPlepTrue()/1000. * cos(GetThetalepTrue());
  }

  double GetEmuGeV() const //GeV
  {
    return GetEmu()/1000.;
  }

  double GetElepTrueGeV() const //GeV
  {
    return GetElepTrue()/1000.;
  }

  int GetTargetNucleon() const {
    return GetInt("mc_targetNucleon");
  }
  
  double GetBjorkenXTrue() const {
    return GetDouble("mc_Bjorkenx");
  }

  double GetBjorkenYTrue() const {
    return GetDouble("mc_Bjorkeny");
  }

  virtual bool IsMinosMatchMuon() const {
    return GetInt("has_interaction_vertex") == 1;
  }  
  
  virtual double GetMuonQP() const {
    return GetDouble((GetAnaToolName() + "_minos_trk_qp").c_str());
  }

  //Some functions to match CCQENuInclusive treatment of DIS weighting. Name matches same Dan area as before.
  virtual double GetTrueExperimentersQ2() const {
    double Enu = GetEnuTrue(); //MeV
    double Emu = GetElepTrue(); //MeV
    double thetaMu = GetThetalepTrue();
    return 4.0*Enu*Emu*pow(sin(thetaMu/2.0),2.0);//MeV^2
  }

  virtual double CalcTrueExperimentersQ2(double Enu, double Emu, double thetaMu) const{
    return 4.0*Enu*Emu*pow(sin(thetaMu/2.0),2.0);//MeV^2
  }

  virtual double GetTrueExperimentersW() const {
    double nuclMass = M_nucleon;
    int struckNucl = GetTargetNucleon();
    if (struckNucl == PDG_n){
      nuclMass=M_n;
    }
    else if (struckNucl == PDG_p){
      nuclMass=M_p;
    }
    double Enu = GetEnuTrue();
    double Emu = GetElepTrue();
    double thetaMu = GetThetalepTrue();
    double Q2 = CalcTrueExperimentersQ2(Enu, Emu, thetaMu);
    return TMath::Sqrt(pow(nuclMass,2) + 2.0*(Enu-Emu)*nuclMass - Q2);
  }

  //Still needed for some systematics to compile, but shouldn't be used for reweighting anymore.
  protected:
  #include "PlotUtils/WeightFunctions.h" // Get*Weight
};

#endif
