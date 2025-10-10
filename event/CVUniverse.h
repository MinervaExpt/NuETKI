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
  static constexpr double M_e = 0.51099895; //MeV
  static constexpr double M_pi = 139.57039; //MeV
  static constexpr double M_nucleon = (1.5*M_n+M_p)/2.5;

  static constexpr int PDG_n = 2112;
  static constexpr int PDG_p = 2212;

  static constexpr double pi = 3.141592653589793;

  //what is this actual number supposed to be? I'm using a hardcoded value of 3.3 right now
  //static constexpr double beam_tilt = 3.2627; //beam tilt downward with respect to z axis, in degrees

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

  //total energy including rest mass in GeV
  double GetElectronEnergy() const {
    return GetVecElem("prong_part_E", 0, 3)/1000.;
  }

  //Copied the overall format here from Ryan/Hang/Sarah's CCNue kinematics calculator, there may be a better way to do this?
  //Also extremely suspect for right now
  //returns in rad
  double GetElectronTheta() const { 
    std::vector<std::vector<double> > electronp = GetVecOfVecDouble("prong_part_E");
    ROOT::Math::XYZVector p(electronp[0][0], electronp[0][1], electronp[0][2]);
    ROOT::Math::RotationX r(-3.3 * (pi / 180.));
    return (r(p)).Theta();
  }

  //mostly for plotting
  double GetElectronThetaDeg() const {
    double rad = GetElectronTheta();
    double deg = rad*180/pi;
    return deg;
  }

  //Electron kinetic energy (T) in GeV
  //Actually I think this is wrong, the E that comes from prong_part_E doesn't have electron mass added...
  //probably cause electron rest mass is basically negligible (0.5MeV, or 0.0005 GeV, or less than the average hit energy...)
  //But either way, makes no sense to subtract it back off if it was never added in the first place
  double GetElectronT() const{
    double leptonT = GetElectronEnergy() - M_e/1000.; //E_e from the function is in GeV, M_e is in MeV
    if (leptonT > 0){
      return leptonT;
    }
    else return 0; 
    //return sqrt(std::max(0, (pow(GetElectronEnergy(), 2) - M_lep_sqr)));
  }
  
  //Electron momentum in GeV
  //And actually same story as above here, this doesn't really make sense. It should just be prong_part_E straight up
  //which again makes sense, for relativistic particles E = p (pretend its squiggly)
  //Doesn't really hurt to do it, but the value is gonna be the same basically.
  double GetElectronP() const{
    double M_lep_sqr = pow(M_e, 2) / pow(10, 6);  //over 10^6 to convert to GeV^2
    double leptonP = pow(GetElectronEnergy(), 2)-M_lep_sqr;
    if (leptonP > 0){
      return sqrt(leptonP);
    }
    else return 0; 
    //return sqrt(std::max(0, (pow(GetElectronEnergy(), 2) - M_lep_sqr)));
  }

  //Transverse momentum of the lepton (electron for me), in GeV
  double GetElectronPt() const{
    return (GetElectronP() * std::sin(GetElectronTheta()));
  }

  //Longitudinal momentum of the lepton (electron for me), in GeV
  double GetElectronPParallel() const{
    return (GetElectronP() * std::cos(GetElectronTheta()));
  }

  //MasterAnaDev_proton_theta, in rad
  double GetProtonTheta() const {
    double protonTheta = GetDouble("MasterAnaDev_proton_theta");
    return protonTheta;

    //TO DO: THESE TWO SETS OF QUANTITIES ARE NOT EXACTLY THE SAME. they almost are, but not quite?
    //Find out why, which one to use, and does it matter? Gonna use MAD one for now

    //double protonTheta1 = GetDouble("MasterAnaDev_proton_theta");
    //double protonPhi1 = GetDouble("MasterAnaDev_proton_phi");

    //double protonTheta2 = (r(protonP_vec)).Theta();
    //double protonPhi2 = (r(protonP_vec)).Phi(); 

  }

  //same thing but in rad, useful for plotting. Also if there's no proton candidate don't convert it
  double GetProtonThetaDeg() const {
    double protonTheta = GetDouble("MasterAnaDev_proton_theta");
    //std::cout << "True (highest) Proton Momentum: " << protonP << std::end; 
    if (protonTheta > -9999){
      return protonTheta * (180/pi);
    }
    else {
      return protonTheta;
    }
  }

  //Kinetic energy (T) of ONLY the primary proton candidate
  double GetProtonT() const {
    double protonT = GetDouble("MasterAnaDev_proton_calib_energy"); //TO DO !! MAKE SURE THIS IS ACTUALLY T AND NOT TOTAL E
    if (protonT>-9999){ //aka is there a proton candidate
      return protonT/1000.;
    }
    else{
      return protonT;
    }
  }

  //Total Kinetic energy (sigma T) of all proton candidates
  double GetSumProtonT() const {
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
  
  //MasterAnaDev_proton_P_fromdEdx, in GeV
  double GetProtonP() const {
    double protonP = GetDouble("MasterAnaDev_proton_P_fromdEdx");
    if (protonP>-9999){ //aka is there a proton candidate
      return protonP/1000.;
    }
    else{
      return protonP;
    }
  }

  //proton transverse momentum, in GeV
  double GetProtonPt() const{
    return (GetProtonP() * std::sin(GetProtonTheta()));
  }

  //proton longitudinal momentum, in GeV
  double GetProtonPParallel() const{
    return (GetProtonP() * std::cos(GetProtonTheta()));
  }

  //Hadronic available energy (no neutrons) in GeV
  double GetEavail() const {
    double E = GetDouble("blob_recoil_E_tracker") + GetDouble("blob_recoil_E_ecal");
    //Check this correction by the electron energy... This is not what Hang does, and in fact it's different for data/MC
    double Eavail = (E * 1.17 - ((0.008 * GetVecElem("prong_part_E", 0, 3)) + 5))/1000.; 
    //std::cout << "Reco Eavail: " << Eavail << std::endl;
    return Eavail;
  }

  //This is the same e-avail as above, but with all tracked proton calibrated kinetic energy subtracted off (
  //this is a completely reco variable by the way, in case I get lost lol
  double GetModifiedEavail() const {
    double Eavail = GetEavail()*1000.; //because my function returns eavail in GeV and all of these branches are in MeV :)))))

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
    return ModEavail/1000.;
  }

  //in GeV, also note I changed the truth version in opt/include/PlotUtils/TruthFunctions.h to match (aka put it in GeV). So if I reinstall or something check that.
  double GetEnu() const {
    double Enu = GetEavail() + GetElectronEnergy();
    //std::cout << "Reco Enu: " << Enu << std::endl;
    return Enu;
  }
  
  //Returns a root XYZVector object containing the lepton (electron) transverse 3 momentum
  ROOT::Math::XYZVector GetElectronPtVec() const{
    std::vector<std::vector<double> > electronP = GetVecOfVecDouble("prong_part_E");
    ROOT::Math::XYZVector electronP_vec(electronP[0][0]/1000., electronP[0][1]/1000., electronP[0][2]/1000.);

    ROOT::Math::RotationX r(-3.3 * (pi / 180.)); //This object represents a slight rotation about the x axis so that the theta we get is wrt to the beam direction and not the z axis, which points down to the ground 3.3 degrees. we can then apply r to the two vectors above
    ROOT::Math::RotationX r2(3.3 * (pi / 180.)); //the reverse rotation to get back to lab coords (3.3 is positive now)

    //perform the rotation, basically transforming to new beam based coord system
    ROOT::Math::XYZVector electronPt_vec = r(electronP_vec);

    //Now set z to zero in the beam frame to get only the transverse components, then rotate back
    electronPt_vec.SetZ(0); 
    electronPt_vec = r2(electronPt_vec);

    return electronPt_vec;
  }
  
  //Returns a root XYZVector object containing the lepton (electron) transverse 3 momentum
  ROOT::Math::XYZVector GetProtonPtVec() const{
    ROOT::Math::XYZVector protonP_vec(GetDouble("MasterAnaDev_proton_Px_fromdEdx")/1000., GetDouble("MasterAnaDev_proton_Py_fromdEdx")/1000., GetDouble("MasterAnaDev_proton_Pz_fromdEdx")/1000.);
    ROOT::Math::RotationX r(-3.3 * (pi / 180.)); 
    ROOT::Math::RotationX r2(3.3 * (pi / 180.));

    ROOT::Math::XYZVector protonPt_vec = r(protonP_vec);
    protonPt_vec.SetZ(0);
    protonPt_vec = r2(protonPt_vec);

    return protonPt_vec;
  }
  
  //Returns a root XYZVector object containing the sum of the proton transverse 3 momentum and the lepton (electron) transverse 3 momentum
  //which is then used to calculate TKI variables
  //Remember: z direction != beam direction so transverse doesn't exactly mean z components are zero, although they should be small
  ROOT::Math::XYZVector GetDeltaPtVec() const{
    ROOT::Math::XYZVector electronPt_vec = GetElectronPtVec();
    ROOT::Math::XYZVector protonPt_vec = GetProtonPtVec();
        
    //ROOT::Math::XYZVector deltaP_vec = electronP_vec + protonP_vec; //sum of the full 3 momenta. Not sure if I need this so commenting it out for now
    ROOT::Math::XYZVector deltaPt_vec = electronPt_vec + protonPt_vec;
    
    //std::cout << "delta P total (kinda useless?): " << sqrt(deltaP_vec.Mag2()) << std::endl;
    //std::cout << "delta Ptx: " << deltaPt_vec.X() << std::endl;
    //std::cout << "delta Pty: " << deltaPt_vec.Y() << std::endl;
    //std::cout << "delta Pt: " << sqrt(deltaPt_vec.Mag2()) << std::endl;
    //std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n" << std::endl;
    return deltaPt_vec;
  }
  
  //delta pT (magnitude of the vector), in GeV
  double GetDeltaPt() const{
    ROOT::Math::XYZVector deltaPt_vec = GetDeltaPtVec();
    //std::cout << "delta Pt RECO: " << sqrt(deltaPt_vec.Mag2()) << std::endl;
    return sqrt(deltaPt_vec.Mag2()); 
  }

  //Magnitude of the component of delta pT orthogonal to electron pT (GeV)
  double GetDeltaPtX() const{
    ROOT::Math::XYZVector electronPt_vec = GetElectronPtVec();
    ROOT::Math::XYZVector deltaPt_vec = GetDeltaPtVec();

    ROOT::Math::RotationX r(-3.3 * (pi / 180.));

    //Rotate both back to beam frame
    ROOT::Math::XYZVector electron_pt_in_beam_frame = r(electronPt_vec);
    ROOT::Math::XYZVector delta_pt_in_beam_frame = r(deltaPt_vec);

    //expression for the scalar rejection, aka perpendicular dot product
    double num = delta_pt_in_beam_frame.Y()*electron_pt_in_beam_frame.X() - delta_pt_in_beam_frame.X()*electron_pt_in_beam_frame.Y();
    double denom = sqrt(electron_pt_in_beam_frame.Mag2());
    return num/denom;
  }

  //Magnitude of the component of  delta pT parallel to electron pT (GeV)
  double GetDeltaPtY() const{
    ROOT::Math::XYZVector electronPt_vec = GetElectronPtVec();
    ROOT::Math::XYZVector deltaPt_vec = GetDeltaPtVec();

    ROOT::Math::RotationX r(-3.3 * (pi / 180.));

    //Rotate both back to beam frame
    ROOT::Math::XYZVector electron_pt_in_beam_frame = r(electronPt_vec);
    ROOT::Math::XYZVector delta_pt_in_beam_frame = r(deltaPt_vec);

    //std::cout << "\nelectron pT vec: ( " << electron_pt_in_beam_frame.X() << " , " << electron_pt_in_beam_frame.Y() << " , " << electron_pt_in_beam_frame.Z() << " ) " << std::endl;
    //std::cout << "delta pT vec: ( " << delta_pt_in_beam_frame.X() << " , " << delta_pt_in_beam_frame.Y() << " , " << delta_pt_in_beam_frame.Z() << " ) " << std::endl;
    //std::cout << "dot product of the two: " << delta_pt_in_beam_frame.Dot(electron_pt_in_beam_frame) << std::endl;
    //std::cout << "magnitude (norm) of e_pT vec: " << sqrt(electron_pt_in_beam_frame.Mag2()) << std::endl;
    
    //expression for the scalar projection
    double num = delta_pt_in_beam_frame.Dot(electron_pt_in_beam_frame);
    double denom = sqrt(electron_pt_in_beam_frame.Mag2());

    //std::cout << "final result: " << num/denom << "\n" << std::endl;
    return -num/denom;
  }
  
  //Alpha_t, the TKI boosting angle. it's the angle between inverted electron pT and delta pT
  double GetAlphaT() const{
    ROOT::Math::XYZVector electronPt_vec = GetElectronPtVec();
    ROOT::Math::XYZVector deltaPt_vec = GetDeltaPtVec();

    //angle between the two vectors
    double numerator = ((-1*electronPt_vec).Dot(deltaPt_vec));
    double denominator = ( sqrt(electronPt_vec.Mag2()) * sqrt(deltaPt_vec.Mag2()) );
    double alpha = std::acos(numerator/denominator);

    //std::cout << "boosting angle alpha: " << alpha * 180/pi << std::endl;
    return alpha * 180/pi;  //return alpha in degrees cause that's how I've set up my bins for now
  }

  //Phi_t 
  double GetPhiT() const{
    ROOT::Math::XYZVector electronPt_vec = GetElectronPtVec();
    ROOT::Math::XYZVector protonPt_vec = GetProtonPtVec();

    //angle between the two vectors 
    double numerator = ((-1*electronPt_vec).Dot(protonPt_vec));
    double denominator = ( sqrt(electronPt_vec.Mag2()) * sqrt(protonPt_vec.Mag2()) );
    double phi = std::acos(numerator/denominator);

    //std::cout << "TKI phi: " << phi * 180/pi << std::endl;
    return phi * 180/pi;  //return in degrees cause that's how I've set up my bins for now
  }
  
  //Delta P parallel, or the sum of longitudinal (wrt beam direction) momentum of proton & electron.
  //Because this is 1D only can just sum the two values I already have
  //Don't need to mess about with vectors, and those functions already return neg & pos w.r.t. beam direction
  double GetSumPParallel() const{    
    return GetElectronPParallel() + GetProtonPParallel();
  }

  //This is the variable reported by Jeffrey & Xianguo, which is NOT simply delta P parallel
  //Delta P Parallel is the sum of lepton + proton parallel, and then this is the difference between that and initial neutrino longitudinal momentum
  //Which is estimated by the magnitude of the recoil momentum of the struck nucleus or something
  //The delta here is technically the same delta as above, but in the the transverse case the neutrino's contribution is 0,
  //so it simplifies to lepton + proton. That's not the case here. 
  double GetDeltaP_L() const{    
    return 0;
  }

  //This variable represents the initial struck nucleon momentum, determined by delta Pt and delta Pl.
  double GetP_n() const{
  }
  
  virtual double GetQ2Reco() const{
    return GetDouble("qsquared_recoil");
  }

  //GetRecoilE is designed to match the NSF validation suite
  virtual double GetRecoilE() const {
    return GetVecElem("recoil_summed_energy", 0);
  }
  
  //Definitely possible mistakes in here...
  virtual double Getq3() const{
    double eavail = GetEavail()/pow(10,3);
    double q2 = GetQ2Reco() / pow(10,6);
    double q3mec = sqrt(eavail*eavail + q2);
    return q3mec;
  }
  
  // ========================================================================
  // Truth kinematics
  // ========================================================================  

  //in GeV
  double GetElectronEnergyTrue() const {
    double out = GetVecElem("mc_primFSLepton", 3)/1000.;
    //std::cout << "True Lepton Energy: " << out << std::endl;
    return out;
  }

  //in rad
  double GetElectronThetaTrue() const {
    std::vector<double> electronp = GetVecDouble("mc_primFSLepton");
    ROOT::Math::XYZVector p(electronp[0], electronp[1], electronp[2]);   
    ROOT::Math::RotationX r(-3.3 * (pi / 180.)); //This is a slight rotation so that the theta we get is wrt to the beam direction and not the z axis, which points down to the ground 3.3 degrees. 
    return (r(p)).Theta();    
  }

  //in deg
  double GetElectronThetaDegTrue() const {
    return GetElectronThetaTrue() * (180./pi);
  }
  
  //True lepton pT (in GeV)
  double GetElectronPtTrue() const {
    double leptonP = sqrt(pow(GetVecElem("mc_primFSLepton", 3), 2) - pow(M_e, 2));
    double out = (leptonP * std::sin(GetElectronThetaTrue()))/1000.;
    //std::cout << "True lepton pT: " << out << std::endl;
    return out;
  }

  //True lepton p parallel (in GeV)
  double GetElectronPParallelTrue() const {
    double leptonP = sqrt(pow(GetVecElem("mc_primFSLepton", 3), 2) - pow(M_e, 2));
    double out = (leptonP * std::cos(GetElectronThetaTrue()))/1000.;
    //std::cout << "True lepton pT: " << out << std::endl;
    return out;
  }

  //gets highest true proton kinetic energy in GeV, -999 (- M_p for now, will fix) if no true protons
  double GetProtonKETrue() const {
    int highestProtonIndex = GetHighestEnergySignalProtonIndex();
    if (highestProtonIndex > -1){
      double out = (GetVecElem("mc_FSPartE", highestProtonIndex) - M_p)/1000.;
      //std::cout << "True (highest) Proton KE: " << out << std::endl;
      return out;
    }
    else {
      return -999; 
    }
  }

  //proton theta (spherical coords, in deg) wrt beam direction, -999 if no true protons
  //in rad
  double GetProtonThetaTrue() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      ROOT::Math::XYZVector p(GetVecElem("mc_FSPartPx", i), GetVecElem("mc_FSPartPy", i), GetVecElem("mc_FSPartPz", i));
      ROOT::Math::RotationX r(-3.3 * (pi / 180.));
      double out = (r(p)).Theta();
      //std::cout << "True (highest energy) proton Theta: " << out << std::endl;
      return out;
    }
    else {
      return -999; 
    }
  }

  //to deg
  double GetProtonThetaDegTrue() const {
    double proton_theta = GetProtonThetaTrue();
    if (proton_theta == -999) {
      return -999;
    } else {
      return proton_theta * (180./pi);
    }
  }


  //Returns momentum of the highest energy true final state proton (in GeV)
  double GetProtonPTrue() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      double protonP = sqrt(pow(GetVecElem("mc_FSPartE", i),2) - pow(M_p, 2));
      //std::cout << "True (highest) Proton Momentum: " << protonP << std::end; 
      return protonP/1000.;
    }
    else {
      return -999;
    }
  }

  //true proton pT (in GeV) of highest energy proton
  double GetProtonPtTrue() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      //returns momentum in GeV, so convert to MeV which is what M_p is in, should be good
      double protonP = GetProtonPTrue()*1000.;
      double protonAngleRad = GetProtonThetaTrue();
      
      double out = (protonP * std::sin(protonAngleRad))/1000.;	
      //std::cout << "True (highest energy) proton pT: " << out << std::endl;
      return out;
    }
    else {
      return -999; 
    }
  }

  //true proton pT (in GeV) of highest energy proton
  double GetProtonPParallelTrue() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      //returns momentum in GeV, so convert to MeV which is what M_p is in, should be good
      double protonP = GetProtonPTrue()*1000.;
      double protonAngleRad = GetProtonThetaTrue();
      
      double out = (protonP * std::cos(protonAngleRad))/1000.;	
      //std::cout << "True (highest energy) proton pT: " << out << std::endl;
      return out;
    }
    else {
      return -999; 
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
    return T_p/1000.;
  }

  //Angle (in deg) between true primary electron and true highest energy proton in rad, -999 if no protons
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
      return -999; 
    }
  }

  //E_avail as defined in truth, basically summing all sources of visible energy
  double GetEavailTrue() const {
    double T_p = 0; //sum of proton kinetic energies
    double T_pi = 0; //sum of charged pion kinetic energies
    double E_pi0 = 0; //sume of neutral pion total energies
    double E_s = 0; //Sum of (strange baryon energy - proton mass)
    double E_sbar = 0; //Sum of (anti baryon energy + proton mass)
    double E_other = 0; //sum of total energy of any other particles, NOT INCLUDING NEUTRONS

    int nFSPart = GetInt("mc_nFSPart");
    std::vector<double> FSPartPx = GetVecDouble("mc_FSPartPx");
    std::vector<double> FSPartPy = GetVecDouble("mc_FSPartPy");
    std::vector<double> FSPartPz = GetVecDouble("mc_FSPartPz");
    std::vector<double> FSPartE = GetVecDouble("mc_FSPartE");
    std::vector<int> FSPartPDG = GetVecInt("mc_FSPartPDG");

    
    for (int i=0; i<nFSPart; i++){
      if (abs(FSPartPDG[i]) == 11 || abs(FSPartPDG[i]) == 13 || FSPartPDG[i] == 2112 || FSPartPDG[i] > 1000000000){
	continue; //don't want to count leptons, neutrons, or nuclear remnants
      } else if (FSPartPDG[i] == 2212){ 
	T_p += FSPartE[i] - M_p; //add only kinetic (not total) energy for protons
      } else if (abs(FSPartPDG[i]) == 211){ 
	T_pi += FSPartE[i] - M_pi; //add only kinetic (not total) energy for charged pions
      } else if (abs(FSPartPDG[i]) == 111){ 
	E_pi0 += FSPartE[i]; //add pi0 total energy
      } else if (FSPartPDG[i] > 2000){
	E_s += FSPartE[i] - M_p; //add total energy - proton mass for strange baryons (why?)
      } else if (FSPartPDG[i] < 2000){
	E_sbar += FSPartE[i] + M_p; //add total energy + proton mass for strange antibaryons (why???)
      } else {
	E_other += FSPartE[i]; //add total energy for anything else (mostly gammas, kaons i think?)
      }
    }

    double E_avail = T_p + T_pi + E_pi0 + E_s + E_sbar + E_other;
    return E_avail/1000.;
  }

  double GetDeltaPtTrue() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      ROOT::Math::XYZVector electronP_vec(GetVecElem("mc_primFSLepton", 0), GetVecElem("mc_primFSLepton", 1), GetVecElem("mc_primFSLepton", 2));
      ROOT::Math::XYZVector protonP_vec(GetVecElem("mc_FSPartPx",i), GetVecElem("mc_FSPartPy",i), GetVecElem("mc_FSPartPz",i));

      ROOT::Math::RotationX r(-3.3 * (pi / 180.)); //rotation into beam frame
      ROOT::Math::RotationX r2(3.3 * (pi / 180.)); //rotation back into lab frame

      ROOT::Math::XYZVector electronPt_vec = r(electronP_vec);
      ROOT::Math::XYZVector protonPt_vec = r(protonP_vec);
      electronPt_vec.SetZ(0);
      protonPt_vec.SetZ(0);

      //do I need to rotate back for this? I don't think so but double check
      ROOT::Math::XYZVector deltaPt_vec = electronPt_vec + protonPt_vec;      
      return sqrt(deltaPt_vec.Mag2())/1000.;
    }
    else {
      return -999; //what do I return for delta pt true if no true protons?
    }
  }

  double GetDeltaPtXTrue() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      ROOT::Math::XYZVector electronP_vec(GetVecElem("mc_primFSLepton", 0), GetVecElem("mc_primFSLepton", 1), GetVecElem("mc_primFSLepton", 2));
      ROOT::Math::XYZVector protonP_vec(GetVecElem("mc_FSPartPx",i), GetVecElem("mc_FSPartPy",i), GetVecElem("mc_FSPartPz",i));

      ROOT::Math::RotationX r(-3.3 * (pi / 180.)); //rotation into beam frame
      ROOT::Math::RotationX r2(3.3 * (pi / 180.)); //rotation back into lab frame

      ROOT::Math::XYZVector electronPt_vec = r(electronP_vec);
      ROOT::Math::XYZVector protonPt_vec = r(protonP_vec);
      electronPt_vec.SetZ(0);
      protonPt_vec.SetZ(0);

      //Still in beam frame
      ROOT::Math::XYZVector deltaPt_vec = electronPt_vec + protonPt_vec;      
      double num = deltaPt_vec.Y()*electronPt_vec.X() - deltaPt_vec.X()*electronPt_vec.Y();
      double denom = sqrt(electronPt_vec.Mag2());
      return (num/denom)/1000.;
    }
    else {
      return -999; //what do I return for delta pt true if no true protons?
    }
  }

  double GetDeltaPtYTrue() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      ROOT::Math::XYZVector electronP_vec(GetVecElem("mc_primFSLepton", 0), GetVecElem("mc_primFSLepton", 1), GetVecElem("mc_primFSLepton", 2));
      ROOT::Math::XYZVector protonP_vec(GetVecElem("mc_FSPartPx",i), GetVecElem("mc_FSPartPy",i), GetVecElem("mc_FSPartPz",i));

      ROOT::Math::RotationX r(-3.3 * (pi / 180.)); //rotation into beam frame
      ROOT::Math::RotationX r2(3.3 * (pi / 180.)); //rotation back into lab frame

      ROOT::Math::XYZVector electronPt_vec = r(electronP_vec);
      ROOT::Math::XYZVector protonPt_vec = r(protonP_vec);
      electronPt_vec.SetZ(0);
      protonPt_vec.SetZ(0);

      //Still in beam frame
      ROOT::Math::XYZVector deltaPt_vec = electronPt_vec + protonPt_vec;      

      double num = deltaPt_vec.Dot(electronPt_vec);
      double denom = sqrt(electronPt_vec.Mag2());
      //std::cout << "TRUE delta PtY: " << (num/denom)/1000. << std::endl;
      return -(num/denom)/1000.;
    }
    else {
      return -999; //what do I return for delta pt true if no true protons?
    }
  }

  double GetAlphaTTrue() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      ROOT::Math::XYZVector electronP_vec(GetVecElem("mc_primFSLepton", 0), GetVecElem("mc_primFSLepton", 1), GetVecElem("mc_primFSLepton", 2));
      ROOT::Math::XYZVector protonP_vec(GetVecElem("mc_FSPartPx",i), GetVecElem("mc_FSPartPy",i), GetVecElem("mc_FSPartPz",i));

      ROOT::Math::RotationX r(-3.3 * (pi / 180.)); //rotation into beam frame
      ROOT::Math::RotationX r2(3.3 * (pi / 180.)); //rotation back into lab frame
      
      ROOT::Math::XYZVector electronPt_vec = r(electronP_vec);
      ROOT::Math::XYZVector protonPt_vec = r(protonP_vec);
      electronPt_vec.SetZ(0);
      protonPt_vec.SetZ(0);

      //Still in beam frame
      ROOT::Math::XYZVector deltaPt_vec = electronPt_vec + protonPt_vec;

      double numerator = ((-1*electronPt_vec).Dot(deltaPt_vec));
      double denominator = ( sqrt(electronPt_vec.Mag2()) * sqrt(deltaPt_vec.Mag2()) );
      double alpha = std::acos(numerator/denominator);
      
      return alpha * 180/pi;
    }
    else {
      return -999;
    }
  }

  double GetPhiTTrue() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      ROOT::Math::XYZVector electronP_vec(GetVecElem("mc_primFSLepton", 0), GetVecElem("mc_primFSLepton", 1), GetVecElem("mc_primFSLepton", 2));
      ROOT::Math::XYZVector protonP_vec(GetVecElem("mc_FSPartPx",i), GetVecElem("mc_FSPartPy",i), GetVecElem("mc_FSPartPz",i));

      ROOT::Math::RotationX r(-3.3 * (pi / 180.)); //rotation into beam frame
      ROOT::Math::RotationX r2(3.3 * (pi / 180.)); //rotation back into lab frame
      
      ROOT::Math::XYZVector electronPt_vec = r(electronP_vec);
      ROOT::Math::XYZVector protonPt_vec = r(protonP_vec);
      electronPt_vec.SetZ(0);
      protonPt_vec.SetZ(0);

      //Still in beam frame
      ROOT::Math::XYZVector deltaPt_vec = electronPt_vec + protonPt_vec;

      double numerator = ((-1*electronPt_vec).Dot(protonPt_vec));
      double denominator = ( sqrt(electronPt_vec.Mag2()) * sqrt(protonPt_vec.Mag2()) );
      double phi = std::acos(numerator/denominator);
      
      return phi * 180/pi;
    }
    else {
      return -999;
    }
  }
  
  // ========================================================================
  // Other reco variables
  // Mostly used for evaluating cuts, mostly straight from the tuple (other than ESC)
  // ========================================================================  

  int GetImprovedNMichel() const {
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
  double GetELepSin2Theta() const {
    //double etheta = GetElectronEnergy() * pow(std::sin(GetElectronTheta()), 2);
    double etheta = GetElectronEnergy() * pow(std::sin(GetElectronTheta()), 2);
    //std::cout << "E_lep*sin(theta)^2: " << etheta << std::endl;
    //std::cout << " " << std::endl;
    return etheta;
  }

  double GetShowerStartZ() const {
    return GetVecElem("prong_axis_vertex", 0, 2);
  }

  double GetProtonStartZ() const {
    return GetDouble("MasterAnaDev_proton_startPointZ");
  }

  double GetProtonEndZ() const {
    return GetDouble("MasterAnaDev_proton_endPointZ");
  }

  double GetDSCalVisE() const {
    //std::cout << "prong_HCALVisE: " << GetVecElem("prong_HCALVisE", 0) << std::endl;
    //std::cout << "prong_ECALVisE: " << GetVecElem("prong_ECALVisE", 0) << std::endl;
    //if (GetVecElem("prong_HCALVisE", 0)==-999 && GetVecElem("prong_ECALVisE", 0)==-999) {	return 0; }
    if (GetVecElem("prong_HCALVisE", 0)<=0 && GetVecElem("prong_ECALVisE", 0)<=0) { return 0; }
    else { return GetVecElem("prong_HCALVisE", 0) / GetVecElem("prong_ECALVisE", 0); }
  }

  double GetODCalVisE() const {
    //std::cout << "prong_ODVisE: " << GetVecElem("prong_ODVisE", 0) << std::endl;
    //std::cout << "prong_SideECALVisE: " << GetVecElem("prong_SideECALVisE", 0) << std::endl;
    if (GetVecElem("prong_ODVisE", 0)<=0 && GetVecElem("prong_SideECALVisE", 0)<=0) { return 0; }
    else { return GetVecElem("prong_ODVisE", 0) / GetVecElem("prong_SideECALVisE", 0); }
  }

  double GetPsi() const {
    return GetDouble("Psi");
  }

  int GetHasTracks() const {
    if (GetInt("n_prongs") > 0) return 1;
    else return 0;
  }

  int GetStartPointVertexMultiplicity() const {
    return GetInt("StartPointVertexMultiplicity");
  }

  int GetHasNoVertexMismatch() const {
    return GetInt("HasNoVertexMismatch");
  }

  int GetVertexTrackMultiplicity() const {
    return GetInt("VertexTrackMultiplicity");
  }

  double GetNonMIPClusFrac() const {
    return GetVecElem("prong_NonMIPClusFrac", 0);
  }

  //For afterpulsing cut
  double GetFirstFireFraction() const {
    return GetVecElem("prong_FirstFireFraction", 0);
  }

  //todo: rename this, its not the highest just the first (electron candidate is always the first prong)
  double GetEMLikeShowerScore() const {
    double score = GetVecElem("prong_part_score", 0);
    //std::cout << "emscore: " << score << std::endl;
    return score;
  }

  double GetTransverseGapScore() const {
    return GetVecElem("prong_TransverseGapScore", 0);
  }

  double GetMeanFrontdEdx() const {
    return GetVecElem("prong_dEdXMeanFrontTracker", 0);
  }

  int GetExitsBack() const {
    return GetInt("HasNoBackExitingTracks");
  }

  //Number of isolated blob prongs
  int GetNIsoBlobs() const {
    return GetInt("nonvtx_iso_blobs_energy_in_prong_sz");
  }

  //Summed energy of all isolated blob prongs
  double GetIsoBlobsEnergy() const {
    return GetDouble("nonvtx_iso_blobs_energy");
  }

  
  ROOT::Math::XYZTVector GetVertex() const {
    ROOT::Math::XYZTVector result;
    result.SetCoordinates(GetVec<double>("vtx").data());
    return result;
  }

  //This returns 0 or 1, 1 if its in it, 0 if its not
  int GetWithinFiducialApothem() const {
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

  virtual int GetTDead() const {
    return GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj");
  }  

  //Needed for response/recoil systematics to work properly after Anezka's 2023 updates
  //Which correspond to reworked particle response branches in p4 and onwards
  //but like, this is target region specific... 
  double ApplyCaloTuning(double calRecoilE) const{
    //Path to MParamFiles calibration files
    std::string pwd = "$MPARAMFILESROOT/data/Calibrations/energy_calib/CalorimetryTunings.txt";
    
    //THESE TWO ARE SUPERRRR DIFFERENT, how do I know which to use??? do I even need to worry about this? But how do I do response systematics without it??
    util::CaloCorrection Nu_Tracker(pwd.c_str(), "NukeCC_Nu_Tracker");
    //util::CaloCorrection Nu_Tracker(pwd.c_str(), "CCNuE_Nu_Tracker");

    return Nu_Tracker.eCorrection(calRecoilE/1000.)*1000.; //MeV
  }

  double GetCalRecoilEnergy() const {
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
    // 7 = nu + e elastic (I should check this, it's like 11 events out of a full 180000 entry tuple...)
    // 8 = MEC (meson exchange current?)
    return GetInt("mc_intType");
  }
  
  virtual bool IsNuEelastic() const {
    return GetInt("mc_intType") == 7;
  }
  
  //This is kind of hacky, might change, but because I'm writing a bunch of proton truth calculators, I don't want to keep finding the 
  //Highest energy proton from mc_FSPartPDG. So this is a function which should only get called by other functions here in CVUniverse.h
  //Just returns the index value (in the final state particles array of the particular event) of the highest energy proton
  //and -999 if no true protons
  int GetHighestEnergySignalProtonIndex() const {
    double highestEnergy = -999;
    int index = -999;
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
	ROOT::Math::RotationX r(-3.3 * (pi / 180.));
	double protonTheta = (r(p)).Theta()*(180/pi); //in degrees
	//if (protonP>450 && protonP<1200 && protonTheta<70){
	if (protonP>450 && protonP<1200 && (protonTheta<70 || protonTheta>110)){ //Testing allowing backwards protons??
	  highestEnergy = energies[i];
	  index = i;
	}
      }
    }
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

  //For background categorization, how many events have protons that just don't land in my signal region cuts?
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

  ROOT::Math::XYZTVector GetTrueVertex() const {
    ROOT::Math::XYZTVector result;
    result.SetCoordinates(GetVec<double>("mc_vtx").data());
    return result;
  }

  virtual int GetCurrent() const { return GetInt("mc_current"); }

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

  //Checks for pions & kaons
  bool GetHasFSMeson() const { 
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    bool hasMeson = false;
    for (int i = 0; i < FSParticles.size(); i++){
      //std::cout << "final state particle: " << i << ": " << FSParticles[i] << std::endl;
      if (abs(FSParticles[i]) == 211 || abs(FSParticles[i]) == 321 || abs(FSParticles[i]) == 311 || abs(FSParticles[i]) == 130 || abs(FSParticles[i]) == 111){
	hasMeson = true;
      }
    }
    //std::cout << "hasMeson: " << hasMeson<< std::endl;
    return hasMeson;
  }

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
