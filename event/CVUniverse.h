// =============================================================================
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

#include <iostream>
#include <TMath.h>

#include "PlotUtils/MinervaUniverse.h"
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
  static constexpr double M_nucleon = (1.5*M_n+M_p)/2.5;

  static constexpr int PDG_n = 2112;
  static constexpr int PDG_p = 2212;

  // ========================================================================
  // Write a "Get" function for all quantities access by your analysis.
  // For composite quantities (e.g. Enu) use a calculator function.
  //
  // In order to properly calculate muon variables and systematics use the
  // various functions defined in MinervaUniverse.
  // E.g. GetPmu, GetEmu, etc.
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

  int GetInteractionType() const {
    return GetInt("mc_intType");
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

  //----------------------Start of Carlos's extra getters----------------------------------
  
  virtual bool IsNuEelastic() const {
    return GetInt("mc_intType") == 7;
  }
  
  //in GeV
  double GetElectronEnergyTrue() const {
    double out = GetVecElem("mc_primFSLepton", 3)/1000.;
    //std::cout << "True Lepton Energy: " << out << std::endl;
    return out;
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
	ROOT::Math::RotationX r(-3.3 * (3.141592 / 180.));
	double protonTheta = (r(p)).Theta()*(180/3.141592); //in degrees
	if (protonP>450 && protonP<1200 && protonTheta<70){
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


  //gets highest true proton kinetic energy in GeV, -999 (- M_p for now, will fix) if no true protons
  double GetTrueProtonKE() const {
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
  double GetTrueProtonAngle() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      ROOT::Math::XYZVector p(GetVecElem("mc_FSPartPx", i), GetVecElem("mc_FSPartPy", i), GetVecElem("mc_FSPartPz", i));
      ROOT::Math::RotationX r(-3.3 * (3.141592 / 180.));
      double out = (r(p)).Theta()*(180/3.141592);
      //std::cout << "True (highest energy) proton Theta: " << out << std::endl;
      return out;
    }
    else {
      return -999; 
    }
  }
  
  //Angle (in deg) between true primary electron and true highest energy proton in rad, -999 if no protons
  double GetTrueProtonToElectronAngle() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      std::vector<double> electronp = GetVecDouble("mc_primFSLepton");

      ROOT::Math::XYZVector electron(electronp[0], electronp[1], electronp[2]);
      ROOT::Math::XYZVector proton(GetVecElem("mc_FSPartPx", i), GetVecElem("mc_FSPartPy", i), GetVecElem("mc_FSPartPz", i));
      
      double dotProd = electron.Dot(proton);

      //no shot this just works first try
      double out = std::acos(dotProd/( sqrt( electron.Mag2() * proton.Mag2() ))) * (180/3.141592);
      //std::cout << "True angle between electron and (highest energy) proton: " << out << std::endl;
      return out;
    }
    else {
      return -999; 
    }
  }
  
  //MasterAnaDev_proton_P_fromdEdx, convert to GeV
  double GetProtonP() const {
    double protonP = GetDouble("MasterAnaDev_proton_P_fromdEdx");
    //std::cout << "True (highest) Proton Momentum: " << protonP << std::end; 
    if (protonP>-9999){
      return protonP/1000.;
    }
    else{
      return protonP;
    }
  }
  
  //MasterAnaDev_proton_theta, then convert to deg
  double GetProtonTheta() const {
    double protonTheta = GetDouble("MasterAnaDev_proton_theta");
    //std::cout << "True (highest) Proton Momentum: " << protonP << std::end; 
    if (protonTheta > -9999){
      return protonTheta * (180/3.141592);
    }
    else {
      return protonTheta;
    }
  }


  //Returns momentum of the highest energy true final state proton (in GeV)
  double GetTrueProtonP() const {
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
  double GetTrueProtonPt() const {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      //returns momentum in GeV, so convert to MeV which is what M_p is in, should be good
      double protonP = GetTrueProtonP()*1000.;
      double protonAngleRad = GetTrueProtonAngle()*(3.141592/180);
      
      double out = (protonP * std::sin(protonAngleRad))/1000.;	
      //std::cout << "True (highest energy) proton pT: " << out << std::endl;
      return out;
    }
    else {
      return -999; 
    }
  }
  
  double GetTrueElectronAngle() const {
    std::vector<double> electronp = GetVecDouble("mc_primFSLepton");
    ROOT::Math::XYZVector p(electronp[0], electronp[1], electronp[2]);   
    ROOT::Math::RotationX r(-3.3 * (3.141592 / 180.));
    return (r(p)).Theta();    
  }

  //True lepton pT (in GeV)
  double GetTrueLeptonPt() const {
    double leptonP = sqrt(pow(GetVecElem("mc_primFSLepton", 3), 2) - pow(M_e, 2));
    double out = (leptonP * std::sin(GetTrueElectronAngle()))/1000.;
    //std::cout << "True lepton pT: " << out << std::endl;
    return out;
  }

  //This is for if I want to have a truth variable, i.e. I want to plot true proton momentum vs true proton angle
  //I need to put something for the reco var, and I don't want to write the function for reco proton momentum and angle yet.
  double GetDummyRecoVar() const {
    return -999; 
  }

  //getter for the branch improved_nmichel, for the michel electron cut. This will probably get fancier. 
  int GetImprovedNMichel() const {
    return GetInt("improved_nmichel");
  }
  
  //Hadronic available energy (no neutrons) in GeV
  double GetEavail() const {
    double E = GetDouble("blob_recoil_E_tracker") + GetDouble("blob_recoil_E_ecal");
    double Eavail = (E * 1.17 - ((0.008 * GetVecElem("prong_part_E", 0, 3)) + 5))/1000.;
    //std::cout << "Reco Eavail: " << Eavail << std::endl;
    return Eavail;
  }

  //This is the same e-avail as above, but with all tracked proton calibrated energy subtracted off (
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

  
  double GetEnu() const {
    double Enu = GetEavail() + GetElectronEnergy();
    //std::cout << "Reco Enu: " << Enu << std::endl;
    return Enu;
  }

  //Fairly confident in this one?
  double GetElectronEnergy() const {
    return GetVecElem("prong_part_E", 0, 3)/1000.;
  }

  //Copied the overall format here from Ryan/Hang/Sarah's CCNue kinematics calculator, there may be a better way to do this?
  //Also extremely suspect for right now
  //returns in rad
  double GetElectronAngle() const { 
    std::vector<std::vector<double> > electronp = GetVecOfVecDouble("prong_part_E");
    ROOT::Math::XYZVector p(electronp[0][0], electronp[0][1], electronp[0][2]);
    ROOT::Math::RotationX r(-3.3 * (3.141592 / 180.));
    //r = ROOT::Math::RotationX(-3.3 * (3.141592 / 180.))(r);
    return (r(p)).Theta();
  }

  double GetProtonESCNodeChi2() const {
    //I really have absolutely no idea where these hardcoded means and sigmas come from, these were in Dan's CCQENu cuts...
    //std::vector<double> means;
    //means.push_back(31.302);
    //means.push_back(11.418);
    //means.push_back(9.769);
    //means.push_back(8.675);
    //means.push_back(7.949);
    //std::vector<double> sigmas;
    //sigmas.push_back(8.997);
    //sigmas.push_back(3.075);
    //sigmas.push_back(2.554);
    //sigmas.push_back(2.484);
    //sigmas.push_back(2.232);
    double z = GetProtonEndZ();
    if (z < -1){ return 999; } //if there's no reco proton, this cut fails
    
    double means[5] = {31.302, 11.418, 9.769, 8.675, 7.949};
    double sigmas[5] = {8.997, 3.075, 2.554, 2.484, 2.232};
    
    double nodeEnergyVal = 0;
    double chi2=0;
    int n_nodes = GetInt("MasterAnaDev_proton_nodes_nodesNormE_sz");
    std::vector<double> nodesNormE = GetVecDouble("MasterAnaDev_proton_nodes_nodesNormE");
    if(n_nodes>5){
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
    //Bit hacky cause in Dan's original CCQENu this was a separate function, but I don't think I'll ever need it otherwise...
    else{
      bool pass;

      //Cut Values based on 22302 (carlos - i assume this is a docdb number? -> yes)
      //Node 0-1
      //Node 2
      //Node 3
      //Node 4
      //Node 5
      //Node 6
      
      double cutval1 = 19;
      double cutval2 = 10;
      double cutval3 = 9;
      double cutval4 = 8;
      double cutval5 = 5;
      
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
      pass = true;
      if(pass) chi2=0;
      else chi2=75;
    }
    //std::cout << "chi2 = " << chi2 << std::endl;
    return chi2;
  }

  
  double GetElectronAngleDeg() const {
    double rad = GetElectronAngle();
    double deg = rad*180/3.141592;
    return deg;
  }
  //p sure this is wrong because it will include invisible energy (from neutrons n shit)
  //Ok for now though, I think that only comes into play when we start looking at migration matrices etc
  //Focus on reco values right now I think so I can match Ryan's selection
  double GetEavailTrue() const {
    return (GetDouble("mc_incomingE") - GetVecElem("mc_primFSLepton", 3))/1000.;
  }

  //Apparently this is already defined somewhere, right now I'm assuming it works like this but should double check that at some point
  //double GetEnuTrue() const {
  //return GetDouble("mc_incomingE")/1000.;
  //}

  //LeptonP in GeV
  double GetLeptonP() const{
    double M_lep_sqr = pow(M_e, 2) / pow(10, 6);  //over 10^6 to convert to GeV^2
    double leptonP = pow(GetElectronEnergy(), 2)-M_lep_sqr;
    if (leptonP > 0){
      return sqrt(leptonP);
    }
    else return 0; 
    //return sqrt(std::max(0, (pow(GetElectronEnergy(), 2) - M_lep_sqr)));
  }

  double GetLeptonPt() const{
    return (GetLeptonP() * std::sin(GetElectronAngle()));
  }

  //This is E_lep * theta_lep^2
  double GetETheta() const {
    //double etheta = GetElectronEnergy() * pow(std::sin(GetElectronAngle()), 2);
    double etheta = GetElectronEnergy() * pow(GetElectronAngle(), 2);
    //std::cout << "E*theta^2: " << etheta << std::endl;
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
    if (GetVecElem("prong_HCALVisE", 0)==0 && GetVecElem("prong_ECALVisE", 0)==0) { return 0; }
    else { return GetVecElem("prong_HCALVisE", 0) / GetVecElem("prong_ECALVisE", 0); }

  }

  double GetODCalVisE() const {
    //std::cout << "prong_ODVisE: " << GetVecElem("prong_ODVisE", 0) << std::endl;
    //std::cout << "prong_SideECALVisE: " << GetVecElem("prong_SideECALVisE", 0) << std::endl;
    if (GetVecElem("prong_ODVisE", 0)==0 && GetVecElem("prong_SideECALVisE", 0)==0) { return 0; }
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

  //For afterpulsing
  double GetFirstFireFraction() const {
    return GetVecElem("prong_FirstFireFraction", 0);
  }

  //todo: rename this
  double GetHighestEMLikeShowerScore() const {
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
  
  ROOT::Math::XYZTVector GetVertex() const
  {
    ROOT::Math::XYZTVector result;
    result.SetCoordinates(GetVec<double>("vtx").data());
    return result;
  }

  ROOT::Math::XYZTVector GetTrueVertex() const
  {
    ROOT::Math::XYZTVector result;
    result.SetCoordinates(GetVec<double>("mc_vtx").data());
    return result;
  }

  virtual int GetTDead() const {
    return GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj");
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
