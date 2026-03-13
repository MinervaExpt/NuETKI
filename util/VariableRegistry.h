#pragma once
#include <vector>
#include <yaml-cpp/yaml.h>
#include "util/Binning.h"
#include "event/CVUniverse.h"


struct VariableDef {
  std::string name;
  std::string title;
  const std::vector<double>* bins;
  double (CVUniverse::*reco)() const;
  double (CVUniverse::*truth)() const;
};

//binning vectors are defined in util/Binning.h
inline const std::vector<VariableDef> ALL_VARIABLES = {
  //Reco variables
  { "E_lep", "E_{e} [GeV]", &electronEnergyBins, &CVUniverse::GetElectronEnergyGeV, &CVUniverse::GetElectronEnergyTrueGeV },
  { "E_avail", "E_{avail} [GeV]", &EavailBins, &CVUniverse::GetEavailGeV, &CVUniverse::GetEavailTrueGeV },
  { "E_nu", "E_{#nu} [GeV]", &electronEnergyBins, &CVUniverse::GetEnuGeV, &CVUniverse::GetEnuTrueGeV },
  { "Lepton_Pt", "p_{T,e} [GeV/c]", &leptonPt_bins, &CVUniverse::GetElectronPtGeV, &CVUniverse::GetElectronPtTrueGeV },
  { "Lepton_Pl", "p_{||,e} [GeV/c]", &Pl_bins, &CVUniverse::GetElectronPParallelGeV, &CVUniverse::GetElectronPParallelTrueGeV},
  { "Theta_lep", "#theta_{e} [deg]", &electronAngleBins, &CVUniverse::GetElectronThetaDeg, &CVUniverse::GetElectronThetaDegTrue},    
  { "Proton_p", "p_{p} [GeV/c]", &protonMomentumBins, &CVUniverse::GetProtonPGeV, &CVUniverse::GetProtonPTrueGeV},
  { "Proton_Pt", "p_{T,p} [GeV/c]", &protonPtBins, &CVUniverse::GetProtonPtGeV, &CVUniverse::GetProtonPtTrueGeV},
  //{ "Proton_Pl", "p_{||,p} [GeV/c]", &Pl_bins, &CVUniverse::GetProtonPParallelGeV, &CVUniverse::GetProtonPParallelTrueGeV}, //this one's pretty much never used
  { "Theta_p", "#theta_{p} [deg]", &protonAngleBins, &CVUniverse::GetProtonThetaDeg, &CVUniverse::GetProtonThetaDegTrue},    
  { "Proton_T", "#T_{p} [GeV]", &T_p_bins, &CVUniverse::GetProtonTGeV, &CVUniverse::GetProtonTTrueGeV}, 

  //tki vars
  { "DeltaPt", "#deltaP_{T} [GeV/c]", &deltaPt_bins, &CVUniverse::GetDeltaPtGeV, &CVUniverse::GetDeltaPtTrueGeV}, //momentum imbalance transverse to neutrino direction
  { "DeltaPtX", "#deltaP_{T,x} [GeV/c]", &deltaPtXBins, &CVUniverse::GetDeltaPtXGeV, &CVUniverse::GetDeltaPtXTrueGeV},
  { "DeltaPtY", "#deltaP_{T,y} [GeV/c]", &deltaPtYBins, &CVUniverse::GetDeltaPtYGeV, &CVUniverse::GetDeltaPtYTrueGeV}, 
  { "DeltaPl", "#deltaP_{l} [GeV/c]", &deltaPlBins, &CVUniverse::GetDeltaPlGeV, &CVUniverse::GetDeltaPlTrueGeV}, //momentum imbalance along neutrino direction
  { "P_n", "#P_{n} [GeV/c]", &PnBins, &CVUniverse::GetPnGeV, &CVUniverse::GetPnTrueGeV}, //struck nucleon momentum
  { "AlphaPt", "#delta#alpha_{T} [deg]", &alphaAngleBins, &CVUniverse::GetAlphaT, &CVUniverse::GetAlphaTTrue}, //boosting angle
  { "PhiPt", "#delta#phi_{T} [deg]", &phiAngleBins, &CVUniverse::GetPhiT, &CVUniverse::GetPhiTTrue}, //acoplanarity

  //true variables (it's for testing, but I'm not really sure if I need this?)
  { "E_lep_true", "True Electron Energy", &electronEnergyBins, &CVUniverse::GetElectronEnergyTrue, &CVUniverse::GetElectronEnergyTrue},
  { "ProtonKE_true", "T_{p} true [GeV}", &KE_bins, &CVUniverse::GetDummyVar, &CVUniverse::GetProtonTTrue},
  { "ProtonP_true", "p_{p} true [GeV/c}", &protonMomentumBins, &CVUniverse::GetDummyVar, &CVUniverse::GetProtonPTrue},
  { "Theta_p_true", "#theta_{p} true [deg}", &protonAngleBins, &CVUniverse::GetDummyVar, &CVUniverse::GetProtonThetaTrue},
  { "Pt_p_true", "p_{T,p} true [GeV/c}", &Pt_bins, &CVUniverse::GetDummyVar, &CVUniverse::GetProtonPtTrue},
  { "Theta_p_e_true", "#theta_{e,p} true [deg}", &protonAngleBins, &CVUniverse::GetDummyVar, &CVUniverse::GetOpeningAngleTrue},
  { "Pt_lep_true", "p_{T,e} true [GeV}", &Pt_bins, &CVUniverse::GetDummyVar, &CVUniverse::GetElectronPtTrue},

  //Cut variables
  //Some of these unfortunately returned ints, so I've changed them to doubles for the sake of making this struct work
  //That works fine for everything except the deadtime cut, for some reason. But I shouldn't really ever need to plot that
  //So I think it's fine?
  { "NoVertexMismatch", "NoVertexMismatch", &binaryCut_bins, &CVUniverse::GetHasNoVertexMismatch, &CVUniverse::GetDummyVar},
  { "VertexZ", "VertexZ [mm]", &VertexZ_bins, &CVUniverse::GetVertexZ, &CVUniverse::GetDummyVar},
  { "InApothem", "InApothem", &binaryCut_bins, &CVUniverse::GetWithinFiducialApothem, &CVUniverse::GetDummyVar},
  { "StartPointVertexMultiplicity", "StartPointVertexMultiplicity", &StartPointVertexMultiplicity_bins, &CVUniverse::GetStartPointVertexMultiplicity, &CVUniverse::GetDummyVar},
  { "Afterpulsing", "Electron prong First Fire fraction", &Afterpulsing_bins, &CVUniverse::GetFirstFireFraction, &CVUniverse::GetDummyVar},
  // { "Deadtime", "N dead discr pairs upstream", &Deadtime_bins, &CVUniverse::GetTDead, &CVUniverse::GetDummyVar}, //will stop code from compiling :(
  { "NoBackExitingTracks", "NoBackExitingTracks", &binaryCut_bins, &CVUniverse::GetExitsBack, &CVUniverse::GetDummyVar},
  { "DSCalVisE", "Ratio of DS HCalVisE/ECalVisE", &DSCalVisE_bins, &CVUniverse::GetDSCalVisE, &CVUniverse::GetDummyVar},
  { "ODCalVisE", "Ratio of OD HCalVisE/ECalVisE", &ODCalVisE_bins, &CVUniverse::GetODCalVisE, &CVUniverse::GetDummyVar},
  { "VertexTrackMultiplicity", "Vertex Track Multiplicity", &VertexTrackMultiplicity_bins, &CVUniverse::GetVertexTrackMultiplicity, &CVUniverse::GetDummyVar},
  { "TransverseGapScore", "Transverse Gap Score", &TransverseGapScore_bins, &CVUniverse::GetTransverseGapScore, &CVUniverse::GetDummyVar},
  { "NonMIPClusFrac", "Non MIP Cluster Fraction", &NonMIPClusFrac_bins, &CVUniverse::GetNonMIPClusFrac, &CVUniverse::GetDummyVar},
  { "EMLikeTrackScore", "EMShower Score", &EMScore_bins, &CVUniverse::GetEMLikeShowerScore, &CVUniverse::GetDummyVar},
  { "MichelCut", "N michels", &NMichel_bins, &CVUniverse::GetImprovedNMichel, &CVUniverse::GetDummyVar},
  { "MeanFrontdEdX", "Mean Front dE/dX [MeV/cm]", &MeanFrontDEDX_bins, &CVUniverse::GetMeanFrontdEdx, &CVUniverse::GetDummyVar},
  // { "E_lep", "E_lep [GeV]", &electronEnergyBins, &CVUniverse::GetElectronEnergy, &CVUniverse::GetDummyVar}, //already in reco section
  { "ModifiedEavailable", "E_avail - sum(proton_E) [GeV]", &Modified_E_avail_bins, &CVUniverse::GetModifiedEavail, &CVUniverse::GetDummyVar},
  { "NIsoBlobs", "# of isolated blobs", &NMichel_bins, &CVUniverse::GetNIsoBlobs, &CVUniverse::GetDummyVar},
  { "ESC", "ESC Proton Node Chi2", &ESCChi2_bins, &CVUniverse::GetProtonESCNodeChi2, &CVUniverse::GetDummyVar},

  //Potential cuts, but not ones I normally do
  { "Psi", "Psi", &Psi_bins, &CVUniverse::GetPsi, &CVUniverse::GetDummyVar},
  { "ProtonMomentum", "Primary Proton momentum [GeV}", &protonMomentumBins, &CVUniverse::GetProtonP, &CVUniverse::GetDummyVar},
  { "ProtonTheta", "{theta}_proton [deg}", &protonAngleBins, &CVUniverse::GetProtonThetaDeg, &CVUniverse::GetDummyVar},
  //{ "LeptonPt", "LeptonPt [GeV}", &Pt_bins, &CVUniverse::GetElectronPt, &CVUniverse::GetDummyVar}, //also already in reco section
  { "Etheta", "E_lep * theta_lep^2", &Etheta_bins, &CVUniverse::GetETheta, &CVUniverse::GetDummyVar }
};
