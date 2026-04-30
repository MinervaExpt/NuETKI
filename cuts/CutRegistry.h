//Registry/Factory of all my cuts and how to construct them.
// This exists to remove a bunch of boilerplate in the event loop
// Basically the event loop will loop through the cuts in the config file,
// and create the cut objects using the constructors here and the arguments from the config

// Author C. Pernas
#pragma once
#include <yaml-cpp/yaml.h> //needed because my cut registry creators have yaml node in the signature, yaml doesn't actually get read here

using CutBase = PlotUtils::Cut<CVUniverse, CCNuEEvent>;

struct CutDef {
  std::string name;
  CutBase* (*creator)(const YAML::Node&);
};

CutBase* MakeHasTracks(const YAML::Node&) {
  return new reco::HasTracks<CVUniverse, CCNuEEvent>();
}

CutBase* MakeNoVertexMismatch(const YAML::Node&) {
  return new reco::NoVertexMismatch<CVUniverse, CCNuEEvent>();
}

CutBase* MakeZRange(const YAML::Node& node) {
  double minZ = node["minZ"].as<double>();
  double maxZ = node["maxZ"].as<double>();
  return new reco::ZRange<CVUniverse, CCNuEEvent>("Vertex in Tracker", minZ, maxZ);
}

CutBase* MakeApothem(const YAML::Node& node) {
  double apothem = node["apothem"].as<double>();
  return new reco::Apothem<CVUniverse, CCNuEEvent>(apothem);
}

CutBase* MakeStartPointVertexMultiplicity(const YAML::Node&) {
  return new reco::StartPointVertexMultiplicity<CVUniverse, CCNuEEvent>();
}

CutBase* MakeAfterpulsing(const YAML::Node& node) {
  double minFirstFireFraction = node["minFirstFireFraction"].as<double>();
  return new reco::Afterpulsing<CVUniverse, CCNuEEvent>(minFirstFireFraction);
}

CutBase* MakeNoDeadtime(const YAML::Node&) {
  return new reco::NoDeadtime<CVUniverse, CCNuEEvent>(1, "Deadtime");
}

CutBase* MakeHasNoBackExitingTracks(const YAML::Node&) {
  return new reco::HasNoBackExitingTracks<CVUniverse, CCNuEEvent>();
}

CutBase* MakeSideExitingMuon(const YAML::Node&) {
  return new reco::SideExitingMuon<CVUniverse, CCNuEEvent>();
}

CutBase* MakeDSCalVisE(const YAML::Node& node) {
  double maxDSCalRatio = node["maxDSCalRatio"].as<double>();
  return new reco::DSCalVisE<CVUniverse, CCNuEEvent>(maxDSCalRatio);
}

CutBase* MakeODCalVisE(const YAML::Node& node) {
  double maxSideCalRatio = node["maxSideCalRatio"].as<double>();
  return new reco::ODCalVisE<CVUniverse, CCNuEEvent>(maxSideCalRatio);
}

CutBase* MakeVertexTrackMultiplicity(const YAML::Node& node) {
  int min = node["minVertexTrackMultiplicity"].as<int>();
  int max = node["maxVertexTrackMultiplicity"].as<int>();
  return new reco::VertexTrackMultiplicity<CVUniverse, CCNuEEvent>(min, max);
}

CutBase* MakeTransverseGapScore(const YAML::Node& node) {
  double minScore = node["minTransverseGapScore"].as<double>();
  return new reco::TransverseGapScore<CVUniverse, CCNuEEvent>(minScore);
}

CutBase* MakeNonMIPClusterFraction(const YAML::Node& node) {
  double minFrac = node["minNonMIPClusFrac"].as<double>();
  return new reco::NonMIPClusterFraction<CVUniverse, CCNuEEvent>(minFrac);
}

CutBase* MakeEMLikeTrackScore(const YAML::Node& node) {
  double minScore = node["minEMTrackScore"].as<double>();
  return new reco::EMLikeTrackScore<CVUniverse, CCNuEEvent>(minScore);
}

CutBase* MakeMichelCut(const YAML::Node&) {
  return new reco::MichelCut<CVUniverse, CCNuEEvent>();
}

CutBase* MakeMeanFrontdEdX(const YAML::Node& node) {
  double minDEDX = node["maxMeanFrontDEDX"].as<double>();
  return new reco::MeanFrontdEdX<CVUniverse, CCNuEEvent>(minDEDX);
}

CutBase* MakeElectronEnergy(const YAML::Node& node) {
  double minEnergy = node["minElectronEnergy"].as<double>();
  return new reco::ElectronEnergy<CVUniverse, CCNuEEvent>(minEnergy);
}

CutBase* MakeModifiedEavailable(const YAML::Node& node) {
  double maxE = node["maxModifiedEAvail"].as<double>();
  return new reco::ModifiedEavailable<CVUniverse, CCNuEEvent>(maxE);
}

CutBase* MakeProtonInEvent(const YAML::Node&) {
  return new reco::ProtonInEvent<CVUniverse, CCNuEEvent>();
}

CutBase* MakeNIsoBlobs(const YAML::Node& node) {
  int maxBlobs = node["maxNIsoBlobs"].as<int>();
  return new reco::NIsoBlobs<CVUniverse, CCNuEEvent>(maxBlobs);
}

CutBase* MakeIsoBlobEnergy(const YAML::Node& node) {
  int maxBlobEnergy = node["maxIsoBlobEnergy"].as<int>();
  return new reco::IsoBlobEnergy<CVUniverse, CCNuEEvent>(maxBlobEnergy);
}

CutBase* MakeUpstreamIsoBlob(const YAML::Node& node) {
  int maxUpstream = node["furthestUpstreamIsoBlobStartZ"].as<int>();
  return new reco::UpstreamIsoBlob<CVUniverse, CCNuEEvent>(maxUpstream);
}

CutBase* MakeESC(const YAML::Node& node) {
  double maxChi2 = node["maxESCChi2"].as<double>();
  return new reco::ESC<CVUniverse, CCNuEEvent>(maxChi2);
}

CutBase* MakeEavailable(const YAML::Node&) {
  return new reco::Eavailable<CVUniverse, CCNuEEvent>();
}

CutBase* MakeProtonMomentum(const YAML::Node&) {
  return new reco::ProtonMomentum<CVUniverse, CCNuEEvent>();
}

CutBase* MakePsi(const YAML::Node&) {
  return new reco::Psi<CVUniverse, CCNuEEvent>();
}

CutBase* MakeProtonTheta(const YAML::Node&) {
  return new reco::ProtonTheta<CVUniverse, CCNuEEvent>();
}

CutBase* MakeElectronPt(const YAML::Node&) {
  return new reco::ElectronPt<CVUniverse, CCNuEEvent>();
}

CutBase* MakeEleptonSin2Theta(const YAML::Node&) {
  return new reco::EleptonSin2Theta<CVUniverse, CCNuEEvent>();
}

inline const std::vector<CutDef> ALL_CUTS = {
  //precuts
  { "HasTracks", &MakeHasTracks },
  { "HasNoBackExitingTracks", &MakeHasNoBackExitingTracks },
  { "SideExitingMuon", &MakeSideExitingMuon },
  { "NoVertexMismatch", &MakeNoVertexMismatch },
  { "ZRange", &MakeZRange },
  { "Apothem", &MakeApothem },
  { "StartPointVertexMultiplicity", &MakeStartPointVertexMultiplicity },
  { "Afterpulsing", &MakeAfterpulsing },
  { "NoDeadtime", &MakeNoDeadtime },
  { "DSCalVisE", &MakeDSCalVisE },
  { "ODCalVisE", &MakeODCalVisE },
  { "VertexTrackMultiplicity", &MakeVertexTrackMultiplicity },
  { "TransverseGapScore", &MakeTransverseGapScore },
  { "NonMIPClusterFraction", &MakeNonMIPClusterFraction },
  { "EMLikeTrackScore", &MakeEMLikeTrackScore },
  { "ElectronEnergy", &MakeElectronEnergy },
  { "ModifiedEavailable", &MakeModifiedEavailable },

  //Sideband cuts
  { "MichelCut", &MakeMichelCut },
  { "MeanFrontdEdX", &MakeMeanFrontdEdX },  

  //QELike Cuts
  { "ProtonInEvent", &MakeProtonInEvent },
  { "NIsoBlobs", &MakeNIsoBlobs },
  { "IsoBlobEnergy", &MakeIsoBlobEnergy },
  { "UpstreamIsoBlob", &MakeUpstreamIsoBlob },
  { "ESC", &MakeESC },

  //Extra cuts, typically unused in my main analysis
  { "Eavailable", &MakeEavailable },
  { "ProtonMomentum", &MakeProtonMomentum },
  { "Psi", &MakePsi },
  { "ProtonTheta", &MakeProtonTheta },
  { "ElectronPt", &MakeElectronPt },
  { "EleptonSin2Theta", &MakeEleptonSin2Theta }
};
