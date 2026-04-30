//Registry/Factory of all my reweighters and how to construct them.
// This exists to remove a bunch of boilerplate in the event loop
// Basically the event loop will loop through the weighters in the config file,
// and instantiate the weighter classes using the constructors here and the arguments from the config

// Author C. Pernas
#pragma once
#include <yaml-cpp/yaml.h> //needed because my cut registry creators have yaml node in the signature, yaml doesn't actually get read here
#include "PlotUtils/Model.h"
#include "PlotUtils/FluxAndCVReweighter.h"
#include "PlotUtils/GENIEReweighter.h"
#include "PlotUtils/LowRecoil2p2hReweighter.h"
#include "PlotUtils/RPAReweighter.h"
#include "PlotUtils/MINOSEfficiencyReweighter.h"
#include "PlotUtils/LowQ2PiReweighter.h"
#include "PlotUtils/AMUDISReweighter.h"
#include "PlotUtils/FSIReweighter.h"

using ReweighterBase = PlotUtils::Reweighter<CVUniverse, CCNuEEvent>;

struct ReweighterDef {
  std::string name;
  ReweighterBase* (*creator)(const YAML::Node&);
  YAML::Node defaultNode = YAML::Node{};
};

//Defaults for GENIE & LowQ2Pi reweight options, for use with pre-defined MnvTunev1 & v2
inline const YAML::Node GENIE_defaults = [](){
  YAML::Node n; n["useNonResPi"] = true;  n["useDeutreriumPionTune"] = false;  return n; }();
inline const YAML::Node LowQ2Pi_defaults = [](){
  YAML::Node n; n["channel"] = "JOINT"; return n; }();

ReweighterBase* MakeFluxAndCV(const YAML::Node& node) {
  return new PlotUtils::FluxAndCVReweighter<CVUniverse, CCNuEEvent>();
}

ReweighterBase* MakeGENIE(const YAML::Node& node) {  
  bool useNonResPi = node["useNonResPi"].as<bool>();
  bool useDeutreriumPionTune = node["useDeutreriumPionTune"].as<bool>();
  return new PlotUtils::GENIEReweighter<CVUniverse, CCNuEEvent>(useNonResPi, useDeutreriumPionTune);
}

ReweighterBase* MakeLowRecoil2p2h(const YAML::Node& node) {
  return new PlotUtils::LowRecoil2p2hReweighter<CVUniverse, CCNuEEvent>();
}

ReweighterBase* MakeMINOSEfficiency(const YAML::Node& node) {
  return new PlotUtils::MINOSEfficiencyReweighter<CVUniverse, CCNuEEvent>();
}

ReweighterBase* MakeRPA(const YAML::Node& node) {
  return new PlotUtils::RPAReweighter<CVUniverse, CCNuEEvent>();
}

ReweighterBase* MakeLowQ2Pi(const YAML::Node& node) {  
  std::string channel = node["channel"].as<std::string>();
  return new PlotUtils::LowQ2PiReweighter<CVUniverse, CCNuEEvent>(channel);
}

ReweighterBase* MakeAMUDIS(const YAML::Node& node) {  
  return new PlotUtils::AMUDISReweighter<CVUniverse, CCNuEEvent>();
}

ReweighterBase* MakeFSI(const YAML::Node& node) {  
  bool useElastic = node["useElastic"].as<bool>();
  bool useAbsorption = node["useAbsorption"].as<bool>();
  return new PlotUtils::FSIReweighter<CVUniverse, CCNuEEvent>(useElastic, useAbsorption);
}

inline const std::vector<ReweighterDef> MnvTunev1_reweights = {
  { "FluxAndCVReweighter", &MakeFluxAndCV },
  { "GENIEReweighter", &MakeGENIE, GENIE_defaults },
  { "LowRecoil2p2hReweighter", &MakeLowRecoil2p2h },
  { "RPAReweighter", &MakeRPA}
};

inline const std::vector<ReweighterDef> MnvTunev2_reweights = {
  { "FluxAndCVReweighter", &MakeFluxAndCV },
  { "GENIEReweighter", &MakeGENIE, GENIE_defaults },
  { "LowRecoil2p2hReweighter", &MakeLowRecoil2p2h },
  { "RPAReweighter", &MakeRPA},
  { "LowQ2PiReweighter", &MakeLowQ2Pi, LowQ2Pi_defaults },
};

inline const std::vector<ReweighterDef> ALL_REWEIGHTS = {
  { "FluxAndCVReweighter", &MakeFluxAndCV },
  { "GENIEReweighter", &MakeGENIE },
  { "LowRecoil2p2hReweighter", &MakeLowRecoil2p2h },
  { "MINOSEfficiencyReweighter", &MakeMINOSEfficiency }, //should never have to use this...
  { "RPAReweighter", &MakeRPA},  
  { "LowQ2PiReweighter", &MakeLowQ2Pi},
  { "AMUDISReweighter", &MakeAMUDIS},
  { "FSIReweighter", &MakeFSI},
};
