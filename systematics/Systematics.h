#ifndef Systematics_h
#define Systematics_h

//==============================================================================
// Get Several standard MINERvA systematics
//==============================================================================

//Includes from this package (NuE_TKI)
#include "event/CVUniverse.h"
#include "ElectronSystematics.h"

//PlotUtils includes
#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/MinosEfficiencySystematics.h"
#include "PlotUtils/MuonSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "PlotUtils/MuonResolutionSystematics.h"
#include "PlotUtils/AngleSystematics.h"
#include "PlotUtils/ResponseSystematics.h"
#include "PlotUtils/MnvHadronReweight.h"
#include "PlotUtils/GeantHadronSystematics.h"
#include "PlotUtils/TargetMassSystematics.h"

typedef std::map<std::string, std::vector<CVUniverse*>> UniverseMap;
UniverseMap GetTestSystematics(PlotUtils::ChainWrapper* chain)
{
  // return map
  UniverseMap error_bands;

  // CV
  error_bands[std::string("cv")].push_back(new CVUniverse(chain));

  //this is here to help me figure out why this affects my analysis...
  //UniverseMap bands_minoseff = PlotUtils::GetMinosEfficiencySystematicsMap<CVUniverse>(chain);
  //error_bands.insert(bands_minoseff.begin(), bands_minoseff.end());
  
  const bool use_ID = true;
  const bool use_OD = true;
  std::string name_tag = "allNonMuonClusters";
  const bool use_neutron = false;
  const bool use_new = false;
  const bool use_proton = true;
  //UniverseMap bands_response = PlotUtils::GetResponseSystematicsMap<CVUniverse>(chain, use_ID, use_OD, name_tag, use_neutron, use_new, use_proton);
  //error_bands.insert(bands_response.begin(), bands_response.end());

  //UniverseMap bands_mass = PlotUtils::GetTargetMassSystematicsMap<CVUniverse>(chain);
  //error_bands.insert(bands_mass.begin(), bands_mass.end());

  //UniverseMap bands_electron = PlotUtils::GetElectronEnergyShiftSystematicsMap<CVUniverse>(chain);
  //error_bands.insert(bands_electron.begin(), bands_electron.end());
  
  //UniverseMap bands_leakage = PlotUtils::GetLeakageSystematicsMap<CVUniverse>(chain);
  //error_bands.insert(bands_leakage.begin(), bands_leakage.end());  
  
  //UniverseMap bands_electron_angle = PlotUtils::GetElectronAngleShiftSystematicsMap<CVUniverse>(chain);
  //error_bands.insert(bands_electron_angle.begin(), bands_electron_angle.end());  

  UniverseMap bands_lowq2 = PlotUtils::GetLowQ2PiSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_lowq2.begin(), bands_lowq2.end());
  
  return error_bands;
}

UniverseMap GetNuETKISystematics(PlotUtils::ChainWrapper* chain)
{
  // return map
  UniverseMap error_bands;

  // CV
  error_bands[std::string("cv")].push_back(new CVUniverse(chain));

  //========================================================================
  // FLUX
  //========================================================================
  UniverseMap bands_flux = PlotUtils::GetFluxSystematicsMap<CVUniverse>(chain, CVUniverse::GetNFluxUniverses());
  error_bands.insert(bands_flux.begin(), bands_flux.end());

  //========================================================================
  // GENIE
  //========================================================================
  // Standard
  UniverseMap bands_genie = PlotUtils::GetGenieSystematicsMap<CVUniverse>(chain); //PlotUtils::GetStandardGenieSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_genie.begin(), bands_genie.end());

  //Ok sweet this works, although I don't actually need this one specifically 
  //error_bands["NormCCRes"].push_back(new PlotUtils::GenieNormCCResUniverse<CVUniverse>(chain, -1.));
  //error_bands["NormCCRes"].push_back(new PlotUtils::GenieNormCCResUniverse<CVUniverse>(chain, 1.));
  
  //========================================================================
  // MnvTunes
  //========================================================================
  // RPA
  UniverseMap bands_rpa = PlotUtils::GetRPASystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_rpa.begin(), bands_rpa.end());

  // 2P2H
  UniverseMap bands_2p2h = PlotUtils::Get2p2hSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_2p2h.begin(), bands_2p2h.end());

  // Low Q2 pion suppression (i think this is a tune that I don't have...
  UniverseMap bands_lowq2 = PlotUtils::GetLowQ2PiSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_lowq2.begin(), bands_lowq2.end());
  
  //========================================================================
  // Muons
  //========================================================================
  // Muon reco in MINERvA -- Catchall systematic for pmu reco in minerva.
  // Lateral-only. Shifts pmu.
  //UniverseMap bands_muon_minerva = PlotUtils::GetMinervaMuonSystematicsMap<CVUniverse>(chain);
  //error_bands.insert(bands_muon_minerva.begin(), bands_muon_minerva.end());

  //UniverseMap bands_electron = PlotUtils::GetMinervaMuonSystematicsMap<CVUniverse>(chain);
  //error_bands.insert(bands_muon_minerva.begin(), bands_muon_minerva.end());

  // Muons in MINOS -- Catchall systematic for wiggle solution -- correlates
  // flux universes and minos muon momentum reco.
  // Lateral AND Vertical systematic. Shifts Pmu and GetFluxAndCVUniverse.
  // Expect a non-zero systematic even when no pmu involved.
  //
  //Carlos - Can I remove these??? Apparently it has flux implications, and the universes do change a bit
  //UniverseMap bands_muon_minos = PlotUtils::GetMinosMuonSystematicsMap<CVUniverse>(chain);
  //error_bands.insert(bands_muon_minos.begin(), bands_muon_minos.end());
  //UniverseMap bands_minoseff = PlotUtils::GetMinosEfficiencySystematicsMap<CVUniverse>(chain);
  //error_bands.insert(bands_minoseff.begin(), bands_minoseff.end());

  //UniverseMap bands_muon_resolution = PlotUtils::GetMuonResolutionSystematicsMap<CVUniverse>(chain);
  //error_bands.insert(bands_muon_resolution.begin(), bands_muon_resolution.end());
  
  UniverseMap bands_geant = PlotUtils::GetGeantHadronSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_geant.begin(), bands_geant.end());

  // Beam angle: these seem to not do anything, so I'll have to write up whatever Hang & Sarah have
  // it's like shower angle something or other
  UniverseMap bands_angle = PlotUtils::GetAngleSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_angle.begin(), bands_angle.end());

  //========================================================================
  // Particle Response Systematics
  //========================================================================
  //Using Anezka's reworked recoil systematics after p4 tuples
  // these are working, but as it stands I don't actually use these exact quantities anywhere in my analysis. See the spreadsheet for details
  const bool use_ID = true;
  const bool use_OD = true;
  std::string name_tag = "allNonMuonClusters";
  const bool use_neutron = false;
  const bool use_new = false;
  const bool use_proton = true;
  UniverseMap bands_response = PlotUtils::GetResponseSystematicsMap<CVUniverse>(chain, use_ID, use_OD, name_tag, use_neutron, use_new, use_proton);
  error_bands.insert(bands_response.begin(), bands_response.end());

  UniverseMap bands_mass = PlotUtils::GetTargetMassSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_mass.begin(), bands_mass.end());

  UniverseMap bands_electron = PlotUtils::GetElectronEnergyShiftSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_electron.begin(), bands_electron.end());
  
  UniverseMap bands_leakage = PlotUtils::GetLeakageSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_leakage.begin(), bands_leakage.end());  
  
  UniverseMap bands_electron_angle = PlotUtils::GetElectronAngleShiftSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_electron_angle.begin(), bands_electron_angle.end());  
  
  return error_bands;
}
//Legacy, this is what was in the tutorial
UniverseMap GetStandardSystematics(PlotUtils::ChainWrapper* chain)
{
  // return map
  UniverseMap error_bands;

  // CV
  error_bands[std::string("cv")].push_back(new CVUniverse(chain));

  //========================================================================
  // FLUX
  //========================================================================
  UniverseMap bands_flux =
      PlotUtils::GetFluxSystematicsMap<CVUniverse>(chain, CVUniverse::GetNFluxUniverses());
  error_bands.insert(bands_flux.begin(), bands_flux.end());

  //========================================================================
  // GENIE
  //========================================================================
  // Standard
  UniverseMap bands_genie =
      PlotUtils::GetGenieSystematicsMap<CVUniverse>(chain); //PlotUtils::GetStandardGenieSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_genie.begin(), bands_genie.end());

  //========================================================================
  // MnvTunes
  //========================================================================
  // RPA
  UniverseMap bands_rpa = PlotUtils::GetRPASystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_rpa.begin(), bands_rpa.end());

  // 2P2H
  UniverseMap bands_2p2h = PlotUtils::Get2p2hSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_2p2h.begin(), bands_2p2h.end());

  //========================================================================
  // Muons
  //========================================================================
  // Muon reco in MINERvA -- Catchall systematic for pmu reco in minerva.
  // Lateral-only. Shifts pmu.
  UniverseMap bands_muon_minerva =
      PlotUtils::GetMinervaMuonSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_muon_minerva.begin(), bands_muon_minerva.end());

  // Muons in MINOS -- Catchall systematic for wiggle solution -- correlates
  // flux universes and minos muon momentum reco.
  // Lateral AND Vertical systematic. Shifts Pmu and GetFluxAndCVUniverse.
  //
  // Expect a non-zero systematic even when no pmu involved.
  UniverseMap bands_muon_minos =
     PlotUtils::GetMinosMuonSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_muon_minos.begin(), bands_muon_minos.end());

  // Vertical only
  UniverseMap bands_minoseff =
      PlotUtils::GetMinosEfficiencySystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_minoseff.begin(), bands_minoseff.end());

  UniverseMap bands_muon_resolution = PlotUtils::GetMuonResolutionSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_muon_resolution.begin(), bands_muon_resolution.end());

  UniverseMap bands_geant = PlotUtils::GetGeantHadronSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_geant.begin(), bands_geant.end());

  // Beam angle
  UniverseMap bands_angle = PlotUtils::GetAngleSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_angle.begin(), bands_angle.end());

  // Hadron inelastics cross sections
  //TODO: There's some special recoil function I need to write for the response systematics to work correctly
  /*UniverseMap bands_response = PlotUtils::GetResponseSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_response.begin(), bands_response.end());*/

  return error_bands;
}

#endif  // Systematics_h
