#ifndef CCNuEEvent_h
#define CCNuEEvent_h

#include "event/CVUniverse.h"
#include "TFile.h"
#include "TTree.h"
#include <bitset>

//This struct definition doesn't have to be HERE specifically, just in some .h file that I import into runEventLoop and my sideband Study objects...
struct OutputTreeManager {
  bool enableOutputTree = false;  
  TFile* outputTreeFile = nullptr;
  TTree* outputTree = nullptr;
  
  //This is so dumb, but basically maps each sideband name to a8 set of entry numbers (which corresponds to i in my main runEventLoop)
  //Avoids the triple event filling thing into my tree but jesus this struct is getting out of hand...
  std::map<std::string, std::unordered_set<Long64_t>> filledEntriesPerTree; 

  const std::vector<std::string> branches_to_keep = {
    // prong variables
    "prong_part_E",
    "prong_ECALCalibE",
    "prong_axis_vertex",
    "prong_HCALVisE",
    "prong_ECALVisE",
    "prong_ODVisE",
    "prong_SideECALVisE",
    "prong_NonMIPClusFrac",
    "prong_FirstFireFraction",
    "prong_part_score",
    "prong_TransverseGapScore",
    "prong_dEdXMeanFrontTracker",
    
    // MasterAnaDev proton variables
    "MasterAnaDev_proton_E_fromdEdx",
    "MasterAnaDev_proton_T_fromdEdx",
    "MasterAnaDev_proton_P_fromdEdx",
    "MasterAnaDev_proton_Px_fromdEdx",
    "MasterAnaDev_proton_Py_fromdEdx",
    "MasterAnaDev_proton_Pz_fromdEdx",
    "MasterAnaDev_proton_theta",
    "MasterAnaDev_proton_phi",
    "MasterAnaDev_proton_calib_energy",
    "MasterAnaDev_proton_startPointZ",
    "MasterAnaDev_proton_endPointZ",
    "MasterAnaDev_proton_nodes_nodesNormE",
    "MasterAnaDev_proton_nodes_nodesNormE_sz",
    "MasterAnaDev_sec_protons_T_fromCalo",
    "MasterAnaDev_hadron_isSideECAL",
    "MasterAnaDev_hadron_isODMatch",
    
    // blob/recoil variables
    "blob_recoil_E_tracker",
    "blob_recoil_E_ecal",
    "blob_recoil_E_hcal",
    "blob_recoil_E_od",
    "blob_recoil_E_nucl",
    "nonvtx_iso_blobs_start_position_z_in_prong",
    "nonvtx_iso_blobs_start_position_z_in_prong_sz",
    "nonvtx_iso_blobs_energy_in_prong_sz",
    "nonvtx_iso_blobs_energy",
    
    // MC truth variables
    "mc_run",
    "mc_subrun",
    "mc_nthEvtInFile",
    "mc_targetZ",
    "mc_targetA",
    "mc_targetNucleon",
    "mc_primFSLepton",
    "mc_initNucVec",
    "mc_Q2",
    "mc_w",
    "mc_FSPartE",
    "mc_FSPartPx",
    "mc_FSPartPy",
    "mc_FSPartPz",
    "mc_FSPartPDG",
    "mc_nFSPart",
    "mc_intType",
    "mc_current",
    "mc_incoming",
    "mc_incomingE",
    "mc_targetNucleon",
    "mc_targetZ",
    "mc_vtx",
    "mc_Bjorkenx",
    "mc_Bjorkeny",
    "mc_fr_nuParentID",
    "mc_fr_nuAncestorID",
    
    // reco event variables
    "vtx",
    "recoil_summed_energy",
    "improved_nmichel",
    "n_prongs",
    "StartPointVertexMultiplicity",
    "HasNoVertexMismatch",
    "VertexTrackMultiplicity",
    "HasNoBackExitingTracks",
    "phys_n_dead_discr_pair_upstream_prim_track_proj",
    
    // other
    "Psi",
    "part_response_total_recoil_passive_allNonMuonClusters_id",
    "part_response_total_recoil_passive_allNonMuonClusters_od",
  };

  //Dummy initialization value, I should never see this value in the output tuple and if I do something got messed up
  int selectionCategory = -60;

  // real selectionCategory values:
  // -999 = signal (green)
  // 0 = NonQE (yellow) bkgd
  // 1 = other NuE CC (red) bkgd
  // 2 = NC Pi0 (pink) bkgd
  // 3 = NuMu CC Pi0 (teal) bkgd
  // -1 = other backgrounds (dark blue)

  //these are to save the results of sideband cuts as an array of ints, matching the order of entries in 
  int n_sideband_cuts;
  TString sideband_names_str;
  std::vector<int> sidebandCutResults;
  // An array of those same bits, to fill the branches
  
  void Init(const std::string& filename, bool enable, TTree& inputTree, const std::vector<std::string>& sideband_names) {
    enableOutputTree = enable;
    if (!enableOutputTree) return;

    n_sideband_cuts = sideband_names.size();
    outputTreeFile = TFile::Open(filename.c_str(), "RECREATE");
    if (!outputTreeFile || outputTreeFile->IsZombie()) {
      throw std::runtime_error("Failed to open output file: " + filename);
    }    

    //keep only branches I use, helps reduce run time & file size
    inputTree.SetBranchStatus("*", 0);
    for (const auto& b : branches_to_keep) {
      inputTree.SetBranchStatus(b.c_str(), 1); 
    }
    inputTree.SetBranchStatus("*", 1);

    sideband_names_str = "";
    for (int i = 0; i < sideband_names.size(); i++) {
      sideband_names_str += sideband_names[i];
      if (i < sideband_names.size()-1) sideband_names_str += ", ";
    }

    sidebandCutResults.resize(n_sideband_cuts, 0);
    outputTree = inputTree.CloneTree(0);
    outputTree->SetName("MasterAnaDev");
    outputTree->Branch("selectionCategory", &selectionCategory, "selectionCategory/I");
    std::string branch_def = "sideband_cut_results[" + std::to_string(n_sideband_cuts) + "]/I";
    //outputTree->Branch("sideband_cut_results", &cutBits, branch_def.c_str());
    outputTree->Branch("sideband_cut_results", sidebandCutResults.data(), branch_def.c_str());

    inputTree.SetBranchStatus("*", 1); //unfortunately I think I need this, would run A LOT faster if I could keep those off though...
  }
  
  void Fill(const std::string& sidebandName, Long64_t entryNumber) {
    if (!enableOutputTree) return;
    auto& filledEntries = filledEntriesPerTree[sidebandName];
    if (filledEntries.count(entryNumber)) return; // Already filled — skip

    outputTree->Fill();
    filledEntries.insert(entryNumber);
  }
  
  void WriteAndClose() {
    if (!enableOutputTree) return;

    if (outputTreeFile) {

      TNamed* sideband_meta = new TNamed("sideband_cut_names", sideband_names_str.Data());
      
      outputTreeFile->cd();
      outputTree->Write();
      sideband_meta->Write();
      outputTreeFile->Close();
      delete outputTreeFile;
      outputTreeFile = nullptr;
      filledEntriesPerTree.clear();
    }
  }
};

extern OutputTreeManager g_OutputTreeManager;

struct CCNuEEvent {
    int m_idx; // Index for Best Michel in nmichels
    double m_bestdist; // in mm 
    std::vector<double> m_best2D; //0: XZ, 1: UZ, 2:VZ   
    double m_best_XZ;
    double m_best_UZ;
    double m_best_VZ;

    Long64_t entryNumber; //entry number in the LoopAndFillMC event loop, need it to make my output trees not have triple entries :(
    //std::vector<Michel*> m_nmichels; //nmatched michels
};
#endif
