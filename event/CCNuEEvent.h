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

  //TTree* outputTree = nullptr;
  std::map<std::string, TTree*> outputTrees; //One for signal region (ie passing all cuts), then one per sideband
  
  //This is so dumb, but basically maps each sideband name to a set of entry numbers (which corresponds to i in my main runEventLoop)
  //Avoids the triple event filling thing into my trees but jesus this struct is getting out of hand...
  std::map<std::string, std::unordered_set<Long64_t>> filledEntriesPerTree; 

  //Dummy initialization value, I should never see this value in the output tuple and if I do something got messed up
  int selectionCategory = -60;

  //This is pretty temporary, and will most likely change frequently. I want to save a tuple of events, including many that don't pass the full selection
  //But save what WOULD have been the results of the cut. I was originally thinking of doing this as a bitset, but since I know I'll only need 6 (or fewer)
  //bits, I guess I'll just save it as 6 additional branches, whatever. 
  static constexpr int numSidebands = 6;
  // Bitset storing cut results
  std::bitset<numSidebands> sidebandCutResults;
  // An array of those same bits, to fill the branches
  //bool cutBits[numSidebands];
  
  void Init(const std::string& filename, bool enable, TTree& inputTree, const std::vector<std::string>& sidebandNames) {
    enableOutputTree = enable;
    if (!enableOutputTree) return;
    
    outputTreeFile = TFile::Open(filename.c_str(), "RECREATE");
    if (!outputTreeFile || outputTreeFile->IsZombie()) {
      throw std::runtime_error("Failed to open output file: " + filename);
    }
    
    for (const auto& name : sidebandNames) {
      TTree *tree = inputTree.CloneTree(0);
      tree->SetName(name.c_str());
      //TTree* tree = new TTree(name.c_str(), name.c_str());
      //TTree* tree = CloneTreeWithNewName(inputTree, name);
      tree->Branch("selectionCategory", &selectionCategory, "selectionCategory/I");
      //This little loop here adds 6 branches to each of the 6 TTrees, showing the status of all 6 of the cuts
      /*
      for (int i = 0; i < numSidebands; ++i) { //skip the first entry, which is the signal region/passes precuts tree, and doesn't need a branch
	std::string branchName = "default";
	if (i==0) { branchName = "passes_em_score_cut"; }
	else if (i==1) { branchName = "passes_michel_cut"; }
	else if (i==2) { branchName = "passes_E_lep_cut"; }
	else if (i==3) { branchName = "passes_Mod_Eavail_cut"; }
	else if (i==4) { branchName = "passes_mean_front_DEDX_cut"; }
	else if (i==5) { branchName = "passes_ESC_proton_cut"; }
	tree->Branch(branchName.c_str(), &cutBits[i]);
      }
      */
      outputTrees[name] = tree;
    }
  }
  
  void Fill(const std::string& sidebandName, Long64_t entryNumber) {
    if (!enableOutputTree || !outputTrees.count(sidebandName)) return;
    auto& filledEntries = filledEntriesPerTree[sidebandName];
    if (filledEntries.count(entryNumber)) return; // Already filled — skip
    /*
    for (int i = 0; i < numSidebands; ++i) {
    cutBits[i] = sidebandCutResults[i];
    }*/
    outputTrees[sidebandName]->Fill();
    filledEntries.insert(entryNumber);
  }
  
  void WriteAndClose() {
    if (!enableOutputTree) return;

    if (outputTreeFile) {
      outputTreeFile->cd();
      for (auto& [name, tree] : outputTrees) {
        tree->Write();
      }
      outputTreeFile->Close();
      delete outputTreeFile;
      outputTreeFile = nullptr;
      outputTrees.clear();
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
