#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include <vector>

void ScanFSParticles() {
  TFile *file = TFile::Open("/exp/minerva/data/users/cpernas/NuE_TKI/MC_TUPLE_WITH_MICHELCUT_Oct_18_2024.root");
  TTree *tree = (TTree*)file->Get("MasterAnaDev");

  //int nEventsWithPDG101 = 0;
  bool printEvent = false;
  int nEventsWithFSpi0 = 0;
  int nEventsWithFSpiPlus = 0;
  int nEventsWithFSpiMinus = 0;
  int nEventsWithoutChargedPions = 0;
  int nEventsWithoutPions = 0;
  int nEventsWithFSProton = 0;
  int nEventsWithoutFSProtonOrChargedPion = 0;
  int totalCount = 0;
  int withoutPiPlusCount = 0;

  
  int mc_incoming;
  int mc_current;
  int selectionCategory;
  int mc_nFSPart;
  int mc_FSPartPDG[1000];
  int improved_nmichel;



  
  tree->SetBranchAddress("mc_incoming", &mc_incoming);
  tree->SetBranchAddress("mc_current", &mc_current);
  tree->SetBranchAddress("selectionCategory", &selectionCategory);
  tree->SetBranchAddress("mc_nFSPart", &mc_nFSPart);
  tree->SetBranchAddress("mc_FSPartPDG", &mc_FSPartPDG);
  tree->SetBranchAddress("improved_nmichel", &improved_nmichel);
  
  std::map<int, int> pdgCounts;

  //loop through events
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (selectionCategory!=1) { continue; } //skip signal events
    totalCount++;
    bool haspi0 = false;
    bool haspiplus = false;
    bool haspiminus = false;
    bool hasproton = false;
    //loop through final state particles
    for (int j = 0; j < mc_nFSPart; j++){
      //std::cout << mc_FSPartPDG[j] << ", ";
      int pdgCode = mc_FSPartPDG[j];
      pdgCounts[pdgCode]++;
      
      if (mc_FSPartPDG[j] == 111 && !haspi0){
	nEventsWithFSpi0++;
	haspi0=true;
      }

      if (mc_FSPartPDG[j] == 211 && !haspiplus && improved_nmichel==0){
	nEventsWithFSpiPlus++;
	haspiplus=true;
      }      

      if (mc_FSPartPDG[j] == -211 && !haspiminus){
	nEventsWithFSpiMinus++;
	haspiminus=true;
      }

      if (mc_FSPartPDG[j] == 2212 && !hasproton){
	nEventsWithFSProton++;
	hasproton=true;
      }

    }
    if (!haspiminus && !haspiplus){
      nEventsWithoutChargedPions++;
      printEvent=true;
    }
    if (!haspiminus && !haspiplus && !haspi0){
      nEventsWithoutPions++;
    }
    if (!haspiminus && !haspiplus && !hasproton){
      nEventsWithoutFSProtonOrChargedPion++;
    }
    if (printEvent){
      std::cout << "-------------- Event #" << i << ", category = " << selectionCategory << " --------------" << std::endl;
      std::cout << "PDG codes of the final state: ";
      for (int j = 0; j < mc_nFSPart; j++){
	std::cout << mc_FSPartPDG[j] << ", ";
      }
      std::cout << std::endl;
      std::cout << std::endl;    
      
    }
    printEvent=false;
  }
  std::cout << "PDG code counts:" << std::endl;
  for (const auto& pair : pdgCounts) {
    std::cout << "PDG Code: " << pair.first << " - Count: " << pair.second << std::endl;
  }
  std::cout << "Total number of background events: " << totalCount << std::endl;
  std::cout << "Total number of background events with final state pi0's: " << nEventsWithFSpi0 << std::endl;
  std::cout << "Total number of background events with final state pi+: " << nEventsWithFSpiPlus << std::endl;
  std::cout << "Total number of background events with final state pi-: " << nEventsWithFSpiMinus << std::endl;
  std::cout << "Total number of background events without any pi+ OR pi-: " << nEventsWithoutChargedPions << std::endl;
  std::cout << "Total number of background events without any pi+ OR pi- OR pi0: " << nEventsWithFSpiMinus << std::endl;
  std::cout << "Total number of background events without any pi+ OR pi- OR proton: " << nEventsWithoutFSProtonOrChargedPion << std::endl;
  std::cout << "Total number of background events with final state protons: " << nEventsWithFSProton << std::endl;
  
  
}
