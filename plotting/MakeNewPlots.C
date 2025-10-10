//This script is basically to make any kind of custom plots I want to make, out of the reduced selected events only tuple that I make from my event loop. 


#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <TH2F.h>
#include <TCanvas.h>

#include <vector>

double CalcEavail(double E_tracker, double E_ecal, double reco_lep_E) {
  double E = E_tracker + E_ecal;
  double Eavail = (E * 1.17 - ((0.008 * reco_lep_E) + 5))/1000.;
  //std::cout << "Reco Eavail: " << Eavail << std::endl;                                                                                                            
  return Eavail;
}

void MakeNewPlots() {
  TFile *file = TFile::Open("/exp/minerva/data/users/cpernas/NuE_TKI/MC_TUPLE_Oct_18_2024.root");
  TTree *tree = (TTree*)file->Get("MasterAnaDev");

  bool printEvent = false;
  int mc_incoming, mc_current, selectionCategory, mc_nFSPart, improved_nmichel;
  int mc_FSPartPDG[1000];  

  double blob_recoil_E_tracker, blob_recoil_E_ecal;
  
  tree->SetBranchAddress("mc_incoming", &mc_incoming);
  tree->SetBranchAddress("mc_current", &mc_current);
  tree->SetBranchAddress("selectionCategory", &selectionCategory);
  tree->SetBranchAddress("mc_nFSPart", &mc_nFSPart);
  tree->SetBranchAddress("mc_FSPartPDG", mc_FSPartPDG);
  tree->SetBranchAddress("blob_recoil_E_tracker", &blob_recoil_E_tracker);
  tree->SetBranchAddress("blob_recoil_E_ecal", &blob_recoil_E_ecal);
  
  TBranch* prongE = tree->GetBranch("prong_part_E");
  prongE->GetEntry(0); //next line complains abt null pointer if I don't load an entry, this gets overwritten @ beginning of loop anyways        
  vector<vector<double>>* prongEarr = (vector<vector<double>>*)(prongE->GetLeaf("prong_part_E")->GetValuePointer());  
  
  std::map<int, int> pdgCounts;
  
  // Create a 2D histogram
  int nBinsX = 40; // Number of bins for X axis, which is elep
  int nBinsY = 40; // Number of bins for Y axis, which is eavail
  TH2F *h1 = new TH2F("h1", "2D Histogram", nBinsX, 2, 10, nBinsY, 0, 2); // reco
  //TH2F *h2 = new TH2F("h1", "2D Histogram", nBinsX, 0, 20, nBinsY, 0, 4); // true
  
  //loop through events
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (selectionCategory!=5) { continue; } //lets me select out which category I wanna look at
    
    double reco_lep_E = prongEarr->data()[0][3];
    double reco_lep_E_GeV = reco_lep_E/1000.;
    double reco_Eavail = CalcEavail(blob_recoil_E_tracker, blob_recoil_E_ecal, reco_lep_E);

    h1->Fill(reco_lep_E_GeV, reco_Eavail);
    //h2->Fill();
  }
  // Draw the histogram
  TCanvas *c1 = new TCanvas("c1", "2D Histogram", 800, 600);
  h1->Draw("colz");
  c1->SaveAs("reco_elep_vs_evail_Numu_CC_pi0.png");
  delete c1;
  delete h1;

  //TCanvas *c2 = new TCanvas("c2", "2D Histogram", 800, 600);
  //h2->Draw("COLZ");
  //c2->SaveAs("true_elep_vs_evail.png");
  //delete c2;
  //delete h2;

  file->Close();
}

