#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"

#include "TH1D.h"
#include "TFile.h"
#include "TKey.h"
#include "TParameter.h"
#include "TCanvas.h"


//c++ includes
#include <iostream>
#include <string.h>
#include <bits/stdc++.h>
#include <vector>
#include <array>
using namespace std;
#include <cmath>

namespace PlotUtils
{
  class MnvH1D;
}

//Using Arrays

bool hasFSProton(int PDGs[183], int nFS_part) {
  bool hasProton = false;
  for (int i = 0; i < nFS_part; i++) {
    int part = abs(PDGs[i]);
    if (part==2212){
      hasProton=true;
    }
  }
  return hasProton;
}
bool hasTrackableFSProton(int PDGs[183], double FS_energies[183], int nFS_part) {
  bool hasTrackableProton = false;
  for (int i = 0; i < nFS_part; i++) {
    int part = abs(PDGs[i]);
    if (part==2212 && FS_energies[i] > 1038 ){
      hasTrackableProton=true;
    }
  }
  return hasTrackableProton;
}


void truthStudies() {
  
  //TChain tree("MasterAnaDev");
  //tree.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110030_Playlist.root");
  ///pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110030_Playlist.root
  
  TFile * myFile = new TFile("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110030_Playlist.root");
  TTree *tree = (TTree*)myFile->Get("MasterAnaDev");

  tree->SetBranchStatus("*", false);
  
  //Truth quantities
  tree->SetBranchStatus("mc_incoming", true);
  tree->SetBranchStatus("mc_nFSPart", true);
  tree->SetBranchStatus("mc_FSPartPDG", true);
  tree->SetBranchStatus("mc_FSPartE", true);
  tree->SetBranchStatus("mc_run", true);
  tree->SetBranchStatus("mc_subrun", true);
  tree->SetBranchStatus("mc_nthEvtInFile", true);

  //MAT reco cuts for reference
  /*  
  preCuts.emplace_back(new reco::HasTracks<CVUniverse, CCNuEEvent>());
  preCuts.emplace_back(new reco::NoVertexMismatch<CVUniverse, CCNuEEvent>());
  preCuts.emplace_back(new reco::HasNoBackExitingTracks<CVUniverse, CCNuEEvent>());
  preCuts.emplace_back(new reco::EMLikeTrackScore<CVUniverse, CCNuEEvent>());

  preCuts.emplace_back(new reco::DSCalVisE<CVUniverse, CCNuEEvent>());
  preCuts.emplace_back(new reco::ODCalVisE<CVUniverse, CCNuEEvent>());
  preCuts.emplace_back(new reco::VertexTrackMultiplicity<CVUniverse, CCNuEEvent>());
  preCuts.emplace_back(new reco::Afterpulsing<CVUniverse, CCNuEEvent>());
  preCuts.emplace_back(new reco::NoDeadtime<CVUniverse, CCNuEEvent>(1, "Deadtime"));

  preCuts.emplace_back(new reco::StartPointVertexMultiplicity<CVUniverse, CCNuEEvent>());
  preCuts.emplace_back(new reco::MeanFrontdEdX<CVUniverse, CCNuEEvent>());
  preCuts.emplace_back(new reco::NonMIPClusterFraction<CVUniverse, CCNuEEvent>());
  preCuts.emplace_back(new reco::TransverseGapScore<CVUniverse, CCNuEEvent>());

  preCuts.emplace_back(new reco::ZRange<CVUniverse, CCNuEEvent>("Tracker", minZ, maxZ));
  preCuts.emplace_back(new reco::Apothem<CVUniverse, CCNuEEvent>(apothem));
  preCuts.emplace_back(new reco::Eavailable<CVUniverse, CCNuEEvent>());

  preCuts.emplace_back(new reco::ProtonEnd<CVUniverse, CCNuEEvent>());
  */

  //Reco branches
  tree->SetBranchStatus("MasterAnaDev_proton_endPointZ", true);
  tree->SetBranchStatus("HasNoBackExitingTracks", true);
  tree->SetBranchStatus("prong_part_score", true);
  tree->SetBranchStatus("prong_HCALVisE", true);
  tree->SetBranchStatus("prong_ECALVisE", true);
  tree->SetBranchStatus("prong_ODVisE", true);
  tree->SetBranchStatus("prong_SideECALVisE", true);
  tree->SetBranchStatus("HasNoVertexMismatch", true);
  tree->SetBranchStatus("VertexTrackMultiplicity", true);
  tree->SetBranchStatus("prong_dEdXMeanFrontTracker", true);
  tree->SetBranchStatus("prong_TransverseGapScore", true);
  //tree->SetBranchStatus("VertexTrackMultiplicity", true);
  //tree->SetBranchStatus("VertexTrackMultiplicity", true);
  //tree->SetBranchStatus("VertexTrackMultiplicity", true);
  //tree->SetBranchStatus("VertexTrackMultiplicity", true);


  int inc, run, subrun, gate, nFS_part, noBackExitingTracks,noVertexMismatch,vertexTrackMultiplicity;
  double recoProton, HCalVisE, ECalVisE,ODVisE,SideCalVisE,frontdEdX,transverseGapScore;
  
  Int_t FS_PDG[183];
  Double_t FS_E[183];

  //vector<int> FS_PDG {0, 1, 2, 3, 4, 5};
  //vector<double> FS_E {27.3, 95.3, 2138.7, 8.9, 17.3, 82.4};
   

  tree->SetBranchAddress("mc_incoming", &inc);
  tree->SetBranchAddress("mc_nFSPart", &nFS_part);
  tree->SetBranchAddress("mc_FSPartPDG", FS_PDG);
  tree->SetBranchAddress("mc_FSPartE", FS_E);
  tree->SetBranchAddress("mc_run", &run);
  tree->SetBranchAddress("mc_subrun", &subrun);
  tree->SetBranchAddress("mc_nthEvtInFile", &gate); //off by one vs arachne gate, I think this number + 1 = arachne gate num

  tree->SetBranchAddress("MasterAnaDev_proton_endPointZ", &recoProton);
  tree->SetBranchAddress("HasNoBackExitingTracks", &noBackExitingTracks);
  tree->SetBranchAddress("prong_HCALVisE", &HCalVisE);
  tree->SetBranchAddress("prong_ECALVisE", &ECalVisE);
  tree->SetBranchAddress("prong_ODVisE", &ODVisE);
  tree->SetBranchAddress("prong_SideECALVisE", &SideCalVisE);
  tree->SetBranchAddress("HasNoVertexMismatch", &noVertexMismatch);
  tree->SetBranchAddress("VertexTrackMultiplicity", &vertexTrackMultiplicity);
  tree->SetBranchAddress("prong_dEdXMeanFrontTracker", &frontdEdX);
  tree->SetBranchAddress("prong_TransverseGapScore", &transverseGapScore);

  int NuMu_RecoProtons = 0;
  int NuMu_TrueProtons = 0;
  int NuE_RecoProtons = 0;
  int NuE_TrueProtons = 0;
  int count = 0;

  cout << "skrt" << endl;
  for (int i = 0; tree->LoadTree(i) >= 0; i++){
    tree->GetEntry(i);

    //for diagnostivcs
    /*
    if (i < 10){
      cout << "size of pdg array: " << sizeof(FS_E)/sizeof(double) << endl;
      int index = 0;
      for (int energy : FS_E){
	cout << "contents of pdg array[" << index << "]: " << energy << endl;
	index++;
	if (index >= nFS_part){ break;}
      }
    }
    */

    //if (inc==14 && recoProton > 0){ NuMu_RecoProtons++;    };
    //if (inc==14 && hasTrackableFSProton(FS_PDG, FS_E, nFS_part)){ NuMu_TrueProtons++; };
    //if (inc==12 && recoProton > 0){ NuE_RecoProtons++; };
    //if (inc==12 && hasTrackableFSProton(FS_PDG, FS_E, nFS_part)){ NuE_TrueProtons++; }
 
    //OK let's try and redo those last two but without the use of mc_incoming... aka the most barebones version of my selection
    double DSCalVisE = HCalVisE/ECalVisE;
    double ODCalVisE = ODVisE/SideCalVisE;  

    //if (noBackExitingTracks==1 && DSCalVisE <= 0.2 && ODCalVisE <= 0.05 && noVertexMismatch==1 && vertexTrackMultiplicity < 6 && frontdEdX < 2.4 && transverseGapScore > 15 && recoProton > 0) { NuE_RecoProtons++; }

    //if (noBackExitingTracks==1 && DSCalVisE <= 0.2 && ODCalVisE <= 0.05 && noVertexMismatch==1 && vertexTrackMultiplicity < 6 && frontdEdX < 2.4 && transverseGapScore > 15 && hasTrackableFSProton(FS_PDG, FS_E, nFS_part)) { NuE_TrueProtons++; }

    
    if (noBackExitingTracks==1 && DSCalVisE <= 0.2 && ODCalVisE <= 0.05 && noVertexMismatch==1 && vertexTrackMultiplicity < 6 && frontdEdX < 2.4 && transverseGapScore > 15) { 
      count++;
      cout << "count : " << count << ", recoProton = " << recoProton << ", nuE_recoprotons = " << NuE_RecoProtons << endl;
      if (recoProton > 0) {
	cout << "WENT IN" << endl;
	NuE_RecoProtons = NuE_RecoProtons + 1; 
      }
      else { cout << "DIDN'T GO IN" << endl; }

      /*
      if (hasTrackableFSProton(FS_PDG, FS_E, nFS_part)) {
	NuE_TrueProtons++;
	}*/
      //double andrea = 12.7+19.6;
    }
    
  

  }
  cout << "numu reco proton events: " << NuMu_RecoProtons << endl;
  cout << "numu true proton events: " << NuMu_TrueProtons << endl;
  cout << "nue reco proton events: " << NuE_RecoProtons << endl;
  cout << "nue true proton events: " << NuE_TrueProtons << endl;

  tree->ResetBranchAddresses();
}

