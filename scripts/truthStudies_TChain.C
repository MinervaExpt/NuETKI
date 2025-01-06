#include "TFile.h"
#include "TH1F.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include <vector>

double CalcApothem(double xval, double yval) {
  double x=abs(xval);
  double y=abs(yval);
  if ( x == 0 or y/x > 1/sqrt(3)) {
    return (y+x/sqrt(3))/2*sqrt(3);
  }
  else {
    return x;
  }
}

void truthStudies_TChain() {

  
  TChain* tree = new TChain("MasterAnaDev");
  //tree->Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110030*.root");
  //tree->Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110031*.root");
  tree->Add("/pnfs/minerva/scratch/users/cpernas/ProtonEfficiencyStudy/run_110040_newShortTracking_test/grid/central_value/minerva/ana/v22r1p1/00/11/00/40/SIM*.root");
  cout << "# of files added: " << tree->GetNtrees() << ", # of events: " << tree->GetEntries() << endl;
  

  //TTree *tree = (TTree*)myFile->Get("MasterAnaDev");
  TTreeReader myReader(tree);

  //True variables
  TTreeReaderValue<Int_t> mc_incoming(myReader, "mc_incoming");
  TTreeReaderValue<Int_t> run(myReader, "mc_run");
  TTreeReaderValue<Int_t> subrun(myReader, "mc_subrun");
  TTreeReaderValue<Int_t> nthEvtInFile(myReader, "mc_nthEvtInFile");
  TTreeReaderArray<Int_t> FS_PDG(myReader, "mc_FSPartPDG");
  TTreeReaderArray<Double_t> FS_E(myReader, "mc_FSPartE");

  //Reco Variables
  TTreeReaderArray<Double_t> vtx(myReader, "vtx");
  TTreeReaderValue<Int_t> nProngs(myReader, "n_prongs");
  TTreeReaderValue<Int_t> noBackExitingTracks(myReader, "HasNoBackExitingTracks");
  TTreeReaderValue<Int_t> noVertexMismatch(myReader, "HasNoVertexMismatch");
  TTreeReaderValue<Int_t> trackMultiplicity(myReader, "VertexTrackMultiplicity");
  TTreeReaderValue<Int_t> vertexMultiplicity(myReader, "StartPointVertexMultiplicity");
  TTreeReaderValue<Double_t> blobRecoilTracker(myReader, "blob_recoil_E_tracker");
  TTreeReaderValue<Double_t> blobRecoilECAL(myReader, "blob_recoil_E_ecal");
  TTreeReaderValue<Int_t> tDead(myReader, "phys_n_dead_discr_pair_upstream_prim_track_proj");
  TTreeReaderArray<Double_t> showerScore(myReader, "prong_part_score");
  TTreeReaderArray<Double_t> firstFireFraction(myReader, "prong_FirstFireFraction");
  TTreeReaderArray<Double_t> HCALVisE(myReader, "prong_HCALVisE");
  TTreeReaderArray<Double_t> ECALVisE(myReader, "prong_ECALVisE");
  TTreeReaderArray<Double_t> ODVisE(myReader, "prong_ODVisE");
  TTreeReaderArray<Double_t> sideECALVisE(myReader, "prong_SideECALVisE");
  TTreeReaderArray<Double_t> meanFrontDEDX(myReader, "prong_dEdXMeanFrontTracker");
  TTreeReaderArray<Double_t> transverseGapScore(myReader, "prong_TransverseGapScore");
  TTreeReaderArray<Double_t> nonMIPClusFrac(myReader, "prong_NonMIPClusFrac");
  TTreeReaderValue<Double_t> MAD_proton_endPointZ(myReader, "MasterAnaDev_proton_endPointZ");
  TTreeReaderArray<Int_t> prong_part_PID(myReader, "prong_part_pid");

  //This is kinda hacky, but by doing it this way I'm basically doing a manual SetBranchAddress, and also don't have to load the whole tree (just the 1 branch)
  //kinda the same way you only load a subset of branches by setting branch status on & off, but not sure if that'll mess w my TTReeReader and this works so 
  // ¯\_(ツ)_/¯
  //tree->GetEntry(0);
  //TBranch* prongE = tree->GetBranch("prong_part_E");
  //prongE->GetEntry(0); //next line complains abt null pointer if I don't load an entry, this gets overwritten @ beginning of loop anyways
  //vector<vector<double>>* prongEarr = (vector<vector<double>>*)(prongE->GetLeaf("prong_part_E")->GetValuePointer());
  
  //TFile* outFile = new TFile("/minerva/data/users/cpernas/NuE_TKI/NuE_Selection_ONLY.root","recreate");
  //TTree *newTree = tree->CloneTree(0);
  //tree->SetBranchStatus("*", false);
  //tree->SetBranchStatus("prong_part_E", true);

  int counter = 0;
  int nTrueProtonsBeforeReco = 0;
  int nRecoProtons = 0;
  int nTrueProtonsAfterReco = 0;
  
  int numuRecoProton = 0;
  int numuTrueProton = 0;
  int nueRecoProton = 0;
  int nueTrueProton = 0;
  
  int unsuccessfulCount = 0;
  int successfulCount = 0;
  
  Int_t i = 0;
  while (myReader.Next()) {
    //if (i > 10000){ break; }
    //if (i%5000==0) {cout << i << endl;}
    
    //tree->GetEntry(i);
    //prongE = tree->GetBranch("prong_part_E");
    //prongE->GetTree()->GetEntry(i);
    //prongE->GetTree();
    //prongE->GetEntry(i);

    tree->GetEntry(i);
    TBranch* prongE = tree->GetBranch("prong_part_E");
    prongE->GetEntry(i); //next line complains abt null pointer if I don't load an entry, this gets overwritten @ beginning of loop anyways
    vector<vector<double>>* prongEarr = (vector<vector<double>>*)(prongE->GetLeaf("prong_part_E")->GetValuePointer());

    double apothem = CalcApothem(vtx[0], vtx[1]);
    double reco_lep_E = prongEarr->data()[0][3];
    double Eavail = ((*blobRecoilTracker + *blobRecoilECAL)*1.17 - ((0.008 * reco_lep_E)+5))/1000.;
    
    double DSCalRatio;
    double ODCalRatio;
    if (HCALVisE[0]==0 && ECALVisE[0]==0) { DSCalRatio = 0; }
    else { DSCalRatio = HCALVisE[0] / ECALVisE[0]; } 
    if (ODVisE[0]==0 && sideECALVisE[0]==0) { ODCalRatio = 0;}
    else { ODCalRatio = ODVisE[0] / sideECALVisE[0]; }
    
    bool hasTrueProton = false;
    for (int j = 0; j < FS_PDG.GetSize(); j++) {
      int part = abs(FS_PDG[j]);
      if (part==2212 && FS_E[j] > 1038 ){
	//if (part==2212){
	hasTrueProton=true;
      }
    }
    
    //if (i==0){ //for investigating specific events
        //cout << "nProngs: " << *nProngs << ", noVertexMismatch: " << *noVertexMismatch << ", noBackExitingTracks: " << *noBackExitingTracks << ", showerScore[0]: " << showerScore[0] << ", DSCalRatio: " << DSCalRatio << ", ODCalRatio: " << ODCalRatio << ", trackMultiplicity: " << *trackMultiplicity << ", firstFireFraction[0]: " << firstFireFraction[0] << ", tDead: " << *tDead << ", vertexMultiplicity: " << *vertexMultiplicity << ", meanFrontdedx[0]: " << meanFrontDEDX[0] << ", nonMIPClusFrac[0]: " << nonMIPClusFrac[0] << ", transverseGapScore[0]: " << transverseGapScore[0] << ", vtx[2]: " << vtx[2] << ", apothem: " << apothem << ", Eavail: " << Eavail << endl;
    //}

    //My "precuts", if you will, basically all cuts to get nu_e events but nothing about the proton yet
    if (*nProngs>0 && *noVertexMismatch==1 && *noBackExitingTracks==1 && showerScore[0]>0.7 && DSCalRatio<=0.2 && ODCalRatio<=0.05 && *trackMultiplicity<6 && firstFireFraction[0]>=0.25 && *tDead<=1 && *vertexMultiplicity==1 && meanFrontDEDX[0]<2.4 && nonMIPClusFrac[0]>0.4 && transverseGapScore[0]>15 && vtx[2]>5980 && vtx[2]<8422&& apothem<850 && Eavail<0.8){
      counter++;
      //tree->GetTree()->GetEntry(i);
      //newTree->Fill();
      
      if (hasTrueProton) { 	
	nTrueProtonsBeforeReco++;  
	if (*MAD_proton_endPointZ<=0){
	  cout << "i: " << i << " - UNSUCCESSFULLY reco'd proton at run: " << *run << ", subrun: " << *subrun << ", nthEvtInFile: " << *nthEvtInFile << endl;
          //cout << "MasterAnaDev_proton_endPointZ: " << *MAD_proton_endPointZ << endl;
	  unsuccessfulCount++;
	}
      }

      if (*MAD_proton_endPointZ>0){ 
	nRecoProtons++; 
	if (hasTrueProton) { 
	  nTrueProtonsAfterReco++;
	  cout << "i: " << i << " - SUCCESSFULLY reco'd proton at run: " << *run << ", subrun: " << *subrun << ", nthEvtInFile: " << *nthEvtInFile << endl;
	  successfulCount++;
	  //cout << "MasterAnaDev_proton_endPointZ: " << *MAD_proton_endPointZ << endl;
	}
      }
      
    }
    

    /*
    if (*mc_incoming == 14 && *MAD_proton_endPointZ > 0){ numuRecoProton++; }
    if (*mc_incoming == 14 && hasTrueProton){ numuTrueProton++; }
    if (*mc_incoming == 12 && *MAD_proton_endPointZ > 0){ nueRecoProton++; }
    if (*mc_incoming == 12 && hasTrueProton){ nueTrueProton++; }
    */

    hasTrueProton = false;
    i++;
  }
  //newTree->Write();
  //outFile->Write();
  //outFile->Close();

  cout << "total # events: " << i << endl;

  //cout << "numu reco proton events: " << numuRecoProton << endl;
  //cout << "numu true proton events: " << numuTrueProton << endl;
  //cout << "nue reco proton events: " << nueRecoProton << endl;
  //cout << "nue true proton events: " << nueTrueProton << endl;
  
  cout << "Unsuccessful count: " << unsuccessfulCount << endl;
  cout << "Successful count: " << successfulCount << endl;

  cout << "# events BEFORE checking for reco protons: " << counter << endl;
  cout << "# events with true protons BEFORE checking for reco protons: " << nTrueProtonsBeforeReco << endl;
  cout << "# events AFTER checking for reco protons: " << nRecoProtons << endl;
  cout << "# events with true protons AFTER checking for reco protons: " << nTrueProtonsAfterReco << endl;
}
