#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TMath.h>
#include <Math/Vector3D.h>
#include <Math/RotationX.h>
#include <iostream>
#include <vector>

void PlotElectronResolution() {
  //TString filename = "/exp/minerva/data/users/cpernas/TUPLE_MC_May_14_2025_me1M_all_cuts_with_Michel_sideband.root";
  TString filename = "/exp/minerva/data/users/cpernas/TUPLE_MC_May_12_2025_NuEOnly_me1M_all_cuts_with_Michel_sideband.root";
  //TString filename = "/pnfs/minerva/persistent/users/cpernas/default_analysis_loc/TUPLE_MC_May_01_2025_me1B_all_cuts_with_Michel_sideband.root";

  TFile* file = TFile::Open(filename);
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file\n";
    return;
  }

  TTree* tree = (TTree*)file->Get("Signal_Region"); // <- Update this
  if (!tree) {
    std::cerr << "Error: TTree not found\n";
    return;
  }

  // electron mass (in MeV)
  const double M_e = 0.51099895;
  
  // Branch variables
  int selectionCategory, mc_run, mc_subrun, mc_nthEvtInFile;
  double trueElectron4Vec[4];
  std::vector<std::vector<double>> *prong_part_E = nullptr;
  
  
  tree->SetBranchAddress("selectionCategory", &selectionCategory);
  tree->SetBranchAddress("mc_run", &mc_run);
  tree->SetBranchAddress("mc_subrun", &mc_subrun);
  tree->SetBranchAddress("mc_nthEvtInFile", &mc_nthEvtInFile);
  tree->SetBranchAddress("mc_primFSLepton", &trueElectron4Vec);
  tree->SetBranchAddress("prong_part_E", &prong_part_E);


  // Create histograms
  TH2F* h1 = new TH2F("h1", "Electron p Resolution;True Electron P (GeV); (Reco - True)/True", 50, 0, 20, 50, -0.5, 0.5);
  TH2F* h2 = new TH2F("h2", "Electron pT Resolution;True Electron pT (GeV); (Reco - True)/True", 50, 0, 2, 50, -0.5, 0.5);
  TH2F* h3 = new TH2F("h3", "Electron Theta Resolution;True Electron Theta (deg); (Reco - True)", 50, 0, 20, 50, -10, 10);

  Long64_t nentries = tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    //if (i>1000) break;
    tree->GetEntry(i);

    if (selectionCategory != -999) continue;
    //std::cout << "\nlooking at event: " << i << "\n" << std::endl;
    std::vector<double> recoElectron4Vec = (*prong_part_E)[0];
    
    //calc true kinematics
    double trueE = trueElectron4Vec[3];
    double trueP = std::sqrt(trueE * trueE - M_e * M_e) / 1000.0; // in GeV
    double px = trueElectron4Vec[0];
    double py = trueElectron4Vec[1];
    double pz = trueElectron4Vec[2];
    ROOT::Math::XYZVector p_vec(px, py, pz);
    ROOT::Math::RotationX r(-3.3 * (TMath::Pi() / 180.)); //represents a rotation about the x axis of 3.3 degrees, aka 3.3 degrees down toward floor to match beam direction
    double trueTheta = (r(p_vec)).Theta();
    double trueThetaDeg = trueTheta * 180.0 / TMath::Pi();
    double truePt = trueP * std::sin(trueTheta);

    //reco kinematics
    double recoE = recoElectron4Vec[3];
    double recoP = sqrt((recoE * recoE - M_e * M_e))/ 1000.0;
    double px2 = recoElectron4Vec[0];
    double py2 = recoElectron4Vec[1];
    double pz2 = recoElectron4Vec[2];
    ROOT::Math::XYZVector p2_vec(px2, py2, pz2);
    double recoTheta = (r(p2_vec)).Theta();
    double recoThetaDeg = recoTheta * 180.0 / TMath::Pi();
    double recoPt = recoP * std::sin(recoTheta);

    
    double resolutionP = (recoP - trueP) / trueP;
    double resolutionPt = (recoPt - truePt) / truePt;
    double resolutionTheta = (recoThetaDeg - trueThetaDeg); //not over theta...? Mike's suggestion
    
    //std::cout << "Filling true P:  " << trueP << ", Resolution: " << resolutionP << std::endl;
    //std::cout << "Filling true Pt:  " << truePt << ", resolutionPt: " << resolutionPt << std::endl;
    //std::cout << "Filling true Theta:  " << trueThetaDeg << ", resolutionTheta: " << resolutionTheta << std::endl;
    //std::cout << "run|subrun|gate: " << mc_run << "|" << mc_subrun << "|" << mc_nthEvtInFile+1 << std::endl;
    /*
    if (fabs(recoThetaDeg - trueThetaDeg)>50){
      std::cout << "big theta diff, reco theta: " << recoThetaDeg << ", true theta: " << trueThetaDeg << ", run|subrun|gate: " << mc_run << "|" << mc_subrun << "|" << mc_nthEvtInFile+1 << std::endl;
      }*/

    h1->Fill(trueP, resolutionP);
    h2->Fill(truePt, resolutionPt);

    h3->Fill(trueThetaDeg, resolutionTheta);
  }

  
  TCanvas* c = new TCanvas("c", "", 800, 600);
  std::vector<std::pair<TH1*, std::string>> plots = {
    //{h1, "electron_p_resolution.png"},
    //{h2, "electron_pt_resolution.png"},
    {h3, "electron_theta_resolution.png"},
  };

    
  for (const auto &pair : plots) {
    c->Clear();
    pair.first->Draw("colz");
    c->Update();
    c->SaveAs(pair.second.c_str());
  }
    
  delete file;
}
