#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TMath.h>
#include <Math/Vector3D.h>
#include <Math/RotationX.h>
#include <iostream>
#include <vector>

void PlotProtonResolution() {
  TString filename = "/exp/minerva/data/users/cpernas/NuE_TKI/MC_June_05_2025_full_FHC.root";

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

  // Proton mass (in MeV)
  const double M_p = 938.272013;
  const int MAX_MCPARTS = 200;
  const int MAX_PROTON_NODES = 200;

  
  // Branch variables
  int selectionCategory, mc_run, mc_subrun, mc_nthEvtInFile;
  int mc_nFSPart;
  double mc_FSPartE[MAX_MCPARTS];
  int mc_FSPartPDG[MAX_MCPARTS];
  double mc_FSPartPx[MAX_MCPARTS];
  double mc_FSPartPy[MAX_MCPARTS];
  double mc_FSPartPz[MAX_MCPARTS];
  double recoPMev;
  double recoTheta;
  int n_protonNodes;
  double protonNodeE[MAX_PROTON_NODES];
  
  
  tree->SetBranchAddress("selectionCategory", &selectionCategory);
  tree->SetBranchAddress("mc_run", &mc_run);
  tree->SetBranchAddress("mc_subrun", &mc_subrun);
  tree->SetBranchAddress("mc_nthEvtInFile", &mc_nthEvtInFile);
  tree->SetBranchAddress("mc_nFSPart", &mc_nFSPart);
  tree->SetBranchAddress("mc_FSPartE", &mc_FSPartE);
  tree->SetBranchAddress("mc_FSPartPDG", &mc_FSPartPDG);
  tree->SetBranchAddress("mc_FSPartPx", &mc_FSPartPx);
  tree->SetBranchAddress("mc_FSPartPy", &mc_FSPartPy);
  tree->SetBranchAddress("mc_FSPartPz", &mc_FSPartPz);
  tree->SetBranchAddress("MasterAnaDev_proton_P_fromdEdx", &recoPMev);
  tree->SetBranchAddress("MasterAnaDev_proton_theta", &recoTheta);
  tree->SetBranchAddress("MasterAnaDev_proton_nodes_nodesNormE_sz", &n_protonNodes);
  tree->SetBranchAddress("MasterAnaDev_proton_nodes_nodesNormE", &protonNodeE);

  // Create histograms
  //momentum bins are like 0.45 to 1.2 (my signal requirements)
  //chi2 bins are like 0 to 40 ish
  //node energy will be like 0 to 50/100 MeV...
  //TH2F* h1 = new TH2F("h1", "Proton p Resolution;True Proton P (GeV); (Reco - True)/True", 50, 0.45, 1.2, 50, -0.5, 0.5);
  //TH2F* h2 = new TH2F("h2", "Proton pT Resolution;True Proton pT (GeV); (Reco - True)/True", 50, 0.45, 1.2, 50, -0.5, 0.5);
  //TH2F* h1 = new TH2F("h1", "Proton p Resolution;dE/dX chi^2; (Reco - True)/True", 30, 0, 30, 50, -0.5, 0.5);
  //TH2F* h2 = new TH2F("h2", "Proton pT Resolution;dE/dX chi^2; (Reco - True)/True", 30, 0, 30, 50, -0.5, 0.5);

  TH2F* h1 = new TH2F("Node 0+1", "Node 0+1;Node energy (MeV); (True - Reco)/True", 20, 0, 50, 20, -0.5, 0.5);
  TH2F* h2 = new TH2F("Node 2", "Node 2;Node energy (MeV); (True - Reco)/True", 20, 0, 20, 20, -0.5, 0.5);
  TH2F* h3 = new TH2F("Node 3", "Node 3;Node energy (MeV); (True - Reco)/True", 20, 0, 20, 20, -0.5, 0.5);
  TH2F* h4 = new TH2F("Node 4", "Node 4;Node energy (MeV); (True - Reco)/True", 20, 0, 20, 20, -0.5, 0.5);
  TH2F* h5 = new TH2F("Node 5", "Node 5;Node energy (MeV); (True - Reco)/True", 20, 0, 20, 20, -0.5, 0.5);
  //TH2F* h2 = new TH2F("h2", "Proton pT Resolution;dE/dX chi^2; (Reco - True)/True", 30, 0, 30, 50, -0.5, 0.5);

  //TH2F* h3 = new TH2F("h3", "Proton Theta Resolution;True Proton Theta (deg); (Reco - True)", 50, 0, 70, 50, -30, 30);
  //TH1F* hchi2 = new TH1F("chi2", "; proton momentum resolution, (reco-true)/true; events", 50, -2, 2);

  Long64_t nentries = tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    //if (i>1000) break;
    tree->GetEntry(i);

    if (selectionCategory != -999) continue;
    //std::cout << "\nlooking at event: " << i << "\n" << std::endl;

    // Get highest-energy signal proton
    double highestE = -999;
    int bestIndex = -1;
    
    for (int j = 0; j < mc_nFSPart; ++j) {
      if (mc_FSPartPDG[j] == 2212 && mc_FSPartE[j] > highestE) {
        // Apply momentum and angle cuts
        double px = mc_FSPartPx[j];
        double py = mc_FSPartPy[j];
        double pz = mc_FSPartPz[j];

        double p = std::sqrt(px * px + py * py + pz * pz);
        //double trueTheta = std::acos(pz / p);
	ROOT::Math::XYZVector p_vec(px, py, pz); //our momentum 3 vector
	ROOT::Math::RotationX r(-3.3 * (TMath::Pi() / 180.)); //represents a rotation about the x axis of 3.3 degrees, aka 3.3 degrees down toward floor to match beam direction
	double trueTheta = (r(p_vec)).Theta();
        double trueThetaDeg = trueTheta * 180.0 / TMath::Pi();
        if (p > 450 && p < 1200 && trueThetaDeg < 70) {
          highestE = mc_FSPartE[j];
          bestIndex = j;
        }
      }
    }
    if (bestIndex < 0 || recoPMev < -9990) continue;

    //Calculate proton end dE/dX chi^2 for end of track if > 5 nodes
    //Otherwise apply individual ESC cuts on all 5 node energies, if it passes them all set chi^2 to 0
    //If not set it to some high value (75 or whatever) (Dan's Hybrid ESC / Chi^2 method)
    //Fit means and Sigmas are taken from Dan at the moment...

    //These values (the fit Bragg Peak fit basically) come from docdb 27170 from Dan, should re-evauate though...
    double means[5] = {31.302, 11.418, 9.769, 8.675, 7.949};
    double sigmas[5] = {8.997, 3.075, 2.554, 2.484, 2.232};
    
    double nodeEnergyVal = 0;
    double chi2=0;   
    if(n_protonNodes>5){
      for(int j=0;j<n_protonNodes;j++){
	if(j==6) break;
	if(j==0) nodeEnergyVal+=protonNodeE[0];
	else if(j==1) nodeEnergyVal+=protonNodeE[1];
	else nodeEnergyVal=protonNodeE[j];
	if(j>=1){
	  chi2+=(nodeEnergyVal-means[j-1])*(nodeEnergyVal-means[j-1])/(sigmas[j-1]*sigmas[j-1]);
	}
      }
    }
    else{
      //chi2=-999;
      bool pass = true;
      //Cut Values based on 22302 (carlos - i assume this is a docdb number? -> yes)
      double cutval1 = 19; //node 0 + node 1
      double cutval2 = 10; //node 2
      double cutval3 = 9;  //node 3
      double cutval4 = 8;  //node 4
      double cutval5 = 5;  //node 5
      //for(int j=0;j<6;j++){       
      //std::cout << "test, proton nodesE[" << j << "]: " << protonNodeE[j] << std::endl;
      //std::cout << " " << std::endl;
      //}
      if(n_protonNodes==0) { pass = false; }///no nodes
      if(n_protonNodes>1) {
	if(protonNodeE[0]+protonNodeE[1] < cutval1) { pass = false; }
      }
      if(n_protonNodes>2) {
	if(protonNodeE[2] < cutval2) { pass = false; }
      }
      if(n_protonNodes>3) {
	if(protonNodeE[3] < cutval3) { pass = false; }
      }
      if(n_protonNodes>4) {
	if(protonNodeE[4] < cutval4) { pass = false; }
      }
      if(n_protonNodes>5) {
	if(protonNodeE[5] < cutval5){ pass = false; }
      }
      //survived loops?
      if(pass) {
	chi2=0;
	//std::cout << "event PASSED direct esc cut" << std::endl;
      }
      else {
	chi2=75;
	//std::cout << "event FAILED direct esc cut" << std::endl;
      }      
      
    }
    //std::cout << "chi2 = " << chi2 << std::endl;

    //calc true kinematics
    double trueE = mc_FSPartE[bestIndex];
    double trueP = std::sqrt(trueE * trueE - M_p * M_p) / 1000.0; // in GeV
    double px = mc_FSPartPx[bestIndex];
    double py = mc_FSPartPy[bestIndex];
    double pz = mc_FSPartPz[bestIndex];
    ROOT::Math::XYZVector p_vec(px, py, pz);
    ROOT::Math::RotationX r(-3.3 * (TMath::Pi() / 180.)); //represents a rotation about the x axis of 3.3 degrees, aka 3.3 degrees down toward floor to match beam direction
    double trueTheta = (r(p_vec)).Theta();
    double trueThetaDeg = trueTheta * 180.0 / TMath::Pi();
    double truePt = trueP * std::sin(trueTheta);

    //reco kinematics
    double recoP = recoPMev / 1000.0;
    double recoPt = recoP * std::sin(recoTheta);
    double recoThetaDeg = recoTheta * 180.0 / TMath::Pi();
    
    double resolutionP = (recoP - trueP) / trueP;
    double resolutionPt = (recoPt - truePt) / truePt;
    double resolutionTheta = (recoThetaDeg - trueThetaDeg); //not over theta...? Mike's suggestion
    
    //std::cout << "Filling true P:  " << trueP << ", Resolution: " << resolutionP << std::endl;
    //std::cout << "Filling true Pt:  " << truePt << ", resolutionPt: " << resolutionPt << std::endl;
    //std::cout << "Filling true Theta:  " << trueThetaDeg << ", resolutionTheta: " << resolutionTheta << std::endl;
    //std::cout << "proton track dE/dX chi^2: " << chi2 << std::endl;
    /*
    if (fabs(recoThetaDeg - trueThetaDeg)>50){
      std::cout << "big theta diff, reco theta: " << recoThetaDeg << ", true theta: " << trueThetaDeg << ", run|subrun|gate: " << mc_run << "|" << mc_subrun << "|" << mc_nthEvtInFile+1 << std::endl;
      }*/

    //h1->Fill(trueP, resolutionP);
    //h2->Fill(truePt, resolutionPt);

    if (n_protonNodes<=5){
      h1->Fill(protonNodeE[0]+protonNodeE[1], -resolutionP);
      h2->Fill(protonNodeE[2], -resolutionP);
      h3->Fill(protonNodeE[3], -resolutionP);
      h4->Fill(protonNodeE[4], -resolutionP);
      h5->Fill(protonNodeE[5], -resolutionP);
    }
    //h3->Fill(trueThetaDeg, resolutionTheta);
    //hchi2->Fill(chi2);
  }

  TCanvas* c = new TCanvas("c", "", 800, 600);
  std::vector<std::pair<TH1*, std::string>> plots = {
    //{h1, "proton_p_resolution.png"},
    //{h2, "proton_pt_resolution.png"},
    //{h3, "proton_theta_resolution.png"},
    {h1, "P_resolution_node0+1.png"},
    {h2, "P_resolution_node2.png"},
    {h3, "P_resolution_node3.png"},
    {h4, "P_resolution_node4.png"},
    {h5, "P_resolution_node5.png"},

    //{hchi2, "proton_chi2.png"},
  };

  int cutvals[5] = {19,10,9,8,5};
  int k = 0;
  for (const auto &pair : plots) {
    c->Clear();
    //c->SetLogz();
    pair.first->Draw("colz");
    
    double ymin = pair.first->GetYaxis()->GetXmin();
    double ymax = pair.first->GetYaxis()->GetXmax();
    double x_val = 0.5;
    TLine* line = new TLine(cutvals[k], ymin, cutvals[k], ymax);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->Draw();

    c->Update();
    c->SaveAs(pair.second.c_str());
    k++;
  }
    
  delete file;
}
