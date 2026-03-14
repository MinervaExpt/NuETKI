//Root macro that loops through an input TTree, preferably one that has already selected my events,
// and recalculates my analysis variables in true and reco space.
// I then make and save resolution plots. 
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TMath.h>
#include <Math/Vector3D.h>
#include <Math/RotationX.h>
#include <iostream>
#include <vector>

const double M_n = 939.56536;  //MeV
const double M_p = 938.272013; //MeV
const double M_e = 0.51099895; //MeV
const double M_pi = 139.57039; //MeV
const double M_nucleon = (1.5*M_n+M_p)/2.5;
const double pi = 3.141592653589793;
const double beam_tilt = 3.343 * (pi / 180);
const bool signal_only = false;

int GetHighestEnergySignalProton(const double* FSPartE, const double* FSPartPx, const double* FSPartPy, const double* FSPartPz, const int* FSPartPDG, int nFSPart) {
  // Get highest-energy signal proton
  double highestE = -9999;
  int bestIndex = -9999;
  
  for (int j = 0; j < nFSPart; ++j) {
    if (FSPartPDG[j] == 2212 && FSPartE[j] > highestE) {
      // Apply momentum and angle cuts
      double px = FSPartPx[j];
      double py = FSPartPy[j];
      double pz = FSPartPz[j];
      
      double p = std::sqrt(px * px + py * py + pz * pz);
      //double trueTheta = std::acos(pz / p);
      ROOT::Math::XYZVector p_vec(px, py, pz); //our momentum 3 vector
      ROOT::Math::RotationX r(-3.3 * (TMath::Pi() / 180.)); //represents a rotation about the x axis of 3.3 degrees, aka 3.3 degrees down toward floor to match beam direction
      double trueTheta = (r(p_vec)).Theta();
      double trueThetaDeg = trueTheta * 180.0 / TMath::Pi();
      if (p > 450 && p < 1200 && trueThetaDeg < 70) {
	highestE = FSPartE[j];
	bestIndex = j;
      }
    }
  }
  return bestIndex;
}

double GetEAvailTrue(const double* FSPartE, const int* FSPartPDG, int nFSPart) {
  double T_p = 0; //sum of proton kinetic energies
  double T_pi = 0; //sum of charged pion kinetic energies
  double E_pi0 = 0; //sume of neutral pion total energies
  double E_s = 0; //Sum of (strange baryon energy - proton mass)
  double E_sbar = 0; //Sum of (anti baryon energy + proton mass)
  double E_other = 0; //sum of total energy of any other particles, NOT INCLUDING NEUTRONS
    
  for (int i=0; i<nFSPart; i++){
    if (abs(FSPartPDG[i]) == 11 || abs(FSPartPDG[i]) == 13 || FSPartPDG[i] == 2112 || FSPartPDG[i] > 1000000000){
      continue; //don't want to count leptons, neutrons, or nuclear remnants
    } else if (FSPartPDG[i] == 2212){
      T_p += FSPartE[i] - M_p; //add only kinetic (not total) energy for protons
    } else if (abs(FSPartPDG[i]) == 211){
      T_pi += FSPartE[i] - M_pi; //add only kinetic (not total) energy for charged pions
    } else if (abs(FSPartPDG[i]) == 111){
      E_pi0 += FSPartE[i]; //add pi0 total energy
    } else if (FSPartPDG[i] > 2000){
      E_s += FSPartE[i] - M_p; //add total energy - proton mass for strange baryons (why?)
    } else if (FSPartPDG[i] < 2000){
      E_sbar += FSPartE[i] + M_p; //add total energy + proton mass for strange antibaryons (why???)
    } else {
      E_other += FSPartE[i]; //add total energy for anything else (mostly gammas, kaons i think?)
    }
  }
  
  double E_avail = T_p + T_pi + E_pi0 + E_s + E_sbar + E_other;
  return E_avail;
}

//These all take vectors as arguments so I can use them for both true and reco calculations
double CalcDeltaPt(const ROOT::Math::XYZVector& lepton_vec, const ROOT::Math::XYZVector& proton_vec) {
  ROOT::Math::XYZVector electronPt_vec = lepton_vec;
  ROOT::Math::XYZVector protonPt_vec = proton_vec;
  electronPt_vec.SetZ(0);
  protonPt_vec.SetZ(0);

  ROOT::Math::XYZVector deltaPt_vec = electronPt_vec + protonPt_vec;
  return deltaPt_vec.R();
}

double CalcDeltaPtX(const ROOT::Math::XYZVector& lepton_vec, const ROOT::Math::XYZVector& proton_vec) {
  ROOT::Math::XYZVector electronPt_vec = lepton_vec;
  ROOT::Math::XYZVector protonPt_vec = proton_vec;
  electronPt_vec.SetZ(0);
  protonPt_vec.SetZ(0);
  ROOT::Math::XYZVector deltaPt_vec = electronPt_vec + protonPt_vec;

  double num = deltaPt_vec.Y()*electronPt_vec.X() - deltaPt_vec.X()*electronPt_vec.Y();
  double denom = electronPt_vec.R();
  return num/denom;
}

double CalcDeltaPtY(const ROOT::Math::XYZVector& lepton_vec, const ROOT::Math::XYZVector& proton_vec) {
  ROOT::Math::XYZVector electronPt_vec = lepton_vec;
  ROOT::Math::XYZVector protonPt_vec = proton_vec;
  electronPt_vec.SetZ(0);
  protonPt_vec.SetZ(0);
  ROOT::Math::XYZVector deltaPt_vec = electronPt_vec + protonPt_vec;

  double num = deltaPt_vec.Dot(electronPt_vec);
  double denom = electronPt_vec.R();
  return -num/denom;
}

double CalcDeltaPl(const ROOT::Math::XYZVector& lepton_vec, const ROOT::Math::XYZVector& proton_vec, double E_lep, double E_proton) {
  ROOT::Math::XYZVector electronPt_vec = lepton_vec;
  ROOT::Math::XYZVector protonPt_vec = proton_vec;
  electronPt_vec.SetZ(0);
  protonPt_vec.SetZ(0);
  ROOT::Math::XYZVector deltaPt_vec = electronPt_vec + protonPt_vec;
  
  double mass_nuke = 11.188*1000;
  double binding_energy = 27.13;
  
  double R = mass_nuke + lepton_vec.Z() + proton_vec.Z() - ( E_lep + E_proton );
  double numerator = pow(mass_nuke - M_p + binding_energy, 2) + pow(deltaPt_vec.R(), 2);
  
  double DeltaPl = 0.5 * R - numerator/(2*R);
  return DeltaPl;

}

double CalcAlphaPt(const ROOT::Math::XYZVector& lepton_vec, const ROOT::Math::XYZVector& proton_vec) {
  ROOT::Math::XYZVector electronPt_vec = lepton_vec;
  ROOT::Math::XYZVector protonPt_vec = proton_vec;
  electronPt_vec.SetZ(0);
  protonPt_vec.SetZ(0);

  ROOT::Math::XYZVector deltaPt_vec = electronPt_vec + protonPt_vec;

  double numerator = ((-1*electronPt_vec).Dot(deltaPt_vec));
  double denominator = ( electronPt_vec.R() * deltaPt_vec.R() );
  double alpha = std::acos(numerator/denominator);
  
  return alpha * 180/pi;
}

double CalcPhiPt(const ROOT::Math::XYZVector& lepton_vec, const ROOT::Math::XYZVector& proton_vec) {
  ROOT::Math::XYZVector electronPt_vec = lepton_vec;
  ROOT::Math::XYZVector protonPt_vec = proton_vec;
  electronPt_vec.SetZ(0);
  protonPt_vec.SetZ(0);

  double numerator = ((-1*electronPt_vec).Dot(protonPt_vec));
  double denominator = ( electronPt_vec.R() * protonPt_vec.R() );
  double phi = std::acos(numerator/denominator);
  
  return phi * 180/pi;
}

void PlotTKIResolution() {
 TString filename = "/exp/minerva/data/users/cpernas/NuE_TKI/selected_tuples/TUPLE_MC_June_12_2025_NuE_Only.root";

  TFile* file = TFile::Open(filename);
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file\n";
    return;
  }

  TTree* tree = (TTree*)file->Get("Signal_Region"); 
  if (!tree) {
    std::cerr << "Error: TTree not found\n";
    return;
  }

  std::vector<std::string> variables = {"E_lep", "E_avail", "E_nu", "lepton_pT", "lepton_pL", "theta_lep", "proton_P", "proton_pT", "proton_theta", "delta_Pt", "delta_PtX", "delta_PtY", "delta_Pl", "P_n", "alpha_Pt", "acoplanarity"};
  std::vector<std::pair<double, double>> axes_limits = {
    {0, 20}, //E_lep
    {0, 1.2}, //E_avail
    {0, 20}, //E_nu
    {0, 1.6}, //lepton_pT
    {0, 20}, //lepton_pL
    {0, 40}, //lepton_theta
    {0, 1.4}, //proton_P
    {0, 1.1}, //proton_pT
    {0, 180}, //proton_theta
    {0, 1.4}, //delta_Pt
    {-2, 2}, //delta_PtX
    {-2, 2}, //delta_PtY
    {-0.5, 1}, //delta_Pl
    {0, 1}, //P_n
    {0, 180}, //alpha_Pt
    {0, 180}, //acoplanarity
  };

  // Create histograms
  std::vector<TH1D*> resolution_hists_1D;
  std::vector<TH2D*> resolution_hists;
  for (int i = 0; i < variables.size(); i++) {
    std::string hist_name_1D = "h_" + variables[i] + "_resolution_1D";
    std::string hist_label_1D = "; " + variables[i] + " resolution; events";
    resolution_hists_1D.push_back(new TH1D(hist_name_1D.c_str(), hist_label_1D.c_str(), 50, -1, 1));
    std::string hist_name = "h_" + variables[i] + "_resolution";
    std::string hist_label = variables[i] + " resolution; true " + variables[i] + "; (reco - true)/true";
    resolution_hists.push_back(new TH2D(hist_name.c_str(), hist_label.c_str(), 50, axes_limits[i].first, axes_limits[i].second, 50, -1, 1));
  }

  std::vector<double> avg_resolutions;
  for (const auto& var : variables) {
    avg_resolutions.push_back(0);
  }

  const int MAX_PRONGS = 10;
  const int MAX_MCPARTS = 200;
  const int MAX_PROTON_NODES = 200;
  
  // Branch variables
  int selectionCategory, mc_run, mc_subrun, mc_nthEvtInFile, mc_nFSPart, n_protonNodes;
  double recoil_E_ecal, recoil_E_tracker, proton_E, proton_P, proton_Px, proton_Py, proton_Pz, proton_theta, mc_incomingE;
  double vtx[4], mc_vtx[4], trueElectron4Vec[4];
  double prong_ECALCalibE[MAX_PRONGS];
  int mc_FSPartPDG[MAX_MCPARTS];
  double mc_FSPartE[MAX_MCPARTS], mc_FSPartPx[MAX_MCPARTS], mc_FSPartPy[MAX_MCPARTS], mc_FSPartPz[MAX_MCPARTS];
  double protonNodeE[MAX_PROTON_NODES];
  std::vector<std::vector<double>> *prong_part_E = nullptr;  

  //-999 = signal, 0 = nonQELike NuE, 1 = otherNuE, 2 = NCPi0, 3 = NuMuCCPi0, -1 = other (check these...)
  tree->SetBranchAddress("selectionCategory", &selectionCategory);
  //reco
  tree->SetBranchAddress("vtx", &vtx);
  tree->SetBranchAddress("MasterAnaDev_proton_E_fromdEdx", &proton_E);
  tree->SetBranchAddress("MasterAnaDev_proton_P_fromdEdx", &proton_P);
  tree->SetBranchAddress("MasterAnaDev_proton_Px_fromdEdx", &proton_Px);
  tree->SetBranchAddress("MasterAnaDev_proton_Py_fromdEdx", &proton_Py);
  tree->SetBranchAddress("MasterAnaDev_proton_Pz_fromdEdx", &proton_Pz);
  tree->SetBranchAddress("MasterAnaDev_proton_theta", &proton_theta);
  tree->SetBranchAddress("MasterAnaDev_proton_nodes_nodesNormE_sz", &n_protonNodes);
  tree->SetBranchAddress("MasterAnaDev_proton_nodes_nodesNormE", &protonNodeE);
  tree->SetBranchAddress("blob_recoil_E_tracker", &recoil_E_tracker);
  tree->SetBranchAddress("blob_recoil_E_ecal", &recoil_E_ecal);
  tree->SetBranchAddress("prong_ECALCalibE", &prong_ECALCalibE);
  tree->SetBranchAddress("prong_part_E", &prong_part_E);
  //truth
  tree->SetBranchAddress("mc_vtx", &mc_vtx);
  tree->SetBranchAddress("mc_run", &mc_run);
  tree->SetBranchAddress("mc_subrun", &mc_subrun);
  tree->SetBranchAddress("mc_nthEvtInFile", &mc_nthEvtInFile);
  tree->SetBranchAddress("mc_incomingE", &mc_incomingE);
  tree->SetBranchAddress("mc_nFSPart", &mc_nFSPart);
  tree->SetBranchAddress("mc_FSPartE", &mc_FSPartE);
  tree->SetBranchAddress("mc_FSPartPDG", &mc_FSPartPDG);
  tree->SetBranchAddress("mc_FSPartPx", &mc_FSPartPx);
  tree->SetBranchAddress("mc_FSPartPy", &mc_FSPartPy);
  tree->SetBranchAddress("mc_FSPartPz", &mc_FSPartPz);
  tree->SetBranchAddress("mc_primFSLepton", &trueElectron4Vec);

  int eventCount = 0;
  Long64_t nentries = tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    //if(i%1000==0) std::cout << i << " / " << nentries << "\r" << std::flush;
    //if (i>1000) break;
    tree->GetEntry(i);
    int protonIndex = GetHighestEnergySignalProton(mc_FSPartE, mc_FSPartPx, mc_FSPartPy, mc_FSPartPz, mc_FSPartPDG, mc_nFSPart);
    if (protonIndex < 0) continue;
    if (selectionCategory != -999 && signal_only) continue;
    //std::cout << "\nlooking at event: " << i << "\n" << std::endl;    

    
    ROOT::Math::RotationX r(-beam_tilt); 
    
    //==============================
    // True kinematics
    //==============================
    double true_E_lep = trueElectron4Vec[3]/1000;
    double true_E_avail = GetEAvailTrue(mc_FSPartE, mc_FSPartPDG, mc_nFSPart)/1000;
    double true_E_nu = true_E_lep + true_E_avail; //this is what I do in my analysis, rather than mc_incomingE
    
    ROOT::Math::XYZVector true_lep_P3D(trueElectron4Vec[0], trueElectron4Vec[1], trueElectron4Vec[2]);
    ROOT::Math::XYZVector true_lep_P3D_beamframe = r(true_lep_P3D);
    double true_lepton_pT = true_lep_P3D_beamframe.Rho()/1000;
    double true_lepton_pL = true_lep_P3D_beamframe.Z()/1000;
    double true_theta_lep = true_lep_P3D_beamframe.Theta()* (180/pi);

    ROOT::Math::XYZVector true_proton_P3D(mc_FSPartPx[protonIndex], mc_FSPartPy[protonIndex], mc_FSPartPz[protonIndex]); 
    ROOT::Math::XYZVector true_proton_P3D_beamframe = r(true_proton_P3D);
    double true_proton_P = true_proton_P3D_beamframe.R()/1000;
    double true_proton_pT = true_proton_P3D_beamframe.Rho()/1000;
    double true_proton_theta = true_proton_P3D_beamframe.Theta()*(180/pi);

    double true_delta_Pt = CalcDeltaPt(true_lep_P3D_beamframe, true_proton_P3D_beamframe)/1000;
    double true_delta_PtX = CalcDeltaPtX(true_lep_P3D_beamframe, true_proton_P3D_beamframe)/1000;
    double true_delta_PtY = CalcDeltaPtY(true_lep_P3D_beamframe, true_proton_P3D_beamframe)/1000;
    double true_delta_Pl = CalcDeltaPl(true_lep_P3D_beamframe, true_proton_P3D_beamframe, trueElectron4Vec[3], mc_FSPartE[protonIndex])/1000;

    double true_P_n = sqrt(pow(true_delta_Pt, 2) + pow(true_delta_Pl, 2));
    double true_alpha_Pt = CalcAlphaPt(true_lep_P3D_beamframe, true_proton_P3D_beamframe); //already returned in deg
    double true_acoplanarity = CalcPhiPt(true_lep_P3D_beamframe, true_proton_P3D_beamframe); //already returned in deg
    
    //==============================
    // Reco kinematics
    //==============================
    double reco_E_lep = ((*prong_part_E)[0][3] + (prong_ECALCalibE[0]*-0.058))/1000;
    double reco_E_avail = (((recoil_E_tracker + recoil_E_ecal)*1.17) - (((*prong_part_E)[0][3] * 0.008) - 5 ))/1000;
    double reco_E_nu = reco_E_lep + reco_E_avail;
    
    //scale the reco lepton 3 momentum by the amount we're shifting the electron energy
    double scale = 1.0 + ((prong_ECALCalibE[0]*-0.058) /(*prong_part_E)[0][3]);
    ROOT::Math::XYZVector reco_lep_P3D((*prong_part_E)[0][0], (*prong_part_E)[0][1], (*prong_part_E)[0][2]);
    reco_lep_P3D *= scale;
    ROOT::Math::XYZVector reco_lep_P3D_beamframe = r(reco_lep_P3D);
    double reco_lepton_pT = reco_lep_P3D_beamframe.Rho()/1000;
    double reco_lepton_pL = reco_lep_P3D_beamframe.Z()/1000;
    double reco_theta_lep = reco_lep_P3D_beamframe.Theta()*(180/pi);

    ROOT::Math::XYZVector reco_proton_P3D(proton_Px, proton_Py, proton_Pz);
    ROOT::Math::XYZVector reco_proton_P3D_beamframe = r(reco_proton_P3D);
    double reco_proton_P = reco_proton_P3D_beamframe.R()/1000;
    double reco_proton_pT = reco_proton_P3D_beamframe.Rho()/1000;
    double reco_proton_theta = reco_proton_P3D_beamframe.Theta()*(180/pi);

    double reco_delta_Pt = CalcDeltaPt(reco_lep_P3D_beamframe, reco_proton_P3D_beamframe)/1000;
    double reco_delta_PtX = CalcDeltaPtX(reco_lep_P3D_beamframe, reco_proton_P3D_beamframe)/1000;
    double reco_delta_PtY = CalcDeltaPtY(reco_lep_P3D_beamframe, reco_proton_P3D_beamframe)/1000;
    double reco_delta_Pl = CalcDeltaPl(reco_lep_P3D_beamframe, reco_proton_P3D_beamframe, reco_E_lep*1000, proton_E)/1000;

    double reco_P_n = sqrt(pow(reco_delta_Pt, 2) + pow(reco_delta_Pl, 2));
    double reco_alpha_Pt = CalcAlphaPt(reco_lep_P3D_beamframe, reco_proton_P3D_beamframe); //already returned in deg
    double reco_acoplanarity = CalcPhiPt(reco_lep_P3D_beamframe, reco_proton_P3D_beamframe); //already returned in deg

    //==============================
    // Resolutions
    //==============================
    double E_lep_resolution = (reco_E_lep - true_E_lep)/true_E_lep;
    double E_avail_resolution = (reco_E_avail - true_E_avail)/true_E_avail;
    double E_nu_resolution = (reco_E_nu - true_E_nu)/true_E_nu;
    double lepton_pT_resolution = (reco_lepton_pT - true_lepton_pT)/true_lepton_pT;
    double lepton_pL_resolution = (reco_lepton_pL - true_lepton_pL)/true_lepton_pL;
    double theta_lep_resolution = (reco_theta_lep - true_theta_lep)/true_theta_lep;
    double proton_P_resolution = (reco_proton_P - true_proton_P)/true_proton_P;
    double proton_pT_resolution = (reco_proton_pT - true_proton_pT)/true_proton_pT;

    double proton_theta_resolution = (reco_proton_theta - true_proton_theta)/true_proton_theta;
    double delta_Pt_resolution = (reco_delta_Pt - true_delta_Pt)/true_delta_Pt;
    double delta_PtX_resolution = (reco_delta_PtX - true_delta_PtX)/true_delta_PtX;
    double delta_PtY_resolution = (reco_delta_PtY - true_delta_PtY)/true_delta_PtY;    
    double delta_Pl_resolution = (reco_delta_Pl - true_delta_Pl)/true_delta_Pl;
    double P_n_resolution = (reco_P_n - true_P_n)/true_P_n;
    double alpha_Pt_resolution = (reco_alpha_Pt - true_alpha_Pt)/true_alpha_Pt;
    double acoplanarity_resolution = (reco_acoplanarity - true_acoplanarity)/true_acoplanarity;

    /*
    std::cout << "================ EVENT " << i << " ================" << std::endl;
    std::cout << "proton_theta: reco = " << reco_proton_theta << ", true = " << true_proton_theta << ", resolution = " << proton_theta_resolution << std::endl;    
    std::cout << "delta_Pt: reco = " << reco_delta_Pt << ", true = " << true_delta_Pt << ", resolution = " << delta_Pt_resolution << std::endl;
    std::cout << "delta_PtX: reco = " << reco_delta_Pt << ", true = " << true_delta_Pt << ", resolution = " << delta_Pt_resolution << std::endl;
    std::cout << "delta_PtY: reco = " << reco_delta_PtY << ", true = " << true_delta_PtY << ", resolution = " << delta_PtY_resolution << std::endl;
    std::cout << "delta_Pl: reco = " << reco_delta_Pl << ", true = " << true_delta_Pl << ", resolution = " << delta_Pl_resolution << std::endl;
    std::cout << "P_n: reco = " << reco_P_n << ", true = " << true_P_n << ", resolution = " << P_n_resolution << std::endl;
    std::cout << "alpha_Pt: reco = " << reco_alpha_Pt << ", true = " << true_alpha_Pt << ", resolution = " << alpha_Pt_resolution << std::endl;
    std::cout << "acoplanarity: reco = " << reco_acoplanarity << ", true = " << true_acoplanarity << ", resolution = " << acoplanarity_resolution << std::endl;
    */
    
    std::vector<double> resolution_values = {
      E_lep_resolution,
      E_avail_resolution,
      E_nu_resolution,
      lepton_pT_resolution,
      lepton_pL_resolution,
      theta_lep_resolution,
      proton_P_resolution,
      proton_pT_resolution,
      proton_theta_resolution,
      delta_Pt_resolution,
      delta_PtX_resolution,
      delta_PtY_resolution,
      delta_Pl_resolution,
      P_n_resolution,
      alpha_Pt_resolution,
      acoplanarity_resolution
    };

    std::vector<double> true_values = {
      true_E_lep,
      true_E_avail,
      true_E_nu,
      true_lepton_pT,
      true_lepton_pL,
      true_theta_lep,
      true_proton_P,
      true_proton_pT,
      true_proton_theta,
      true_delta_Pt,
      true_delta_PtX,
      true_delta_PtY,
      true_delta_Pl,
      true_P_n,
      true_alpha_Pt,
      true_acoplanarity
    };
    
    for (int i = 0; i < resolution_hists.size(); i++) {
      resolution_hists_1D[i]->Fill(resolution_values[i]);
      resolution_hists[i]->Fill(true_values[i], resolution_values[i]);
      avg_resolutions[i] += abs(resolution_values[i]);
    }    
    eventCount++;
  }

  std::cout << " ----------------- FINISHED EVENT LOOP ----------------- " << std::endl;
  for (int i = 0; i < variables.size(); i++) {
    std::cout << "Avg " << variables[i] << " resolution = " << avg_resolutions[i]/eventCount << std::endl;
  }

  TCanvas* c = new TCanvas("c", "", 800, 600);
  for (int i = 0; i < variables.size(); i++) {
    c->Clear();
    //c->SetLogz();
    resolution_hists_1D[i]->Draw("hist");
    c->Update();
    std::string outName_1D = variables[i] + "_resolution_1D.png";
    c->SaveAs(outName_1D.c_str());    

    c->Clear();
    //c->SetLogz();
    resolution_hists[i]->Draw("colz");
    c->Update();
    std::string outName = variables[i] + "_resolution.png";
    c->SaveAs(outName.c_str());    
  }

  /*
  auto leg = new TLegend(0.65,0.70,0.88,0.88);
  leg->AddEntry(h_resolution1, "No ESC/Chi2 cut", "1");
  leg->AddEntry(h_resolution2, "Looser ESC/Chi2 cut", "1");
  leg->AddEntry(h_resolution3, "Original ESC/Chi2 cut", "1");
  leg->Draw();
  c->SetLogy(1);
  c->Update();
  */
  
  //delete file;
  }
