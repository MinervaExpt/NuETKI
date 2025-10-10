//A script that takes an input MAD tuple, and plots a dataMC stack the same way MnvPlotter does
//Input file is just hard coded in, I need to add more variables anyways.
//
//As for output, uhh it kinda depends? Could be plots, could be print output, not sure yet... probably will be an evolving thing
//Just implement both and control with a boolean?

//includes
#include "Math/Vector3D.h"

//some constants
constexpr double M_p = 938.272013; //MeV
constexpr double M_e = 0.51099895; //MeV
constexpr double beam_tilt = -3.3; //degrees
constexpr double pi = 3.141592653589793;

//some hardcode config
//constexpr bool printOutput = true;
constexpr bool printOutput = false;
constexpr bool makePlots = true;
//constexpr bool makePlots = false;

//variable calculators
double CalcDSCALVisE(double HCalVisE, double ECalVisE) {
  if (HCalVisE<=0 && ECalVisE<=0) { return 0; }
  else { return HCalVisE / ECalVisE; }
}

double CalcODCALVisE(double ODVisE, double SideECalVisE) {
  if (ODVisE<=0 && SideECalVisE<=0) { return 0; }
  else { return ODVisE / SideECalVisE; }
}

double CalcEAvail(double E_recoil_tracker, double E_recoil_ecal, double E_electron) { //in GeV
  double E_recoil = E_recoil_tracker + E_recoil_ecal;
  double E_avail = (E_recoil * 1.17 - ((0.008 * E_electron) + 5))/1000.; 
  return E_avail;
}

double CalcModEAvail(double E_avail, double E_prim_proton, const double* E_sec_proton, int nSecProtons) { //in GeV
  if (E_prim_proton < -1) { E_prim_proton= 0; } //Making sure the event has a reco proton candidate

  double sum_secondary_proton_energy = 0;
  for (int i = 0; i < nSecProtons; i++){
    sum_secondary_proton_energy += E_sec_proton[i]; 
  }
  
  double ModEavail = E_avail*1000. - E_prim_proton - sum_secondary_proton_energy;
  return ModEavail/1000.;
}

double CalcESCChi2(int n_nodes, const double* proton_nodes_E) {
    double nodeEnergyVal = 0;
    double chi2=0;

    if(n_nodes>5){ //use chi^2 method, not ESC node method
      double means[5] = {31.302, 11.418, 9.769, 8.675, 7.949};
      double sigmas[5] = {8.997, 3.075, 2.554, 2.484, 2.232};
      for(int i=0;i<n_nodes;i++){
        if(i==6) break;
        if(i==0) nodeEnergyVal+=proton_nodes_E[0];
        else if(i==1) nodeEnergyVal+=proton_nodes_E[1];
        else nodeEnergyVal=proton_nodes_E[i];
        if(i>=1){
          chi2+=(nodeEnergyVal-means[i-1])*(nodeEnergyVal-means[i-1])/(sigmas[i-1]*sigmas[i-1]);
        }
      }
    }
    
    else{ //Has 5 or fewer nodes, so use ESC node by node method
      bool pass = true;
      double cutval1 = 19;
      double cutval2 = 10;
      double cutval3 = 9;
      double cutval4 = 8;
      double cutval5 = 5; 
      if(n_nodes==0) { pass = false; }
      if(n_nodes>1) {
        if(proton_nodes_E[0]+proton_nodes_E[1] < cutval1) { pass = false; }
      }
      if(n_nodes>2) {
        if(proton_nodes_E[2] < cutval2) { pass = false; }
      }
      if(n_nodes>3) {
        if(proton_nodes_E[3] < cutval3) { pass = false; }
      }
      if(n_nodes>4) {
        if(proton_nodes_E[4] < cutval4) { pass = false; }
      }
      if(n_nodes>5) {
        if(proton_nodes_E[5] < cutval5){ pass = false; }
      }
      if(pass) chi2=0;
      else chi2=75; 
    }
    return chi2;
}

//wrt beam axis, takes a 3 vector, doesn't care about particle species
double CalcTheta(double px, double py, double pz) {
  ROOT::Math::XYZVector p(px, py, pz);
  ROOT::Math::RotationX r(-3.3 * (pi / 180.)); //represents rotation of 3.3 deg downward to match beam angle
  double theta_1 = (r(p)).Theta();

  //testing manual calculation just to convince myself its the same
  double p_vec[3] = {px, py, pz};
  double beam_tilt_unit_vec[3] = {0, sin(-3.3 * (pi / 180.)), cos(-3.3 * (pi / 180.))};
  double num = p_vec[0]*beam_tilt_unit_vec[0] + p_vec[1]*beam_tilt_unit_vec[1] + p_vec[2]*beam_tilt_unit_vec[2];
  double denom = sqrt(pow(p_vec[0], 2) + pow(p_vec[1], 2) + pow(p_vec[2], 2));
  double theta_2 = acos(num / denom);
  
  return theta_1;
}


//returns a root XYZVector of the transverse momentum vector only (completely orthogonal to beam dir)
//note this has nonzero z component because neutrino direction is 3.3 degrees offset from the z direction
ROOT::Math::XYZVector CalcPtVec(double px, double py, double pz) {
  ROOT::Math::XYZVector p_vec(px, py, pz);

  ROOT::Math::RotationX r(beam_tilt * (pi/180.)); //rotation downwards, into beam direction
  ROOT::Math::RotationX r2(-1*beam_tilt * (pi/180.)); //rotation back up, into lab frame

  ROOT::Math::XYZVector pt_vec = r(p_vec);
  //std::cout << "3 MOMENTUM IN BEAM COORDINATES: ( " << pt_vec.X() << ", " << pt_vec.Y() << ", " << pt_vec.Z() << " ) " << std::endl;
  pt_vec.SetZ(0);
  pt_vec = r2(pt_vec);

  return pt_vec;
}

//Defined as the component of delta Pt orthogonal to electron Pt
//AKA the scalar rejection, or the magnitude of the orthogonal vector
double CalcDeltaPtx(ROOT::Math::XYZVector& electron_pt, ROOT::Math::XYZVector& delta_pt) {
  //First rotate both vectors back into the transverse plane
  ROOT::Math::RotationX r(beam_tilt * (pi/180.)); //rotation downwards, into beam direction

  ROOT::Math::XYZVector electron_pt_in_beam_frame = r(electron_pt);
  ROOT::Math::XYZVector delta_pt_in_beam_frame = r(delta_pt);

  //expression for the scalar rejection, aka perpendicular dot product
  double num = delta_pt_in_beam_frame.Y()*electron_pt_in_beam_frame.X() - delta_pt_in_beam_frame.X()*electron_pt_in_beam_frame.Y();
  double denom = sqrt(electron_pt_in_beam_frame.Mag2());
  double delta_ptx = num / denom;

  return delta_ptx;
}

//Defined as the component of delta Pt parallel to electron Pt
//AKA the scalar projection, or the magnitude of the projected vector
double CalcDeltaPty(ROOT::Math::XYZVector& electron_pt, ROOT::Math::XYZVector& delta_pt) {
  //First rotate both vectors back into the transverse plane
  ROOT::Math::RotationX r(beam_tilt * (pi/180.)); //rotation downwards, into beam direction

  ROOT::Math::XYZVector electron_pt_in_beam_frame = r(electron_pt);
  ROOT::Math::XYZVector delta_pt_in_beam_frame = r(delta_pt);

  //expression for the scalar projection
  double num = delta_pt_in_beam_frame.Dot(electron_pt_in_beam_frame);
  //double num = delta_pt_in_beam_frame.X()*electron_pt_in_beam_frame.X() + delta_pt_in_beam_frame.Y()*electron_pt_in_beam_frame.Y();
  double denom = sqrt(electron_pt_in_beam_frame.Mag2());
  double delta_pty = num / denom;

  return delta_pty;
}

//Angle between the (negative) lepton_pt vec and delta pt vec
double CalcAlphaPt(ROOT::Math::XYZVector& electron_pt, ROOT::Math::XYZVector& delta_pt) {
  double num = ((-1*electron_pt).Dot(delta_pt));
  double denom = ( sqrt(electron_pt.Mag2()) * sqrt(delta_pt.Mag2()) );
  double alpha = std::acos(num/denom);
  
  //std::cout << "boosting angle alpha: " << alpha * 180/pi << std::endl;
  return alpha * 180/pi;  //return alpha in degrees cause that's how I've set up my bins for now
}

//Angle between the (negative) lepton_pt vec and the proton pt vec
double CalcPhiPt(ROOT::Math::XYZVector& electron_pt, ROOT::Math::XYZVector& proton_pt) {
  double num = ((-1*electron_pt).Dot(proton_pt));
  double denom = ( sqrt(electron_pt.Mag2()) * sqrt(proton_pt.Mag2()) );
  double phi = std::acos(num/denom);
  
  //std::cout << "boosting angle alpha: " << alpha * 180/pi << std::endl;
  return phi * 180/pi;  //return alpha in degrees cause that's how I've set up my bins for now
}

void ReadSelectedMADTuple() {
  //TString filename = "/exp/minerva/data/users/cpernas/TUPLE_MC_May_01_2025_me1B_all_cuts_with_Michel_sideband.root";
  //TString filename = "/exp/minerva/data/users/cpernas/TUPLE_MC_May_14_2025_me1M_all_cuts_with_Michel_sideband.root";
  //TString filename = "/exp/minerva/data/users/cpernas/TUPLE_MC_May_12_2025_NuEOnly_me1M_all_cuts_with_Michel_sideband.root";
  //TString filename = "/exp/minerva/data/users/cpernas/TUPLE_MC_June_05_2025_me1B.root";
  TString filename = "/exp/minerva/app/users/cpernas/TUPLE_MC_June_12_2025_VertexTrackMultiplicity.root";

  TString treename = "Signal_Region";

  int nEventsWithNoIsoBlobs = 0;

  
  TFile *file = TFile::Open(filename);
  if (!file || file->IsZombie()) {
    std::cerr << "Error: cannot open file " << filename << std::endl;
    return;
  }

  TTree *tree = (TTree*)file->Get(treename);
  if (!tree) {
    std::cerr << "Error: cannot find tree " << treename << " in file." << std::endl;
    return;
  }

  const int kMaxProtonNodes = 200;
  const int kMaxProtons = 100; // must be >= max size ever used
  const int maxProngs = 10; // must be >= max size ever used

  //ints
  int hasNoVertexMismatch, startPointVertexMultiplicity, vertexTrackMultiplicity, hasNoBackExitingTracks, nDeadDiscrPairsUpstream, nMichels, nProtonTrackNodes, nSecProtons, nIsoBlobs, selectionCategory, mc_run, mc_subrun, mc_nthEvtInFile, n_prongs;
  double blob_recoil_E_tracker, blob_recoil_E_ecal, Psi, protonTheta, protonEnergy, protonP, protonPx, protonPy, protonPz, nonvtx_iso_blobs_energy;

  double prong_firstFireFraction[maxProngs], prong_ECALVisE[maxProngs], prong_HCALVisE[maxProngs], prong_ODVisE[maxProngs], prong_SideECALVisE[maxProngs], prong_TransverseGapScore[maxProngs], prong_NonMIPClusFrac[maxProngs], prong_dEdXMeanFrontTracker[maxProngs], prong_part_score[maxProngs];
  double nonvtx_iso_blobs_distance_in_prong[maxProngs], nonvtx_iso_blobs_energy_in_prong[maxProngs], nonvtx_iso_blobs_start_position_z_in_prong[maxProngs];
  double proton_nodes_E[kMaxProtonNodes];
  double sec_protons_T_fromdEdx[kMaxProtons];
  double vtx[4];
  double electron_4p_vec[4];

  //truth branches
  double truth_proton_px[20], truth_proton_py[20], truth_proton_pz[20], truth_proton_E[20], truth_proton_theta_wrtbeam[20];
  
  //Reco Branches
  tree->SetBranchAddress("HasNoVertexMismatch", &hasNoVertexMismatch);
  tree->SetBranchAddress("StartPointVertexMultiplicity", &startPointVertexMultiplicity);
  tree->SetBranchAddress("VertexTrackMultiplicity", &vertexTrackMultiplicity);
  tree->SetBranchAddress("HasNoBackExitingTracks", &hasNoBackExitingTracks);
  tree->SetBranchAddress("phys_n_dead_discr_pair_upstream_prim_track_proj", &nDeadDiscrPairsUpstream);
  tree->SetBranchAddress("prong_FirstFireFraction", &prong_firstFireFraction);
  tree->SetBranchAddress("prong_ECALVisE", &prong_ECALVisE);
  tree->SetBranchAddress("prong_HCALVisE", &prong_HCALVisE);
  tree->SetBranchAddress("prong_ODVisE", &prong_ODVisE);
  tree->SetBranchAddress("prong_SideECALVisE", &prong_SideECALVisE);
  tree->SetBranchAddress("prong_TransverseGapScore", &prong_TransverseGapScore);
  tree->SetBranchAddress("prong_NonMIPClusFrac", &prong_NonMIPClusFrac);
  tree->SetBranchAddress("prong_dEdXMeanFrontTracker", &prong_dEdXMeanFrontTracker);
  tree->SetBranchAddress("prong_part_score", &prong_part_score);
  tree->SetBranchAddress("blob_recoil_E_tracker", &blob_recoil_E_tracker);
  tree->SetBranchAddress("blob_recoil_E_ecal", &blob_recoil_E_ecal);
  tree->SetBranchAddress("improved_nmichel", &nMichels);
  tree->SetBranchAddress("MasterAnaDev_proton_nodes_nodesNormE_sz", &nProtonTrackNodes);
  tree->SetBranchAddress("MasterAnaDev_proton_nodes_nodesNormE", &proton_nodes_E);
  tree->SetBranchAddress("MasterAnaDev_sec_protons_T_fromCalo_sz", &nSecProtons);
  tree->SetBranchAddress("MasterAnaDev_proton_theta", &protonTheta);
  tree->SetBranchAddress("MasterAnaDev_proton_calib_energy", &protonEnergy);
  tree->SetBranchAddress("MasterAnaDev_proton_P_fromdEdx", &protonP);
  tree->SetBranchAddress("MasterAnaDev_proton_Px_fromdEdx", &protonPx);
  tree->SetBranchAddress("MasterAnaDev_proton_Py_fromdEdx", &protonPy);
  tree->SetBranchAddress("MasterAnaDev_proton_Pz_fromdEdx", &protonPz);
  tree->SetBranchAddress("MasterAnaDev_sec_protons_T_fromdEdx", &sec_protons_T_fromdEdx);
  //check if this is always the same as &nonvtx_iso_blobs_energy and &nonvtx_iso_blobs_energy_all
  tree->SetBranchAddress("nonvtx_iso_blobs_energy", &nonvtx_iso_blobs_energy);
  tree->SetBranchAddress("nonvtx_iso_blobs_energy_in_prong_sz", &nIsoBlobs);
  tree->SetBranchAddress("nonvtx_iso_blobs_energy_in_prong", &nonvtx_iso_blobs_energy_in_prong);
  tree->SetBranchAddress("nonvtx_iso_blobs_distance_in_prong", &nonvtx_iso_blobs_distance_in_prong);
  tree->SetBranchAddress("nonvtx_iso_blobs_start_position_z_in_prong", &nonvtx_iso_blobs_start_position_z_in_prong);
  tree->SetBranchAddress("n_prongs", &n_prongs);
  
  //tree->SetBranchAddress("prong_part_E", &prong_part_E);
  tree->SetBranchAddress("MasterAnaDev_leptonE", &electron_4p_vec);
  tree->SetBranchAddress("vtx", &vtx);
  
  //Truth Branches
  tree->SetBranchAddress("selectionCategory", &selectionCategory);
  tree->SetBranchAddress("mc_run", &mc_run);
  tree->SetBranchAddress("mc_subrun", &mc_subrun);
  tree->SetBranchAddress("mc_nthEvtInFile", &mc_nthEvtInFile);
  tree->SetBranchAddress("truth_proton_px", &truth_proton_px);
  tree->SetBranchAddress("truth_proton_py", &truth_proton_py);
  tree->SetBranchAddress("truth_proton_pz", &truth_proton_pz);
  tree->SetBranchAddress("truth_proton_E", &truth_proton_E);
  tree->SetBranchAddress("truth_proton_theta_wrtbeam", &truth_proton_theta_wrtbeam);

  TString xaxis_label = "Blob Start Z - vtx Z (mm)";
  //double xbins[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2};
  //int nbins = sizeof(xbins)/sizeof(xbins[0])-1;

  double xbinslow = -2000;
  double xbinshigh= 2000;
  int nbins = 40;
  
  //Create individual histos for stack
  //TH1F* h_sig  = new TH1F("h_sig",  "Signal;" + xaxis_label + ";Events", nbins, xbins);
  //TH1F* h_bkg1 = new TH1F("h_bkg1", "NuECC nonQELike;" + xaxis_label + ";Events", nbins, xbins);
  //TH1F* h_bkg2 = new TH1F("h_bkg2", "other NuECC;" + xaxis_label + ";Events", nbins, xbins);
  //TH1F* h_bkg3 = new TH1F("h_bkg3", "NC Pi0;" + xaxis_label + ";Events", nbins, xbins);
  //TH1F* h_bkg4 = new TH1F("h_bkg4", "NuMu CC Pi0;" + xaxis_label + ";Events", nbins, xbins);
  //TH1F* h_bkg5 = new TH1F("h_bkg5", "Other;" + xaxis_label + ";Events", nbins, xbins);

  TH1F* h_sig  = new TH1F("h_sig",  "Signal;" + xaxis_label + ";Events", nbins, xbinslow, xbinshigh);
  TH1F* h_bkg1 = new TH1F("h_bkg1", "NuECC nonQELike;" + xaxis_label + ";Events", nbins, xbinslow, xbinshigh);
  TH1F* h_bkg2 = new TH1F("h_bkg2", "other NuECC;" + xaxis_label + ";Events", nbins, xbinslow, xbinshigh);
  TH1F* h_bkg3 = new TH1F("h_bkg3", "NC Pi0;" + xaxis_label + ";Events", nbins, xbinslow, xbinshigh);
  TH1F* h_bkg4 = new TH1F("h_bkg4", "NuMu CC Pi0;" + xaxis_label + ";Events", nbins, xbinslow, xbinshigh);
  TH1F* h_bkg5 = new TH1F("h_bkg5", "Other;" + xaxis_label + ";Events", nbins, xbinslow, xbinshigh);

  TH2F* h_nBlobs_vs_totalBlobEnergy = new TH2F("h_nBlobs_vs_totalBlobEnergy", "nBlobs vs total blob energy", 5, 1, 6, 20, 0, 800);
  TH2F* h_BlobEnergy_vs_distToVertex = new TH2F("h_BlobEnergy_vs_distToVertex", "blob energy vs blob distance to vertex",20,0,500,20,0,2000);

  //Colors
  h_sig->SetFillColor(416);   // green
  h_bkg1->SetFillColor(5);    // yellowish
  h_bkg2->SetFillColor(2);    // red
  h_bkg3->SetFillColor(6);    // magenta
  h_bkg4->SetFillColor(7);    // cyan
  h_bkg5->SetFillColor(4);    // blue

  //Loop over tree
  Long64_t nentries = tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    tree->GetEntry(i);
    
    //if (i > 100) { break; }
    //if (selectionCategory != -999) { continue; } 
    //if (nIsoBlobs == 0) { continue; } 
    bool printEvent = false;
    
    //Calculating some cut values
    double DSCalVisE = CalcDSCALVisE(prong_HCALVisE[0], prong_ECALVisE[0]);
    double ODCalVisE = CalcODCALVisE(prong_ODVisE[0], prong_SideECALVisE[0]);

    double E_avail = CalcEAvail(blob_recoil_E_tracker, blob_recoil_E_ecal, electron_4p_vec[3]); //e_tracker, e_ecal, e_lep
    double modified_E_avail = CalcModEAvail(E_avail, protonEnergy, sec_protons_T_fromdEdx, nSecProtons); //e_avail, e_primary_prot, E_sec_protons, n_sec_protons

    double ESCChi2 = CalcESCChi2(nProtonTrackNodes, proton_nodes_E);

    double furthestUpstreamIsoBlobStartZ = -9999;
    for (int k=0; k < nIsoBlobs; k++){
      if ((nonvtx_iso_blobs_start_position_z_in_prong[k] < furthestUpstreamIsoBlobStartZ) || furthestUpstreamIsoBlobStartZ<-9000){
	furthestUpstreamIsoBlobStartZ = nonvtx_iso_blobs_start_position_z_in_prong[k];
      }
    }

    double isoBlobStartZ_to_vtxZ = furthestUpstreamIsoBlobStartZ - vtx[2];
    
    //fun kinematic stuff
    double proton_3p_vec[3] = {protonPx, protonPy, protonPz}; 

    double electronTheta = CalcTheta(electron_4p_vec[0], electron_4p_vec[1], electron_4p_vec[2]); //also MasterAnaDev_muon_theta? sus...
    double protonTheta = CalcTheta(proton_3p_vec[0], proton_3p_vec[1], proton_3p_vec[2]); //also MasterAnaDev_proton_theta and theta_fromdEdx

    ROOT::Math::XYZVector electron_pt = CalcPtVec(electron_4p_vec[0], electron_4p_vec[1], electron_4p_vec[2]);
    ROOT::Math::XYZVector proton_pt = CalcPtVec(proton_3p_vec[0], proton_3p_vec[1], proton_3p_vec[2]);
    ROOT::Math::XYZVector delta_pt = electron_pt + proton_pt;

    double delta_ptx = CalcDeltaPtx(electron_pt, delta_pt);
    double delta_pty = CalcDeltaPty(electron_pt, delta_pt);
    double alpha_pt = CalcAlphaPt(electron_pt, delta_pt);
    double phi_pt = CalcPhiPt(electron_pt, proton_pt);

    //if (selectionCategory == -999 || modified_E_avail > 0.5) { continue; } 
    //if (selectionCategory == -999 && isoBlobStartZ_to_vtxZ < 0) { printEvent = true; }
    
    if (printOutput){
    //if (printEvent){

      std::cout << " " << std::endl;
      std::cout << "-------------------- EVENT: " << mc_run << " | " << mc_subrun << " | " << mc_nthEvtInFile+1 << " , selectionCategory = " << selectionCategory << ", i= = " << i << " --------------------" << std::endl;
      //std::cout << mc_run << " | " << mc_subrun << " | " << mc_nthEvtInFile+1 << std::endl;
      
      //Cut vars
      //std::cout << "HasNoVertexMismatch" << hasNoVertexMismatch << std::endl;
      //std::cout << "vertex: ( " << vtx[0] << ", " << vtx[1] << ", " << vtx[2] << ", " << vtx[3] << " ) " << std::endl;
      //std::cout << "StartPointVertexMultiplicity: " << startPointVertexMultiplicity << std::endl;
      //std::cout << "FirstFireFraction (afterpulsing): " << prong_firstFireFraction[0] << std::endl;
      //std::cout << "nDeadDiscrPairsUpstream (deadtime): " << nDeadDiscrPairsUpstream << std::endl;
      //std::cout << "HasNoBackExitingTracks: " << hasNoBackExitingTracks << std::endl;
      //std::cout << "DSCalVisE: " << DSCalVisE << std::endl;
      //std::cout << "ODCalVisE: " << ODCalVisE << std::endl;
      //std::cout << "VertexTrackMultiplicity: " << vertexTrackMultiplicity << std::endl;
      //std::cout << "TransverseGapScore: " << prong_TransverseGapScore[0] << std::endl;
      //std::cout << "NonMIPClusFrac: " << prong_NonMIPClusFrac[0] << std::endl;
      //std::cout << "prong_part_score[0]: " << prong_part_score[0] << std::endl;
      
      std::cout << "nMichels: " << nMichels << std::endl;
      std::cout << "Mean Front dEdX: " << prong_dEdXMeanFrontTracker[0] << std::endl;
      std::cout << "Modified E_avail: " << modified_E_avail << std::endl;
      std::cout << "Hybrid ESC/Chi^2: " << ESCChi2 << std::endl;
      std::cout << "# of isolated blobs: " << nIsoBlobs << std::endl;
      std::cout << "Total energy of isolated blobs: " << nonvtx_iso_blobs_energy << std::endl;
      std::cout << "Furthest upstream iso blob start z - vtx Z: " << isoBlobStartZ_to_vtxZ << std::endl;

      std::cout << " n prongs: " << n_prongs << std::endl;
      
      //Kinematic Variables and other interesting "output" vars      
      //std::cout << "Available Energy: " << E_avail << std::endl;
      //std::cout << "Neutrino Energy estimate: " << E_avail+electron_4p_vec[3] << std::endl;
      
      std::cout << "Electron 4 momentum: ( " << electron_4p_vec[0] << ", " << electron_4p_vec[1] << ", " << electron_4p_vec[2] << ", " << electron_4p_vec[3] << " ) " << std::endl;
      std::cout << "Electron Theta: " << electronTheta*180./pi << std::endl;
      //std::cout << "Electron_Pt: ( " << electron_pt.X() << ", " << electron_pt.Y() << ", " << electron_pt.Z() << " ) - Magnitude = " << sqrt(electron_pt.Mag2()) << std::endl;
      
      std::cout << "Proton momentum: ( " << proton_3p_vec[0] << ", " << proton_3p_vec[1] << ", " << proton_3p_vec[2] << " ) - total = " << protonP << std::endl;
      std::cout << "Proton Theta: " << protonTheta*180./pi << std::endl;
      //std::cout << "Proton_Pt: ( " << proton_pt.X() << ", " << proton_pt.Y() << ", " << proton_pt.Z() << " ) - Magnitude = " << sqrt(proton_pt.Mag2()) << std::endl;
      
      //std::cout << "Delta_Pt: ( " << delta_pt.X() << ", " << delta_pt.Y() << ", " << delta_pt.Z() << " ) - Magnitude = " << sqrt(delta_pt.Mag2()) << std::endl;
      //std::cout << "Delta_Ptx: " << delta_ptx << std::endl;
      //std::cout << "Delta_Pty: " << delta_pty << std::endl;
      //std::cout << "Alpha_Pt: " << alpha_pt << std::endl;
      //std::cout << "Phi_Pt: " << phi_pt << std::endl;

      //Truth info
      //for (int j; j++; j<20){
      //if (truth_proton_E[j] < 0) { break; } //the array is filled with -1s after we run out of protons, so skip
	//std::cout << "True proton #" << j+1 << " 4-momentum: ( " << truth_proton_px[j] << ", " << truth_proton_py[j] << ", " << truth_proton_pz[j] << ", " truth_proton_E[j] << " ) " << std::endl;
	//std::cout << "True proton #" << j+1 " theta wrt beam: " << truth_proton_theta_wrtbeam[j] << std::endl;
      //}
    }

    // Fill appropriate histo(s)
    if (nIsoBlobs>0){

      if (selectionCategory == -999)
      h_sig->Fill(isoBlobStartZ_to_vtxZ);
    else if (selectionCategory == 0)
      h_bkg1->Fill(isoBlobStartZ_to_vtxZ);
    else if (selectionCategory == 1)
      h_bkg2->Fill(isoBlobStartZ_to_vtxZ);
    else if (selectionCategory == 2)
      h_bkg3->Fill(isoBlobStartZ_to_vtxZ);
    else if (selectionCategory == 3)
      h_bkg4->Fill(isoBlobStartZ_to_vtxZ);
    else if (selectionCategory == -1)
      h_bkg5->Fill(isoBlobStartZ_to_vtxZ);

      //h_nBlobs_vs_totalBlobEnergy->Fill(nIsoBlobs, nonvtx_iso_blobs_energy);
      //for (int k=0; k < nIsoBlobs; k++){      
      //h_BlobEnergy_vs_distToVertex->Fill(nonvtx_iso_blobs_energy_in_prong[k], nonvtx_iso_blobs_distance_in_prong[k]);
      //}

    }
    
  } // end loop over tree

  //std::cout << "number of events with no isolated blobs: " << nEventsWithNoIsoBlobs << std::endl;
  if (makePlots){
    //Make & plot signal/background stack
    THStack* hs = new THStack("hs", ";" + xaxis_label + ";Events");
    hs->Add(h_bkg5);
    hs->Add(h_bkg4);
    hs->Add(h_bkg3);
    hs->Add(h_bkg2);
    hs->Add(h_bkg1);
    hs->Add(h_sig);
  
    auto c1 = new TCanvas("c1", "Stacked histos", 800, 600);
    hs->Draw("hist");
    c1->BuildLegend(0.65, 0.75, 0.88, 0.88);
    c1->SaveAs("blobZ-vtxZ_stack.png");  
    
    //c1->Clear();
    //h_nBlobs_vs_totalBlobEnergy->Draw("colz");
    //c1->Update();
    //c1->SaveAs("nBlobs_vs_totalBlobEnergy.png");

    //c1->Clear();
    //h_BlobEnergy_vs_distToVertex->Draw("colz");
    //c1->Update();
    //c1->SaveAs("blobEnergy_vs_blobDistToVtx.png");    
  }
}
