//Another version of PlotMCStackFromTree for simpler variables
//By simple I mean variables that can be pulled or calculated in line directly from the Draw() command
//eg mc_incomingE , as opposed to something more complicated like Modified E_avail (requires looping over protons, etc)
void PlotMCStackFromTree() {
  // ==== CONFIGURATION ====
  TString filename = "/pnfs/minerva/persistent/users/cpernas/default_analysis_loc/TUPLE_MC_May_01_2025_me1B_all_cuts_with_Michel_sideband.root";
  TString treename = "Signal_Region"; 
  //TString branch   = "nonvtx_iso_blobs_start_position_z_in_prong_sz";
  //TString xaxis_label = "N isolated blobs";
  //TString branch   = "nonvtx_iso_blobs_energy";
  //TString xaxis_label = "nonvtx iso blobs energy [MeV]";
  //TString branch   = "event_tracks_energy_sz";
  //TString xaxis_label = "total number of tracks in event";
  
  TString branch   = "multiplicity";
  TString xaxis_label = "multiplicity";

  int nbins = 2; 
  
  // ==== Open File and Get Tree ====
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

  // ==== Canvas ====
  auto c1 = new TCanvas("c1", "Stacked histos", 800, 600);

  //off screen draw to get xmin and xmax for the stack
  tree->Draw(branch, "", "goff");
  double xmin = TMath::MinElement(tree->GetSelectedRows(), tree->GetV1());
  double xmax = TMath::MaxElement(tree->GetSelectedRows(), tree->GetV1());  
  //double xmin = 0;
  //double xmax = 3000;
  
  // ==== Histograms ====
  TH1F* h_sig  = new TH1F("h_sig",  "Signal;" + xaxis_label + ";Events", nbins, xmin, xmax);
  TH1F* h_bkg1 = new TH1F("h_bkg1", "NuECC nonQELike;" + xaxis_label + ";Events", nbins, xmin, xmax);
  TH1F* h_bkg2 = new TH1F("h_bkg2", "other NuECC;" + xaxis_label + ";Events", nbins, xmin, xmax);
  TH1F* h_bkg3 = new TH1F("h_bkg3", "NC Pi0;" + xaxis_label + ";Events", nbins, xmin, xmax);
  TH1F* h_bkg4 = new TH1F("h_bkg4", "NuMu CC Pi0;" + xaxis_label + ";Events", nbins, xmin, xmax);
  TH1F* h_bkg5 = new TH1F("h_bkg5", "Other;" + xaxis_label + ";Events", nbins, xmin, xmax);

  // ==== Colors ====
  h_sig->SetFillColor(416);   // green
  h_bkg1->SetFillColor(5);    // yellowish
  h_bkg2->SetFillColor(2);    // red
  h_bkg3->SetFillColor(6);    // magenta
  h_bkg4->SetFillColor(7);    // cyan
  h_bkg5->SetFillColor(4);    // blue

  // ==== Fill histograms ====
  TString draw_expr = branch + " >> ";

  tree->Draw(draw_expr + "h_sig",  "selectionCategory == -999", "goff");
  tree->Draw(draw_expr + "h_bkg1", "selectionCategory == 0",    "goff");
  tree->Draw(draw_expr + "h_bkg2", "selectionCategory == 1",    "goff");
  tree->Draw(draw_expr + "h_bkg3", "selectionCategory == 2",    "goff");
  tree->Draw(draw_expr + "h_bkg4", "selectionCategory == 3",    "goff");
  tree->Draw(draw_expr + "h_bkg5", "selectionCategory == -1",   "goff");

  // ==== Stack ====
  THStack* hs = new THStack("hs", branch + ";" + xaxis_label + ";Events");
  hs->Add(h_bkg5);
  hs->Add(h_bkg4);
  hs->Add(h_bkg3);
  hs->Add(h_bkg2);
  hs->Add(h_bkg1);
  hs->Add(h_sig);

  hs->Draw("hist");
  c1->BuildLegend(0.65, 0.75, 0.88, 0.88);
}
