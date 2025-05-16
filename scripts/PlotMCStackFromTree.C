//A script that takes an input MAD tuple, and plots a dataMC stack the same way MnvPlotter does
//Input is just hard coded in, I need to add more variables anyways.
//TO DO: Is there a way to get this to actually USE MnvPlotter? Is it even worth trying to figure out?
void PlotMCStackFromTree() {
  TString filename = "/exp/minerva/data/users/cpernas/TUPLE_MC_May_14_2025_me1M_all_cuts_with_Michel_sideband.root";
  TString treename = "Signal_Region";

  TString xaxis_label = "Modified Available Energy (GeV)";

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
  
  const int kMaxProtons = 100; // must be >= max size ever used
  
  double blob_recoil_E_tracker = 0, blob_recoil_E_ecal = 0;
  std::vector<std::vector<double>>* prong_part_E = nullptr;
  double MasterAnaDev_proton_calib_energy = 0;
  double sec_protons_T_fromdEdx[kMaxProtons];
  int sec_protons_T_fromdEdx_sz;
  int selectionCategory = 0;

  tree->SetBranchAddress("blob_recoil_E_tracker", &blob_recoil_E_tracker);
  tree->SetBranchAddress("blob_recoil_E_ecal", &blob_recoil_E_ecal);
  tree->SetBranchAddress("prong_part_E", &prong_part_E);
  tree->SetBranchAddress("MasterAnaDev_proton_calib_energy", &MasterAnaDev_proton_calib_energy);
  tree->SetBranchAddress("MasterAnaDev_sec_protons_T_fromdEdx", sec_protons_T_fromdEdx);
  tree->SetBranchAddress("MasterAnaDev_sec_protons_T_fromdEdx_sz", &sec_protons_T_fromdEdx_sz);
  tree->SetBranchAddress("selectionCategory", &selectionCategory);


  double xbins[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2};
  int nbins = sizeof(xbins)/sizeof(xbins[0])-1;
  // ==== Histograms ====
  TH1F* h_sig  = new TH1F("h_sig",  "Signal;" + xaxis_label + ";Events", nbins, xbins);
  TH1F* h_bkg1 = new TH1F("h_bkg1", "NuECC nonQELike;" + xaxis_label + ";Events", nbins, xbins);
  TH1F* h_bkg2 = new TH1F("h_bkg2", "other NuECC;" + xaxis_label + ";Events", nbins, xbins);
  TH1F* h_bkg3 = new TH1F("h_bkg3", "NC Pi0;" + xaxis_label + ";Events", nbins, xbins);
  TH1F* h_bkg4 = new TH1F("h_bkg4", "NuMu CC Pi0;" + xaxis_label + ";Events", nbins, xbins);
  TH1F* h_bkg5 = new TH1F("h_bkg5", "Other;" + xaxis_label + ";Events", nbins, xbins);

  // ==== Colors ====
  h_sig->SetFillColor(416);   // green
  h_bkg1->SetFillColor(5);    // yellowish
  h_bkg2->SetFillColor(2);    // red
  h_bkg3->SetFillColor(6);    // magenta
  h_bkg4->SetFillColor(7);    // cyan
  h_bkg5->SetFillColor(4);    // blue

  // ==== Fill histograms ====
  Long64_t nentries = tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    tree->GetEntry(i);
    
    // Compute Eavail
    float prongE = 0;
    if (prong_part_E && !prong_part_E->empty() && (*prong_part_E)[0].size() > 3)
      prongE = (*prong_part_E)[0][3];
    
    float E = blob_recoil_E_tracker + blob_recoil_E_ecal;
    float Eavail = (E * 1.17 - (0.008 * prongE + 5)) / 1000.; // GeV
    
    // Modified Eavail
    float primary_proton_E = MasterAnaDev_proton_calib_energy;
    if (primary_proton_E < -1) primary_proton_E = 0;
    
    double sum_secondary_proton_energy = 0.0;
    for (int i = 0; i < sec_protons_T_fromdEdx_sz; ++i) {
      sum_secondary_proton_energy += sec_protons_T_fromdEdx[i];
    }


        float ModEavail = (Eavail * 1000. - primary_proton_E - sum_secondary_proton_energy) / 1000.; // GeV

        // Fill appropriate histo
        if (selectionCategory == -999)
            h_sig->Fill(ModEavail);
        else if (selectionCategory == 0)
            h_bkg1->Fill(ModEavail);
        else if (selectionCategory == 1)
            h_bkg2->Fill(ModEavail);
        else if (selectionCategory == 2)
            h_bkg3->Fill(ModEavail);
        else if (selectionCategory == 3)
            h_bkg4->Fill(ModEavail);
        else if (selectionCategory == -1)
            h_bkg5->Fill(ModEavail);
    }


  // ==== Stack ====
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
  c1->SaveAs("ModifiedEavail_stack.png");
}
