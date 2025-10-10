void FitNodeEnergyHistogram() {
  // Load your file
  TFile *f = TFile::Open("your_file.root");
  TTree *t = (TTree*)f->Get("your_tree");

  // Set up the histogram
  TH1D *hNode2 = new TH1D("hNode2", "Energy at Node[2];Energy [MeV];Events", 50, 0, 100);

  // Fill the histogram with node[2] energy
  // Replace 'proton_nodesNormE[2]' with your actual branch access pattern
  t->Draw("proton_nodesNormE[2] >> hNode2", "selection_criteria_here", "goff");

  // Optional: draw before fitting
  hNode2->Draw();

  // Fit a Gaussian
  hNode2->Fit("gaus", "Q"); // "Q" = quiet mode, no print

  // Get the fit
  TF1 *fit = hNode2->GetFunction("gaus");
  double mean  = fit->GetParameter(1);
  double sigma = fit->GetParameter(2);

  std::cout << "Node[2] Mean Energy: " << mean << " MeV" << std::endl;
  std::cout << "Node[2] Sigma:       " << sigma << " MeV" << std::endl;

  // Optional: save the plot
  TCanvas *c = new TCanvas();
  hNode2->Draw("E");
  c->SaveAs("node2_fit.png");
}
