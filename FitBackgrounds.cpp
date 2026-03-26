// FitBackgrounds.cpp
// Takes as input root files from runEventLoop containing signal samples and sideband samples of meanfront dE/dX and >= 1 michel events.
// Calculates scale factors (with multiple options how), creates plots of scale factors, scaled sidebands, scaled signal
// And saves the scaled signal region to a new root file, to be used in ExtractCrossSection.cpp
//
// Usage (example):
//   FitBackgrounds data.root mc.root 1000 Lepton_Pt
//   input data histos, input mc histos, regularization strength lambda, variable name

//For testing, forces all output scale factors to 1. DON'T LEAVE THIS ON BY ACCIDENT
bool set_scale_factors_to_1 = false;

#include "util/GetIngredient.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvVertErrorBand.h"
#include "PlotUtils/MnvPlotter.h"

//ROOT includes
#include "TH1D.h"
#include "TFile.h"
#include "TKey.h"
#include "TParameter.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>

//c++ includes
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

struct ScaleFactors {
  PlotUtils::MnvH1D* meanFrontBkg_mnvhist;
  PlotUtils::MnvH1D* meanFrontSig_mnvhist;
  PlotUtils::MnvH1D* michelBkg_mnvhist;
  PlotUtils::MnvH1D* michelSig_mnvhist;
};

// data_hists: sidebands in order [signalRegion, meanFrontSB, michelSB] (same as your code)
// mc_hists: vector per category: mc_hists[cat][sideband_index]
// mcScale: overall mc normalization
// lambda: regularization strength (0 => no regularization -> exact per-bin solution)
// regularizeSignal: whether to regularize signal scale-factors too (default false)
ScaleFactors ExtractScaleFactors(
    const std::vector<PlotUtils::MnvH1D*>& data_hists,
    const std::vector<std::vector<PlotUtils::MnvH1D*>>& mc_hists,
    double mcScale,
    double lambda,
    bool regularizeSignal = false)
{
  ScaleFactors sf; //create and prep MnvH1Ds for holding scale factors, unique set per sideband and per universe
  sf.meanFrontBkg_mnvhist = dynamic_cast<PlotUtils::MnvH1D*>(mc_hists[0][0]->Clone((mc_hists[0][0]->GetName()+std::string("_clone")).c_str()));
  sf.meanFrontSig_mnvhist = dynamic_cast<PlotUtils::MnvH1D*>(mc_hists[0][0]->Clone((mc_hists[0][0]->GetName()+std::string("_clone")).c_str()));
  sf.michelBkg_mnvhist = dynamic_cast<PlotUtils::MnvH1D*>(mc_hists[0][0]->Clone((mc_hists[0][0]->GetName()+std::string("_clone")).c_str()));
  sf.michelSig_mnvhist = dynamic_cast<PlotUtils::MnvH1D*>(mc_hists[0][0]->Clone((mc_hists[0][0]->GetName()+std::string("_clone")).c_str()));

  sf.meanFrontBkg_mnvhist->Reset("ICES");
  sf.meanFrontSig_mnvhist->Reset("ICES");
  sf.michelBkg_mnvhist->Reset("ICES");
  sf.michelSig_mnvhist->Reset("ICES");
  
  const int nbins = data_hists[0]->GetNbinsX();
  // number of unknowns: for each bin we keep [b(i), s(i)], per fit
  const int nUnknowns = nbins * 2;
  
  auto solve_with_reg = [&](bool forMeanFront) {    
    const int neq_data = 2 * nbins;
    const int nRegBkg = std::max(0, nbins - 1);
    const int nRegSig = (regularizeSignal ? std::max(0, nbins - 1) : 0);
    const int nRegRows = static_cast<int>(lambda > 0 ? (nRegBkg + nRegSig) : 0);

    const int nRows = neq_data + nRegRows;

    std::vector<std::string> vertErrorBandNames = mc_hists[0][0]->GetVertErrorBandNames();
    vertErrorBandNames.push_back("cv"); //add the cv to the list since it doesn't get returned, this way don't have to write a separate loop
    
    size_t band_index = 0;
    //std::cout << "======================= BEGINNING LOOP OVER ERROR BANDS ========================" << std::endl;
    for (const auto& bandName : vertErrorBandNames){//Loop over error bands
      bool isCV = (bandName == "cv");
      PlotUtils::MnvVertErrorBand* band;
      int nHists;
      if (isCV) { nHists = 1; }
      else {
	band = mc_hists[0][0]->GetVertErrorBand( bandName );
	nHists = band->GetNHists();
      }
      //std::cout << "Currently on band #" << band_index << ", which is: " << bandName << " and contains " << nHists << " universes/histograms. " << std::endl;
      for (int universe_index = 0; universe_index < nHists; universe_index++){ //Loop over universes within that error band (normally 2, although flux has 100)
	//std::cout << "---------- universe " << universe_index << " in band " << bandName << ": now looping through its bins. ----------" << std::endl;
	//if (bandName == "GENIE_FrAbs_N") { std::cout << "---------- universe " << universe_index << " in band " << bandName << ": now looping through its bins. ----------" << std::endl; }
	
	TMatrixD A(nRows, nUnknowns); // zero-initialized
	TVectorD d(nRows);            // RHS
	
	for (int ib = 1; ib <= nbins; ++ib) {
	  int row0 = (ib - 1) * 2;
	  double data_signal_region = data_hists[0]->GetBinContent(ib);
	  double data_sb = forMeanFront ? data_hists[1]->GetBinContent(ib) : data_hists[2]->GetBinContent(ib);
	  
	  double d_s, d_sb, other_s, other_sb;
	  double mc_sig_s, mc_sig_sb, mc_bkg_s, mc_bkg_sb;

	  if (isCV){
	    mc_sig_s = mc_hists[0][0]->GetBinContent(ib) * mcScale; // signal in signal region
	    mc_sig_sb = mc_hists[0][ forMeanFront ? 1 : 2 ]->GetBinContent(ib) * mcScale; // signal in corresponding sideband
	  
	    if (forMeanFront) {
	      other_s = (mc_hists[1][0]->GetBinContent(ib) + mc_hists[2][0]->GetBinContent(ib) + mc_hists[5][0]->GetBinContent(ib)) * mcScale;
	      other_sb= (mc_hists[1][1]->GetBinContent(ib) + mc_hists[2][1]->GetBinContent(ib) + mc_hists[5][1]->GetBinContent(ib)) * mcScale;
	      mc_bkg_s  = (mc_hists[3][0]->GetBinContent(ib) + mc_hists[4][0]->GetBinContent(ib)) * mcScale;
	      mc_bkg_sb = (mc_hists[3][1]->GetBinContent(ib) + mc_hists[4][1]->GetBinContent(ib)) * mcScale;
	    } else {
	      other_s = (mc_hists[2][0]->GetBinContent(ib) + mc_hists[3][0]->GetBinContent(ib) + mc_hists[4][0]->GetBinContent(ib) + mc_hists[5][0]->GetBinContent(ib)) * mcScale;
	      other_sb= (mc_hists[2][2]->GetBinContent(ib) + mc_hists[3][2]->GetBinContent(ib) + mc_hists[4][2]->GetBinContent(ib) + mc_hists[5][2]->GetBinContent(ib)) * mcScale;
	      mc_bkg_s  = (mc_hists[1][0]->GetBinContent(ib)) * mcScale;
	      mc_bkg_sb = (mc_hists[1][2]->GetBinContent(ib)) * mcScale; // for michel, index 2 sideband	    
	    }
	  }
	  else { //not cv, have to grab specifically by errorband and universe. Dear lord this is ugly
	    mc_sig_s = mc_hists[0][0]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib) * mcScale; // signal in signal region
	    mc_sig_sb = mc_hists[0][ forMeanFront ? 1 : 2 ]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib) * mcScale; // signal in corresponding sideband
	  
	    if (forMeanFront) {
	      other_s = (mc_hists[1][0]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib) + mc_hists[2][0]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib) + mc_hists[5][0]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib)) * mcScale;
	      other_sb= (mc_hists[1][1]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib) + mc_hists[2][1]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib) + mc_hists[5][1]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib)) * mcScale;
	      mc_bkg_s  = (mc_hists[3][0]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib) + mc_hists[4][0]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib)) * mcScale;
	      mc_bkg_sb = (mc_hists[3][1]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib) + mc_hists[4][1]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib)) * mcScale;
	    } else {
	      other_s = (mc_hists[2][0]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib) + mc_hists[3][0]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib) + mc_hists[4][0]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib) + mc_hists[5][0]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib)) * mcScale;
	      other_sb= (mc_hists[2][2]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib) + mc_hists[3][2]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib) + mc_hists[4][2]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib) + mc_hists[5][2]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib)) * mcScale;
	      mc_bkg_s  = (mc_hists[1][0]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib)) * mcScale;
	      mc_bkg_sb = (mc_hists[1][2]->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib)) * mcScale; // for michel, index 2 sideband	    
	    }
	  }

	  d_s = data_signal_region - other_s;
	  d_sb = data_sb - other_sb;
	    
	  int col_b = (ib - 1) * 2 + 0; // background position for bin i
	  int col_s = (ib - 1) * 2 + 1; // signal position for bin i
	  
	  A(row0, col_b) = mc_bkg_s;
	  A(row0, col_s) = mc_sig_s;
	  d(row0)       = d_s;
	  
	  A(row0+1, col_b) = mc_bkg_sb;
	  A(row0+1, col_s) = mc_sig_sb;
	  d(row0+1)       = d_sb;
	} // end bin loop for data eqs
	
	// Add regularization rows after the data rows
	int regRowStart = neq_data;
	int r = 0;
	if (lambda > 0) {
	  // background diffs
	  for (int ib = 1; ib <= nbins - 1; ++ib) {
	    int row = regRowStart + r++;
	    double w = std::sqrt(lambda);
	    int col_b_i   = (ib - 1) * 2 + 0;
	    int col_b_ip1 = (ib    ) * 2 + 0;
	    A(row, col_b_i)   =  w;
	    A(row, col_b_ip1) = -w;
	    d(row) = 0.0;
	  }
	  // optional signal diffs
	  if (regularizeSignal) {
	    for (int ib = 1; ib <= nbins - 1; ++ib) {
	      int row = regRowStart + r++;
	      double w = std::sqrt(lambda);
	      int col_s_i   = (ib - 1) * 2 + 1;
	      int col_s_ip1 = (ib    ) * 2 + 1;
	      A(row, col_s_i)   =  w;
	      A(row, col_s_ip1) = -w;
	      d(row) = 0.0;
	    }
	  }
	}
	
	// Solve least-squares: minimize ||A x - d|| using SVD
	TDecompSVD svd(A);
	Bool_t ok;
	TVectorD rhs = d;            // copy RHS because Solve modifies it
	TVectorD x = svd.Solve(rhs, ok); // vector of length nUnknowns: [b1,s1,b2,s2,...]
	if (!ok) {
	  std::cerr << "Warning: SVD solve failed (singular). Returning unity scales.\n";
	  TVectorD one(nUnknowns);
	  for (int i = 0; i < nUnknowns; ++i) one(i) = 1.0;
	  //return one;
	}
	//Manually set everything to 1 for testing, DON'T LEAVE THIS ON BY ACCIDENT
	if (set_scale_factors_to_1){
	  for (int i = 0; i < nUnknowns; ++i) x(i) = 1.0;
	}
	//Write output to my struct's MnvH1Ds, per universe.
	if (isCV){
	  for (int ib = 0; ib < nbins; ++ib) { //double check that this is right, setBinContent(0, ...) does the underflow bin I think...
	    if (forMeanFront){
	      sf.meanFrontBkg_mnvhist->SetBinContent(ib+1, x(2*ib + 0));
	      sf.meanFrontSig_mnvhist->SetBinContent(ib+1, x(2*ib + 1));
	    } else {
	      sf.michelBkg_mnvhist->SetBinContent(ib+1, x(2*ib + 0));
	      sf.michelSig_mnvhist->SetBinContent(ib+1, x(2*ib + 1));
	    }
	  }
	}
	else { //not the CV, so all other universes
	  for (int ib = 0; ib < nbins; ++ib) { //double check that this is right, setBinContent(0, ...) does the underflow bin I think...
	    if (forMeanFront){
	      sf.meanFrontBkg_mnvhist->GetVertErrorBand(bandName)->GetHist(universe_index)->SetBinContent(ib+1, x(2*ib + 0));
	      sf.meanFrontSig_mnvhist->GetVertErrorBand(bandName)->GetHist(universe_index)->SetBinContent(ib+1, x(2*ib + 1));
	    } else {
	      sf.michelBkg_mnvhist->GetVertErrorBand(bandName)->GetHist(universe_index)->SetBinContent(ib+1, x(2*ib + 0));
	      sf.michelSig_mnvhist->GetVertErrorBand(bandName)->GetHist(universe_index)->SetBinContent(ib+1, x(2*ib + 1));
	    } 
	  }  //end bin loop for output writing
	}  //else statement for nonCV universes
      }  //end loop over universes within an error band
      band_index++;
    }  //end loop over error bands
  }; // end lambda solve_with_reg

  // run meanFront fit
  solve_with_reg(true);
  // run michel fit
  solve_with_reg(false);

  return sf;
} 


void saveSFPlot(PlotUtils::MnvH1D* mnvhist, const std::string& filename) {
  std::unique_ptr<TCanvas> c(new TCanvas(("c_" + filename).c_str(), "", 800, 600));
  c->SetGrid();
  
  // Make a local clone but let the canvas own it
  TH1D cvhist = mnvhist->GetCVHistoWithError();
  TH1D* hist = &cvhist;
  
  hist->SetLineWidth(3);
  hist->SetLineColor(kBlue + 1);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1.0);
  hist->SetMarkerColor(kBlue + 1);
  hist->GetYaxis()->SetTitle("Scale Factor");
  
  hist->Draw("E1");
  gPad->SetTicks();
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  
  c->SaveAs((filename + ".png").c_str());    
};

void saveStackPlot(PlotUtils::MnvH1D* data, const std::vector<PlotUtils::MnvH1D*>& mc_scaled,
		   const std::string& outName, const std::string& titleSuffix, double dataPOT, double mcPOT) {
    // set titles (use python's ordering & labels)
    std::vector<std::string> labels = {
      "signal (nu_e QELike + proton)",
      "nu_e nonQE (has FS mesons)",
      "Other nu_eCC",
      "NC with pi0",
      "nu_mu CC with pi0",
      "other"
    };
    PlotUtils::MnvPlotter plotter;
    plotter.legend_text_size = 0.015;
    plotter.data_line_width = 2;
    plotter.data_marker_size = 1.5;
    
    // MC category colors 
    const std::vector<int> mcColors = {4, 7, 6, 2, 5, 416};
    int arr_int[6];
    for (size_t i=0;i<mcColors.size();++i) arr_int[i]=mcColors[i];
    int* arr = arr_int;
    const double mcScale = dataPOT / mcPOT;

    for (size_t c=0;c<mc_scaled.size();++c) {
      mc_scaled[c]->SetTitle(labels[c].c_str());
      mc_scaled[c]->SetLineColor(kBlack);
      mc_scaled[c]->SetFillColor(mcColors[c]);
      mc_scaled[c]->SetLineWidth(1);
    }
    plotter.mc_line_width = 2;
    data->SetTitle("data");
    // Create TObjArray in reverse order so the stack looks like python (signal on top)
    TObjArray array;
    array.SetOwner(false);
    for (int k = (int)mc_scaled.size()-1; k >= 0; --k) { array.Add(mc_scaled[k]); }

    std::unique_ptr<TCanvas> c(new TCanvas("c", "", 1200, 900));
    //plotter.DrawDataStackedMC(data_hists[s], &arr, nullptr, mcScale, "TR", "Data", 1001, data_hists[0]->GetXaxis()->GetTitle(), "N events");
    plotter.DrawDataStackedMC(data, &array, arr, mcScale, "TR", "Data", 1001, data->GetXaxis()->GetTitle(), "N events");
    plotter.AddPOTNormBox(dataPOT, mcPOT, 0.3, 0.85);
    c->SaveAs(outName.c_str());
  };

// Copy all top-level keys from inputFilePath whose name contains 'prefix'
// except those in skipNames. Also copy any TParameter with "POT" in the name
// and explicitly try to copy "POTUsed".
//
void CopyObjectsWithPrefix(TFile* inFile,
                                   TFile* outFile,
                                   const std::vector<std::string>& skipNames,
                                   const std::string& prefix)
{
  const char* inputFilePath = "dummy";
  if (!inFile || inFile->IsZombie()) {
    std::cerr << "ERROR: could not open input file: " << inputFilePath << std::endl;
    return;
  }

  // pointer to list of keys in the input file
  TList* keys = inFile->GetListOfKeys();
  if (!keys) {
    std::cerr << "ERROR: input file has no keys: " << inputFilePath << std::endl;
    return;
  }

  const Int_t nkeys = keys->GetSize();
  outFile->cd();

  for (Int_t i = 0; i < nkeys; ++i) {
    TKey* key = dynamic_cast<TKey*>(keys->At(i));
    if (!key) continue;

    std::string kname = key->GetName();

    // Only consider keys with the requested prefix
    if (kname.find(prefix) == std::string::npos) continue;

    // Skip names explicitly listed in skipNames
    if (std::find(skipNames.begin(), skipNames.end(), kname) != skipNames.end()) {
      // std::cout << "Skipping modified object: " << kname << std::endl;
      continue;
    }

    // Read the object (allocates an object)
    TObject* obj = key->ReadObj();
    if (!obj) {
      std::cerr << "Warning: cannot read object " << kname << " from " << inputFilePath << std::endl;
      continue;
    }

    // Write into output file with same key name; overwrite if present
    outFile->cd();
    obj->Write(kname.c_str(), TObject::kOverwrite);
  }

  // Try to copy POTUsed explicitly if present
  TKey* kPot = dynamic_cast<TKey*>(keys->FindObject("POTUsed"));
  if (kPot) {
    TObject* potObj = kPot->ReadObj();
    if (potObj) {
      outFile->cd();
      potObj->Write("POTUsed", TObject::kOverwrite);
    }
  } else {
    // Fallback: copy any TParameter-like keys containing "POT"
    for (Int_t i = 0; i < nkeys; ++i) {
      TKey* key = dynamic_cast<TKey*>(keys->At(i));
      if (!key) continue;
      std::string kname = key->GetName();
      if (kname.find("POT") == std::string::npos) continue;

      // Optionally check class: copy only TParameter or numeric parameters
      const char* cls = key->GetClassName();
      if (cls && (std::string(cls).find("TParameter") != std::string::npos || std::string(cls).find("TObjString") != std::string::npos)) {
        TObject* obj = key->ReadObj();
        if (!obj) continue;
        outFile->cd();
        obj->Write(kname.c_str(), TObject::kOverwrite);
      }
    }
  }
}

int main(int argc, char** argv) {
  TH1::AddDirectory(kFALSE); // avoid ownership issues with ROOT directories
  gROOT->SetBatch(kTRUE);
  
  if (argc < 2) {
    std::cerr << "USAGE: " << argv[0] << " <data.root> <mc.root> [variableName=DeltaPt] [lambda=1000]\n";
    return 1;
  }

  const char* dataPath = argv[1];
  const char* mcPath   = argv[2];
  double lambda = std::stod(argv[3]);
  std::string varName = "DeltaPt";
  if (argc >= 5) varName = argv[4]; //if a variable name is provided, use it, otherwise default to DeltaPt

  std::cout << "varName = " << varName << std::endl;
  TFile* dataFile = TFile::Open(dataPath, "READ");
  if (!dataFile || dataFile->IsZombie()) {
    std::cerr << "Failed to open data file: " << dataPath << std::endl;
    return 2;
  }
  TFile* mcFile = TFile::Open(mcPath, "READ");
  if (!mcFile || mcFile->IsZombie()) {
    std::cerr << "Failed to open MC file: " << mcPath << std::endl;
    return 3;
  }

  // POT scaling
  double mcPOT = 1.0, dataPOT = 1.0;
  try {
    auto mp = util::GetIngredient<TParameter<double>>(*mcFile, "POTUsed");
    if (mp) mcPOT = mp->GetVal();
  } catch (...) {
    if (mcFile->Get("POTUsed")) {
      TParameter<double>* p = nullptr;
      mcFile->GetObject("POTUsed", p);
      if (p) mcPOT = p->GetVal();
    }
  }
  try {
    auto dp = util::GetIngredient<TParameter<double>>(*dataFile, "POTUsed");
    if (dp) dataPOT = dp->GetVal();
  } catch (...) {
    if (dataFile->Get("POTUsed")) {
      TParameter<double>* p = nullptr;
      dataFile->GetObject("POTUsed", p);
      if (p) dataPOT = p->GetVal();
    }
  }
  const double mcScale = dataPOT / mcPOT;
  std::cout << "mc POT scale = " << mcScale << "  (dataPOT=" << dataPOT << ", mcPOT=" << mcPOT << ")\n";

  const std::vector<std::string> bkgdCategoryNames = {
    "selected_signal_reco", "background_NuECC_with_pions", "background_Other_NueCC",
    "background_NC_pi0", "background_CC_Numu_pi0", "background_Other"
  };
  const std::vector<std::string> sidebands = {"_", "_MeanFrontDEDXSB_", "_MichelSB_"};

  // Load histograms (MnvH1D) for data & MC
  std::vector<PlotUtils::MnvH1D*> data_hists;
  std::vector<std::vector<PlotUtils::MnvH1D*>> mc_hists(6); // 6 categories, each has 3 sidebands. First index is background category(0 to 5), second index is sideband region (0=sig region, 1=meanFrontdEdX, 2=michel)

  for (size_t s = 0; s < sidebands.size(); ++s) {
    std::string dataName = varName + sidebands[s] + "data";
    PlotUtils::MnvH1D* d = nullptr;
    dataFile->GetObject(dataName.c_str(), d);
    if (!d) {
      std::cerr << "ERROR: data hist " << dataName << " not found in data file\n";
      return 10;
    }
    data_hists.push_back(d);

    for (size_t c = 0; c < bkgdCategoryNames.size(); ++c) {
      std::string mname = varName + sidebands[s] + bkgdCategoryNames[c];
      PlotUtils::MnvH1D* m = nullptr;
      mcFile->GetObject(mname.c_str(), m);
      if (!m) {
        std::cerr << "ERROR: mc histogram " << mname << " not found in mc file\n";
        return 11;
      }
      mc_hists[c].push_back(m);
    }
  }

  ScaleFactors sfs = ExtractScaleFactors(data_hists, mc_hists, mcScale, lambda);
  // Save each scale factor histogram
  saveSFPlot(sfs.meanFrontBkg_mnvhist, "meanFront_bkg_scale_factors");
  saveSFPlot(sfs.meanFrontSig_mnvhist, "meanFront_sig_scale_factors");
  saveSFPlot(sfs.michelBkg_mnvhist,  "michel_bkg_scale_factors");
  saveSFPlot(sfs.michelSig_mnvhist,  "michel_sig_scale_factors");
 
  // --- lil lambda function to build scaled versions for a given sideband index s:
  // For s==0 (signal region): leave signal unscaled, apply michel SF to 1, and apply meanFront SF to 3&4
  // For s==1 (MeanFront): apply meanFront bkg SF to categories 3 & 4, and sig SF to cat 0
  // For s==2 (Michel): apply michel bkg SF to category 1, and sig SF to cat 0
  auto applyScaleFactors = [&](size_t s) -> std::vector<PlotUtils::MnvH1D*> {
    std::vector<PlotUtils::MnvH1D*> scaled(bkgdCategoryNames.size(), nullptr);

    for (size_t c = 0; c < bkgdCategoryNames.size(); ++c) {
      PlotUtils::MnvH1D* base = mc_hists[c][s];
      PlotUtils::MnvH1D* clone = dynamic_cast<PlotUtils::MnvH1D*>(base->Clone((std::string(base->GetName())+std::string("_scaled")).c_str()));

      //Loop through the categories, apply correct scale factors for region and category
      //categories 2 (other NuECC, red) and 5 (other, dark blue) remain completely unscaled, at least for now
     
      if (c == 0) { //Signal category (green), only scale these for the sidebands (s=1&2) to make validation plots.
	if (s == 1){ clone->Multiply(clone, sfs.meanFrontSig_mnvhist); } 
	else if (s == 2){ clone->Multiply(clone, sfs.michelSig_mnvhist); }
      }
      if (c == 1) { //NonQELike yellow category, scale by michel scale factors
	if (s == 0 || s == 2){ clone->Multiply(clone, sfs.michelBkg_mnvhist); }
      }
      if (c == 3 || c == 4) { //NC Pi0 (purple) and NumuCC Pi0 (teal), scale by mean front dE/dX scale factors
	if (s == 0 || s == 1){ clone->Multiply(clone, sfs.meanFrontBkg_mnvhist); }
      }

      scaled[c] = clone;
      }
    return scaled;
  };
      
  auto finalSignalScaled = applyScaleFactors(0); //signal region, signal is NOT scaled
  auto meanScaled = applyScaleFactors(1); //mean front region, with signal also scaled
  auto michelScaled = applyScaleFactors(2); //michel region, with signal also scaled

  // meanFront SB plot
  saveStackPlot(data_hists[1], meanScaled, (varName + "_MeanFrontSB_scaled.png"), "MeanFront SB scaled", dataPOT, mcPOT);
  // michel SB plot 
  saveStackPlot(data_hists[2], michelScaled, (varName + "_MichelSB_scaled.png"), "Michel SB scaled", dataPOT, mcPOT);
  // final signal region
  saveStackPlot(data_hists[0], finalSignalScaled, (varName + "_SignalRegion_finalScaled.png"), "Signal region final scaled", dataPOT, mcPOT);

  // ---------------------------
  // Write output root file containing the signal-region scaled MC histograms with original names.
  // give em the same names as the originals so ExtractCrossSection works as intended
  // ---------------------------
  TFile* outFile = TFile::Open("scaled_mc.root", "RECREATE");
  if (!outFile || outFile->IsZombie()) {
    std::cerr << "ERROR: couldn't open scaled_mc.root for writing\n";
  } else {
    std::vector<std::string> modifiedNames = {
      varName + "_selected_signal_reco",
      varName + "_background_NuECC_with_pions",
      varName + "_background_Other_NueCC",
      varName + "_background_NC_pi0",
      varName + "_background_CC_Numu_pi0",
      varName + "_background_Other"
    };
    //outFile->cd();
    // Copy all other objects from the mc input with prefix varName
    // (preserving key names), and copy POTUsed/TParameter(POT*) as well.
    CopyObjectsWithPrefix(mcFile, outFile, modifiedNames, varName);

    for (size_t c = 0; c < bkgdCategoryNames.size(); ++c) {
      std::string origName = mc_hists[c][0]->GetName(); // original name in MC file
      finalSignalScaled[c]->SetName(origName.c_str());
      finalSignalScaled[c]->Write(); //gets deleted at the very end
    }    

    sfs.meanFrontBkg_mnvhist->SetName("MeanFrontDEDXSB_bkg_scale_factors");
    sfs.meanFrontSig_mnvhist->SetName("MeanFrontDEDXSB_sig_scale_factors");
    sfs.michelBkg_mnvhist->SetName("MichelSB_bkg_scale_factors");
    sfs.michelSig_mnvhist->SetName("MichelSB_sig_scale_factors");

    sfs.meanFrontBkg_mnvhist->Write();
    sfs.meanFrontSig_mnvhist->Write();
    sfs.michelBkg_mnvhist->Write();
    sfs.michelSig_mnvhist->Write();
        
    std::cout << "Writing complete\n";
    // flush and close
    outFile->Close();
    delete outFile;
    std::cout << "outFile closed and deleted successfully\n";
  }

  delete dataFile;
  delete mcFile;
  
  gROOT->GetListOfFunctions()->Delete();
  std::cout << "Done.\n";
  return 0;
}
