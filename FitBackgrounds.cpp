// FitBackgrounds_fromPython_fixed.cpp
// Converted from your scale_factors_test.py and adjusted to:
//  - compute scale factors identical to the python script
//  - make only three plots: meanFront sideband, michel sideband, and final signal region
//  - write signal-region scaled MC histograms (with same names) into scaled_mc.root
//
// Usage (example):
//   root -l -b -q 'FitBackgrounds_fromPython_fixed.cpp("data.root","mc.root",3,"DeltaPt")'
//
// Notes:
//  - Requires PlotUtils (MnvH1D, MnvPlotter) available to ROOT/CINT (ROOT_INCLUDE_PATH etc).
//  - Compiles in cling / ROOT JIT (works as a ROOT macro). If you prefer a compiled binary, integrate
//    into your build system and link the same libraries.

#include "util/GetIngredient.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#include "MinervaUnfold/MnvUnfold.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"
#pragma GCC diagnostic pop

#include "TH1D.h"
#include "TFile.h"
#include "TKey.h"
#include "TClass.h"
#include "TObject.h"
#include "TDirectory.h"
#include "TParameter.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TList.h"
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>

#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

struct ScaleFactors {
  std::vector<double> meanFrontBkg;
  std::vector<double> meanFrontSig;
  std::vector<double> michelBkg;
  std::vector<double> michelSig;
};

ScaleFactors ExtractScaleFactors(
    const std::vector<PlotUtils::MnvH1D*>& data_hists,
    const std::vector<std::vector<PlotUtils::MnvH1D*>>& mc_hists,
    double mcScale)
{
  ScaleFactors sf;

  const int nbins = data_hists[0]->GetNbinsX();
  sf.meanFrontBkg.resize(nbins, 1.0);
  sf.meanFrontSig.resize(nbins, 1.0);
  sf.michelBkg.resize(nbins, 1.0);
  sf.michelSig.resize(nbins, 1.0);

  for (int ib = 1; ib <= nbins; ++ib) {
    double data_signal_region = data_hists[0]->GetBinContent(ib);
    double data_meanFront_sb = data_hists[1]->GetBinContent(ib);
    double data_michel_sb = data_hists[2]->GetBinContent(ib);

    // MeanFront sideband calculation
    double d_s = data_signal_region - (
      (mc_hists[1][0]->GetBinContent(ib) + mc_hists[2][0]->GetBinContent(ib) + mc_hists[5][0]->GetBinContent(ib)) * mcScale
    );
    double d_sb = data_meanFront_sb - (
      (mc_hists[1][1]->GetBinContent(ib) + mc_hists[2][1]->GetBinContent(ib) + mc_hists[5][1]->GetBinContent(ib)) * mcScale
    );

    double mc_sig_s = mc_hists[0][0]->GetBinContent(ib) * mcScale;
    double mc_sig_sb = mc_hists[0][1]->GetBinContent(ib) * mcScale;
    double mc_bkg_s  = (mc_hists[3][0]->GetBinContent(ib) + mc_hists[4][0]->GetBinContent(ib)) * mcScale;
    double mc_bkg_sb = (mc_hists[3][1]->GetBinContent(ib) + mc_hists[4][1]->GetBinContent(ib)) * mcScale;

    double delta = mc_sig_s * mc_bkg_sb - mc_sig_sb * mc_bkg_s;
    double bkg_scale_1 = 1.0, sig_scale_1 = 1.0;
    if (std::fabs(delta) > 1e-12) {
      bkg_scale_1 = ( mc_sig_sb * (-d_s) + mc_sig_s * d_sb ) / delta;
      sig_scale_1 = ( d_s * mc_bkg_sb - d_sb * mc_bkg_s ) / delta;
    } else {
      bkg_scale_1 = 1.0;
      sig_scale_1 = 1.0;
    }
    sf.meanFrontBkg[ib-1] = bkg_scale_1;
    sf.meanFrontSig[ib-1] = sig_scale_1;

    // Michel sideband calculation
    d_s = data_signal_region - (
      (mc_hists[2][0]->GetBinContent(ib) + mc_hists[3][0]->GetBinContent(ib) + mc_hists[4][0]->GetBinContent(ib) + mc_hists[5][0]->GetBinContent(ib)) * mcScale
    );
    d_sb = data_michel_sb - (
      (mc_hists[2][2]->GetBinContent(ib) + mc_hists[3][2]->GetBinContent(ib) + mc_hists[4][2]->GetBinContent(ib) + mc_hists[5][2]->GetBinContent(ib)) * mcScale
    );

    mc_sig_s = mc_hists[0][0]->GetBinContent(ib) * mcScale;
    mc_sig_sb = mc_hists[0][2]->GetBinContent(ib) * mcScale;
    mc_bkg_s  = (mc_hists[1][0]->GetBinContent(ib)) * mcScale;
    mc_bkg_sb = (mc_hists[1][2]->GetBinContent(ib)) * mcScale;

    delta = mc_sig_s * mc_bkg_sb - mc_sig_sb * mc_bkg_s;
    double bkg_scale_2 = 1.0, sig_scale_2 = 1.0;
    if (std::fabs(delta) > 1e-12) {
      bkg_scale_2 = ( mc_sig_sb * (-d_s) + mc_sig_s * d_sb ) / delta;
      sig_scale_2 = ( d_s * mc_bkg_sb - d_sb * mc_bkg_s ) / delta;
    } else {
      bkg_scale_2 = 1.0;
      sig_scale_2 = 1.0;
    }
    sf.michelBkg[ib-1] = bkg_scale_2;
    sf.michelSig[ib-1] = sig_scale_2;
  } // end bin loop  

  return sf;
}

// data_hists: sidebands in order [signalRegion, meanFrontSB, michelSB] (same as your code)
// mc_hists: vector per category: mc_hists[cat][sideband_index]
// mcScale: overall mc normalization
// lambda: regularization strength (0 => no regularization -> exact per-bin solution)
// regularizeSignal: whether to regularize signal scale-factors too (default false)
ScaleFactors ExtractScaleFactorsWithReg(
    const std::vector<PlotUtils::MnvH1D*>& data_hists,
    const std::vector<std::vector<PlotUtils::MnvH1D*>>& mc_hists,
    double mcScale,
    double lambda,
    bool regularizeSignal = false)
{
  const int nbins = data_hists[0]->GetNbinsX();
  // number of unknowns: for each bin we keep [b(i), s(i)], per fit
  const int nUnknowns = nbins * 2;

  // Number of data equations: 2 per bin (signal region and that sideband's equation)
  // We'll have two independent sets (meanFront and michel). We'll build two separate fits:
  //  - meanFront fit gives meanFront bkg/sig vectors
  //  - michel fit gives michel bkg/sig vectors
  // Implement as one function reused twice (so we avoid duplicated code).
  auto solve_with_reg = [&](bool forMeanFront) {
    // Build A and d for all bins stacked. Each bin has 2 equations -> total 2*nbins rows.
    // Then add regularization rows: (nbins-1) rows for background diffs; optionally also for signal diffs.

    const int neq_data = 2 * nbins;
    const int nRegBkg = std::max(0, nbins - 1);
    const int nRegSig = (regularizeSignal ? std::max(0, nbins - 1) : 0);
    const int nRegRows = static_cast<int>(lambda > 0 ? (nRegBkg + nRegSig) : 0);

    const int nRows = neq_data + nRegRows;

    TMatrixD A(nRows, nUnknowns); // zero-initialized
    TVectorD d(nRows);            // RHS

    // fill data equations
    for (int ib = 1; ib <= nbins; ++ib) {
      int row0 = (ib - 1) * 2;
      // compute d_s and d_sb as in your code but adjusted for which sideband
      double data_signal_region = data_hists[0]->GetBinContent(ib);
      double data_meanFront_sb   = data_hists[1]->GetBinContent(ib);
      double data_michel_sb      = data_hists[2]->GetBinContent(ib);

      // For meanFront fit: use meanFront sideband (sideband index 1)
      // For michel fit: use michel sideband (sideband index 2)
      double data_sb = forMeanFront ? data_meanFront_sb : data_michel_sb;

      // compute "d" as data minus other (fixed) MC categories *mcScale
      // you must mirror your python logic: which MC components are "other" (not scaled)
      // I'll use the same combinations you had in your working code:
      // -- for meanFront: other backgrounds in the d_s / d_sb subtraction were categories 1,2,5 (these were excluded from the 2x2 solve)
      // -- for michel: other backgrounds were categories 2,3,4,5 per your earlier code. Adjust this block if your exact grouping differs.

      double d_s = 0;
      double d_sb = 0;

      if (forMeanFront) {
        // in your earlier meanFront calculation:
        // d_s = data_signal_region - (mc_hists[1][0] + mc_hists[2][0] + mc_hists[5][0]) * mcScale
        // d_sb= data_meanFront_sb   - (mc_hists[1][1] + mc_hists[2][1] + mc_hists[5][1]) * mcScale
        double other_s = (mc_hists[1][0]->GetBinContent(ib) + mc_hists[2][0]->GetBinContent(ib) + mc_hists[5][0]->GetBinContent(ib)) * mcScale;
        double other_sb= (mc_hists[1][1]->GetBinContent(ib) + mc_hists[2][1]->GetBinContent(ib) + mc_hists[5][1]->GetBinContent(ib)) * mcScale;
        d_s = data_signal_region - other_s;
        d_sb = data_sb - other_sb;
        // mc_sig_s, mc_sig_sb, mc_bkg_s, mc_bkg_sb as in your code
      } else {
        // michel sideband calculation earlier:
        // d_s = data_signal_region - (mc_hists[2][0] + mc_hists[3][0] + mc_hists[4][0] + mc_hists[5][0]) * mcScale
        // d_sb= data_michel_sb - (mc_hists[2][2]+mc_hists[3][2]+mc_hists[4][2]+mc_hists[5][2]) * mcScale
        double other_s = (mc_hists[2][0]->GetBinContent(ib) + mc_hists[3][0]->GetBinContent(ib) + mc_hists[4][0]->GetBinContent(ib) + mc_hists[5][0]->GetBinContent(ib)) * mcScale;
        double other_sb= (mc_hists[2][2]->GetBinContent(ib) + mc_hists[3][2]->GetBinContent(ib) + mc_hists[4][2]->GetBinContent(ib) + mc_hists[5][2]->GetBinContent(ib)) * mcScale;
        d_s = data_signal_region - other_s;
        d_sb= data_sb - other_sb;
      }

      // mc components that multiply s and b in the 2x2 system:
      double mc_sig_s = mc_hists[0][0]->GetBinContent(ib) * mcScale; // signal in signal region
      double mc_sig_sb = mc_hists[0][ forMeanFront ? 1 : 2 ]->GetBinContent(ib) * mcScale; // signal in corresponding sideband

      double mc_bkg_s, mc_bkg_sb;
      if (forMeanFront) {
        // meanFront bkg is categories 3+4? In your earlier code meanFront bkg used categories 3 & 4 (I think)
        mc_bkg_s  = (mc_hists[3][0]->GetBinContent(ib) + mc_hists[4][0]->GetBinContent(ib)) * mcScale;
        mc_bkg_sb = (mc_hists[3][1]->GetBinContent(ib) + mc_hists[4][1]->GetBinContent(ib)) * mcScale;
      } else {
        // michel bkg used category 1 only in your earlier code
        mc_bkg_s  = (mc_hists[1][0]->GetBinContent(ib)) * mcScale;
        mc_bkg_sb = (mc_hists[1][ forMeanFront ? 1 : 2 ]->GetBinContent(ib)) * mcScale; // for michel, index 2 sideband
      }

      // row row0: mc_sig_s * s + mc_bkg_s * b = d_s
      // row row0+1: mc_sig_sb * s + mc_bkg_sb * b = d_sb
      // But remember our unknown order is [b(i), s(i)] i.e. b then s
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
    TVectorD x = svd.Solve(rhs, ok);
    if (!ok) {
      std::cerr << "Warning: SVD solve failed (singular). Returning unity scales.\n";
      TVectorD one(nUnknowns);
      for (int i = 0; i < nUnknowns; ++i) one(i) = 1.0;
      return one;
    }
    return x; // vector of length nUnknowns: [b1,s1,b2,s2,...]
  }; // end lambda solve_with_reg

  // run meanFront fit
  TVectorD x_mean = solve_with_reg(true);
  // run michel fit
  TVectorD x_michel = solve_with_reg(false);

  // unpack into output vectors
  ScaleFactors out;
  out.meanFrontBkg.resize(nbins);
  out.meanFrontSig.resize(nbins);
  out.michelBkg.resize(nbins);
  out.michelSig.resize(nbins);

  for (int ib = 0; ib < nbins; ++ib) {
    out.meanFrontBkg[ib] = x_mean(2*ib + 0);
    out.meanFrontSig[ib] = x_mean(2*ib + 1);
    out.michelBkg[ib]    = x_michel(2*ib + 0);
    out.michelSig[ib]    = x_michel(2*ib + 1);
  }

  return out;
} 

// Copy all top-level keys from inputFilePath whose name contains 'prefix'
// except those in skipNames. Also copy any TParameter with "POT" in the name
// and explicitly try to copy "POTUsed".
//
// outFile must already be open for writing (and will be the destination).
void CopyObjectsWithPrefix(const char* inputFilePath,
                                   TFile* outFile,
                                   const std::vector<std::string>& skipNames,
                                   const std::string& prefix = "DeltaPt_")
{
  TFile* in = TFile::Open(inputFilePath, "READ");
  if (!in || in->IsZombie()) {
    std::cerr << "ERROR: could not open input file: " << inputFilePath << std::endl;
    return;
  }

  // pointer to list of keys in the input file
  TList* keys = in->GetListOfKeys();
  if (!keys) {
    std::cerr << "ERROR: input file has no keys: " << inputFilePath << std::endl;
    in->Close();
    delete in;
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

    // Clean up the object created by ReadObj()
    delete obj;
  }

  // Try to copy POTUsed explicitly if present
  TKey* kPot = dynamic_cast<TKey*>(keys->FindObject("POTUsed"));
  if (kPot) {
    TObject* potObj = kPot->ReadObj();
    if (potObj) {
      outFile->cd();
      potObj->Write("POTUsed", TObject::kOverwrite);
      delete potObj;
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
        delete obj;
      } else {
        // If you want to copy everything with "POT" regardless of type, uncomment:
        // TObject* obj = key->ReadObj(); outFile->cd(); obj->Write(kname.c_str(), TObject::kOverwrite); delete obj;
      }
    }
  }

  in->Close();
  delete in;
}


int main(int argc, char** argv) {
#ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable();
#endif

  TH1::AddDirectory(kFALSE); // avoid ownership issues with ROOT directories

  if (argc < 3) {
    std::cerr << "USAGE: " << argv[0] << " <data.root> <mc.root> [rebinFactor=3] [variableName=DeltaPt]\n";
    return 1;
  }

  const char* dataPath = argv[1];
  const char* mcPath   = argv[2];
  int rebinFactor = 3;
  std::string varName = "DeltaPt";
  if (argc >= 4) rebinFactor = std::stoi(argv[3]);
  if (argc >= 5) varName = argv[4];

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

  PlotUtils::MnvPlotter plotter;
  plotter.legend_text_size = 0.015;
  plotter.data_line_width = 2;
  plotter.data_marker_size = 1.5;
  plotter.draw_normalized_to_bin_width = 0;

  // MC colors (same order used in Python)
  const std::vector<int> mcColors = {4, 7, 6, 2, 5, 416};
  int arr_int[6];
  for (size_t i=0;i<mcColors.size();++i) arr_int[i]=mcColors[i];
  int* arr = arr_int;

  // categories and sidebands exactly as in Python
  const std::vector<std::string> bkgdCategoryNames = {
    "selected_signal_reco", "background_NuECC_with_pions", "background_Other_NueCC",
    "background_NC_pi0", "background_CC_Numu_pi0", "background_Other"
  };
  const std::vector<std::string> sidebands = {"_", "_MeanFrontDEDXSB_", "_MichelSB_"};

  // Load histograms (MnvH1D) for data & MC
  std::vector<PlotUtils::MnvH1D*> data_hists;
  std::vector<std::vector<PlotUtils::MnvH1D*>> mc_hists(6); // 6 categories, each has 3 sidebands

  for (size_t s = 0; s < sidebands.size(); ++s) {
    std::string dataName = varName + sidebands[s] + "data";
    PlotUtils::MnvH1D* d = nullptr;
    dataFile->GetObject(dataName.c_str(), d);
    if (!d) {
      std::cerr << "ERROR: data hist " << dataName << " not found in data file\n";
      return 10;
    }
    // Make clones / rebin and detach from file
    PlotUtils::MnvH1D* dclone = dynamic_cast<PlotUtils::MnvH1D*>(d->Clone((d->GetName()+std::string("_clone")).c_str()));
    dclone->SetDirectory(nullptr);
    if (rebinFactor > 1) {
      // MnvH1D::Rebin returns TH1*, but calling Rebin on MnvH1D is supported and returns a TH1*;
      // to keep MnvH1D behaviour, we rebin the clone object in-place using Rebin(rebinFactor)
      dclone->Rebin(rebinFactor);
    }
    data_hists.push_back(dclone);

    for (size_t c = 0; c < bkgdCategoryNames.size(); ++c) {
      std::string mname = varName + sidebands[s] + bkgdCategoryNames[c];
      PlotUtils::MnvH1D* m = nullptr;
      mcFile->GetObject(mname.c_str(), m);
      if (!m) {
        std::cerr << "ERROR: mc histogram " << mname << " not found in mc file\n";
        return 11;
      }
      PlotUtils::MnvH1D* mclone = dynamic_cast<PlotUtils::MnvH1D*>(m->Clone((m->GetName())));
      mclone->SetDirectory(nullptr);
      if (rebinFactor > 1) mclone->Rebin(rebinFactor);
      mc_hists[c].push_back(mclone);
    }
  }

  // sanity: nbins same
  const int nbins = data_hists[0]->GetNbinsX();
  for (auto* h : data_hists) {
    if (h->GetNbinsX() != nbins) {
      std::cerr << "ERROR: inconsistent binning among data histos\n";
      return 12;
    }
  }

  double lambda = 1000;
  //ScaleFactors sfs = ExtractScaleFactors(data_hists, mc_hists, mcScale);
  ScaleFactors sfs = ExtractScaleFactorsWithReg(data_hists, mc_hists, mcScale, lambda);
  //ScaleFactors sfs;
  //sfs.meanFrontBkg.resize(nbins, 1.0);
  //sfs.meanFrontSig.resize(nbins, 1.0);
  //sfs.michelBkg.resize(nbins, 1.0);
  //sfs.michelSig.resize(nbins, 1.0);

  // Make scale factor histograms
  auto makeSFhist = [&](const std::string& name, const std::vector<double>& vals) -> PlotUtils::MnvH1D* {
    PlotUtils::MnvH1D* h = dynamic_cast<PlotUtils::MnvH1D*>(data_hists[0]->Clone(name.c_str()));
    h->Reset("ICES");
    for (int ib = 1; ib <= nbins; ++ib) h->SetBinContent(ib, vals[ib-1]);
    h->SetTitle(name.c_str());
    h->SetDirectory(nullptr);
    return h;
  };
  auto mean_bkg_h = makeSFhist("meanFrontBkg_hist", sfs.meanFrontBkg);
  auto mean_sig_h = makeSFhist("meanFrontSig_hist", sfs.meanFrontSig);
  auto mic_bkg_h  = makeSFhist("michelBkg_hist", sfs.michelBkg);
  auto mic_sig_h  = makeSFhist("michelSig_hist", sfs.michelSig);

  auto saveSFPlot = [&](PlotUtils::MnvH1D* hist, const std::string& filename) {
    TCanvas* c = new TCanvas(("c_" + filename).c_str(), "", 800, 600);
    c->SetGrid();
    
    hist->SetLineWidth(3);
    hist->SetLineColor(kBlue+1);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(1.0);
    hist->SetMarkerColor(kBlue+1);
    
    hist->GetYaxis()->SetTitle("Scale Factor");
    hist->Draw("E1");
    
    gPad->SetTicks();
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    
    c->SaveAs((filename + ".png").c_str());
    delete c;
  };

  // Save each scale factor histogram
  saveSFPlot(mean_bkg_h, "meanFront_bkg_scale_factors");
  saveSFPlot(mean_sig_h, "meanFront_sig_scale_factors");
  saveSFPlot(mic_bkg_h,  "michel_bkg_scale_factors");
  saveSFPlot(mic_sig_h,  "michel_sig_scale_factors");
 
  // ---------------------------
  // Create scaled histograms for two sidebands (for plotting correctness) and final signal region.
  // We'll follow the python application logic exactly.
  // ---------------------------

  // --- Helper to build scaled versions for a given sideband index s:
  // For s==1 (MeanFront): apply meanFront bkg SF to categories 3 & 4, and sig SF to cat 0
  // For s==2 (Michel): apply michel bkg SF to category 1, and sig SF to cat 0
  auto makeScaledForSideband = [&](size_t s) -> std::vector<PlotUtils::MnvH1D*> {
    std::vector<PlotUtils::MnvH1D*> scaled(bkgdCategoryNames.size(), nullptr);
    for (size_t c = 0; c < bkgdCategoryNames.size(); ++c) {
      PlotUtils::MnvH1D* base = mc_hists[c][s];
      PlotUtils::MnvH1D* clone = dynamic_cast<PlotUtils::MnvH1D*>(base->Clone((std::string(base->GetName())+std::string("_scaled")).c_str()));
      clone->SetDirectory(nullptr);
      clone->Reset("ICES");
      for (int ib = 1; ib <= nbins; ++ib) {
        double oldc = base->GetBinContent(ib);
        double olde = base->GetBinError(ib);
        double newc = oldc;
        double newe = olde;
        if (s == 1) { // MeanFront
          if (c == 3 || c == 4) {
            double sf = sfs.meanFrontBkg[ib-1];
            newc = oldc * sf;
            newe = olde * fabs(sf);
          }
          // signal (c==0) scaled by sig SF in the "sig+bg scaled" version below also; here we produce a single
          // scaled MC set where we also apply the signal SF so the plot will show the full scaled MC = (sig*sigSF + bkg*bkgSF + ...)
          if (c == 0) {
            double sf = sfs.meanFrontSig[ib-1];
            newc = oldc * sf;
            newe = olde * fabs(sf);
          }
        } else if (s == 2) { // Michel
          if (c == 1) {
            double sf = sfs.michelBkg[ib-1];
            newc = oldc * sf;
            newe = olde * fabs(sf);
          }
          if (c == 0) {
            double sf = sfs.michelSig[ib-1];
            newc = oldc * sf;
            newe = olde * fabs(sf);
          }
        }
        clone->SetBinContent(ib, newc);
        clone->SetBinError(ib, newe);
      }
      scaled[c] = clone;
    }
    return scaled;
  };

  // make meanFront scaled MC (signal+background scaled for that sideband)
  auto meanScaled = makeScaledForSideband(1);
  // make michel scaled MC (signal+background scaled for that sideband)
  auto michelScaled = makeScaledForSideband(2);

  // ---------------------------
  // Make the final signal-region scaled MC (what we save).
  // Python logic for the full/combined scaling (i==0 case in python):
  //   - For full scaled signal-region: scale category 1 (NueCCPion) with MICHEL bkg SF,
  //     scale categories 3 and 4 (NCPi0 and CCnumuPi0) with MEANFRONT bkg SF,
  //     leave signal (cat 0) UNCHANGED.
  // ---------------------------
  std::vector<PlotUtils::MnvH1D*> finalSignalScaled(bkgdCategoryNames.size(), nullptr);
  for (size_t c = 0; c < bkgdCategoryNames.size(); ++c) {
    PlotUtils::MnvH1D* base = mc_hists[c][0];
    PlotUtils::MnvH1D* clone = dynamic_cast<PlotUtils::MnvH1D*>(base->Clone((std::string(base->GetName())+std::string("_finalScaled")).c_str()));
    clone->SetDirectory(nullptr);
    clone->Reset("ICES");
    for (int ib = 1; ib <= nbins; ++ib) {
      double oldc = base->GetBinContent(ib);
      double olde = base->GetBinError(ib);
      double newc = oldc;
      double newe = olde;
      if (c == 1) { // category 1 scaled by MICHEL bkg SF
        double sf = sfs.michelBkg[ib-1];
        newc = oldc * sf;
        newe = olde * fabs(sf);
      } else if (c == 3 || c == 4) { // categories 3 & 4 scaled by MEANFRONT bkg SF
        double sf = sfs.meanFrontBkg[ib-1];
        newc = oldc * sf;
        newe = olde * fabs(sf);
      } else {
        // signal (c==0) and others left alone
        newc = oldc;
        newe = olde;
      }
      clone->SetBinContent(ib, newc);
      clone->SetBinError(ib, newe);
    }
    finalSignalScaled[c] = clone;
  }

  // ---------------------------
  // Plotting:
  //  * meanFront sideband plot: data_hists[1] vs meanScaled (stacked; colors in arr)
  //  * michel sideband plot:    data_hists[2] vs michelScaled
  //  * final signal region plot: data_hists[0] vs finalSignalScaled (signal unscaled, backgrounds scaled)
  // ---------------------------

  auto makeStackAndPlot = [&](PlotUtils::MnvH1D* data, const std::vector<PlotUtils::MnvH1D*>& mc_scaled,
                             const std::string& outName, const std::string& titleSuffix) {
    // set titles (use python's ordering & labels)
    std::vector<std::string> labels = {
      "signal (nu_e QELike + proton)",
      "nu_e nonQE (has FS mesons)",
      "Other nu_eCC",
      "NC with pi0",
      "nu_mu CC with pi0",
      "other"
    };
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
    for (int k = (int)mc_scaled.size()-1; k >= 0; --k) array.Add(mc_scaled[k]);

    TCanvas c("c", titleSuffix.c_str(), 1200, 900);
    //plotter.DrawDataStackedMC(data_hists[s], &arr, nullptr, mcScale, "TR", "Data", 1001, data_hists[0]->GetXaxis()->GetTitle(), "N events");
    plotter.DrawDataStackedMC(data, &array, arr, mcScale, "TR", "Data", 1001, data->GetXaxis()->GetTitle(), "N events");
    std::string normalization = std::string("Normalized to ") + std::to_string(dataPOT) + " data POT";
    plotter.AddPlotLabel(normalization.c_str(), 0.2, 0.95, 0.03);
    c.SaveAs(outName.c_str());
  };

  // meanFront SB (index 1)
  makeStackAndPlot(data_hists[1], meanScaled, (varName + "_MeanFrontSB_scaled.png"), "MeanFront SB scaled");

  // michel SB (index 2)
  makeStackAndPlot(data_hists[2], michelScaled, (varName + "_MichelSB_scaled.png"), "Michel SB scaled");

  // final signal region: data_hists[0] vs finalSignalScaled
  makeStackAndPlot(data_hists[0], finalSignalScaled, (varName + "_SignalRegion_finalScaled.png"), "Signal region final scaled");

  // ---------------------------
  // Write output root file containing the signal-region scaled MC histograms with original names.
  // The user requested these to have the same names as originals for downstream usage.
  // ---------------------------
  TFile* outFile = TFile::Open("scaled_mc.root", "RECREATE");
  if (!outFile || outFile->IsZombie()) {
    std::cerr << "ERROR: couldn't open scaled_mc.root for writing\n";
  } else {
    // for each category, write the scaled histogram but give it the original name (mc_hists[c][0]->GetName())
    std::vector<std::string> modifiedNames = {
      "DeltaPt_selected_signal_reco",
      "DeltaPt_background_NuECC_with_pions",
      "DeltaPt_background_Other_NueCC",
      "DeltaPt_background_NC_pi0",
      "DeltaPt_background_CC_Numu_pi0",
      "DeltaPt_background_Other"
      // ... add any other names you've modified/written ...
    };
    // Copy all other objects from the mc input with prefix "DeltaPt_"
    // (preserving key names), and copy POTUsed/TParameter(POT*) as well.
    CopyObjectsWithPrefix(mcPath, outFile, modifiedNames, "DeltaPt_");

    for (size_t c = 0; c < bkgdCategoryNames.size(); ++c) {
      // clone finalSignalScaled[c] and set name to the original MC histogram name for signal-region
      std::string origName = mc_hists[c][0]->GetName(); // original name in MC file
      PlotUtils::MnvH1D* toWrite = dynamic_cast<PlotUtils::MnvH1D*>(finalSignalScaled[c]->Clone((origName + std::string("_out")).c_str()));
      toWrite->SetDirectory(outFile);
      toWrite->SetName(origName.c_str()); // overwrite name to original
      outFile->cd();
      toWrite->Write();
      delete toWrite;
    }
    
    // flush and close
    outFile->Write();
    outFile->Close();
    delete outFile;
    //outFile->Close();
    std::cout << "Wrote scaled MC signal-region histograms to scaled_mc.root\n";
  }

  // clean up
  for (auto* h : data_hists) delete h;
  for (size_t c=0;c<mc_hists.size();++c) {
    for (auto* h : mc_hists[c]) delete h;
  }
  for (auto* p : meanScaled) delete p;
  for (auto* p : michelScaled) delete p;
  for (auto* p : finalSignalScaled) delete p;
  delete mean_bkg_h;
  delete mean_sig_h;
  delete mic_bkg_h;
  delete mic_sig_h;

  std::cout << "Done.\n";
  return 0;
}
