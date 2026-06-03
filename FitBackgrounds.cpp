// FitBackgrounds.cpp
// Takes as input root files from runEventLoop containing signal samples and sideband samples of meanfront dE/dX and >= 1 michel events.
// Calculates scale factors (with multiple options how), creates plots of scale factors, scaled sidebands, scaled signal
// And saves the scaled signal region to a new root file, to be used in ExtractCrossSection.cpp
//
// Usage (example):
//   FitBackgrounds data.root mc.root 1000 Lepton_Pt
//   input data histos, input mc histos, regularization strength lambda, variable name

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
#include "TColor.h"
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
#include <TMath.h>
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

//c++ includes
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

//For testing, forces all output scale factors to 1. DON'T LEAVE THIS ON BY ACCIDENT
const bool set_scale_factors_to_1 = false;
double mcScale = 1.0;  // will be set in main
std::string varName = "DeltaPt"; //default value that gets overwritten in main

//const std::vector<std::string> bkgdCategoryNames = {"selected_signal_reco", "background_NuECC_with_pions", "background_Other_NueCC", "background_NC_pi0", "background_CC_Numu_pi0", "background_Other"};
const std::vector<std::string> bkgdCategoryNames = {"selected_signal_reco", "background_NuECC_nonQELike_single_pi_plus", "background_NuECC_nonQELike_single_pi_zero", "background_NuECC_nonQELike_single_pi_minus", "background_NuECC_nonQELike_Npi", "background_Other_NueCC", "background_NC_pi0", "background_CC_Numu_pi0", "background_Other"};
const std::vector<std::string> sidebands = {"_", "_MeanFrontDEDXSB_", "_MichelSB_"};

//indices of the different background categories, in the past I had only 1 nonQELike bkg category
// so if I wanna fit those histograms then I need to use nonQELikeIdx = {1}, pi0Idx = {3, 4}, and fixedIdx = {0, 2, 5}
//const std::vector<int> nonQELikeIdx = {1};
//const std::vector<int> pi0Idx = {3, 4};
//const std::vector<int> fixedIdx = {0,2,5}; // signal + otherNueCC + other
//const std::vector<int> fixedIdxNoSignal = {2,5}; // otherNueCC + other, for the separate fits where I also calc signal scale factors
const std::vector<int> nonQELikeIdx = {1,2,3,4};
const std::vector<int> pi0Idx = {6,7};
const std::vector<int> fixedIdx = {0,5,8}; // signal + otherNueCC + other
const std::vector<int> fixedIdxNoSignal = {5,8}; // otherNueCC + other, for the separate fits where I also calc signal scale factors

//Struct to contain the scale factors for each background in each universe, saved as an MnvH1D
struct ScaleFactors {
  PlotUtils::MnvH1D* meanFrontBkg_mnvhist;
  PlotUtils::MnvH1D* meanFrontSig_mnvhist;
  PlotUtils::MnvH1D* michelBkg_mnvhist;
  PlotUtils::MnvH1D* michelSig_mnvhist;

  void Init(const PlotUtils::MnvH1D* mc_hist) {
    auto clone = [&](const std::string& suffix) {
      auto* h = dynamic_cast<PlotUtils::MnvH1D*>(mc_hist->Clone((varName + suffix).c_str()));
      h->Reset("ICES");
      return h;
    };
    meanFrontBkg_mnvhist = clone("_meanFrontBkg");
    meanFrontSig_mnvhist = clone("_meanFrontSig");
    michelBkg_mnvhist    = clone("_michelBkg");
    michelSig_mnvhist    = clone("_michelSig");
  }
  
  // sets one bin on either CV or a specific universe
  void SetBin(PlotUtils::MnvH1D* hist, int ib, double val,
	      bool isCV, const std::string& bandName, int universe_index) {
    if (isCV) hist->SetBinContent(ib, val);
    else hist->GetVertErrorBand(bandName)->GetHist(universe_index)->SetBinContent(ib, val);
  }
  
  // vals_per_bin: function that maps bin index (1-based) to {meanFrontBkg, meanFrontSig, michelBkg, michelSig}
  void WriteOutput(int nbins, bool isCV, const std::string& bandName, int universe_index,
		   std::function<std::array<double,4>(int)> vals_per_bin) {
    for (int ib = 1; ib <= nbins; ib++) {
      auto [mfBkg, mfSig, miBkg, miSig] = vals_per_bin(ib);
      SetBin(meanFrontBkg_mnvhist, ib, mfBkg, isCV, bandName, universe_index);
      SetBin(meanFrontSig_mnvhist, ib, mfSig, isCV, bandName, universe_index);
      SetBin(michelBkg_mnvhist,    ib, miBkg, isCV, bandName, universe_index);
      SetBin(michelSig_mnvhist,    ib, miSig, isCV, bandName, universe_index);
    }
  }  

};

struct SidebandData {
  // bin contents for each region, indexed by [category][region][bin]
  // regions: 0=signal, 1=meanFront, 2=michel
  // categories: same as mc_hists indices
  std::vector<std::vector<std::vector<double>>> mc;  
  std::vector<std::vector<double>> data;  // [region][bin]
  bool isCV;
  int nbins;
  
  //loop through provided histograms and fill bin contents (for a single universe)
  void fillSidebandData(const std::vector<PlotUtils::MnvH1D*>& data_hists, const std::vector<std::vector<PlotUtils::MnvH1D*>>& mc_hists, std::string bandName, int universe_index) {

    //Helper for getting bin content from CV vs one of the error bands
    auto getBin = [&](PlotUtils::MnvH1D* h, int ib) -> double {
      if (isCV) return h->GetBinContent(ib) * mcScale;
      return h->GetVertErrorBand(bandName)->GetHist(universe_index)->GetBinContent(ib) * mcScale;
    };
    
    nbins = data_hists[0]->GetNbinsX();
    int nCategories = mc_hists.size();    // 9
    int nRegions    = mc_hists[0].size(); // 3
    mc.assign(nCategories, std::vector<std::vector<double>>(nRegions, std::vector<double>(nbins + 1, 0.0)));
    
    for (int cat = 0; cat < nCategories; cat++) {
      for (int reg = 0; reg < nRegions; reg++) {
	for (int ib = 1; ib <= nbins; ib++) {
	  mc[cat][reg][ib] = getBin(mc_hists[cat][reg], ib);
	}
      }
    }
    data.assign(nRegions, std::vector<double>(nbins + 1, 0.0));
    for (int reg = 0; reg < nRegions; reg++) {
      for (int ib = 1; ib <= nbins; ib++) {
	data[reg][ib] = data_hists[reg]->GetBinContent(ib);
      }
    }    
  }

};

//Helper to sum bins of the bkg categories being considered, as in, combining all the pi0 bkgs or the NonQELike bkgs
double sumBins(std::vector<int> indices, int region, int ib, std::vector<std::vector<std::vector<double>>> mc){
  double sum = 0;
  for (int idx : indices) sum += mc[idx][region][ib];
  return sum;
 };	

//calculates a chi^2 value between the data (obsreved) and the MC (predicted), using a single overall scale factor per bkg
// this function is what gets passed to Minuit2 and minimized, using the scale factors as parameters
double Chi2NormOnly(double alpha_nonQELike, double alpha_pi0, 
                    const SidebandData& sb) {
    double chi2 = 0;
    
    for (int region = 0; region < 3; region++) {
        for (int ib = 1; ib <= sb.nbins; ib++) {
            // sum fixed backgrounds (assume they're correct, subtract from data without any scaling)
            double fixed = 0;
            for (int idx : fixedIdx) fixed += sb.mc[idx][region][ib];
            
            // predicted = scaled nonQELike + scaled pi0 + fixed
            double nonQELike_pred = alpha_nonQELike * sumBins(nonQELikeIdx, region, ib, sb.mc);
            double pi0_pred       = alpha_pi0       * sumBins(pi0Idx,       region, ib, sb.mc);
            double predicted = nonQELike_pred + pi0_pred + fixed;
            
            double observed = sb.data[region][ib];
            double sigma2 = observed > 0 ? observed : 1.0;  // Poisson uncertainty: sigma^2 = N
            
            chi2 += pow(observed - predicted, 2) / sigma2;
        }
    }
    return chi2;
}

ScaleFactors ExtractScaleFactors_chi2_simultaneous_normOnly(
    const std::vector<PlotUtils::MnvH1D*>& data_hists,
    const std::vector<std::vector<PlotUtils::MnvH1D*>>& mc_hists)
{
  ScaleFactors sf;  
  sf.Init( mc_hists[0][0] );

  std::vector<std::string> vertErrorBandNames = mc_hists[0][0]->GetVertErrorBandNames();
  vertErrorBandNames.push_back("cv"); //add the cv to the list since it doesn't get returned, this way don't have to write a separate loop
  
  for (const auto& bandName : vertErrorBandNames){//Loop over error bands
    bool isCV = (bandName == "cv");
    PlotUtils::MnvVertErrorBand* band;
    int nHists;
    if (isCV) { nHists = 1; }
    else {
      band = mc_hists[0][0]->GetVertErrorBand( bandName );
      nHists = band->GetNHists();
    }
    for (int universe_index = 0; universe_index < nHists; universe_index++){ //Loop over universes within that error band (normally 2, although flux has 100)

      SidebandData sb;
      sb.isCV = (bandName == "cv");
      sb.fillSidebandData(data_hists, mc_hists, bandName, universe_index);
      
      ROOT::Math::Functor fcn([=](const double* p) { return Chi2NormOnly(p[0], p[1], sb); }, 2);      

      std::unique_ptr<ROOT::Math::Minimizer> min(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
      min->SetFunction(fcn);
      min->SetVariable(0, "alpha_nonQELike", 1.0, 0.01);  // start at 1, step 0.01
      min->SetVariable(1, "alpha_pi0",       1.0, 0.01);
      min->SetVariableLimits(0, 0.0, 5.0);  // physical constraint: non-negative
      min->SetVariableLimits(1, 0.0, 5.0);
      min->Minimize();

      double alpha_nonQELike = min->X()[0];
      double alpha_pi0       = min->X()[1];
      double err_nonQELike   = min->Errors()[0];  // from Hessian
      double err_pi0         = min->Errors()[1];  // these are statistical only uncertainties on the minuit fit itself (per universe), not sure yet if I'll use em for anything

      sf.WriteOutput(sb.nbins, isCV, bandName, universe_index, [&](int ib) {  return std::array<double,4>{ alpha_pi0, 1.0, alpha_nonQELike, 1.0 }; });
    } //end universe loop within error band
  } //end error band loop
  return sf;
}

double Chi2NormOnly_separate(double alpha, double const_scale_factor, const SidebandData& sb, bool forMeanFront) {
  double chi2 = 0;  
  //since this method does the two fits separately, each separate fit considers only signal region + relevant sideband, loop thru twice
  for (int i = 0; i < 2; i++) {
    int region;
    if (i==0) region = 0; //signal region is always the first one
    else{
      region = forMeanFront ? 1 : 2; //meanfront index = 1, michel index = 2
    }
    for (int ib = 1; ib <= sb.nbins; ib++) {
      // sum fixed backgrounds (taken as correct, subtracted from data)
      double fixed = sumBins(fixedIdx, region, ib, sb.mc);

      double nonQELike_pred = 0;
      double pi0_pred       = 0;
      // predicted = scaled nonQELike + scaled pi0 + fixed
      // alpha is the varying scale factor, const_scale_factor is the constant scale factor applied to the bkg not under consideration
      if (forMeanFront){ 
	nonQELike_pred = const_scale_factor * sumBins(nonQELikeIdx, region, ib, sb.mc);
	pi0_pred       = alpha * sumBins(pi0Idx, region, ib, sb.mc);
      } else { 
	nonQELike_pred = alpha * sumBins(nonQELikeIdx, region, ib, sb.mc);
	pi0_pred       = const_scale_factor * sumBins(pi0Idx, region, ib, sb.mc);
      }
      double predicted = nonQELike_pred + pi0_pred + fixed;
      
      double observed = sb.data[region][ib];
      double sigma2 = observed > 0 ? observed : 1.0;  // Poisson: sigma^2 = N
      
      chi2 += pow(observed - predicted, 2) / sigma2;
    }
  }
  return chi2;
}

ScaleFactors ExtractScaleFactors_chi2_separate_normOnly( const std::vector<PlotUtils::MnvH1D*>& data_hists, const std::vector<std::vector<PlotUtils::MnvH1D*>>& mc_hists){
  ScaleFactors sf;  
  sf.Init( mc_hists[0][0] );

  std::vector<std::string> vertErrorBandNames = mc_hists[0][0]->GetVertErrorBandNames();
  vertErrorBandNames.push_back("cv"); //add the cv to the list since it doesn't get returned, this way don't have to write a separate loop
  
  for (const auto& bandName : vertErrorBandNames){//Loop over error bands
    bool isCV = (bandName == "cv");
    PlotUtils::MnvVertErrorBand* band;
    int nHists;
    if (isCV) { nHists = 1; }
    else {
      band = mc_hists[0][0]->GetVertErrorBand( bandName );
      nHists = band->GetNHists();
    }
    for (int universe_index = 0; universe_index < nHists; universe_index++){ //Loop over universes within that error band (normally 2, although flux has 100)

      SidebandData sb;
      sb.isCV = (bandName == "cv");
      sb.fillSidebandData(data_hists, mc_hists, bandName, universe_index);

      double chi2_MF, chi2_Michel;
      double alpha_nonQELike = 1;
      double alpha_pi0 = 1;
      //run 2 fits, one for pi0 scale factor (using meanfront sb + signal region), and one for nonQELike scale factor (using michel + signal region)
      for (int i=1; i<31; i++){
	bool forMeanFront = static_cast<bool>(i%2);
	double const_scale_factor;
	if (forMeanFront) const_scale_factor = alpha_nonQELike;
	else const_scale_factor = alpha_pi0;
	
	ROOT::Math::Functor fcn([=](const double* p) { return Chi2NormOnly_separate(p[0], const_scale_factor, sb, forMeanFront); }, 1);
	
	std::unique_ptr<ROOT::Math::Minimizer> min(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
	min->SetFunction(fcn);
	min->SetVariable(0, "alpha", 1.0, 0.01);  // start at 1, step 0.01
	min->SetVariableLimits(0, 0.0, 5.0);  // physical constraint: non-negative
	min->Minimize();
	
	double alpha = min->X()[0];
	double err = min->Errors()[0];  // from Hessian
	double chi2 = min->MinValue();
	if (forMeanFront) {
	  alpha_pi0 = alpha;
	  chi2_MF = min->MinValue();
	} else {
	  alpha_nonQELike = alpha;
	  chi2_Michel = min->MinValue();
	}
      }
      sf.WriteOutput(sb.nbins, isCV, bandName, universe_index, [&](int ib) {  return std::array<double,4>{ alpha_pi0, 1.0, alpha_nonQELike, 1.0 }; });
    } //end universe loop within error band
  } //end error band loop
  return sf;
}

ScaleFactors ExtractScaleFactors_simultaneous_binByBin(
    const std::vector<PlotUtils::MnvH1D*>& data_hists,
    const std::vector<std::vector<PlotUtils::MnvH1D*>>& mc_hists,
    long long lambda)
{
  ScaleFactors sf;  
  sf.Init(mc_hists[0][0]);
  
  std::vector<std::string> vertErrorBandNames = mc_hists[0][0]->GetVertErrorBandNames();
  vertErrorBandNames.push_back("cv"); //add the cv to the list since it doesn't get returned, this way don't have to write a separate loop
  
  for (const auto& bandName : vertErrorBandNames){//Loop over error bands
    bool isCV = (bandName == "cv");
    PlotUtils::MnvVertErrorBand* band;
    int nHists;
    if (isCV) { nHists = 1; }
    else {
      band = mc_hists[0][0]->GetVertErrorBand( bandName );
      nHists = band->GetNHists();
    }
    for (int universe_index = 0; universe_index < nHists; universe_index++){ //Loop over universes within that error band (normally 2, although flux has 100)

      SidebandData sb;
      sb.isCV = isCV;
      sb.fillSidebandData(data_hists, mc_hists, bandName, universe_index);
      
      // Build the chi2-weighted system Ax = d
      // unknowns x = [alpha_nonQELike_1, alpha_pi0_1, alpha_nonQELike_2, alpha_pi0_2, ...]
      int nUnknowns = 2 * sb.nbins;
      int nEquations = 3 * sb.nbins;  // 3 regions * nbins, all simultaneously
      int nRegRows = static_cast<int>(lambda > 0 ? (2*(sb.nbins-1)) : 0); //regularization rows, need nbins-1 rows to link neighboring bins of each scale factor, and we have 2 sets of scale factors
      int nRows = nEquations + nRegRows;
      TMatrixD A(nRows, nUnknowns);
      TVectorD d(nRows);
      //2 bins, region = 1, ib = 1
      for (int region = 0; region < 3; region++) {
	for (int ib = 1; ib <= sb.nbins; ib++) {
	  int row = region * sb.nbins + (ib-1);
	  int col_nonQELike = (ib-1) * 2 + 0;
	  int col_pi0       = (ib-1) * 2 + 1;
	  
	  double observed = sb.data[region][ib];
	  double sigma = observed > 0 ? sqrt(observed) : 1.0;
	  
	  //get my mc event counts by fitting category (either fixed, takes a pi0 scale factor, or takes a nonQELike scale factor)
	  double fixed        = sumBins(fixedIdx,     region, ib, sb.mc);
	  double nonQELike_mc = sumBins(nonQELikeIdx, region, ib, sb.mc);
	  double pi0_mc       = sumBins(pi0Idx,       region, ib, sb.mc);
	  
	  // divide everything by sigma for chi2 weighting
	  A(row, col_nonQELike) = nonQELike_mc / sigma;
	  A(row, col_pi0)       = pi0_mc / sigma;
	  d(row)                = (observed - fixed) / sigma;
	}
      }

      // add smoothness/regularization rows which couple neighboring scale factors together and penalize differences
      if (lambda > 0) {
	int regRowStart = nEquations;
	int r = 0;
	for (int ib = 1; ib <= sb.nbins - 1; ++ib) {
	  int row_nonQELike = regRowStart + r++;
	  int row_pi0       = regRowStart + r++;
	  double w = std::sqrt(lambda);

	  //nonQELike regularization rows
	  int col_nonQELike_pos = (ib - 1) * 2;
	  int col_nonQELike_neg = (ib    ) * 2;
	  A(row_nonQELike, col_nonQELike_pos) =  w;
	  A(row_nonQELike, col_nonQELike_neg) = -w;
	  d(row_nonQELike) = 0.0;

	  //pi0 regularization rows
	  int col_pi0_pos = (ib - 1) * 2 + 1;
	  int col_pi0_neg = (ib    ) * 2 + 1;
	  A(row_pi0, col_pi0_pos) =  w;
	  A(row_pi0, col_pi0_neg) = -w;
	  d(row_pi0) = 0.0;
	  
	}
      }
      
      TDecompSVD svd(A);
      //svd.SetTol(1e-3);
      Bool_t ok;
      TVectorD rhs = d;
      TVectorD x = svd.Solve(rhs, ok);
      if (isCV) {
	std::cout << "solution vector x = [";
	for (int j=0; j<x.GetNrows(); j++){
	  if (j==x.GetNrows()-1) std::cout << x(j) << "]" << std::endl;
	  else std::cout << x(j) << ", ";
	}
      }
      sf.WriteOutput(sb.nbins, isCV, bandName, universe_index, [&](int ib) {
	return std::array<double,4>{ x(2*(ib-1)+1), 1.0, x(2*(ib-1)+0), 1.0 };
      });

    } //end universe loop within error band
  } //end error band loop
  return sf;
}

// data_hists: sidebands in order [signalRegion, meanFrontSB, michelSB] (same as your code)
// mc_hists: vector per category: mc_hists[cat][sideband_index]
// lambda: regularization strength (0 => no regularization -> exact per-bin solution)
// regularizeSignal: whether to regularize signal scale-factors too (default false)
ScaleFactors ExtractScaleFactors_separate_binByBin(
    const std::vector<PlotUtils::MnvH1D*>& data_hists,
    const std::vector<std::vector<PlotUtils::MnvH1D*>>& mc_hists,
    long long lambda)
{
  ScaleFactors sf; //create and prep MnvH1Ds for holding scale factors, unique set per sideband and per universe
  sf.Init( mc_hists[0][0] );
  
  std::vector<std::string> vertErrorBandNames = mc_hists[0][0]->GetVertErrorBandNames();
  vertErrorBandNames.push_back("cv"); //add the cv to the list since it doesn't get returned, this way don't have to write a separate loop
    
  for (const auto& bandName : vertErrorBandNames){//Loop over error bands
    bool isCV = (bandName == "cv");
    PlotUtils::MnvVertErrorBand* band;
    int nHists;
    if (isCV) { nHists = 1; }
    else {
      band = mc_hists[0][0]->GetVertErrorBand( bandName );
      nHists = band->GetNHists();
    }
    for (int universe_index = 0; universe_index < nHists; universe_index++){ //Loop over universes within that error band (normally 2, although flux has 100)
      
      SidebandData sb;
      sb.isCV = isCV;
      sb.fillSidebandData(data_hists, mc_hists, bandName, universe_index);

      int nRows = 2*sb.nbins;
      int nUnknowns = 1*sb.nbins;
      int nRegRows = static_cast<int>(lambda > 0 ? (sb.nbins-1) : 0); 
      int region;
      std::vector<double> alpha_nonQELike(sb.nbins, 1.0);
      std::vector<double> alpha_pi0(sb.nbins, 1.0);
      for (int i=1; i<31; i++){ //number of back and forth iterations...
	bool forMeanFront = static_cast<bool>(i%2); //want this to flip flop every iteration...
	region = forMeanFront ? 1 : 2;

	TMatrixD A(nRows, nUnknowns); // zero-initialized
	TVectorD d(nRows);            // RHS
	for (int ib = 1; ib <= sb.nbins; ++ib) {
	  double predicted_SR; //mc prediction of ONLY the relevant mc contribution in signal region
	  double predicted_SB; //same thing in sideband region (these things are pi0 bkgs and mean front dedx respectively, or nonQELike bkgs and michel SB respectively)
	  double observed_SR = sb.data[0][ib];
	  double observed_SB = sb.data[region][ib];
	  double sigma_SR = observed_SR > 0 ? sqrt(observed_SR) : 1.0;
	  double sigma_SB = observed_SR > 0 ? sqrt(observed_SR) : 1.0;

          //get my mc event counts by fitting category (either fixed, takes a pi0 scale factor, or takes a nonQELike scale factor)
          double fixed_SR     = sumBins(fixedIdx,     0, ib, sb.mc);
          double nonQELike_SR = sumBins(nonQELikeIdx, 0, ib, sb.mc);
          double pi0_SR       = sumBins(pi0Idx,       0, ib, sb.mc);

          double fixed_SB     = sumBins(fixedIdx,     region, ib, sb.mc);
          double nonQELike_SB = sumBins(nonQELikeIdx, region, ib, sb.mc);
          double pi0_SB       = sumBins(pi0Idx,       region, ib, sb.mc);

	  if (forMeanFront) {
	    predicted_SR = pi0_SR;
	    predicted_SB = pi0_SB;
	    observed_SR = observed_SR - (alpha_nonQELike[ib-1]*nonQELike_SR + fixed_SR);
	    observed_SB = observed_SB - (alpha_nonQELike[ib-1]*nonQELike_SB + fixed_SB);
	  } else {
	    predicted_SR = nonQELike_SR;
	    predicted_SB = nonQELike_SB;
	    observed_SR = observed_SR - (alpha_pi0[ib-1]*pi0_SR + fixed_SR);
	    observed_SB = observed_SB - (alpha_pi0[ib-1]*pi0_SB + fixed_SB);
	  }

	  int row_SR = (ib-1)*2;
	  int row_SB = (ib-1)*2 + 1;
	  A(row_SR, ib-1) = predicted_SR / sigma_SR;
	  A(row_SB, ib-1) = predicted_SB / sigma_SB ;
	  d(row_SR) = observed_SR / sigma_SR;
	  d(row_SB) = observed_SB / sigma_SB;
	  
	} // end bin loop for data eqs
	
	// Minimize chi2: ||A x - d|| using SVD
	TDecompSVD svd(A);
	svd.SetTol(1e-3);
	Bool_t ok;
	TVectorD rhs = d;            // copy RHS because Solve modifies it
	TVectorD x = svd.Solve(rhs, ok); // vector of length nUnknowns: [b1,s1,b2,s2,...]
	if (forMeanFront) {
	  for (int jb=0; jb < sb.nbins; jb++){
	    alpha_pi0[jb] = x(jb);
	  }
	} else {
	  for (int jb=0; jb < sb.nbins; jb++){
	    alpha_nonQELike[jb] = x(jb);
	  }
	}
      }
      sf.WriteOutput(sb.nbins, isCV, bandName, universe_index, [&](int ib) {
	return std::array<double,4>{ alpha_pi0[ib-1], 1.0, alpha_nonQELike[ib-1], 1.0 };
      });

      /*
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
	*/
    }  //end loop over universes within an error band
  }//end loop over error bands
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
  // set titles (use my python function's ordering & labels)

  std::vector<std::string> labels;
  std::vector<int> mcColors;
  if (bkgdCategoryNames.size() == 6) {
    labels = {"signal (nu_e QELike + proton)", "nu_e nonQE (has FS mesons)", "Other nu_eCC", "NC with pi0", "nu_mu CC with pi0", "other"};
    mcColors = {4, 7, 6, 2, 5, 416};
  } else if (bkgdCategoryNames.size() == 9) {
    labels = { "Signal", "Single #pi^{+}", "Single #pi^{-}", "Single #pi^{0}", "N#pi", "Other #nu_{e}CC", "NC with #pi^{0}", "#nu_{#mu}CC with #pi^{0}", "other"};
    mcColors = { TColor::GetColor("#0000FF"), TColor::GetColor("#00FFFF"), TColor::GetColor("#FF00FF"), TColor::GetColor("#FF0000"), TColor::GetColor("#DAA520"), TColor::GetColor("#FFD700"), TColor::GetColor("#FFFF00"), TColor::GetColor("#FFFACD"), kGreen};
  }
  
  PlotUtils::MnvPlotter plotter;
  plotter.legend_text_size = 0.015;
  plotter.data_line_width = 2;
  plotter.data_marker_size = 1.5;
  
  int arr_int[bkgdCategoryNames.size()];
  for (size_t i=0;i<mcColors.size();++i) arr_int[i]=mcColors[i];
  int* arr = arr_int;
  
  for (size_t c=0;c<mc_scaled.size();++c) {
    mc_scaled[c]->SetTitle(labels[c].c_str());
    mc_scaled[c]->SetLineColor(kBlack);
    mc_scaled[c]->SetFillColor(mcColors[c]);
    mc_scaled[c]->SetLineWidth(1);
  }
  plotter.mc_line_width = 2;
  data->SetTitle("data");
  // Create TObjArray in reverse order so the stack looks like the python version (signal on top)
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

    // only consider keys with the requested prefix
    if (kname.find(prefix) == std::string::npos) continue;

    // skip names explicitly listed in skipNames
    // these are the ones we scaled, modified, and saved after fitting. don't wanna copy the old versions over to the new file
    if (std::find(skipNames.begin(), skipNames.end(), kname) != skipNames.end()) {
      // std::cout << "Skipping modified object: " << kname << std::endl;
      continue;
    }

    // read the object (allocates an object)
    TObject* obj = key->ReadObj();
    if (!obj) {
      std::cerr << "Warning: cannot read object " << kname << " from " << inputFilePath << std::endl;
      continue;
    }

    // write into output file with same key name; overwrite it if present
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
  //int method = std::stod(argv[3]);
  long long lambda = std::stod(argv[3]);
  std::cout << "lambda = " << lambda << std::endl;
  int method = 3;
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

  // get POT scaling
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
  mcScale = dataPOT / mcPOT;
  std::cout << "mc POT scale = " << mcScale << "  (dataPOT=" << dataPOT << ", mcPOT=" << mcPOT << ")\n";


  // load histograms (MnvH1D) for data & MC
  std::vector<PlotUtils::MnvH1D*> data_hists;
  std::vector<std::vector<PlotUtils::MnvH1D*>> mc_hists( bkgdCategoryNames.size() );
  
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

  std::cout << "Now trying to extract scale factors..." << std::endl;
  ScaleFactors sfs;
  std::string methodName;
  if (method==0) {
    methodName = "simultaneous_norm_only";
    sfs = ExtractScaleFactors_chi2_simultaneous_normOnly(data_hists, mc_hists);
  }
  else if (method==1) {
    methodName = "iterated_norm_only";
    sfs = ExtractScaleFactors_chi2_separate_normOnly(data_hists, mc_hists);
  }
  else if (method==2) {
    methodName = "iterated_bin_by_bin";    
    sfs = ExtractScaleFactors_separate_binByBin(data_hists, mc_hists, lambda);
  }
  else if (method==3) {
    methodName = "simultaneous_bin_by_bin";    
    sfs = ExtractScaleFactors_simultaneous_binByBin(data_hists, mc_hists, lambda);
  }
  std::cout << "Succeeded" << std::endl; 

  int nbins = data_hists[0]->GetNbinsX();
  /*
  std::cout << "pi0 scale factors      : [";
  for (int i=0; i<nbins; i++){
    if (i==(nbins-1)) std::cout << sfs.meanFrontBkg_mnvhist->GetBinContent(i+1) << "]" << std::endl;
    else std::cout << sfs.meanFrontBkg_mnvhist->GetBinContent(i+1) << ", ";
  }
  std::cout << "nonQELike scale factors: [";
  for (int i=0; i<nbins; i++){
    if (i==(nbins-1)) std::cout << sfs.michelBkg_mnvhist->GetBinContent(i+1) << "]" << std::endl;
    else std::cout << sfs.michelBkg_mnvhist->GetBinContent(i+1) << ", ";
  }
  */
  
  // Save each scale factor histogram
  saveSFPlot(sfs.meanFrontBkg_mnvhist, varName + "_" + methodName + "_meanFront_bkg_scale_factors");
  //saveSFPlot(sfs.meanFrontSig_mnvhist, varName + methodName + "meanFront_sig_scale_factors");
  saveSFPlot(sfs.michelBkg_mnvhist, varName + "_" + methodName + "_michel_bkg_scale_factors");
  //saveSFPlot(sfs.michelSig_mnvhist, varName + methodName + "michel_sig_scale_factors");
 
  // --- lil lambda function to build scaled versions of my distributions for a given sideband index s:
  // For s==0 (signal region): leave signal unscaled, apply michel SF to nonQELike bkgs, and apply meanFront SF to pi0 bkgs
  // For s==1 (MeanFront): apply meanFront bkg SF to pi0 bkgs, and sig SF to signal if applicable
  // For s==2 (Michel): apply michel bkg SF to nonQELike bkgs, and sig SF to signal if applicable
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
      if (std::find(nonQELikeIdx.begin(), nonQELikeIdx.end(), c) != nonQELikeIdx.end()) { //NonQELike yellow categories, scaled by michel scale factors
	clone->Multiply(clone, sfs.michelBkg_mnvhist);
      }
      if (std::find(pi0Idx.begin(), pi0Idx.end(), c) != pi0Idx.end()) { //NC Pi0 (purple) and NumuCC Pi0 (teal), scale by mean front dE/dX scale factors
	clone->Multiply(clone, sfs.meanFrontBkg_mnvhist);
      }

      scaled[c] = clone;
    }
    return scaled;
  };

  //apply scale factors returns a vector of mnvh1ds, one for each category
  auto finalSignalScaled = applyScaleFactors(0); 
  auto meanScaled = applyScaleFactors(1); 
  auto michelScaled = applyScaleFactors(2); 

  // meanFront SB plot
  saveStackPlot(data_hists[1], meanScaled, (varName + "_MeanFrontSB_scaled.png"), "MeanFront SB scaled", dataPOT, mcPOT);
  // michel SB plot 
  saveStackPlot(data_hists[2], michelScaled, (varName + "_MichelSB_scaled.png"), "Michel SB scaled", dataPOT, mcPOT);
  // final signal region
  saveStackPlot(data_hists[0], finalSignalScaled, (varName + "_SignalRegion_finalScaled.png"), "Signal region final scaled", dataPOT, mcPOT);

  //Now calculate a chi2 for each scaled region (CV only), and a total chi2 which is the sum of all of them.
  double chi2_SR = 0;
  double chi2_MF = 0;
  double chi2_Michel = 0;
  for (int region = 0; region<3; region++){
    for (int ib=1; ib <= nbins; ib++){
      double observed = data_hists[region]->GetBinContent(ib);

      std::vector<PlotUtils::MnvH1D*> mc_region;
      if (region==0) mc_region = finalSignalScaled;
      else if (region==1) mc_region = meanScaled;
      else if (region==2) mc_region = michelScaled;
      double predicted = 0;
      for (auto cat: mc_region){
	predicted += cat->GetBinContent(ib);
      }
      predicted *= mcScale;

      double chi2_bin = 0.0;
      if (predicted == 0.0 && observed == 0.0) {
	chi2_bin = 0.0;
      } else if (observed == 0.0) {
	chi2_bin = predicted * predicted;  // sigma^2 = 1 fallback
      } else if (predicted == 0.0) {
	chi2_bin = observed;  // (obs-0)^2/obs = obs, print warning
	std::cerr << "Warning: zero prediction with " << observed 
		  << " observed events in bin " << ib << std::endl;
      } else {
	chi2_bin = pow(observed - predicted, 2) / observed;
      }
      if (region==0) chi2_SR += chi2_bin;
      else if (region==1) chi2_MF += chi2_bin;
      else if (region==2) chi2_Michel += chi2_bin;
    }
  }

  int ndof;
  if ( (method==0) | (method==1) ) { ndof = 3*nbins - 2; } //for the normOnly fit 
  else if ( (method==2) | (method==3) ) { ndof = nbins; }

  double chi2_total_per_ndof = (chi2_SR + chi2_MF + chi2_Michel) / ndof;
  std::cout << "Total chi2 per ndof = " << chi2_total_per_ndof << std::endl;
  std::cout << "  chi2 contribution from signal region = " << chi2_SR << std::endl;
  std::cout << "  chi2 contribution from meanFront SB  = " << chi2_MF << std::endl;
  std::cout << "  chi2 contribution from Michel SB     = " << chi2_Michel << std::endl;

  //Append chi^2 results to a csv, useful for evaluating a bunch of fits at once
  std::ofstream csv("chi2_results.csv", std::ios::app);  // append mode
  csv << varName << ", " << methodName << ", " << chi2_total_per_ndof << ", "
      << chi2_SR << ", " << chi2_MF << ", " << chi2_Michel << "\n";

  // ---------------------------
  // Write output root file containing the signal-region scaled MC histograms with original names.
  // give em the same names as the originals so ExtractCrossSection works as intended
  // ---------------------------
  //TFile* outFile = TFile::Open("scaled_mc.root", "RECREATE"); 
  TFile* outFile = TFile::Open(("scaled_" + varName + "_" + methodName + "_mc.root").c_str(), "RECREATE");
  if (!outFile || outFile->IsZombie()) {
    std::cerr << "ERROR: couldn't open scaled_mc.root for writing\n";
  } else {
    std::vector<std::string> modifiedNames = {};
    for (auto cat: bkgdCategoryNames){
      modifiedNames.push_back(varName + "_" + cat);
    }
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
        
    // flush and close
    outFile->Close();
    delete outFile;
  }

  delete dataFile;
  delete mcFile;
  
  gROOT->GetListOfFunctions()->Delete();
  std::cout << "Done.\n";
  return 0;
}
