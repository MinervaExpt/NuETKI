//Quick little script that loops through all universes in an MnvH1D and prints out bin contents for each one.
//I notice that in some cases my universes aren't actually different than the CV, is that a problem??

#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvVertErrorBand.h"
#include "TFile.h"
#include "TH1D.h"

using namespace std;

//void loopUniverses(const char* filepath = "/exp/minerva/data/users/cpernas/NuE_TKI/Nov_07_per_universe_scaling/mc.root",
//		   const char* histName = "DeltaPt_data"){
void loopUniverses(){
  //const char* filepath = "/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/test/MC_Jan_28.root";
  const char* filepath = "/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/MC_Feb_05.root";
  const char* histName = "DeltaPt_data";

  
  //has to be an mc file
  TFile* inFile = TFile::Open(filepath, "READ");
  if (!inFile || inFile->IsZombie()) {
    std::cerr << "Failed to open data file: " << filepath << std::endl;
    return;
  }

  PlotUtils::MnvH1D* hist = nullptr;
  inFile->GetObject(histName, hist);
  if (!hist || hist->IsZombie()) {
    std::cerr << "Failed to get hist: " << histName << std::endl;
    return;
  }

  
  int nbins = hist->GetNbinsX();
  vector<double> cv_bin_contents(nbins);
  for (int i=1; i<=nbins; i++){
    cv_bin_contents[i-1] = hist->GetBinContent(i);
  }

  vector<string> unchanged_errorbands; //list of universes that aren't changed wrt CV
  int total_universes = 0;
  bool sameAsCV;
  vector<double> universe_bin_contents(nbins);
  vector<string> errorBandNames = hist->GetVertErrorBandNames();
  for (auto& bandName : errorBandNames){
    PlotUtils::MnvVertErrorBand* band = hist->GetVertErrorBand(bandName);
    vector<TH1D*> universes = band->GetHists();
    int i = 0;
    sameAsCV = false;
    for (auto& universe : universes){
      //cout << "================ Band " << bandName << " Universe " << i << " ================" << endl;
      //cout << "( ";
      for (int j=1; j<=nbins; j++){
	universe_bin_contents[j-1] = universe->GetBinContent(j);
	if (j==nbins){
	  //cout << universe_bin_contents[j-1] << " )" << endl;
	} else {
	  //cout << universe_bin_contents[j-1] << ", ";
	}	
      }
      if (universe_bin_contents == cv_bin_contents){
	sameAsCV = true;
      }
      i++;
      total_universes++;
      //cout << endl;
    }
    if (sameAsCV){
      unchanged_errorbands.push_back(bandName);
    }
  }

  //finally print out cv bins at the very end to compare
  cout << "================ CV ==================" << endl;
  cout << "( ";
  for (int i=0; i<nbins; i++){
    if (i==(nbins-1)){
      cout << cv_bin_contents[i] << " )" << endl;
    } else {
      cout << cv_bin_contents[i] << ", ";
    }
  }

  cout << "Done, MnvH1D contains " << errorBandNames.size() << " error bands and " << total_universes << " total universes." << endl;
  cout << "There were " << unchanged_errorbands.size() <<" error bands with identical bin contents to the CV:" << endl;
  for (auto& name : unchanged_errorbands) { cout << name << ", "; }
  cout << endl;
}
