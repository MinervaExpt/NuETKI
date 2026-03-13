//Usage: root -l -b -q 'SumHistogramsFromFiles.C("MCHists_Feb_25.txt", "mc_merged.root")'

//Script to merge output of my event loop job histograms
//Sums POT and entries from histograms of various playlists, but have to be careful with:

// fiducial_nucleons -> definitely DON'T sum these, they're all the same so just take one

// reweighted_flux_integrated -> this ones tricky, the units for the flux histo are per POT so they stay
// relatively constant across playlists, but not exactly...
// me1A, me1B, me1C, me1D, me1E, me1F: reweightedflux_integrated->Integral() = 5.0354872e-05
//                   me1G, me1L, me1M: reweightedflux_integrated->Integral() = 5.0100729e-05
//                   me1N, me1O, me1P: reweightedflux_integrated->Integral() = 4.9288562e-05
// I DO NOT KNOW if this is the correct thing to do, but I will take the weighted average (by POT) of each of these.
// they're pretty close so it shouldn't make a huge difference? I can also check how GENIEXSecExtract does it I guess
// since you can just run it on all the playlists at once

//Minerva includes
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"

#include <TFile.h>
#include <TSystem.h>
#include <TH1.h>
#include <TKey.h>
#include <TClass.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

void SumHistogramsFromFiles(const char* filelistPath = "filelist.txt", const char* outputPath = "summedOutput.root") {
  std::vector<std::string> inputFiles;
  std::ifstream infile(filelistPath);
  std::string line;

  // Read input file paths
  while (std::getline(infile, line)) {
    if (!line.empty()) inputFiles.push_back(line);
  }

  if (inputFiles.empty()) {
    std::cerr << "Error: No input files found in " << filelistPath << std::endl;
    return;
  }

  std::cout << "Found " << inputFiles.size() << " files to process.\n";

  // Open first file to get list of histogram names
  TFile* firstFile = TFile::Open(inputFiles[0].c_str());
  if (!firstFile || firstFile->IsZombie()) {
    std::cerr << "Error opening file: " << inputFiles[0] << std::endl;
    return;
  }

  std::vector<std::string> histNames;
  TIter next(firstFile->GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)next())) {
    TObject* obj = key->ReadObj();
    //if (obj->IsA() == PlotUtils::MnvH1D::Class()) {
    //if (obj->InheritsFrom(TH1::Class())) { //this allows the MnvH2Ds as well
    histNames.push_back(key->GetName());
    std::cout << "key name: " << key->GetName() << ", obj name: " << obj->GetName() << ", type: (" << key->GetClassName() << ")" << std::endl;
      //}
    delete obj;
  }
  firstFile->Close();
  std::cout << histNames.size() << " keys found" << std::endl;
  //for (const auto& name : histNames) std::cout << name << std::endl;
  std::cout << std::endl;

  // Maps to hold summed outputs
  std::map<std::string, PlotUtils::MnvH1D*> summedHists;
  std::map<std::string, PlotUtils::MnvH2D*> summedHists2D;
  std::map<std::string, double> summedParams;

  std::cout << inputFiles.size() << " files found. " << std::endl;
  bool isNotFirstFile = false;
  // Loop through each file
  for (const auto& filepath : inputFiles) {
    TFile* file = TFile::Open(filepath.c_str());
    std::cout << "- FILE - " << filepath << std::endl;
    if (!file || file->IsZombie()) {
      std::cerr << "Warning: Skipping file " << filepath << std::endl;
      continue;
    }
    double filePOT;
    for (const auto& name : histNames) {
      TObject* obj = nullptr;
      file->GetObject(name.c_str(), obj);
      if (!obj) {
	std::cerr << "Warning: Histogram " << name << " not found in file " << filepath << std::endl;
	continue;
      }
      if ((name.find("fiducial_nucleons") != std::string::npos) && isNotFirstFile){ continue; } //copy # of nucleons for the first file, then don't sum anymore
      else if (name.find("reweightedflux_integrated") == std::string::npos){ //sum everything that's not n_nucleons or flux integrated (cause its per POT)
	//PlotUtils::MnvH1D* hist = nullptr;
	//file->GetObject(name.c_str(), hist);
	//for the MnvH1Ds
	if (auto* h1 = dynamic_cast<PlotUtils::MnvH1D*>(obj)) {
	  if (summedHists.find(name) == summedHists.end()) {
	    h1->SetDirectory(nullptr);
	    auto* clone = static_cast<PlotUtils::MnvH1D*>(h1->Clone());
	    clone->SetDirectory(nullptr);
	    clone->SetName(name.c_str());
	    summedHists[name] = clone;
	  } else {
	    summedHists[name]->Add(h1);
	  }
	} //For the MnvH2Ds
	else if (auto* h2 = dynamic_cast<PlotUtils::MnvH2D*>(obj)) {
	  if (summedHists2D.find(name) == summedHists2D.end()) {
	    h2->SetDirectory(nullptr);
	    auto* clone = static_cast<PlotUtils::MnvH2D*>(h2->Clone());
	    clone->SetDirectory(nullptr);
	    clone->SetName(name.c_str());
	    summedHists2D[name] = clone;
	  } else {
	    summedHists2D[name]->Add(h2);
	  }
	} else if (auto* param = dynamic_cast<TParameter<double>*>(obj)) {
	  if (name.find("POTUsed") != std::string::npos){ filePOT = param->GetVal(); } //keep track of this file's POT for the flux histogram
	  summedParams[name] += param->GetVal();
	}
      }
      else if (name.find("reweightedflux_integrated") != std::string::npos){ //weighted average of reweightedflux_integrated by POT
	if (auto* h1 = dynamic_cast<PlotUtils::MnvH1D*>(obj)) {
	  if (summedHists.find(name) == summedHists.end()) {
	    h1->SetDirectory(nullptr);
	    auto* clone = static_cast<PlotUtils::MnvH1D*>(h1->Clone());
	    clone->SetDirectory(nullptr);
	    clone->SetName(name.c_str());
	    clone->Scale(filePOT); //we want to sum the flux hist * POT, which we'll then divide by total POT
	    summedHists[name] = clone;
	  } else {
	    h1->Scale(filePOT);
	    summedHists[name]->Add(h1);
	  }
	}
      }
      else { std::cerr << "Warning: Object \"" << name << "\" is neither MnvH1D or MnvH2D or TParameter<double>\n"; }
      //std::cout << "hist: " << name << ", # entries: " << hist->GetEntries() << std::endl;
    }
    isNotFirstFile = true;
    file->cd();
    file->Clear();
    file->Close();
    delete file;
  }
  // Write summed histograms to output file
  TFile* outFile = TFile::Open(outputPath, "RECREATE");
  for (const auto& entry : summedHists) {
    //divide the flux histogram sums by total POT
    if (entry.first.find("reweightedflux_integrated") != std::string::npos){
      entry.second->Scale(1./summedParams["POTUsed"]);
    }
    entry.second->Write();
  }
  for (const auto& entry : summedHists2D) {
    entry.second->Write();
  }
  for (const auto& [entry, value] : summedParams) {
    TParameter<double> param(entry.c_str(), value);
    param.Write();
  }
  outFile->Close();
  //for (auto& pair : summedHists) {
  //delete pair.second;
  //}

  std::cout << "Summed histograms written to: " << outputPath << std::endl;
}
