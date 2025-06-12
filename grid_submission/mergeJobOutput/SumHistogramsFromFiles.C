//Usage: root -l -b -q 'SumHistogramsFromFiles.C("mcHistNames.txt", "mc_sum.root")'

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
  std::cout << histNames.size() << " Histograms found" << std::endl;
  //for (const auto& name : histNames) std::cout << name << std::endl;
  std::cout << std::endl;

  // Maps to hold summed outputs
  std::map<std::string, PlotUtils::MnvH1D*> summedHists;
  std::map<std::string, PlotUtils::MnvH2D*> summedHists2D;
  std::map<std::string, double> summedParams;
  
  // Loop through each file
  for (const auto& filepath : inputFiles) {
    TFile* file = TFile::Open(filepath.c_str());
    std::cout << "- FILE - " << filepath << std::endl;
    
    if (!file || file->IsZombie()) {
      std::cerr << "Warning: Skipping file " << filepath << std::endl;
      continue;
    }

    for (const auto& name : histNames) {
      //PlotUtils::MnvH1D* hist = nullptr;
      //file->GetObject(name.c_str(), hist);
      TObject* obj = nullptr;
      file->GetObject(name.c_str(), obj);
      if (!obj) {
	std::cerr << "Warning: Histogram " << name << " not found in file " << filepath << std::endl;
	continue;
      }
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
      //For the MnvH2Ds
      } else if (auto* h2 = dynamic_cast<PlotUtils::MnvH2D*>(obj)) {
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
	summedParams[name] += param->GetVal();
      }

      else { std::cerr << "Warning: Object \"" << name << "\" is neither MnvH1D or MnvH2D or TParameter<double>\n"; }
      //std::cout << "hist: " << name << ", # entries: " << hist->GetEntries() << std::endl;
    
    }
    file->cd();
    file->Clear();
    file->Close();
    delete file;
  }

  // Write summed histograms to output file
  TFile* outFile = TFile::Open(outputPath, "RECREATE");
  for (const auto& entry : summedHists) {
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
