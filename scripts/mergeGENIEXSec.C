// Script I whipped up to merge individual playlists of GENIEXSecLooper cross sections
// because cross section = n events / normalization, and GENIEXSecExtract saves the event rates,
// you can get a normalization (flux * targets) per playlist just by dividing the event rate by the cross section
// which should be completely flat.
// Then do sum all the event rate histograms from each playlist, and divide by the SUM of the scalars obtained from the above division.
//
// author - Carlos P.
// usage: root -b -q 'mergeGENIEXSec.C("genie_files.txt", "GENIEXSECEXTRACT_all_playlists.root")'
// genie_files.txt = text file with paths to the input genie xsecs
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

void mergeGENIEXSec(
  const char* GENIEXSec_filelist,
  const char* outName)
{
  TH1::AddDirectory(false);
  
  std::vector<std::string> GENIEXSec_files;
  std::ifstream GENIE_infile(GENIEXSec_filelist);
  std::string line;

  std::vector<std::string> vars = {"E_lep","E_avail","E_nu","Lepton_Pt","Lepton_Pl","Theta_lep","Proton_p","Proton_Pt","Theta_p","Proton_T","DeltaPt","DeltaPtX","DeltaPtY","DeltaPl","P_n","AlphaPt","PhiPt"};
  //std::vector<std::string> vars = {"Lepton_Pt"};
  //populate input file names to vectors of strings
  while (std::getline(GENIE_infile, line)) {
    if (!line.empty()) GENIEXSec_files.push_back(line);
  }

  //merge the cross sections by summing (evRate/(fluxIntegral*pot))
  //then divide by bin width and n targets
  TFile* fout = new TFile(outName, "RECREATE");
  std::map <string, PlotUtils::MnvH1D*> evRate_sums;
  std::map <string, PlotUtils::MnvH1D*> norm_hists;
  std::map <string, double> total_norms;
  for (int i=0; i < GENIEXSec_files.size(); i++){

    TFile* genie_file = TFile::Open(GENIEXSec_files[i].c_str());
    if (!genie_file || genie_file->IsZombie()) {
      std::cerr << "Error opening file " << GENIEXSec_files[i] << "\n";
      return;
    }

    for (const auto& var : vars){
      PlotUtils::MnvH1D* xsec = (PlotUtils::MnvH1D*)genie_file->Get((var+"_xsec").c_str());
      PlotUtils::MnvH1D* evRate = (PlotUtils::MnvH1D*)genie_file->Get((var+"_xsec_evRate").c_str());

      //Undo bin width normalization for xsec, since evRate isn't bin width normalized
      for (int j=1; j < xsec->GetNbinsX()+1; j++){
	double bin_width = xsec->GetBinWidth(j);
	double bin_content = xsec->GetBinContent(j);
	xsec->SetBinContent(j, bin_content*bin_width);
      }

      PlotUtils::MnvH1D* norm_hist = evRate->Clone((var+"_normalization").c_str());
      norm_hist->SetDirectory(nullptr);
      norm_hist->Divide(norm_hist, xsec); //the result of this should be completely flat, other than bins with no events (which are zero)
      //this little loop is just to ensure we get a proper bin value which isn't 0,
      // some vars start with 0 bin counts (E_lep), some end with it, some have some in the middle etc. For some reason GetMaximum() doesn't work.
      // I was getting the maximum negative value of a 32 bit float, so something got messed up. 
      double norm = 0;
      for (int j=1; j < norm_hist->GetNbinsX()+1; j++){
	norm = norm_hist->GetBinContent(j);
	if (norm != 0){ break; }
      }
      auto it = evRate_sums.find(var);
      if (it == evRate_sums.end()){
	PlotUtils::MnvH1D* evRate_sum = evRate->Clone((var+"_xsec_evRate").c_str());
	evRate_sum->SetDirectory(nullptr);
	evRate_sums.insert({var, evRate_sum});
	//norm_hists.insert({var, norm_hist});
	total_norms.insert({var, norm});
      }
      else {
	evRate_sums[var]->Add(evRate);
	//norm_hists[var]->Add(norm_hist);
	total_norms[var] += norm;
      }
    }
    //genie_file->cd();
    //genie_file->Clear();
    genie_file->Close();
    //delete genie_file;
  }

  //Now use the summed event rates and total normalizations to extract cross sections
  // then normalize by bin width, and write to output file
  for (const auto& var : vars){
    PlotUtils::MnvH1D* xsec = evRate_sums[var]->Clone("cross_section");
    xsec->SetName((var+"_xsec").c_str());
    xsec->Scale(1.0/total_norms[var]);
    xsec->Scale(1.0, "width");
    //for (int j=1; j < xsec->GetNbinsX()+1; j++){
    //double bin_width = xsec->GetBinWidth(j);
    //double bin_content = xsec->GetBinContent(j);
    //xsec->SetBinContent(j, bin_content/bin_width);
    //}
    fout->cd();
    evRate_sums[var]->Write();
    //norm_hists[var]->Write();
    xsec->Write();
  }

  fout->Write();
  fout->Close();
  delete fout;
  std::cout << "Success" << std::endl;
}

