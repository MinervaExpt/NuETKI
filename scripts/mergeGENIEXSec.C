// Script I whipped up to merge individual playlists of GENIEXSecLooper cross sections
// this was pretty tricky and I decided to do it using the event rates from the genie file
// and flux histos from the corresponding mc files run on the same playlist.
// Then I basically just re-extract the cross section after merging all of the event rates
//
// author - Carlos P.
// usage: root -l -b -q 'mergeGENIEXSec.C("genie_files.txt", "mc_files.txt", "GENIEXSECEXTRACT_all_playlists.root")'
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
  const char* MC_filelist,
  const char* outName)
{
  TH1::AddDirectory(false);
  
  std::vector<std::string> GENIEXSec_files;
  std::vector<std::string> mc_files;
  std::ifstream GENIE_infile(GENIEXSec_filelist);
  std::ifstream mc_infile(MC_filelist);
  std::string line;

  std::vector<std::string> vars = {"E_lep","E_avail","E_nu","Lepton_Pt","Lepton_Pl","Theta_lep","Proton_p","Proton_Pt","Theta_p","Proton_T","DeltaPt","DeltaPtX","DeltaPtY","DeltaPl","P_n","AlphaPt","PhiPt"};
  //std::vector<std::string> vars = {"Lepton_Pt"};
  //populate input file names to vectors of strings
  while (std::getline(GENIE_infile, line)) {
    if (!line.empty()) GENIEXSec_files.push_back(line);
  }
  while (std::getline(mc_infile, line)) {
    if (!line.empty()) mc_files.push_back(line);
  }
  
  //loop through mc files first, get necessary ingredients
  // flux histogram, pot, also n_nucleons (same for all of them)    
  std::map<std::string, std::vector<PlotUtils::MnvH1D*>> flux_histos; //each variable has 12 flux histograms, one per playlist
  std::vector<double> pots;
  double n_nucleons = 0;

  for (const auto& file : mc_files){
    TFile* mc_file = TFile::Open(file.c_str());
    //get pot
    TParameter<double>* pot = nullptr;
    mc_file->GetObject("POTUsed", pot);
    pots.push_back( pot->GetVal());

    //get n nucleons
    if (n_nucleons == 0){
      TParameter<double>* nucleons = nullptr;
      mc_file->GetObject("DeltaPt_fiducial_nucleons", nucleons); //they're all the same, so can just use 1
      n_nucleons = nucleons->GetVal();
      std::cout << "updating n_nucleons to " << n_nucleons << std::endl;
    }

    //get flux histograms, there are 12*17 (n_playlists*n_variables) of them
    for (const auto& var : vars){
      std::string flux_hist_name = var + std::string("_reweightedflux_integrated");
      PlotUtils::MnvH1D* flux_hist = nullptr;
      mc_file->GetObject(flux_hist_name.c_str(), flux_hist);
      flux_hist->SetDirectory(nullptr);
      flux_histos[var].push_back(flux_hist);
    }
    mc_file->Close();
  }

  //merge the cross sections by summing (evRate/(fluxIntegral*pot))
  //then divide by bin width and n targets
  TFile* fout = new TFile(outName, "RECREATE");
  for (const auto& var : vars){
    std::cout << "merging xsec for var " << var << std::endl;
    PlotUtils::MnvH1D* evRate_sum = nullptr;
    PlotUtils::MnvH1D* hXSec = nullptr;
    double pot_sum = 0;
    for (int i = 0; i < GENIEXSec_files.size(); ++i)
      {
	pot_sum += pots[i];
	TFile* genie_file = TFile::Open(GENIEXSec_files[i].c_str());
	std::string evRate_name = var + std::string("_xsec_evRate");
	PlotUtils::MnvH1D* flux_hist = flux_histos[var][i];	

	TObject* obj = nullptr;
	genie_file->GetObject(evRate_name.c_str(), obj);
	if (auto* evRate = dynamic_cast<PlotUtils::MnvH1D*>(obj)) {
	  //PlotUtils::MnvH1D* evRate = (PlotUtils::MnvH1D*)genie_file->Get(evRate_name.c_str());

	  if (!evRate_sum){
	    evRate->SetDirectory(nullptr);
	    auto* clone = static_cast<PlotUtils::MnvH1D*>(evRate->Clone());
	    clone->SetDirectory(nullptr);

	    evRate_sum = static_cast<PlotUtils::MnvH1D*>(clone->Clone());;
	    evRate_sum->SetDirectory(nullptr);

	    clone->PopVertErrorBand("GENIE");
	    clone->Divide(clone, flux_hist);
	    clone->Scale(1.e4); //flux histogram is in units of m^-2 but want to report in cm^-2
	    hXSec = clone;
	    hXSec->SetDirectory(nullptr);
	  }
	  else{
	    evRate_sum->Add(evRate);

	    evRate->PopVertErrorBand("GENIE");
	    evRate->Divide(evRate, flux_hist);
	    evRate->Scale(1.e4); //flux histogram is in units of m^-2 but want to report in cm^-2
	    hXSec->Add(evRate);
	  }
	}

	//double fluxIntegral = flux_hist->Integral();
	//Phi_total += fluxIntegral * pots[i];

	genie_file->cd();
	genie_file->Clear();
	genie_file->Close();
	delete genie_file;
      }
    fout->cd();
    evRate_sum->Write();
    //std::cout << "pot_sum: " << pot_sum << std::endl;    
    //PlotUtils::MnvH1D* hXSec = (PlotUtils::MnvH1D*)evRate_sum->Clone("hXSec");  
    std::string xsec_name = var + "_xsec";
    hXSec->SetName(xsec_name.c_str());
    double globalScale = 1.0 / (n_nucleons*pot_sum);
    hXSec->Scale(1./(n_nucleons*pot_sum));
    hXSec->Scale(1.0, "width");
    /*
    for (int b = 1; b <= hXSec->GetNbinsX(); ++b)
      {
	//double N = evRate_sum->GetBinContent(b);
	double N = hXSec->GetBinContent(b);
	double width = hXSec->GetBinWidth(b);	
	double sigma = N / (n_nucleons * width * pot_sum); //cross section in that bin
	//double sigma = N / (n_nucleons * width); //cross section in that bin
	if (b==3){
	  std::cout << "n_nucleons: " << n_nucleons << std::endl;
	  std::cout << "bin width: " << width << std::endl;
	}
	hXSec->SetBinContent(b, sigma);
	}*/
    
    
    //std::cout << "writing output histogram" << std::endl;
    fout->cd();
    hXSec->Write();
    delete hXSec;
    delete evRate_sum;
  }
  std::cout << "closing output file" << std::endl;
  fout->Write();
  fout->Close();
  delete fout;
}

