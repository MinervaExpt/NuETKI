#include "PlotUtils/MnvH1D.h"

#include <TFile.h>
#include <vector>
void doClosureTest(
  const char* xsec_file_name,
  const char* genieXSec_file_name,
  const char* var)
{
  TH1::AddDirectory(false);

  TFile* xsec_file = TFile::Open(xsec_file_name);
  TFile* genie_file = TFile::Open(genieXSec_file_name);

  std::string mc_xsec_name = "simulatedCrossSection"; //ExtractCrossSection saves the variable name in the file name, not on the hist itself
  std::string genie_xsec_name = std::string(var) + "_xsec";
  std::string hist_name = std::string(var) + " Closure Test";
  std::string out_name = std::string(var) + "_closure.png";
  
  PlotUtils::MnvH1D* mc_xsec = (PlotUtils::MnvH1D*)xsec_file->Get(mc_xsec_name.c_str());
  PlotUtils::MnvH1D* genie_xsec = (PlotUtils::MnvH1D*)genie_file->Get(genie_xsec_name.c_str());

  std::vector<std::string> error_band_names = mc_xsec->GetVertErrorBandNames();

  //don't need error bands for closure test, but they mess with the division if the denominator doesn't have all the ones in the numerator
  for (auto& band : error_band_names){
    mc_xsec->PopVertErrorBand(band);
  }
  mc_xsec->Divide(mc_xsec, genie_xsec);
  mc_xsec->SetName(hist_name.c_str());
  
  TCanvas *c1 = new TCanvas("c1", "Closure Test", 800, 600);
  mc_xsec->Draw("hist");
  c1->SaveAs(out_name.c_str());
  
  xsec_file->Close();
  genie_file->Close();
}
