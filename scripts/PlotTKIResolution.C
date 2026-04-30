//Root macro that loops through an input TTree, preferably one that has already selected my events,
// and recalculates my analysis variables in true and reco space.
// I then make and save resolution plots. 
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TMath.h>
#include <Math/Vector3D.h>
#include <Math/RotationX.h>
#include <iostream>
#include <vector>

#include "event/CVUniverse.h"
#include "event/CCNuEEvent.h"

#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/Model.h"
#include "PlotUtils/FluxAndCVReweighter.h"
#include "PlotUtils/GENIEReweighter.h"
#include "PlotUtils/LowRecoil2p2hReweighter.h"
#include "PlotUtils/RPAReweighter.h"
#include "PlotUtils/LowQ2PiReweighter.h"
#include "PlotUtils/FSIReweighter.h"

const double M_n = 939.56536;  //MeV
const double M_p = 938.272013; //MeV
const double M_e = 0.51099895; //MeV
const double M_pi = 139.57039; //MeV
const double M_nucleon = (1.5*M_n+M_p)/2.5;
const double pi = 3.141592653589793;
const double beam_tilt = 0.05887; //in rads, about 3.373 deg

const std::map<int, std::string> runNumberToPlaylist = {
  {110000, "minervame1A"},
  {111000, "minervame1B"},
  {111030, "minervame1C"},
  {111100, "minervame1D"},
  {111325, "minervame1E"},
  {111490, "minervame1F"},
  {110150, "minervame1G"},
  {113000, "minervame1L"},
  {113020, "minervame1M"},
  {113270, "minervame1N"},
  {113375, "minervame1O"},
  {112000, "minervame1P"},
};

void PlotTKIResolution() {
  std::vector<std::string> variables = {
    "E_lep",
    "E_avail",
    "E_nu",
    "Lepton_Pt",
    "Lepton_pl",
    "Theta_lep",
    "proton_p",
    "proton_Pt",
    "Theta_p",
    "DeltaPt",
    "DeltaPtX",
    "DeltaPtY",
    "DeltaPl",
    "P_n",
    "AlphaPt",
    "PhiPt",
  };

  //for the x axes of each variable
  std::vector<std::pair<double, double>> axes_limits = { 
    {0, 20}, //E_lep
    {0, 1.2}, //E_avail
    {0, 20}, //E_nu
    {0, 1.6}, //lepton_pT
    {0, 20}, //lepton_pL
    {0, 40}, //lepton_theta
    {0, 1.4}, //proton_P
    {0, 1.1}, //proton_pT
    {0, 180}, //proton_theta
    {0, 1.4}, //delta_Pt
    {-1.5, 1.5}, //delta_PtX
    {-1.5, 1.5}, //delta_PtY
    {-0.5, 1}, //delta_Pl
    {0, 1}, //P_n
    {0, 180}, //alpha_Pt
    {0, 180}, //acoplanarity
  };

  //for the reco - true differences, these are in units of GeV or degrees
  std::vector<std::pair<double, double>> axes_limits_diff = { 
    {-1, 1}, //E_lep
    {-0.5, 0.5}, //E_avail
    {-1, 1}, //E_nu
    {-0.5, 0.5}, //lepton_pT
    {-1, 1}, //lepton_pL
    {-10, 10}, //lepton_theta
    {-0.25, 0.25}, //proton_P
    {-0.25, 0.25}, //proton_pT
    {-30, 30}, //proton_theta
    {-0.3, 0.3}, //delta_Pt
    {-0.3, 0.3}, //delta_PtX
    {-0.3, 0.3}, //delta_PtY
    {-0.2, 0.2}, //delta_Pl
    {-0.3, 0.3}, //P_n
    {-30, 30}, //alpha_Pt
    {-20, 20}, //acoplanarity
  };

  //binning for migration matrices
  std::vector<std::vector<double>> migration_binnings = { 
    {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,20}, //E_lep
    {0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.8, 1, 1.2}, //E_avail
    {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,20}, //E_nu
    {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5}, //lepton_pT
    {0, 3, 4.5, 5.5, 6.5, 8, 12, 30}, //lepton_pL
    {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40}, //lepton_theta
    {0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375}, //proton_P
    {0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.1}, //proton_pT
    {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180}, //proton_theta
    {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.4}, //delta_Pt
    {-2.0, -0.7, -0.4, -0.2, 0.0, 0.2, 0.4, 0.7, 2.0}, //delta_PtX
    {-2.0, -1.1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 1.3}, //delta_PtY
    {-0.25, 0.00, 0.10, 0.20, 0.30, 0.40, 0.60}, //delta_Pl
    {0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 1.0}, //P_n
    {0,20,40,60,80,100,120,140,160,180}, //alpha_Pt
    {0,10,30,50,80,120,180}, //acoplanarity
  };

    // Create histograms
  //std::vector<TH1D*> resolution_hists_1D;
  std::vector<TH2D*> resolution_hists;
  std::vector<TH2D*> diff_hists;
  std::vector<TH2D*> migration_hists;
  for (int i = 0; i < variables.size(); i++) {
    //std::string hist_name_1D = "h_" + variables[i] + "_resolution_1D";
    //std::string hist_label_1D = "; " + variables[i] + " resolution; events";
    //resolution_hists_1D.push_back(new TH1D(hist_name_1D.c_str(), hist_label_1D.c_str(), 50, -1, 1));

    std::string hist_name = "h_" + variables[i] + "_resolution";
    std::string hist_label = variables[i] + " resolution; true " + variables[i] + "; (reco - true)/true";
    resolution_hists.push_back(new TH2D(hist_name.c_str(), hist_label.c_str(), 50, axes_limits[i].first, axes_limits[i].second, 50, -0.5, 0.5));

    std::string hist_name_diff = "h_" + variables[i] + "_diff";
    std::string hist_label_diff;
    std::string migration_name = "h_" + variables[i] + "_migration";
    std::string migration_label;
    auto& bins = migration_binnings[i];
    if ( variables[i]=="Theta_lep" || variables[i]=="Theta_p" || variables[i]=="AlphaPt" || variables[i]=="PhiPt") {
      hist_label_diff = variables[i] + " reco - true diff; true " + variables[i] + " [degrees]; reco - true, [degrees]";
      migration_label = variables[i] + " migration; Reco " + variables[i] + " [degrees]; True " + variables[i] + " [degrees]";
    } else {
      hist_label_diff = variables[i] + " reco - true diff; true " + variables[i] + " [GeV]; reco - true, [GeV]";
      migration_label = variables[i] + " migration; Reco " + variables[i] + " [GeV]; True " + variables[i] + " [GeV]";
    }
    diff_hists.push_back(new TH2D(hist_name_diff.c_str(), hist_label_diff.c_str(), 50, axes_limits[i].first, axes_limits[i].second, 50, axes_limits_diff[i].first, axes_limits_diff[i].second));
    migration_hists.push_back(new TH2D(migration_name.c_str(), migration_label.c_str(), bins.size()-1, bins.data(), bins.size()-1, bins.data()));
  }

  TH2D* tejin_plot = new TH2D("dptx_vs_dpty", "dptx vs dpty; dptx [GeV]; dpty [GeV]", 50, -0.25, 0.25, 50, -0.25, 0.25);
  std::vector<double> avg_resolutions;
  for (const auto& var : variables) {
    avg_resolutions.push_back(0);
  }

  //Get my reweighters  
  std::vector<std::unique_ptr<PlotUtils::Reweighter<CVUniverse, CCNuEEvent>>> MyModel;
  MyModel.emplace_back(new PlotUtils::FluxAndCVReweighter<CVUniverse, CCNuEEvent>());
  MyModel.emplace_back(new PlotUtils::GENIEReweighter<CVUniverse, CCNuEEvent>(true, false));
  MyModel.emplace_back(new PlotUtils::LowRecoil2p2hReweighter<CVUniverse, CCNuEEvent>());
  MyModel.emplace_back(new PlotUtils::RPAReweighter<CVUniverse, CCNuEEvent>());
  MyModel.emplace_back(new PlotUtils::LowQ2PiReweighter<CVUniverse, CCNuEEvent>("JOINT"));
  MyModel.emplace_back(new PlotUtils::FSIReweighter<CVUniverse, CCNuEEvent>(true, false));  
  PlotUtils::Model<CVUniverse, CCNuEEvent> model(std::move(MyModel));

  //make my chainwrapper and CVUniverse, set my options
  PlotUtils::MinervaUniverse::SetNuEConstraint(true); //the neutrino-electron flux scattering constraint?? So should leave true?
  //PlotUtils::MinervaUniverse::SetPlaylist(options.m_plist_string);
  PlotUtils::MinervaUniverse::SetPlaylist("minervame1a");
  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(12); //Apparently this is: Analysis neutrino identity -- used for flux weights
  PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
  PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false); //to do: investigate
  PlotUtils::MinervaUniverse::RPAMaterials(false); //Sets whether you want to use different nuclei for RPA reweight, with false it just uses carbon which is fine for me                  
  ChainWrapper* chw = makeChainWrapperPtr("custom_tuple_path.txt", "MasterAnaDev"); //mc_file list is string
  CVUniverse* cvUniv = new CVUniverse(chw);

  std::string currentPlaylist = "minervame1A";
  int eventCount = 0;
  const Long64_t nEntries = chw->GetEntries();
  std::cout << "n entries = " << nEntries << std::endl;
  for (Long64_t i = 0; i < nEntries; ++i) {
    //if(i%1000==0) std::cout << i << " / " << nEntries << "\r" << std::flush;
    //if (i>100) break;
    CCNuEEvent cvEvent;
    cvUniv->SetEntry(i);
    int selectionCategory = cvUniv->GetInt("selectionCategory");
    std::vector<int> sideband_cut_results = cvUniv->GetVecInt("sideband_cut_results");

    //if (abs(mc_incoming) != 12) continue;
    //if (protonIndex < 0) continue;
    //if (!(sideband_cut_results[0] == 1 && sideband_cut_results[1] == 1 && sideband_cut_results[2] == 1)) continue;
    if (!(sideband_cut_results[0] == 1 && sideband_cut_results[1] == 1)) continue;
    if (selectionCategory != -999) continue;
    //if (mc_intType != 1) continue;

    //Some logic to set playlist to the right value whenever we cross a run number threshold, since my tuples are multi-playlist
    //required for flux weights to be set right
    int eventRun = cvUniv->GetInt("mc_run");
    auto found = runNumberToPlaylist.upper_bound(eventRun);
    if (found != runNumberToPlaylist.begin()) {
        std::string eventPlaylist = std::prev(found)->second;
        if (eventPlaylist != currentPlaylist) {
            currentPlaylist = eventPlaylist;
            PlotUtils::MinervaUniverse::SetPlaylist(currentPlaylist);
        }
    }    

    model.SetEntry(*cvUniv, cvEvent);
    const double cvWeight = model.GetWeight(*cvUniv, cvEvent);
    //const double cvWeight = 1;
    //std::cout << "\nlooking at event: " << i << "\n" << std::endl;    
    //std::cout << "signal event: " << cvUniv->GetInt("mc_run") << " | " << cvUniv->GetInt("mc_subrun") << " | " << cvUniv->GetInt("mc_nthEvtInFile")+1 << std::endl;
    
    ROOT::Math::RotationX r(-beam_tilt); 
    
    //==============================
    // True kinematics
    //==============================
    double true_E_lep = cvUniv->GetElectronEnergyTrueGeV();
    double true_E_avail = cvUniv->GetEavailTrueGeV();
    double true_E_nu = cvUniv->GetEnuTrueGeV();
    
    double true_lepton_pT = cvUniv->GetElectronPtTrueGeV();
    double true_lepton_pL = cvUniv->GetElectronPParallelTrueGeV();
    double true_lepton_theta = cvUniv->GetElectronThetaDegTrue();

    double true_proton_P = cvUniv->GetProtonPTrueGeV();
    double true_proton_pT = cvUniv->GetProtonPtTrueGeV();
    double true_proton_theta = cvUniv->GetProtonThetaDegTrue();

    double true_delta_Pt = cvUniv->GetDeltaPtTrueGeV();
    double true_delta_PtX = cvUniv->GetDeltaPtXTrueGeV();
    double true_delta_PtY = cvUniv->GetDeltaPtYTrueGeV();
    double true_delta_Pl = cvUniv->GetDeltaPlTrueGeV();
    double true_P_n = cvUniv->GetPnTrueGeV();
    double true_alpha_Pt = cvUniv->GetAlphaTTrue();
    double true_acoplanarity = cvUniv->GetPhiTTrue();
    
    //==============================
    // Reco kinematics
    //==============================
    double reco_E_lep = cvUniv->GetElectronEnergyGeV();
    double reco_E_avail = cvUniv->GetEavailGeV();
    double reco_E_nu = cvUniv->GetEnuGeV();
    
    double reco_lepton_pT = cvUniv->GetElectronPtGeV();
    double reco_lepton_pL = cvUniv->GetElectronPParallelGeV();
    double reco_lepton_theta = cvUniv->GetElectronThetaDeg();

    double reco_proton_P = cvUniv->GetProtonPGeV();
    double reco_proton_pT = cvUniv->GetProtonPtGeV();
    double reco_proton_theta = cvUniv->GetProtonThetaDeg();

    double reco_delta_Pt = cvUniv->GetDeltaPtGeV();
    double reco_delta_PtX = cvUniv->GetDeltaPtXGeV();
    double reco_delta_PtY = cvUniv->GetDeltaPtYGeV();
    double reco_delta_Pl = cvUniv->GetDeltaPlGeV();
    double reco_P_n = cvUniv->GetPnGeV();
    double reco_alpha_Pt = cvUniv->GetAlphaT();
    double reco_acoplanarity = cvUniv->GetPhiT();

    //now that I have reco & true, fill this vec of vecs with the output values
    //format is {reco, true, resolution = (reco-true)/true, diff = reco-true} for each variable. 
    std::vector<std::vector<double>> output_values = {
      {reco_E_lep, true_E_lep, (reco_E_lep - true_E_lep)/true_E_lep, (reco_E_lep - true_E_lep)},
      {reco_E_avail, true_E_avail, (reco_E_avail - true_E_avail)/true_E_avail, (reco_E_avail - true_E_avail)},
      {reco_E_nu, true_E_nu, (reco_E_nu - true_E_nu)/true_E_nu, (reco_E_nu - true_E_nu)},
      {reco_lepton_pT, true_lepton_pT, (reco_lepton_pT - true_lepton_pT)/true_lepton_pT, (reco_lepton_pT - true_lepton_pT)},
      {reco_lepton_pL, true_lepton_pL, (reco_lepton_pL - true_lepton_pL)/true_lepton_pL, (reco_lepton_pL - true_lepton_pL)},
      {reco_lepton_theta, true_lepton_theta, (reco_lepton_theta - true_lepton_theta)/true_lepton_theta, (reco_lepton_theta - true_lepton_theta)},
      {reco_proton_P, true_proton_P, (reco_proton_P - true_proton_P)/true_proton_P, (reco_proton_P - true_proton_P)},
      {reco_proton_pT, true_proton_pT, (reco_proton_pT - true_proton_pT)/true_proton_pT, (reco_proton_pT - true_proton_pT)},
      {reco_proton_theta, true_proton_theta, (reco_proton_theta - true_proton_theta)/true_proton_theta, (reco_proton_theta - true_proton_theta)},
      {reco_delta_Pt, true_delta_Pt, (reco_delta_Pt - true_delta_Pt)/true_delta_Pt, (reco_delta_Pt - true_delta_Pt)},
      {reco_delta_PtX, true_delta_PtX, (reco_delta_PtX - true_delta_PtX)/true_delta_PtX, (reco_delta_PtX - true_delta_PtX)},
      {reco_delta_PtY, true_delta_PtY, (reco_delta_PtY - true_delta_PtY)/true_delta_PtY, (reco_delta_PtY - true_delta_PtY)},
      {reco_delta_Pl, true_delta_Pl, (reco_delta_Pl - true_delta_Pl)/true_delta_Pl, (reco_delta_Pl - true_delta_Pl)},
      {reco_P_n, true_P_n, (reco_P_n - true_P_n)/true_P_n, (reco_P_n - true_P_n)},
      {reco_alpha_Pt, true_alpha_Pt, (reco_alpha_Pt - true_alpha_Pt)/true_alpha_Pt, (reco_alpha_Pt - true_alpha_Pt)},
      {reco_acoplanarity, true_acoplanarity, (reco_acoplanarity - true_acoplanarity)/true_acoplanarity, (reco_acoplanarity - true_acoplanarity)}
    };

    //tejins fancy 2d contour plot of dpty and dptx using truth
    tejin_plot->Fill(output_values[10][1], output_values[11][1]);
    for (int i = 0; i < resolution_hists.size(); i++) {
      //resolution_hists_1D[i]->Fill(resolution_values[i]);
      resolution_hists[i]->Fill(output_values[i][1], output_values[i][2], cvWeight);
      diff_hists[i]->Fill(output_values[i][1], output_values[i][3], cvWeight);
      migration_hists[i]->Fill(output_values[i][0], output_values[i][1], cvWeight);
	
      avg_resolutions[i] += output_values[i][2];
    }
    eventCount++;
  }

  std::cout << " ----------------- FINISHED EVENT LOOP ----------------- " << std::endl;
  std::cout << "events filled: " << eventCount << std::endl;
  /*
  for (int i = 0; i < variables.size(); i++) {
    std::cout << "Avg " << variables[i] << " resolution = " << avg_resolutions[i]/eventCount << std::endl;
    }
*/

  TFile *f = new TFile("AlphaPt_Migration.root", "RECREATE");
  migration_hists[14]->Write();
  f->Close();

  PlotUtils::MnvPlotter plotter;
  TCanvas* c = new TCanvas("c", "", 800, 600);
  gErrorIgnoreLevel = kWarning;
  //tejin_plot->Draw("colz");
  //c->SaveAs("tejin_plot_QE_only.png");
  for (int i = 0; i < variables.size(); i++) {
    /*
    c->Clear();
    //c->SetLogz();
    resolution_hists_1D[i]->Draw("hist");
    c->Update();
    std::string outName_1D = variables[i] + "_resolution_1D.png";
    c->SaveAs(outName_1D.c_str());    
    */

    c->Clear();
    //c->SetLogz();
    resolution_hists[i]->Draw("colz");
    c->Update();
    std::string outName = "resolutions/" + variables[i] + "_resolution.png";
    c->SaveAs(outName.c_str());    

    
    c->Clear();
    //c->SetLogz();
    diff_hists[i]->Draw("colz");
    c->Update();
    std::string outName2 = "diffs/" + variables[i] + "_diff.png";
    c->SaveAs(outName2.c_str());    
    

    c->Clear();
    //c->SetLogz();
    //migration_hists[i]->Draw("colz"); //-> replace with mnv plotter row norm version...
    plotter.DrawNormalizedMigrationHistogram(migration_hists[i]);
    c->Update();
    std::string outName3 = "migrations/" + variables[i] + "_migration.png";
    c->SaveAs(outName3.c_str());      }

  /*
  auto leg = new TLegend(0.65,0.70,0.88,0.88);
  leg->AddEntry(h_resolution1, "No ESC/Chi2 cut", "1");
  leg->AddEntry(h_resolution2, "Looser ESC/Chi2 cut", "1");
  leg->AddEntry(h_resolution3, "Original ESC/Chi2 cut", "1");
  leg->Draw();
  c->SetLogy(1);
  c->Update();
  */
  
  //delete file;
  }
