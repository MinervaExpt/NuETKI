#include "GENIEXSecExtract/XSecLooper.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>

#include <cstdlib>
typedef unsigned int uint;

class MinModDepCCQEXSec : public XSec
{
public:
  MinModDepCCQEXSec(const char* name)
    :XSec(name)
  {
  };

bool isNuEQELikeSignal( PlotUtils::ChainWrapper& chw, int entry )
{
  int mc_incoming              = static_cast<int>(chw.GetValue("mc_incoming",entry)); //should be 12 or -12
  int current                  = static_cast<int>(chw.GetValue("mc_current",entry));
  std::vector<double> electron_4vec(4);
  for (int i = 0; i < 4; ++i) { electron_4vec[i] = chw.GetValue("mc_primFSLepton", entry, i); }

  int nFS = static_cast<int>(chw.GetValue("mc_nFSPart",entry));

  std::vector<int> FSPartPDG(nFS);
  std::vector<double> FSPartE(nFS), FSPartPx(nFS), FSPartPy(nFS), FSPartPz(nFS);
  for (int i = 0; i < nFS; ++i) { 
    FSPartPDG[i] = static_cast<int>(chw.GetValue("mc_FSPartPDG",entry, i));
    FSPartE[i]   = chw.GetValue("mc_FSPartE", entry, i);
    FSPartPx[i]  = chw.GetValue("mc_FSPartPx", entry, i);
    FSPartPy[i]  = chw.GetValue("mc_FSPartPy", entry, i);
    FSPartPz[i]  = chw.GetValue("mc_FSPartPz", entry, i);
  }
  
  bool hasSignalProton = false;
  bool hasMeson = false;
  bool hasPhoton = false;
  for (int i=0; i < FSPartPDG.size(); i++){
    if (FSPartPDG[i] == 2212){
      double protonP = sqrt(pow(FSPartE[i],2) - pow(938.272013, 2));

      double numi_beam_angle_rad = -0.05887;
      double pyprime = -1.0*sin(numi_beam_angle_rad)*FSPartPz[i] + cos(numi_beam_angle_rad)*FSPartPy[i];
      double pzprime =  1.0*cos(numi_beam_angle_rad)*FSPartPz[i] + sin(numi_beam_angle_rad)*FSPartPy[i];
      double pSquare = pow(FSPartPx[i],2) + pow(FSPartPy[i],2) + pow(FSPartPz[i],2);
      double protonTheta = acos( pzprime / sqrt(pSquare) );
      protonTheta *= 180./3.14159;

      if (protonP>450 && protonP<1200 && (protonTheta<70 || protonTheta>110)){
	hasSignalProton = true;
      }
    }
    else if (abs(FSPartPDG[i]) == 211 || abs(FSPartPDG[i]) == 321 || abs(FSPartPDG[i]) == 311 || abs(FSPartPDG[i]) == 130 || abs(FSPartPDG[i]) == 111){
      hasMeson = true;
    }
    else if (abs(FSPartPDG[i]) == 22 && FSPartE[i] > 10){
      hasPhoton = true;
    }
  }

  if(abs(mc_incoming)==12 && current==1 && hasSignalProton && !hasMeson && !hasPhoton && electron_4vec[3]>2500) { return true; }
  return false;

}
  // Override this method from the base class to decide what events to
  // include in this selection
  virtual bool passesCuts(PlotUtils::ChainWrapper& chw, int entry)
  {
    //if((int)chw.GetValue("mc_incoming", entry)!=14) return false;
    //if((int)chw.GetValue("mc_current", entry)!=1) return false;
    if(!isNuEQELikeSignal  ( chw, entry ) ) return false;
    
    return true;
  }
};

int main(const int argc, const char** argv)
{
  //Read a playlist file from the command line
  if(argc != 2)
  {
    std::cerr << "Expected exactly 1 command line argument, but got " << argc - 1 << ".\n\n"
              << "USAGE: runXSecLooper <MCPlaylist.txt>\n\n"
              << "MCPlaylist.txt shall contain one .root file per line that has a Truth tree in it.\n"
              << "This program returns 0 when it suceeds.  It produces a .root file with GENIEXSECEXTRACT in its name.\n";
    return 1;
  }

  const std::string playlistFile = argv[1]; //argv[0] is the name of the executable

  // Create the XSecLooper and tell it the input files
  // Inputs should be the merged ntuples:
  XSecLooper loop(playlistFile.c_str());

  // Tell the XSecLooper which neutrino type we're considering (mandatory)
  loop.setNuPDG(12);

  // Setting the number of Universes in the GENIE error band (default 100, put 0 if you do not want to include the universes)
  loop.setNumUniv(0); 
  loop.setFiducial(5980, 8422);

  // Add the differential cross sections
  double lepton_E_edges[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,20};
  int lepton_E_nbins = (sizeof(lepton_E_edges) / sizeof(lepton_E_edges[0])) - 1;
  double E_avail_edges[] = {0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.8, 1, 1.2};
  int E_avail_nbins = (sizeof(E_avail_edges) / sizeof(E_avail_edges[0])) - 1;
  double E_nu_edges[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,20};
  int E_nu_nbins = (sizeof(E_nu_edges) / sizeof(E_nu_edges[0])) - 1;
  double lepton_pt_edges[] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
  int lepton_pt_nbins =  (sizeof(lepton_pt_edges) / sizeof(lepton_pt_edges[0])) - 1;
  double lepton_pz_edges[] = {0, 3, 4.5, 5.5, 6.5, 8, 12, 30};
  int lepton_pz_nbins = (sizeof(lepton_pz_edges) / sizeof(lepton_pz_edges[0])) - 1;
  double lepton_theta_edges[] = {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40};
  int lepton_theta_nbins = (sizeof(lepton_theta_edges) / sizeof(lepton_theta_edges[0])) - 1;

  double proton_p_edges[] = {0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375};
  int proton_p_nbins = (sizeof(proton_p_edges) / sizeof(proton_p_edges[0])) - 1;
  double proton_pt_edges[] = {0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.1};
  int proton_pt_nbins = (sizeof(proton_pt_edges) / sizeof(proton_pt_edges[0])) - 1;
  double proton_theta_edges[] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180};
  int proton_theta_nbins = (sizeof(proton_theta_edges) / sizeof(proton_theta_edges[0])) - 1;
  double proton_KE_edges[] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
  int proton_KE_nbins = (sizeof(proton_KE_edges) / sizeof(proton_KE_edges[0])) - 1;

  
  double delta_pt_edges[] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.4};
  int delta_pt_nbins = (sizeof(delta_pt_edges) / sizeof(delta_pt_edges[0])) - 1;
  double delta_ptx_edges[] = {-2.0, -0.7, -0.4, -0.2, 0.0, 0.2, 0.4, 0.7, 2.0};
  int delta_ptx_nbins = (sizeof(delta_ptx_edges) / sizeof(delta_ptx_edges[0])) - 1;
  double delta_pty_edges[] = {-2.0, -1.1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 1.3};
  int delta_pty_nbins = (sizeof(delta_pty_edges) / sizeof(delta_pty_edges[0])) - 1;
  double delta_pl_edges[] = {-0.25, 0.00, 0.10, 0.20, 0.30, 0.40, 0.60};
  int delta_pl_nbins = (sizeof(delta_pl_edges) / sizeof(delta_pl_edges[0])) - 1;
  double neutron_momentum_edges[] = {0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 1.0};
  int neutron_momentum_nbins = (sizeof(neutron_momentum_edges) / sizeof(neutron_momentum_edges[0])) - 1;
  double alpha_edges[] = {0,20,40,60,80,100,120,140,160,180};
  int alpha_nbins =  (sizeof(alpha_edges) / sizeof(alpha_edges[0])) - 1; 
  double phi_edges[] = {0,10,30,50,80,120,180};
  int phi_nbins =  (sizeof(phi_edges) / sizeof(phi_edges[0])) - 1; 
  
  //=========================================
  // Lepton variables ( + E_avail and E_nu )
  //=========================================
  MinModDepCCQEXSec* ds_dElep = new MinModDepCCQEXSec("E_lep");
  ds_dElep->setBinEdges(lepton_E_nbins, lepton_E_edges);
  ds_dElep->setVariable(XSec::kELep);
  ds_dElep->setIsFluxIntegrated(true);
  ds_dElep->setDimension(1);
  ds_dElep->setFluxIntLimits(0.0, 100.0);
  ds_dElep->setNormalizationType(XSec::kPerNucleon);  
  ds_dElep->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dElep);

  MinModDepCCQEXSec* ds_dEavail = new MinModDepCCQEXSec("E_avail");
  ds_dEavail->setBinEdges(E_avail_nbins, E_avail_edges);
  ds_dEavail->setVariable(XSec::kEAvail);
  ds_dEavail->setIsFluxIntegrated(true);
  ds_dEavail->setDimension(1);
  ds_dEavail->setFluxIntLimits(0.0, 100.0);
  ds_dEavail->setNormalizationType(XSec::kPerNucleon);  
  ds_dEavail->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dEavail);

  MinModDepCCQEXSec* ds_dEnu = new MinModDepCCQEXSec("E_nu");
  ds_dEnu->setBinEdges(E_nu_nbins, E_nu_edges);
  ds_dEnu->setVariable(XSec::kENu);
  ds_dEnu->setIsFluxIntegrated(true);
  ds_dEnu->setDimension(1);
  ds_dEnu->setFluxIntLimits(0.0, 100.0);
  ds_dEnu->setNormalizationType(XSec::kPerNucleon);  
  ds_dEnu->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dEnu);    

  MinModDepCCQEXSec* ds_dpT = new MinModDepCCQEXSec("Lepton_Pt");
  ds_dpT->setBinEdges(lepton_pt_nbins, lepton_pt_edges);
  ds_dpT->setVariable(XSec::kPTLep);
  ds_dpT->setIsFluxIntegrated(true);
  ds_dpT->setDimension(1);
  ds_dpT->setFluxIntLimits(0.0, 100.0);
  ds_dpT->setNormalizationType(XSec::kPerNucleon);  
  ds_dpT->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT);   
 
  MinModDepCCQEXSec* ds_dpZ = new MinModDepCCQEXSec("Lepton_Pl");
  ds_dpZ->setBinEdges(lepton_pz_nbins, lepton_pz_edges);
  ds_dpZ->setVariable(XSec::kPZLep);
  ds_dpZ->setIsFluxIntegrated(true);
  ds_dpZ->setDimension(1);
  ds_dpZ->setFluxIntLimits(0.0, 100.0);
  ds_dpZ->setNormalizationType(XSec::kPerNucleon);  
  ds_dpZ->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpZ);

  MinModDepCCQEXSec* ds_dThetaLep = new MinModDepCCQEXSec("Theta_lep");
  ds_dThetaLep->setBinEdges(lepton_theta_nbins, lepton_theta_edges);
  ds_dThetaLep->setVariable(XSec::kThetaLep);
  ds_dThetaLep->setIsFluxIntegrated(true);
  ds_dThetaLep->setDimension(1);
  ds_dThetaLep->setFluxIntLimits(0.0, 100.0);
  ds_dThetaLep->setNormalizationType(XSec::kPerNucleon);  
  ds_dThetaLep->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dThetaLep);    

  //=========================================
  // Proton Variables
  //=========================================
  MinModDepCCQEXSec* ds_dprotonP = new MinModDepCCQEXSec("Proton_p");
  ds_dprotonP->setBinEdges(proton_p_nbins, proton_p_edges);
  ds_dprotonP->setVariable(XSec::kTprotonmomentum);
  ds_dprotonP->setIsFluxIntegrated(true);
  ds_dprotonP->setDimension(1);
  ds_dprotonP->setFluxIntLimits(0.0, 100.0);
  ds_dprotonP->setNormalizationType(XSec::kPerNucleon);  
  ds_dprotonP->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dprotonP);    

  MinModDepCCQEXSec* ds_dprotonPt = new MinModDepCCQEXSec("Proton_Pt");
  ds_dprotonPt->setBinEdges(proton_pt_nbins, proton_pt_edges);
  ds_dprotonPt->setVariable(XSec::kTprotonPt);
  ds_dprotonPt->setIsFluxIntegrated(true);
  ds_dprotonPt->setDimension(1);
  ds_dprotonPt->setFluxIntLimits(0.0, 100.0);
  ds_dprotonPt->setNormalizationType(XSec::kPerNucleon);  
  ds_dprotonPt->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dprotonPt);    

  MinModDepCCQEXSec* ds_dprotonTheta = new MinModDepCCQEXSec("Theta_p");
  ds_dprotonTheta->setBinEdges(proton_theta_nbins, proton_theta_edges);
  ds_dprotonTheta->setVariable(XSec::kTprotontheta);
  ds_dprotonTheta->setIsFluxIntegrated(true);
  ds_dprotonTheta->setDimension(1);
  ds_dprotonTheta->setFluxIntLimits(0.0, 100.0);
  ds_dprotonTheta->setNormalizationType(XSec::kPerNucleon);  
  ds_dprotonTheta->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dprotonTheta);    

  MinModDepCCQEXSec* ds_dprotonKE = new MinModDepCCQEXSec("Proton_T");
  ds_dprotonKE->setBinEdges(proton_KE_nbins, proton_KE_edges);
  ds_dprotonKE->setVariable(XSec::kTprotonKE);
  ds_dprotonKE->setIsFluxIntegrated(true);
  ds_dprotonKE->setDimension(1);
  ds_dprotonKE->setFluxIntLimits(0.0, 100.0);
  ds_dprotonKE->setNormalizationType(XSec::kPerNucleon);  
  ds_dprotonKE->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dprotonKE);   
  
  //=========================================
  // Single TKI variables
  //=========================================
  MinModDepCCQEXSec* ds_ddpT = new MinModDepCCQEXSec("DeltaPt");
  ds_ddpT->setBinEdges(delta_pt_nbins, delta_pt_edges);
  ds_ddpT->setVariable(XSec::kTdpt);
  ds_ddpT->setIsFluxIntegrated(true);
  ds_ddpT->setDimension(1);
  ds_ddpT->setFluxIntLimits(0.0, 100.0);
  ds_ddpT->setNormalizationType(XSec::kPerNucleon);  
  ds_ddpT->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_ddpT);

  MinModDepCCQEXSec* ds_ddpTX = new MinModDepCCQEXSec("DeltaPtX");
  ds_ddpTX->setBinEdges(delta_ptx_nbins, delta_ptx_edges);
  ds_ddpTX->setVariable(XSec::kTdptx);
  ds_ddpTX->setIsFluxIntegrated(true);
  ds_ddpTX->setDimension(1);
  ds_ddpTX->setFluxIntLimits(0.0, 100.0);
  ds_ddpTX->setNormalizationType(XSec::kPerNucleon);  
  ds_ddpTX->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_ddpTX);   

  MinModDepCCQEXSec* ds_ddpTY = new MinModDepCCQEXSec("DeltaPtY");
  ds_ddpTY->setBinEdges(delta_pty_nbins, delta_pty_edges);
  ds_ddpTY->setVariable(XSec::kTdpty);
  ds_ddpTY->setIsFluxIntegrated(true);
  ds_ddpTY->setDimension(1);
  ds_ddpTY->setFluxIntLimits(0.0, 100.0);
  ds_ddpTY->setNormalizationType(XSec::kPerNucleon);  
  ds_ddpTY->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_ddpTY);

  MinModDepCCQEXSec* ds_ddpL = new MinModDepCCQEXSec("DeltaPl");
  ds_ddpL->setBinEdges(delta_pl_nbins, delta_pl_edges);
  ds_ddpL->setVariable(XSec::kTdpl);
  ds_ddpL->setIsFluxIntegrated(true);
  ds_ddpL->setDimension(1);
  ds_ddpL->setFluxIntLimits(0.0, 100.0);
  ds_ddpL->setNormalizationType(XSec::kPerNucleon);  
  ds_ddpL->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_ddpL);     

  MinModDepCCQEXSec* ds_dnp = new MinModDepCCQEXSec("P_n"); //struck neutron momentum
  ds_dnp->setBinEdges(neutron_momentum_nbins, neutron_momentum_edges);
  ds_dnp->setVariable(XSec::kTneutronmomentum);
  ds_dnp->setIsFluxIntegrated(true);
  ds_dnp->setDimension(1);
  ds_dnp->setFluxIntLimits(0.0, 100.0);
  ds_dnp->setNormalizationType(XSec::kPerNucleon);  
  ds_dnp->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dnp);   
  
  MinModDepCCQEXSec* ds_dAlphapT = new MinModDepCCQEXSec("AlphaPt"); //boosting angle
  ds_dAlphapT->setBinEdges(alpha_nbins, alpha_edges);
  ds_dAlphapT->setVariable(XSec::kTdalphat);
  ds_dAlphapT->setIsFluxIntegrated(true);
  ds_dAlphapT->setDimension(1);
  ds_dAlphapT->setFluxIntLimits(0.0, 100.0);
  ds_dAlphapT->setNormalizationType(XSec::kPerNucleon);  
  ds_dAlphapT->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dAlphapT);

  MinModDepCCQEXSec* ds_dPhiT = new MinModDepCCQEXSec("PhiPt"); //acoplanarity
  ds_dPhiT->setBinEdges(phi_nbins, phi_edges);
  ds_dPhiT->setVariable(XSec::kTdphit);
  ds_dPhiT->setIsFluxIntegrated(true);
  ds_dPhiT->setDimension(1);
  ds_dPhiT->setFluxIntLimits(0.0, 100.0);
  ds_dPhiT->setNormalizationType(XSec::kPerNucleon);  
  ds_dPhiT->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dPhiT);
  
  loop.runLoop();

  // Get the output histograms and save them to file
  std::cout << "CARLOS TEST" << std::endl;  
  std::cout << playlistFile << std::endl;
  std::cout << playlistFile.rfind("/")+1 << std::endl;
  std::cout << playlistFile.find(".") << std::endl;
  std::cout << playlistFile.substr(playlistFile.rfind("/")+1, playlistFile.find(".")) << std::endl;
  std::string geniefilename =  "GENIEXSECEXTRACT_" + playlistFile.substr(playlistFile.rfind("/")+1, playlistFile.find(".")) + ".root";
  TFile fout(geniefilename.c_str(), "RECREATE");
  for(uint i=0; i<loop.getXSecs().size(); ++i)
  {
    loop.getXSecs()[i]->getXSecHist()->Write();
    loop.getXSecs()[i]->getEvRateHist()->Write();
  }

  return 0;
}
