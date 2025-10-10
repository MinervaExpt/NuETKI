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

bool isNuEQELikeSignal( ChainWrapper& chw, int entry )
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

  if(abs(mc_incoming)==12 && current==1 && hasSignalProton && !hasMeson && !hasPhoton && electron_4vec[3]>2500) return true;
  return false;

}
  // Override this method from the base class to decide what events to
  // include in this selection
  virtual bool passesCuts(ChainWrapper& chw, int entry)
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

  // Add the differential cross section dsigma/ds_dpT
  double pt_edges[] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};
  int pt_nbins = 15; 
  double alpha_edges[] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180};
  int alpha_nbins = 18; 
  double SumTp_edges[] = {0,100,200,300,400,500,600,700,800,900,1000};
  int SumTp_nbins = 10; 
  
  // Flux-integrated over the range 0.0 to 100.0 GeV
  MinModDepCCQEXSec* ds_dpT = new MinModDepCCQEXSec("pT");
  ds_dpT->setBinEdges(pt_nbins, pt_edges);
  ds_dpT->setVariable(XSec::kPTLep);
  ds_dpT->setIsFluxIntegrated(true);
  ds_dpT->setDimension(1);
  ds_dpT->setFluxIntLimits(0.0, 100.0);
  ds_dpT->setNormalizationType(XSec::kPerNucleon);  
  ds_dpT->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT);
  
  
  MinModDepCCQEXSec* ds_dDeltapT = new MinModDepCCQEXSec("DeltapT");
  ds_dDeltapT->setBinEdges(pt_nbins, pt_edges);
  ds_dDeltapT->setVariable(XSec::kTdpt);
  ds_dDeltapT->setIsFluxIntegrated(true);
  ds_dDeltapT->setDimension(1);
  ds_dDeltapT->setFluxIntLimits(0.0, 100.0);
  ds_dDeltapT->setNormalizationType(XSec::kPerNucleon);  
  ds_dDeltapT->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dDeltapT);

  MinModDepCCQEXSec* ds_dAlphapT = new MinModDepCCQEXSec("AlphapT");
  ds_dAlphapT->setBinEdges(alpha_nbins, alpha_edges);
  ds_dAlphapT->setVariable(XSec::kTdalphat);
  ds_dAlphapT->setIsFluxIntegrated(true);
  ds_dAlphapT->setDimension(1);
  ds_dAlphapT->setFluxIntLimits(0.0, 100.0);
  ds_dAlphapT->setNormalizationType(XSec::kPerNucleon);  
  ds_dAlphapT->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dAlphapT);

  MinModDepCCQEXSec* ds_SumT_p = new MinModDepCCQEXSec("SumTp");
  ds_SumT_p->setBinEdges(SumTp_nbins, SumTp_edges);
  ds_SumT_p->setVariable(XSec::kSumTP);
  ds_SumT_p->setIsFluxIntegrated(true);
  ds_SumT_p->setDimension(1);
  ds_SumT_p->setFluxIntLimits(0.0, 100.0);
  ds_SumT_p->setNormalizationType(XSec::kPerNucleon);  
  ds_SumT_p->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_SumT_p);
  
  loop.runLoop();

  // Get the output histograms and save them to file
  string geniefilename =  "GENIEXSECEXTRACT_" + playlistFile.substr(playlistFile.rfind("/")+1, playlistFile.find(".")) + ".root";
  TFile fout(geniefilename.c_str(), "RECREATE");
  for(uint i=0; i<loop.getXSecs().size(); ++i)
  {
    loop.getXSecs()[i]->getXSecHist()->Write();
    loop.getXSecs()[i]->getEvRateHist()->Write();
  }

  return 0;
}
