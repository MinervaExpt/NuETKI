//Whether or not to save a TTree of my selected events (adds a lot of runtime if true)
bool write_tree = false; //true;
int n_sideband_cuts; //have to do this because the number of sideband cuts in Cutter is private and not accessible from here :/

#define USAGE					\
"\n*** USAGE ***\n"\
"runEventLoop <dataPlaylist.txt> <mcPlaylist.txt> <configFile.yaml>\n\n"\
"*** Explanation ***\n"\
"Reduce MasterAnaDev AnaTuples to event selection histograms to extract a\n"\
"CCQE-Like electron neutrino cross section, adapted from the 2021 MINERvA 101 tutorial.\n\n"\
"*** The Input Files ***\n"\
"Playlist files are plaintext files with 1 file name per line.  Filenames may be\n"\
"xrootd URLs or refer to the local filesystem.  The first playlist file's\n"\
"entries will be treated like data, and the second playlist's entries must\n"\
"have the \"Truth\" tree to use for calculating the efficiency denominator.\n"\
"The yaml config file is used to control output variables, cuts, etc. without\n"\
"needing to recompile.\n\n"\
"*** Output ***\n"\
"Produces a two files, a Data.root and an MC.root file, with\n"\
"all histograms needed for the ExtractCrossSection program also built by this\n"\
"package.  You'll need a .rootlogon.C that loads ROOT object definitions from\n"\
"PlotUtils to access systematics information from these files.\n\n"\
"*** Environment Variables ***\n"\
"Setting up this package appends to PATH and LD_LIBRARY_PATH.  PLOTUTILSROOT,\n"\
"MPARAMFILESROOT, and MPARAMFILES must be set according to the setup scripts in\n"\
"those packages for systematics and flux reweighters to function.\n"\
"If MNV101_SKIP_SYST is defined at all, output histograms will have no error bands.\n"\
"This is useful for debugging the CV and running warping studies.\n\n"\
"*** Return Codes ***\n"\
"0 indicates success.  All histograms are valid only in this case.  Any other\n"\
"return code indicates that histograms should not be used.  Error messages\n"\
"about what went wrong will be printed to stderr.  So, they'll end up in your\n"\
"terminal, but you can separate them from everything else with something like:\n"\
"\"runEventLoop data.txt mc.txt 2> errors.txt\"\n"

enum ErrorCodes
{
  success = 0,
  badCmdLine = 1,
  badInputFile = 2,
  badFileRead = 3,
  badOutputFile = 4
};

//PlotUtils includes
//No junk from PlotUtils please!  I already
//know that MnvH1D does horrible horrible things.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"

//Includes from this package
#include "event/CVUniverse.h"
#include "event/CCNuEEvent.h"
#include "systematics/Systematics.h"
#include "cuts/MaxPzMu.h"
#include "util/Variable.h"
#include "util/Variable2D.h"
#include "util/VariableRegistry.h"
#include "util/ReweighterRegistry.h"
#include "util/GetFluxIntegral.h"
#include "util/GetPlaylist.h"
#include "util/Binning.h"
#include "cuts/SignalDefinition.h"
#include "cuts/q3RecoCut.h"
#include "cuts/NuETKICuts.h"
#include "cuts/NuETKISignal.h"
#include "cuts/CutRegistry.h"

//Studies and sidebands
#include "studies/Study.h"
#include "studies/MeanFrontDEDXSideband.h"
#include "studies/MichelSideband.h"
#include "studies/NIsoClusSideband.h"
#include "studies/NPiSideband.h"

//PlotUtils includes
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/MacroUtil.h"
#include "PlotUtils/MnvPlotter.h"
//#include "PlotUtils/CCInclusiveCuts.h"
//#include "PlotUtils/CCInclusiveSignal.h"  //taking these out to replace with my own set of cuts, NuETKI
#include "PlotUtils/CrashOnROOTMessage.h" //Sets up ROOT's debug callbacks by itself
#include "PlotUtils/Cutter.h"
#include "PlotUtils/Model.h"
#include "PlotUtils/FluxAndCVReweighter.h"
#include "PlotUtils/GENIEReweighter.h"
#include "PlotUtils/LowRecoil2p2hReweighter.h"
#include "PlotUtils/RPAReweighter.h"
#include "PlotUtils/MINOSEfficiencyReweighter.h"
#include "PlotUtils/LowQ2PiReweighter.h"
#include "PlotUtils/AMUDISReweighter.h"
#include "PlotUtils/FSIReweighter.h"
#include "PlotUtils/TargetUtils.h"
#pragma GCC diagnostic pop

//ROOT includes
#include "TParameter.h"
#include "Math/Vector3D.h"
#include "Math/AxisAngle.h"

//c++ includes
#include <iostream>
#include <cstdlib> //getenv()
#include <ctime>
#include <sstream>
#include <iomanip>
#include <chrono>

//yaml
#include <yaml-cpp/yaml.h>


OutputTreeManager g_OutputTreeManager; //defining my output tree manager (struct declaration and full def in CCNuEEvent.h)

//added these for grid time tracking
auto getTimeStamp() {
    return std::chrono::high_resolution_clock::now();
}
auto printElapsed(const std::string& label, 
                  std::chrono::high_resolution_clock::time_point start,
                  std::chrono::high_resolution_clock::time_point end) {
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << label << ": " << elapsed << " ms\n";
}

//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFillEventSelection(
    PlotUtils::ChainWrapper* chain,
    std::map<std::string, std::vector<CVUniverse*> > error_bands,
    std::vector<Variable*> vars,
    std::vector<Variable2D*> vars2D,
    std::vector<Study*> studies,
    PlotUtils::Cutter<CVUniverse, CCNuEEvent>& michelcuts,
    PlotUtils::Model<CVUniverse, CCNuEEvent>& model)
{
  assert(!error_bands["cv"].empty() && "\"cv\" error band is empty!  Can't set Model weight.");
  auto& cvUniv = error_bands["cv"].front();

  std::cout << "Starting MC reco loop...\n";
  
  const int nEntries = chain->GetEntries();
  //chain->GetChain()->SetBranchStatus("*", 0);
  for (Long64_t i=0; i<nEntries; ++i)
  {
    //if(i%1000==0) std::cout << i << " / " << nEntries << "\r" << std::flush;

    CCNuEEvent cvEvent;
    cvUniv->SetEntry(i);
    model.SetEntry(*cvUniv, cvEvent);
    const double cvWeight = model.GetWeight(*cvUniv, cvEvent);

    //=========================================
    // Systematics loop(s)
    //=========================================
    for (auto band : error_bands)
    {
      std::vector<CVUniverse*> error_band_universes = band.second;
      //int count = 0; //only used for testing/printing universe # within a band
      for (auto universe : error_band_universes)
      {
        CCNuEEvent myevent; // make sure your event is inside the error band loop. 
        // Tell the Event which entry in the TChain it's looking at
        universe->SetEntry(i);
	myevent.primaryProtonIndex = universe->GetHighestEnergySignalProtonIndex();
        // This is where you would Access/create a Michel
        //weight is ignored in isMCSelected() for all but the CV Universe.
	//if ths event doesnt pass selection cuts, skip

	//std::cout << "band/universe: " << band.first << "/" << count <<  ", Model Weight = " << cvWeight << std::endl;
	//count++;

	//isMCSelected returns a bitset which sideband cuts are passed. It is all 1s other than specifically the sideband cuts, ie we lose info about specific precuts (which is fine). 
	//So if filling sidebands I don't want to check !results.all() because the sideband result will always be 1
	//Instead I'll make a bitmask (basically a comparison bitset which I'll use to check all bits except the first one, which can be whatever
	//and then I'll use that first cut result to fill my sideband (if 0, goes into sideband, if 1, goes into selected events)
	std::bitset<64> cut_results = michelcuts.isMCSelected(*universe, myevent, cvWeight);  

	std::bitset<64> mask;
	mask.set(); //set all bits in mask to 1
	mask &= ~std::bitset<64>((1ULL << n_sideband_cuts) - 1); //very fancy bit manipulation from chat gpt, just sets the first (or more accurately, least significant) n bits to 0

	if ((cut_results & mask) != mask) continue; //skips events if it doesn't pass all cuts OTHER than sideband cuts

	//save the results of the sideband cuts only as a vector, to pass to my output tree manager
	std::vector<int> sideband_results;
	for (int i = 0; i < n_sideband_cuts; i++) {
	  sideband_results.push_back((int)cut_results[i]);
	}

	if (write_tree){
	  universe->GetTree()->GetTree()->GetEntry(i); //how is this different than universe->SetEntry(i)? not sure but it definitely is...
	  myevent.entryNumber = i;
	}
	g_OutputTreeManager.sidebandCutResults = sideband_results;
	g_OutputTreeManager.cv_weight = cvWeight;
	
	const double weight = model.GetWeight(*universe, myevent); //Only calculate the per-universe weight for events that will actually use it.	
	const bool isSignal = michelcuts.isSignal(*universe, weight);

	int bkgd_ID;
	int sig_ID;
	if (isSignal){
	  g_OutputTreeManager.selectionCategory = -999;
	  sig_ID = 4; //sig_ID == 4 -> NuE CCQELike other, mostly (all?) coherent events in my signal sample
	  if (universe->GetInteractionType()==1){ sig_ID = 0; } //NuE CCQELike QE
	  else if (universe->GetInteractionType()==2){ sig_ID = 1; } //NuE CCQELike RES
	  else if (universe->GetInteractionType()==3){ sig_ID = 2; } //NuE CCQELike DIS
	  else if (universe->GetInteractionType()==8){ sig_ID = 3; } //NuE CCQELike 2p2h
	}
	else {
	  bkgd_ID = -1; //this is "other" selectionCategory
	  //NueCC with final state mesons (mostly pions)
	  //if (abs(universe->GetTruthNuPDG())==12 && universe->GetCurrent()==1 && universe->GetHasFSMeson()) { bkgd_ID = 0; }
	  if (abs(universe->GetTruthNuPDG())==12 && universe->GetCurrent()==1 && universe->GetHasSinglePiPlus()) { bkgd_ID = 0; }
	  else if (abs(universe->GetTruthNuPDG())==12 && universe->GetCurrent()==1 && universe->GetHasSinglePiMinus()) { bkgd_ID = 1; }
	  else if (abs(universe->GetTruthNuPDG())==12 && universe->GetCurrent()==1 && universe->GetHasSinglePiZero()) { bkgd_ID = 2; }
	  else if (abs(universe->GetTruthNuPDG())==12 && universe->GetCurrent()==1 && universe->GetHasMultiplePions()) { bkgd_ID = 3; }

	  //Other NueCC. mostly just consists of completely signal-like events but protons are outside of kinematic reqs
	  else if (abs(universe->GetTruthNuPDG())==12 && universe->GetCurrent()==1) { bkgd_ID = 4; }
	  //Other NC pi0, all "other" NC i get should have pi0 so if this is small and the "other" selectionCategory is large that's a problem
	  else if (universe->GetCurrent()==2 && universe->GetHasFSPi0()==1) { bkgd_ID = 5; }
	  //same logic here with cc numu pi0
	  else if (abs(universe->GetTruthNuPDG())==14 && universe->GetCurrent()==1 && universe->GetHasFSPi0()==1) { bkgd_ID = 6; }
	  
	  //If I need to use Hang's bkgd categories (cc numu, NC Coh, other NC, nu + e) use these (and also use the corresponding Variable.h)	  
	  //int bkgd_ID = -1;
	  //if (abs(universe->GetTruthNuPDG())==14 && universe->GetCurrent()==1){ bkgd_ID = 0; } //cc numu
	  //else if (universe->GetCurrent()==2 && universe->GetInteractionType()==4){ bkgd_ID = 1; } //NC (current=2) and Coh (inttype = 4)
	  //else if (universe->GetCurrent()==2){ bkgd_ID = 2; } //other NC
	  //else if (universe->GetInteractionType()==7){ bkgd_ID = 3; } //nu + e elastic	 
	  
	  g_OutputTreeManager.selectionCategory = bkgd_ID;
	}
	//Here we fill sidebands
	for(auto& study: studies){ study->SelectedMC(*universe, myevent, weight); }
	g_OutputTreeManager.Fill("MasterAnaDev", myevent.entryNumber);

	if ((cut_results & ~mask) != ~mask) continue; //Check the remaining sideband cuts, everything that passes this has now passed ALL cuts

	//print out events if I want to
	if (band.first == "cv" && isSignal) {
	  std::cout << "signal event: " << universe->GetInt("mc_run") << " | " << universe->GetInt("mc_subrun") << " | " << universe->GetInt("mc_nthEvtInFile")+1 << ", total cvWeight = " << cvWeight << std::endl;
	}
	
        for(auto& var: vars) var->selectedMCReco->FillUniverse(universe, var->GetRecoValue(*universe), weight); //"Fake data" for closure 
		
        if(isSignal)
        {
          //for(auto& study: studies) study->SelectedSignal(*universe, myevent, weight);
          for(auto& var: vars)
          {
            //Cross section components
            var->efficiencyNumerator->FillUniverse(universe, var->GetTrueValue(*universe), weight);
            var->migration->FillUniverse(universe, var->GetRecoValue(*universe), var->GetTrueValue(*universe), weight);
            var->selectedSignalReco->FillUniverse(universe, var->GetRecoValue(*universe), weight); //Efficiency numerator in reco variables.  Useful for warping studies.
	    //My additional signal breakdown, not used for any cross section calculation but useful for plotting
	    (*var->m_signalBreakdown)[sig_ID].FillUniverse(universe, var->GetRecoValue(*universe), weight);
          }

          for(auto& var: vars2D)
          {
            var->efficiencyNumerator->FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
          }
        }
	else{ //The actual logic for which category an event gets put into happens above (so that we can save it in the tuple & sidebands)
          for(auto& var: vars) (*var->m_backgroundHists)[bkgd_ID].FillUniverse(universe, var->GetRecoValue(*universe), weight);
          for(auto& var: vars2D) (*var->m_backgroundHists)[bkgd_ID].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
	}
      } // End band's universe loop
    } // End Band loop
  } //End entries loop
  
  std::cout << "Finished MC reco loop.\n";
}

void LoopAndFillData( PlotUtils::ChainWrapper* data,
			        std::vector<CVUniverse*> data_band,
				std::vector<Variable*> vars,
                                std::vector<Variable2D*> vars2D,
                                std::vector<Study*> studies,
				PlotUtils::Cutter<CVUniverse, CCNuEEvent>& michelcuts)

{
  std::cout << "Starting data loop...\n";
  const int nEntries = data->GetEntries();
  for (int i=0; i<data->GetEntries(); ++i) {
    for (auto universe : data_band) {
      universe->SetEntry(i);
      //if(i%1000==0) std::cout << i << " / " << nEntries << "\r" << std::flush;
      CCNuEEvent myevent;
      
      //if (!michelcuts.isDataSelected(*universe, myevent).all()) continue;
      std::bitset<64> cut_results = michelcuts.isDataSelected(*universe, myevent);  
      std::bitset<64> mask;
      mask.set(); //set all bits in mask to 1
      mask &= ~std::bitset<64>((1ULL << n_sideband_cuts) - 1); //very fancy bit manipulation from chat gpt, just sets the first (or more accurately, least significant) n bits to 0
      if ((cut_results & mask) != mask) continue; //skips events if it doesn't pass all cuts OTHER than sideband cuts

      //fill sidebands
      for(auto& study: studies) study->SelectedData(*universe, myevent, 1); 

      if ((cut_results & ~mask) != ~mask) continue; //Check the remaining sideband cuts, everything that passes this has now passed ALL cuts       

      for(auto& var: vars)
      {
        var->dataHist->FillUniverse(universe, var->GetRecoValue(*universe, myevent.m_idx), 1);
      }

      for(auto& var: vars2D)
      {
        var->dataHist->FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), 1);
      }
    }
  }
  std::cout << "Finished data loop.\n";
}

void LoopAndFillEffDenom( PlotUtils::ChainWrapper* truth,
    				std::map<std::string, std::vector<CVUniverse*> > truth_bands,
    				std::vector<Variable*> vars,
                                std::vector<Variable2D*> vars2D,
    				PlotUtils::Cutter<CVUniverse, CCNuEEvent>& michelcuts,
                                PlotUtils::Model<CVUniverse, CCNuEEvent>& model)
{
  assert(!truth_bands["cv"].empty() && "\"cv\" error band is empty!  Could not set Model entry.");
  auto& cvUniv = truth_bands["cv"].front();

  std::cout << "Starting efficiency denominator loop...\n";
  const int nEntries = truth->GetEntries();
  for (int i=0; i<nEntries; ++i)
  {
    //if(i%1000==0) std::cout << i << " / " << nEntries << "\r" << std::flush;

    CCNuEEvent cvEvent;
    cvUniv->SetEntry(i);
    model.SetEntry(*cvUniv, cvEvent);
    const double cvWeight = model.GetWeight(*cvUniv, cvEvent);

    //=========================================
    // Systematics loop(s)
    //=========================================
    for (auto band : truth_bands)
    {
      std::vector<CVUniverse*> truth_band_universes = band.second;
      for (auto universe : truth_band_universes)
      {
        CCNuEEvent myevent; //Only used to keep the Model happy

        // Tell the Event which entry in the TChain it's looking at
        universe->SetEntry(i);

        if (!michelcuts.isEfficiencyDenom(*universe, cvWeight)) continue; //Weight is ignored for isEfficiencyDenom() in all but the CV universe 
        const double weight = model.GetWeight(*universe, myevent); //Only calculate the weight for events that will use it

        //Fill efficiency denominator now: 
        for(auto var: vars)
        {
	  //if (band.first == "cv") { std::cout << "Event i=" << i << ": " << universe->GetInt("mc_run") << " | " << universe->GetInt("mc_subrun") << " | " << universe->GetInt("mc_nthEvtInFile") << ", val = " << var->GetTrueValue(*universe) << ", cvweight = " << cvWeight << std::endl; }
	  var->efficiencyDenominator->FillUniverse(universe, var->GetTrueValue(*universe), weight);
        }

        for(auto var: vars2D)
        {
          var->efficiencyDenominator->FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
        }
      }
    }
  }
  std::cout << "Finished efficiency denominator loop.\n";
}

//Returns false if recoTreeName could not be inferred
bool inferRecoTreeNameAndCheckTreeNames(const std::string& mcPlaylistName, const std::string& dataPlaylistName, std::string& recoTreeName)
{
  const std::vector<std::string> knownTreeNames = {"Truth", "Meta"};
  bool areFilesOK = false;

  std::ifstream playlist(mcPlaylistName);
  std::string firstFile = "";
  playlist >> firstFile;
  auto testFile = TFile::Open(firstFile.c_str());
  if(!testFile)
  {
    std::cerr << "Failed to open the first MC file at " << firstFile << "\n";
    return false;
  }

  //Does the MC playlist have the Truth tree?  This is needed for the efficiency denominator.
  const auto truthTree = testFile->Get("Truth");
  if(truthTree == nullptr || !truthTree->IsA()->InheritsFrom(TClass::GetClass("TTree")))
  {
    std::cerr << "Could not find the \"Truth\" tree in MC file named " << firstFile << "\n";
    return false;
  }

  //Figure out what the reco tree name is
  for(auto key: *testFile->GetListOfKeys())
  {
    if(static_cast<TKey*>(key)->ReadObj()->IsA()->InheritsFrom(TClass::GetClass("TTree"))
       && std::find(knownTreeNames.begin(), knownTreeNames.end(), key->GetName()) == knownTreeNames.end())
    {
      recoTreeName = key->GetName();
      areFilesOK = true;
    }
  }
  delete testFile;
  testFile = nullptr;

  //Make sure the data playlist's first file has the same reco tree
  playlist.open(dataPlaylistName);
  playlist >> firstFile;
  testFile = TFile::Open(firstFile.c_str());
  if(!testFile)
  {
    std::cerr << "Failed to open the first data file at " << firstFile << "\n";
    return false;
  }

  const auto recoTree = testFile->Get(recoTreeName.c_str());
  if(recoTree == nullptr || !recoTree->IsA()->InheritsFrom(TClass::GetClass("TTree")))
  {
    std::cerr << "Could not find the \"" << recoTreeName << "\" tree in data file named " << firstFile << "\n";
    return false;
  }

  return areFilesOK;
}

//==============================================================================
// Main
//==============================================================================
int main(const int argc, const char** argv)
{
  TH1::AddDirectory(false);

  std::chrono::high_resolution_clock::time_point t_start, t_mc_start, t_eff_start, t_data_start, t_data_end, t_end;
  t_start = getTimeStamp();
  
  //Validate input.
  //I expect a data playlist file name and an MC playlist file name, and an optional 3rd config file (otherwise use default)
  const int nArgsExpected = 2;
  if(argc < nArgsExpected + 1) //argc is the size of argv.  I check for number of arguments + 1 because
                                //argv[0] is always the path to the executable. and the 3rd config option is optional
  {
    std::cerr << "Expected " << nArgsExpected << " arguments, but got " << argc - 1 << "\n" << USAGE << "\n";
    return badCmdLine;
  }
  
  //One playlist must contain only MC files, and the other must contain only data files.
  //Only checking the first file in each playlist because opening each file an extra time
  //remotely (e.g. through xrootd) can get expensive.
  
  const std::string mc_file_list = argv[2],
                    data_file_list = argv[1];

  //Check that necessary TTrees exist in the first file of mc_file_list and data_file_list
  std::string reco_tree_name;
  if(!inferRecoTreeNameAndCheckTreeNames(mc_file_list, data_file_list, reco_tree_name))
  {
    std::cerr << "Failed to find required trees in MC playlist " << mc_file_list << " and/or data playlist " << data_file_list << ".\n" << USAGE << "\n";
    return badInputFile;
  }

  //Get date & time for output file names
  std::stringstream ss_1, ss_2, ss_3;
  std::time_t t = std::time(nullptr);
  std::tm tm = *std::localtime(&t);    
  //ss_1 << "MC_" << std::put_time(&tm, "%b_%d_%H%M") << ".root";
  ss_1 << "MC_" << std::put_time(&tm, "%b_%d") << ".root";
  //ss_2 << "Data_" << std::put_time(&tm, "%b_%d_%H%M") << ".root";
  ss_2 << "Data_" << std::put_time(&tm, "%b_%d") << ".root";
  ss_3 << "TUPLE_MC_" << std::put_time(&tm, "%b_%d") << ".root";
  std::string mc_out_filename = ss_1.str();
  std::string data_out_filename = ss_2.str();
  std::string mc_tuple_out_filename = ss_3.str();


  //MacroUtil member variables are:
  //PlotUtils::ChainWrapper* m_data, m_mc, m_truth
  //string m_plist_string
  //double m_data_pot, m_mc_pot
  //And then it also has some getters: GetDataEntries(), GetMCEntries(), GetTruthEntries(), PrintMacroConfiguration(std::string macro_name = "")
  PlotUtils::MacroUtil options(reco_tree_name, mc_file_list, data_file_list, "minervame1A", true);
  options.m_plist_string = util::GetPlaylist(*options.m_mc, true); //TODO: Put GetPlaylist into PlotUtils::MacroUtil

  //Load YAML config options
  std::string cfgFile = "/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/config/analysis.yaml"; //default config
  if (argc > 3) cfgFile = argv[3]; //override default config if one is provided
  YAML::Node config = YAML::LoadFile(cfgFile);
  auto variables_cfg = config["variables"];  //list of individual variables to write xsec ingredients for
  auto precuts_cfg = config["precuts"];      //list of precuts to apply, checked first
  auto sbcuts_cfg = config["sideband_cuts"]; //list of sideband cuts to apply, checked after filling sidebands
  auto sidebands_cfg = config["sidebands"];  //list of sideband regions to write
  auto model_cfg = config["CVModel"];        //list of reweighters to apply to the CV
  auto syst_cfg = config["systematics"];     //true or false (for now), either add all of them or none
  auto tree_cfg = config["write_tree"];      //true or false, whether or not to write and save output TTree of selected events (and sb events). Adds a lot of runtime.
  if (tree_cfg) { write_tree = tree_cfg.as<bool>(); }
    
  // You're required to make some decisions
  PlotUtils::MinervaUniverse::SetNuEConstraint(true); //the neutrino-electron flux scattering constraint?? So should leave true?
  PlotUtils::MinervaUniverse::SetPlaylist(options.m_plist_string); 
  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(12); //Apparently this is: Analysis neutrino identity -- used for flux weights
  PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
  PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false); //to do: investigate
  PlotUtils::MinervaUniverse::RPAMaterials(false); //Sets whether you want to use different nuclei for RPA reweight, with false it just uses carbon which is fine for me

  //Now that we've defined what a cross section is, decide which sample and model
  //we're extracting a cross section for.

  //cuts defined here
  PlotUtils::Cutter<CVUniverse, CCNuEEvent>::reco_t sideband_cuts, preCuts;
  PlotUtils::Cutter<CVUniverse, CCNuEEvent>::truth_t signalDefinition, phaseSpace;

  //Loop through cut names in the config file, find corresponding entry in cut registry (defined in cuts/NuETKICuts.h),
  // construct it using parameters in the config, and finally add it to the list of cuts to be evaluated
  for (auto it = precuts_cfg.begin(); it != precuts_cfg.end(); ++it) {
    std::string cutName = it->first.as<std::string>();
    YAML::Node cutNode = it->second;
    
    if (!cutNode["enabled"] || !cutNode["enabled"].as<bool>())
      continue;
    
    // find matching registry entry
    auto found = std::find_if( ALL_CUTS.begin(), ALL_CUTS.end(), [&](const CutDef& def){ return def.name == cutName; });
    
    if (found == ALL_CUTS.end()) {
      std::cout << "Unknown cut in config: " << cutName << std::endl;
      continue;
    }
    preCuts.emplace_back(found->creator(cutNode));
  }
  
  //Sideband cuts, same exact thing as the precuts but add the cut to sideband_cuts instead of preCuts
  //also save the names so I can have a branch next to the cut results with the corresponding names to help keep track
  std::vector<std::string> sideband_cut_names; 
  for (auto it = sbcuts_cfg.begin(); it != sbcuts_cfg.end(); ++it) {
    std::string cutName = it->first.as<std::string>();
    YAML::Node cutNode = it->second;
    
    if (!cutNode["enabled"] || !cutNode["enabled"].as<bool>())
      continue;
    
    // find matching registry entry
    auto found = std::find_if( ALL_CUTS.begin(), ALL_CUTS.end(), [&](const CutDef& def){ return def.name == cutName; });
    
    if (found == ALL_CUTS.end()) {
      std::cout << "Unknown cut in config: " << cutName << std::endl;
      continue;
    }
    sideband_cuts.emplace_back(found->creator(cutNode));
    sideband_cut_names.emplace_back(cutName);
  }
  n_sideband_cuts = sideband_cut_names.size();
  
  const double minZ = precuts_cfg["ZRange"]["minZ"].as<double>();
  const double maxZ = precuts_cfg["ZRange"]["maxZ"].as<double>();
  const double apothem = precuts_cfg["Apothem"]["apothem"].as<double>(); //these are values in mm pulled from config file
  //Hang also has a "hits nucleus" requirement in his signal definition, whats the best way to implement that...?
  //decided not to use configs for my signal definition since it basically never changes
  signalDefinition.emplace_back(new truth::IsNue<CVUniverse>());
  signalDefinition.emplace_back(new truth::IsCC<CVUniverse>());
  signalDefinition.emplace_back(new truth::TrueElectronEnergy<CVUniverse>(precuts_cfg["ElectronEnergy"]["minElectronEnergy"].as<double>()));
  signalDefinition.emplace_back(new truth::HasSignalProton<CVUniverse>());
  signalDefinition.emplace_back(new truth::HasNoMeson<CVUniverse>());
  signalDefinition.emplace_back(new truth::HasNoPhoton<CVUniverse>());

  phaseSpace.emplace_back(new truth::ZRange<CVUniverse>("Tracker", minZ, maxZ));
  phaseSpace.emplace_back(new truth::Apothem<CVUniverse>(apothem));

  PlotUtils::Cutter<CVUniverse, CCNuEEvent> mycuts(std::move(preCuts), std::move(sideband_cuts) , std::move(signalDefinition),std::move(phaseSpace));

  //Define CV model using config inputs
  std::vector<std::unique_ptr<PlotUtils::Reweighter<CVUniverse, CCNuEEvent>>> MyModel;  
  //Check if I ask for a specific MnvTune defined in util/ReweighterRegistry.h, if so add those reweighters (defined in util/ReweighterRegistry.h)
  if (model_cfg.IsScalar()) {
    std::string tuneName = model_cfg.as<std::string>();
    const std::vector<ReweighterDef>* tuneReweights = nullptr;
    if      (tuneName == "MnvTunev1") tuneReweights = &MnvTunev1_reweights;
    else if (tuneName == "MnvTunev2") tuneReweights = &MnvTunev2_reweights;
    else {
      std::cout << "Unknown tune in config: " << tuneName << std::endl;
      // handle error — throw, return, etc.
    }
    if (tuneReweights) {
      for (const auto& def : *tuneReweights) {
	MyModel.emplace_back(def.creator(def.defaultNode));
      }
    }
  } else if (model_cfg.IsMap()) { //if not search through all known reweighters and add the ones in config individually
    for (auto it = model_cfg.begin(); it != model_cfg.end(); ++it) {
      std::string reweighterName = it->first.as<std::string>();
      YAML::Node reweighterNode = it->second;      
      if (!reweighterNode["enabled"] || !reweighterNode["enabled"].as<bool>())
	continue;      
      // find matching registry entry
      auto found = std::find_if( ALL_REWEIGHTS.begin(), ALL_REWEIGHTS.end(), [&](const ReweighterDef& def){ return def.name == reweighterName; });      
      if (found == ALL_REWEIGHTS.end()) {
	std::cout << "Unknown reweighter in config: " << reweighterName << std::endl;
	continue;
      }
      MyModel.emplace_back(found->creator(reweighterNode));
    }
  }
  PlotUtils::Model<CVUniverse, CCNuEEvent> model(std::move(MyModel));

  // Make a map of systematic universes
  // Leave out systematics when making validation histograms
  const bool doSystematics = syst_cfg.as<bool>(); //whether or not to do systematics is controlled by yaml config file
                                                  //to do: control which systematics get included? might not really be necessary/useful though
  if(!doSystematics){
    std::cout << "Skipping systematics (except 1 flux universe)\n";
    PlotUtils::MinervaUniverse::SetNFluxUniverses(2); //Necessary to get Flux integral later...  Doesn't work with just 1 flux universe though because _that_ triggers "spread errors".
  }

  std::map< std::string, std::vector<CVUniverse*> > error_bands;
  //defined in NuE_TKI/systematics/Systematics.h
  if(doSystematics) error_bands = GetNuETKISystematics(options.m_mc);
  //if(doSystematics) error_bands = GetTestSystematics(options.m_mc);
  else{
    std::map<std::string, std::vector<CVUniverse*> > band_flux = PlotUtils::GetFluxSystematicsMap<CVUniverse>(options.m_mc, CVUniverse::GetNFluxUniverses());
    error_bands.insert(band_flux.begin(), band_flux.end()); //Necessary to get flux integral later...
  }
  error_bands["cv"] = {new CVUniverse(options.m_mc)};
  std::map< std::string, std::vector<CVUniverse*> > truth_bands;
  if(doSystematics) { truth_bands = GetNuETKISystematics(options.m_truth); }
  //ExtractCrossSection throws an error when you run without systematics
  //Because it can't divide effNum/effDenom, as effDenom doesn't have a flux error band. So just adding this here. 
  if(!doSystematics) {
    std::map<std::string, std::vector<CVUniverse*> > band_flux = PlotUtils::GetFluxSystematicsMap<CVUniverse>(options.m_truth, CVUniverse::GetNFluxUniverses());
    truth_bands.insert(band_flux.begin(), band_flux.end());
  }
  truth_bands["cv"] = {new CVUniverse(options.m_truth)};

  //Variables to plot as histograms are created here
  //Loops through the variable registry in util/VariableRegistry.h, yaml config file determines which ones are filled and saved
  //binning defined in util/Binning.h
  std::vector<Variable*> vars;
  for (const auto& def : ALL_VARIABLES) {
    if (variables_cfg &&
        variables_cfg[def.name] &&
        variables_cfg[def.name].as<bool>()) {

      vars.push_back(new Variable(
        def.name,
        def.title,
        *def.bins,
        def.reco,
        def.truth
      ));
    }
  }
  
  //I don't have any 2D variables, so no config loop needed
  std::vector<Variable2D*> vars2D;

  CVUniverse* data_universe = new CVUniverse(options.m_data);
  std::vector<CVUniverse*> data_band = {data_universe};
  std::map<std::string, std::vector<CVUniverse*> > data_error_bands;
  data_error_bands["cv"] = data_band;
  
  //The actual sideband histograms are created here, they are "Studies".
  //Useful to have these as configs because sometimes for testing I don't need to fill sidebands
  std::vector<Study*> data_studies;
  std::vector<Study*> studies;
  if (sidebands_cfg["MeanFrontDEDXSB"] && sidebands_cfg["MeanFrontDEDXSB"].as<bool>()) {
    data_studies.push_back(new MeanFrontDEDXSideband(vars, error_bands, truth_bands, data_band));
    studies.push_back(new MeanFrontDEDXSideband(vars, error_bands, truth_bands, data_band));
  }
  if (sidebands_cfg["MichelSB"] && sidebands_cfg["MichelSB"].as<bool>()) {
    data_studies.push_back(new MichelSideband(vars, error_bands, truth_bands, data_band));
    studies.push_back(new MichelSideband(vars, error_bands, truth_bands, data_band));
  }
  if (sidebands_cfg["NIsoClusSB"] && sidebands_cfg["NIsoClusSB"].as<bool>()) {
    data_studies.push_back(new NIsoClusSideband(vars, error_bands, truth_bands, data_band));
    studies.push_back(new NIsoClusSideband(vars, error_bands, truth_bands, data_band));
  }
  if (sidebands_cfg["NPiSB"] && sidebands_cfg["NPiSB"].as<bool>()) {
    data_studies.push_back(new NPiSideband(vars, error_bands, truth_bands, data_band));
    studies.push_back(new NPiSideband(vars, error_bands, truth_bands, data_band));
  }
  
  for(auto& var: vars) var->InitializeMCHists(error_bands, truth_bands);
  for(auto& var: vars) var->InitializeDATAHists(data_band);
  
  for(auto& var: vars2D) var->InitializeMCHists(error_bands, truth_bands);
  for(auto& var: vars2D) var->InitializeDATAHists(data_band);
 
  //Initialize my output tree stuff. if write_tree is false I still get an instance of the class, but it doesn't really do much or add any complexity/time
  g_OutputTreeManager.Init(mc_tuple_out_filename, write_tree, *options.m_mc->GetChain()->GetTree(), sideband_cut_names);

  // Loop entries and fill
  try
  {
    
    CVUniverse::SetTruth(false);
    t_mc_start = getTimeStamp();
    LoopAndFillEventSelection(options.m_mc, error_bands, vars, vars2D, studies, mycuts, model);
    CVUniverse::SetTruth(true);
    t_eff_start = getTimeStamp();
    LoopAndFillEffDenom(options.m_truth, truth_bands, vars, vars2D, mycuts, model);

    std::cout << "---------------- YAML config ----------------" << std::endl;
    std::cout << config << std::endl;
    //std::cout << model_cfg << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
    options.PrintMacroConfiguration(argv[0]);
    std::cout << "MC cut summary:\n" << mycuts << "\n";
    mycuts.resetStats();

    CVUniverse::SetTruth(false);
    t_data_start = getTimeStamp();    
    LoopAndFillData(options.m_data, data_band, vars, vars2D, data_studies, mycuts);
    std::cout << "Data cut summary:\n" << mycuts << "\n";
    t_data_end = getTimeStamp();

    
    //Write MC results
    TFile* mcOutDir = TFile::Open(mc_out_filename.c_str(), "RECREATE");

    if(!mcOutDir)
    {
      std::cerr << "Failed to open a file named " << mc_out_filename << " in the current directory for writing histograms.\n";
      return badOutputFile;
    }

    for(auto& study: studies) study->SaveOrDrawMC(*mcOutDir);
    //for(auto& study: studies) study->SaveOrDraw(*mcOutDir);
    for(auto& var: vars) var->WriteMC(*mcOutDir);
    for(auto& var: vars2D) var->Write(*mcOutDir);

    
    //Protons On Target
    auto mcPOT = new TParameter<double>("POTUsed", options.m_mc_pot);
    mcPOT->Write();

    PlotUtils::TargetUtils targetInfo;
    assert(error_bands["cv"].size() == 1 && "List of error bands must contain a universe named \"cv\" for the flux integral.");

    for(const auto& var: vars)
    {
      //Flux integral only if systematics are being done (temporary solution)
      
      auto flux_hist = util::GetFluxIntegral(*error_bands["cv"].front(), var->efficiencyNumerator->hist);
      flux_hist->Write((var->GetName() + "_reweightedflux_integrated").c_str());
      //Always use MC number of nucleons for cross section
      auto nNucleons = new TParameter<double>((var->GetName() + "_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons(minZ, maxZ, true, apothem));
      nNucleons->Write();
    }
    
    //Put this here after we finish writing out all the info to the MC histo file
    g_OutputTreeManager.WriteAndClose();
    
    //Write data results
    
    TFile* dataOutDir = TFile::Open(data_out_filename.c_str(), "RECREATE");
    if(!dataOutDir)
    {
      std::cerr << "Failed to open a file named " << data_out_filename << " in the current directory for writing histograms.\n";
      return badOutputFile;
    }
    
    for(auto& study: data_studies) study->SaveOrDrawData(*dataOutDir);
    for(auto& var: vars) var->WriteData(*dataOutDir);
    

    //Protons On Target
    auto dataPOT = new TParameter<double>("POTUsed", options.m_data_pot);
    dataPOT->Write();
    
    std::cout << "Success" << std::endl;
  }
  catch(const ROOT::exception& e)
  {
    std::cerr << "Ending on a ROOT error message.  No histograms will be produced.\n"
              << "If the message talks about \"TNetXNGFile\", this could be a problem with dCache.  The message is:\n"
              << e.what() << "\n" << USAGE << "\n";
    return badFileRead;
  }
  t_end = getTimeStamp();
  
  std::cout << "--------- grid job time benchmarking ---------\n";
  printElapsed("   time to set up & open input", t_start, t_mc_start);
  printElapsed("              time for MC loop", t_mc_start, t_eff_start);
  printElapsed("       time for eff denom loop", t_eff_start, t_data_start);
  printElapsed("            time for data loop", t_data_start, t_data_end);
  printElapsed("time to write and close output", t_data_end, t_end);
  
  return success;
}
