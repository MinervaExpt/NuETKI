#define MC_OUT_FILE_NAME "MC_June_12_2025_NuE_Only.root"
#define DATA_OUT_FILE_NAME "Data_June_12_2025_NuE_Only.root"
#define MC_TUPLE_OUT_FILE_NAME "TUPLE_MC_June_12_2025_NuE_Only.root"
//#define MC_OUT_FILE_NAME "efficiencyTesting.root"
//#define DATA_OUT_FILE_NAME "dummyData.root"

bool write_tree = true;
//bool write_tree = false;

#define USAGE					\
"\n*** USAGE ***\n"\
"runEventLoop <dataPlaylist.txt> <mcPlaylist.txt>\n\n"\
"*** Explanation ***\n"\
"Reduce MasterAnaDev AnaTuples to event selection histograms to extract a\n"\
"single-differential inclusive cross section for the 2021 MINERvA 101 tutorial.\n\n"\
"*** The Input Files ***\n"\
"Playlist files are plaintext files with 1 file name per line.  Filenames may be\n"\
"xrootd URLs or refer to the local filesystem.  The first playlist file's\n"\
"entries will be treated like data, and the second playlist's entries must\n"\
"have the \"Truth\" tree to use for calculating the efficiency denominator.\n\n"\
"*** Output ***\n"\
"Produces a two files, " MC_OUT_FILE_NAME " and " DATA_OUT_FILE_NAME ", with\n"\
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
#include "event/MichelEvent.h"
#include "systematics/Systematics.h"
#include "cuts/MaxPzMu.h"
#include "util/Variable.h"
#include "util/Variable2D.h"
#include "util/GetFluxIntegral.h"
#include "util/GetPlaylist.h"
#include "cuts/SignalDefinition.h"
#include "cuts/q3RecoCut.h"
#include "cuts/NuETKICuts.h"
#include "cuts/NuETKISignal.h"
//#include "Binning.h" //TODO: Fix me

//Studies and sidebands
#include "studies/Study.h"
//#include "studies/EMShowerScoreSB.h"
#include "studies/MichelSideband.h"
//#include "studies/E_lep_SB.h"
//#include "studies/ModEAvailSB.h"
//#include "studies/MeanFrontDEDX_SB.h"
//#include "studies/ESC_SB.h"

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
#include "PlotUtils/TargetUtils.h"
#pragma GCC diagnostic pop

//ROOT includes
#include "TParameter.h"

//c++ includes
#include <iostream>
#include <cstdlib> //getenv()
#include <ctime>

OutputTreeManager g_OutputTreeManager; //defining my output tree manager (struct declaration and full def in MichelEvent.h

//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFillEventSelection(
    PlotUtils::ChainWrapper* chain,
    std::map<std::string, std::vector<CVUniverse*> > error_bands,
    std::vector<Variable*> vars,
    std::vector<Variable2D*> vars2D,
    std::vector<Study*> studies,
    PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts,
    PlotUtils::Model<CVUniverse, MichelEvent>& model)
{
  assert(!error_bands["cv"].empty() && "\"cv\" error band is empty!  Can't set Model weight.");
  auto& cvUniv = error_bands["cv"].front();

  std::cout << "Starting MC reco loop...\n";
  
  const int nEntries = chain->GetEntries();
  
  for (Long64_t i=0; i<nEntries; ++i)
  {
    if(i%1000==0) std::cout << i << " / " << nEntries << "\r" << std::flush;
    
    MichelEvent cvEvent;
    cvUniv->SetEntry(i);
    model.SetEntry(*cvUniv, cvEvent);
    const double cvWeight = model.GetWeight(*cvUniv, cvEvent);

    //=========================================
    // Systematics loop(s)
    //=========================================
    for (auto band : error_bands)
    {
      std::vector<CVUniverse*> error_band_universes = band.second;
      for (auto universe : error_band_universes)
      {
        MichelEvent myevent; // make sure your event is inside the error band loop. 
    
        // Tell the Event which entry in the TChain it's looking at
        universe->SetEntry(i);
        // This is where you would Access/create a Michel
        //weight is ignored in isMCSelected() for all but the CV Universe.
	//if ths event doesnt pass selection cuts, skip

	//isMCSelected returns a bitset which sideband cuts are passed. It is all 1s other than specifically the sideband cuts, ie we lose info about specific precuts (which is fine). 
	//So if filling sidebands I don't want to check !results.all() because the sideband result will always be 1
	//Instead I'll make a bitmask (basically a comparison bitset which I'll use to check all bits except the first one, which can be whatever
	//and then I'll use that first cut result to fill my sideband (if 0, goes into sideband, if 1, goes into selected events)
	std::bitset<64> cut_results = michelcuts.isMCSelected(*universe, myevent, cvWeight);  
	constexpr int n_sidebands = 1; //just setting this manually, it's mostly for testing anyways, it won't change frequently or anything
	std::bitset<64> mask;
	mask.set(); //set all bits in mask to 1
	//mask.reset(0); //set first bit to 0 (first bit representing the first sideband cut, I can do more of these when I add more sidebands)
	mask &= ~std::bitset<64>((1ULL << n_sidebands) - 1); //very fancy bit manipulation from chat gpt, just sets the first (or more accurately, least significant) n bits to 0

	//if (!cut_results.all()) continue; //skips events that don't pass ALL CUTS
	if ((cut_results & mask) != mask) continue; //skips events if it doesn't pass all cuts OTHER than sideband cuts
	std::bitset<n_sidebands> sideband_results(cut_results.to_ullong() & ((1ULL << n_sidebands) - 1)); //also just for testing mostly, this just isolates the results of those n cuts and just drops the massive set of 1s that always precedes it, mostly so I can print it / write it with more readability
	if (write_tree){
	  universe->GetTree()->GetTree()->GetEntry(i); //how is this different than universe->SetEntry(i)? not sure but it definitely is...
	  myevent.entryNumber = i;
	}
	//g_OutputTreeManager.sidebandCutResults = sideband_results;
	
	const double weight = model.GetWeight(*universe, myevent); //Only calculate the per-universe weight for events that will actually use it.	
	const bool isSignal = michelcuts.isSignal(*universe, weight);

	int bkgd_ID;
	if (isSignal){ g_OutputTreeManager.selectionCategory = -999;	}
	else {
	  bkgd_ID = -1; //this is "other" selectionCategory
	  //NueCC with final state mesons (mostly pions)
	  if (abs(universe->GetTruthNuPDG())==12 && universe->GetCurrent()==1 && universe->GetHasFSMeson()) { bkgd_ID = 0; }
	  //Other NueCC. So this will be NON QELike nuecc, will include all the stuff with mesons
	  else if (abs(universe->GetTruthNuPDG())==12 && universe->GetCurrent()==1) { bkgd_ID = 1; }
	  //Other NC pi0, all "other" NC i get should have pi0 so if this is small and the "other" selectionCategory is large that's a problem
	  else if (universe->GetCurrent()==2 && universe->GetHasFSPi0()==1) { bkgd_ID = 2; }
	  //same logic here with cc numu pi0
	  else if (abs(universe->GetTruthNuPDG())==14 && universe->GetCurrent()==1 && universe->GetHasFSPi0()==1) { bkgd_ID = 3; }
	  
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
	
	if ((cut_results & ~mask) != ~mask) continue; //Check the remaining sideband cuts, everything that passes this has now passed ALL cuts

	//print out events if I want to
	/*
	if (isSignal){
	  std::cout << "Selected event: " << universe->GetInt("mc_run") << " | " << universe->GetInt("mc_subrun") << " | " << universe->GetInt("mc_nthEvtInFile")+1 << " , SIGNAL EVENT" << std::endl;
	} else {
	  std::cout << "Selected event: " << universe->GetInt("mc_run") << " | " << universe->GetInt("mc_subrun") << " | " << universe->GetInt("mc_nthEvtInFile")+1 << " , bkgd_ID = " << bkgd_ID << std::endl;
	}
	*/
	
	g_OutputTreeManager.Fill("Signal_Region", myevent.entryNumber); 

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
          }

          for(auto& var: vars2D)
          {
            var->efficiencyNumerator->FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
          }
        }
	else{ //The actual logic for which category an event gets put into happens above (so that we can save it in the tuple & sidebands)
          for(auto& var: vars) (*var->m_backgroundHists)[bkgd_ID].FillUniverse(universe, var->GetRecoValue(*universe), weight);
	  //for testing my proton efficiency in true vars 2d plot
          for(auto& var: vars2D) (*var->m_backgroundHists)[bkgd_ID].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
          //for(auto& var: vars2D) (*var->m_backgroundHists)[bkgd_ID].FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
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
				PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts)

{
  std::cout << "Starting data loop...\n";
  const int nEntries = data->GetEntries();
  for (int i=0; i<data->GetEntries(); ++i) {
    for (auto universe : data_band) {
      universe->SetEntry(i);
      if(i%1000==0) std::cout << i << " / " << nEntries << "\r" << std::flush;
      MichelEvent myevent;
      
      //if (!michelcuts.isDataSelected(*universe, myevent).all()) continue;
      std::bitset<64> cut_results = michelcuts.isDataSelected(*universe, myevent);  
      constexpr int n_sidebands = 1; //just setting this manually, it's mostly for testing anyways, it won't change frequently or anything
      std::bitset<64> mask;
      mask.set(); //set all bits in mask to 1
      mask &= ~std::bitset<64>((1ULL << n_sidebands) - 1); //very fancy bit manipulation from chat gpt, just sets the first (or more accurately, least significant) n bits to 0
      if ((cut_results & mask) != mask) continue; //skips events if it doesn't pass all cuts OTHER than sideband cuts
      std::bitset<n_sidebands> sideband_results(cut_results.to_ullong() & ((1ULL << n_sidebands) - 1)); //also just for testing mostly, this just isolates the results of those n cuts and just drops the massive set of 1s that always precedes it, mostly so I can print it / write it with more readability

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
    				PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts,
                                PlotUtils::Model<CVUniverse, MichelEvent>& model)
{
  assert(!truth_bands["cv"].empty() && "\"cv\" error band is empty!  Could not set Model entry.");
  auto& cvUniv = truth_bands["cv"].front();

  std::cout << "Starting efficiency denominator loop...\n";
  const int nEntries = truth->GetEntries();
  for (int i=0; i<nEntries; ++i)
  {
    if(i%1000==0) std::cout << i << " / " << nEntries << "\r" << std::flush;

    MichelEvent cvEvent;
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
        MichelEvent myevent; //Only used to keep the Model happy

        // Tell the Event which entry in the TChain it's looking at
        universe->SetEntry(i);

        if (!michelcuts.isEfficiencyDenom(*universe, cvWeight)) continue; //Weight is ignored for isEfficiencyDenom() in all but the CV universe 
        const double weight = model.GetWeight(*universe, myevent); //Only calculate the weight for events that will use it

        //Fill efficiency denominator now: 
        for(auto var: vars)
        {
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

  //Validate input.
  //I expect a data playlist file name and an MC playlist file name which is exactly 2 arguments.
  const int nArgsExpected = 2;
  if(argc != nArgsExpected + 1) //argc is the size of argv.  I check for number of arguments + 1 because
                                //argv[0] is always the path to the executable.
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

  const bool doCCQENuValidation = (reco_tree_name == "CCQENu"); //Enables extra histograms and might influence which systematics I use if the reco tree name is CCQENu...

  //const bool is_grid = false; //TODO: Are we going to put this back?  Gonzalo needs it iirc.
  PlotUtils::MacroUtil options(reco_tree_name, mc_file_list, data_file_list, "minervame1A", true);
  options.m_plist_string = util::GetPlaylist(*options.m_mc, true); //TODO: Put GetPlaylist into PlotUtils::MacroUtil

  // You're required to make some decisions
  PlotUtils::MinervaUniverse::SetNuEConstraint(true); //what is this? I think it's the neutrino-electron flux scattering constraint?? So should leave true?
  PlotUtils::MinervaUniverse::SetPlaylist(options.m_plist_string); //TODO: Infer this from the files somehow?
  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(12); //Apparently this is: Analysis neutrino identity -- used for flux weights
  PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
  PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false);

  //Now that we've defined what a cross section is, decide which sample and model
  //we're extracting a cross section for.

  //CUTS DEFINED HERE! reco_t and truth_t are just vectors of Cut objects (from MAT), emplace_back just adds something to the end of the vector
  //So I guess reco:: and truth:: are just types of Cuts? -> they are namespaces in CCInclusiveCuts.h and CCInclusiveSignal.h respectively...
  PlotUtils::Cutter<CVUniverse, MichelEvent>::reco_t sidebands, preCuts;
  PlotUtils::Cutter<CVUniverse, MichelEvent>::truth_t signalDefinition, phaseSpace;

  //So I guess the question is, what are these reco:: and truth:: namespaces? because cuts are defined as members of those...
  const double minZ = 5980, maxZ = 8422, apothem = 850; //All in mm
  const double minFirstFireFraction = 0.25; //For Afterpulsing quality cut
  const double maxDSCalRatio = 0.2; //vis E in DS HCAL / ECAL
  const double maxSideCalRatio = 0.05; //vis E in Side HCAL / ECAL (kinda useless tbh)
  const int minVertexTrackMultiplicity = 0; //Amount of tracks coming directly out of the vertex
  const int maxVertexTrackMultiplicity = 2; //cut is inclusive (so will include events with 0 and 2 if that's what its set to)
  const double minTransverseGapScore = 15;
  const double minNonMIPClusFrac = 0.7; 

  const double minEMTrackScore = 0.7;
  const double maxMeanFrontDEDX = 2.4; //in MeV/cm
  const double minElectronEnergy = 2.5; //in GeV
  const double maxModifiedEAvail = 1; //in GeV
  const double maxESCChi2 = 20;

  const int maxNIsoBlobs = 2; //cut is exclusive
  const double maxIsoBlobEnergy = 250; //in MeV
  const double furthestUpstreamIsoBlobStartZ = -50; //in mm, aka do not allow any blobs where blob start Z - vtx Z is lower than x
  //ie do not allow any events with blobs that start more than -x mm upstream of vtx
  
  //"precuts"
  preCuts.emplace_back(new reco::HasTracks<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::NoVertexMismatch<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::ZRange<CVUniverse, MichelEvent>("Vertex in Tracker", minZ, maxZ));
  preCuts.emplace_back(new reco::Apothem<CVUniverse, MichelEvent>(apothem));
  preCuts.emplace_back(new reco::StartPointVertexMultiplicity<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::Afterpulsing<CVUniverse, MichelEvent>(minFirstFireFraction));
  preCuts.emplace_back(new reco::NoDeadtime<CVUniverse, MichelEvent>(1, "Deadtime"));
  preCuts.emplace_back(new reco::HasNoBackExitingTracks<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::DSCalVisE<CVUniverse, MichelEvent>(maxDSCalRatio));
  preCuts.emplace_back(new reco::ODCalVisE<CVUniverse, MichelEvent>(maxSideCalRatio));
  preCuts.emplace_back(new reco::VertexTrackMultiplicity<CVUniverse, MichelEvent>(minVertexTrackMultiplicity, maxVertexTrackMultiplicity));
  preCuts.emplace_back(new reco::TransverseGapScore<CVUniverse, MichelEvent>(minTransverseGapScore));
  preCuts.emplace_back(new reco::NonMIPClusterFraction<CVUniverse, MichelEvent>(minNonMIPClusFrac));
  
  //not technically sidebands, but maybe more sophisticated/specific cuts
  preCuts.emplace_back(new reco::EMLikeTrackScore<CVUniverse, MichelEvent>(minEMTrackScore));
  //preCuts.emplace_back(new reco::MichelCut<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::MeanFrontdEdX<CVUniverse, MichelEvent>(maxMeanFrontDEDX));
  preCuts.emplace_back(new reco::ElectronEnergy<CVUniverse, MichelEvent>(minElectronEnergy));
  preCuts.emplace_back(new reco::ModifiedEavailable<CVUniverse, MichelEvent>(maxModifiedEAvail));

  //This is reduntant, but separating them helps to see the effect of requiring a proton track vs the actual chi^2 / esc node cut effects
  preCuts.emplace_back(new reco::ProtonInEvent<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::NIsoBlobs<CVUniverse, MichelEvent>(maxNIsoBlobs));
  preCuts.emplace_back(new reco::IsoBlobEnergy<CVUniverse, MichelEvent>(maxIsoBlobEnergy));
  preCuts.emplace_back(new reco::UpstreamIsoBlob<CVUniverse, MichelEvent>(furthestUpstreamIsoBlobStartZ));
  preCuts.emplace_back(new reco::ESC<CVUniverse, MichelEvent>(maxESCChi2));
  
  //older / obsolete cuts
  //preCuts.emplace_back(new reco::Eavailable<CVUniverse, MichelEvent>());
  //preCuts.emplace_back(new reco::ProtonMomentum<CVUniverse, MichelEvent>());
  //preCuts.emplace_back(new reco::Psi<CVUniverse, MichelEvent>());
  //preCuts.emplace_back(new reco::ProtonTheta<CVUniverse, MichelEvent>());
  //preCuts.emplace_back(new reco::ElectronPt<CVUniverse, MichelEvent>());
  //preCuts.emplace_back(new reco::EleptonSin2Theta<CVUniverse, MichelEvent>());
  //preCuts.emplace_back(new reco::ProtonInEvent<CVUniverse, MichelEvent>());  

  //Hangs cuts exactly, when I want to run that selection
  /*
  preCuts.emplace_back(new reco::HasTracks<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::HasNoBackExitingTracks<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::EMLikeTrackScore<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::DSCalVisE<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::ODCalVisE<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::VertexTrackMultiplicity<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::Afterpulsing<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::NoDeadtime<CVUniverse, MichelEvent>(1, "Deadtime"));
  preCuts.emplace_back(new reco::StartPointVertexMultiplicity<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::MeanFrontdEdX<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::NonMIPClusterFraction<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::TransverseGapScore<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::ZRange<CVUniverse, MichelEvent>("Tracker", minZ, maxZ));
  preCuts.emplace_back(new reco::Apothem<CVUniverse, MichelEvent>(apothem));
  preCuts.emplace_back(new reco::Eavailable<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new reco::ElectronEnergy<CVUniverse, MichelEvent>());
  */

  //Sidebands
  //sidebands.emplace_back(new reco::MeanFrontdEdXAbove<CVUniverse, MichelEvent>());
  //sidebands.emplace_back(new reco::EMLikeTrackScore<CVUniverse, MichelEvent>());
  sidebands.emplace_back(new reco::MichelCut<CVUniverse, MichelEvent>());
  //sidebands.emplace_back(new reco::ElectronEnergy<CVUniverse, MichelEvent>());
  //sidebands.emplace_back(new reco::ModifiedEavailable<CVUniverse, MichelEvent>());
  //sidebands.emplace_back(new reco::MeanFrontdEdX<CVUniverse, MichelEvent>());
  //sidebands.emplace_back(new reco::ESC<CVUniverse, MichelEvent>());
  
  //Hang also has a "hits nucleus" requirement in his signal definition, whats the best way to implement that...?
  signalDefinition.emplace_back(new truth::IsNue<CVUniverse>());
  signalDefinition.emplace_back(new truth::IsCC<CVUniverse>());
  signalDefinition.emplace_back(new truth::TrueElectronEnergy<CVUniverse>(minElectronEnergy));
  
  signalDefinition.emplace_back(new truth::HasSignalProton<CVUniverse>());
  signalDefinition.emplace_back(new truth::HasNoMeson<CVUniverse>());
  signalDefinition.emplace_back(new truth::HasNoPhoton<CVUniverse>());

  phaseSpace.emplace_back(new truth::ZRange<CVUniverse>("Tracker", minZ, maxZ));
  phaseSpace.emplace_back(new truth::Apothem<CVUniverse>(apothem));

  PlotUtils::Cutter<CVUniverse, MichelEvent> mycuts(std::move(preCuts), std::move(sidebands) , std::move(signalDefinition),std::move(phaseSpace));

  std::vector<std::unique_ptr<PlotUtils::Reweighter<CVUniverse, MichelEvent>>> MnvTunev1;

  //these cause a problem for some reason?? or one of them does... RPA I think actually
  
  //MnvTunev1.emplace_back(new PlotUtils::FluxAndCVReweighter<CVUniverse, MichelEvent>());
  //MnvTunev1.emplace_back(new PlotUtils::GENIEReweighter<CVUniverse, MichelEvent>(true, false));
  //MnvTunev1.emplace_back(new PlotUtils::LowRecoil2p2hReweighter<CVUniverse, MichelEvent>());
  //MnvTunev1.emplace_back(new PlotUtils::MINOSEfficiencyReweighter<CVUniverse, MichelEvent>());

  //MnvTunev1.emplace_back(new PlotUtils::RPAReweighter<CVUniverse, MichelEvent>());
  
  
  PlotUtils::Model<CVUniverse, MichelEvent> model(std::move(MnvTunev1));

  // Make a map of systematic universes
  // Leave out systematics when making validation histograms
  //const bool doSystematics = (getenv("MNV101_SKIP_SYST") == nullptr); 
  const bool doSystematics = false; //hard coding to false for now, we'll get to those eventually (carlos p, nov 2023)
  if(!doSystematics){
    std::cout << "Skipping systematics (except 1 flux universe) because environment variable MNV101_SKIP_SYST is set.\n";
    PlotUtils::MinervaUniverse::SetNFluxUniverses(2); //Necessary to get Flux integral later...  Doesn't work with just 1 flux universe though because _that_ triggers "spread errors".
  }

  std::map< std::string, std::vector<CVUniverse*> > error_bands;
  if(doSystematics) error_bands = GetStandardSystematics(options.m_mc);
  else{
    std::map<std::string, std::vector<CVUniverse*> > band_flux = PlotUtils::GetFluxSystematicsMap<CVUniverse>(options.m_mc, CVUniverse::GetNFluxUniverses());
    error_bands.insert(band_flux.begin(), band_flux.end()); //Necessary to get flux integral later...
  }
  error_bands["cv"] = {new CVUniverse(options.m_mc)};
  std::map< std::string, std::vector<CVUniverse*> > truth_bands;
  if(doSystematics) truth_bands = GetStandardSystematics(options.m_truth);
  truth_bands["cv"] = {new CVUniverse(options.m_truth)};
  
  //Carlos: I should clean this up at some point
  std::vector<double> dansPTBins = {0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1, 1.25, 1.5, 2.5, 4.5},
                      dansPzBins = {1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20, 40, 60},
                      //protonMomentumBins = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0},
                      protonMomentumBins = {0.0,0.08,0.16,0.24,0.32,0.40,0.48,0.56,0.64,0.72,0.8,0.88,0.96,1.04,1.12,1.2,1.28,1.36,1.44,1.52,1.6},
                      //Pt_bins = {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2},
		      Pt_bins = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5},
		      KE_bins = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0, 1.25, 1.5},
		      EavailBins = {0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.8, 1, 1.2},
		      electronAngleBins = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40},
		      protonAngleBins = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180},
                      phiAngleBins = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100},
                      deltaPtXBins = {-2, -1.8, -1.6, -1.4, -1.2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2},
                      electronEnergyBins = {0.0,0.5,1.0,1.5,2,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5,13.0,13.5,14.0,14.5,20};
                      //electronEnergyBins = {0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,9,10,12.5,15,17.5,20};
  
  std::vector<Variable*> vars = {
    //Reco variables
    new Variable("E_nu", "E_{#nu} [GeV]", electronEnergyBins, &CVUniverse::GetEnu, &CVUniverse::GetEnuTrue),
    new Variable("E_avail", "E_{avail} [GeV]", EavailBins, &CVUniverse::GetEavail, &CVUniverse::GetEavailTrue),
    new Variable("Lepton_Pt", "p_{T,e} [GeV/c]", Pt_bins, &CVUniverse::GetElectronPt, &CVUniverse::GetElectronPtTrue),
    new Variable("E_lep", "E_{e} [GeV]", electronEnergyBins, &CVUniverse::GetElectronEnergy, &CVUniverse::GetElepTrueGeV),
    new Variable("Theta_lep", "#theta_{e} [deg]", electronAngleBins, &CVUniverse::GetElectronThetaDeg, &CVUniverse::GetEnuTrue),    
    new Variable("RecoProtonP", "P_{p} [GeV/c]", protonMomentumBins, &CVUniverse::GetProtonP, &CVUniverse::GetProtonPTrue),    
    new Variable("RecoProtonTheta", "#theta_{p}", protonAngleBins, &CVUniverse::GetProtonThetaDeg, &CVUniverse::GetProtonThetaTrue),    
    //tki vars start here
    new Variable("DeltaPt", "#deltaP_{T} [GeV/c]", Pt_bins, &CVUniverse::GetDeltaPt, &CVUniverse::GetDummyVar),
    new Variable("DeltaPtX", "#deltaP_{T,x} [GeV/c]", deltaPtXBins, &CVUniverse::GetDeltaPtX, &CVUniverse::GetDummyVar),//these 2 are wrong, for sure
    new Variable("DeltaPtY", "#deltaP_{T,y} [GeV/c]", deltaPtXBins, &CVUniverse::GetDeltaPtY, &CVUniverse::GetDummyVar),//
    new Variable("AlphaPt", "#delta#alpha_{T} [deg]", protonAngleBins, &CVUniverse::GetAlphaT, &CVUniverse::GetDummyVar),
    new Variable("PhiPt", "#delta#phi_{T} [deg]", phiAngleBins, &CVUniverse::GetPhiT, &CVUniverse::GetDummyVar),
    
    //true variables
    //new Variable("E_lep", "True Electron Energy", electronEnergyBins, &CVUniverse::GetElectronEnergyTrue, &CVUniverse::GetElectronEnergyTrue),
    //new Variable("ProtonKE", "T_{p} true [GeV]", KE_bins, &CVUniverse::GetDummyVar, &CVUniverse::GetProtonKETrue),
    //new Variable("ProtonP", "p_{p} true [GeV/c]", protonMomentumBins, &CVUniverse::GetDummyVar, &CVUniverse::GetProtonPTrue),
    //new Variable("Theta_p", "#theta_{p} true [deg]", protonAngleBins, &CVUniverse::GetDummyVar, &CVUniverse::GetProtonThetaTrue),
    //new Variable("Pt_p", "p_{T,p} true [GeV/c]", Pt_bins, &CVUniverse::GetDummyVar, &CVUniverse::GetProtonPtTrue),
    //new Variable("Theta_p_e", "#theta_{e,p} true [deg]", protonAngleBins, &CVUniverse::GetDummyVar, &CVUniverse::GetOpeningAngleTrue),
    //new Variable("Pt_lep", "p_{T,e} true [GeV]", Pt_bins, &CVUniverse::GetDummyVar, &CVUniverse::GetElectronPtTrue)
  };
  
/*
  //bins for cut plots
  std::vector<double> binaryCut_bins = {-0.5,0.5,1.5},
    VertexZ_bins,
    //VertexApothem_bins, //just gonna use binary bins, not an easy way to show this other it is or it isn't in the hexagon
    StartPointVertexMultiplicity_bins = {-1,0,1,2,3,4,5,6,7,8},
    Afterpulsing_bins = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1},
    Deadtime_bins = {0,1,2,3,4,5,6},
    DSCalVisE_bins = {0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25},
    ODCalVisE_bins = {0,0.0025,0.005,0.0075,0.01,0.0125,0.015,0.0175,0.02,0.0225,0.025,0.0275,0.03,0.0325,0.035,0.0375,0.04,0.0425,0.045,0.0475,0.05,0.0525,0.055,0.0575,0.06},
    VertexTrackMultiplicity_bins = {0,1,2,3,4,5,6},
    //TransverseGapScore_bins = {0,2,4,6,8,10,12,14,15,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50},
    TransverseGapScore_bins = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100},
    NonMIPClusFrac_bins = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1},
    //EMScore_bins = {0.499,0.549,0.599,0.649,0.699,0.749,0.799,0.849,0.899,0.949,0.999,1.049},
    EMScore_bins = {0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05},
    NMichel_bins = {0, 1, 2, 3, 4, 5},
    MeanFrontDEDX_bins = {0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.2,4.4,4.6,4.8,5,5.2,5.4,5.6,5.8,6,6.2,6.4,6.6,6.8,7,8,9,10},
    Modified_E_avail_bins = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2},
    ESCChi2_bins = {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40},
    Psi_bins = {0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5},
    Etheta_bins = {0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.015,0.02,0.025};
    //E_lep_Sin2Theta_bins = {0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.015,0.02,0.025,0.03}; 
    //Etheta_bins = {0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1};

  //5980, maxZ = 8422,
  const double VertexZBinWidth = 100; //mm
  for(int whichBin = 0; whichBin < 27 + 1; ++whichBin) VertexZ_bins.push_back(5880 + (VertexZBinWidth * whichBin));

  
  //std::vector<Variable*> vars = {new Variable("VertexTrackMultiplicity", "Vertex Track Multiplicity", VertexTrackMultiplicity_bins, &CVUniverse::GetVertexTrackMultiplicity, &CVUniverse::GetDummyVar)};
  
  //Cut variables (for when I want to plot the things I'm cutting on
  std::vector<Variable*> vars = {
    new Variable("NoVertexMismatch", "NoVertexMismatch", binaryCut_bins, &CVUniverse::GetHasNoVertexMismatch, &CVUniverse::GetDummyVar),
    new Variable("VertexZ", "VertexZ [mm]", VertexZ_bins, &CVUniverse::GetVertexZ, &CVUniverse::GetDummyVar),
    new Variable("InApothem", "InApothem", binaryCut_bins, &CVUniverse::GetWithinFiducialApothem, &CVUniverse::GetDummyVar),
    new Variable("StartPointVertexMultiplicity", "StartPointVertexMultiplicity", StartPointVertexMultiplicity_bins, &CVUniverse::GetStartPointVertexMultiplicity, &CVUniverse::GetDummyVar),
    new Variable("Afterpulsing", "Electron prong First Fire fraction", Afterpulsing_bins, &CVUniverse::GetFirstFireFraction, &CVUniverse::GetDummyVar),
    new Variable("Deadtime", "N dead discr pairs upstream", Deadtime_bins, &CVUniverse::GetTDead, &CVUniverse::GetDummyVar),
    new Variable("NoBackExitingTracks", "NoBackExitingTracks", binaryCut_bins, &CVUniverse::GetExitsBack, &CVUniverse::GetDummyVar),
    new Variable("DSCalVisE", "Ratio of DS HCalVisE/ECalVisE", DSCalVisE_bins, &CVUniverse::GetDSCalVisE, &CVUniverse::GetDummyVar),
    new Variable("ODCalVisE", "Ratio of OD HCalVisE/ECalVisE", ODCalVisE_bins, &CVUniverse::GetODCalVisE, &CVUniverse::GetDummyVar),
    new Variable("VertexTrackMultiplicity", "Vertex Track Multiplicity", VertexTrackMultiplicity_bins, &CVUniverse::GetVertexTrackMultiplicity, &CVUniverse::GetDummyVar),
    new Variable("TransverseGapScore", "Transverse Gap Score", TransverseGapScore_bins, &CVUniverse::GetTransverseGapScore, &CVUniverse::GetDummyVar),
    new Variable("NonMIPClusFrac", "Non MIP Cluster Fraction", NonMIPClusFrac_bins, &CVUniverse::GetNonMIPClusFrac, &CVUniverse::GetDummyVar),
    new Variable("EMScore", "EMShower Score", EMScore_bins, &CVUniverse::GetEMLikeShowerScore, &CVUniverse::GetDummyVar),
    new Variable("NMichels", "N michels", NMichel_bins, &CVUniverse::GetImprovedNMichel, &CVUniverse::GetDummyVar),
    new Variable("MeanFrontDEDX", "Mean Front dE/dX [MeV/cm]", MeanFrontDEDX_bins, &CVUniverse::GetMeanFrontdEdx, &CVUniverse::GetDummyVar),
    new Variable("E_lep", "E_lep [GeV]", electronEnergyBins, &CVUniverse::GetElectronEnergy, &CVUniverse::GetDummyVar),
    //new Variable("E_avail", "E_avail", Modified_E_avail_bins, &CVUniverse::GetEavail, &CVUniverse::GetDummyVar),
    new Variable("Modified_E_avail", "E_avail - sum(proton_E) [GeV]", Modified_E_avail_bins, &CVUniverse::GetModifiedEavail, &CVUniverse::GetDummyVar),
    new Variable("ESCChi2", "ESC Proton Node Chi2", ESCChi2_bins, &CVUniverse::GetProtonESCNodeChi2, &CVUniverse::GetDummyVar),

    //Potential cuts, but not ones I actually do. Can make all of these at once in the same job, unlike the others? well maybe not the proton ones, since ESC kinda messes that up
    new Variable("Psi", "Psi", Psi_bins, &CVUniverse::GetPsi, &CVUniverse::GetDummyVar),
    new Variable("ProtonMomentum", "Primary Proton momentum [GeV]", protonMomentumBins, &CVUniverse::GetProtonP, &CVUniverse::GetDummyVar),
    new Variable("ProtonTheta", "{theta}_proton [deg]", protonAngleBins, &CVUniverse::GetProtonThetaDeg, &CVUniverse::GetDummyVar),
    new Variable("LeptonPt", "LeptonPt [GeV]", Pt_bins, &CVUniverse::GetElectronPt, &CVUniverse::GetDummyVar),
    new Variable("Etheta", "E_lep * theta_lep^2", Etheta_bins, &CVUniverse::GetETheta, &CVUniverse::GetDummyVar)
    };*/
  
  std::vector<Variable2D*> vars2D;
  //2D proton efficiency plot, should be true proton angle in x and true proton momentum in y.
  //Jefferey has bins of size 5 degrees from 0 to 180 for angle
  //and bins of 0.04 GeV from 0 to 1.6 GeV
  //Let's just double these for now, because my statistics are significantly worse... 
  //vars2D.push_back(new Variable2D(*vars[7], *vars[6])); 

  CVUniverse* data_universe = new CVUniverse(options.m_data);
  std::vector<CVUniverse*> data_band = {data_universe};
  std::map<std::string, std::vector<CVUniverse*> > data_error_bands;
  data_error_bands["cv"] = data_band;
  
  std::vector<Study*> data_studies;
  data_studies.push_back(new MichelSideband(vars, error_bands, truth_bands, data_band));

  
  std::vector<Study*> studies;
  //studies.push_back(new EMShowerScoreSB(vars, error_bands, truth_bands, data_band));
  studies.push_back(new MichelSideband(vars, error_bands, truth_bands, data_band));
  //studies.push_back(new E_lep_SB(vars, error_bands, truth_bands, data_band));
  //studies.push_back(new ModEAvailSB(vars, error_bands, truth_bands, data_band));
  //studies.push_back(new MeanFrontDEDX_SB(vars, error_bands, truth_bands, data_band));
  //studies.push_back(new ESC_SB(vars, error_bands, truth_bands, data_band));

  for(auto& var: vars) var->InitializeMCHists(error_bands, truth_bands);
  for(auto& var: vars) var->InitializeDATAHists(data_band);
  
  for(auto& var: vars2D) var->InitializeMCHists(error_bands, truth_bands);
  for(auto& var: vars2D) var->InitializeDATAHists(data_band);

  //This is for my saved output trees, just make sure it matches the elements of "studies"
  //std::vector<std::string> sidebandNames = {"Signal_Region"};
  std::vector<std::string> sidebandNames = {"Signal_Region", "Michel_Sideband"};
  //std::vector<std::string> sidebandNames = {"passes_precuts"};

  //Initialize my output tree stuff. if write_tree is false I still get an instance of the class, but it doesn't really do much or add any complexity/time  
  g_OutputTreeManager.Init(MC_TUPLE_OUT_FILE_NAME, write_tree, *options.m_mc->GetChain()->GetTree(), sidebandNames);

  // Loop entries and fill
  try
  {
    
    CVUniverse::SetTruth(false);
    LoopAndFillEventSelection(options.m_mc, error_bands, vars, vars2D, studies, mycuts, model);
    CVUniverse::SetTruth(true);
    LoopAndFillEffDenom(options.m_truth, truth_bands, vars, vars2D, mycuts, model);
    options.PrintMacroConfiguration(argv[0]);
    std::cout << "MC cut summary:\n" << mycuts << "\n";
    mycuts.resetStats();

    CVUniverse::SetTruth(false);
    LoopAndFillData(options.m_data, data_band, vars, vars2D, data_studies, mycuts);
    std::cout << "Data cut summary:\n" << mycuts << "\n";
    
    //Write MC results
    TFile* mcOutDir = TFile::Open(MC_OUT_FILE_NAME, "RECREATE");

    if(!mcOutDir)
    {
      std::cerr << "Failed to open a file named " << MC_OUT_FILE_NAME << " in the current directory for writing histograms.\n";
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
      util::GetFluxIntegral(*error_bands["cv"].front(), var->efficiencyNumerator->hist)->Write((var->GetName() + "_reweightedflux_integrated").c_str());
      //Always use MC number of nucleons for cross section
      auto nNucleons = new TParameter<double>((var->GetName() + "_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons(minZ, maxZ, true, apothem));
      nNucleons->Write();
    }
    
    //Put this here after we finish writing out all the info to the MC histo file
    g_OutputTreeManager.WriteAndClose();
    
    //Write data results
    
    TFile* dataOutDir = TFile::Open(DATA_OUT_FILE_NAME, "RECREATE");
    if(!dataOutDir)
    {
      std::cerr << "Failed to open a file named " << DATA_OUT_FILE_NAME << " in the current directory for writing histograms.\n";
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

  return success;
}
