//studies includes
#include "studies/Study.h"

//Mehreen's includes
#include "event/CCNuEEvent.h"
#include "util/Categorized.h"
#include "event/CVUniverse.h"

//c++ includes
#include <functional> //for std::function

class PerEventVarByGENIELabel: public Study
{
  public:
    //PerMichelVarByGENIELabel fills a histogram with 1 entry per Michel with some variable calculated from that Michel.  Your function will get to see the CVUniverse, the CCNuEEvent (= reconstructed Michels), and which Michel it's looping over.
    using reco_t = std::function<double(const CVUniverse&, const CCNuEEvent&)>;

    PerEventVarByGENIELabel(reco_t reco, const std::string& varName, const std::string& varUnits, const int nBins, const double minBin, const double maxBin, const std::map<std::string, std::vector<CVUniverse*>>& univs): Study(), fReco(reco)
    {
      std::map<int, std::string> GENIELabels = {{1, "QE"},
                                                {8, "2p2h"},
                                                {2, "RES"},
                                                {3, "DIS"}};
      m_VarToGENIELabel = new util::Categorized<HIST, int>(varName, varName + " [" + varUnits + "]", GENIELabels, nBins, minBin, maxBin, univs);
    }

    void SaveOrDraw(TDirectory& outDir)
    {
       outDir.cd();
       m_VarToGENIELabel->visit([](HIST& wrapper)
                                {
                                  wrapper.SyncCVHistos();
                                  wrapper.hist->Write();
                                });

       //TODO: You could do plotting here
    }

  private:
    using HIST = PlotUtils::HistWrapper<CVUniverse>;

    reco_t fReco;

    util::Categorized<HIST, int>* m_VarToGENIELabel;

    //Overriding base class functions
    //Do nothing for now...  Good place for data comparisons in the future. 
    void fillSelected(const CVUniverse& univ, const CCNuEEvent& evt, const double weight) {}

    //All of your plots happen here so far.
    void fillSelectedSignal(const CVUniverse& univ, const CCNuEEvent& evt, const double weight)
    {
        
        (*m_VarToGENIELabel)[univ.GetInteractionType()].FillUniverse(&univ, fReco(univ, evt), weight);
    }

    //Do nothing for now...  Good place for efficiency denominators in the future.
    void fillTruthSignal(const CVUniverse& univ, const CCNuEEvent& evt, const double weight) {}
};
