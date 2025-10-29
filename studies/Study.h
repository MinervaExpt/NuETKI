#ifndef STUDY_H
#define STUDY_H

#include "event/CCNuEEvent.h"

//ROOT includes
#include "TDirectory.h"

class CCNuEEvent;
class CVUniverse;

class Study
{
  public:
    Study() {} //TODO: Any base constructor needed?

    //not sure if splitting up fillSelectedMC and fillSelectedData is the right move, david doesn't do it like this...
    //just seems like, since the event loop already handles the MC and data loops separately I can just call the
    //study->SelectedMC or study->SelectedData separately in each loop, rather than having the same
    //study->Selected() call in each one and having this function figure out which loop its being called from...

    //The only problem with this is that I'll have to redefine the functions in the base Study.h class, because the event loop itself
    //uses pointers to the base class and not the specific sideband derived class... But this would be the same problem for every sideband/Study
    //so might be OK?      
    void SelectedMC(const CVUniverse& univ, const CCNuEEvent& evt, const double weight)
    {
      fillSelectedMC(univ, evt, weight);
    }

    void SelectedData(const CVUniverse& univ, const CCNuEEvent& evt, const double weight)
    {
      fillSelectedData(univ, evt, weight);
    }
    
    void SelectedSignal(const CVUniverse& univ, const CCNuEEvent& evt, const double weight)
    {
      fillSelectedSignal(univ, evt, weight);
    }

    void TruthSignal(const CVUniverse& univ, const CCNuEEvent& evt, const double weight)
    {
      fillTruthSignal(univ, evt, weight);
    }

    //Find Andrew if you need to know how to overload functions for drawing.
    //Only need this when you write a new Study.
    virtual void SaveOrDraw(TFile& outDir) = 0;
    virtual void SaveOrDrawData(TFile& outFile) = 0;
    virtual void SaveOrDrawMC(TFile& outFile) = 0;

  
  private:
    using Hist = PlotUtils::HistWrapper<CVUniverse>;

    virtual void fillSelectedMC(const CVUniverse& univ, const CCNuEEvent& evt, const double weight) = 0;
    virtual void fillSelectedData(const CVUniverse& univ, const CCNuEEvent& evt, const double weight) = 0;
    virtual void fillSelectedSignal(const CVUniverse& univ, const CCNuEEvent& evt, const double weight) = 0;
    virtual void fillTruthSignal(const CVUniverse& univ, const CCNuEEvent& evt, const double weight) = 0;
};

#endif //STUDY_H
