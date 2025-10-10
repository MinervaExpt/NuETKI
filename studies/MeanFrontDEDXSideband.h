#include "studies/Study.h"
#include "event/CVUniverse.h"
#include "event/MichelEvent.h"
#include "util/Variable.h"
#include "PlotUtils/Cutter.h"

#include "TDirectory.h"

//class MichelEvent;
//class CVUniverse;

class MeanFrontDEDXSideband: public Study
{
  private:
    std::vector<Variable*> fVars;
    std::vector<Variable2D*> fVars2D;
  
  public:
  //what args do I need for my sideband constructor??
    MeanFrontDEDXSideband(std::vector<Variable*> vars,
		 std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
		 std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands,
		   std::vector<CVUniverse*>& data_error_bands): Study() {

      //Making copies of all my existing variables but appending MeanFrontDEDXSB to distinguish which ones are the sideband
      for (auto& var : vars){
	fVars.push_back(new Variable((var->GetName()+"_MeanFrontDEDXSB").c_str(), var->GetAxisLabel(), var->GetBinVec(), var->GetRecoFunc(), var->GetTrueFunc()));
      }

      //Initialize SB histos to be filled and written
      for (auto& var : fVars){
	var->InitializeMCHists(mc_error_bands, truth_error_bands);
	var->InitializeDATAHists(data_error_bands);    
      }
	//fVars2D = vars2D;
    }

  //not sure if splitting up fillSelectedMC and fillSelectedData is the right move, david doesn't do it like this...
  //just seems like, since the event loop already handles the MC and data loops separately I can just call the
  //study->SelectedMC or study->SelectedData separately in each loop, rather than having the same
  //study->Selected() call in each one and having this function figure out which loop its being called from...

  //The only problem with this is that I'll have to redefine the functions in the base Study.h class, because the event loop itself
  //uses pointers to the base class and not the specific sideband derived class... But this would be the same problem for every sideband/Study
  //so might be OK?
    void SelectedMC(const CVUniverse& univ, const MichelEvent& evt, const double weight)
    {
      fillSelectedMC(univ, evt, weight);
    }

    void SelectedData(const CVUniverse& univ, const MichelEvent& evt, const double weight)
    {
      fillSelectedData(univ, evt, weight);
    }
    
    void SelectedSignal(const CVUniverse& univ, const MichelEvent& evt, const double weight)
    {
      fillSelectedSignal(univ, evt, weight);
    }

    void TruthSignal(const CVUniverse& univ, const MichelEvent& evt, const double weight)
    {
      fillTruthSignal(univ, evt, weight);
    }

    //Find Andrew if you need to know how to overload functions for drawing.
    //Only need this when you write a new Study.
    void SaveOrDraw(TFile& outDir) {
      //for (auto& var : fVars) var->WriteMC(outDir);
      return;
    }
    void SaveOrDrawMC(TFile& outDir) {
      for (auto& var : fVars) var->WriteMC(outDir);
      //return;
    }
    void SaveOrDrawData(TFile& outDir) {      
      for (auto& var : fVars) var->WriteData(outDir);
      //return;
    }    
  
  private:
    using Hist = PlotUtils::HistWrapper<CVUniverse>;

      void fillSelectedMC(const CVUniverse& univ, const MichelEvent& evt, const double weight) {
      if (univ.GetMeanFrontdEdx() > 2.4){
	g_OutputTreeManager.Fill("MeanFrontDEDX_Sideband", evt.entryNumber); //add this event to the michel sb tree of my output tree, selectionCategory is already set earlier
	for (auto& var: fVars) var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(univ), weight);
	//const bool isSignal = false; //figure out how to check this later, I don't have michelcuts to do it for me :/
	const auto vertex = univ.GetTrueVertex();
        const bool inApothem = (fabs(vertex.y()) < (-1./sqrt(3.))*fabs(vertex.x()) + 2.*850/sqrt(3.)) && (fabs(vertex.x()) < 850);

	//This checks all of my signal definition requirements (tbh I can also use my selectioncategory from g_outputManager for this, plus it uses michelcuts.issignal
	if (abs(univ.GetTruthNuPDG())==12 && univ.GetCurrent()==1 && univ.GetElectronEnergyTrue()>=2.5 && univ.GetHasSignalFSProton()==1 && univ.GetHasFSMeson()==0 && univ.GetHasFSPhoton()==0 && univ.GetTrueVertex().z() >= 5980 && univ.GetTrueVertex().z() <= 8422 && inApothem) {
	//if (isSignal) { 
	  for (auto& var: fVars){
	    var->efficiencyNumerator->FillUniverse(&univ, var->GetTrueValue(univ), weight);
            var->migration->FillUniverse(&univ, var->GetRecoValue(univ), var->GetTrueValue(univ), weight);
            var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValue(univ), weight);
	  }
	  //for(auto& var: fVars2D) var->efficiencyNumerator->FillUniverse(&univ, var->GetTrueValueX(univ) var->GetTrueValueY(univ), weight);
  	}
	else{
	  int bkgd_ID = -1;
          if (abs(univ.GetTruthNuPDG())==12 && univ.GetCurrent()==1 && univ.GetHasFSMeson()) { bkgd_ID = 0; }
          else if (abs(univ.GetTruthNuPDG())==12 && univ.GetCurrent()==1) { bkgd_ID = 1; }
          else if (univ.GetCurrent()==2 && univ.GetHasFSPi0()==1) { bkgd_ID = 2; }
          else if (abs(univ.GetTruthNuPDG())==14 && univ.GetCurrent()==1 && univ.GetHasFSPi0()==1) { bkgd_ID = 3; }
	  for(auto& var: fVars) (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	  //for(auto& var: fVars2D) (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	}
      }
    }


    void fillSelectedData(const CVUniverse& univ, const MichelEvent& evt, const double weight) {
      for (auto& var: fVars){
	if (univ.GetMeanFrontdEdx() > 2.4){
	  var->dataHist->FillUniverse(&univ, var->GetRecoValue(univ), 1); //should never be weighting data
        }
      }
    }

  
    void fillSelectedSignal(const CVUniverse& univ, const MichelEvent& evt, const double weight) {
      for (auto& var: fVars){
	if (univ.GetMeanFrontdEdx() > 2.4){
	  var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValue(univ), weight);
        }
      }

    }

    void fillTruthSignal(const CVUniverse& univ, const MichelEvent& evt, const double weight) {
      //std::cout << "carlos testing fillTruthSignal\n" << std::endl;
      for (auto& var: fVars){
	if (univ.GetMeanFrontdEdx() > 2.4){
	  var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValue(univ), weight);
        }
      }
    }
};
