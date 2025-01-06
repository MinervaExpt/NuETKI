//studies includes
#include "studies/Study.h"

//David's includes
#include "event/NeutCands.h"
#include "util/NeutronVariable.h"
#include "util/EventVariable.h"
#include "util/Variable.h"
#include "util/GetBackgroundID.h"
#include "event/CVUniverse.h"

class EMSideBands: public Study
{
  private:
    std::vector<EventVariable*> fEVars_NBlobs;
    std::vector<EventVariable*> fEVars_ENHit;
    std::vector<EventVariable*> fEVars_Sel;
    std::vector<Variable*> fVars_NBlobs;
    std::vector<Variable*> fVars_ENHit;

  public:
    EMSideBands(std::vector<Variable*> vars,
		     std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
		     std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands,
		     std::vector<CVUniverse*>& data_error_bands): Study()
    {
      for (auto& var : vars){
	fVars_NBlobs.push_back(new Variable(false,(var->GetName()+"_NBlobsSB").c_str(), var->GetAxisLabel(), var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
	fVars_ENHit.push_back(new Variable(false,(var->GetName()+"_ENHitSB").c_str(), var->GetAxisLabel(), var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
      }

      std::vector<double> myBlobNBins;
      const double myBlobNBinWidth = 1.;
      for (int whichBin=0; whichBin < 11; ++whichBin) myBlobNBins.push_back(myBlobNBinWidth * whichBin);

      std::vector<double> myBlobEBins;
      const double myBlobEBinWidth = 10.;
      for (int whichBin=0; whichBin < 51; ++whichBin) myBlobEBins.push_back(myBlobEBinWidth * whichBin);

      std::vector<double> myBlobNHitBins;
      const double myBlobNHitBinWidth = 5.;
      for (int whichBin=0; whichBin < 11; ++whichBin) myBlobNHitBins.push_back(myBlobNHitBinWidth * whichBin);

      std::vector<double> myBlobERatBins;
      const double myBlobERatBinWidth = 1.;
      for (int whichBin=0; whichBin < 51; ++whichBin) myBlobERatBins.push_back(myBlobERatBinWidth * whichBin);

      fEVars_Sel = {
	new EventVariable("EMnBlobs","No.", myBlobNBins, &NeutronEvent::GetEMNBlobs),
	new EventVariable("EMBlobE","E [MeV]", myBlobEBins, &NeutronEvent::GetEMBlobE),
	new EventVariable("EMBlobNHit","No.", myBlobNHitBins, &NeutronEvent::GetEMBlobNHits),
	new EventVariable("EMBlobENHitRatio","E/No. [MeV]", myBlobERatBins, &NeutronEvent::GetEMBlobENHitRatio),
	new EventVariable("MaxEFSNeutKE","KE [MeV]", myBlobEBins, &NeutronEvent::GetMaxFSNeutronKE),
      };

      fEVars_NBlobs = {
	new EventVariable("EMnBlobs_NBlobsSB","No.", myBlobNBins, &NeutronEvent::GetEMNBlobs),
	new EventVariable("EMBlobE_NBlobsSB","E [MeV]", myBlobEBins, &NeutronEvent::GetEMBlobE),
	new EventVariable("EMBlobNHit_NBlobsSB","No.", myBlobNHitBins, &NeutronEvent::GetEMBlobNHits),
	new EventVariable("EMBlobENHitRatio_NBlobsSB","E/No. [MeV]", myBlobERatBins, &NeutronEvent::GetEMBlobENHitRatio),
	new EventVariable("MaxEFSNeutKE_NBlobsSB","KE [MeV]", myBlobEBins, &NeutronEvent::GetMaxFSNeutronKE),
      };

      fEVars_ENHit = {
	new EventVariable("EMnBlobs_ENHitSB","No.", myBlobNBins, &NeutronEvent::GetEMNBlobs),
	new EventVariable("EMBlobE_ENHitSB","E [MeV]", myBlobEBins, &NeutronEvent::GetEMBlobE),
	new EventVariable("EMBlobNHit_ENHitSB","No.", myBlobNHitBins, &NeutronEvent::GetEMBlobNHits),
	new EventVariable("EMBlobENHitRatio_ENHitSB","E/No. [MeV]", myBlobERatBins, &NeutronEvent::GetEMBlobENHitRatio),
	new EventVariable("MaxEFSNeutKE_ENHitSB","KE [MeV]", myBlobEBins, &NeutronEvent::GetMaxFSNeutronKE),
      };

      for(auto& var: fVars_NBlobs) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fVars_NBlobs) var->InitializeDATAHists(data_error_bands);

      for(auto& var: fEVars_NBlobs) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fEVars_NBlobs) var->InitializeDATAHists(data_error_bands);

      for(auto& var: fVars_ENHit) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fVars_ENHit) var->InitializeDATAHists(data_error_bands);

      for(auto& var: fEVars_ENHit) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fEVars_ENHit) var->InitializeDATAHists(data_error_bands);

      for(auto& var: fEVars_Sel) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fEVars_Sel) var->InitializeDATAHists(data_error_bands);
    }

    void SaveOrDraw(TDirectory& outDir){
      return;
    }

    void SaveOrDrawMC(TFile& outFile)
    {
      for (auto& var : fVars_NBlobs) var->WriteMC(outFile);
      for (auto& var : fVars_ENHit) var->WriteMC(outFile);
      for (auto& var : fEVars_NBlobs) var->WriteMC(outFile);
      for (auto& var : fEVars_ENHit) var->WriteMC(outFile);
      for (auto& var : fEVars_Sel) var->WriteMC(outFile);
    }

    void SaveOrDrawData(TFile& outFile)
    {
      for (auto& var : fVars_NBlobs) var->WriteData(outFile);
      for (auto& var : fVars_ENHit) var->WriteData(outFile);
      for (auto& var : fEVars_NBlobs) var->WriteData(outFile);
      for (auto& var : fEVars_ENHit) var->WriteData(outFile);
      for (auto& var : fEVars_Sel) var->WriteData(outFile);
    }

  private:
    //Overriding base class functions
    //Do nothing for now...  Good place for data comparisons in the future. 
    void fillSelected(const CVUniverse& univ, const NeutronEvent& evt, const double weight) {
      std::bitset<64> SBStat = evt.GetSideBandStat();
      int intType = evt.GetIntType();
      int tgtType = evt.GetTgtZ();
      int leadBlobType = evt.GetLeadingNeutCand().GetPDGBin();

      if (evt.IsMC()){
	if (evt.IsSignal()){
	  if (SBStat[0] == 0){
	    for (auto& var: fVars_NBlobs){
	      var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(univ), weight); //"Fake data" for closure                                                                          
	      var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValue(univ), weight);
	      (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	      (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	      //(*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    }

	    for (auto& var: fEVars_NBlobs){
	      var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(evt), weight); //"Fake data" for closure                                                                          
	      var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValue(evt), weight);	      
	      (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	      (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	      //(*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	    }
	  }

	  else if (SBStat[1] == 0){
	    for (auto& var: fVars_ENHit){
	      var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(univ), weight); //"Fake data" for closure                                                                          
	      var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValue(univ), weight);
	      (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	      (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	      //(*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    }

	    for (auto& var: fEVars_ENHit){
	      var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(evt), weight); //"Fake data" for closure                                                                          
	      var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValue(evt), weight);	      
	      (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	      (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	      //(*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	    }
	  }

	  else{
	    for (auto& var: fEVars_Sel){
	      var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(evt), weight); //"Fake data" for closure                                                                          
	      var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValue(evt), weight);	      
	      (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	      (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	      //(*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	    }
	  }
	}

	else{
	  int bkgd_ID = -1;	  
	  bkgd_ID = util::GetBackgroundID(univ);

	  if (SBStat[0] == 0){
	    for(auto& var: fVars_NBlobs){
	      var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(univ), weight); //"Fake data" for closure                                                                          
	      (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    //Various breakdowns of selected backgrounds
	      (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	      (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	      //(*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    }

	    for(auto& var: fEVars_NBlobs){
	      var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(evt), weight); //"Fake data" for closure                                                                          
	      (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	    //Various breakdowns of selected backgrounds
	      (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	      (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	      //(*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	    }
	  }

	  else if (SBStat[1] == 0){
	    for(auto& var: fVars_ENHit){
	      var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(univ), weight); //"Fake data" for closure                                                                          
	      (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    //Various breakdowns of selected backgrounds
	      (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	      (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	      //(*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    }

	    for(auto& var: fEVars_ENHit){
	      var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(evt), weight); //"Fake data" for closure                                                                          
	      (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	    //Various breakdowns of selected backgrounds
	      (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	      (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	      //(*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	    }
	  }

	  else{
	    for(auto& var: fEVars_Sel){
	      var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(evt), weight); //"Fake data" for closure                                                                          
	      (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	    //Various breakdowns of selected backgrounds
	      (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	      (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	      //(*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(evt), weight);
	    }
	  }
	}
      }

      else{
	if (SBStat[0] == 0){
	  for (auto& var : fVars_NBlobs){
	    var->dataHist->FillUniverse(&univ, var->GetRecoValue(univ), 1);
	  }

	  for (auto& var : fEVars_NBlobs){
	    var->dataHist->FillUniverse(&univ, var->GetRecoValue(evt), 1);
	  }
	}

	else if (SBStat[1] == 0){
	  for (auto& var : fVars_ENHit){
	    var->dataHist->FillUniverse(&univ, var->GetRecoValue(univ), 1);
	  }

	  for (auto& var : fEVars_ENHit){
	    var->dataHist->FillUniverse(&univ, var->GetRecoValue(evt), 1);
	  }
	}

	else{
	  for (auto& var : fEVars_Sel){
	    var->dataHist->FillUniverse(&univ, var->GetRecoValue(evt), 1);
	  }
	}
      }
      return;
    }

    //All of your plots happen here so far.
    void fillSelectedSignal(const CVUniverse& univ, const NeutronEvent& evt, const double weight)
    {
      return;
      //(*m_VarToGENIELabel)[univ.GetInteractionType()].FillUniverse(&univ, fReco(univ, evt), weight);
    }

    //Do nothing for now...  Good place for efficiency denominators in the future.
    void fillTruthSignal(const CVUniverse& univ, const NeutronEvent& evt, const double weight) { return; }
};
