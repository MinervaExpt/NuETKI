//studies includes
#include "studies/Study.h"

//David's includes
#include "util/Variable.h"
#include "util/Variable2D.h"
//#include "util/GetBackgroundID.h"
#include "util/Categorized.h"
#include "event/CVUniverse.h"

//Renamed to PreRecoil, includes all target breakdowns needed to fit plastic backgrounds etc.
class PreRecoil: public Study
{
  private:
    std::vector<Variable*> fVars;
    std::vector<util::Categorized<Variable, int>*> fVars_US_ByTgt;//US Plastic inside target
    std::vector<util::Categorized<Variable, int>*> fVars_DS_ByTgt;//DS Plastic inside target
    std::vector<util::Categorized<Variable, int>*> fVars_US_Post_ByTgt;//See above US defined object, but before recoil cut
    std::vector<util::Categorized<Variable, int>*> fVars_DS_Post_ByTgt;//See immediately above, but DS Plastic
    std::vector<util::Categorized<Variable, int>*> fVars_ByTgt;
    util::Categorized<Variable, int>* fVtxTest_US_Post_ByTgt;
    util::Categorized<Variable, int>* fVtxTest_US_ByTgt;
    util::Categorized<Variable, int>* fVtxTest_DS_Post_ByTgt;
    util::Categorized<Variable, int>* fVtxTest_DS_ByTgt;
    util::Categorized<Variable2D, int>* fVtxTest2D_US_Post_ByTgt;
    util::Categorized<Variable2D, int>* fVtxTest2D_US_ByTgt;
    util::Categorized<Variable2D, int>* fVtxTest2D_DS_Post_ByTgt;
    util::Categorized<Variable2D, int>* fVtxTest2D_DS_ByTgt;
    std::map<int, util::Categorized<Variable,int>*> fVtx_US_ByTgt;
    std::map<int, util::Categorized<Variable,int>*> fVtx_DS_ByTgt;
    //std::map<int, util::Categorized<Variable2D,int>*> fVtx2D_US_ByTgt;
    //std::map<int, util::Categorized<Variable2D,int>*> fVtx2D_DS_ByTgt;
    std::vector<Variable2D*> fVars2D;
    std::vector<util::Categorized<Variable2D, int>*> fVars2D_ByTgt;
    std::vector<util::Categorized<Variable2D, int>*> fVars2D_US_ByTgt;
    std::vector<util::Categorized<Variable2D, int>*> fVars2D_DS_ByTgt;
    std::vector<util::Categorized<Variable2D, int>*> fVars2D_US_Post_ByTgt;
    std::vector<util::Categorized<Variable2D, int>*> fVars2D_DS_Post_ByTgt;
    std::map<int,Variable*> fRecoilBinned;
    std::map<int, util::Categorized<Variable, int>*> fRecoilBinned_ByTgt;
    std::map<int, util::Categorized<Variable, int>*> fRecoilBinned_US_ByTgt;
    std::map<int, util::Categorized<Variable, int>*> fRecoilBinned_DS_ByTgt;
    std::map<int,std::map<int, util::Categorized<Variable,int>*>> fVtxBinned_US_ByTgt;
    std::map<int,std::map<int, util::Categorized<Variable,int>*>> fVtxBinned_DS_ByTgt;

    TString fFVregionName;
    bool fSplitRecoil;
    bool fDoNeutronCuts;
    bool fDoVtx;
    int fTgtID;

  public:
    PreRecoil(std::vector<Variable*> vars,
		     std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
		     std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands,
	      std::vector<CVUniverse*>& data_error_bands, bool splitRecoil, bool doNeutronCuts, TString FVregionName, int TgtID, bool doVtx): Study(), fSplitRecoil(splitRecoil), fDoNeutronCuts(doNeutronCuts), fFVregionName(FVregionName), fTgtID(TgtID), fDoVtx(doVtx)
    {

      if (fSplitRecoil){
	//Constructing recoil variable in bins of pt/pz. Just make a map of it all.
	std::vector<double> myRecoilBins;
	const double myRecoilBinWidth = 1.0/50.;
	for(int whichBin = 0; whichBin < 51; ++whichBin) myRecoilBins.push_back(myRecoilBinWidth * whichBin);

	int nBins = 14;//Hard-coded from bins from Amit's analysis... Not appropriate binning given my cuts, but going to get the structure going first.
	int lowBin = -1;//Default return value. hard-coded as well.

	for (int iBin=lowBin; iBin < nBins; ++iBin){
	  std::string binName = std::to_string(iBin);
	  if (iBin == lowBin) binName = "lost";
	  if (fFVregionName.Contains("Target")){
	    if (!fDoVtx){
	      fRecoilBinned_ByTgt[iBin] = new util::Categorized<Variable, int>("", "ByTgt", false, ("recoilE_bin_"+binName).c_str(),"Recoil E [GeV]", util::TgtCodeList[fTgtID], myRecoilBins, &CVUniverse::GetDANRecoilEnergyGeV);
	      fRecoilBinned_US_ByTgt[iBin] = new util::Categorized<Variable, int>("", "ByTgt", false, ("recoilE_bin_"+binName+"_InnerUSPlastic").c_str(),"Recoil E [GeV]", util::TgtCodeList[fTgtID], myRecoilBins, &CVUniverse::GetDANRecoilEnergyGeV);
	      fRecoilBinned_DS_ByTgt[iBin] = new util::Categorized<Variable, int>("", "ByTgt", false, ("recoilE_bin_"+binName+"_InnerDSPlastic").c_str(),"Recoil E [GeV]", util::TgtCodeList[fTgtID], myRecoilBins, &CVUniverse::GetDANRecoilEnergyGeV);
	      fRecoilBinned_ByTgt[iBin]->visit([&mc_error_bands, &truth_error_bands](Variable& var)
					       {
						 var.InitializeMCHists(mc_error_bands, truth_error_bands);
					       });
	      fRecoilBinned_ByTgt[iBin]->visit([&data_error_bands](Variable& var)
					       {
						 var.InitializeDATAHists(data_error_bands);
					       });
	      fRecoilBinned_US_ByTgt[iBin]->visit([&mc_error_bands, &truth_error_bands](Variable& var)
						  {
						    var.InitializeMCHists(mc_error_bands, truth_error_bands);
						  });
	      fRecoilBinned_US_ByTgt[iBin]->visit([&data_error_bands](Variable& var)
						  {
						    var.InitializeDATAHists(data_error_bands);
						  });
	      fRecoilBinned_DS_ByTgt[iBin]->visit([&mc_error_bands, &truth_error_bands](Variable& var)
						  {
						    var.InitializeMCHists(mc_error_bands, truth_error_bands);
						  });
	      fRecoilBinned_DS_ByTgt[iBin]->visit([&data_error_bands](Variable& var)
						  {
						    var.InitializeDATAHists(data_error_bands);
						  });
	    }
	  }
	  else{
	    fRecoilBinned[iBin] = new Variable(false,("recoilE_bin_"+binName).c_str(), "Recoil E [GeV]", myRecoilBins, &CVUniverse::GetDANRecoilEnergyGeV);
	    fRecoilBinned[iBin]->InitializeMCHists(mc_error_bands, truth_error_bands);
	    fRecoilBinned[iBin]->InitializeDATAHists(data_error_bands);
	  }
	}
      }

      std::vector<double> coarseVtxBinsUS;//Modified to be a number of planes instead of vertex Z values.
      for (int iBin=0; iBin<11;++iBin){
	coarseVtxBinsUS.push_back((double)iBin-10.5);
      }

      std::vector<double> coarseVtxBinsDS;//Modified to be a number of planes instead of vertex Z values.
      for (int iBin=11; iBin<22;++iBin){
	coarseVtxBinsDS.push_back((double)iBin-10.5);
      }

      fVars_ByTgt = {};
      fVtx_US_ByTgt = {};
      fVtx_DS_ByTgt = {};
      for (auto& var : vars){
	fVars.push_back(new Variable(false, (var->GetName()+"_PreRecoilCut").c_str(), var->GetAxisLabel(), var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
	if (fFVregionName.Contains("Target")){
	  fVars.back()->SetFillVar(false);
	  TString nameCheck = var->GetName();
	  if (nameCheck=="vtxZ"){
	    if (fSplitRecoil){
	      if (fDoVtx){
		int nBins = 14;//Hard-coded from bins from Amit's analysis... Not appropriate binning given my cuts, but going to get the structure going first.
		int lowBin = -1;//Default return value. hard-coded as well.
	      
		for (int iBin=lowBin; iBin < nBins; ++iBin){
		  std::string binName = std::to_string(iBin);
		  if (iBin == lowBin) binName = "lost";
		  for (auto& tgt: util::TgtCodeList[fTgtID]){
		    fVtxBinned_US_ByTgt[iBin][tgt.first]=new util::Categorized<Variable, int>(("ByTgt_"+tgt.second).c_str(), "US_ByType", false, (var->GetName()+"_bin_"+binName+"_ByTgt_"+tgt.second+"_PreRecoilCut").c_str(),var->GetAxisLabel().c_str(), util::TgtTypeList2,coarseVtxBinsUS,var->GetRecoFunc(),var->GetTrueFunc());
		    fVtxBinned_US_ByTgt[iBin][tgt.first]->visit([&mc_error_bands, &truth_error_bands](Variable& vtx)
								{
								  vtx.InitializeMCHists(mc_error_bands, truth_error_bands);
								});
		    fVtxBinned_US_ByTgt[iBin][tgt.first]->visit([&data_error_bands](Variable& vtx)
								{
								  vtx.InitializeDATAHists(data_error_bands);
								});
		    fVtxBinned_DS_ByTgt[iBin][tgt.first]=new util::Categorized<Variable, int>(("ByTgt_"+tgt.second).c_str(), "DS_ByType", false, (var->GetName()+"_bin_"+binName+"_ByTgt_"+tgt.second+"_PreRecoilCut").c_str(),var->GetAxisLabel().c_str(), util::TgtTypeList2,coarseVtxBinsDS,var->GetRecoFunc(),var->GetTrueFunc());
		    fVtxBinned_DS_ByTgt[iBin][tgt.first]->visit([&mc_error_bands, &truth_error_bands](Variable& vtx)
								{
								  vtx.InitializeMCHists(mc_error_bands, truth_error_bands);
								});
		    fVtxBinned_DS_ByTgt[iBin][tgt.first]->visit([&data_error_bands](Variable& vtx)
								{
								  vtx.InitializeDATAHists(data_error_bands);
								});
		  }
		}
	      }
	    }
	    else {
	      fVtxTest_US_Post_ByTgt = new util::Categorized<Variable, int>(var->GetDirectoryName(), "ByTgt", false, (var->GetName()+"_OuterUSPlastic").c_str(),var->GetAxisLabel().c_str(), util::TgtCodeList[fTgtID],coarseVtxBinsUS,var->GetRecoFunc(),var->GetTrueFunc());
	      fVtxTest_US_ByTgt = new util::Categorized<Variable, int>(var->GetDirectoryName(), "ByTgt", false, (var->GetName()+"_OuterUSPlastic_PreRecoilCut").c_str(),var->GetAxisLabel().c_str(), util::TgtCodeList[fTgtID],coarseVtxBinsUS,var->GetRecoFunc(),var->GetTrueFunc());
	      fVtxTest_DS_Post_ByTgt = new util::Categorized<Variable, int>(var->GetDirectoryName(), "ByTgt", false, (var->GetName()+"_OuterDSPlastic").c_str(),var->GetAxisLabel().c_str(), util::TgtCodeList[fTgtID],coarseVtxBinsDS,var->GetRecoFunc(),var->GetTrueFunc());
	      fVtxTest_DS_ByTgt = new util::Categorized<Variable, int>(var->GetDirectoryName(), "ByTgt", false, (var->GetName()+"_OuterDSPlastic_PreRecoilCut").c_str(),var->GetAxisLabel().c_str(), util::TgtCodeList[fTgtID],coarseVtxBinsDS,var->GetRecoFunc(),var->GetTrueFunc());

	      fVtxTest2D_US_Post_ByTgt = new util::Categorized<Variable2D, int>(var->GetDirectoryName(), "ByTgt", false, (var->GetName()+"_v_pT_OuterUSPlastic").c_str(), util::TgtCodeList[fTgtID], *fVars[0], (*fVtxTest_US_Post_ByTgt)[0]);
	      fVtxTest2D_US_ByTgt = new util::Categorized<Variable2D, int>(var->GetDirectoryName(), "ByTgt", false, (var->GetName()+"_v_pT_OuterUSPlastic_PreRecoilCut").c_str(), util::TgtCodeList[fTgtID], *fVars[0], (*fVtxTest_US_ByTgt)[0]);
	      fVtxTest2D_DS_Post_ByTgt = new util::Categorized<Variable2D, int>(var->GetDirectoryName(), "ByTgt", false, (var->GetName()+"_v_pT_OuterDSPlastic").c_str(), util::TgtCodeList[fTgtID], *fVars[0], (*fVtxTest_DS_Post_ByTgt)[0]);
	      fVtxTest2D_DS_ByTgt = new util::Categorized<Variable2D, int>(var->GetDirectoryName(), "ByTgt", false, (var->GetName()+"_v_pT_OuterDSPlastic_PreRecoilCut").c_str(), util::TgtCodeList[fTgtID], *fVars[0], (*fVtxTest_DS_ByTgt)[0]);

	      for (auto& tgt: util::TgtCodeList[fTgtID]){
		fVtx_US_ByTgt[tgt.first]=new util::Categorized<Variable, int>(("ByTgt_"+tgt.second).c_str(), "US_ByType", false, (var->GetName()+"_ByTgt_"+tgt.second+"_PreRecoilCut").c_str(),var->GetAxisLabel().c_str(), util::TgtTypeList,coarseVtxBinsUS,var->GetRecoFunc(),var->GetTrueFunc());
		fVtx_DS_ByTgt[tgt.first]=new util::Categorized<Variable, int>(("ByTgt_"+tgt.second).c_str(), "DS_ByType", false, (var->GetName()+"_ByTgt_"+tgt.second+"_PreRecoilCut").c_str(),var->GetAxisLabel().c_str(), util::TgtTypeList,coarseVtxBinsDS,var->GetRecoFunc(),var->GetTrueFunc());
		//fVtx2D_US_ByTgt[tgt.first]=new util::Categorized<Variable2D, int>(("ByTgt_"+tgt.second).c_str(), "US_ByType", false, (var->GetName()+"2D_ByTgt_"+tgt.second+"_PreRecoilCut").c_str(),*fVars[0],*fVtx_US_ByTgt[tgt.first]);
		//fVtx2D_DS_ByTgt[tgt.first]=new util::Categorized<Variable2D, int>(("ByTgt_"+tgt.second).c_str(), "DS_ByType", false, (var->GetName()+"2D_ByTgt_"+tgt.second+"_PreRecoilCut").c_str(),*fVars[0],*fVtx_DS_ByTgt[tgt.first]);
	      }
	    }
	  }
	  else if ((nameCheck == "pTmu" || nameCheck == "recoilE") && fDoNeutronCuts){
	    fVars_ByTgt.push_back(new util::Categorized<Variable, int>(var->GetDirectoryName(), "ByTgt", false, (var->GetName()+"_PreRecoilCut").c_str(),var->GetAxisLabel().c_str(), util::TgtCodeList[fTgtID],var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
	    fVars_US_ByTgt.push_back(new util::Categorized<Variable, int>(var->GetDirectoryName(), "ByTgt", false, (var->GetName()+"_InnerUSPlastic_PreRecoilCut").c_str(),var->GetAxisLabel().c_str(), util::TgtCodeList[fTgtID],var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
	    fVars_DS_ByTgt.push_back(new util::Categorized<Variable, int>(var->GetDirectoryName(), "ByTgt", false, (var->GetName()+"_InnerDSPlastic_PreRecoilCut").c_str(),var->GetAxisLabel().c_str(), util::TgtCodeList[fTgtID],var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
	    fVars_US_Post_ByTgt.push_back(new util::Categorized<Variable, int>(var->GetDirectoryName(), "ByTgt", false, (var->GetName()+"_InnerUSPlastic").c_str(),var->GetAxisLabel().c_str(), util::TgtCodeList[fTgtID],var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
	    fVars_DS_Post_ByTgt.push_back(new util::Categorized<Variable, int>(var->GetDirectoryName(), "ByTgt", false, (var->GetName()+"_InnerDSPlastic").c_str(),var->GetAxisLabel().c_str(), util::TgtCodeList[fTgtID],var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
	  }
	  else if (fDoNeutronCuts){
	    fVars_ByTgt.push_back(new util::Categorized<Variable, int>(var->GetDirectoryName(), "ByTgt", false, (var->GetName()+"_PreRecoilCut").c_str(),var->GetAxisLabel().c_str(), util::TgtCodeList[fTgtID],var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
	  }
	}
      }
      
      fVars2D = {};
      fVars2D_ByTgt = {};
      //Look to main code for other 2D variable initializations.
      if (!fFVregionName.Contains("Target")) fVars2D.push_back(new Variable2D(false,"recoil_v_pT_PreRecoilCut",*vars[0],*vars[vars.size()-3]));
      if (!fDoNeutronCuts){
	for (auto& var: fVars) var->SetFillVar(false);
	if (fFVregionName.Contains("Target")){
	  //Skips these if the only plots we want are the vertex plots.
	  if (!fDoVtx){
	    fVars2D_ByTgt.push_back(new util::Categorized<Variable2D, int>("", "ByTgt", true, "pmu2D_PreRecoilCut", util::TgtCodeList[fTgtID], *fVars[1], *fVars[0]));
	    fVars2D_US_ByTgt.push_back(new util::Categorized<Variable2D, int>("", "ByTgt", true, "pmu2D_InnerUSPlastic_PreRecoilCut", util::TgtCodeList[fTgtID], *fVars[1], *fVars[0]));
	    fVars2D_DS_ByTgt.push_back(new util::Categorized<Variable2D, int>("", "ByTgt", true, "pmu2D_InnerDSPlastic_PreRecoilCut", util::TgtCodeList[fTgtID], *fVars[1], *fVars[0]));
	    fVars2D_US_Post_ByTgt.push_back(new util::Categorized<Variable2D, int>("", "ByTgt", true, "pmu2D_InnerUSPlastic", util::TgtCodeList[fTgtID], *fVars[1], *fVars[0]));
	    fVars2D_DS_Post_ByTgt.push_back(new util::Categorized<Variable2D, int>("", "ByTgt", true, "pmu2D_InnerDSPlastic", util::TgtCodeList[fTgtID], *fVars[1], *fVars[0]));
	  }
	}
	else{
	  fVars2D.push_back(new Variable2D(true,"pmu2D_PreRecoilCut",*fVars[1],*fVars[0]));
	}
      }

      for(auto& var: fVars) if (var->IsFill()) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fVars) if (var->IsFill()) var->InitializeDATAHists(data_error_bands);

      for(auto& tgt: fVars_ByTgt){
	tgt->visit([&mc_error_bands, &truth_error_bands](Variable& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fVars_ByTgt){
	tgt->visit([&data_error_bands](Variable& var)
		   {
		     var.InitializeDATAHists(data_error_bands);
		   });
      }

      for(auto& tgt: fVars_US_ByTgt){
	tgt->visit([&mc_error_bands, &truth_error_bands](Variable& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fVars_US_ByTgt){
	tgt->visit([&data_error_bands](Variable& var)
		   {
		     var.InitializeDATAHists(data_error_bands);
		   });
      }

      for(auto& tgt: fVars_DS_ByTgt){
	tgt->visit([&mc_error_bands, &truth_error_bands](Variable& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fVars_DS_ByTgt){
	tgt->visit([&data_error_bands](Variable& var)
		   {
		     var.InitializeDATAHists(data_error_bands);
		   });
      }

      for(auto& tgt: fVars_US_Post_ByTgt){
	tgt->visit([&mc_error_bands, &truth_error_bands](Variable& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fVars_US_Post_ByTgt){
	tgt->visit([&data_error_bands](Variable& var)
		   {
		     var.InitializeDATAHists(data_error_bands);
		   });
      }

      for(auto& tgt: fVars_DS_Post_ByTgt){
	tgt->visit([&mc_error_bands, &truth_error_bands](Variable& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fVars_DS_Post_ByTgt){
	tgt->visit([&data_error_bands](Variable& var)
		   {
		     var.InitializeDATAHists(data_error_bands);
		   });
      }

      if (fFVregionName.Contains("Target")){
	fVtxTest_US_Post_ByTgt->visit([&mc_error_bands, &truth_error_bands](Variable& var){
	    var.InitializeMCHists(mc_error_bands, truth_error_bands);
	  });
	fVtxTest_US_Post_ByTgt->visit([&data_error_bands](Variable& var){
	    var.InitializeDATAHists(data_error_bands);
	  });
	fVtxTest2D_US_Post_ByTgt->visit([&mc_error_bands, &truth_error_bands](Variable2D& var){
	    var.InitializeMCHists(mc_error_bands, truth_error_bands);
	  });
	fVtxTest2D_US_Post_ByTgt->visit([&data_error_bands](Variable2D& var){
	    var.InitializeDATAHists(data_error_bands);
	  });
	
	fVtxTest_US_ByTgt->visit([&mc_error_bands, &truth_error_bands](Variable& var){
	    var.InitializeMCHists(mc_error_bands, truth_error_bands);
	  });
	fVtxTest_US_ByTgt->visit([&data_error_bands](Variable& var){
	    var.InitializeDATAHists(data_error_bands);
	  });
	fVtxTest2D_US_ByTgt->visit([&mc_error_bands, &truth_error_bands](Variable2D& var){
	    var.InitializeMCHists(mc_error_bands, truth_error_bands);
	  });
	fVtxTest2D_US_ByTgt->visit([&data_error_bands](Variable2D& var){
	    var.InitializeDATAHists(data_error_bands);
	  });
	
	fVtxTest_DS_Post_ByTgt->visit([&mc_error_bands, &truth_error_bands](Variable& var){
	    var.InitializeMCHists(mc_error_bands, truth_error_bands);
	  });
	fVtxTest_DS_Post_ByTgt->visit([&data_error_bands](Variable& var){
	    var.InitializeDATAHists(data_error_bands);
	  });
	fVtxTest2D_DS_Post_ByTgt->visit([&mc_error_bands, &truth_error_bands](Variable2D& var){
	    var.InitializeMCHists(mc_error_bands, truth_error_bands);
	  });
	fVtxTest2D_DS_Post_ByTgt->visit([&data_error_bands](Variable2D& var){
	    var.InitializeDATAHists(data_error_bands);
	  });
	
	fVtxTest_DS_ByTgt->visit([&mc_error_bands, &truth_error_bands](Variable& var){
	    var.InitializeMCHists(mc_error_bands, truth_error_bands);
	  });
	fVtxTest_DS_ByTgt->visit([&data_error_bands](Variable& var){
	    var.InitializeDATAHists(data_error_bands);
	  });
	fVtxTest2D_DS_ByTgt->visit([&mc_error_bands, &truth_error_bands](Variable2D& var){
	    var.InitializeMCHists(mc_error_bands, truth_error_bands);
	  });
	fVtxTest2D_DS_ByTgt->visit([&data_error_bands](Variable2D& var){
	    var.InitializeDATAHists(data_error_bands);
	  });
      }

      for(auto& tgt: fVtx_US_ByTgt){
	tgt.second->visit([&mc_error_bands, &truth_error_bands](Variable& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fVtx_US_ByTgt){
	tgt.second->visit([&data_error_bands](Variable& var)
		   {
		     var.InitializeDATAHists(data_error_bands);
		   });
      }

      for(auto& tgt: fVtx_DS_ByTgt){
	tgt.second->visit([&mc_error_bands, &truth_error_bands](Variable& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fVtx_DS_ByTgt){
	tgt.second->visit([&data_error_bands](Variable& var)
		   {
		     var.InitializeDATAHists(data_error_bands);
		   });
      }

      for(auto& var: fVars2D) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fVars2D) var->InitializeDATAHists(data_error_bands);

      for(auto& tgt: fVars2D_ByTgt){
	tgt->visit([&mc_error_bands, &truth_error_bands](Variable2D& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fVars2D_ByTgt){
	tgt->visit([&data_error_bands](Variable2D& var)
		   {
		     var.InitializeDATAHists(data_error_bands);
		   });
      }
      for(auto& tgt: fVars2D_US_ByTgt){
	tgt->visit([&mc_error_bands, &truth_error_bands](Variable2D& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fVars2D_US_ByTgt){
	tgt->visit([&data_error_bands](Variable2D& var)
		   {
		     var.InitializeDATAHists(data_error_bands);
		   });
      }
      for(auto& tgt: fVars2D_DS_ByTgt){
	tgt->visit([&mc_error_bands, &truth_error_bands](Variable2D& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fVars2D_DS_ByTgt){
	tgt->visit([&data_error_bands](Variable2D& var)
		   {
		     var.InitializeDATAHists(data_error_bands);
		   });
      }
      for(auto& tgt: fVars2D_US_Post_ByTgt){
	tgt->visit([&mc_error_bands, &truth_error_bands](Variable2D& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fVars2D_US_Post_ByTgt){
	tgt->visit([&data_error_bands](Variable2D& var)
		   {
		     var.InitializeDATAHists(data_error_bands);
		   });
      }
      for(auto& tgt: fVars2D_DS_Post_ByTgt){
	tgt->visit([&mc_error_bands, &truth_error_bands](Variable2D& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fVars2D_DS_Post_ByTgt){
	tgt->visit([&data_error_bands](Variable2D& var)
		   {
		     var.InitializeDATAHists(data_error_bands);
		   });
      }
    }

    void SaveOrDraw(TDirectory& outDir){
      return;
    }

    void SaveOrDrawMC(TFile& outFile)
    {
      for (auto& recoil : fRecoilBinned) recoil.second->WriteMC(outFile);
      for (auto& recoil : fRecoilBinned_ByTgt){
	recoil.second->visit([&outFile](Variable& var)
			     {
			       var.WriteMC(outFile);
			     });
      }
      for (auto& recoil : fRecoilBinned_US_ByTgt){
	recoil.second->visit([&outFile](Variable& var)
			     {
			       var.WriteMC(outFile);
			     });
      }
      for (auto& recoil : fRecoilBinned_DS_ByTgt){
	recoil.second->visit([&outFile](Variable& var)
			     {
			       var.WriteMC(outFile);
			     });
      }
      for (auto& var : fVars) if (var->IsFill()) var->WriteMC(outFile);
      for(auto& tgt: fVars_ByTgt){
	tgt->visit([&outFile](Variable& var)
		   {
		     var.WriteMC(outFile);
		   });
      }
      for(auto& tgt: fVars_US_ByTgt){
	tgt->visit([&outFile](Variable& var)
		   {
		     var.WriteMC(outFile);
		   });
      }
      for(auto& tgt: fVars_DS_ByTgt){
	tgt->visit([&outFile](Variable& var)
		   {
		     var.WriteMC(outFile);
		   });
      }
      for(auto& tgt: fVars_US_Post_ByTgt){
	tgt->visit([&outFile](Variable& var)
		   {
		     var.WriteMC(outFile);
		   });
      }
      for(auto& tgt: fVars_DS_Post_ByTgt){
	tgt->visit([&outFile](Variable& var)
		   {
		     var.WriteMC(outFile);
		   });
      }

      if (fFVregionName.Contains("Target")){
	fVtxTest_US_Post_ByTgt->visit([&outFile](Variable& var){
	    var.WriteMC(outFile);
	  });
	fVtxTest2D_US_Post_ByTgt->visit([&outFile](Variable2D& var){
	    var.WriteMC(outFile);
	  });
	
	fVtxTest_US_ByTgt->visit([&outFile](Variable& var){
	    var.WriteMC(outFile);
	  });
	fVtxTest2D_US_ByTgt->visit([&outFile](Variable2D& var){
	    var.WriteMC(outFile);
	  });
	
	fVtxTest_DS_Post_ByTgt->visit([&outFile](Variable& var){
	    var.WriteMC(outFile);
	  });
	fVtxTest2D_DS_Post_ByTgt->visit([&outFile](Variable2D& var){
	    var.WriteMC(outFile);
	  });
	
	fVtxTest_DS_ByTgt->visit([&outFile](Variable& var){
	    var.WriteMC(outFile);
	  });
	fVtxTest2D_DS_ByTgt->visit([&outFile](Variable2D& var){
	    var.WriteMC(outFile);
	  });
      }

      for(auto& tgt: fVtx_US_ByTgt){
	tgt.second->visit([&outFile](Variable& var)
		   {
		     var.WriteMC(outFile);
		   });
      }
      for(auto& tgt: fVtx_DS_ByTgt){
	tgt.second->visit([&outFile](Variable& var)
		   {
		     var.WriteMC(outFile);
		   });
      }
      for(auto& recoil: fVtxBinned_US_ByTgt){
	for(auto& tgt: recoil.second){
	  tgt.second->visit([&outFile](Variable& var)
			    {
			      var.WriteMC(outFile);
			    });
	}
      }
      for(auto& recoil: fVtxBinned_DS_ByTgt){
	for(auto& tgt: recoil.second){
	  tgt.second->visit([&outFile](Variable& var)
			    {
			      var.WriteMC(outFile);
			    });
	}
      }

      for (auto& var : fVars2D) var->WriteMC(outFile);
      for(auto& tgt: fVars2D_ByTgt){
	tgt->visit([&outFile](Variable2D& var)
		   {
		     var.WriteMC(outFile);
		   });
      }
      for(auto& tgt: fVars2D_US_ByTgt){
	tgt->visit([&outFile](Variable2D& var)
		   {
		     var.WriteMC(outFile);
		   });
      }
      for(auto& tgt: fVars2D_DS_ByTgt){
	tgt->visit([&outFile](Variable2D& var)
		   {
		     var.WriteMC(outFile);
		   });
      }
      for(auto& tgt: fVars2D_US_Post_ByTgt){
	tgt->visit([&outFile](Variable2D& var)
		   {
		     var.WriteMC(outFile);
		   });
      }
      for(auto& tgt: fVars2D_DS_Post_ByTgt){
	tgt->visit([&outFile](Variable2D& var)
		   {
		     var.WriteMC(outFile);
		   });
      }
    }

    void SaveOrDrawData(TFile& outFile)
    {
      for (auto& recoil : fRecoilBinned) recoil.second->WriteData(outFile);
      for (auto& recoil : fRecoilBinned_ByTgt){
	recoil.second->visit([&outFile](Variable& var)
			     {
			       var.WriteData(outFile);
			     });
      }
      for (auto& recoil : fRecoilBinned_US_ByTgt){
	recoil.second->visit([&outFile](Variable& var)
			     {
			       var.WriteData(outFile);
			     });
      }
      for (auto& recoil : fRecoilBinned_DS_ByTgt){
	recoil.second->visit([&outFile](Variable& var)
			     {
			       var.WriteData(outFile);
			     });
      }
      for (auto& var : fVars) if (var->IsFill()) var->WriteData(outFile);
      for(auto& tgt: fVars_ByTgt){
	tgt->visit([&outFile](Variable& var)
		   {
		     var.WriteData(outFile);
		   });
      }

      if (fFVregionName.Contains("Target")){
	fVtxTest_US_Post_ByTgt->visit([&outFile](Variable& var){
	    var.WriteData(outFile);
	  });
	fVtxTest2D_US_Post_ByTgt->visit([&outFile](Variable2D& var){
	    var.WriteData(outFile);
	  });
	
	fVtxTest_US_ByTgt->visit([&outFile](Variable& var){
	    var.WriteData(outFile);
	  });
	fVtxTest2D_US_ByTgt->visit([&outFile](Variable2D& var){
	  var.WriteData(outFile);
	  });
	
	fVtxTest_DS_Post_ByTgt->visit([&outFile](Variable& var){
	    var.WriteData(outFile);
	  });
	fVtxTest2D_DS_Post_ByTgt->visit([&outFile](Variable2D& var){
	    var.WriteData(outFile);
	  });
	
	fVtxTest_DS_ByTgt->visit([&outFile](Variable& var){
	    var.WriteData(outFile);
	  });
	fVtxTest2D_DS_ByTgt->visit([&outFile](Variable2D& var){
	    var.WriteData(outFile);
	  });
      }

      for(auto& tgt: fVtx_US_ByTgt){
	tgt.second->visit([&outFile](Variable& var)
		   {
		     var.WriteData(outFile);
		   });
      }
      for(auto& tgt: fVtx_DS_ByTgt){
	tgt.second->visit([&outFile](Variable& var)
		   {
		     var.WriteData(outFile);
		   });
      }
      for(auto& recoil: fVtxBinned_US_ByTgt){
	for(auto& tgt: recoil.second){
	  tgt.second->visit([&outFile](Variable& var)
			    {
			      var.WriteData(outFile);
			    });
	}
      }
      for(auto& recoil: fVtxBinned_DS_ByTgt){
	for(auto& tgt: recoil.second){
	  tgt.second->visit([&outFile](Variable& var)
			    {
			      var.WriteData(outFile);
			    });
	}
      }

      for (auto& var : fVars2D) var->WriteData(outFile);
      for(auto& tgt: fVars2D_ByTgt){
	tgt->visit([&outFile](Variable2D& var)
		   {
		     var.WriteData(outFile);
		   });
      }
      for(auto& tgt: fVars2D_US_ByTgt){
	tgt->visit([&outFile](Variable2D& var)
		   {
		     var.WriteData(outFile);
		   });
      }
      for(auto& tgt: fVars2D_DS_ByTgt){
	tgt->visit([&outFile](Variable2D& var)
		   {
		     var.WriteData(outFile);
		   });
      }
      for(auto& tgt: fVars2D_US_Post_ByTgt){
	tgt->visit([&outFile](Variable2D& var)
		   {
		     var.WriteData(outFile);
		   });
      }
      for(auto& tgt: fVars2D_DS_Post_ByTgt){
	tgt->visit([&outFile](Variable2D& var)
		   {
		     var.WriteData(outFile);
		   });
      }
    }
    
  private:
    //Overriding base class functions
    //Do nothing for now...  Good place for data comparisons in the future. 
    void fillSelected(const CVUniverse& univ, const NeutronEvent& evt, const double weight) {
      std::bitset<64> SBStat = evt.GetSideBandStat();
      int intType = evt.GetIntType();
      int intCode = evt.GetIntCode();
      int tgtType = evt.GetTgtZ();
      int leadBlobType = evt.GetLeadingNeutCand().GetPDGBin();
      int iBin = evt.GetBinPTPZ();
      const bool isTgts = (fFVregionName.Contains("Target")) ? true: false;
      std::vector<double> muonMom = {univ.GetMuon4V().X(),univ.GetMuon4V().Y(),univ.GetMuon4V().Z()};

      //Could code this into the sideband definition, but easier to test the idea by hard-coding here.
      //bool isSideband = univ.GetDANRecoilEnergyGeV() >= 0.5 && univ.GetDANRecoilEnergyGeV() <= 1.0;
      bool isSideband = (!SBStat.all());

      // At some point just put this into the event structure itself.
      std::vector<double> vtx = univ.GetVtx();
      double vtx_x = vtx.at(0);
      double vtx_y = vtx.at(1);
      double vtx_z = vtx.at(2);

      int tgtID = util::GetRecoTargetZ(vtx_x,vtx_y,vtx_z);
      int tgtCode = util::GetRecoTargetCode(vtx_x,vtx_y,vtx_z,muonMom);

      //Adding this in to change the checks that are done to get the right bkgd id to match what's in runEventLoop... this is due to a change in utilizng the true tgtCode to match instead of just the material.
      double mc_vtx_x = -999.;
      double mc_vtx_y = -999.;
      double mc_vtx_z = -999.;
      int tgtZ = -999;
      if (evt.IsMC()){
	std::vector<double> mc_vtx = univ.GetTrueVtx();
	mc_vtx_x = mc_vtx.at(0);
	mc_vtx_y = mc_vtx.at(1);
	mc_vtx_z = mc_vtx.at(2);
	tgtZ = univ.GetTargetZ();
      }

      int trueTgtCode = util::GetTrueTgtCode(tgtZ, mc_vtx_x, mc_vtx_y, mc_vtx_z);

      int USTgt = -1;
      int DSTgt = -1;
      //std::cout << fFVregionName << std::endl;
      if (tgtCode == -1  && fFVregionName.Contains("Target")){
	//std::cout << fFVregionName << " changing USTgt and DSTgt." << std::endl;
	USTgt = util::GetUSTgtCode(tgtID,vtx_x,vtx_y,vtx_z,muonMom,fTgtID);
	DSTgt = util::GetDSTgtCode(tgtID,vtx_x,vtx_y,vtx_z,muonMom,fTgtID);
	/*Debugging removed
	if (tgtID >= 10 && USTgt == -999 && tgtType == 1){
	  std::cout << "TGT ID: " << tgtID << std::endl;
	  std::cout << "US Tgt Code: " << USTgt << std::endl;
	  std::cout << "DS Tgt Code: " << DSTgt << std::endl;
	  std::cout << "Event Weight: " << weight << std::endl;
	  if(evt.IsSignal()) std::cout << "Signal" << std::endl;
	  else std::cout << "BKG" << std::endl;
	  if (tgtID == 10) std::cout << "N Planes US: " << util::GetNPlanesUSOfTarget(1111,vtx_z) << std::endl;
	  else if (tgtID == 21) std::cout << "N Planes US: " << util::GetNPlanesUSOfTarget(3333,vtx_z) << std::endl;
	  else if (tgtID == 32) std::cout << "N Planes US: " << util::GetNPlanesUSOfTarget(3333,vtx_z) << std::endl;
	  else if (tgtID == 63) std::cout << "N Planes US: " << util::GetNPlanesUSOfTarget(6666,vtx_z) << std::endl;
	  else if (tgtID == 46) std::cout << "N Planes US: " << util::GetNPlanesUSOfTarget(4444,vtx_z) << std::endl;
	  else if (tgtID == 54) std::cout << "N Planes US: " << util::GetNPlanesUSOfTarget(5555,vtx_z) << std::endl;
	  std::cout << "" << std::endl;
	}
	else if ((tgtID == 0 || tgtID > 11) && DSTgt == -999 && tgtType == 1){
	  std::cout << "TGT ID: " << tgtID << std::endl;
	  std::cout << "US Tgt Code: " << USTgt << std::endl;
	  std::cout << "DS Tgt Code: " << DSTgt << std::endl;
	  std::cout << "Event Weight: " << weight << std::endl;
	  if(evt.IsSignal()) std::cout << "Signal" << std::endl;
	  else std::cout << "BKG" << std::endl;

	  if(tgtID==21) std::cout << "N Planes DS: " << util::GetNPlanesDSOfTarget(1111,vtx_z) << std::endl;
	  else if(tgtID==32) std::cout << "N Planes DS: " << util::GetNPlanesDSOfTarget(2222,vtx_z) << std::endl;
	  else if(tgtID==63) std::cout << "N Planes DS: " << util::GetNPlanesDSOfTarget(3333,vtx_z) << std::endl;
	  else if(tgtID==46) std::cout << "N Planes DS: " << util::GetNPlanesDSOfTarget(6666,vtx_z) << std::endl;
	  else if(tgtID==54) std::cout << "N Planes DS: " << util::GetNPlanesDSOfTarget(4444,vtx_z) << std::endl;
	  else if(tgtID==0) std::cout << "N Planes DS: " << util::GetNPlanesDSOfTarget(5555,vtx_z) << std::endl;
	  std::cout << "" << std::endl;
	}
	*/
      }
      
      if (evt.IsMC()){
	if (evt.IsSignal()){
	  //Could put this inside the isSideband, but this gives flexibility of fitting still that I like.
	  if (fSplitRecoil){
	    if (fRecoilBinned.size() > 0){
	      fRecoilBinned[iBin]->selectedMCReco->FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight); //"Fake data" for closure
	      fRecoilBinned[iBin]->selectedSignalReco->FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	      if (fRecoilBinned[iBin]->IsBroken()){
		(*fRecoilBinned[iBin]->m_SigIntTypeHists)[intType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
		(*fRecoilBinned[iBin]->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
		//(*fRecoilBinned[iBin]->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	      }
	    }

	    if (tgtCode != -1 && fRecoilBinned_ByTgt.size() > 0){
	      (*fRecoilBinned_ByTgt[iBin])[tgtCode].selectedMCReco->FillUniverse(&univ, (*fRecoilBinned_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
	      (*fRecoilBinned_ByTgt[iBin])[tgtCode].selectedSignalReco->FillUniverse(&univ, (*fRecoilBinned_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
	      if ((*fRecoilBinned_ByTgt[iBin])[tgtCode].IsBroken()){
		(*(*fRecoilBinned_ByTgt[iBin])[tgtCode].m_SigIntTypeHists)[intType].FillUniverse(&univ, (*fRecoilBinned_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
		(*(*fRecoilBinned_ByTgt[iBin])[tgtCode].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*fRecoilBinned_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
		//(*(*fRecoilBinned_ByTgt[iBin])[tgtCode].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*fRecoilBinned_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
	      }
	    }
	  }

	  if (isSideband){
	    for (auto& var: fVars){
	      if (var->IsFill()){
		var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(univ), weight); //"Fake data" for closure
		var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValue(univ), weight);
		if (var->IsBroken()){
		  (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
		  (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
		  //(*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
		}
	      }
	    }

	    for (auto& var: fVars2D){
	      var->selectedMCReco->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight); //"Fake data" for closure
	      var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	      (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	      (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	      //(*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    }

	    //Only fill for non-interstitial plastics
	    if(tgtCode != -1){
	      for (auto& var: fVars_ByTgt){
		(*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
		(*var)[tgtCode].selectedSignalReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		if ((*var)[tgtCode].IsBroken()){
		  (*(*var)[tgtCode].m_SigIntTypeHists)[intType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		  (*(*var)[tgtCode].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		  //(*(*var)[tgtCode].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		}
	      }

	      for (auto& var: fVars2D_ByTgt){
		(*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight); //"Fake data" for closure
		(*var)[tgtCode].selectedSignalReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
		(*(*var)[tgtCode].m_SigIntTypeHists)[intType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
		(*(*var)[tgtCode].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
		//(*(*var)[tgtCode].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
	      }
	    }
	    else if (tgtCode == -1){
	      if(USTgt > 0){
		(*fVtxTest_US_ByTgt)[USTgt].selectedMCReco->FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_ByTgt)[USTgt].GetRecoValue(univ)), weight); //"Fake data" for closure
		(*fVtxTest_US_ByTgt)[USTgt].selectedSignalReco->FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_ByTgt)[USTgt].GetRecoValue(univ)), weight);
		if ((*fVtxTest_US_ByTgt)[USTgt].IsBroken()){
		  (*(*fVtxTest_US_ByTgt)[USTgt].m_SigIntTypeHists)[intType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_ByTgt)[USTgt].GetRecoValue(univ)), weight);
		  (*(*fVtxTest_US_ByTgt)[USTgt].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_ByTgt)[USTgt].GetRecoValue(univ)), weight);
		  //(*(*fVtxTest_US_ByTgt)[USTgt].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_ByTgt)[USTgt].GetRecoValue(univ)), weight);
		}

		(*fVtxTest2D_US_ByTgt)[USTgt].selectedMCReco->FillUniverse(&univ, (*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueY(univ)), weight); //"Fake data" for closure
		(*fVtxTest2D_US_ByTgt)[USTgt].selectedSignalReco->FillUniverse(&univ, (*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueY(univ)), weight);
		(*(*fVtxTest2D_US_ByTgt)[USTgt].m_SigIntTypeHists)[intType].FillUniverse(&univ, (*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueY(univ)), weight);
		(*(*fVtxTest2D_US_ByTgt)[USTgt].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueY(univ)), weight);
		//(*(*fVtxTest2D_US_ByTgt)[USTgt].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueY(univ)), weight);
	      }

	      if(DSTgt > 0){
		(*fVtxTest_DS_ByTgt)[DSTgt].selectedMCReco->FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_ByTgt)[DSTgt].GetRecoValue(univ)), weight); //"Fake data" for closure
		(*fVtxTest_DS_ByTgt)[DSTgt].selectedSignalReco->FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_ByTgt)[DSTgt].GetRecoValue(univ)), weight);
		if ((*fVtxTest_DS_ByTgt)[DSTgt].IsBroken()){
		  (*(*fVtxTest_DS_ByTgt)[DSTgt].m_SigIntTypeHists)[intType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_ByTgt)[DSTgt].GetRecoValue(univ)), weight);
		  (*(*fVtxTest_DS_ByTgt)[DSTgt].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_ByTgt)[DSTgt].GetRecoValue(univ)), weight);
		  //(*(*fVtxTest_DS_ByTgt)[DSTgt].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_ByTgt)[DSTgt].GetRecoValue(univ)), weight);
		}

		(*fVtxTest2D_DS_ByTgt)[DSTgt].selectedMCReco->FillUniverse(&univ, (*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueY(univ)), weight); //"Fake data" for closure
		(*fVtxTest2D_DS_ByTgt)[DSTgt].selectedSignalReco->FillUniverse(&univ, (*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueY(univ)), weight);
		(*(*fVtxTest2D_DS_ByTgt)[DSTgt].m_SigIntTypeHists)[intType].FillUniverse(&univ, (*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueY(univ)), weight);
		(*(*fVtxTest2D_DS_ByTgt)[DSTgt].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueY(univ)), weight);
		//(*(*fVtxTest2D_DS_ByTgt)[DSTgt].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueY(univ)), weight);
	      }
	    }
	  }
	  else if (tgtCode == -1){
	    if(USTgt > 0){
	      (*fVtxTest_US_Post_ByTgt)[USTgt].selectedMCReco->FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_Post_ByTgt)[USTgt].GetRecoValue(univ)), weight); //"Fake data" for closure
	      (*fVtxTest_US_Post_ByTgt)[USTgt].selectedSignalReco->FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_Post_ByTgt)[USTgt].GetRecoValue(univ)), weight);
	      if ((*fVtxTest_US_Post_ByTgt)[USTgt].IsBroken()){
		(*(*fVtxTest_US_Post_ByTgt)[USTgt].m_SigIntTypeHists)[intType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_Post_ByTgt)[USTgt].GetRecoValue(univ)), weight);
		(*(*fVtxTest_US_Post_ByTgt)[USTgt].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_Post_ByTgt)[USTgt].GetRecoValue(univ)), weight);
		//(*(*fVtxTest_US_Post_ByTgt)[USTgt].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_Post_ByTgt)[USTgt].GetRecoValue(univ)), weight);
	      }

	      (*fVtxTest2D_US_Post_ByTgt)[USTgt].selectedMCReco->FillUniverse(&univ, (*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueY(univ)), weight); //"Fake data" for closure
	      (*fVtxTest2D_US_Post_ByTgt)[USTgt].selectedSignalReco->FillUniverse(&univ, (*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueY(univ)), weight);
	      (*(*fVtxTest2D_US_Post_ByTgt)[USTgt].m_SigIntTypeHists)[intType].FillUniverse(&univ, (*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueY(univ)), weight);
	      (*(*fVtxTest2D_US_Post_ByTgt)[USTgt].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueY(univ)), weight);
	      //(*(*fVtxTest2D_US_Post_ByTgt)[USTgt].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueY(univ)), weight);

	      /* Debugging removed
		 std::cout << "IsSignal" << std::endl;
		 std::cout << "Material: " << tgtType << std::endl;
		 std::cout << "USTgt: " << USTgt << std::endl;
		 std::cout << "Vertex Z: " << vtx_z << std::endl;
		 std::cout << "Variable Vtx Z: " << (*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ) << std::endl;
		 std::cout << "NPlanes US for Target: " << util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)) << std::endl;
		 std::cout << "" << std::endl;
	      */
	      if (fVtx_US_ByTgt.size() > 0){
		(*(fVtx_US_ByTgt[USTgt]))[tgtType].selectedMCReco->FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), weight); //"Fake data" for closure
		(*(fVtx_US_ByTgt[USTgt]))[tgtType].selectedSignalReco->FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), weight);
		if ((*(fVtx_US_ByTgt[USTgt]))[tgtType].IsBroken()){
		  (*(*(fVtx_US_ByTgt[USTgt]))[tgtType].m_SigIntTypeHists)[intType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), weight);
		  (*(*(fVtx_US_ByTgt[USTgt]))[tgtType].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), weight);
		  //(*(*(fVtx_US_ByTgt[USTgt]))[tgtType].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), weight);
		}	
	      }
	      if (fSplitRecoil && fVtxBinned_US_ByTgt.size() > 0){
		(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].selectedMCReco->FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].GetRecoValue(univ)), weight); //"Fake data" for closure
		(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].selectedSignalReco->FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].GetRecoValue(univ)), weight);
		if ((*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].IsBroken()){
		  (*(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].m_SigIntTypeHists)[intType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].GetRecoValue(univ)), weight);
		  (*(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].GetRecoValue(univ)), weight);
		  //(*(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].GetRecoValue(univ)), weight);
		}
	      }
	    }

	    if(DSTgt > 0){
	      (*fVtxTest_DS_Post_ByTgt)[DSTgt].selectedMCReco->FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_Post_ByTgt)[DSTgt].GetRecoValue(univ)), weight); //"Fake data" for closure
	      (*fVtxTest_DS_Post_ByTgt)[DSTgt].selectedSignalReco->FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_Post_ByTgt)[DSTgt].GetRecoValue(univ)), weight);
	      if ((*fVtxTest_DS_Post_ByTgt)[DSTgt].IsBroken()){
		(*(*fVtxTest_DS_Post_ByTgt)[DSTgt].m_SigIntTypeHists)[intType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_Post_ByTgt)[DSTgt].GetRecoValue(univ)), weight);
		(*(*fVtxTest_DS_Post_ByTgt)[DSTgt].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_Post_ByTgt)[DSTgt].GetRecoValue(univ)), weight);
		//(*(*fVtxTest_DS_Post_ByTgt)[DSTgt].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_Post_ByTgt)[DSTgt].GetRecoValue(univ)), weight);
	      }

	      (*fVtxTest2D_DS_Post_ByTgt)[DSTgt].selectedMCReco->FillUniverse(&univ, (*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueY(univ)), weight); //"Fake data" for closure
	      (*fVtxTest2D_DS_Post_ByTgt)[DSTgt].selectedSignalReco->FillUniverse(&univ, (*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueY(univ)), weight);
	      (*(*fVtxTest2D_DS_Post_ByTgt)[DSTgt].m_SigIntTypeHists)[intType].FillUniverse(&univ, (*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueY(univ)), weight);
	      (*(*fVtxTest2D_DS_Post_ByTgt)[DSTgt].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueY(univ)), weight);
	      //(*(*fVtxTest2D_DS_Post_ByTgt)[DSTgt].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueY(univ)), weight);

	      /*Debugging removed
		std::cout << "IsSignal" << std::endl;
		std::cout << "Material: " << tgtType << std::endl;
		std::cout << "DSTgt: " << DSTgt << std::endl;
		std::cout << "Vertex Z: " << vtx_z << std::endl;
		std::cout << "Variable Vtx Z: " << (*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ) << std::endl;
		std::cout << "NPlanes DS for Target: " << util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)) << std::endl;
		std::cout << "" << std::endl;
	      */
	      if (fVtx_DS_ByTgt.size() > 0){
		(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].selectedMCReco->FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), weight); //"Fake data" for closure
		(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].selectedSignalReco->FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), weight);
		if ((*(fVtx_DS_ByTgt[DSTgt]))[tgtType].IsBroken()){
		  (*(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].m_SigIntTypeHists)[intType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), weight);
		  (*(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), weight);
		  //(*(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), weight);
		}	    
	      }
	      if (fSplitRecoil && fVtxBinned_DS_ByTgt.size() > 0){
		(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].selectedMCReco->FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].GetRecoValue(univ)), weight); //"Fake data" for closure
		(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].selectedSignalReco->FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].GetRecoValue(univ)), weight);
		if ((*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].IsBroken()){
		  (*(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].m_SigIntTypeHists)[intType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].GetRecoValue(univ)), weight);
		  (*(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].GetRecoValue(univ)), weight);
		  //(*(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].GetRecoValue(univ)), weight);
		}
	      }
	    }
	  }
	}

	else{
	  int bkgd_ID = -1;
	  if (isTgts){
	    bkgd_ID = 44;
	    if (tgtType == 200 || tgtType == 300){
	      if (evt.IsFSSignal()){
		if(tgtCode != -1){
		  if (tgtType == 200){
		    if (isSideband){
		      for (auto& var: fVars_US_ByTgt){
			(*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
			(*var)[tgtCode].selectedSignalReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			if ((*var)[tgtCode].IsBroken()){
			  (*(*var)[tgtCode].m_SigIntTypeHists)[intCode].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			  (*(*var)[tgtCode].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			  //(*(*var)[tgtCode].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			}
		      }
		      for (auto& var: fVars2D_US_ByTgt){
			(*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight); //"Fake data" for closure
			(*var)[tgtCode].selectedSignalReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			(*(*var)[tgtCode].m_SigIntTypeHists)[intType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			(*(*var)[tgtCode].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			//(*(*var)[tgtCode].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
		      }
		      if (fSplitRecoil && fRecoilBinned_US_ByTgt.size() > 0){
			(*fRecoilBinned_US_ByTgt[iBin])[tgtCode].selectedMCReco->FillUniverse(&univ, (*fRecoilBinned_US_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
			(*fRecoilBinned_US_ByTgt[iBin])[tgtCode].selectedSignalReco->FillUniverse(&univ, (*fRecoilBinned_US_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
			if ((*fRecoilBinned_US_ByTgt[iBin])[tgtCode].IsBroken()){
			  (*(*fRecoilBinned_US_ByTgt[iBin])[tgtCode].m_SigIntTypeHists)[intType].FillUniverse(&univ, (*fRecoilBinned_US_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
			  (*(*fRecoilBinned_US_ByTgt[iBin])[tgtCode].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*fRecoilBinned_US_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
			  //(*(*fRecoilBinned_US_ByTgt[iBin])[tgtCode].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*fRecoilBinned_US_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
			}
		      }
		    }
		    else if (SBStat.all()){
		      for (auto& var: fVars_US_Post_ByTgt){
			(*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
			(*var)[tgtCode].selectedSignalReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			if ((*var)[tgtCode].IsBroken()){
			  (*(*var)[tgtCode].m_SigIntTypeHists)[intCode].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			  (*(*var)[tgtCode].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			  //(*(*var)[tgtCode].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			}
		      }
		      for (auto& var: fVars2D_US_Post_ByTgt){
			(*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight); //"Fake data" for closure
			(*var)[tgtCode].selectedSignalReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			(*(*var)[tgtCode].m_SigIntTypeHists)[intType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			(*(*var)[tgtCode].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			//(*(*var)[tgtCode].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
		      }
		    }
		  }
		  else{
		    if (isSideband){
		      for (auto& var: fVars_DS_ByTgt){
			(*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
			(*var)[tgtCode].selectedSignalReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			if ((*var)[tgtCode].IsBroken()){
			  (*(*var)[tgtCode].m_SigIntTypeHists)[intCode].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			  (*(*var)[tgtCode].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			  //(*(*var)[tgtCode].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			}
		      }
		      for (auto& var: fVars2D_DS_ByTgt){
			(*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight); //"Fake data" for closure
			(*var)[tgtCode].selectedSignalReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			(*(*var)[tgtCode].m_SigIntTypeHists)[intType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			(*(*var)[tgtCode].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			//(*(*var)[tgtCode].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
		      }
		      if (fSplitRecoil && fRecoilBinned_DS_ByTgt.size() > 0){
			(*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].selectedMCReco->FillUniverse(&univ, (*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
			(*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].selectedSignalReco->FillUniverse(&univ, (*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
			if ((*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].IsBroken()){
			  (*(*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].m_SigIntTypeHists)[intType].FillUniverse(&univ, (*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
			  (*(*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
			  //(*(*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
			}
		      }
		    }
		    else if (SBStat.all()){
		      for (auto& var: fVars_DS_Post_ByTgt){
			(*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
			(*var)[tgtCode].selectedSignalReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			if ((*var)[tgtCode].IsBroken()){
			  (*(*var)[tgtCode].m_SigIntTypeHists)[intCode].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			  (*(*var)[tgtCode].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			  //(*(*var)[tgtCode].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			}
		      }
		      for (auto& var: fVars2D_DS_Post_ByTgt){
			(*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight); //"Fake data" for closure
			(*var)[tgtCode].selectedSignalReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			(*(*var)[tgtCode].m_SigIntTypeHists)[intType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			(*(*var)[tgtCode].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			//(*(*var)[tgtCode].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
		      }
		    }
		  }
		}
	      }
	      else{
		bkgd_ID = util::GetBackgroundID(univ);
		if(tgtCode != -1){
		  if (tgtType == 200){
		    if (isSideband){
		      for(auto& var: fVars_US_ByTgt){
			(*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
			(*(*var)[tgtCode].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			(*(*var)[tgtCode].m_BkgIntTypeHists)[intCode].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			if ((*var)[tgtCode].IsBroken()){
			  (*(*var)[tgtCode].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			  //(*(*var)[tgtCode].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			}
		      }
		      for(auto& var: fVars2D_US_ByTgt){
			(*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight); //"Fake data" for closure
			(*(*var)[tgtCode].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			(*(*var)[tgtCode].m_BkgIntTypeHists)[intType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			(*(*var)[tgtCode].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			//(*(*var)[tgtCode].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
		      }
		      if (fSplitRecoil && fRecoilBinned_US_ByTgt.size() > 0){
			(*fRecoilBinned_US_ByTgt[iBin])[tgtCode].selectedMCReco->FillUniverse(&univ, (*fRecoilBinned_US_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
			(*(*fRecoilBinned_US_ByTgt[iBin])[tgtCode].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*fRecoilBinned_US_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
			(*(*fRecoilBinned_US_ByTgt[iBin])[tgtCode].m_BkgIntTypeHists)[intType].FillUniverse(&univ, (*fRecoilBinned_US_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
			if ((*fRecoilBinned_US_ByTgt[iBin])[tgtCode].IsBroken()){
			  (*(*fRecoilBinned_US_ByTgt[iBin])[tgtCode].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*fRecoilBinned_US_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
			  //(*(*fRecoilBinned_US_ByTgt[iBin])[tgtCode].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*fRecoilBinned_US_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
			}
		      }
		    }
		    else if (SBStat.all()){
		      for(auto& var: fVars_US_Post_ByTgt){
			(*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
			(*(*var)[tgtCode].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			(*(*var)[tgtCode].m_BkgIntTypeHists)[intCode].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			if ((*var)[tgtCode].IsBroken()){
			  (*(*var)[tgtCode].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			  //(*(*var)[tgtCode].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			}
		      }
		      for(auto& var: fVars2D_US_Post_ByTgt){
			(*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight); //"Fake data" for closure
			(*(*var)[tgtCode].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			(*(*var)[tgtCode].m_BkgIntTypeHists)[intType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			(*(*var)[tgtCode].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			//(*(*var)[tgtCode].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
		      }
		    }
		  }
		  else{
		    if (isSideband){
		      for(auto& var: fVars_DS_ByTgt){
			(*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
			(*(*var)[tgtCode].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			(*(*var)[tgtCode].m_BkgIntTypeHists)[intCode].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			if ((*var)[tgtCode].IsBroken()){
			  (*(*var)[tgtCode].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			  //(*(*var)[tgtCode].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			}
		      }
		      for(auto& var: fVars2D_DS_ByTgt){
			(*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight); //"Fake data" for closure
			(*(*var)[tgtCode].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			(*(*var)[tgtCode].m_BkgIntTypeHists)[intType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			(*(*var)[tgtCode].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			//(*(*var)[tgtCode].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
		      }
		      if (fSplitRecoil && fRecoilBinned_DS_ByTgt.size() > 0){
			(*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].selectedMCReco->FillUniverse(&univ, (*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
			(*(*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
			(*(*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].m_BkgIntTypeHists)[intType].FillUniverse(&univ, (*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
			if ((*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].IsBroken()){
			  (*(*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
			  //(*(*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*fRecoilBinned_DS_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
			}
		      }
		    }
		    else if (SBStat.all()){
		      for(auto& var: fVars_DS_Post_ByTgt){
			(*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
			(*(*var)[tgtCode].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			(*(*var)[tgtCode].m_BkgIntTypeHists)[intCode].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			if ((*var)[tgtCode].IsBroken()){
			  (*(*var)[tgtCode].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			  //(*(*var)[tgtCode].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
			}
		      }
		      for(auto& var: fVars2D_DS_Post_ByTgt){
			(*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight); //"Fake data" for closure
			(*(*var)[tgtCode].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			(*(*var)[tgtCode].m_BkgIntTypeHists)[intType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			(*(*var)[tgtCode].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
			//(*(*var)[tgtCode].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
		      }
		    }
		  }
		}
	      }
	      bkgd_ID = intType;
	    }
	    else if (util::CorrectTargetMaterial(tgtCode,trueTgtCode) || (tgtCode == -1 && tgtCode == trueTgtCode)) bkgd_ID = util::GetBackgroundID(univ);//Modify this line so that tgtCode == -1 has a wrong nucleus component.
	    if (bkgd_ID == 44) intType = bkgd_ID;
	  }
	  else bkgd_ID = util::GetBackgroundID(univ);

	  if (fSplitRecoil){
	    if (fRecoilBinned.size() > 0){
	      fRecoilBinned[iBin]->selectedMCReco->FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight); //"Fake data" for closure
	      (*fRecoilBinned[iBin]->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	      (*fRecoilBinned[iBin]->m_BkgIntTypeHists)[intType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	      if (fRecoilBinned[iBin]->IsBroken()){
		(*fRecoilBinned[iBin]->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
		//(*fRecoilBinned[iBin]->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	      }
	    }

	    if (tgtCode != -1 && fRecoilBinned_ByTgt.size() > 0){
	      (*fRecoilBinned_ByTgt[iBin])[tgtCode].selectedMCReco->FillUniverse(&univ, (*fRecoilBinned_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
	      (*(*fRecoilBinned_ByTgt[iBin])[tgtCode].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*fRecoilBinned_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
	      (*(*fRecoilBinned_ByTgt[iBin])[tgtCode].m_BkgIntTypeHists)[intType].FillUniverse(&univ, (*fRecoilBinned_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
	      if ((*fRecoilBinned_ByTgt[iBin])[tgtCode].IsBroken()){
		(*(*fRecoilBinned_ByTgt[iBin])[tgtCode].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*fRecoilBinned_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
		//(*(*fRecoilBinned_ByTgt[iBin])[tgtCode].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*fRecoilBinned_ByTgt[iBin])[tgtCode].GetRecoValue(univ), weight);
	      }
	    }
	  }

	  if (isSideband){
	    for(auto& var: fVars){
	      if (var->IsFill()){
		var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(univ), weight); //"Fake data" for closure
		(*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValue(univ), weight);
		(*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
		if (var->IsBroken()){
		  (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
		  //(*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
		}
	      }
	    }
	  
	    for(auto& var: fVars2D){
	      var->selectedMCReco->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight); //"Fake data" for closure
	      (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	      (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	      (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	      //(*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    }

	    //Only fill for non-interstitial plastics
	    if(tgtCode != -1){
	      for(auto& var: fVars_ByTgt){
		(*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
		(*(*var)[tgtCode].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		(*(*var)[tgtCode].m_BkgIntTypeHists)[intType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		if ((*var)[tgtCode].IsBroken()){
		  (*(*var)[tgtCode].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		  //(*(*var)[tgtCode].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		}
	      }

	      for(auto& var: fVars2D_ByTgt){
		(*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight); //"Fake data" for closure
		(*(*var)[tgtCode].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
		(*(*var)[tgtCode].m_BkgIntTypeHists)[intType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
		(*(*var)[tgtCode].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
		//(*(*var)[tgtCode].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), weight);
	      }
	    }
	    else if (tgtCode == -1){
	      if (USTgt > 0){
		(*fVtxTest_US_ByTgt)[USTgt].selectedMCReco->FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_ByTgt)[USTgt].GetRecoValue(univ)), weight); //"Fake data" for closure
		(*(*fVtxTest_US_ByTgt)[USTgt].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_ByTgt)[USTgt].GetRecoValue(univ)), weight);
		(*(*fVtxTest_US_ByTgt)[USTgt].m_BkgIntTypeHists)[intType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_ByTgt)[USTgt].GetRecoValue(univ)), weight);
		if ((*fVtxTest_US_ByTgt)[USTgt].IsBroken()){
		  (*(*fVtxTest_US_ByTgt)[USTgt].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_ByTgt)[USTgt].GetRecoValue(univ)), weight);
		  //(*(*fVtxTest_US_ByTgt)[USTgt].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_ByTgt)[USTgt].GetRecoValue(univ)), weight);
		}

		(*fVtxTest2D_US_ByTgt)[USTgt].selectedMCReco->FillUniverse(&univ, (*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueY(univ)), weight); //"Fake data" for closure
		(*(*fVtxTest2D_US_ByTgt)[USTgt].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueY(univ)), weight);
		(*(*fVtxTest2D_US_ByTgt)[USTgt].m_BkgIntTypeHists)[intType].FillUniverse(&univ, (*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueY(univ)), weight);
		(*(*fVtxTest2D_US_ByTgt)[USTgt].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueY(univ)), weight);
		//(*(*fVtxTest2D_US_ByTgt)[USTgt].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueY(univ)), weight);
	      }

	      if (DSTgt > 0){
		(*fVtxTest_DS_ByTgt)[DSTgt].selectedMCReco->FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_ByTgt)[DSTgt].GetRecoValue(univ)), weight); //"Fake data" for closure
		(*(*fVtxTest_DS_ByTgt)[DSTgt].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_ByTgt)[DSTgt].GetRecoValue(univ)), weight);
		(*(*fVtxTest_DS_ByTgt)[DSTgt].m_BkgIntTypeHists)[intType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_ByTgt)[DSTgt].GetRecoValue(univ)), weight);
		if ((*fVtxTest_DS_ByTgt)[DSTgt].IsBroken()){
		  (*(*fVtxTest_DS_ByTgt)[DSTgt].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_ByTgt)[DSTgt].GetRecoValue(univ)), weight);
		  //(*(*fVtxTest_DS_ByTgt)[DSTgt].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_ByTgt)[DSTgt].GetRecoValue(univ)), weight);
		}

		(*fVtxTest2D_DS_ByTgt)[DSTgt].selectedMCReco->FillUniverse(&univ, (*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueY(univ)), weight); //"Fake data" for closure
		(*(*fVtxTest2D_DS_ByTgt)[DSTgt].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueY(univ)), weight);
		(*(*fVtxTest2D_DS_ByTgt)[DSTgt].m_BkgIntTypeHists)[intType].FillUniverse(&univ, (*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueY(univ)), weight);
		(*(*fVtxTest2D_DS_ByTgt)[DSTgt].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueY(univ)), weight);
		//(*(*fVtxTest2D_DS_ByTgt)[DSTgt].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueY(univ)), weight);
	      }
	    }
	  }
	  else if (tgtCode == -1){
	    if(USTgt > 0){
	      (*fVtxTest_US_Post_ByTgt)[USTgt].selectedMCReco->FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_Post_ByTgt)[USTgt].GetRecoValue(univ)), weight); //"Fake data" for closure
	      (*(*fVtxTest_US_Post_ByTgt)[USTgt].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_Post_ByTgt)[USTgt].GetRecoValue(univ)), weight);
	      (*(*fVtxTest_US_Post_ByTgt)[USTgt].m_BkgIntTypeHists)[intType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_Post_ByTgt)[USTgt].GetRecoValue(univ)), weight);
	      if ((*fVtxTest_US_Post_ByTgt)[USTgt].IsBroken()){
		(*(*fVtxTest_US_Post_ByTgt)[USTgt].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_Post_ByTgt)[USTgt].GetRecoValue(univ)), weight);
		//(*(*fVtxTest_US_Post_ByTgt)[USTgt].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_Post_ByTgt)[USTgt].GetRecoValue(univ)), weight);
	      }

	      (*fVtxTest2D_US_Post_ByTgt)[USTgt].selectedMCReco->FillUniverse(&univ, (*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueY(univ)), weight); //"Fake data" for closure
	      (*(*fVtxTest2D_US_Post_ByTgt)[USTgt].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueY(univ)), weight);
	      (*(*fVtxTest2D_US_Post_ByTgt)[USTgt].m_BkgIntTypeHists)[intType].FillUniverse(&univ, (*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueY(univ)), weight);
	      (*(*fVtxTest2D_US_Post_ByTgt)[USTgt].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueY(univ)), weight);
	      //(*(*fVtxTest2D_US_Post_ByTgt)[USTgt].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueY(univ)), weight);

	      /* Debugging removed
		 std::cout << "IsBKG" << std::endl;
		 std::cout << "Material: " << tgtType << std::endl;
		 std::cout << "USTgt: " << USTgt << std::endl;
		 std::cout << "Vertex Z: " << vtx_z << std::endl;
		 std::cout << "Variable Vtx Z: " << (*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ) << std::endl;
		 std::cout << "NPlanes US for Target: " << util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)) << std::endl;
		 std::cout << "" << std::endl;
	      */
	      if (fVtx_US_ByTgt.size() > 0){
		(*(fVtx_US_ByTgt[USTgt]))[tgtType].selectedMCReco->FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), weight); //"Fake data" for closure
		(*(*(fVtx_US_ByTgt[USTgt]))[tgtType].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), weight);
		(*(*(fVtx_US_ByTgt[USTgt]))[tgtType].m_BkgIntTypeHists)[intType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), weight);
		if ((*(fVtx_US_ByTgt[USTgt]))[tgtType].IsBroken()){
		  (*(*(fVtx_US_ByTgt[USTgt]))[tgtType].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), weight);
		  //(*(*(fVtx_US_ByTgt[USTgt]))[tgtType].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), weight);
		}
	      }
	      if (fSplitRecoil && fVtxBinned_US_ByTgt.size() > 0){
		(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].selectedMCReco->FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].GetRecoValue(univ)), weight); //"Fake data" for closure
		(*(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].GetRecoValue(univ)), weight);
		(*(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].m_BkgIntTypeHists)[intType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].GetRecoValue(univ)), weight);
		if ((*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].IsBroken()){
		  (*(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].GetRecoValue(univ)), weight);
		  //(*(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].GetRecoValue(univ)), weight);
		}
	      }
	    }
	    if(DSTgt > 0){
	      (*fVtxTest_DS_Post_ByTgt)[DSTgt].selectedMCReco->FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_Post_ByTgt)[DSTgt].GetRecoValue(univ)), weight); //"Fake data" for closure
	      (*(*fVtxTest_DS_Post_ByTgt)[DSTgt].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_Post_ByTgt)[DSTgt].GetRecoValue(univ)), weight);
	      (*(*fVtxTest_DS_Post_ByTgt)[DSTgt].m_BkgIntTypeHists)[intType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_Post_ByTgt)[DSTgt].GetRecoValue(univ)), weight);
	      if ((*fVtxTest_DS_Post_ByTgt)[DSTgt].IsBroken()){
		(*(*fVtxTest_DS_Post_ByTgt)[DSTgt].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_Post_ByTgt)[DSTgt].GetRecoValue(univ)), weight);
		//(*(*fVtxTest_DS_Post_ByTgt)[DSTgt].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_Post_ByTgt)[DSTgt].GetRecoValue(univ)), weight);
	      }

	      (*fVtxTest2D_DS_Post_ByTgt)[DSTgt].selectedMCReco->FillUniverse(&univ, (*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueY(univ)), weight); //"Fake data" for closure
	      (*(*fVtxTest2D_DS_Post_ByTgt)[DSTgt].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueY(univ)), weight);
	      (*(*fVtxTest2D_DS_Post_ByTgt)[DSTgt].m_BkgIntTypeHists)[intType].FillUniverse(&univ, (*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueY(univ)), weight);
	      (*(*fVtxTest2D_DS_Post_ByTgt)[DSTgt].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueY(univ)), weight);
	      //(*(*fVtxTest2D_DS_Post_ByTgt)[DSTgt].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueY(univ)), weight);

	      /*Debugging removed
		std::cout << "IsBKG" << std::endl;
		std::cout << "Material: " << tgtType << std::endl;
		std::cout << "DSTgt: " << DSTgt << std::endl;
		std::cout << "Vertex Z: " << vtx_z << std::endl;
		std::cout << "Variable Vtx Z: " << (*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ) << std::endl;
		std::cout << "NPlanes DS for Target: " << util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)) << std::endl;
		std::cout << "" << std::endl;
	      */
	      if (fVtx_DS_ByTgt.size() > 0){
		(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].selectedMCReco->FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), weight); //"Fake data" for closure
		(*(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), weight);
		(*(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].m_BkgIntTypeHists)[intType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), weight);
		if ((*(fVtx_DS_ByTgt[DSTgt]))[tgtType].IsBroken()){
		  (*(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), weight);
		  //(*(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), weight);
		}
	      }
	      if (fSplitRecoil && fVtxBinned_DS_ByTgt.size() > 0){
		(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].selectedMCReco->FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].GetRecoValue(univ)), weight); //"Fake data" for closure
		(*(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].GetRecoValue(univ)), weight);
		(*(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].m_BkgIntTypeHists)[intType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].GetRecoValue(univ)), weight);
		if ((*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].IsBroken()){
		  (*(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].GetRecoValue(univ)), weight);
		  //(*(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].GetRecoValue(univ)), weight);
		}
	      }
	    }
	  }
	}
      }
      
      else{

	if (fSplitRecoil){
	  if (fRecoilBinned.size() > 0){
	    fRecoilBinned[iBin]->dataHist->FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), 1);
	  }
	  if (tgtCode != -1 && fRecoilBinned_ByTgt.size() > 0) (*fRecoilBinned_ByTgt[iBin])[tgtCode].dataHist->FillUniverse(&univ, (*fRecoilBinned_ByTgt[iBin])[tgtCode].GetRecoValue(univ), 1);
	}

	if (isSideband){
	  for (auto& var : fVars){
	    if (var->IsFill()) var->dataHist->FillUniverse(&univ, var->GetRecoValue(univ), 1);
	  }

	  for (auto& var : fVars2D){
	    var->dataHist->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), 1);
	  }
	  
	  //Only fill for non-interstitial plastics
	  if(tgtCode != -1){
	    for (auto& var : fVars_ByTgt){
	      (*var)[tgtCode].dataHist->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), 1);
	    }

	    for (auto& var : fVars2D_ByTgt){
	      (*var)[tgtCode].dataHist->FillUniverse(&univ, (*var)[tgtCode].GetRecoValueX(univ), (*var)[tgtCode].GetRecoValueY(univ), 1);
	    }
	  }
	  
	  else if (tgtCode == -1){
	    if(USTgt > 0){
	      (*fVtxTest_US_ByTgt)[USTgt].dataHist->FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_ByTgt)[USTgt].GetRecoValue(univ)), 1);
	      (*fVtxTest2D_US_ByTgt)[USTgt].dataHist->FillUniverse(&univ, (*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_ByTgt)[USTgt].GetRecoValueY(univ)), 1); //"Fake data" for closure
	    }
	    if(DSTgt > 0){
	      (*fVtxTest_DS_ByTgt)[DSTgt].dataHist->FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_ByTgt)[DSTgt].GetRecoValue(univ)), 1);
	      (*fVtxTest2D_DS_ByTgt)[DSTgt].dataHist->FillUniverse(&univ, (*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_ByTgt)[DSTgt].GetRecoValueY(univ)), 1); //"Fake data" for closure
	    }
	  }
	}

	// REMOVE FILLING THE PLASTIC DATA. It's useless... except it's what I need to fit... wheeeeeeeeee.... so I un-removed... but this set me back a minute... oopsie...
	else if (tgtCode == -1){
	  if(USTgt > 0){
	    (*fVtxTest_US_Post_ByTgt)[USTgt].dataHist->FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest_US_Post_ByTgt)[USTgt].GetRecoValue(univ)), 1);

	    (*fVtxTest2D_US_Post_ByTgt)[USTgt].dataHist->FillUniverse(&univ, (*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueX(univ), util::GetNPlanesUSOfTarget(USTgt,(*fVtxTest2D_US_Post_ByTgt)[USTgt].GetRecoValueY(univ)), 1); //"Fake data" for closure

	    /* removed debugging
	       std::cout << "IsData" << std::endl;
	       std::cout << "Material: " << tgtType << std::endl;
	       std::cout << "USTgt: " << USTgt << std::endl;
	       std::cout << "Vertex Z: " << vtx_z << std::endl;
	       std::cout << "Variable Vtx Z: " << (*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ) << std::endl;
	       std::cout << "NPlanes US for Target: " << util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)) << std::endl;
	       std::cout << "" << std::endl;
	    */
	    if (fVtx_US_ByTgt.size() > 0){
	      (*(fVtx_US_ByTgt[USTgt]))[tgtType].dataHist->FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), 1);
	    }
	    if (fSplitRecoil && fVtxBinned_US_ByTgt.size() > 0) (*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].dataHist->FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtxBinned_US_ByTgt[iBin][USTgt]))[tgtType].GetRecoValue(univ)), 1);
	  }
	  if(DSTgt > 0){
	    (*fVtxTest_DS_Post_ByTgt)[DSTgt].dataHist->FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest_DS_Post_ByTgt)[DSTgt].GetRecoValue(univ)), 1);

	    (*fVtxTest2D_DS_Post_ByTgt)[DSTgt].dataHist->FillUniverse(&univ, (*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueX(univ), util::GetNPlanesDSOfTarget(DSTgt,(*fVtxTest2D_DS_Post_ByTgt)[DSTgt].GetRecoValueY(univ)), 1); //"Fake data" for closure

	    /* removed debugging
	       std::cout << "IsData" << std::endl;
	       std::cout << "Material: " << tgtType << std::endl;
	       std::cout << "DSTgt: " << DSTgt << std::endl;
	       std::cout << "Vertex Z: " << vtx_z << std::endl;
	       std::cout << "Variable Vtx Z: " << (*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ) << std::endl;
	       std::cout << "NPlanes DS for Target: " << util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)) << std::endl;
	       std::cout << "" << std::endl;
	    */
	    if (fVtx_DS_ByTgt.size() > 0){
	      (*(fVtx_DS_ByTgt[DSTgt]))[tgtType].dataHist->FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), 1);
	    }
	    if (fSplitRecoil && fVtxBinned_DS_ByTgt.size() > 0) (*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].dataHist->FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtxBinned_DS_ByTgt[iBin][DSTgt]))[tgtType].GetRecoValue(univ)), 1);	  
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

//Actually requires the sideband to be satisfied. Written as though the recoil cut is the 0th sideband cut
class RecoilSB: public Study
{
  private:
    std::vector<Variable*> fVars;
    std::vector<Variable2D*> fVars2D;
    std::map<int,Variable*> fRecoilBinned;
    bool fSplitRecoil;

  public:
    RecoilSB(std::vector<Variable*> vars,
		     std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
		     std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands,
	     std::vector<CVUniverse*>& data_error_bands, bool splitRecoil): Study(), fSplitRecoil(splitRecoil)
    {

      if (fSplitRecoil){
	//Constructing recoil variable in bins of pt/pz. Just make a map of it all.
	std::vector<double> myRecoilBins;
	const double myRecoilBinWidth = 1.0/50.;
	for(int whichBin = 0; whichBin < 51; ++whichBin) myRecoilBins.push_back(myRecoilBinWidth * whichBin);

	int nBins = 14;//Hard-coded from bins from Amit's analysis... Not appropriate binning given my cuts, but going to get the structure going first.
	int lowBin = -1;//Default return value. hard-coded as well.

	for (int iBin=lowBin; iBin < nBins; ++iBin){
	  std::string binName = std::to_string(iBin);
	  if (iBin == lowBin) binName = "lost";
	  fRecoilBinned[iBin] = new Variable(false,("recoilE_RecoilSB_bin_"+binName).c_str(), "Recoil E [GeV]", myRecoilBins, &CVUniverse::GetDANRecoilEnergyGeV);
	  fRecoilBinned[iBin]->InitializeMCHists(mc_error_bands, truth_error_bands);
	  fRecoilBinned[iBin]->InitializeDATAHists(data_error_bands);
	}
      }

      for (auto& var : vars){
	fVars.push_back(new Variable(false,(var->GetName()+"_RecoilSB").c_str(), var->GetAxisLabel(), var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
      }

      fVars2D = {
	new Variable2D(false, *fVars[4],*fVars[3]),//recoil v. Q2     
	new Variable2D(false, *fVars[fVars.size()-1],*fVars[0]),//pT v. recoilQ2Bin     
      };

      for(auto& var: fVars) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fVars) var->InitializeDATAHists(data_error_bands);

      for(auto& var: fVars2D) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fVars2D) var->InitializeDATAHists(data_error_bands);
    }

    void SaveOrDraw(TDirectory& outDir){
      return;
    }

    void SaveOrDrawMC(TFile& outFile)
    {
      for (auto& recoil : fRecoilBinned) recoil.second->WriteMC(outFile);
      for (auto& var : fVars) var->WriteMC(outFile);
      for (auto& var : fVars2D) var->WriteMC(outFile);
    }

    void SaveOrDrawData(TFile& outFile)
    {
      for (auto& recoil : fRecoilBinned) recoil.second->WriteData(outFile);
      for (auto& var : fVars) var->WriteData(outFile);
    }
    
  private:
    //Overriding base class functions
    //Do nothing for now...  Good place for data comparisons in the future. 
    void fillSelected(const CVUniverse& univ, const NeutronEvent& evt, const double weight) {
      std::bitset<64> SBStat = evt.GetSideBandStat();
      int intType = evt.GetIntType();
      int tgtType = evt.GetTgtZ();
      int leadBlobType = evt.GetLeadingNeutCand().GetPDGBin();
      int iBin = evt.GetBinPTPZ();

      if (SBStat[0] == 1) return;
      
      if (evt.IsMC()){
	if (evt.IsSignal()){
	  if (fSplitRecoil){
	    fRecoilBinned[iBin]->selectedMCReco->FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight); //"Fake data" for closure
	    fRecoilBinned[iBin]->selectedSignalReco->FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    (*fRecoilBinned[iBin]->m_SigIntTypeHists)[intType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    (*fRecoilBinned[iBin]->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    //(*fRecoilBinned[iBin]->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	  }

	  for (auto& var: fVars){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(univ), weight); //"Fake data" for closure
	    var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    //(*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	  }

	  for (auto& var: fVars2D){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight); //"Fake data" for closure
	    var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    //(*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	  }
	}
	
	else{
	  int bkgd_ID = -1;	  
	  bkgd_ID = util::GetBackgroundID(univ);
	  
	  if (fSplitRecoil){
	    fRecoilBinned[iBin]->selectedMCReco->FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight); //"Fake data" for closure
	    (*fRecoilBinned[iBin]->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    (*fRecoilBinned[iBin]->m_BkgIntTypeHists)[intType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    (*fRecoilBinned[iBin]->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    //(*fRecoilBinned[iBin]->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	  }

	  for(auto& var: fVars){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(univ), weight); //"Fake data" for closure
	    (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    //(*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	  }

	  for(auto& var: fVars2D){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight); //"Fake data" for closure
	    (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    //(*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	  }
	}
      }
      
      else{
	if (fSplitRecoil) fRecoilBinned[iBin]->dataHist->FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), 1);

	for (auto& var : fVars){
	  var->dataHist->FillUniverse(&univ, var->GetRecoValue(univ), 1);
	}

	for (auto& var : fVars2D){
	  var->dataHist->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), 1);
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
