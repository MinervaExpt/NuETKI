//==============================================================================
//In this file several inclusive cuts are defined.
//
//Each cut is a class that inherits from PlotUtils::Cut.
//
//At minimum, cuts must have a name and they must override the checkCut
//function.
//
//You cut also takes two template parameters:
//* UNIVERSE template parameter: Cuts are built on Universe objects. Plug
//in your CVUniverse as a template parameter so the cut can access your
//branches and do so correctly within a systematic universe.
//
// * EVENT template parameter: sometimes cuts need to do more than just return
//a bool when you call the checkCut function. The event object is a
//user-specifiable object that can hold onto information that is learned within
//checkCut. E.g. A HasMichel cut determines which tracks are potential pion
//tracks. See CCPionCuts.h for a fleshed out example making use of EVENT.
//==============================================================================

//Original Author: Andrew Olivier aolivier@ur.rochester.edu

#include "PlotUtils/Cut.h"
#include <sstream>
#ifndef BEN_CCINCLUSIVECUTS_H
#define BEN_CCINCLUSIVECUTS_H

#ifdef __GCCXML__
#define override
namespace std {
  std::string to_string(double x) {
    std::stringstream ss;
    ss << x;
    return ss.str();
  }
}
#endif

namespace reco
{
  //============================================================================
  //Example 1: The simplest cut example. Just derive from Cut<> base class.
  //============================================================================
  //completely unused for my analysis (Carlos)
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class HasMINOSMatch: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
      // Constructor
      HasMINOSMatch(): PlotUtils::Cut<UNIVERSE, EVENT>("Has MINOS Match") {}

    private:
      // THE cut function
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
      {
        // Call a CVUniverse member function to make the cut
        return univ.IsMinosMatchMuon();
      }
  };

  //============================================================================
  //Example 2: Many cuts require variables to be above, below, or at some value.
  //To simply help avoid typos and reduce a little typing, use Minimum,
  //Maximum, and IsSame helper templates.
  //In this case the muon energy is required to be at least X, at minimum.
  //The specific value of X gets set when you instantiate this cut (in
  //getCCInclusiveCuts).
  //============================================================================
#ifndef __GCCXML__
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  using MuonEnergyMin = PlotUtils::Minimum<UNIVERSE, double, &UNIVERSE::GetEmu, EVENT>;
#endif
  //============================================================================
  //Example 3: The first, commented-out dead time cut inherits directly from Ben's
  //CVUniverse class. (In form, it resembles case 1, by the way.) While this is
  //valid, it means that no one else can use this cut (because they don't also
  //use my CVUniverse).
  //============================================================================
  /*class NoDeadtime: public Cut<CVUniverse, detail::empty>
  {
    public:
      NoDeadtime(const int nDead = 1): Cut("Dead Discriminators"), m_nDeadAllowed(nDead)
      {
      }

    private:
      const int m_nDeadAllowed; //Number of dead discriminators allowed

      bool checkCut(const CVUniverse& univ, const detail::empty&) const override
      {
        return univ.GetTDead() < m_nDeadAllowed;
      }
  };*/

  //============================================================================
  //Example 3 contd: It's better to make the cut a template so that any analyzer's
  //CVUniverse will work with it.
  //Furthermore, look how much space we save by using the Maximum helper
  //template.
  //============================================================================
#ifndef __GCCXML__
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  using NoDeadtime = PlotUtils::Maximum<UNIVERSE, int, &UNIVERSE::GetTDead, EVENT>;
#endif

  //Eavail Cut
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class Eavailable: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    // Constructor                                                                                                                                                  
    Eavailable(): PlotUtils::Cut<UNIVERSE, EVENT>("Available energy < 2 GeV") {}

    private:
    // THE cut function                                                                                                                                             
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      // Call a CVUniverse member function to make the cut                                                                                                   
      //      std::cout <<"blob_recoil_E_tracker: "
      //double Eavail = (univ.GetDouble("blob_recoil_E_tracker") + univ.GetDouble("blob_recoil_E_ecal")) * 1.17 - (0.008 * univ.GetDouble("prong_part_E") + 5);
      const double Eavail = univ.GetEavail();
      return Eavail < 2;
      //return true;
    }
  };

  //Modified Eavail Cut , which is Eavail minus sum of all tracked proton energies
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class ModifiedEavailable: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    // Constructor                                                                                                                                                  
    ModifiedEavailable(): PlotUtils::Cut<UNIVERSE, EVENT>("Modified Available energy < 2 GeV") {}

    private:
    // THE cut function                                                                                                                                             
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      // Call a CVUniverse member function to make the cut                                                                                                   
      //      std::cout <<"blob_recoil_E_tracker: "
      //double Eavail = (univ.GetDouble("blob_recoil_E_tracker") + univ.GetDouble("blob_recoil_E_ecal")) * 1.17 - (0.008 * univ.GetDouble("prong_part_E") + 5);
      const double ModEavail = univ.GetModifiedEavail();
      return ModEavail < 2;
      //return true;
    }
  };

  //LeptonPt , gonna leave this one for later since it's not a branch straight up, it's some calculated quantity. double check ryan's to see how he does this cut
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class LeptonPt: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    LeptonPt(): PlotUtils::Cut<UNIVERSE, EVENT>("Electron transverse momentum") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return (univ.GetLeptonPt() > 0.2 && univ.GetLeptonPt() < 1.6);
    }
  };

  //Electron energy, kevin suggested i hit a SUPER aggressive cut just to isolate my events. screw it lets do it
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class ElectronEnergy: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    ElectronEnergy(): PlotUtils::Cut<UNIVERSE, EVENT>("Electron energy > 2.5 GeV") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetElectronEnergy() > 2.5;
    }
  };

  //EleptonSin2Theta
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class EleptonSin2Theta: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    EleptonSin2Theta(): PlotUtils::Cut<UNIVERSE, EVENT>("E_lep * Theta_lep^2 > 0.003") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return (univ.GetETheta() > 0.003);
    }
  };

  //Proton starts near electron prong start point
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class ZDifference: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    ZDifference(): PlotUtils::Cut<UNIVERSE, EVENT>("ProtonStartZ - ShowerStartZ") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      double zdif = univ.GetProtonStartZ() - univ.GetShowerStartZ();
      //std::cout << "ProtonStart: " << univ.GetProtonStartZ() << std::endl;
      //std::cout << "ShowerStart: " << univ.GetShowerStartZ() << std::endl;
      return zdif >= -30 && zdif <= 50;
    }
  };

  //MasterAnaDev_proton_* branches filled, altho this one specifically checks MasterAnaDev_proton_endPointZ
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class ProtonInEvent: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    ProtonInEvent(): PlotUtils::Cut<UNIVERSE, EVENT>("Has Reco Primary Proton") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      //std::cout << "protonEndZ: " << univ.GetProtonEndZ() << std::endl;
      double protonEndZ = univ.GetProtonEndZ();
      //return protonEndZ < 8500 && protonEndZ > -10;

      //This is effectively just, do we have a primary proton in the MADtuple.
      return protonEndZ > -1;
    }
  };

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class ProtonMomentum: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    ProtonMomentum(): PlotUtils::Cut<UNIVERSE, EVENT>("Reco proton momentum b/w 0.45 and 1.2 GeV/c") {}
    
    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      double protonP = univ.GetProtonP();
      return (protonP>0.45 && protonP<1.2);
    }
  };

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class ProtonAngle: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    ProtonAngle(): PlotUtils::Cut<UNIVERSE, EVENT>("Reco proton angle < 70 deg") {}
    
    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      double protonTheta = univ.GetProtonTheta();
      return (protonTheta < 70);
    }
  };

  //HasTracks
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class HasTracks: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    HasTracks(): PlotUtils::Cut<UNIVERSE, EVENT>("Has Tracks") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return (univ.GetHasTracks()==1);
    }
  };

  //HasNoBackExitingTracks
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class HasNoBackExitingTracks: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    HasNoBackExitingTracks(): PlotUtils::Cut<UNIVERSE, EVENT>("No back exiting tracks") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return (univ.GetExitsBack() == 1);
    }
  };

  //VertexZ is in tracker
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class ZRange: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
      ZRange(const std::string& name, const double zMin, const double zMax): PlotUtils::Cut<UNIVERSE, EVENT>(name), fMin(zMin), fMax(zMax)
      {
      }

    private:
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
      {
        return univ.GetVertexZ() >= fMin && univ.GetVertexZ() <= fMax;
      }

      const double fMin;
      const double fMax;
  };

  //Apothem
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class Apothem: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    Apothem(const double apothem): PlotUtils::Cut<UNIVERSE, EVENT>(std::string("Apothem ") + std::to_string(apothem)), fApothem(apothem), fSlope(-1./sqrt(3.))//A regular hexagon has angles of 2*M_PI/3, so I can find this is 1/tan(M_PI/3.)
      {
      }

    private:
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
      {
        const ROOT::Math::XYZTVector vertex = univ.GetVertex();
        return (fabs(vertex.y()) < fSlope*fabs(vertex.x()) + 2.*fApothem/sqrt(3.))
               && (fabs(vertex.x()) < fApothem);
      }

      const double fApothem;
      const double fSlope; 
  };


  //EMLikeTrackScore
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class EMLikeTrackScore: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    EMLikeTrackScore(): PlotUtils::Cut<UNIVERSE, EVENT>("EM Shower candidate score") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      /*      bool passesCut = false;
      std::vector<double> prong_part_scores = univ.GetVecDouble("prong_part_score");
      for (int i = 0; i < prong_part_scores.size(); i++){
	if (prong_part_scores[i] > 0.7) passesCut = true;
      }
      */
      return (univ.GetHighestEMLikeShowerScore() > 0.7);
    }
  };

  //DSCalVisE
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class DSCalVisE: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    DSCalVisE(): PlotUtils::Cut<UNIVERSE, EVENT>("DSCalVisE") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return (univ.GetDSCalVisE() <= 0.2);
    }
  };

  //ODCalVisE
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class ODCalVisE: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    ODCalVisE(): PlotUtils::Cut<UNIVERSE, EVENT>("ODCalVisE") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return (univ.GetODCalVisE() <= 0.05);
    }
  };

  //Afterpulsing
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class Afterpulsing: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    Afterpulsing(): PlotUtils::Cut<UNIVERSE, EVENT>("Afterpulsing") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return (univ.GetFirstFireFraction() >= 0.25);
    }
  };

  //NonMIPClusterFraction
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class NonMIPClusterFraction: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    NonMIPClusterFraction(): PlotUtils::Cut<UNIVERSE, EVENT>("Prong Non MIP cluster fraction") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return (univ.GetNonMIPClusFrac() > 0.7); //testing 0.7, was 0.4
    }
  };

  //NoVertexMismatch
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class NoVertexMismatch: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    NoVertexMismatch(): PlotUtils::Cut<UNIVERSE, EVENT>("No Vertex Mismatch") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return (univ.GetHasNoVertexMismatch() == 1);
    }
  };

  //VertexTrackMultiplicity
  //can probably tune this one... ?
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class VertexTrackMultiplicity: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    VertexTrackMultiplicity(): PlotUtils::Cut<UNIVERSE, EVENT>("Vertex Track multiplicity") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      //return (univ.GetVertexTrackMultiplicity() > 2 && univ.GetVertexTrackMultiplicity() < 6);
      return (univ.GetVertexTrackMultiplicity() <= 2); //was at < 6, testing <= 2 but this removes multi tracked proton events so idk...

    }
  };

  //StartPointVertexMultiplicity
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class StartPointVertexMultiplicity: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    StartPointVertexMultiplicity(): PlotUtils::Cut<UNIVERSE, EVENT>("Start point vertex multiplicity") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return (univ.GetStartPointVertexMultiplicity() == 1);
    }
  };


  //TransverseGapScore
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class TransverseGapScore: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    TransverseGapScore(): PlotUtils::Cut<UNIVERSE, EVENT>("Transverse gap score > 15") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetTransverseGapScore() > 15;
    }
  };

  //MeanFrontdEdX
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class MeanFrontdEdX: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    MeanFrontdEdX(): PlotUtils::Cut<UNIVERSE, EVENT>("Mean front shower dE/dX < 2.4 MeV/cm") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetMeanFrontdEdx() < 2.4;
    }
  };

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class MeanFrontdEdXAbove: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    MeanFrontdEdXAbove(): PlotUtils::Cut<UNIVERSE, EVENT>("Mean front shower dE/dX >= 2.4 MeV/cm (sideband)") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetMeanFrontdEdx() >= 2.4;
    }
  };

  
  //Extra Energy Ratio, aka Psi
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class Psi: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    Psi(): PlotUtils::Cut<UNIVERSE, EVENT>("Psi (Extra Energy ratio) < 0.1") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetPsi() < 0.1;
    }
  };

  //Michel electron cut, very rudimentary rn. 
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class HasNoMichel: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    HasNoMichel(): PlotUtils::Cut<UNIVERSE, EVENT>("improved_nmichel = 0") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetImprovedNMichel()==0;
    }
  };

  //Michel electron cut flipped, for sideband
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class HasMichel: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    HasMichel(): PlotUtils::Cut<UNIVERSE, EVENT>("improved_nmichel > 0") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetImprovedNMichel()>0;
    }
  };


  //Elastically Scattered Contained (ESC) proton cut
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class ESC: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    ESC(): PlotUtils::Cut<UNIVERSE, EVENT>("ESC Proton Cut (chi2 < 10)") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetProtonESCNodeChi2() < 10; //please don't leave this hardcoded carlos
    }
  };


#ifndef __GCCXML__
  //============================================================================
  // This function instantiates each of the above cuts and adds them to a
  // container, over which we'll loop during our event selection to apply the
  // cuts.
  // 
  // The return type for this function is a `cuts_t<UNIVERSE, EVENT>`, which is
  // shorthand for std::vector<std::unique_ptr<PlotUtils::Cut<UNIVERSE, EVENT>>>;
  //============================================================================
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  PlotUtils::cuts_t<UNIVERSE, EVENT> GetCCInclusiveCuts()
  {
    PlotUtils::cuts_t<UNIVERSE, EVENT> inclusive_cuts;

    //inclusive_cuts.emplace_back(new HasMINOSMatch<UNIVERSE, EVENT>());
    //inclusive_cuts.emplace_back(new MuonEnergyMin<UNIVERSE, EVENT>(2e3, "Emu"));
    inclusive_cuts.emplace_back(new NoDeadtime<UNIVERSE, EVENT>(1, "Deadtime"));
    //inclusive_cuts.emplace_back(new IsNeutrino<UNIVERSE, EVENT>());

    return inclusive_cuts;
  }

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  PlotUtils::cuts_t<UNIVERSE, EVENT> GetCCInclusive2DCuts()
  { 
    PlotUtils::cuts_t<UNIVERSE, EVENT> inclusive_cuts;
    inclusive_cuts.emplace_back(new ZRange<UNIVERSE, EVENT>("Tracker", 5980, 8422));
    inclusive_cuts.emplace_back(new Apothem<UNIVERSE, EVENT>(850.));
    //inclusive_cuts.emplace_back(new MaxMuonAngle<UNIVERSE, EVENT>(20.));
    //inclusive_cuts.emplace_back(new HasMINOSMatch<UNIVERSE, EVENT>());
    inclusive_cuts.emplace_back(new NoDeadtime<UNIVERSE, EVENT>(1, "Deadtime"));
    //inclusive_cuts.emplace_back(new IsNeutrino<UNIVERSE, EVENT>());

    return inclusive_cuts;
  }
#endif
  //TODO: MnvGENIEv1, nu-e constraint, 50 flux universes
  //TODO: Binning: [ 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20, 40, 60]
  //               [ 0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1, 1.25, 1.5, 2.5, 4.5]

}
#ifdef __GCCXML__
#undef override
#endif

#endif //BEN_CCINCLUSIVECUTS_H
