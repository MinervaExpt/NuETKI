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

namespace utils {
  std::string to_string_trimmed(double x, int precision = 2) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(precision) << x;
    std::string result = ss.str();

    result.erase(result.find_last_not_of('0') + 1); //trims trailing zeroes

    if (!result.empty() && result.back() == '.') { //trims dot if that's the last char
        result.pop_back();
    }

    return result;
  }
}
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

  //Modified Eavail Cut , which is Eavail minus sum of all candidate proton energies
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class ModifiedEavailable: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    // Constructor                                                                                                                                                  
    ModifiedEavailable(const double max): PlotUtils::Cut<UNIVERSE, EVENT>("Non Proton Available energy < " + utils::to_string_trimmed(max) + " GeV"), fMax(max) {}

    private:
    // THE cut function                                                                                                                                             
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetModifiedEavail() < fMax;
    }
    const double fMax;
  };

  //ElectronPt , gonna leave this one for later since it's not a branch straight up, it's some calculated quantity. double check ryan's to see how he does this cut
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class ElectronPt: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    ElectronPt(): PlotUtils::Cut<UNIVERSE, EVENT>("Electron transverse momentum") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return (univ.GetElectronPt() > 0.2 && univ.GetElectronPt() < 1.6);
    }
  };

  //Electron energy, kevin suggested i hit a SUPER aggressive cut just to isolate my events. screw it lets do it
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class ElectronEnergy: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    ElectronEnergy(const double min): PlotUtils::Cut<UNIVERSE, EVENT>("Lepton energy > " + utils::to_string_trimmed(min) + " GeV"), fMin(min) {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetElectronEnergy() > fMin;
    }
    const double fMin;
  };

  //EleptonSin2Theta
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class EleptonSin2Theta: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    EleptonSin2Theta(): PlotUtils::Cut<UNIVERSE, EVENT>("E_lep * Theta_lep^2 > 0.003") {} //No sin involved at all??? double check ryan's...

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
  class ProtonTheta: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    ProtonTheta(): PlotUtils::Cut<UNIVERSE, EVENT>("Reco proton angle < 70 deg") {}
    
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
    Apothem(const double apothem): PlotUtils::Cut<UNIVERSE, EVENT>(std::string("Apothem ") + utils::to_string_trimmed(apothem)), fApothem(apothem), fSlope(-1./sqrt(3.))//A regular hexagon has angles of 2*M_PI/3, so I can find this is 1/tan(M_PI/3.)
      {
      }

    private:
      bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
      {
        const ROOT::Math::XYZTVector vertex = univ.GetVertex();
        return (fabs(vertex.y()) < fSlope*fabs(vertex.x()) + 2.*fApothem/sqrt(3.)) && (fabs(vertex.x()) < fApothem);
	//return univ.GetWithinFiducialApothem();
      }

      const double fApothem;
      const double fSlope; 
  };


  //EMLikeTrackScore
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class EMLikeTrackScore: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    EMLikeTrackScore(const double min): PlotUtils::Cut<UNIVERSE, EVENT>("EM Shower candidate score > " + utils::to_string_trimmed(min)), fMin(min) {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return (univ.GetEMLikeShowerScore() > fMin);
    }
    const double fMin;
  };

  //DSCalVisE
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class DSCalVisE: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    DSCalVisE(const double max): PlotUtils::Cut<UNIVERSE, EVENT>("DSCalVisE Ratio < " + utils::to_string_trimmed(max)), fMax(max) {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return (univ.GetDSCalVisE() <= fMax);
    }
    const double fMax;
  };

  //ODCalVisE
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class ODCalVisE: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    ODCalVisE(const double max): PlotUtils::Cut<UNIVERSE, EVENT>("ODCalVisE Ratio < " + utils::to_string_trimmed(max)), fMax(max) {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return (univ.GetODCalVisE() <= fMax);
    }
    const double fMax;
  };

  //Afterpulsing
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class Afterpulsing: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    Afterpulsing(const double min): PlotUtils::Cut<UNIVERSE, EVENT>("Afterpulsing"), fMin(min) {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return (univ.GetFirstFireFraction() >= fMin);
    }

    const double fMin;
  };

  //NonMIPClusterFraction
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class NonMIPClusterFraction: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    NonMIPClusterFraction(const double min): PlotUtils::Cut<UNIVERSE, EVENT>("Lepton Non MIP cluster fraction > " + utils::to_string_trimmed(min)), fMin(min) {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return (univ.GetNonMIPClusFrac() > fMin); //testing 0.7, was 0.4
    }
    const double fMin;
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
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class VertexTrackMultiplicity: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    VertexTrackMultiplicity(const int min, const int max): PlotUtils::Cut<UNIVERSE, EVENT>("Vertex Track multiplicity: " + utils::to_string_trimmed(min) + " <= n <= " + utils::to_string_trimmed(max)), fMin(min), fMax(max) {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return (univ.GetVertexTrackMultiplicity() >= fMin && univ.GetVertexTrackMultiplicity() <= fMax);
    }
    const int fMin;
    const int fMax;
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
    TransverseGapScore(const double min): PlotUtils::Cut<UNIVERSE, EVENT>("Lepton Transverse gap score > " + utils::to_string_trimmed(min)), fMin(min) {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetTransverseGapScore() > fMin;
    }
    const double fMin;
  };

  //MeanFrontdEdX
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class MeanFrontdEdX: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    MeanFrontdEdX(const double max): PlotUtils::Cut<UNIVERSE, EVENT>("Mean front shower dE/dX < " + utils::to_string_trimmed(max) +" MeV/cm"), fMax(max) {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetMeanFrontdEdx() < fMax;
    }
    const double fMax;
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

  //Michel electron cut, also used for michel sideband
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class MichelCut: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    MichelCut(): PlotUtils::Cut<UNIVERSE, EVENT>("Has No Michel") {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetImprovedNMichel()==0;
    }
  };

  //Number of nonvtx isolated blobs cut
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class NIsoBlobs: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    NIsoBlobs(const double max): PlotUtils::Cut<UNIVERSE, EVENT>("# of nonvtx iso blobs < " + utils::to_string_trimmed(max)), fMax(max) {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetNIsoBlobs() < fMax;
    }
    const double fMax;
  };

  //Cut on total energy of nonvtx isolated blobs
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class IsoBlobEnergy: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    IsoBlobEnergy(const double max): PlotUtils::Cut<UNIVERSE, EVENT>("Total energy of nonvtx iso blobs < " + utils::to_string_trimmed(max) + " MeV"), fMax(max) {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetIsoBlobsEnergy() < fMax;
    }
    const double fMax;
  };

  //Cut on any isolated blob more than some distance upstream of the neutrino vertex
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class UpstreamIsoBlob: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    UpstreamIsoBlob(const double min): PlotUtils::Cut<UNIVERSE, EVENT>("No Iso blobs more than " + utils::to_string_trimmed(-1*min) + " mm upstream of neutrino vertex"), fMin(min) {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetBlobZDiffToVtxZ() > fMin;
    }
    const double fMin;
  };

  //Elastically Scattered Contained (ESC) proton cut
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class ESC: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
    public:
    ESC(const double max): PlotUtils::Cut<UNIVERSE, EVENT>("Hybrid Chi^2/ESC Proton Cut (chi2 < " + utils::to_string_trimmed(max) + ")"), fMax(max) {}

    private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetProtonESCNodeChi2() < fMax;
    }
    const double fMax;
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
