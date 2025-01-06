#!/usr/bin/python

#USAGE: backgroundStack.py <dataFile.root> <mcFile.root>

#from ROOT import *
import ROOT
from ROOT import PlotUtils

import sys
from array import array
import math

def CalcApothem(x,y):
    x=abs(x)
    y=abs(y)
    if ( x == 0 or y/x > 1/math.sqrt(3)):
        return (y+x/math.sqrt(3))/2*math.sqrt(3)
    else:
        return x

#Determines if an event is signal based on truth quantities
#Right now this means incoming electron neutrino (currently excludes antineutrino), charged current, and with true vtx within the tracker & apothem cut
def isTrueSignal(event, true_apothem):
    if event.mc_vtx[2] > 5980 and event.mc_vtx[2] < 8422 and true_apothem < 850 and event.mc_incoming==12 and event.mc_current==1:
        return True
    else:
        return False


#basic, self explanatory cuts that I don't really need to look at specific histograms for (zrange in tracker, within 850cm apothem and deadtime cut
def passesPreCuts(event, apothem):
    if event.vtx[2] > 5980 and event.vtx[2] < 8422 and apothem < 850 and event.phys_n_dead_discr_pair_upstream_prim_track_proj==0:
        return True
    else:
        return False
    
def EventLoopToMakeHistos():
    daChain = ROOT.TChain("MasterAnaDev")
    
#    daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110030_Playlist.root")
#    daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110031_Playlist.root")
#    daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110032_Playlist.root")
#    daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110033_Playlist.root")
#    daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110034_Playlist.root")
#    daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110035_Playlist.root")
    daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110036_Playlist.root")
    daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110037_Playlist.root")
    daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110038_Playlist.root")
    daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110039_Playlist.root")
    daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110040_Playlist.root")
    #/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110038_Playlist.root
    #/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110039_Playlist.root
    #/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110040_Playlist.root


    daChain.SetBranchStatus("*", False)

    daChain.SetBranchStatus("vtx", True)
    daChain.SetBranchStatus("phys_n_dead_discr_pair_upstream_prim_track_proj", True)
    daChain.SetBranchStatus("blob_recoil_E_tracker", True)
    daChain.SetBranchStatus("blob_recoil_E_ecal", True)
    daChain.SetBranchStatus("prong_part_E", True)
    daChain.SetBranchStatus("HasNoBackExitingTracks", True)
    daChain.SetBranchStatus("prong_part_score", True)
    daChain.SetBranchStatus("prong_TransverseGapScore", True)
    daChain.SetBranchStatus("prong_dEdXMeanFrontTracker", True)

    daChain.SetBranchStatus("mc_vtx", True)
    daChain.SetBranchStatus("mc_incoming", True)
    daChain.SetBranchStatus("mc_current", True)

    #Initialize Histograms, there's a lot... :/
    #for every cut (before & after), I have a histogram in reco E_nu for all selected events, and only selected signal events. 
    #Combined backgrounds are then selected - signal
    #Then, for each cut (again before & after) I also have a histogram in the VARIABLE that we're cutting on, one each for selected and signal. 
    #total is 28 histograms

    energy_bins=[0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,9,10,12.5,15,17.5,20]
    n_energy_bins=len(energy_bins)-1

    emscore_bins=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
    n_emscore_bins=len(emscore_bins)-1

    gap_score_bins=[10 * i for i in range(20)]
    n_gap_score_bins=len(gap_score_bins)-1

    dEdX_bins=[0.2 * i for i in range(26)]
    n_dEdX_bins=len(dEdX_bins)-1

    #selected events at various cut stages (reco E_nu)
    reco_E_nu_selected_PRE_EAVAIL_CUT = ROOT.TH1D( 'reco_E_nu_selected_PRE_EAVAIL_CUT', 'Reco E_nu', n_energy_bins, array('d', energy_bins))
    reco_E_nu_selected_POST_EAVAIL_CUT = ROOT.TH1D( 'reco_E_nu_selected_POST_EAVAIL_CUT', 'Reco E_nu', n_energy_bins, array('d', energy_bins))
    reco_E_nu_selected_POST_BACK_TRACK_REJECTION = ROOT.TH1D( 'reco_E_nu_selected_POST_BACK_TRACK_REJECTION', 'Reco E_nu', n_energy_bins, array('d', energy_bins))
    reco_E_nu_selected_POST_SHOWER_SCORE_CUT = ROOT.TH1D( 'reco_E_nu_selected_POST_SHOWER_SCORE_CUT', 'Reco E_nu', n_energy_bins, array('d', energy_bins))
    reco_E_nu_selected_POST_TRANSVERSE_GAP_CUT = ROOT.TH1D( 'reco_E_nu_selected_POST_TRANSVERSE_GAP_CUT', 'Reco E_nu', n_energy_bins, array('d', energy_bins))
    reco_E_nu_selected_POST_FRONTdEdX_CUT = ROOT.TH1D( 'reco_E_nu_selected_POST_FRONTdEdX_CUT', 'Reco E_nu', n_energy_bins, array('d', energy_bins))

    #true signal events (which were also selected) at various cut stages (reco E_nu)
    reco_E_nu_signal_selected_PRE_EAVAIL_CUT = ROOT.TH1D( 'reco_E_nu_signal_selected_PRE_EAVAIL_CUT', 'Reco E_nu', n_energy_bins, array('d', energy_bins))
    reco_E_nu_signal_selected_POST_EAVAIL_CUT = ROOT.TH1D( 'reco_E_nu_signal_selected_POST_EAVAIL_CUT', 'Reco E_nu', n_energy_bins, array('d', energy_bins))
    reco_E_nu_signal_selected_POST_BACK_TRACK_REJECTION = ROOT.TH1D( 'reco_E_nu_signal_selected_POST_BACK_TRACK_REJECTION', 'Reco E_nu', n_energy_bins, array('d', energy_bins))
    reco_E_nu_signal_selected_POST_SHOWER_SCORE_CUT = ROOT.TH1D( 'reco_E_nu_signal_selected_POST_SHOWER_SCORE_CUT', 'Reco E_nu', n_energy_bins, array('d', energy_bins))
    reco_E_nu_signal_selected_POST_TRANSVERSE_GAP_CUT = ROOT.TH1D( 'reco_E_nu_signal_selected_POST_TRANSVERSE_GAP_CUT', 'Reco E_nu', n_energy_bins, array('d', energy_bins))
    reco_E_nu_signal_selected_POST_FRONTdEdX_CUT = ROOT.TH1D( 'reco_E_nu_signal_selected_POST_FRONTdEdX_CUT', 'Reco E_nu', n_energy_bins, array('d', energy_bins))

    #E_avail
    reco_E_avail_selected_PRE_EAVAIL_CUT = ROOT.TH1D( 'reco_E_Avail_selected_PRE_EAVAIL_CUT', 'Reco E_avail', n_energy_bins, array('d', energy_bins))
    reco_E_avail_selected_POST_EAVAIL_CUT = ROOT.TH1D( 'reco_E_Avail_selected_POST_EAVAIL_CUT', 'Reco E_avail', n_energy_bins, array('d', energy_bins))

    reco_E_avail_signal_selected_PRE_EAVAIL_CUT = ROOT.TH1D( 'reco_E_avail_signal_selected_PRE_EAVAIL_CUT', 'Reco E_nu', n_energy_bins, array('d', energy_bins))
    reco_E_avail_signal_selected_POST_EAVAIL_CUT = ROOT.TH1D( 'reco_E_avail_signal_selected_POST_EAVAIL_CUT', 'Reco E_nu', n_energy_bins, array('d', energy_bins))

    #EMShower_score
    EMShower_score_selected_PRE_SHOWER_SCORE_CUT = ROOT.TH1D( 'EMShower_score_selected_PRE_SHOWER_SCORE_CUT', 'Highest EMLike Shower Score', n_emscore_bins, array('d', emscore_bins))
    EMShower_score_selected_POST_SHOWER_SCORE_CUT = ROOT.TH1D( 'EMShower_score_selected_POST_SHOWER_SCORE_CUT', 'Highest EMLike Shower Score', n_emscore_bins, array('d', emscore_bins))

    EMShower_score_signal_selected_PRE_SHOWER_SCORE_CUT = ROOT.TH1D( 'EMShower_score_signal_selected_PRE_SHOWER_SCORE_CUT', 'Highest EMLike Shower Score', n_emscore_bins, array('d', emscore_bins))
    EMShower_score_signal_selected_POST_SHOWER_SCORE_CUT = ROOT.TH1D( 'EMShower_score_signal_selected_POST_SHOWER_SCORE_CUT', 'Highest EMLike Shower Score', n_emscore_bins, array('d', emscore_bins))

    #transverse_gap_score
    transverse_gap_score_selected_PRE_GAP_CUT = ROOT.TH1D( 'transverse_gap_score_selected_PRE_GAP_CUT', 'Transverse Gap Score', n_gap_score_bins, array('d', gap_score_bins))
    transverse_gap_score_selected_POST_GAP_CUT = ROOT.TH1D( 'transverse_gap_score_selected_POST_GAP_CUT', 'Transverse Gap Score', n_gap_score_bins, array('d', gap_score_bins))

    transverse_gap_score_signal_selected_PRE_GAP_CUT = ROOT.TH1D( 'transverse_gap_score_signal_selected_PRE_GAP_CUT', 'Transverse Gap Score', n_gap_score_bins, array('d', gap_score_bins))
    transverse_gap_score_signal_selected_POST_GAP_CUT = ROOT.TH1D( 'transverse_gap_score_signal_selected_POST_GAP_CUT', 'Transverse Gap Score', n_gap_score_bins, array('d', gap_score_bins))


    #frontdEdX
    frontdEdX_selected_PRE_FRONTdEdX_CUT = ROOT.TH1D( 'frontdEdX_selected_PRE_FRONTdEdX_CUT', 'Mean Front dE/dx', n_dEdX_bins, array('d', dEdX_bins))
    frontdEdX_selected_POST_FRONTdEdX_CUT = ROOT.TH1D( 'frontdEdX_selected_POST_FRONTdEdX_CUT', 'Mean Front dE/dx', n_dEdX_bins, array('d', dEdX_bins))

    frontdEdX_signal_selected_PRE_FRONTdEdX_CUT = ROOT.TH1D( 'frontdEdX_signal_selected_PRE_FRONTdEdX_CUT', 'Mean Front dE/dx', n_dEdX_bins, array('d', dEdX_bins))
    frontdEdX_signal_selected_POST_FRONTdEdX_CUT = ROOT.TH1D( 'frontdEdX_signal_selected_POST_FRONTdEdX_CUT', 'Mean Front dE/dx', n_dEdX_bins, array('d', dEdX_bins))

    #Ok now I'm going to start looping through, one cut at a time, and filling histos. 

    for i, event in enumerate(daChain):
        #limit to 1k events for testing purposes
        #if i>1000:
        #    break

        #print statements to make sure its still running
        if i%1000==0:
            print "Currently on event #", i
        #Calculating some cut variables
        apothem = CalcApothem(event.vtx[0], event.vtx[1])
        true_apothem = CalcApothem(event.mc_vtx[0], event.mc_vtx[1])

        E_avail = ((event.blob_recoil_E_tracker + event.blob_recoil_E_ecal)*1.17 - (0.008*event.prong_part_E[0][3]+5))/1000 #in GeV, so divide by 1e3
        
        reco_E_nu = E_avail + event.prong_part_E[0][3]/1000
        
        highest_EMScore = 0
        for score in event.prong_part_score:
            if score > highest_EMScore:
                highest_EMScore = score

        #First just the pre cuts, this would be classified as PRE_EAVAIL
        if passesPreCuts(event, apothem):
            reco_E_nu_selected_PRE_EAVAIL_CUT.Fill(reco_E_nu)
            reco_E_avail_selected_PRE_EAVAIL_CUT.Fill(E_avail)
            if isTrueSignal(event, true_apothem):
                reco_E_nu_signal_selected_PRE_EAVAIL_CUT.Fill(reco_E_nu)
                reco_E_avail_signal_selected_PRE_EAVAIL_CUT.Fill(E_avail)

        #POST_EAVAIL (also pre back exiting tracks but im not plotting that
        if passesPreCuts(event, apothem) and E_avail < 1.5:
            reco_E_nu_selected_POST_EAVAIL_CUT.Fill(reco_E_nu)
            reco_E_avail_selected_POST_EAVAIL_CUT.Fill(E_avail)
            if isTrueSignal(event, true_apothem):
                reco_E_nu_signal_selected_POST_EAVAIL_CUT.Fill(reco_E_nu)
                reco_E_avail_signal_selected_POST_EAVAIL_CUT.Fill(E_avail)

        #POST_BACK_TRACK_REJECTION and also PRE_SHOWER_SCORE_CUT
        if passesPreCuts(event, apothem) and E_avail < 1.5 and event.HasNoBackExitingTracks==1:
            reco_E_nu_selected_POST_BACK_TRACK_REJECTION.Fill(reco_E_nu)
            EMShower_score_selected_PRE_SHOWER_SCORE_CUT.Fill(highest_EMScore)
            if isTrueSignal(event, true_apothem):
                reco_E_nu_signal_selected_POST_BACK_TRACK_REJECTION.Fill(reco_E_nu)
                EMShower_score_signal_selected_PRE_SHOWER_SCORE_CUT.Fill(highest_EMScore)
                
        #POST_SHOWER_SCORE_CUT and also PRE_GAP_CUT
        if passesPreCuts(event, apothem) and E_avail < 1.5 and event.HasNoBackExitingTracks==1 and highest_EMScore > 0.7:
            reco_E_nu_selected_POST_SHOWER_SCORE_CUT.Fill(reco_E_nu)
            EMShower_score_selected_POST_SHOWER_SCORE_CUT.Fill(highest_EMScore)
            transverse_gap_score_selected_PRE_GAP_CUT.Fill(event.prong_TransverseGapScore[0])
            if isTrueSignal(event, true_apothem):
                reco_E_nu_signal_selected_POST_SHOWER_SCORE_CUT.Fill(reco_E_nu)
                EMShower_score_signal_selected_POST_SHOWER_SCORE_CUT.Fill(highest_EMScore)
                transverse_gap_score_signal_selected_PRE_GAP_CUT.Fill(event.prong_TransverseGapScore[0])
                 
        #POST_GAP_CUT and also PRE_FRONTdEdX_CUT
        if passesPreCuts(event, apothem) and E_avail < 1.5 and event.HasNoBackExitingTracks==1 and highest_EMScore > 0.7 and event.prong_TransverseGapScore[0] > 15:
            reco_E_nu_selected_POST_TRANSVERSE_GAP_CUT.Fill(reco_E_nu)
            transverse_gap_score_selected_POST_GAP_CUT.Fill(event.prong_TransverseGapScore[0])
            frontdEdX_selected_PRE_FRONTdEdX_CUT.Fill(event.prong_dEdXMeanFrontTracker[0])
            if isTrueSignal(event, true_apothem):
                reco_E_nu_signal_selected_POST_TRANSVERSE_GAP_CUT.Fill(reco_E_nu)
                transverse_gap_score_signal_selected_POST_GAP_CUT.Fill(event.prong_TransverseGapScore[0])
                frontdEdX_signal_selected_PRE_FRONTdEdX_CUT.Fill(event.prong_dEdXMeanFrontTracker[0])

        #POST_FRONTdEdX_CUT
        if passesPreCuts(event, apothem) and E_avail < 1.5 and event.HasNoBackExitingTracks==1 and highest_EMScore > 0.7 and event.prong_TransverseGapScore[0] > 15 and event.prong_dEdXMeanFrontTracker[0] < 2.4:
            reco_E_nu_selected_POST_FRONTdEdX_CUT.Fill(reco_E_nu)
            frontdEdX_selected_POST_FRONTdEdX_CUT.Fill(event.prong_dEdXMeanFrontTracker[0])
            if isTrueSignal(event, true_apothem):
                reco_E_nu_signal_selected_POST_FRONTdEdX_CUT.Fill(reco_E_nu)
                frontdEdX_signal_selected_POST_FRONTdEdX_CUT.Fill(event.prong_dEdXMeanFrontTracker[0])

    #Now, using selected - signal, going to make 14 more histos as bkgd
    reco_E_nu_background_PRE_EAVAIL_CUT = ROOT.TH1D.Clone(reco_E_nu_selected_PRE_EAVAIL_CUT)
    reco_E_nu_background_PRE_EAVAIL_CUT.Add(reco_E_nu_signal_selected_PRE_EAVAIL_CUT, -1)
    reco_E_nu_background_PRE_EAVAIL_CUT.SetName("reco_E_nu_background_PRE_EAVAIL_CUT")

    reco_E_nu_background_POST_EAVAIL_CUT = ROOT.TH1D.Clone(reco_E_nu_selected_POST_EAVAIL_CUT)
    reco_E_nu_background_POST_EAVAIL_CUT.Add(reco_E_nu_signal_selected_POST_EAVAIL_CUT, -1)
    reco_E_nu_background_POST_EAVAIL_CUT.SetName("reco_E_nu_background_PRE_EAVAIL_CUT")

    reco_E_nu_background_POST_BACK_TRACK_REJECTION = ROOT.TH1D.Clone(reco_E_nu_selected_POST_BACK_TRACK_REJECTION)
    reco_E_nu_background_POST_BACK_TRACK_REJECTION.Add(reco_E_nu_signal_selected_POST_BACK_TRACK_REJECTION, -1)
    reco_E_nu_background_POST_BACK_TRACK_REJECTION.SetName("reco_E_nu_background_POST_BACK_TRACK_REJECTION")

    reco_E_nu_background_POST_SHOWER_SCORE_CUT = ROOT.TH1D.Clone(reco_E_nu_selected_POST_SHOWER_SCORE_CUT)
    reco_E_nu_background_POST_SHOWER_SCORE_CUT.Add(reco_E_nu_signal_selected_POST_SHOWER_SCORE_CUT, -1)
    reco_E_nu_background_POST_SHOWER_SCORE_CUT.SetName("reco_E_nu_background_POST_SHOWER_SCORE_CUT")

    reco_E_nu_background_POST_TRANSVERSE_GAP_CUT = ROOT.TH1D.Clone(reco_E_nu_selected_POST_TRANSVERSE_GAP_CUT)
    reco_E_nu_background_POST_TRANSVERSE_GAP_CUT.Add(reco_E_nu_signal_selected_POST_TRANSVERSE_GAP_CUT, -1)
    reco_E_nu_background_POST_TRANSVERSE_GAP_CUT.SetName("reco_E_nu_background_POST_TRANSVERSE_GAP_CUT")

    reco_E_nu_background_POST_FRONTdEdX_CUT = ROOT.TH1D.Clone(reco_E_nu_selected_POST_FRONTdEdX_CUT)
    reco_E_nu_background_POST_FRONTdEdX_CUT.Add(reco_E_nu_signal_selected_POST_FRONTdEdX_CUT, -1)
    reco_E_nu_background_POST_FRONTdEdX_CUT.SetName("reco_E_nu_background_POST_FRONTdEdX_CUT")

    reco_E_avail_background_PRE_EAVAIL_CUT = ROOT.TH1D.Clone(reco_E_avail_selected_PRE_EAVAIL_CUT)
    reco_E_avail_background_PRE_EAVAIL_CUT.Add(reco_E_avail_signal_selected_PRE_EAVAIL_CUT, -1)
    reco_E_avail_background_PRE_EAVAIL_CUT.SetName("reco_E_avail_background_PRE_EAVAIL_CUT")

    reco_E_avail_background_POST_EAVAIL_CUT = ROOT.TH1D.Clone(reco_E_avail_selected_POST_EAVAIL_CUT)
    reco_E_avail_background_POST_EAVAIL_CUT.Add(reco_E_avail_signal_selected_POST_EAVAIL_CUT, -1)
    reco_E_avail_background_POST_EAVAIL_CUT.SetName("reco_E_avail_background_POST_EAVAIL_CUT")

    EMShower_score_background_PRE_SHOWER_SCORE_CUT = ROOT.TH1D.Clone(EMShower_score_selected_PRE_SHOWER_SCORE_CUT)
    EMShower_score_background_PRE_SHOWER_SCORE_CUT.Add(EMShower_score_signal_selected_PRE_SHOWER_SCORE_CUT, -1)
    EMShower_score_background_PRE_SHOWER_SCORE_CUT.SetName("EMShower_score_background_PRE_SHOWER_SCORE_CUT")

    EMShower_score_background_POST_SHOWER_SCORE_CUT = ROOT.TH1D.Clone(EMShower_score_selected_POST_SHOWER_SCORE_CUT)
    EMShower_score_background_POST_SHOWER_SCORE_CUT.Add(EMShower_score_signal_selected_POST_SHOWER_SCORE_CUT, -1)
    EMShower_score_background_POST_SHOWER_SCORE_CUT.SetName("EMShower_score_background_POST_SHOWER_SCORE_CUT")

    transverse_gap_score_background_PRE_GAP_CUT = ROOT.TH1D.Clone(transverse_gap_score_selected_PRE_GAP_CUT)
    transverse_gap_score_background_PRE_GAP_CUT.Add(transverse_gap_score_signal_selected_PRE_GAP_CUT, -1)
    transverse_gap_score_background_PRE_GAP_CUT.SetName("transverse_gap_score_background_PRE_GAP_CUT")

    transverse_gap_score_background_POST_GAP_CUT = ROOT.TH1D.Clone(transverse_gap_score_selected_POST_GAP_CUT)
    transverse_gap_score_background_POST_GAP_CUT.Add(transverse_gap_score_signal_selected_POST_GAP_CUT, -1)
    transverse_gap_score_background_POST_GAP_CUT.SetName("transverse_gap_score_background_POST_GAP_CUT")

    frontdEdX_background_PRE_FRONTdEdX_CUT = ROOT.TH1D.Clone(frontdEdX_selected_PRE_FRONTdEdX_CUT)
    frontdEdX_background_PRE_FRONTdEdX_CUT.Add(frontdEdX_signal_selected_PRE_FRONTdEdX_CUT, -1)
    frontdEdX_background_PRE_FRONTdEdX_CUT.SetName("frontdEdX_background_PRE_FRONTdEdX_CUT")

    frontdEdX_background_POST_FRONTdEdX_CUT = ROOT.TH1D.Clone(frontdEdX_selected_POST_FRONTdEdX_CUT)
    frontdEdX_background_POST_FRONTdEdX_CUT.Add(frontdEdX_signal_selected_POST_FRONTdEdX_CUT, -1)
    frontdEdX_background_POST_FRONTdEdX_CUT.SetName("frontdEdX_background_POST_FRONTdEdX_CUT")

    #k now that all that nonsense is over with, going to save these to a root file for future use
    outFile = ROOT.TFile('CutSummaryHistos.root', 'RECREATE')

    reco_E_nu_selected_PRE_EAVAIL_CUT.Write()
    reco_E_nu_selected_POST_EAVAIL_CUT.Write()
    reco_E_nu_selected_POST_BACK_TRACK_REJECTION.Write()
    reco_E_nu_selected_POST_SHOWER_SCORE_CUT.Write()
    reco_E_nu_selected_POST_TRANSVERSE_GAP_CUT.Write()
    reco_E_nu_selected_POST_FRONTdEdX_CUT.Write()
    
    reco_E_nu_signal_selected_PRE_EAVAIL_CUT.Write()
    reco_E_nu_signal_selected_POST_EAVAIL_CUT.Write()
    reco_E_nu_signal_selected_POST_BACK_TRACK_REJECTION.Write()
    reco_E_nu_signal_selected_POST_SHOWER_SCORE_CUT.Write()
    reco_E_nu_signal_selected_POST_TRANSVERSE_GAP_CUT.Write()
    reco_E_nu_signal_selected_POST_FRONTdEdX_CUT.Write()

    reco_E_nu_background_PRE_EAVAIL_CUT.Write()
    reco_E_nu_background_POST_EAVAIL_CUT.Write()
    reco_E_nu_background_POST_BACK_TRACK_REJECTION.Write()
    reco_E_nu_background_POST_SHOWER_SCORE_CUT.Write()
    reco_E_nu_background_POST_TRANSVERSE_GAP_CUT.Write()
    reco_E_nu_background_POST_FRONTdEdX_CUT.Write()
    
    reco_E_avail_selected_PRE_EAVAIL_CUT.Write()
    reco_E_avail_selected_POST_EAVAIL_CUT.Write()
    reco_E_avail_signal_selected_PRE_EAVAIL_CUT.Write()
    reco_E_avail_signal_selected_POST_EAVAIL_CUT.Write()
    reco_E_avail_background_PRE_EAVAIL_CUT.Write()
    reco_E_avail_background_POST_EAVAIL_CUT.Write()
        
    EMShower_score_selected_PRE_SHOWER_SCORE_CUT.Write()
    EMShower_score_selected_POST_SHOWER_SCORE_CUT.Write()
    EMShower_score_signal_selected_PRE_SHOWER_SCORE_CUT.Write()
    EMShower_score_signal_selected_POST_SHOWER_SCORE_CUT.Write()
    EMShower_score_background_PRE_SHOWER_SCORE_CUT.Write()
    EMShower_score_background_POST_SHOWER_SCORE_CUT.Write()
    
    transverse_gap_score_selected_PRE_GAP_CUT.Write()
    transverse_gap_score_selected_POST_GAP_CUT.Write()
    transverse_gap_score_signal_selected_PRE_GAP_CUT.Write()
    transverse_gap_score_signal_selected_POST_GAP_CUT.Write()
    transverse_gap_score_background_PRE_GAP_CUT.Write()
    transverse_gap_score_background_POST_GAP_CUT.Write()

    frontdEdX_selected_PRE_FRONTdEdX_CUT.Write()
    frontdEdX_selected_POST_FRONTdEdX_CUT.Write()
    frontdEdX_signal_selected_PRE_FRONTdEdX_CUT.Write()
    frontdEdX_signal_selected_POST_FRONTdEdX_CUT.Write()
    frontdEdX_background_PRE_FRONTdEdX_CUT.Write()
    frontdEdX_background_POST_FRONTdEdX_CUT.Write()
    
    outFile.Close()

selected_histos = ["reco_E_nu_selected_PRE_EAVAIL_CUT", "reco_E_nu_selected_POST_EAVAIL_CUT", "reco_E_nu_selected_POST_BACK_TRACK_REJECTION", "reco_E_nu_selected_POST_SHOWER_SCORE_CUT", "reco_E_nu_selected_POST_TRANSVERSE_GAP_CUT", "reco_E_nu_selected_POST_FRONTdEdX_CUT", "reco_E_avail_selected_PRE_EAVAIL_CUT"]
signal_histos = []

#EventLoopToMakeHistos()

#TH1.AddDirectory(False)

mcFile = ROOT.TFile.Open(sys.argv[1])
#mcFile = ROOT.TFile.Open("CutSummaryHistos.root")
plotter = PlotUtils.MnvPlotter()
plotter.ApplyStyle(PlotUtils.kCCNuPionIncStyle)

#mcPOT = mcFile.Get("POTUsed").GetVal()
#dataPOT = dataFile.Get("POTUsed").GetVal()

#mcScale = dataPOT/mcPOT

#I have no idea what I'm doing
#plotter.DrawDataMC(dataHist, mcHist, mcScale, "TR")

#I need to pass this guy a TObjArray of two histograms -> one that is pTmu_selected_signal_reco, and one that is all of the backgrounds
#for all the backgrounds I can either add up the 3 background histograms (pTmu_background_Wrongsign+pTmu_background_NC+pTmu_background_other)
#or I can do all minus signal (pTmu_data - pTmu_selected_signal_reco)

selected = mcFile.Get('pTmu_data')
#true_signal = mcFile.Get('pTmu_selected_signal_reco')

#selected = mcFile.Get('E_avail_data')
#true_signal = mcFile.Get('E_avail_selected_signal_reco')

#selected = mcFile.Get('')
#backgrounds = PlotUtils.MnvH1D.Clone(selected)
#backgrounds.Add(true_signal, -1)

#true_signal.SetFillColor(kGreen)
#backgrounds.SetFillColor(kRed)

signal = mcFile.Get('E_avail_selected_signal_reco')
bkg1 = mcFile.Get('E_avail_background_NueCC_0pi_0p')
bkg2 = mcFile.Get('E_avail_background_Other_NueCC')
bkg3 = mcFile.Get('E_avail_background_NumuCC')
bkg4 = mcFile.Get('E_avail_background_NC_with_pi0')
bkg5 = mcFile.Get('E_avail_background_Other')


#selected.SetTitle('selected')
#true_signal.SetTitle('signal')
#backgrounds.SetTitle('background')

array = ROOT.TObjArray()

canvASS = ROOT.TCanvas( 'canvASS', 'background', 200, 10, 700, 500 )
canvASS.cd()

#array.Add(true_signal)
#array.Add(backgrounds)
array.Add(signal)
array.Add(bkg1)
array.Add(bkg2)
array.Add(bkg3)
array.Add(bkg4)
array.Add(bkg5)

plotter.DrawStackedMC(array, 1, 'TR', 2, 1, 3001, "", "")
#plotter.WritePreliminary()

