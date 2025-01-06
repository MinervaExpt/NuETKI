#!/usr/bin/python

#USAGE: backgroundStack.py <dataFile.root> <mcFile.root>

#in progress edit to plotCuts.py trying to make it so I dont have to manually initialize like 8000 histograms, currently unfinished im on like line 99
#tbh might not even finish it cause I should move back to the actual MAT event loop at some point. 
#from ROOT import *

#to make it not display the canvases as it draws and saves them, saves a bunch of time
import sys
sys.argv.append( '-b' )

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

#proton pdg code
#mesons to exclude: pi0 (111), pi+ (211), pi- (
def isTrueSignalWithProton(event, true_apothem):
    hasProton = False
    hasMeson = False
    for pdg in event.mc_FSPartPDG:
        part = abs(pdg)
        if pdg==2212:
            hasProton = True
        #Checking for mesons here, just pions & kaons
        elif part==211 or part==111 or part==321 or part==311 or part==130:
            hasMeson = True
    if event.mc_vtx[2] > 5980 and event.mc_vtx[2] < 8422 and true_apothem < 850 and abs(event.mc_incoming)==12 and event.mc_current==1 and hasProton and not hasMeson:
        return True
    else:
        return False

def isTrueSignalWithoutProton(event, true_apothem):
    hasProton = False
    hasMeson = False
    for pdg in event.mc_FSPartPDG:
        part = abs(pdg)
        if pdg==2212:
            hasProton = True
        #Checking for mesons here, just pions & kaons                 
        elif part==211 or part==111 or part==321 or part==311 or part==130:
            hasMeson = True
    if event.mc_vtx[2] > 5980 and event.mc_vtx[2] < 8422 and true_apothem < 850 and abs(event.mc_incoming)==12 and event.mc_current==1 and not hasProton and not hasMeson:
        return True
    else:
        return False

def isTrueOtherNueCC(event, true_apothem):
    hasProton = False
    hasMeson = False
    for pdg in event.mc_FSPartPDG:
        part = abs(pdg)
        if pdg==2212:
            hasProton = True
        #Checking for mesons here, just pions & kaons                 
        elif part==211 or part==111 or part==321 or part==311 or part==130:
            hasMeson = True
    if event.mc_vtx[2] > 5980 and event.mc_vtx[2] < 8422 and true_apothem < 850 and abs(event.mc_incoming)==12 and event.mc_current==1 and hasMeson:
        return True
    else:
        return False


def isTrueNumuCC(event, true_apothem):
    if event.mc_vtx[2] > 5980 and event.mc_vtx[2] < 8422 and true_apothem < 850 and abs(event.mc_incoming)==14 and event.mc_current==1:
        return True
    else:
        return False

def isTrueNCpi0(event, true_apothem):
    hasPi0 = False
    for pdg in event.mc_FSPartPDG:
        part = abs(pdg)
        if part==111:
            hasPi0= True
    if event.mc_vtx[2] > 5980 and event.mc_vtx[2] < 8422 and true_apothem < 850 and event.mc_current==2 and hasPi0:
        return True
    else:
        return False


#basic, self explanatory cuts that I don't really need to look at specific histograms for (zrange in tracker, within 850cm apothem and deadtime cut
def passesPreCuts(event, apothem):
    if event.vtx[2] > 5980 and event.vtx[2] < 8422 and apothem < 850 and event.phys_n_dead_discr_pair_upstream_prim_track_proj==0:
        return True
    else:
        return False
    
def EventLoopToMakeHistos(selectedHistoNames):
    daChain = ROOT.TChain("MasterAnaDev")
    
#    daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110030_Playlist.root")
#    daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110031_Playlist.root")
#    daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110032_Playlist.root")
#    daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110033_Playlist.root")
#    daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110034_Playlist.root")
#    daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110035_Playlist.root")
    daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110036_Playlist.root")
    #daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110037_Playlist.root")
    #daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110038_Playlist.root")
    #daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110039_Playlist.root")
    #daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110040_Playlist.root")

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
    daChain.SetBranchStatus("mc_FSPartPDG", True)
    

    #Initialize Histograms, there's a lot... :/

    energy_bins=[0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,9,10,12.5,15,17.5,20]
    n_energy_bins=len(energy_bins)-1

    emscore_bins=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
    n_emscore_bins=len(emscore_bins)-1

    gap_score_bins=[10 * i for i in range(20)]
    n_gap_score_bins=len(gap_score_bins)-1

    dEdX_bins=[0.2 * i for i in range(26)]
    n_dEdX_bins=len(dEdX_bins)-1

    #lists of all of the names of each histogram
    signalWithProtonHistoNames=[]
    signalWithOutProtonHistoNames=[]
    otherNueCCHistoNames=[]
    numuCCHistoNames=[]
    NCpi0HistoNames=[]
    otherBackgroundHistoNames=[]

    for i in selectedHistoNames:
        signalWithProtonHistoNames.append(i.replace("selected", "mesonless_nue_CC_with_proton"))
        signalWithOutProtonHistoNames.append(i.replace("selected", "mesonless_nue_CC_without_proton"))
        otherNueCCHistoNames.append(i.replace("selected", "other_nue_CC"))
        numuCCHistoNames.append(i.replace("selected", "numu_CC"))
        NCpi0HistoNames.append(i.replace("selected", "NC_pi0"))
        otherBackgroundHistoNames.append(i.replace("selected", "other_background"))

    #python lists of the histogram objects themselves, all seven of these (for a single variable at a single cut point) will make up 1 plot when stacked
    selectedHistos=[]
    signalWithProtonHistos=[]
    signalWithOutProtonHistos=[]
    otherNueCCHistos=[]
    numuCCHistos=[]
    NCpi0Histos=[]
    otherBackgroundHistos=[]

    for i in range(len(selectedHistoNames)): 
        #figure out which variable we're plotting (for binning & labelling), then initialize all necessary histos
        if "reco_E_nu" in selectedHistoNames[i]:
            label = "Reco E_nu (GeV)"
            nbins = n_energy_bins
            binning = energy_bins
        elif "reco_E_avail" in selectedHistoNames[i]:
            label = "Reco E_avail (GeV)"
            nbins = n_energy_bins
            binning = energy_bins
        elif "EMShower_score" in selectedHistoNames[i]:
            label = "EMLike Shower Score"
            nbins = n_emscore_bins
            binning = emscore_bins
        elif "transverse_gap_score" in selectedHistoNames[i]:
            label = "Transverse Gap Score (mm?)"
            nbins = n_gap_score_bins
            binning = gap_score_bins
        elif "frontdEdX" in selectedHistoNames[i]:
            label = "Mean Front dE/dX (MeV/cm ??)"
            nbins = n_dEdX_bins
            binning = dEdX_bins

        selectedHistos.append(ROOT.TH1D( selectedHistoNames[i], label, nbins, array('d', binning)))
        signalWithProtonHistos.append(ROOT.TH1D( signalWithProtonHistoNames[i], label, nbins, array('d', binning)))
        signalWithOutProtonHistos.append(ROOT.TH1D( signalWithOutProtonHistoNames[i], label, nbins, array('d', binning)))
        otherNueCCHistos.append(ROOT.TH1D( otherNueCCHistoNames[i], label, nbins, array('d', binning)))
        numuCCHistos.append(ROOT.TH1D( numuCCHistoNames[i], label, nbins, array('d', binning)))
        NCpi0Histos.append(ROOT.TH1D( NCpi0HistoNames[i], label, nbins, array('d', binning)))
        otherBackgroundHistos.append(ROOT.TH1D( otherBackgroundHistoNames[i], label, nbins, array('d', binning)))

    #Ok now I'm going to start looping through the events, do one cut at a time, and filling histos. 

    for i, event in enumerate(daChain):
        #limit to 1k events for testing purposes
        #if i>1000:
            #break

        #print statements to make sure its still running
        if i%1000==0:
            print "Currently on event #", i
        #Calculating some cut variables
        apothem = CalcApothem(event.vtx[0], event.vtx[1])
        true_apothem = CalcApothem(event.mc_vtx[0], event.mc_vtx[1])

        reco_E_avail = ((event.blob_recoil_E_tracker + event.blob_recoil_E_ecal)*1.17 - (0.008*event.prong_part_E[0][3]+5))/1000 #in GeV, so divide by 1e3
        
        reco_E_nu = reco_E_avail + event.prong_part_E[0][3]/1000
        
        #First just the pre cuts, this would be classified as PRE_EAVAIL
        if passesPreCuts(event, apothem):
            selectedHistos[0].Fill(reco_E_nu)
            selectedHistos[6].Fill(reco_E_avail)
            if isTrueSignalWithProton(event, true_apothem):
                signalWithProtonHistos[0].Fill(reco_E_nu)
                signalWithProtonHistos[6].Fill(reco_E_avail)
            elif isTrueSignalWithoutProton(event, true_apothem):
                signalWithOutProtonHistos[0].Fill(reco_E_nu)
                signalWithOutProtonHistos[6].Fill(reco_E_avail)
            elif isTrueOtherNueCC(event, true_apothem):
                otherNueCCHistos[0].Fill(reco_E_nu)
                otherNueCCHistos[6].Fill(reco_E_avail)
            elif isTrueNumuCC(event, true_apothem):
                numuCCHistos[0].Fill(reco_E_nu)
                numuCCHistos[6].Fill(reco_E_avail)
            elif isTrueNCpi0(event, true_apothem):
                NCpi0Histos[0].Fill(reco_E_nu)
                NCpi0Histos[6].Fill(reco_E_avail)
            else:  #category for any other background not already covered
                if event.mc_vtx[2] > 5980 and event.mc_vtx[2] < 8422 and true_apothem < 850 and abs(event.mc_incoming)==14 and event.mc_current==1:
                    print "looks like i didn't cover all of the other nu_e-CC events in other nu_e-CC :(((("
                otherBackgroundHistos[0].Fill(reco_E_nu)
                otherBackgroundHistos[6].Fill(reco_E_avail)


        #POST_EAVAIL (also pre back exiting tracks but im not plotting that
        #if passesPreCuts(event, apothem) and E_avail < 1.5:
            #reco_E_nu_selected_POST_EAVAIL_CUT.Fill(reco_E_nu)
            #reco_E_avail_selected_POST_EAVAIL_CUT.Fill(E_avail)
            #if isTrueSignal(event, true_apothem):
                #reco_E_nu_signal_selected_POST_EAVAIL_CUT.Fill(reco_E_nu)
                #reco_E_avail_signal_selected_POST_EAVAIL_CUT.Fill(E_avail)

        if passesPreCuts(event, apothem) and reco_E_avail < 1.5:
            selectedHistos[1].Fill(reco_E_nu)
            selectedHistos[7].Fill(reco_E_avail)
            if isTrueSignalWithProton(event, true_apothem):
                signalWithProtonHistos[1].Fill(reco_E_nu)
                signalWithProtonHistos[7].Fill(reco_E_avail)
            elif isTrueSignalWithoutProton(event, true_apothem):
                signalWithOutProtonHistos[1].Fill(reco_E_nu)
                signalWithOutProtonHistos[7].Fill(reco_E_avail)
            elif isTrueOtherNueCC(event, true_apothem):
                otherNueCCHistos[1].Fill(reco_E_nu)
                otherNueCCHistos[7].Fill(reco_E_avail)
            elif isTrueNumuCC(event, true_apothem):
                numuCCHistos[1].Fill(reco_E_nu)
                numuCCHistos[7].Fill(reco_E_avail)
            elif isTrueNCpi0(event, true_apothem):
                NCpi0Histos[1].Fill(reco_E_nu)
                NCpi0Histos[7].Fill(reco_E_avail)
            else:  #category for any other background not already covered
                if event.mc_vtx[2] > 5980 and event.mc_vtx[2] < 8422 and true_apothem < 850 and abs(event.mc_incoming)==14 and event.mc_current==1:
                    print "looks like i didn't cover all of the other nu_e-CC events in other nu_e-CC :(((("
                otherBackgroundHistos[1].Fill(reco_E_nu)
                otherBackgroundHistos[7].Fill(reco_E_avail)


        #POST_BACK_TRACK_REJECTION and also PRE_SHOWER_SCORE_CUT
        #if passesPreCuts(event, apothem) and E_avail < 1.5 and event.HasNoBackExitingTracks==1:
            #reco_E_nu_selected_POST_BACK_TRACK_REJECTION.Fill(reco_E_nu)
            #EMShower_score_selected_PRE_SHOWER_SCORE_CUT.Fill(event.prong_part_score[0])
            #if isTrueSignal(event, true_apothem):
                #reco_E_nu_signal_selected_POST_BACK_TRACK_REJECTION.Fill(reco_E_nu)
                #EMShower_score_signal_selected_PRE_SHOWER_SCORE_CUT.Fill(event.prong_part_score[0])
                
        if passesPreCuts(event, apothem) and reco_E_avail < 1.5 and event.HasNoBackExitingTracks==1:
            selectedHistos[2].Fill(reco_E_nu)
            selectedHistos[8].Fill(event.prong_part_score[0])
            if isTrueSignalWithProton(event, true_apothem):
                signalWithProtonHistos[2].Fill(reco_E_nu)
                signalWithProtonHistos[8].Fill(event.prong_part_score[0])
            elif isTrueSignalWithoutProton(event, true_apothem):
                signalWithOutProtonHistos[2].Fill(reco_E_nu)
                signalWithOutProtonHistos[8].Fill(event.prong_part_score[0])
            elif isTrueOtherNueCC(event, true_apothem):
                otherNueCCHistos[2].Fill(reco_E_nu)
                otherNueCCHistos[8].Fill(event.prong_part_score[0])
            elif isTrueNumuCC(event, true_apothem):
                numuCCHistos[2].Fill(reco_E_nu)
                numuCCHistos[8].Fill(event.prong_part_score[0])
            elif isTrueNCpi0(event, true_apothem):
                NCpi0Histos[2].Fill(reco_E_nu)
                NCpi0Histos[8].Fill(event.prong_part_score[0])
            else:  #category for any other background not already covered
                if event.mc_vtx[2] > 5980 and event.mc_vtx[2] < 8422 and true_apothem < 850 and abs(event.mc_incoming)==14 and event.mc_current==1:
                    print "looks like i didn't cover all of the other nu_e-CC events in other nu_e-CC :(((("
                otherBackgroundHistos[2].Fill(reco_E_nu)
                otherBackgroundHistos[8].Fill(event.prong_part_score[0])


        #POST_SHOWER_SCORE_CUT and also PRE_GAP_CUT
        #if passesPreCuts(event, apothem) and E_avail < 1.5 and event.HasNoBackExitingTracks==1 and event.prong_part_score[0] > 0.7:
            #reco_E_nu_selected_POST_SHOWER_SCORE_CUT.Fill(reco_E_nu)
            #EMShower_score_selected_POST_SHOWER_SCORE_CUT.Fill(event.prong_part_score[0])
            #transverse_gap_score_selected_PRE_GAP_CUT.Fill(event.prong_TransverseGapScore[0])
            #if isTrueSignal(event, true_apothem):
                #reco_E_nu_signal_selected_POST_SHOWER_SCORE_CUT.Fill(reco_E_nu)
                #EMShower_score_signal_selected_POST_SHOWER_SCORE_CUT.Fill(event.prong_part_score[0])
                #transverse_gap_score_signal_selected_PRE_GAP_CUT.Fill(event.prong_TransverseGapScore[0])

        if passesPreCuts(event, apothem) and reco_E_avail < 1.5 and event.HasNoBackExitingTracks==1 and event.prong_part_score[0] > 0.7:
            selectedHistos[3].Fill(reco_E_nu)
            selectedHistos[9].Fill(event.prong_part_score[0])
            selectedHistos[10].Fill(event.prong_TransverseGapScore[0])
            if isTrueSignalWithProton(event, true_apothem):
                signalWithProtonHistos[3].Fill(reco_E_nu)
                signalWithProtonHistos[9].Fill(event.prong_part_score[0])
                signalWithProtonHistos[10].Fill(event.prong_TransverseGapScore[0])
            elif isTrueSignalWithoutProton(event, true_apothem):
                signalWithOutProtonHistos[3].Fill(reco_E_nu)
                signalWithOutProtonHistos[9].Fill(event.prong_part_score[0])
                signalWithOutProtonHistos[10].Fill(event.prong_TransverseGapScore[0])
            elif isTrueOtherNueCC(event, true_apothem):
                otherNueCCHistos[3].Fill(reco_E_nu)
                otherNueCCHistos[9].Fill(event.prong_part_score[0])
                otherNueCCHistos[10].Fill(event.prong_TransverseGapScore[0])
            elif isTrueNumuCC(event, true_apothem):
                numuCCHistos[3].Fill(reco_E_nu)
                numuCCHistos[9].Fill(event.prong_part_score[0])
                numuCCHistos[10].Fill(event.prong_TransverseGapScore[0])
            elif isTrueNCpi0(event, true_apothem):
                NCpi0Histos[3].Fill(reco_E_nu)
                NCpi0Histos[9].Fill(event.prong_part_score[0])
                NCpi0Histos[10].Fill(event.prong_TransverseGapScore[0])
            else:  #category for any other background not already covered
                if event.mc_vtx[2] > 5980 and event.mc_vtx[2] < 8422 and true_apothem < 850 and abs(event.mc_incoming)==14 and event.mc_current==1:
                    print "looks like i didn't cover all of the other nu_e-CC events in other nu_e-CC :(((("
                otherBackgroundHistos[3].Fill(reco_E_nu)
                otherBackgroundHistos[9].Fill(event.prong_part_score[0])
                otherBackgroundHistos[10].Fill(event.prong_TransverseGapScore[0])


        #POST_GAP_CUT and also PRE_FRONTdEdX_CUT
        #if passesPreCuts(event, apothem) and E_avail < 1.5 and event.HasNoBackExitingTracks==1 and event.prong_part_score[0] > 0.7 and event.prong_TransverseGapScore[0] > 15:
            #reco_E_nu_selected_POST_TRANSVERSE_GAP_CUT.Fill(reco_E_nu)
            #transverse_gap_score_selected_POST_GAP_CUT.Fill(event.prong_TransverseGapScore[0])
            #frontdEdX_selected_PRE_FRONTdEdX_CUT.Fill(event.prong_dEdXMeanFrontTracker[0])
            #if isTrueSignal(event, true_apothem):
                #reco_E_nu_signal_selected_POST_TRANSVERSE_GAP_CUT.Fill(reco_E_nu)
                #transverse_gap_score_signal_selected_POST_GAP_CUT.Fill(event.prong_TransverseGapScore[0])
                #frontdEdX_signal_selected_PRE_FRONTdEdX_CUT.Fill(event.prong_dEdXMeanFrontTracker[0])

        if passesPreCuts(event, apothem) and reco_E_avail < 1.5 and event.HasNoBackExitingTracks==1 and event.prong_part_score[0] > 0.7 and event.prong_TransverseGapScore[0] > 15:
            selectedHistos[4].Fill(reco_E_nu)
            selectedHistos[12].Fill(event.prong_dEdXMeanFrontTracker[0])
            selectedHistos[11].Fill(event.prong_TransverseGapScore[0])
            if isTrueSignalWithProton(event, true_apothem):
                signalWithProtonHistos[4].Fill(reco_E_nu)
                signalWithProtonHistos[12].Fill(event.prong_dEdXMeanFrontTracker[0])
                signalWithProtonHistos[11].Fill(event.prong_TransverseGapScore[0])
            elif isTrueSignalWithoutProton(event, true_apothem):
                signalWithOutProtonHistos[4].Fill(reco_E_nu)
                signalWithOutProtonHistos[12].Fill(event.prong_dEdXMeanFrontTracker[0])
                signalWithOutProtonHistos[11].Fill(event.prong_TransverseGapScore[0])
            elif isTrueOtherNueCC(event, true_apothem):
                otherNueCCHistos[4].Fill(reco_E_nu)
                otherNueCCHistos[12].Fill(event.prong_dEdXMeanFrontTracker[0])
                otherNueCCHistos[11].Fill(event.prong_TransverseGapScore[0])
            elif isTrueNumuCC(event, true_apothem):
                numuCCHistos[4].Fill(reco_E_nu)
                numuCCHistos[12].Fill(event.prong_dEdXMeanFrontTracker[0])
                numuCCHistos[11].Fill(event.prong_TransverseGapScore[0])
            elif isTrueNCpi0(event, true_apothem):
                NCpi0Histos[4].Fill(reco_E_nu)
                NCpi0Histos[12].Fill(event.prong_dEdXMeanFrontTracker[0])
                NCpi0Histos[11].Fill(event.prong_TransverseGapScore[0])
            else:  #category for any other background not already covered
                if event.mc_vtx[2] > 5980 and event.mc_vtx[2] < 8422 and true_apothem < 850 and abs(event.mc_incoming)==14 and event.mc_current==1:
                    print "looks like i didn't cover all of the other nu_e-CC events in other nu_e-CC :(((("
                otherBackgroundHistos[4].Fill(reco_E_nu)
                otherBackgroundHistos[12].Fill(event.prong_dEdXMeanFrontTracker[0])
                otherBackgroundHistos[11].Fill(event.prong_TransverseGapScore[0])


        #POST_FRONTdEdX_CUT
        #if passesPreCuts(event, apothem) and E_avail < 1.5 and event.HasNoBackExitingTracks==1 and event.prong_part_score[0] > 0.7 and event.prong_TransverseGapScore[0] > 15 and event.prong_dEdXMeanFrontTracker[0] < 2.4:
            #reco_E_nu_selected_POST_FRONTdEdX_CUT.Fill(reco_E_nu)
            #frontdEdX_selected_POST_FRONTdEdX_CUT.Fill(event.prong_dEdXMeanFrontTracker[0])
            #if isTrueSignal(event, true_apothem):
                #reco_E_nu_signal_selected_POST_FRONTdEdX_CUT.Fill(reco_E_nu)
                #frontdEdX_signal_selected_POST_FRONTdEdX_CUT.Fill(event.prong_dEdXMeanFrontTracker[0])

        if passesPreCuts(event, apothem) and reco_E_avail < 1.5 and event.HasNoBackExitingTracks==1 and event.prong_part_score[0] > 0.7 and event.prong_TransverseGapScore[0] > 15 and event.prong_dEdXMeanFrontTracker[0] < 2.4:
            selectedHistos[5].Fill(reco_E_nu)
            selectedHistos[13].Fill(event.prong_dEdXMeanFrontTracker[0])
            if isTrueSignalWithProton(event, true_apothem):
                signalWithProtonHistos[5].Fill(reco_E_nu)
                signalWithProtonHistos[13].Fill(event.prong_dEdXMeanFrontTracker[0])
            elif isTrueSignalWithoutProton(event, true_apothem):
                signalWithOutProtonHistos[5].Fill(reco_E_nu)
                signalWithOutProtonHistos[13].Fill(event.prong_dEdXMeanFrontTracker[0])
            elif isTrueOtherNueCC(event, true_apothem):
                otherNueCCHistos[5].Fill(reco_E_nu)
                otherNueCCHistos[13].Fill(event.prong_dEdXMeanFrontTracker[0])
            elif isTrueNumuCC(event, true_apothem):
                numuCCHistos[5].Fill(reco_E_nu)
                numuCCHistos[13].Fill(event.prong_dEdXMeanFrontTracker[0])
            elif isTrueNCpi0(event, true_apothem):
                NCpi0Histos[5].Fill(reco_E_nu)
                NCpi0Histos[13].Fill(event.prong_dEdXMeanFrontTracker[0])
            else:  #category for any other background not already covered
                if event.mc_vtx[2] > 5980 and event.mc_vtx[2] < 8422 and true_apothem < 850 and abs(event.mc_incoming)==14 and event.mc_current==1:
                    print "looks like i didn't cover all of the other nu_e-CC events in other nu_e-CC :(((("
                otherBackgroundHistos[5].Fill(reco_E_nu)
                otherBackgroundHistos[13].Fill(event.prong_dEdXMeanFrontTracker[0])


    #k now that all that nonsense is over with, going to save these to a root file for future use
    outFile = ROOT.TFile('CutSummaryHistos.root', 'RECREATE')

    for i in range(len(selectedHistoNames)): 
        selectedHistos[i].Write()
        signalWithProtonHistos[i].Write()
        signalWithOutProtonHistos[i].Write()
        otherNueCCHistos[i].Write()
        numuCCHistos[i].Write()
        NCpi0Histos[i].Write()
        otherBackgroundHistos[i].Write()
    
    outFile.Close()

#names in the root file of the histograms we're interested, at the moment only affects plotting but i should change to also have the event loop use this 
#itll cut down the number of lines by a ton
selectedHistoNames = ["reco_E_nu_selected_PRE_EAVAIL_CUT", "reco_E_nu_selected_POST_EAVAIL_CUT", "reco_E_nu_selected_POST_BACK_TRACK_REJECTION", "reco_E_nu_selected_POST_SHOWER_SCORE_CUT", "reco_E_nu_selected_POST_TRANSVERSE_GAP_CUT", "reco_E_nu_selected_POST_FRONTdEdX_CUT", "reco_E_avail_selected_PRE_EAVAIL_CUT", "reco_E_avail_selected_POST_EAVAIL_CUT", "EMShower_score_selected_PRE_SHOWER_SCORE_CUT", "EMShower_score_selected_POST_SHOWER_SCORE_CUT","transverse_gap_score_selected_PRE_GAP_CUT", "transverse_gap_score_selected_POST_GAP_CUT","frontdEdX_selected_PRE_FRONTdEdX_CUT", "frontdEdX_selected_POST_FRONTdEdX_CUT"]

signalWithProtonHistoNames=[]
signalWithOutProtonHistoNames=[]
otherNueCCHistoNames=[]
numuCCHistoNames=[]
NCpi0HistoNames=[]
otherBackgroundHistoNames=[]

for i in selectedHistoNames:
    signalWithProtonHistoNames.append(i.replace("selected", "mesonless_nue_CC_with_proton"))
    signalWithOutProtonHistoNames.append(i.replace("selected", "mesonless_nue_CC_without_proton"))
    otherNueCCHistoNames.append(i.replace("selected", "other_nue_CC"))
    numuCCHistoNames.append(i.replace("selected", "numu_CC"))
    NCpi0HistoNames.append(i.replace("selected", "NC_pi0"))
    otherBackgroundHistoNames.append(i.replace("selected", "other_background"))


EventLoopToMakeHistos(selectedHistoNames)

#TH1.AddDirectory(False)
print "done writing"
#dataFile = TFile.Open(sys.argv[1])
#mcFile = TFile.Open(sys.argv[2])
mcFile = ROOT.TFile.Open("CutSummaryHistos.root")
plotter = PlotUtils.MnvPlotter()
#plotter.ApplyStyle(PlotUtils.kCCNuPionIncStyle)

#mcPOT = mcFile.Get("POTUsed").GetVal()
#dataPOT = dataFile.Get("POTUsed").GetVal()

#mcScale = dataPOT/mcPOT

#I have no idea what I'm doing
#plotter.DrawDataMC(dataHist, mcHist, mcScale, "TR")

#I need to pass this guy a TObjArray of two histograms -> one that is pTmu_selected_signal_reco, and one that is all of the backgrounds
#for all the backgrounds I can either add up the 3 background histograms (pTmu_background_Wrongsign+pTmu_background_NC+pTmu_background_other)
#or I can do all minus signal (pTmu_data - pTmu_selected_signal_reco)

for i in range(len(selectedHistoNames)):
    signalWithProton = mcFile.Get(signalWithProtonHistoNames[i])
    signalWithoutProton = mcFile.Get(signalWithOutProtonHistoNames[i])
    otherNueCC = mcFile.Get(otherNueCCHistoNames[i])
    numuCC = mcFile.Get(numuCCHistoNames[i])
    NCPi0 = mcFile.Get(NCpi0HistoNames[i])
    otherBkgd = mcFile.Get(otherBackgroundHistoNames[i])

    signalWithProton.SetTitle('mesonless nue-CC w proton')
    signalWithoutProton.SetTitle('mesonless nue-CC w/o proton')
    otherNueCC.SetTitle('other nue-CC')
    numuCC.SetTitle('numu-CC')
    NCPi0.SetTitle('NC with pi0')
    otherBkgd.SetTitle('other')

    #cause my dumb ass made em TH1D's instead of MnvH1D's
    signalWithProton = PlotUtils.MnvH1D(signalWithProton)
    signalWithoutProton = PlotUtils.MnvH1D(signalWithoutProton)
    otherNueCC = PlotUtils.MnvH1D(otherNueCC)
    numuCC = PlotUtils.MnvH1D(numuCC)
    NCPi0 = PlotUtils.MnvH1D(NCPi0)
    otherBkgd = PlotUtils.MnvH1D(otherBkgd)

    array = ROOT.TObjArray()
    array.Add(otherBkgd)
    array.Add(NCPi0)
    array.Add(numuCC)
    array.Add(otherNueCC)
    array.Add(signalWithoutProton)
    array.Add(signalWithProton)
        
    outName=selectedHistoNames[i].replace("_selected_", "_")
    canvas = ROOT.TCanvas( 'canvas', outName, 0, 0, 2000, 1600 )

    #to get variable name for the x axis
    head, sep, tail = selectedHistoNames[i].partition('_selected')

    #arguments for stackedMC array are:
    #DrawStackedMC(mcHists, mcScale, legend position, base color, color offset, fill style, xaxislabel, yaxislabel)
    plotter.DrawStackedMC(array, 1.0, "TR", 2, 1, 3001, head, "N events")
    
    outName = outName + ".png"
    canvas.SaveAs(outName)
    canvas.Delete()
#plotter.WritePreliminary()


