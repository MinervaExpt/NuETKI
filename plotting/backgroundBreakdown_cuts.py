#!/usr/bin/python

#USAGE: backgroundStack.py <dataFile.root> <mcFile.root>

#in progress edit to plotCuts.py trying to make it so I dont have to manually initialize like 8000 histograms, currently unfinished im on like line 99
#tbh might not even finish it cause I should move back to the actual MAT event loop at some point. 
#from ROOT import *


import ROOT
from ROOT import PlotUtils
from array import array
import math
import ctypes
import sys

#to make it not display the canvases as it draws and saves them, saves a bunch of time
ROOT.gROOT.SetBatch(True)

drawData = False
#drawData = True

#cutsToDraw = ["NoVertexMismatch","VertexZ","InApothem","StartPointVertexMultiplicity","Afterpulsing","Deadtime","NoBackExitingTracks","DSCalVisE","ODCalVisE","VertexTrackMultiplicity","TransverseGapScore","NonMIPClusFrac","EMScore","NMichels","MeanFrontDEDX","E_lep","Modified_E_avail","ESCChi2","Psi", "ProtonMomentum","ProtonTheta","LeptonPt","Etheta"]
#cutsToDraw = ["Afterpulsing", "DSCalVisE", "E_lep", "EMLikeTrackScore", "ESC", "MeanFrontdEdX", "MichelCut", "ModifiedEavailable", "NIsoBlobs", "NonMIPClusFrac", "ODCalVisE", "StartPointVertexMultiplicity", "TransverseGapScore", "VertexTrackMultiplicity", "Psi", "Lepton_Pt", "Etheta", "ProtonMomentum", "ProtonTheta"]
cutsToDraw = ["ESC", "MeanFrontdEdX", "ModifiedEavailable"]

for cut in cutsToDraw:

    #names in the root file of the histograms we're interested in, just put the signal and it'll find the other bkgd categories from that
    signalHistoNames = [cut+"_selected_signal_reco"] 
    print(signalHistoNames)

    #can't think of a better way to do this... I have to put the cut arrows in a different spot for each plot, so here we go:
    
    if cut=="NoVertexMismatch":
        cutArrowLoc = 1 #where the cut is on x axis
        cutArrowDirection = "R" #"L" or "R"
        cutArrowHeight = 45 #height on the histogram (in terms of y axis units, aka n events)
        cutArrowLength = 0.2 #cut arrow horizontal length on the part that sticks out, in terms of x axis units
    if cut=="VertexZ":
        cutArrowLoc = 5980 #where the cut is on x axis
        cutArrowDirection = "R" #"L" or "R"
        cutArrowHeight = 98 #height on the histogram (in terms of y axis units, aka n events)
        cutArrowLength = 200 #cut arrow horizontal length on the part that sticks out, in terms of x axis units
    if cut=="InApothem":
        cutArrowLoc = 1 #where the cut is on x axis
        cutArrowDirection = "R" #"L" or "R"
        cutArrowHeight = 600 #height on the histogram (in terms of y axis units, aka n events)
        cutArrowLength = 0.2 #cut arrow horizontal length on the part that sticks out, in terms of x axis units
    if cut=="StartPointVertexMultiplicity":
        cutArrowLoc = 1 #where the cut is on x axis
        cutArrowDirection = "L" #"L" or "R"
        cutArrowHeight = 600 #height on the histogram (in terms of y axis units, aka n events)
        cutArrowLength = 0.5 #cut arrow horizontal length on the part that sticks out, in terms of x axis units
    if cut=="ESC":
        cutArrowLoc = 20 #where the cut is on x axis
        cutArrowDirection = "L" #"L" or "R"
        cutArrowHeight = 200 #height on the histogram (in terms of y axis units, aka n events)
        cutArrowLength = 2 #cut arrow horizontal length on the part that sticks out, in terms of x axis units
    if cut=="MeanFrontdEdX":
        cutArrowLoc = 2.4 #where the cut is on x axis
        cutArrowDirection = "L" #"L" or "R"
        cutArrowHeight = 150 #height on the histogram (in terms of y axis units, aka n events)
        cutArrowLength = 0.3 #cut arrow horizontal length on the part that sticks out, in terms of x axis units
    if cut=="ModifiedEavailable":
        cutArrowLoc = 600 #where the cut is on x axis
        cutArrowDirection = "L" #"L" or "R"
        cutArrowHeight = 150 #height on the histogram (in terms of y axis units, aka n events)
        cutArrowLength = 50 #cut arrow horizontal length on the part that sticks out, in terms of x axis units    
    if cut=="EMLikeTrackScore":
        cutArrowLoc = 0.7 #where the cut is on x axis
        cutArrowDirection = "R" #"L" or "R"
        cutArrowHeight = 100 #height on the histogram (in terms of y axis units, aka n events)
        cutArrowLength = 0.1 #cut arrow horizontal length on the part that sticks out, in terms of x axis units
    """
    if cut="NoVertexMismatch":
        cutArrowLoc = 1 #where the cut is on x axis
        cutArrowDirection = "R" #"L" or "R"
        cutArrowHeight = 45 #height on the histogram (in terms of y axis units, aka n events)
        cutArrowLength = 0.2 #cut arrow horizontal length on the part that sticks out, in terms of x axis units
    if cut="NoVertexMismatch":
        cutArrowLoc = 1 #where the cut is on x axis
        cutArrowDirection = "R" #"L" or "R"
        cutArrowHeight = 45 #height on the histogram (in terms of y axis units, aka n events)
        cutArrowLength = 0.2 #cut arrow horizontal length on the part that sticks out, in terms of x axis units
    if cut="NoVertexMismatch":
        cutArrowLoc = 1 #where the cut is on x axis
        cutArrowDirection = "R" #"L" or "R"
        cutArrowHeight = 45 #height on the histogram (in terms of y axis units, aka n events)
        cutArrowLength = 0.2 #cut arrow horizontal length on the part that sticks out, in terms of x axis units
    if cut="NoVertexMismatch":
        cutArrowLoc = 1 #where the cut is on x axis
        cutArrowDirection = "R" #"L" or "R"
        cutArrowHeight = 45 #height on the histogram (in terms of y axis units, aka n events)
        cutArrowLength = 0.2 #cut arrow horizontal length on the part that sticks out, in terms of x axis units
    if cut="NoVertexMismatch":
        cutArrowLoc = 1 #where the cut is on x axis
        cutArrowDirection = "R" #"L" or "R"
        cutArrowHeight = 45 #height on the histogram (in terms of y axis units, aka n events)
        cutArrowLength = 0.2 #cut arrow horizontal length on the part that sticks out, in terms of x axis units
    if cut="NoVertexMismatch":
        cutArrowLoc = 1 #where the cut is on x axis
        cutArrowDirection = "R" #"L" or "R"
        cutArrowHeight = 45 #height on the histogram (in terms of y axis units, aka n events)
        cutArrowLength = 0.2 #cut arrow horizontal length on the part that sticks out, in terms of x axis units
    if cut="NoVertexMismatch":
        cutArrowLoc = 1 #where the cut is on x axis
        cutArrowDirection = "R" #"L" or "R"
        cutArrowHeight = 45 #height on the histogram (in terms of y axis units, aka n events)
        cutArrowLength = 0.2 #cut arrow horizontal length on the part that sticks out, in terms of x axis units
    """

    
    NueCC0PiHistoNames=[]
    otherNueCCHistoNames=[]
    NCPi0HistoNames=[]
    CCnumuPi0HistoNames=[]
    otherBackgroundHistoNames=[]
    
    dataHistoNames=[]
    
    for i in signalHistoNames:
        NueCC0PiHistoNames.append(i.replace("selected_signal_reco", "background_NuECC_with_pions"))
        otherNueCCHistoNames.append(i.replace("selected_signal_reco", "background_Other_NueCC"))
        NCPi0HistoNames.append(i.replace("selected_signal_reco", "background_NC_pi0"))
        CCnumuPi0HistoNames.append(i.replace("selected_signal_reco", "background_CC_Numu_pi0"))
        otherBackgroundHistoNames.append(i.replace("selected_signal_reco", "background_Other"))
        
        dataHistoNames.append(i.replace("selected_signal_reco", "data"))
        

    #if drawData:
    #dataFile = ROOT.TFile.Open(sys.argv[1])
    #mcFile = ROOT.TFile.Open(sys.argv[2])
    if cut == "Psi" or cut == "Lepton_Pt" or cut == "Etheta" or cut == "ProtonMomentum" or cut == "ProtonTheta":
        dataFile = ROOT.TFile.Open("/exp/minerva/data/users/cpernas/NuE_TKI/N-1_Feb_22/ExtraCuts/Data_Feb_23.root")
        mcFile = ROOT.TFile.Open("/exp/minerva/data/users/cpernas/NuE_TKI/N-1_Feb_22/ExtraCuts/MC_Feb_23.root")
    else:        
        dataFile = ROOT.TFile.Open("/exp/minerva/data/users/cpernas/NuE_TKI/N-1_Feb_22/"+cut+"/MC_Feb_23.root")
        mcFile = ROOT.TFile.Open("/exp/minerva/data/users/cpernas/NuE_TKI/N-1_Feb_22/"+cut+"/MC_Feb_23.root")
    print(dataFile)
    print(mcFile)


    plotter = PlotUtils.MnvPlotter()
    #plotter.ApplyStyle(PlotUtils.kCCNuPionIncStyle)
    plotter.legend_text_size = 0.018
    plotter.legend_offset_x = 0.07
    if drawData:
        plotter.data_line_width = 2
        plotter.data_marker_size = 2
    else:
        plotter.data_line_width = 0
        plotter.data_marker_size = 0

    #python to cpp nonsense idk
    mcColors = [4, 7, 6, 2, 5, 416]
    #kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
    #kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
    #kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900

    
    arr =(ctypes.c_int * len(mcColors))(*mcColors)



    mcPOT = mcFile.Get("POTUsed").GetVal()
    if drawData:
        dataPOT = dataFile.Get("POTUsed").GetVal()
        mcScale = dataPOT/mcPOT
    else:
        mcScale = 1

    #I need to pass this guy a TObjArray of two histograms -> one that is pTmu_selected_signal_reco, and one that is all of the backgrounds
    #for all the backgrounds I can either add up the 3 background histograms (pTmu_background_Wrongsign+pTmu_background_NC+pTmu_background_other)
    #or I can do all minus signal (pTmu_data - pTmu_selected_signal_reco)

    for i in range(len(signalHistoNames)):
        signal = mcFile.Get(signalHistoNames[i])
        NueCC0Pi = mcFile.Get(NueCC0PiHistoNames[i])
        otherNueCC = mcFile.Get(otherNueCCHistoNames[i])
        NCPi0 = mcFile.Get(NCPi0HistoNames[i])
        CCnumuPi0 = mcFile.Get(CCnumuPi0HistoNames[i])
        otherBkgd = mcFile.Get(otherBackgroundHistoNames[i])

        signal.SetTitle('signal (nu_e QELike with proton)')
        NueCC0Pi.SetTitle('nu_e nonQE (has FS mesons)')
        otherNueCC.SetTitle('Other nu_eCC')
        NCPi0.SetTitle('NC with pi0')
        CCnumuPi0.SetTitle('nu_mu CC with pi0')
        otherBkgd.SetTitle('other')

        array = ROOT.TObjArray()
        array.Add(otherBkgd)
        array.Add(CCnumuPi0)
        array.Add(NCPi0)
        array.Add(otherNueCC)
        array.Add(NueCC0Pi)
        array.Add(signal)
        
        outName=signalHistoNames[i].replace("_selected_signal_reco", "")
        canvas = ROOT.TCanvas( 'canvas', outName, 0, 0, 2000, 1600 )

        
        #to get variable name for the x axis
        head, sep, tail = signalHistoNames[i].partition('_selected')            
        
        #else:
        data = dataFile.Get(dataHistoNames[i])
        bins = data.GetXaxis()
        
        for i in range(bins.GetNbins()):
            data.SetBinContent(i, 0)

        #plotter.axis_minimum=0.001
            
        #plotter.SetYRangeUser(0, 1000)
        #plotter.DrawDataStackedMC(data, array, arr, mcScale, "TR", "Data", 1001, head, "N events")
        if cut=="ModifiedEavailable":
            plotter.DrawDataStackedMC(data, array, arr, mcScale, "TR", "Data", 1001, "EAvail - Proton Energy [MeV]", "N events")
            
        else:               
            plotter.DrawDataStackedMC(data, array, arr, mcScale, "TR", "Data", 1001, signal.GetXaxis().GetTitle(), "N events")

        
        #void MnvPlotter::AddCutArrow(const double cut_location, const double y1, const double y2, const double arrow_length, const std::string& arrow_direction)
        #i think y1 is bottom of line, y2 is top? in the y axis, not pixels or anything. So if my stuff goes up to 500 events, do 500 or smth
        #arrow length is the horizontal part sticking out
        plotter.arrow_line_color = 415
        plotter.AddCutArrow(cutArrowLoc, 0, cutArrowHeight, cutArrowLength, cutArrowDirection)
        #plotter.AddCutArrow(cutArrowLoc, 0, cutArrowHeight, cutArrowLength, cutArrowDirection) #some cuts will need two cut arrows, if i cut both low and high

        #canvas.SetLogy()
        plotter.WritePreliminary("TC", 0.035, 0.01, -0.14, True)


        outName = outName + ".png"
        canvas.SaveAs(outName)
        canvas.Delete()

    dataFile.Close()
    mcFile.Close()
