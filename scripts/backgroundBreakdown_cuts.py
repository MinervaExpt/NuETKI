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
import ctypes


drawData = False
#drawData = True

#cutsToDraw = ["NoVertexMismatch","VertexZ","InApothem","StartPointVertexMultiplicity","Afterpulsing","Deadtime","NoBackExitingTracks","DSCalVisE","ODCalVisE","VertexTrackMultiplicity","TransverseGapScore","NonMIPClusFrac","EMScore","NMichels","MeanFrontDEDX","E_lep","Modified_E_avail","ESCChi2","nonProton_extraCuts","proton_extraCuts"]

#cutsToDraw = ["NoVertexMismatch","VertexZ","InApothem","StartPointVertexMultiplicity","Afterpulsing","Deadtime","NoBackExitingTracks","DSCalVisE","ODCalVisE","VertexTrackMultiplicity","TransverseGapScore","NonMIPClusFrac","EMScore","NMichels","MeanFrontDEDX","E_lep","Modified_E_avail","ESCChi2","Psi", "ProtonMomentum","ProtonTheta","LeptonPt","Etheta"]
#cutsToDraw = ["ESCChi2","Psi", "ProtonMomentum","ProtonTheta","LeptonPt","Etheta"]
cutsToDraw = ["VertexTrackMultiplicity"]

for cut in cutsToDraw:

    #names in the root file of the histograms we're interested in, just put the signal and it'll find the other bkgd categories from that
    signalHistoNames = [cut+"_selected_signal_reco"] 
    print(signalHistoNames)

    #can't think of a better way to do this... I have to put the cut arrows in a different spot for each plot, so here we go:
    """
    if cut="NoVertexMismatch":
        cutArrowLoc = 1 #where the cut is on x axis
        cutArrowDirection = "R" #"L" or "R"
        cutArrowHeight = 45 #height on the histogram (in terms of y axis units, aka n events)
        cutArrowLength = 0.2 #cut arrow horizontal length on the part that sticks out, in terms of x axis units
    if cut="VertexZ":
        cutArrowLoc = 5980 #where the cut is on x axis
        cutArrowDirection = "R" #"L" or "R"
        cutArrowHeight = 98 #height on the histogram (in terms of y axis units, aka n events)
        cutArrowLength = 200 #cut arrow horizontal length on the part that sticks out, in terms of x axis units
    if cut="InApothem":
        cutArrowLoc = 1 #where the cut is on x axis
        cutArrowDirection = "R" #"L" or "R"
        cutArrowHeight = 600 #height on the histogram (in terms of y axis units, aka n events)
        cutArrowLength = 0.2 #cut arrow horizontal length on the part that sticks out, in terms of x axis units
    if cut="StartPointVertexMultiplicity":
        cutArrowLoc = 1 #where the cut is on x axis
        cutArrowDirection = "L" #"L" or "R"
        cutArrowHeight = 600 #height on the histogram (in terms of y axis units, aka n events)
        cutArrowLength = 0.5 #cut arrow horizontal length on the part that sticks out, in terms of x axis units
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
        

    #print(NueCC0PiHistoNames)
    #print(otherNueCCHistoNames)
    #print(NCPi0HistoNames)
    #print(CCnumuPi0HistoNames)
    #print(otherBackgroundHistoNames.append)
    #print(dataHistoNames)
    
    #if drawData:
    #dataFile = ROOT.TFile.Open(sys.argv[1])
    #mcFile = ROOT.TFile.Open(sys.argv[2])
    if cut == "Psi" or cut == "LeptonPt" or cut == "Etheta":
        dataFile = ROOT.TFile.Open("/pnfs/minerva/persistent/users/cpernas/default_analysis_loc/Data_May_01_2025_me1M_me1M_Modified_nonProton_extraCuts.root")
        mcFile = ROOT.TFile.Open("/pnfs/minerva/persistent/users/cpernas/default_analysis_loc/MC_May_01_2025_me1M_nonProton_extraCuts.root")
    elif cut == "ProtonMomentum" or cut == "ProtonTheta":
        dataFile = ROOT.TFile.Open("/pnfs/minerva/persistent/users/cpernas/default_analysis_loc/Data_May_01_2025_me1M_me1M_Modified_proton_extraCuts.root")
        mcFile = ROOT.TFile.Open("/pnfs/minerva/persistent/users/cpernas/default_analysis_loc/MC_May_01_2025_me1M_proton_extraCuts.root")
    else:        
        #dataFile = ROOT.TFile.Open("/pnfs/minerva/persistent/users/cpernas/default_analysis_loc/Data_May_01_2025_me1M_me1M_Modified_"+cut+".root")
        #mcFile = ROOT.TFile.Open("/pnfs/minerva/persistent/users/cpernas/default_analysis_loc/MC_May_01_2025_me1M_"+cut+".root")
        dataFile = ROOT.TFile.Open("/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/Data_June_12_2025_"+cut+".root")
        mcFile = ROOT.TFile.Open("/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/MC_June_12_2025_"+cut+".root")
#Data_May_01_2025_me1M_me1M_Modified_ESCChi2.root
    print(dataFile)
    print(mcFile)
    #dataFile = ROOT.TFile.Open("/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/playlists/Data_Jan_23_2025_EAvail_WITHCUT.root")
    #mcFile = ROOT.TFile.Open("/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/playlists/MC_Jan_23_2025_EAvail_WITHCUT.root")


    plotter = PlotUtils.MnvPlotter()

    #plotter.ApplyStyle(PlotUtils.kCCNuPionIncStyle)
    plotter.legend_text_size = 0.02
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

        #cause my dumb ass made em TH1D's instead of MnvH1D's
        #signalWithProton = PlotUtils.MnvH1D(signalWithProton)
        #withoutProton = PlotUtils.MnvH1D(withoutProton)
        #otherNueCC = PlotUtils.MnvH1D(otherNueCC)
        #NCCoh = PlotUtils.MnvH1D(NCCoh)
        #OtherNC = PlotUtils.MnvH1D(OtherNC)
        #otherBkgd = PlotUtils.MnvH1D(otherBkgd)

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
        
        #if head == "E_avail" or head == "E_lep" or head == "E_nu":
        #head+=" [GeV]"
        #elif head == "Lepton_Pt":
        #head+=" [GeV/c]"
        #elif head == "Theta_lep":
        #head+= " [deg]"

        #some plotter options

        #arguments for stackedMC array are:
        #DrawStackedMC(mcHists, mcScale, legend position, base color, color offset, fill style, xaxislabel, yaxislabel)
        
        #else:
        data = dataFile.Get(dataHistoNames[i])
        bins = data.GetXaxis()
        
        for i in range(bins.GetNbins()):
            data.SetBinContent(i, 0)

        #plotter.axis_minimum=0.001
            
        #data = ROOT.TH1D("test","test", signal.GetNbinsX(), 
        #plotter.DrawStackedMC(array, 1.0, "TR", 2, 1, 1001, head, "N events")

        #plotter.SetYRangeUser(0, 1000)
        #plotter.DrawDataStackedMC(data, array, arr, mcScale, "TR", "Data", 1001, head, "N events")
        plotter.DrawDataStackedMC(data, array, arr, mcScale, "TR", "Data", 1001, signal.GetXaxis().GetTitle(), "N events")

        
        #void MnvPlotter::AddCutArrow(const double cut_location, const double y1, const double y2, const double arrow_length, const std::string& arrow_direction)
        #i think y1 is bottom of line, y2 is top? in the y axis, not pixels or anything. So if my stuff goes up to 500 events, do 500 or smth
        #arrow length is the horizontal part sticking out
        #plotter.arrow_line_color = 415
        #plotter.AddCutArrow(cutArrowLoc, 0, cutArrowHeight, cutArrowLength, cutArrowDirection)
        #plotter.AddCutArrow(cutArrowLoc, 0, cutArrowHeight, cutArrowLength, cutArrowDirection) #some cuts will need two cut arrows, if i cut both low and high

    
        #canvas.SetLogy()

        outName = outName + ".png"
        canvas.SaveAs(outName)
        
        canvas.Delete()
        #plotter.WritePreliminary()

    dataFile.Close()
    mcFile.Close()
