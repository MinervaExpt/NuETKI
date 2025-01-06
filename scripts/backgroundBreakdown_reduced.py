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


#drawData = False
drawData = True

#names in the root file of the histograms we're interested in, just put the signal and it'll find the other bkgd categories from that
signalHistoNames = ["E_nu_selected_signal_reco", "E_avail_selected_signal_reco", "Lepton_Pt_selected_signal_reco", "E_lep_selected_signal_reco", "Theta_lep_selected_signal_reco", "RecoProtonP_selected_signal_reco", "RecoProtonTheta_selected_signal_reco"]
#signalHistoNames = ["ProtonP_selected_signal_reco", "Theta_p_selected_signal_reco"]


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

#print signalHistoNames
#print withoutProtonHistoNames
#print nueelHistoNames
#print NCCohHistoNames
#print NCPi0HistoNames
#print otherBackgroundHistoNames

#if drawData:
dataFile = ROOT.TFile.Open(sys.argv[1])
mcFile = ROOT.TFile.Open(sys.argv[2])

plotter = PlotUtils.MnvPlotter()

#plotter.ApplyStyle(PlotUtils.kCCNuPionIncStyle)
plotter.legend_text_size = 0.015
if drawData:
    plotter.data_line_width = 2
    plotter.data_marker_size = 2
else:
    plotter.data_line_width = 0
    plotter.data_marker_size = 0

#python to cpp nonsense idk
mcColors = [4, 7, 6, 2, 5, 416]
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

    signal.SetTitle('signal (nu_e QELike + proton)')
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
    #NCPi0 = PlotUtils.MnvH1D(NCPi0)
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
    
    if drawData:
        data = dataFile.Get(dataHistoNames[i])
        data.SetTitle('data')
        plotter.DrawDataStackedMC(data, array, arr, mcScale, "TR", "Data", 1001, head, "N events")
        #normalization = "Normalized to 8.98e+19 data POT"
        #plotter.AddPlotLabel(normalization, 0.2, 0.95, 0.03)

    else:
        data = dataFile.Get(dataHistoNames[i])
        bins = data.GetXaxis()
        for i in range(bins.GetNbins()):
            data.SetBinContent(i, 0)
        #data = ROOT.TH1D("test","test", signal.GetNbinsX(), 
        #plotter.DrawStackedMC(array, 1.0, "TR", 2, 1, 3001, head, "N events")
        plotter.DrawDataStackedMC(data, array, arr, mcScale, "TR", "Data", 1001, head, "N events")

    
    outName = outName + ".png"
    canvas.SaveAs(outName)
    canvas.Delete()
#plotter.WritePreliminary()


