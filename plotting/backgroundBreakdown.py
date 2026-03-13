#!/usr/bin/python

#USAGE: backgroundStack.py <dataFile.root> <mcFile.root>

#in progress edit to plotCuts.py trying to make it so I dont have to manually initialize like 8000 histograms, currently unfinished im on like line 99
#tbh might not even finish it cause I should move back to the actual MAT event loop at some point. 
#from ROOT import *


import ROOT
from ROOT import PlotUtils

import sys
from array import array
import math
import ctypes

#to make it not display the canvases as it draws and saves them, saves a bunch of time
ROOT.gROOT.SetBatch(True)

#drawData = False
drawData = True

varNames = ["E_lep", "E_avail", "E_nu", "Lepton_Pt", "Lepton_Pl", "Theta_lep", "Proton_p", "Proton_Pt", "Theta_p", "Proton_T","DeltaPt","DeltaPtX","DeltaPtY","DeltaPl","P_n","AlphaPt","PhiPt"]
#signalHistoNames = ["Lepton_Pt"
#varNames = ["E_lep_MichelSB", "E_avail_MichelSB", "E_nu_MichelSB", "Lepton_Pt_MichelSB", "Lepton_Pl_MichelSB", "Theta_lep_MichelSB", "Proton_p_MichelSB", "Proton_Pt_MichelSB", "Theta_p_MichelSB", "Proton_T_MichelSB","DeltaPt_MichelSB","DeltaPtX_MichelSB","DeltaPtY_MichelSB","DeltaPl_MichelSB","P_n_MichelSB","AlphaPt_MichelSB","PhiPt_MichelSB"]
#varNames = ["E_lep_MeanFrontDEDXSB", "E_avail_MeanFrontDEDXSB", "E_nu_MeanFrontDEDXSB", "Lepton_Pt_MeanFrontDEDXSB", "Lepton_Pl_MeanFrontDEDXSB", "Theta_lep_MeanFrontDEDXSB", "Proton_p_MeanFrontDEDXSB", "Proton_Pt_MeanFrontDEDXSB", "Theta_p_MeanFrontDEDXSB", "Proton_T_MeanFrontDEDXSB","DeltaPt_MeanFrontDEDXSB","DeltaPtX_MeanFrontDEDXSB","DeltaPtY_MeanFrontDEDXSB","DeltaPl_MeanFrontDEDXSB","P_n_MeanFrontDEDXSB","AlphaPt_MeanFrontDEDXSB","PhiPt_MeanFrontDEDXSB"]

signalHistoNames = [i + "_selected_signal_reco" for i in varNames]
NueCC0PiHistoNames = [i + "_background_NuECC_with_pions" for i in varNames]
otherNueCCHistoNames = [i + "_background_Other_NueCC" for i in varNames]
NCPi0HistoNames = [i + "_background_NC_pi0" for i in varNames]
CCnumuPi0HistoNames = [i + "_background_CC_Numu_pi0" for i in varNames]
otherBackgroundHistoNames = [i + "_background_Other" for i in varNames]
dataHistoNames = [i + "_data" for i in varNames]

#print(varNames)

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

    print(varNames[i])

    signal.SetTitle('signal (nu_e QELike + proton)')
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
    
    outName = varNames[i]
    canvas = ROOT.TCanvas( 'canvas', "test", 0, 0, 2000, 1600 )

    #some plotter options
    #arguments for stackedMC array are:
    #DrawStackedMC(mcHists, mcScale, legend position, base color, color offset, fill style, xaxislabel, yaxislabel)

    #this lets me fix my y axis scale if I'm trying to make comparison plots for something
    #plotter.axis_maximum = 1000;

    if drawData:
        data = dataFile.Get(dataHistoNames[i])
        #data = PlotUtils.MnvH1D(data)
        #print(dataHistoNames[i])
        data.SetTitle('data')
        plotter.DrawDataStackedMC(data, array, arr, mcScale, "TR", "Data", 1001, signal.GetXaxis().GetTitle(), "N events")
        #plotter.DrawDataStackedMC(data, array, arr, mcScale, "TR", "Data", 1001, head, "N events")
        plotter.WriteNorm("Normalized to ", "TL", 0.03, -0.05, 0, dataPOT)

    else:
        #data = dataFile.Get(dataHistoNames[i])
        data = mcFile.Get(dataHistoNames[i])
        bins = data.GetXaxis()
        #bins = signal.GetXaxis()
        for i in range(bins.GetNbins()):
            data.SetBinContent(i, 0)
        #data = ROOT.TH1D("test","test", signal.GetNbinsX(), 
        #plotter.DrawStackedMC(array, 1.0, "TR", 2, 1, 3001, head, "N events")
        plotter.DrawDataStackedMC(data, array, arr, mcScale, "TR", "Data", 1001, signal.GetXaxis().GetTitle(), "N events")

    plotter.WritePreliminary("TL", 0.03, 0, 0, True)
    outName = outName + ".png"
    canvas.SaveAs(outName)
    canvas.Delete()


