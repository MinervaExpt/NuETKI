#!/usr/bin/python

#USAGE: backgroundBreakdown.py <dataFile.root> <mcFile.root>

#from ROOT import *

#TO DO - Add (either in this script or as a new script) data/MC ratios on the bottom, the way other minerva ppl do. Is that a function in MnvPlotter? gotta check

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

varNames1 = ["E_lep", "E_avail", "E_nu", "Lepton_Pt", "Lepton_Pl", "Theta_lep", "Proton_p", "Proton_Pt", "Theta_p", "Proton_T","DeltaPt","DeltaPtX","DeltaPtY","DeltaPl","P_n","AlphaPt","PhiPt"]
varNames2 = ["E_lep_MichelSB", "E_avail_MichelSB", "E_nu_MichelSB", "Lepton_Pt_MichelSB", "Lepton_Pl_MichelSB", "Theta_lep_MichelSB", "Proton_p_MichelSB", "Proton_Pt_MichelSB", "Theta_p_MichelSB", "Proton_T_MichelSB","DeltaPt_MichelSB","DeltaPtX_MichelSB","DeltaPtY_MichelSB","DeltaPl_MichelSB","P_n_MichelSB","AlphaPt_MichelSB","PhiPt_MichelSB"]
varNames3 = ["E_lep_MeanFrontDEDXSB", "E_avail_MeanFrontDEDXSB", "E_nu_MeanFrontDEDXSB", "Lepton_Pt_MeanFrontDEDXSB", "Lepton_Pl_MeanFrontDEDXSB", "Theta_lep_MeanFrontDEDXSB", "Proton_p_MeanFrontDEDXSB", "Proton_Pt_MeanFrontDEDXSB", "Theta_p_MeanFrontDEDXSB", "Proton_T_MeanFrontDEDXSB","DeltaPt_MeanFrontDEDXSB","DeltaPtX_MeanFrontDEDXSB","DeltaPtY_MeanFrontDEDXSB","DeltaPl_MeanFrontDEDXSB","P_n_MeanFrontDEDXSB","AlphaPt_MeanFrontDEDXSB","PhiPt_MeanFrontDEDXSB"]

varNames = ["DeltaPt"];
#varNames = ["Lepton_Pt", "Lepton_Pt_MichelSB", "Lepton_Pt_MeanFrontDEDXSB"]
#varNames = varNames1 + varNames2 + varNames3

# (name in root file, legend label, map key, plotting color)
#CATEGORIES = [
#    ("_selected_signal_reco",        "signal (nu_e CCQELike + proton)", "all_signal"), 
#    ("_background_NuECC_with_pions", "nu_e nonQELike (has FS mesons)",  "NueCC0Pi"),
#    ("_background_Other_NueCC",      "Other nu_eCC",                    "otherNueCC"),
#    ("_background_NC_pi0",           "NC with pi0",                     "NCPi0"),
#    ("_background_CC_Numu_pi0",      "nu_mu CC with pi0",               "CCnumuPi0"),
#    ("_background_Other",            "other",                           "otherBkgd"),
#]
CATEGORIES = [
    ("_signal_QE",                                  "Sig. + QE",                "signalQE",      ROOT.TColor.GetColor("#90EE90")), #all signal is different shades of green
    ("_signal_RES",                                 "Sig. + RES",               "signalRES",     ROOT.TColor.GetColor("#4CBB17")),
    ("_signal_DIS",                                 "Sig. + DIS",               "signal2p2h",    ROOT.TColor.GetColor("#008000")),
    ("_signal_2p2h",                                "Sig. + 2p2h",              "signalDIS",     ROOT.TColor.GetColor("#005500")),
    ("_signal_Other",                               "Sig. + Other",             "signalOther",   ROOT.TColor.GetColor("#CCFF00")),
    ("_background_NuECC_nonQELike_single_pi_plus",  "Single #pi^{+}",           "singlePiPlus",  ROOT.TColor.GetColor("#FFFF99")), #all nonQELike in different shades of yellow
    ("_background_NuECC_nonQELike_single_pi_minus", "Single #pi^{-}",           "singlePiMinus", ROOT.TColor.GetColor("#FFD700")), 
    ("_background_NuECC_nonQELike_single_pi_zero",  "Single #pi^{0}",           "singlePiZero",  ROOT.TColor.GetColor("#FFA500")),
    ("_background_NuECC_nonQELike_Npi",             "N#pi",                     "multiPi",       ROOT.TColor.GetColor("#FF8C00")),
    ("_background_Other_NueCC",                     "Other #nu_{e}CC",          "otherNueCC",    ROOT.TColor.GetColor("#FF0000")), 
    ("_background_NC_pi0",                          "NC with #pi^{0}",          "NCPi0",         ROOT.TColor.GetColor("#FF00FF")),
    ("_background_CC_Numu_pi0",                     "#nu_{#mu}CC with #pi^{0}", "CCnumuPi0",     ROOT.TColor.GetColor("#00FFFF")),
    ("_background_Other",                           "other",                    "otherBkgd",     ROOT.TColor.GetColor("#0000FF")),
]

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

#build my array of colors
#mcColors = [ROOT.kBlue, ROOT.kCyan, ROOT.kMagenta, ROOT.kRed, ROOT.kYellow, ROOT.kGreen]
mcColors = []
for _, _, _, color in reversed(CATEGORIES):
    mcColors.append(color)
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

for varName in varNames:
    signal = mcFile.Get(varName + "_selected_signal_reco")

    histos = {}
    for suffix, title, key, color in CATEGORIES:
        h = mcFile.Get(varName + suffix)
        h.SetTitle(title)
        histos[key] = h

    array = ROOT.TObjArray()
    for _, _, key, _ in reversed(CATEGORIES):
        array.Add(histos[key])
    
    canvas = ROOT.TCanvas( 'canvas', "test", 0, 0, 2000, 1600 )

    #some plotter options
    #arguments for stackedMC array are:
    #DrawStackedMC(mcHists, mcScale, legend position, base color, color offset, fill style, xaxislabel, yaxislabel)

    #this lets me fix my y axis scale if I'm trying to make comparison plots for something
    #plotter.axis_maximum = 1000;

    if drawData:
        data = dataFile.Get(varName + "_data")
        #data = PlotUtils.MnvH1D(data)
        data.SetTitle('data')
        plotter.DrawDataStackedMC(data, array, arr, mcScale, "TR", "Data", 1001, signal.GetXaxis().GetTitle(), "N events")
        #plotter.DrawDataStackedMC(data, array, arr, mcScale, "TR", "Data", 1001, head, "N events")
        plotter.WriteNorm("Normalized to ", "TL", 0.03, -0.05, 0, dataPOT)

    else:
        data = mcFile.Get(varName + "_data")
        bins = data.GetXaxis()
        #bins = signal.GetXaxis()
        for i in range(bins.GetNbins()):
            data.SetBinContent(i, 0)
        #data = ROOT.TH1D("test","test", signal.GetNbinsX(), 
        #plotter.DrawStackedMC(array, 1.0, "TR", 2, 1, 3001, head, "N events")
        plotter.DrawDataStackedMC(data, array, arr, mcScale, "TR", "Data", 1001, signal.GetXaxis().GetTitle(), "N events")

    plotter.WritePreliminary("TL", 0.03, 0, 0, True)
    outName = varName + ".png"
    canvas.SaveAs(outName)
    canvas.Delete()
