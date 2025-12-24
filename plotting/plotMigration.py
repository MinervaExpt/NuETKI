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


signalHistoNames = ["Psi_selected_signal_reco"]



#mcFile = ROOT.TFile.Open("/pnfs/minerva/persistent/users/cpernas/default_analysis_loc/MC_Oct_25_2024_Psi.root")
f = ROOT.TFile.Open("/exp/minerva/data/users/cpernas/NuE_TKI/Nov_07_per_universe_scaling/DeltaPt/maxReg/scaled_mc.root")

matrix = f.Get("DeltaPt_migration")
plotter = PlotUtils.MnvPlotter()

#plotter.ApplyStyle(PlotUtils.kCCNuPionIncStyle)
# plotter.legend_text_size = 0.02

# plotter.data_line_width = 0
# plotter.data_marker_size = 0

    #  void DrawNormalizedMigrationHistogram(
    #         const TH2D* h_migration,
    #         const bool drawAsMatrix = false,
    #         const bool coarseContours = false,
    #         const bool includeFlows = true,
    #         const bool no_text = false
    #         );
    
canvas = ROOT.TCanvas( 'canvas', 'migration', 0, 0, 800, 600 )
plotter.DrawNormalizedMigrationHistogram(matrix)
# outName = outName + "_signal.png"
canvas.SaveAs("DeltaPt_migration.png")
canvas.Delete()

#python to cpp nonsense idk
# mcColors = [416]
#kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
#kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
#kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900

