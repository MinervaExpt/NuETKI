#!/usr/bin/python

#USAGE: plotMigration.py <mcFile.root> <varName>

# quick python script to make a row normalized migration matrix for the given var and file

import ROOT
from ROOT import PlotUtils

import sys
from array import array
import math
import ctypes

ROOT.gROOT.SetBatch(True)


fName = sys.argv[1]
var = sys.argv[2]
f = ROOT.TFile.Open(fName)

matrix = f.Get(var+"_migration")
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
canvas.SaveAs(var+"_migration.png")
canvas.Delete()
