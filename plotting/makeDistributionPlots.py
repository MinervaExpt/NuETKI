#!/usr/bin/python

#USAGE: makeErrorSummaryPlots.py <file.root> <histToPlot>

import sys
import ROOT
from ROOT import PlotUtils
sys.argv.append( '-b' )
ROOT.gROOT.SetBatch(True)


#mcFile = ROOT.TFile.Open(sys.argv[1])
#hist = mcFile.Get(argv[2])
mcFile = ROOT.TFile.Open("/exp/minerva/data/users/cpernas/NuE_TKI/systematics_test_Oct_16/me1L/MC_Oct_16.root")
hist = mcFile.Get("DeltaPt_data")

plotter = PlotUtils.MnvPlotter()
plotter.ApplyStyle(PlotUtils.kCCQENuEStyle)
plotter.axis_maximum = 0.4;

canvas = ROOT.TCanvas( 'canvas', "test", 0, 0, 2000, 1600 )

plotter.DrawDataMC()

#Plotter functions
#DrawStackedMC()
#DrawDataStackedMC()
#DrawDataStackedMCWithErrorBand()
#DrawDataMCRatio()


#So i guess to get those neat plots with distributions on top and ratios under, I gotta make it manually???
