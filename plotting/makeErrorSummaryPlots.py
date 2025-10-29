#!/usr/bin/python

#USAGE: makeErrorSummaryPlots.py <file.root> <histToPlot>

import sys
import ROOT
from ROOT import PlotUtils
sys.argv.append( '-b' )
ROOT.gROOT.SetBatch(True)


#mcFile = ROOT.TFile.Open(sys.argv[1])
#hist = mcFile.Get(argv[2])
mcFile = ROOT.TFile.Open("/exp/minerva/data/users/cpernas/NuE_TKI/Oct_17_syst_test/merged/mc.root")
hist = mcFile.Get("DeltaPt_data")

plotter = PlotUtils.MnvPlotter()
plotter.ApplyStyle(PlotUtils.kCCQENuEStyle)
plotter.axis_maximum = 0.4;

canvas = ROOT.TCanvas( 'canvas', "test", 0, 0, 2000, 1600 )

#DrawErrorSummary( hist, legPos, includeStat, solidLinesOnly, ignoreThreshold, covAreaNormalize, errorGroupName, asfrac, Ytitle, ignoreUngrouped, histDrawOption  

#Main categories
plotter.DrawErrorSummary(hist)
canvas.Print("uncertaintySummary.png")

#outName = outName + ".png"
#canvas.SaveAs(outName)
#canvas.Delete()

#Category breakdowns
#plotter.DrawErrorSummary(hist, "TR", true, true, 1e-5, false, "Cross Section Models")

#plotter.DrawErrorSummary(hist, "TR", true, true, 1e-5, false, "FSI Models")

#plotter.DrawErrorSummary(hist, "TR", true, true, 1e-5, false, "Minerva Tunes")

#plotter.DrawErrorSummary(hist, "TR", true, true, 1e-5, false, "Detector Model")
