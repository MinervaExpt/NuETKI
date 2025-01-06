#!/usr/bin/python

#USAGE: rebin.py <Histos.root>

#Util to rebin histograms that were output by runEventLoop, useful if I ran a huge amount of files and then realize my binning sucks. 
#Probably a file that's going to be edited and messed with often. 

#to make it not display the canvases as it draws and saves them, saves a bunch of time
import ROOT
from ROOT import PlotUtils

import sys
from array import array
import math

#I'm going to try and make it rebin ALL histograms of a certain variable, don't see any reason why protonKE_signal and protonKE_bkgd1 should have different binning
varToRebin = "ProtonKE"
newBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.5,1.75,2.0,2.5,3.0]

inFile = ROOT.TFile.Open(sys.argv[1])


#plotter = PlotUtils.MnvPlotter()

#oldHist = inFile.Get("ProtonKE_data")
oldHist = inFile.Get("ProtonKE_efficiency_numerator") 
oldHist2 = inFile.Get("ProtonKE_efficiency_denominator")
#inFile.Close()

#oldHist.SetTitle("really")

print(type(oldHist))

#newHist = oldHist.Clone()
newHist = oldHist.GetCVHistoWithStatError()
newHist3 = oldHist2.GetCVHistoWithStatError()


newBins2 = array('d', newBins)
newHist2 = newHist.Rebin(18, "newHist2", newBins2)
newHist4 = newHist3.Rebin(18, "newHist4", newBins2)

outFile = ROOT.TFile.Open("outTest.root", "RECREATE")

newHist2.Write()
newHist4.Write()
outFile.Close()



