#!/usr/bin/python

#USAGE: backgroundStack.py <mcFile1.root> <mcFile2.root>

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

#names in the root file of the histograms we're interested in, just put the signal and it'll find the other bkgd categories from that
signalHistoNames = ["E_nu_selected_signal_reco", "E_avail_selected_signal_reco", "Lepton_Pt_selected_signal_reco", "E_lep_selected_signal_reco", "Theta_lep_selected_signal_reco"]

numuCCHistoNames=[]
NCCohHistoNames=[]
otherNCHistoNames=[]
nueElasticHistoNames=[]
otherBackgroundHistoNames=[]

dataHistoNames=[]

for i in signalHistoNames:
    numuCCHistoNames.append(i.replace("selected_signal_reco", "background_CC_numu"))
    NCCohHistoNames.append(i.replace("selected_signal_reco", "background_NC_Coh"))
    otherNCHistoNames.append(i.replace("selected_signal_reco", "background_Other_NC"))
    nueElasticHistoNames.append(i.replace("selected_signal_reco", "background_nu___e_Elastic"))
    otherBackgroundHistoNames.append(i.replace("selected_signal_reco", "background_Other"))

    dataHistoNames.append(i.replace("selected_signal_reco", "data"))

mcFile1 = ROOT.TFile.Open(sys.argv[1])
mcFile2 = ROOT.TFile.Open(sys.argv[2])

plotter = PlotUtils.MnvPlotter()

#plotter.ApplyStyle(PlotUtils.kCCNuPionIncStyle)
plotter.draw_normalized_to_bin_width = False
plotter.legend_text_size = 0.02
plotter.axis_draw_grid_y = True
if drawData:
    plotter.data_line_width = 2
    plotter.data_marker_size = 2
else:
    plotter.data_line_width = 0
    plotter.data_marker_size = 0

#kBlack, kPink+1, kOrange+2, kBlue-7, kYellow, kOrange
mcColors = [1, 901, 802, 593, 400, 800]
arr =(ctypes.c_int * len(mcColors))(*mcColors)



mcPOT1 = mcFile1.Get("POTUsed").GetVal()
mcPOT2 = mcFile2.Get("POTUsed").GetVal()

POT_ratio = mcPOT2/mcPOT1 #multiply the numerator by this number
print("pot ratio: ", POT_ratio)

#I need to pass this guy a TObjArray of two histograms -> one that is pTmu_selected_signal_reco, and one that is all of the backgrounds
#for all the backgrounds I can either add up the 3 background histograms (pTmu_background_Wrongsign+pTmu_background_NC+pTmu_background_other)
#or I can do all minus signal (pTmu_data - pTmu_selected_signal_reco)

for i in range(len(signalHistoNames)):
    _signal1 = mcFile1.Get(signalHistoNames[i])
    _numuCC1 = mcFile1.Get(numuCCHistoNames[i])
    _NCCoh1 = mcFile1.Get(NCCohHistoNames[i])
    _otherNC1 = mcFile1.Get(otherNCHistoNames[i])
    _nueElastic1 = mcFile1.Get(nueElasticHistoNames[i])
    _otherBkgd1 = mcFile1.Get(otherBackgroundHistoNames[i])

    _signal2 = mcFile2.Get(signalHistoNames[i])
    _numuCC2 = mcFile2.Get(numuCCHistoNames[i])
    _NCCoh2 = mcFile2.Get(NCCohHistoNames[i])
    _otherNC2 = mcFile2.Get(otherNCHistoNames[i])
    _nueElastic2 = mcFile2.Get(nueElasticHistoNames[i])
    _otherBkgd2 = mcFile2.Get(otherBackgroundHistoNames[i])

    signal1 = _signal1.GetCVHistoWithStatError()
    numuCC1 = _numuCC1.GetCVHistoWithStatError()
    NCCoh1 = _NCCoh1.GetCVHistoWithStatError()
    otherNC1 = _otherNC1.GetCVHistoWithStatError()
    nueElastic1 = _nueElastic1.GetCVHistoWithStatError()
    otherBkgd1 = _otherBkgd1.GetCVHistoWithStatError()

    signal2 = _signal2.GetCVHistoWithStatError()
    numuCC2 = _numuCC2.GetCVHistoWithStatError()
    NCCoh2 = _NCCoh2.GetCVHistoWithStatError()
    otherNC2 = _otherNC2.GetCVHistoWithStatError()
    nueElastic2 = _nueElastic2.GetCVHistoWithStatError()
    otherBkgd2 = _otherBkgd2.GetCVHistoWithStatError()

    #signal1.SetTitle('cc nue')
    #numuCC1.SetTitle('cc numu')
    #NCCoh1.SetTitle('NC Coh')
    #otherNC1.SetTitle('Other NC')
    #nueElastic1.SetTitle('nu + e Elastic')
    #otherBkgd1.SetTitle('other')

    signal1.Add(numuCC1)
    signal1.Add(NCCoh1)
    signal1.Add(otherNC1)
    signal1.Add(nueElastic1)
    signal1.Add(otherBkgd1)

    signal2.Add(numuCC2)
    signal2.Add(NCCoh2)
    signal2.Add(otherNC2)
    signal2.Add(nueElastic2)
    signal2.Add(otherBkgd2)


    #array = ROOT.TObjArray()
    #array.Add(otherBkgd)
    #array.Add(nueElastic)
    #array.Add(otherNC)
    #array.Add(NCCoh)
    #array.Add(numuCC)
    #array.Add(signal)

    #scale numerator to denom POT
    signal1.Scale(POT_ratio)
    
    signal1.Divide(signal2)

    outName=signalHistoNames[i].replace("_selected_signal_reco", "")
    canvas = ROOT.TCanvas( 'canvas', outName, 0, 0, 2000, 1600 )

    #to get variable name for the x axis
    head, sep, tail = signalHistoNames[i].partition('_selected')

    xAxis = signal1.GetXaxis()
    signal1.GetYaxis().SetRangeUser(0.6,1.4)    

    signal1.Draw()

    lowX = xAxis.GetBinLowEdge(xAxis.GetFirst())
    highX = xAxis.GetBinUpEdge(xAxis.GetLast())
    line = ROOT.TLine()
    line.SetLineStyle(2);
    line.SetLineWidth(3);
    line.SetLineColor(36);
    line.DrawLine(lowX, 1., highX, 1.)
    
    
    outName = outName + "_Ratio.png"
    canvas.SaveAs(outName)
    canvas.Delete()
#plotter.WritePreliminary()


