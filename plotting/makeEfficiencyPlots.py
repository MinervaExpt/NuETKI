#!/usr/bin/python

#USAGE: python makeEfficiencyPlots.py <mcHistos.root>

#Makes efficiency plots for the variables listed below and saves them as pngs

#to make it not display the canvases as it draws and saves them, saves a bunch of time
import sys
sys.argv.append( '-b' )

import ROOT
from ROOT import PlotUtils

import sys
from array import array
import math

#names in the root file of the histograms we're interested in, just put the num and it'll find the denom
#numeratorNames = ["E_e_efficiency_numerator", "ProtonKE_efficiency_numerator", "Theta_p_efficiency_numerator", "Pt_p_efficiency_numerator", "Theta_p_e_efficiency_numerator", "Pt_lep_efficiency_numerator"]
numeratorNames = ["ProtonP_efficiency_numerator", "Theta_p_efficiency_numerator","Theta_p_ProtonP_efficiency_numerator"]

denominatorNames = []
for i in numeratorNames:
    denominatorNames.append(i.replace("numerator", "denominator"))

mcFile = ROOT.TFile.Open(sys.argv[1])
plotter = PlotUtils.MnvPlotter()

#util = PlotUtils.HistogramUtils();

for i in range(len(numeratorNames)):
    num = mcFile.Get(numeratorNames[i])
    denom = mcFile.Get(denominatorNames[i])

    numerator = num.GetCVHistoWithStatError()
    denominator = denom.GetCVHistoWithStatError()
    

    numerator.SetTitle('efficiency')

    numerator.Divide(denominator)
    
    #this is super hacky lol
    #array = ROOT.TObjArray()
    #array.Add(numerator)

    outName=numeratorNames[i].replace("_efficiency_numerator", "_efficiency")
    canvas = ROOT.TCanvas( 'canvas', outName, 0, 0, 2000, 1600 )

    #to get variable name for the x axis
    head, sep, tail = numeratorNames[i].partition('_efficiency')

    #arguments for stackedMC array are:
    #DrawStackedMC(mcHists, mcScale, legend position, base color, color offset, fill style, xaxislabel, yaxislabel)
    #plotter.DrawStackedMC(array, 1.0, "TR", 2, 1, 3001, head, "Efficiency")

    #Manually titling my x axes, will change back:
    
    if "E_e_ef" in numeratorNames[i]:        
        numerator.GetXaxis().SetTitle("True Lepton Energy")
    elif "ProtonKE_ef" in numeratorNames[i]:
        numerator.GetXaxis().SetTitle("True Proton KE")
    elif "Theta_p_ef" in numeratorNames[i]:
        numerator.GetXaxis().SetTitle("True Proton Theta")
    elif "Pt_p_ef" in numeratorNames[i]:
        numerator.GetXaxis().SetTitle("True Proton pT")
    elif "Theta_p_e_ef" in numeratorNames[i]:
        numerator.GetXaxis().SetTitle("True opening angle between proton & lepton")
    elif "Pt_lep_ef" in numeratorNames[i]:
        numerator.GetXaxis().SetTitle("True Lepton pT")
    
    numerator.GetYaxis().SetTitle("Efficiency")
    numerator.SetFillColorAlpha(ROOT.kBlue, 0.5)
    numerator.SetLineWidth(3)
    numerator.Draw()
    
    """
    if "ProtonP_ef" in numeratorNames[i]:
        numerator.GetXaxis().SetTitle("Proton Momentum")
        numerator.Draw()
    if "Theta_p_ef" in numeratorNames[i]:
        numerator.GetXaxis().SetTitle("Proton Theta")
        numerator.Draw()
    elif "ProtonP_Theta_p_e_ef" in numeratorNames[i]:
        numerator.GetXaxis().SetTitle("Proton Momentum")
        numerator.GetYaxis().SetTitle("Opening angle between proton & lepton")
        numerator.GetZaxis().SetTitle("Nevents")
        numerator.Draw("colz")
    """

    #denominator.Draw("hist")

    outName = outName + ".png"
    canvas.SaveAs(outName)
    canvas.Delete()
#plotter.WritePreliminaraxy()


