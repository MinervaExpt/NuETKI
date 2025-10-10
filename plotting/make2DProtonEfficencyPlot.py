#!/usr/bin/python

#USAGE: python make2DProtonEfficiencyPlot.py <SelectionWithProtonCut.root> <SelectionWithoutProtonCut.root>

#Makes a 2D, MAD proton reco efficiency plot by using all selected events with proton cut, divided by all selected events without proton cut.
#For this plot, I don't care about signal in the sense of nue ccqe, I care about signal in the sense of has a true proton.
#Which is already required by the fact that these plots are in true proton vars.
#Any event which does not have a true proton is sitting at -999 for both values

#to make it not display the canvases as it draws and saves them, saves a bunch of time
import sys
sys.argv.append( '-b' )

import ROOT
from ROOT import PlotUtils

import sys
from array import array
import math

#names in the root file of the histograms we're interested in, should be all the bkgd categories plus efficiency numerator plot
signalName = "Theta_p_ProtonP_efficiency_numerator"
bkgdNames = ["Theta_p_ProtonP_by_BKG_Label_Bkg_Wrong_Sign", "Theta_p_ProtonP_by_BKG_Label_NC_Bkg","Theta_p_ProtonP_by_BKG_Label_Other"]

numeratorFile = ROOT.TFile.Open(sys.argv[1])
denomFile = ROOT.TFile.Open(sys.argv[2])
plotter = PlotUtils.MnvPlotter()

#util = PlotUtils.HistogramUtils();

#Gets the signal MnvH2D with and without proton cut
numer = numeratorFile.Get(signalName)
denomin = denomFile.Get(signalName)

#Gets their underlying TH2D's
numerator = numer.GetCVHistoWithStatError()
print("numerator -1: ",numerator.Integral(0,18,0,20))
denominator = denomin.GetCVHistoWithStatError()
print("denominator -1: ", denominator.Integral(0,18,0,20))

#Now loop through background histos and add them to the signal category (remember we want all selected events, signal + all bkgds)
for i in range(len(bkgdNames)): 

    #Gets the MnvH2D's
    num = numeratorFile.Get(bkgdNames[i])
    denom = denomFile.Get(bkgdNames[i])

    #Gets the underlying TH2D's
    numm = num.GetCVHistoWithStatError()
    denomm = denom.GetCVHistoWithStatError()
    print("numerator ",i, ": ",  numerator.Integral(0,18,0,20))
    print("denominator ",i,": ", denominator.Integral(0,18,0,20))


    numerator.Add(numm)
    denominator.Add(denomm)
    print("numerator -1 after addition ",i, ": ",  numerator.Integral(0,18,0,20))
    print("denominator -1 after addition ",i,": ", denominator.Integral(0,18,0,20))


numerator.SetTitle('2D Proton Reco Efficiency')

#numerator.Divide(denominator)
    
#Removing events below 450 MeV proton momentum, these cannot possibly be true reco'd protons.
#These are events with reco'd final protons and true final state protons, but the momentum is too low for those 
#final state protons to be reco'd. So, the "reco'd proton" is probably something else, mis-reconstructed
#These will be removed once we add in an ESC type cut or something most likely

#for i in range(18): #x bins 
#    for j in range(20): #y bins
#        if (j < 6):
#            numerator.SetBinContent(i, j, 0)

#this is super hacky lol
#array = ROOT.TObjArray()
#array.Add(numerator)

outName=signalName.replace("_efficiency_numerator", "_efficiency")
canvas = ROOT.TCanvas( 'canvas', outName, 0, 0, 2000, 1600 )

numerator.GetXaxis().SetTitle("True Proton Angle [deg]")
numerator.GetYaxis().SetTitle("True Proton Momentum [GeV]")


numerator.SetFillColorAlpha(ROOT.kBlue, 0.5)
numerator.SetLineWidth(3)
numerator.Draw("colz")

#denominator.Draw("hist")
    
outName = outName + ".png"
canvas.SaveAs(outName)
canvas.Delete()
#plotter.WritePreliminaraxy()


