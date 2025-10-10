#!/usr/bin/python

#to make it not display the canvases as it draws and saves them, saves a bunch of time
import sys

import ROOT
from ROOT import PlotUtils
sys.argv.append( '-b' )
ROOT.gROOT.SetBatch(True)

import sys
from array import array
import math
import ctypes


dataFile = ROOT.TFile.Open(sys.argv[1])
mcFile = ROOT.TFile.Open(sys.argv[2])
#dataFile = ROOT.TFile.Open("/exp/minerva/data/users/cpernas/NuE_TKI/Data_July_12_2025_1.5_no_syst.root")
#mcFile = ROOT.TFile.Open("/exp/minerva/data/users/cpernas/NuE_TKI/MC_July_12_2025_1.5_no_syst.root")

plotter = PlotUtils.MnvPlotter()

#plotter.ApplyStyle(PlotUtils.kCCNuPionIncStyle)
plotter.legend_text_size = 0.015
plotter.data_line_width = 2
plotter.data_marker_size = 2
plotter.draw_normalized_to_bin_width = 0

#python to cpp nonsense idk
mcColors = [4, 7, 6, 2, 5, 416]
arr =(ctypes.c_int * len(mcColors))(*mcColors)

rebinFactor = 3 #if I wanna rebin more coarsely, this sets new nbins to old_nbins/rebinFactor

mcPOT = mcFile.Get("POTUsed").GetVal()
dataPOT = dataFile.Get("POTUsed").GetVal()
mcScale = dataPOT/mcPOT

print("mc POT scale = ", mcScale)

bkgdCategoryNames = ["selected_signal_reco", "background_NuECC_with_pions", "background_Other_NueCC", "background_NC_pi0", "background_CC_Numu_pi0", "background_Other"]

sidebands = ["_", "_MeanFrontDEDXSB_", "_MichelSB_"]
data_hists = []
mc_hists = [ [], [], [], [], [], [] ] #each embedded list is a bkgd category in same order as bkgd category names, and will have 3 (the # of sidebands) entries

#eventually I'm gonna need another loop of the actual variables... eep
for index, region in enumerate(sidebands): #gonna just count signal as one of these sidebands, keep that in mind
    data_hists.append(dataFile.Get("DeltaPt" + region + "data"))
    data_hists[index].SetDirectory(0)
    rebinned_data = data_hists[index].Rebin(rebinFactor, data_hists[index].GetName())
    data_hists[index] = rebinned_data
    #get all of the mc histos separated by bkgd category
    for j, category in enumerate(bkgdCategoryNames):
        mc_hists[j].append(mcFile.Get("DeltaPt" + region + category))
        mc_hists[j][index].SetDirectory(0)
        rebinned_MC = mc_hists[j][index].Rebin(rebinFactor, mc_hists[j][index].GetName())
        mc_hists[j][index] = rebinned_MC

nbins = data_hists[0].GetNbinsX() #they should all have the same # of bins...
print("# OF BINS = ", nbins)
meanFrontSB_bkg_scale_factors=[]
meanFrontSB_sig_scale_factors=[]
michelSB_bkg_scale_factors=[]
michelSB_sig_scale_factors=[]

for i in range(1, nbins+1):   # ROOT bins start at 1
    data_signal_region = data_hists[0].GetBinContent(i)
    data_meanFront_sb = data_hists[1].GetBinContent(i)
    data_michel_sb = data_hists[2].GetBinContent(i)

    mc_sig_michel_sb = mc_hists[0][2].GetBinContent(i) * mcScale

    mc_michel_bkg_signal_region = mc_hists[1][0].GetBinContent(i) * mcScale
    #mc_michel_bkg_signal_region = ( mc_hists[1][0].GetBinContent(i) + mc_hists[2][0].GetBinContent(i) ) * mcScale #to include or not include the red background...  
    mc_michel_bkg_michel_sb = mc_hists[1][2].GetBinContent(i) * mcScale
    #mc_michel_bkg_michel_sb = ( mc_hists[1][2].GetBinContent(i) + mc_hists[2][2].GetBinContent(i) ) * mcScale #to include or not include the red background...  

    #For each sideband, subtract off from the data value the bkg contributions I plan to hold fixed, then solve for the rest using the same eq

    #calc scale factors for the mean front dE/dX sideband
    d_s = data_signal_region - (mc_hists[1][0].GetBinContent(i)+mc_hists[2][0].GetBinContent(i)+mc_hists[5][0].GetBinContent(i))*mcScale
    d_sb = data_meanFront_sb - (mc_hists[1][1].GetBinContent(i)+mc_hists[2][1].GetBinContent(i)+mc_hists[5][1].GetBinContent(i))*mcScale
    
    mc_sig_s = mc_hists[0][0].GetBinContent(i) * mcScale
    mc_sig_sb = mc_hists[0][1].GetBinContent(i) * mcScale
    mc_bkg_s = ( mc_hists[3][0].GetBinContent(i) + mc_hists[4][0].GetBinContent(i) ) * mcScale
    mc_bkg_sb = ( mc_hists[3][1].GetBinContent(i) + mc_hists[4][1].GetBinContent(i) ) * mcScale

    delta = mc_sig_s*mc_bkg_sb - mc_sig_sb*mc_bkg_s
    bkg_scale_1 = ( mc_sig_sb*(-d_s) + mc_sig_s*d_sb ) / delta
    sig_scale_1 = ( d_s*mc_bkg_sb - d_sb*mc_bkg_s ) / delta

    #calc scale factors for the Michel sideband
    d_s = data_signal_region - (mc_hists[2][0].GetBinContent(i)+mc_hists[3][0].GetBinContent(i)+mc_hists[4][0].GetBinContent(i)+mc_hists[5][0].GetBinContent(i))*mcScale
    d_sb  = data_michel_sb - (mc_hists[2][2].GetBinContent(i)+mc_hists[3][2].GetBinContent(i)+mc_hists[4][2].GetBinContent(i)+mc_hists[5][2].GetBinContent(i))*mcScale

    mc_sig_s = mc_hists[0][0].GetBinContent(i) * mcScale
    mc_sig_sb = mc_hists[0][2].GetBinContent(i) * mcScale
    mc_bkg_s = ( mc_hists[1][0].GetBinContent(i) ) * mcScale
    mc_bkg_sb = ( mc_hists[1][2].GetBinContent(i) ) * mcScale

    delta = mc_sig_s*mc_bkg_sb - mc_sig_sb*mc_bkg_s
    bkg_scale_2 = ( mc_sig_sb*(-d_s) + mc_sig_s*d_sb ) / delta
    sig_scale_2 = ( d_s*mc_bkg_sb - d_sb*mc_bkg_s ) / delta

    meanFrontSB_bkg_scale_factors.append(bkg_scale_1)
    meanFrontSB_sig_scale_factors.append(sig_scale_1)

    michelSB_bkg_scale_factors.append(bkg_scale_2)
    michelSB_sig_scale_factors.append(sig_scale_2)

    print(f"MeanFront, Bin {i:2d}: bkg scale={bkg_scale_1:.2f}, sig_scale={sig_scale_1:.2f}")
    print(f"Michel, Bin {i:2d}: bkg scale={bkg_scale_2:.2f}, sig_scale={sig_scale_2:.2f}\n")

    print(f"Data in Signal Region = {data_signal_region:.2f}")
    print(f"Data in MeanFront SB  = {data_meanFront_sb:.2f}\n")

    unscaled_total_mc_meanfront = (mc_hists[0][1].GetBinContent(i) + mc_hists[1][1].GetBinContent(i) + mc_hists[2][1].GetBinContent(i) + mc_hists[3][1].GetBinContent(i) + mc_hists[4][1].GetBinContent(i) + mc_hists[5][1].GetBinContent(i))*mcScale
    scaled_total_mc_meanfront = (sig_scale_1*mc_hists[0][1].GetBinContent(i) + mc_hists[1][1].GetBinContent(i) + mc_hists[2][1].GetBinContent(i) + bkg_scale_1*mc_hists[3][1].GetBinContent(i) + bkg_scale_1*mc_hists[4][1].GetBinContent(i) + mc_hists[5][1].GetBinContent(i))*mcScale


    print(f"unscaled mc : {unscaled_total_mc_meanfront:.2f}")
    print(f"unscaled signal   : {(mc_hists[0][1].GetBinContent(i)*mcScale):.2f}")
    print(f"unscaled nonQE    : {(mc_hists[1][1].GetBinContent(i)*mcScale):.2f}")
    print(f"unscaled otherNue : {(mc_hists[2][1].GetBinContent(i)*mcScale):.2f}")
    print(f"unscaled NCPi0    : {(mc_hists[3][1].GetBinContent(i)*mcScale):.2f}")
    print(f"unscaled CCnumuPi0: {(mc_hists[4][1].GetBinContent(i)*mcScale):.2f}")
    print(f"unscaled otherbkgd: {(mc_hists[5][1].GetBinContent(i)*mcScale):.2f}\n")

    print(f"scaled mc : {scaled_total_mc_meanfront:.2f}")
    print(f"scaled signal   : {(mc_hists[0][1].GetBinContent(i)*sig_scale_1*mcScale):.2f}")
    print(f"scaled nonQE    : {(mc_hists[1][1].GetBinContent(i)*mcScale):.2f}")
    print(f"scaled otherNue : {(mc_hists[2][1].GetBinContent(i)*mcScale):.2f}")
    print(f"scaled NCPi0    : {(mc_hists[3][1].GetBinContent(i)*bkg_scale_1*mcScale):.2f}")
    print(f"scaled CCnumuPi0: {(mc_hists[4][1].GetBinContent(i)*bkg_scale_1*mcScale):.2f}")
    print(f"scaled otherbkgd: {(mc_hists[5][1].GetBinContent(i)*mcScale):.2f}\n")

    
    print("-------------")
print("\n------------------------------------------------\n")

#------------
# Plotting
#------------

#Histograms of scale factors
meanFrontSB_bkg_scale_factors_hist = data_hists[0].Clone("meanFrontSB_bkg_scale_factors_hist")
meanFrontSB_sig_scale_factors_hist = data_hists[0].Clone("meanFrontSB_sig_scale_factors_hist")
meanFrontSB_bkg_scale_factors_hist.SetTitle("Background Scale Factors")
meanFrontSB_sig_scale_factors_hist.SetTitle("Signal Scale Factors")
meanFrontSB_bkg_scale_factors_hist.Reset("ICES")
meanFrontSB_sig_scale_factors_hist.Reset("ICES")

michelSB_bkg_scale_factors_hist = data_hists[0].Clone("michelSB_bkg_scale_factors_hist")
michelSB_sig_scale_factors_hist = data_hists[0].Clone("michelSB_sig_scale_factors_hist")
michelSB_bkg_scale_factors_hist.SetTitle("Background Scale Factors")
michelSB_sig_scale_factors_hist.SetTitle("Signal Scale Factors")
michelSB_bkg_scale_factors_hist.Reset("ICES")
michelSB_sig_scale_factors_hist.Reset("ICES")
for i in range(1, nbins+1):
    meanFrontSB_bkg_scale_factors_hist.SetBinContent(i, meanFrontSB_bkg_scale_factors[i-1])
    meanFrontSB_sig_scale_factors_hist.SetBinContent(i, meanFrontSB_sig_scale_factors[i-1])

    michelSB_bkg_scale_factors_hist.SetBinContent(i, michelSB_bkg_scale_factors[i-1])
    michelSB_sig_scale_factors_hist.SetBinContent(i, michelSB_sig_scale_factors[i-1])

canvas = ROOT.TCanvas( 'canvas', "test", 0, 0, 2000, 1600 )
meanFrontSB_bkg_scale_factors_hist.Draw("E")
canvas.SaveAs("meanFrontSB_bkg_scale_factors.png")
canvas.Delete()
canvas = ROOT.TCanvas( 'canvas', "test", 0, 0, 2000, 1600 )
meanFrontSB_sig_scale_factors_hist.Draw("E")
canvas.SaveAs("meanFrontSB_sig_scale_factors.png")
canvas.Delete()
canvas = ROOT.TCanvas( 'canvas', "test", 0, 0, 2000, 1600 )
michelSB_bkg_scale_factors_hist.Draw("E")
canvas.SaveAs("michelSB_bkg_scale_factors.png")
canvas.Delete()
canvas = ROOT.TCanvas( 'canvas', "test", 0, 0, 2000, 1600 )
michelSB_sig_scale_factors_hist.Draw("E")
canvas.SaveAs("michelSB_sig_scale_factors.png")
canvas.Delete()


#Applying the scale factors and plotting respective hists...
#plotting unscaled, both scaled, and bkg scaled only for each sideband
for i in range(1, len(sidebands)): 
    #super hacky way of figuring out which scale factors to use...
    if i == 1: #meanFrontSB
        sig_scale_factors = meanFrontSB_sig_scale_factors
        bkg_scale_factors = meanFrontSB_bkg_scale_factors
    elif i == 2: #michelSB
        sig_scale_factors = michelSB_sig_scale_factors
        bkg_scale_factors = michelSB_bkg_scale_factors            
    else: #signal region
        sig_scale_factors = [1, 1, 1, 1, 1]
        bkg_scale_factors = [1, 1, 1, 1, 1]
            
    #plot unscaled
    mc_hists[0][i].SetTitle('signal (nu_e QELike + proton)')
    mc_hists[1][i].SetTitle('nu_e nonQE (has FS mesons)')
    mc_hists[2][i].SetTitle('Other nu_eCC')
    mc_hists[3][i].SetTitle('NC with pi0')
    mc_hists[4][i].SetTitle('nu_mu CC with pi0')
    mc_hists[5][i].SetTitle('other')
    data_hists[i].SetTitle('data')
    
    array = ROOT.TObjArray()
    for k in range(5, -1, -1):
        array.Add(mc_hists[k][i])
        
    canvas = ROOT.TCanvas( 'canvas', "test", 0, 0, 2000, 1600 )
    
    plotter.DrawDataStackedMC(data_hists[i], array, arr, mcScale, "TR", "Data", 1001, mc_hists[0][0].GetXaxis().GetTitle(), "N events")
    normalization = f"Normalized to {dataPOT:.3} data POT"
    plotter.AddPlotLabel(normalization, 0.2, 0.95, 0.03)
    
    outName = "DeltaPt" + sidebands[i] + ".png"
    canvas.SaveAs(outName)
    canvas.Delete()

    NueCCPion = mc_hists[1][i].Clone()
    NueCCPion.Reset("ICES")
    otherNueCC = mc_hists[2][i].Clone()
    otherNueCC.Reset("ICES")
    NCPi0 = mc_hists[3][i].Clone()
    NCPi0.Reset("ICES")
    CCnumuPi0 = mc_hists[4][i].Clone()
    CCnumuPi0.Reset("ICES")
    otherBkgd = mc_hists[5][i].Clone()
    otherBkgd.Reset("ICES")
    
    #apply bkg scale factors only, first
    for ibin in range(1, nbins+1):
        #then apply bkg scale factors and plot those 
        #NueCCPion, bkg1
        bkg1_content = mc_hists[1][i].GetBinContent(ibin)
        bkg1_error   = mc_hists[1][i].GetBinError(ibin)
        if i == 2: #only apply for michel sb, or if we're doing the whole thing
            NueCCPion.SetBinContent(ibin, bkg1_content * bkg_scale_factors[ibin-1])
            NueCCPion.SetBinError(ibin, bkg1_error * bkg_scale_factors[ibin-1])
        else:
            NueCCPion.SetBinContent(ibin, bkg1_content)
            NueCCPion.SetBinError(ibin, bkg1_error)
        #print(f"UNSCALED NueCCPion count in bin {ibin}: {bkg1_content*mcScale:.2f}")
        #print(f"SCALED NueCCPion count in bin {ibin}: {NueCCPion.GetBinContent(ibin)*mcScale:.2f}")

        #otherNueCC, bkg2 
        bkg2_content = mc_hists[2][i].GetBinContent(ibin)
        bkg2_error   = mc_hists[2][i].GetBinError(ibin)
        #gonna test not scaling this one yet...
        otherNueCC.SetBinContent(ibin, bkg2_content)
        otherNueCC.SetBinError(ibin, bkg2_error)
        #print(f"UNSCALED otherNueCC count in bin {ibin}: {bkg2_content*mcScale:.2f}")
        #print(f"SCALED otherNueCC count in bin {ibin}: {otherNueCC.GetBinContent(ibin)*mcScale:.2f}")

        #NCPi0, bkg3
        bkg3_content = mc_hists[3][i].GetBinContent(ibin)
        bkg3_error   = mc_hists[3][i].GetBinError(ibin)
        if i == 1 or i == 0: #only apply for michel sb, or if we're doing the whole thing
            NCPi0.SetBinContent(ibin, bkg3_content * bkg_scale_factors[ibin-1])
            NCPi0.SetBinError(ibin, bkg3_error * bkg_scale_factors[ibin-1])
        else:
            NCPi0.SetBinContent(ibin, bkg3_content)
            NCPi0.SetBinError(ibin, bkg3_error)
        #print(f"UNSCALED NCPi0 count in bin {ibin}: {bkg3_content*mcScale:.2f}")
        #print(f"SCALED NCPi0 count in bin {ibin}: {NCPi0.GetBinContent(ibin)*mcScale:.2f}")

        #CCnumuPi0, bkg4
        bkg4_content = mc_hists[4][i].GetBinContent(ibin)
        bkg4_error   = mc_hists[4][i].GetBinError(ibin)
        if i == 1 or i == 0: #only apply for michel sb, or if we're doing the whole thing
            CCnumuPi0.SetBinContent(ibin, bkg4_content * bkg_scale_factors[ibin-1])
            CCnumuPi0.SetBinError(ibin, bkg4_error * bkg_scale_factors[ibin-1])
        else:
            CCnumuPi0.SetBinContent(ibin, bkg4_content)
            CCnumuPi0.SetBinError(ibin, bkg4_error)            
        #print(f"UNSCALED CCnumuPi0 count in bin {ibin}: {bkg4_content*mcScale:.2f}")
        #print(f"SCALED CCnumuPi0 count in bin {ibin}: {CCnumuPi0.GetBinContent(ibin)*mcScale:.2f}")

        #otherBkgd, bkg5
        bkg5_content = mc_hists[5][i].GetBinContent(ibin)
        bkg5_error   = mc_hists[5][i].GetBinError(ibin)
        otherBkgd.SetBinContent(ibin, bkg5_content)
        otherBkgd.SetBinError(ibin, bkg5_error)
        #print(f"UNSCALED otherBkgd count in bin {ibin}: {bkg5_content*mcScale:.2f}")
        #print(f"SCALED otherBkgd count in bin {ibin}: {otherBkgd.GetBinContent(ibin)*mcScale:.2f}")

        #print(f"UNSCALED MC total for bin {ibin} = {(sig_content + bkg1_content + bkg2_content + bkg3_content + bkg4_content + bkg5_content)*mcScale:.2f}")
        #print(f"SCALED MC total for bin {ibin} = {(signal.GetBinContent(ibin) + NueCCPion.GetBinContent(ibin) + otherNueCC.GetBinContent(ibin) + NCPi0.GetBinContent(ibin) + CCnumuPi0.GetBinContent(ibin) + otherBkgd.GetBinContent(ibin))*mcScale:.2f}\n")

    array = ROOT.TObjArray()
    array.Add(otherBkgd)
    array.Add(CCnumuPi0)
    array.Add(NCPi0)
    array.Add(otherNueCC)
    array.Add(NueCCPion)
    array.Add(mc_hists[0][i])
    
    canvas = ROOT.TCanvas( 'canvas', "test", 0, 0, 2000, 1600 )
    
    plotter.DrawDataStackedMC(data_hists[i], array, arr, mcScale, "TR", "Data", 1001, mc_hists[0][0].GetXaxis().GetTitle(), "N events")
    normalization = f"Normalized to {dataPOT:.3} data POT"
    plotter.AddPlotLabel(normalization, 0.2, 0.95, 0.03)
    
    outName = "DeltaPt" + sidebands[i] + "bkgd_scaled_ONLY.png"
    canvas.SaveAs(outName)
    canvas.Delete()

    signal = mc_hists[0][i].Clone()
    signal.Reset("ICES")

    #then apply sig scale factors and plot those
    for ibin in range(1, nbins+1):        
        sig_content = mc_hists[0][i].GetBinContent(ibin)
        sig_error   = mc_hists[0][i].GetBinError(ibin)
        signal.SetBinContent(ibin, sig_content * sig_scale_factors[ibin-1])
        signal.SetBinError(ibin, sig_error * sig_scale_factors[ibin-1])
        #print(f"UNSCALED signal count in bin {ibin}: {sig_content*mcScale:.2f}")
        #print(f"SCALED signal count in bin {ibin}: {signal.GetBinContent(ibin)*mcScale:.2f}")

    array = ROOT.TObjArray()
    array.Add(otherBkgd)
    array.Add(CCnumuPi0)
    array.Add(NCPi0)
    array.Add(otherNueCC)
    array.Add(NueCCPion)
    array.Add(signal)
    
    canvas = ROOT.TCanvas( 'canvas', "test", 0, 0, 2000, 1600 )
    
    plotter.DrawDataStackedMC(data_hists[i], array, arr, mcScale, "TR", "Data", 1001, signal.GetXaxis().GetTitle(), "N events")
    normalization = f"Normalized to {dataPOT:.3} data POT"
    plotter.AddPlotLabel(normalization, 0.2, 0.95, 0.03)
    
    outName = "DeltaPt" + sidebands[i] + "bkgd_AND_sig_scaled.png"
    canvas.SaveAs(outName)
    canvas.Delete()
        
#Now plot signal region, completely unscaled, and then FULLY scaled using both sets of scale factors

# ------ PLOTTING SIGNAL REGION, UNSCALED ----------
mc_hists[0][0].SetTitle('signal (nu_e QELike + proton)')
mc_hists[1][0].SetTitle('nu_e nonQE (has FS mesons)')
mc_hists[2][0].SetTitle('Other nu_eCC')
mc_hists[3][0].SetTitle('NC with pi0')
mc_hists[4][0].SetTitle('nu_mu CC with pi0')
mc_hists[5][0].SetTitle('other')
data_hists[0].SetTitle('data')

array = ROOT.TObjArray()
for k in range(5, -1, -1):
    array.Add(mc_hists[k][0])

canvas = ROOT.TCanvas( 'canvas', "test", 0, 0, 2000, 1600 )

plotter.DrawDataStackedMC(data_hists[0], array, arr, mcScale, "TR", "Data", 1001, mc_hists[0][0].GetXaxis().GetTitle(), "N events")
normalization = f"Normalized to {dataPOT:.3} data POT"
plotter.AddPlotLabel(normalization, 0.2, 0.95, 0.03)
            
outName = "DeltaPt" + ".png"
canvas.SaveAs(outName)
canvas.Delete()

# --------------- PLOTTING SIGNAL REGION, BACKGROUND AND SIGNAL SCALED TO MEANFRONT -------------


signal = mc_hists[0][0].Clone()
signal.Reset("ICES")
NueCCPion = mc_hists[1][0].Clone()
#NueCCPion.Reset("ICES")
otherNueCC = mc_hists[2][0].Clone()
#otherNueCC.Reset("ICES")
NCPi0 = mc_hists[3][0].Clone()
NCPi0.Reset("ICES")
CCnumuPi0 = mc_hists[4][0].Clone()
CCnumuPi0.Reset("ICES")
otherBkgd = mc_hists[5][0].Clone()
#otherBkgd.Reset("ICES")

for ibin in range(1, nbins+1):
    #then apply bkg scale factors and plot those 
    #signal
    sig_content = mc_hists[0][0].GetBinContent(ibin)
    sig_error   = mc_hists[0][0].GetBinError(ibin)
    signal.SetBinContent(ibin, sig_content * meanFrontSB_sig_scale_factors[ibin-1])
    signal.SetBinError(ibin, sig_error * meanFrontSB_sig_scale_factors[ibin-1])

    #NCPi0, bkg3
    bkg3_content = mc_hists[3][0].GetBinContent(ibin)
    bkg3_error   = mc_hists[3][0].GetBinError(ibin)
    NCPi0.SetBinContent(ibin, bkg3_content * meanFrontSB_bkg_scale_factors[ibin-1])
    NCPi0.SetBinError(ibin, bkg3_error * meanFrontSB_bkg_scale_factors[ibin-1])
    
    #CCnumuPi0, bkg4
    bkg4_content = mc_hists[4][0].GetBinContent(ibin)
    bkg4_error   = mc_hists[4][0].GetBinError(ibin)
    CCnumuPi0.SetBinContent(ibin, bkg4_content * meanFrontSB_bkg_scale_factors[ibin-1])
    CCnumuPi0.SetBinError(ibin, bkg4_error * meanFrontSB_bkg_scale_factors[ibin-1])    

array = ROOT.TObjArray()
array.Add(otherBkgd)
array.Add(CCnumuPi0)
array.Add(NCPi0)
array.Add(otherNueCC)
array.Add(NueCCPion)
array.Add(signal)

canvas = ROOT.TCanvas( 'canvas', "test", 0, 0, 2000, 1600 )

plotter.DrawDataStackedMC(data_hists[0], array, arr, mcScale, "TR", "Data", 1001, signal.GetXaxis().GetTitle(), "N events")
normalization = f"Normalized to {dataPOT:.3} data POT"
plotter.AddPlotLabel(normalization, 0.2, 0.95, 0.03)
            
outName = "DeltaPt_SIG_MEANFRONT_SCALED.png"
canvas.SaveAs(outName)
canvas.Delete()

#------------- PLOTTING SIGNAL REGION, BACKGROUND AND SIGNAL SCALED TO MICHEL -----------
signal = mc_hists[0][0].Clone()
signal.Reset("ICES")
NueCCPion = mc_hists[1][0].Clone()
NueCCPion.Reset("ICES")
otherNueCC = mc_hists[2][0].Clone()
otherNueCC.Reset("ICES")
NCPi0 = mc_hists[3][0].Clone()
#NCPi0.Reset("ICES")
CCnumuPi0 = mc_hists[4][0].Clone()
#CCnumuPi0.Reset("ICES")
otherBkgd = mc_hists[5][0].Clone()
#otherBkgd.Reset("ICES")

for ibin in range(1, nbins+1):
    #then apply bkg scale factors and plot those 
    #signal
    sig_content = mc_hists[0][0].GetBinContent(ibin)
    sig_error   = mc_hists[0][0].GetBinError(ibin)
    signal.SetBinContent(ibin, sig_content * michelSB_sig_scale_factors[ibin-1])
    signal.SetBinError(ibin, sig_error * michelSB_sig_scale_factors[ibin-1])
    
    #NueCCPion, bkg1
    bkg1_content = mc_hists[1][0].GetBinContent(ibin)
    bkg1_error   = mc_hists[1][0].GetBinError(ibin)
    NueCCPion.SetBinContent(ibin, bkg1_content * michelSB_bkg_scale_factors[ibin-1])
    NueCCPion.SetBinError(ibin, bkg1_error * michelSB_bkg_scale_factors[ibin-1])

    #otherNueCC, bkg2 
    bkg2_content = mc_hists[2][0].GetBinContent(ibin)
    bkg2_error   = mc_hists[2][0].GetBinError(ibin)
    #gonna test not scaling this one yet...
    otherNueCC.SetBinContent(ibin, bkg2_content)
    otherNueCC.SetBinError(ibin, bkg2_error)
    
array = ROOT.TObjArray()
array.Add(otherBkgd)
array.Add(CCnumuPi0)
array.Add(NCPi0)
array.Add(otherNueCC)
array.Add(NueCCPion)
array.Add(signal)

canvas = ROOT.TCanvas( 'canvas', "test", 0, 0, 2000, 1600 )

plotter.DrawDataStackedMC(data_hists[0], array, arr, mcScale, "TR", "Data", 1001, signal.GetXaxis().GetTitle(), "N events")
normalization = f"Normalized to {dataPOT:.3} data POT"
plotter.AddPlotLabel(normalization, 0.2, 0.95, 0.03)
            
outName = "DeltaPt_SIG_MICHEL_SCALED.png"
canvas.SaveAs(outName)
canvas.Delete()

#------------ PLOTTING FINAL, FULLY SCALED SIGNAL REGION WITH SIGNAL UNSCALED AND BACKGROUNDS WITH THEIR RESPECTIVE SCALE FACTORS

signal = mc_hists[0][0].Clone()
#signal.Reset("ICES")
NueCCPion = mc_hists[1][0].Clone()
NueCCPion.Reset("ICES")
otherNueCC = mc_hists[2][0].Clone()
otherNueCC.Reset("ICES")
NCPi0 = mc_hists[3][0].Clone()
NCPi0.Reset("ICES")
CCnumuPi0 = mc_hists[4][0].Clone()
CCnumuPi0.Reset("ICES")
otherBkgd = mc_hists[5][0].Clone()
#otherBkgd.Reset("ICES")

for ibin in range(1, nbins+1):
    #then apply bkg scale factors and plot those 
    #NueCCPion, bkg1
    bkg1_content = mc_hists[1][0].GetBinContent(ibin)
    bkg1_error   = mc_hists[1][0].GetBinError(ibin)
    NueCCPion.SetBinContent(ibin, bkg1_content * michelSB_bkg_scale_factors[ibin-1])
    NueCCPion.SetBinError(ibin, bkg1_error * michelSB_bkg_scale_factors[ibin-1])

    #otherNueCC, bkg2 
    bkg2_content = mc_hists[2][0].GetBinContent(ibin)
    bkg2_error   = mc_hists[2][0].GetBinError(ibin)
    #gonna test not scaling this one yet...
    otherNueCC.SetBinContent(ibin, bkg2_content)
    otherNueCC.SetBinError(ibin, bkg2_error)
    
    #NCPi0, bkg3
    bkg3_content = mc_hists[3][0].GetBinContent(ibin)
    bkg3_error   = mc_hists[3][0].GetBinError(ibin)
    NCPi0.SetBinContent(ibin, bkg3_content * meanFrontSB_bkg_scale_factors[ibin-1])
    NCPi0.SetBinError(ibin, bkg3_error * meanFrontSB_bkg_scale_factors[ibin-1])
    
    #CCnumuPi0, bkg4
    bkg4_content = mc_hists[4][0].GetBinContent(ibin)
    bkg4_error   = mc_hists[4][0].GetBinError(ibin)
    CCnumuPi0.SetBinContent(ibin, bkg4_content * meanFrontSB_bkg_scale_factors[ibin-1])
    CCnumuPi0.SetBinError(ibin, bkg4_error * meanFrontSB_bkg_scale_factors[ibin-1])    

array = ROOT.TObjArray()
array.Add(otherBkgd)
array.Add(CCnumuPi0)
array.Add(NCPi0)
array.Add(otherNueCC)
array.Add(NueCCPion)
array.Add(signal)

canvas = ROOT.TCanvas( 'canvas', "test", 0, 0, 2000, 1600 )

plotter.DrawDataStackedMC(data_hists[0], array, arr, mcScale, "TR", "Data", 1001, signal.GetXaxis().GetTitle(), "N events")
normalization = f"Normalized to {dataPOT:.3} data POT"
plotter.AddPlotLabel(normalization, 0.2, 0.95, 0.03)
            
outName = "DeltaPt" + "_FULLY_SCALED" + ".png"
canvas.SaveAs(outName)
canvas.Delete()
#print(f"bkg scale factors: {bkg_scale_factors=}")    
#print(f"sig scale factors: {sig_scale_factors=}")    
    

