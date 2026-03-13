#!/usr/bin/python

#USAGE: python makeEfficiencyPurityPlots.py

#Script to use N-1 plots to scan through efficiency & purity and also plot some other values:
# significance = signal / sqrt(signal + bkg)
# efficiency * purity
# bkg fraction = 1 - purity = bkg / (signal + bkg)

#right now plotting significance at the same time as the others is broken cause its on a different y scale
# but if you plot it by itself it works fine... maybe ill get around to fixing it, maybe not. 


import ROOT
from ROOT import PlotUtils
import sys
import array
import math
import ctypes

ROOT.gROOT.SetBatch(True)

cutsToDraw = ["Afterpulsing", "DSCalVisE", "E_lep", "EMLikeTrackScore", "ESC", "MeanFrontdEdX", "MichelCut", "ModifiedEavailable", "NIsoBlobs", "NonMIPClusFrac", "ODCalVisE", "TransverseGapScore"]
#cutsToDraw = ["E_lep"]
#cutsToDraw = ["NIsoBlobs", "DSCalVisE"]

for cut in cutsToDraw:
    signalHistoNames = [cut+"_selected_signal_reco"] 
    print(signalHistoNames)
    
    dataHistoNames=[]    
    for i in signalHistoNames:
        dataHistoNames.append(i.replace("selected_signal_reco", "data"))
        
        
    if cut == "Psi" or cut == "Lepton_Pt" or cut == "Etheta" or cut == "ProtonMomentum" or cut == "ProtonTheta":
        mcFile = ROOT.TFile.Open("/exp/minerva/data/users/cpernas/NuE_TKI/N-1_Feb_22/ExtraCuts/MC_Feb_23.root")
    else:        
        mcFile = ROOT.TFile.Open("/exp/minerva/data/users/cpernas/NuE_TKI/N-1_Feb_22/"+cut+"/MC_Feb_23.root")
    print(mcFile)

    for hist in range(len(signalHistoNames)):
        total = mcFile.Get(dataHistoNames[hist])
        signal = mcFile.Get(signalHistoNames[hist])
        nbins = signal.GetNbinsX()

        #calculate efficiency denominator (for this plot at least)
        signal_events_total = signal.Integral(1, nbins)
        cuts = []
        efficiency = []
        purity = []
        efficiency_x_purity = []
        significance = []
        bkg_fraction = []
        
        #greater than cuts
        if cut=="Afterpulsing" or cut=="E_lep" or cut=="EMLikeTrackScore" or cut=="NonMIPClusFrac" or cut=="TransverseGapScore":
            for i in range(1, nbins+1):
                print(i)
                cut_value = signal.GetBinLowEdge(i)
                
                #integrate from current bin i to end
                signal_events_selected = signal.Integral(i, nbins)
                all_events_selected = total.Integral(i, nbins)
                
                if all_events_selected <=0:
                    continue
                
                eff = signal_events_selected / signal_events_total
                pur = signal_events_selected / all_events_selected
                sig = signal_events_selected / math.sqrt(all_events_selected)
                bkg_frac = (all_events_selected - signal_events_selected)/all_events_selected
                
                #print("Evaluating cut at {0}, efficiency = {1:.3f} ; purity = {2:.3f} ; efficiency*purity = {3:.3f}".format(cut_value, eff, pur, eff*pur))
                cuts.append(cut_value)
                efficiency.append(eff)
                purity.append(pur)            
                efficiency_x_purity.append( eff * pur )        
                significance.append(sig)
                bkg_fraction.append(bkg_frac)
        elif cut=="DSCalVisE" or cut=="ODCalVisE" or cut=="ESC" or cut=="MeanFrontdEdX" or cut=="MichelCut" or cut=="ModifiedEavailable" or cut=="NIsoBlobs":
        #less than cuts
            for i in range(nbins, 0, -1):
                print(i)
                cut_value = signal.GetBinLowEdge(i)
                
                #integrate from current bin 1 to current bin i
                signal_events_selected = signal.Integral(1, i)
                all_events_selected = total.Integral(1, i)
                
                if all_events_selected <=0:
                    continue
                
                eff = signal_events_selected / signal_events_total
                pur = signal_events_selected / all_events_selected
                sig = signal_events_selected / math.sqrt(all_events_selected)
                bkg_frac = (all_events_selected - signal_events_selected)/all_events_selected
                
                #print("Evaluating cut at {0}, efficiency = {1:.3f} ; purity = {2:.3f} ; efficiency*purity = {3:.3f}".format(cut_value, eff, pur, eff*pur))
                cuts.append(cut_value)
                efficiency.append(eff)
                purity.append(pur)            
                efficiency_x_purity.append( eff * pur )
                significance.append(sig)
                bkg_fraction.append(bkg_frac)

        cuts_arr = array.array('d', cuts)
        eff_arr = array.array('d',efficiency)
        pur_arr = array.array('d', purity)
        prod_arr = array.array('d', efficiency_x_purity)
        sig_arr = array.array('d', significance)
        bkg_frac_arr = array.array('d', bkg_fraction)
        
        g_eff = ROOT.TGraph(len(cuts_arr), cuts_arr, eff_arr)
        g_pur = ROOT.TGraph(len(cuts_arr), cuts_arr, pur_arr)
        g_prod = ROOT.TGraph(len(cuts_arr), cuts_arr, prod_arr)
        g_sig =  ROOT.TGraph(len(cuts_arr), cuts_arr, sig_arr)
        g_bkg_frac = ROOT.TGraph(len(cuts_arr), cuts_arr, bkg_frac_arr)

        #to get variable name for the x axis
        head, sep, tail = signalHistoNames[hist].partition('_selected')
        outName = signalHistoNames[hist].replace("_selected_signal_reco", "")
        canvas = ROOT.TCanvas( 'canvas', outName, 0, 0, 2000, 1600 )                

        g_eff.SetLineColor(ROOT.kBlue)
        g_pur.SetLineColor(ROOT.kRed)
        g_prod.SetLineColor(ROOT.kBlack)
        g_sig.SetLineColor(ROOT.kMagenta+2)
        g_bkg_frac.SetLineColor(ROOT.kGreen+2)

        g_eff.SetLineWidth(3)
        g_pur.SetLineWidth(3)
        g_prod.SetLineWidth(3)
        g_sig.SetLineWidth(3)
        g_bkg_frac.SetLineWidth(3)
        
        g_eff.SetTitle("Efficiency and Purity Scan")
        g_eff.GetXaxis().SetTitle(head)
        g_eff.GetYaxis().SetTitle("Fraction")

        g_eff.Draw("AL")
        g_eff.GetYaxis().SetRangeUser(0.0, 1.0)
        
        g_pur.Draw("L same")
        g_prod.Draw("L same")
        #g_sig.Draw("AL")
        g_bkg_frac.Draw("L same")
        
        legend = ROOT.TLegend(0.60, 0.20, 0.85, 0.40)  # (x1,y1,x2,y2) in NDC
        legend.AddEntry(g_eff,  "Efficiency", "l")
        legend.AddEntry(g_pur,  "Purity", "l")
        legend.AddEntry(g_prod, "Efficiency #times Purity", "l")
        legend.AddEntry(g_sig, "Significance", "l")
        legend.AddEntry(g_bkg_frac, "Bkg Fraction", "l")
        
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.Draw()
        
        canvas.Update()
        outName = outName + ".png"
        canvas.SaveAs(outName)
        
        canvas.Delete()

    mcFile.Close()
