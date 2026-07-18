import ROOT
from ROOT import PlotUtils
import sys

chi2SummaryDir = "Chi2_Iteration_Dists"
chi2SummaryName = "h_chi2_modelData_trueData_iter_chi2"
medianHistName = "h_median_chi2_modelData_trueData_iter_chi2"
meanChi2ProfileName = "m_avg_chi2_modelData_trueData_iter_chi2_truncated"
ratioDir = "Ratio_Unfolded_True_Histograms"
ratio0Name = "Ratio_modelData_trueData_Stat_0_Iter_"

lineWidth = 3
iterChosen = float(sys.argv[2])
matTag = sys.argv[3]+" "+sys.argv[4]

baseName = sys.argv[1]

if baseName.find(".root") != -1:
  print ("GOTTA NAME YOUR OUTPUT FILE")
  exit

can = ROOT.TCanvas("chi2")

for fileName in sys.argv[5:]:

  ROOT.gROOT.SetBatch(True)
  
  myFile = ROOT.TFile.Open(fileName)

  #Try to infer a useful universe name from the file name
  univName = fileName[fileName.find("MnvTune") + len("MnvTune"):(fileName.find("FVregion")-1)]
  #if fileName.find("SuSA") != -1: #The SuSA warp is a stand-alone CV, so it needs special treatment
    #univName = "SuSA"

  outName = baseName+"_Warp"+univName

  spread = myFile.Get(chi2SummaryDir).Get(chi2SummaryName)

  #Infer number of degrees of freedom from y axis title
  axisTitle = spread.GetYaxis().GetTitle()
  yNDF = int(axisTitle[axisTitle.find("ndf=") + 4:axisTitle.find(")")])

  spread.SetTitle("MC: " + univName)
  spread.GetXaxis().SetTitle("No. of Iterations")
  spread.SetTitleOffset(0.785, "X")
  spread.SetTitleSize(0.05,"X")
  spread.SetTitleOffset(0.835, "Y")
  spread.SetTitleSize(0.05,"XY")
  spread.Draw("colz")
  can.Print(outName + "_fullRange.png")
  can.Print(outName + "_fullRange.pdf")
  can.Print(outName + "_fullRange.C")
  spread.GetXaxis().SetRangeUser(1.0,30.0)
  spread.GetYaxis().SetRangeUser(0.0,70.0)
  spread.Draw("colz")

  profile = myFile.Get(chi2SummaryDir).Get(meanChi2ProfileName)
  profile.SetTitle("Mean Chi2")
  profile.SetLineWidth(lineWidth)
  profile.SetLineColor(ROOT.kBlue)
  profile.SetMarkerStyle(0)
  profile.Draw("SAME")

  median = myFile.Get(chi2SummaryDir).Get(medianHistName)
  median.SetTitle("Median Chi2")
  median.SetLineWidth(lineWidth)
  median.SetLineColor(ROOT.kBlack)
  median.Draw("HIST SAME")

  #Draw lines at number of degrees of freedom and 2x NDF
  #ndfLine = ROOT.TLine(1, yNDF, spread.GetXaxis().GetXmax(), yNDF)
  ndfLine = ROOT.TLine(1, yNDF, 30, yNDF)
  ndfLine.SetLineWidth(lineWidth)
  ndfLine.SetLineStyle(ROOT.kDashed)
  ndfLine.Draw()

  #doubleNDFLine = ROOT.TLine(1, 2*yNDF, spread.GetXaxis().GetXmax(), 2*yNDF)
  doubleNDFLine = ROOT.TLine(1, 2*yNDF, 30, 2*yNDF)
  doubleNDFLine.SetLineColor(ROOT.kRed)
  doubleNDFLine.SetLineWidth(lineWidth)
  doubleNDFLine.SetLineStyle(ROOT.kDashed)
  doubleNDFLine.Draw()

  #Draw a line at the chosen number of iterations.
  iterLine = ROOT.TLine(iterChosen + 0.5, 0, iterChosen + 0.5, 70.0)
  iterLine.SetLineWidth(lineWidth)
  iterLine.SetLineStyle(ROOT.kDotted)
  iterLine.Draw()

  #Make a custom legend because I don't want to include the 2D histogram
  #while I must Draw() it first to set the right axis limits.
  leg = ROOT.TLegend(0.6, 0.6, 0.9, 0.9)
  leg.AddEntry(profile)
  leg.AddEntry(median)
  leg.AddEntry(ndfLine, "Number of Bins", "l")
  leg.AddEntry(doubleNDFLine, "2x Number of Bins", "l")
  leg.AddEntry(iterLine, str(iterChosen) + " iterations", "l")
  leg.Draw()

  #Add a tag for the name
  nameTagTitle = "#it{"+matTag+"}"
  nameTag = ROOT.TLatex(0.65,0.525,nameTagTitle)
  nameTag.SetNDC()
  nameTag.SetTextColor(ROOT.kRed)
  nameTag.SetTextFont(43)
  nameTag.SetTextFont(32)
  nameTag.SetLineWidth(2)
  nameTag.Draw()

  WIPTag = ROOT.TLatex(0.175,0.83,"#it{MINER#nuA Work in Progress}")
  WIPTag.SetNDC()
  WIPTag.SetTextColor(ROOT.kRed)
  WIPTag.SetTextFont(43)
  WIPTag.SetTextFont(32)
  WIPTag.SetLineWidth(2)
  WIPTag.Draw()

  can.Print(outName + "_ZoomedIn.png")
  can.Print(outName + "_ZoomedIn.pdf")
  can.Print(outName + "_ZoomedIn.C")
  
  iterations = [str(it) for it in range(1,11)]+[str(10*it) for it in range(2,3)]+[str(50*it) for it in range(1,3)]
  colorsList = [ROOT.kBlack,ROOT.kRed,ROOT.kRed-6,ROOT.kGreen,ROOT.kGreen-6,ROOT.kBlue,ROOT.kBlue-6,ROOT.kOrange,ROOT.kOrange-6,ROOT.kTeal,ROOT.kTeal-6,ROOT.kViolet,ROOT.kViolet-6]
  leg2 = ROOT.TLegend(0.6, 0.6, 0.9, 0.9)

  hists=[]
  counter=0
  for it in iterations:
    hist = myFile.Get(ratioDir).Get(ratio0Name+it)
    hist.SetTitle("Unfolded Data/True Data No Stat. Variation")
    hist.GetYaxis().SetRangeUser(1-1e-9,1+1e-9)
    hist.GetYaxis().SetTitle("Unfolded Data/True Data")
    hist.GetXaxis().SetTitle("Muon Transverse Momentum [GeV/c]")
    hist.SetLineColor(colorsList[counter])
    hists.append(hist)
    counter+=1

  ROOT.gStyle.SetOptStat(0)
  counter=0
  for hist in hists:
    if counter==0:
      hist.Draw()
    else:
      hist.Draw("same")
    leg2.AddEntry(hist,"Iteration "+iterations[counter])
    counter+=1
  leg2.Draw()
  nameTag.Draw()
  WIPTag.Draw()
  can.Print(outName+"_Ratio_Stat_0.png")
  can.Print(outName+"_Ratio_Stat_0.pdf")
  can.Print(outName+"_Ratio_Stat_0.C")
