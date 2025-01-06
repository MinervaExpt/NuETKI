import sys
sys.argv.append( '-b' )
import ROOT
from ROOT import PlotUtils

import sys
from array import array
import math
import operator


daChain = ROOT.TChain("MasterAnaDev")
daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110036_Playlist.root")


daChain.SetBranchStatus("*", False)

daChain.SetBranchStatus("mc_FSPartPDG", True)
daChain.SetBranchStatus("mc_nFSPart", True)

#actually nvm these, im going to do this with a dictionary where the key is the partPDG code (maybe as a string), and the value is the # of occurrences. seems easier
particlePDGs = {}

for i, event in enumerate(daChain):
#    if i > 1000:
#        break
    for particle in event.mc_FSPartPDG:
        PDGcode = str(particle)
        if PDGcode not in particlePDGs.keys():
            particlePDGs[PDGcode] = 1
        else:
            particlePDGs[PDGcode] += 1

sorted_pdgs = sorted(particlePDGs.items(), key=operator.itemgetter(1))
#print sorted_pdgs
for i in xrange(len(sorted_pdgs)-1, 0, -1):
    print sorted_pdgs[i][0] + ": ", sorted_pdgs[i][1]
