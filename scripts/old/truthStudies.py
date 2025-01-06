import ROOT
import math


def CalcApothem(x,y):
    x=abs(x)
    y=abs(y)
    if ( x == 0 or y/x > 1/math.sqrt(3)):
        return (y+x/math.sqrt(3))/2*math.sqrt(3)
    else:
        return x

#Determines if an event is signal based on truth quantities
#Right now this means incoming electron neutrino (currently excludes antineutrino), charged current, and with true vtx within the tracker & apothem cut
def isTrueSignal(event, true_apothem):
    if event.mc_vtx[2] > 5980 and event.mc_vtx[2] < 8422 and true_apothem < 850 and event.mc_incoming==12 and event.mc_current==1:
        return True
    else:
        return False

def hasTrackableFSProton(event):
    hasTrackableProton = False
    for i, pdg in enumerate(event.mc_FSPartPDG):
        part = abs(pdg)
        if part==2212 and event.mc_FSPartE[i]>1038: #1038 - 938 -> 100 MeV of kinetic energy. Just a guess, tune this.
            hasTrackableProton=True

    return hasTrackableProton
    
def hasFSProton(event):
    hasProton = False
    for pdg in event.mc_FSPartPDG:
        part = abs(pdg)
        if part==2212:
            hasProton=True
    return hasProton


def hasFSMeson(event):
    hasMeson = False
    for pdg in (event.mc_FSPartPDG):
        part = abs(pdg)
        if part==211 or part==111 or part==321 or part==311 or part==130:
            hasMeson = True
    return hasMeson

#basic, self explanatory cuts that I don't really need to look at specific histograms for (zrange in tracker, within 850cm apothem and deadtime cut
def passesPreCuts(event, apothem):
    if event.vtx[2] > 5980 and event.vtx[2] < 8422 and apothem < 850 and event.phys_n_dead_discr_pair_upstream_prim_track_proj==0:
        return True
    else:
        return False



#mc_file = ROOT.TFile( filename, 'read' )
#mc_tree = mc_file.Get('NukeCC')
#metaTree = mc_file.Get("Meta")


#Chaining together all subruns for one run of the special sample, no NukeCC selection -> SHOULD HAVE nuE's
print "hope this works? "
#daChain = ROOT.TChain("Truth")
daChain = ROOT.TChain("MasterAnaDev")

daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110030_Playlist.root")
#daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110031_Playlist.root")
#daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110032_Playlist.root")
#daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110033_Playlist.root")
#daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110034_Playlist.root")
#daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110035_Playlist.root")
#daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110036_Playlist.root")
#daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110037_Playlist.root")
#daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110038_Playlist.root")
#daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110039_Playlist.root")
#daChain.Add("/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110040_Playlist.root")
#/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110038_Playlist.root
#/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110039_Playlist.root
#/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110040_Playlist.root

print "Added ", daChain.GetNtrees(), " trees, and ", daChain.GetEntries(), " total events"

#entries = mc_tree.GetEntries()
#metaTree.GetEntry(0)
#print metaTree.POT_Used, " POT_Used, ", entries, " events total"


daChain.SetBranchStatus("*", False)

#Truth quantities
daChain.SetBranchStatus("mc_vtx", True)
daChain.SetBranchStatus("mc_incoming", True)
daChain.SetBranchStatus("mc_targetZ", True)
daChain.SetBranchStatus("mc_current", True)
daChain.SetBranchStatus("mc_FSPartPDG", True)
daChain.SetBranchStatus("mc_FSPartE", True)
daChain.SetBranchStatus("mc_intType", True)
daChain.SetBranchStatus("mc_Q2", True)
daChain.SetBranchStatus("mc_w", True)
daChain.SetBranchStatus("mc_incomingE", True)
daChain.SetBranchStatus("mc_run", True)
daChain.SetBranchStatus("mc_subrun", True)
daChain.SetBranchStatus("mc_nthEvtInFile", True)

#reco branches
daChain.SetBranchStatus("MasterAnaDev_proton_endPointZ", True)
daChain.SetBranchStatus("HasNoBackExitingTracks", True)




NuMu_RecoProtons = 0
NuMu_TrueProtons = 0
NuE_RecoProtons = 0
NuE_TrueProtons = 0
RecoNue_RecoProtons = 0

for i, event in enumerate(daChain):
    #apothem = CalcApothem(event.mc_vtx[0], event.mc_vtx[1])
    if i%5000 == 0:
        print i

    if event.mc_incoming==14 and event.MasterAnaDev_proton_endPointZ > 0:
        NuMu_RecoProtons+=1
    if event.mc_incoming==14 and hasTrackableFSProton(event):
        NuMu_TrueProtons+=1
    if event.mc_incoming==12 and event.MasterAnaDev_proton_endPointZ > 0:
        NuE_RecoProtons+=1
        #print "true nu_e, reco'd proton - run: ", event.mc_run, ", subrun: ", event.mc_subrun, " nthEvtInFile: ", event.mc_nthEvtInFile
    #if event.HasNoBackExitingTracks==1 and event.MasterAnaDev_proton_endPointZ > 0 and event.mc_incoming==12:
    #RecoNue_RecoProtons+=1
    if event.mc_incoming==12 and hasTrackableFSProton(event):
        NuE_TrueProtons+=1

    #if i>2000:
        #break


print("numu reco proton events: ", NuMu_RecoProtons)
print("numu true proton events: ", NuMu_TrueProtons)
print("true nue reco proton events: ", NuE_RecoProtons)
#print("\"reco\" nue reco proton events: ", RecoNue_RecoProtons)
print("true nue true proton events: ", NuE_TrueProtons)
