#sample command to run:

#python SubmitJobsToGrid.py --stage eventLoop --playlist minervame1A --basedir /exp/minerva/app/users/cpernas/MAT_AL9/ --outdir /pnfs/minerva/persistent/users/cpernas/default_analysis_loc

# other options: --pdg, --
import os,sys,time
from optparse import OptionParser
import datetime

###########################################################
#  Script: JobSubmission python script for submitting
#          analysis jobs to the grid
###########################################################

# Write the command you used to run your analysis
def writeTransWarp(mywrapper,outdir):
    mywrapper.write(". runTransWarp_NEUTRINO_nu.sh lowQ2PiTuneOn ")
    mywrapper.write("\n")
    mywrapper.write("ifdh cp ./*.root "+outdir)

def writeEventLoop(mywrapper,outdir,playlist,pdg,filename):
    #mywrapper.write("./Event_Selection_Tracker_anyMode . "+playlist+" "+pdg+" "+str(filename))

    #Determine playlist from input...
    mywrapper.write("runEventLoop ${CONDOR_DIR_INPUT}/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/playlists/MAD_data_"+playlist+".txt ${CONDOR_DIR_INPUT}/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/playlists/MAD_mc_"+playlist+".txt |& tee CutSummary.txt")

    
    #one run of data and MC each, for testing
    #mywrapper.write("runEventLoop ${CONDOR_DIR_INPUT}/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/playlists/test_data.txt ${CONDOR_DIR_INPUT}/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/playlists/test_mc.txt |& tee CutSummary.txt")
    
    #same as above just different runs, for cross checking
    #mywrapper.write("runEventLoop ${CONDOR_DIR_INPUT}/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/playlists/xrd_paths/p6/data_test_2.txt ${CONDOR_DIR_INPUT}/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/playlists/xrd_paths/p6/mc_test_2.txt")

    #one quarter me1A, using xrd
    #mywrapper.write("runEventLoop ${CONDOR_DIR_INPUT}/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/playlists/p6/xrd_paths/data_one_fourth_me1A.txt ${CONDOR_DIR_INPUT}/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/playlists/p6/xrd_paths/mc_one_fourth_me1A.txt |& tee CutSummary.txt")
    
    #Full me1A, using xrd
    #mywrapper.write("runEventLoop ${CONDOR_DIR_INPUT}/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/playlists/xrd_paths/p6/mad_data_minervame1A.txt ${CONDOR_DIR_INPUT}/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/playlists/xrd_paths/p6/mad_mc_minervame1A.txt |& tee CutSummary.txt")

    #full me1A, direct pnfs paths with no xrd
    #mywrapper.write("runEventLoop ${CONDOR_DIR_INPUT}/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/playlists/p6/test_data.txt ${CONDOR_DIR_INPUT}/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/playlists/p6/test_MC.txt |& tee CutSummary.txt")

    
    mywrapper.write("\n")
    mywrapper.write("echo \"carlos - pwd, then printing out contents of current dir\"\n")
    mywrapper.write("pwd\n")
    mywrapper.write("ls -l\n")
    mywrapper.write("echo \"carlos - now trying to copy output to pnfs\"\n")
    mywrapper.write("ifdh cp ./*.root "+outdir)
    #Carlos added so I could get my cut summary
    mywrapper.write("\n")
    mywrapper.write("ifdh cp ./CutSummary.txt "+outdir)

def writeMigration(mywrapper,outdir,playlist,pdg,filename):
    mywrapper.write("./Migration_Tracker_anyMode . "+playlist+" "+pdg+" "+str(filename))
    mywrapper.write("\n")
    #mywrapper.write("./M . 6A\n")
    mywrapper.write("ifdh cp ./*.root "+outdir)

def writeEfficiency(mywrapper,outdir,playlist,pdg,filename):
    mywrapper.write("./Efficiency_Tracker_anyMode . "+playlist+" "+pdg+" "+str(filename))
    mywrapper.write("\n")
    mywrapper.write("ifdh cp ./*.root "+outdir)

def writeAll(mywrapper,outdir):
    mywrapper.write("./ES . 6J\n")
    mywrapper.write("ifdh cp ./*.root "+outdir)
    mywrapper.write("./M . 6J\n")
    mywrapper.write("ifdh cp ./*.root "+outdir)
    mywrapper.write("./Eff . 6J\n")
    mywrapper.write("ifdh cp ./*.root "+outdir)
 
def closureTest(mywrapper,outdir):
    mywrapper.write("./c\n")
    mywrapper.write("ifdh cp ./*.root "+outdir)

def writeSetups(mywrapper,basedir):
    # create topdir + source minerva products
    topdir = "$CONDOR_DIR_INPUT" + basedir
    mywrapper.write("cd "+topdir+"\n")
   # mywrapper.write("source /cvmfs/minerva.opensciencegrid.org/minerva/setup/setup_minerva_products.sh\n")
    # set environment variables

    #mywrapper.write("source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh\n")
    #mywrapper.write("spack load fife-utils@3.7.4\n")
    #mywrapper.write("spack cmake@3.27.7\n")
    #mywrapper.write("spack root@6.28.12\n")    
    #mywrapper.write("spack load gcc@12.2.0\n")


    mywrapper.write("source /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh\n")
    mywrapper.write("spack load root@6.28.12 arch=linux-almalinux9-x86_64_v3\n")
    mywrapper.write("spack load cmake@3.27.9%gcc@11.4.1 arch=linux-almalinux9-x86_64_v3\n")
    mywrapper.write("spack load ifdhc@2.8.0%gcc@12.2.0 arch=linux-almalinux9-x86_64_v3\n")
    #mywrapper.write("spack load gcc@12.2.0 arch=linux-almalinux9-x86_64_v3\n")


    #Do I need these? setup.sh and setupGRID.sh are VERY similar scripts and I run them both... why? -carlos
    #mywrapper.write("cd "+topdir+"\n")
    #mywrapper.write("cd opt/bin/\n")
    #mywrapper.write("source setup.sh\n")

    #Maria had these in her setup, because it looks like she has some additional compiled cxx stuff that she runs...
    #mywrapper.write("export NUKECCSRCROOT=${CONDOR_DIR_INPUT}/exp/minerva/app/users/cpernas/MAT_AL9/GitRepo-Analysis/NUKECCSRC\n")
    #mywrapper.write("export NUKECCSRC_ANA=${NUKECCSRCROOT}/ana_common\n")
    #mywrapper.write("export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CONDOR_DIR_INPUT}/exp/minerva/app/users/cpernas/MAT_AL9/GitRepo-Analysis/NUKECCSRC/ana_common/src\n")

    # source bash script
    mywrapper.write("cd "+topdir+"opt/bin\n")
    mywrapper.write("source setupGRID.sh\n")
    mywrapper.write("export LD_LIBRARY_PATH=${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}\n")
    #mywrapper.write("export PYTHONPATH=/exp/minerva/app/users/cpernas/MAT_AL9/MAT-MINERvA/python:/exp/minerva/app/users/cpernas/MAT_AL9/MAT-MINERvA/python/PlotUtils\n")
    #mywrapper.write("source /cvmfs/larsoft.opensciencegrid.org/spack-packages/opt/spack/linux-almalinux9-x86_64_v2/gcc-12.2.0/root-6.28.12-sfwfmqorvxttrxgfrfhoq5kplou2pddd/bin/thisroot.sh\n")
    # get to directory to run code
    mywrapper.write("cd "+topdir+"NuE_TKI\n")
    mywrapper.write("pwd\n")
    mywrapper.write("ls -a\n")
    #mywrapper.write("cd GitRepo-Analysis/ana/gridSubmission/\n")
    #mywrapper.write("echo Rint.Logon: ./rootlogon_grid.C > ./.rootrc\n")

# Define function to create tarball
def createTarball(outdir,tag,basedir):
    found = os.path.isfile("%s/myareatar_%s.tar.gz"%(outdir,tag))
    if(not found):
        cmd = "tar -czf /exp/minerva/app/users/$USER/myareatar_%s.tar.gz %s"%(tag,basedir)
        print("Making tar",cmd)
        os.system(cmd)
        cmd2 = "cp /exp/minerva/app/users/$USER/myareatar_%s.tar.gz %s/"%(tag,outdir)
        print("Copying tar",cmd2)
        os.system(cmd2)
# Define function to unpack tarball
def writeTarballProceedure(mywrapper,tag,basedir):
    print("I will be making the tarball upacking with this version") 
    print("Path is",basedir)
    #mywrapper.write("source /cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/setup.sh\n")
    mywrapper.write("export XRD_NETWORKSTACK=IPv4\n")
    #mywrapper.write("export UPS_OVERRIDE='-H Linux64bit+3.10-2.17'\n")
    mywrapper.write("cd $CONDOR_DIR_INPUT\n")
    mywrapper.write("tar -xvzf myareatar_%s.tar.gz\n"%tag)
    mywrapper.write("cd $CONDOR_DIR_INPUT\n")
    mywrapper.write("export MINERVA_PREFIX=${CONDOR_DIR_INPUT}/exp/minerva/app/users/cpernas/MAT_AL9/opt\n")
    writeSetups(mywrapper,basedir)

def writeOptions(parser):
    print("Now write options to the parser") 
    # Directory to write output
    parser.add_option('--outdir', dest='outdir', help='Directory to write output to', default = "/pnfs/minerva/persistent/users/"+_user_+"/default_analysis_loc/")
    parser.add_option('--basedir', dest='basedir', help='Base directory for making tarball', default = "NONE")
    parser.add_option('--stage', dest='stage', help='Process type', default="NONE")
    parser.add_option('--sample', dest='sample', help='Sample type', default="NONE")
    parser.add_option('--playlist', dest='playlist', help='Playlist type', default="NONE")
    ##########################################################################
    #  Options for making tar files....Basically you can make tarfiles 
    #######################################################################
    parser.add_option('--tag', dest='tag', help="Tag your release",default="tag_")
    parser.add_option('--mail',dest='mail',help="Want mail after job completion or failure",default=False,action="store_true")
    parser.add_option('--sametar',dest='sametar',help="Recycle the same tar file for jobsubmission",default=False,action="store_true")
    parser.add_option('--tarfilename',dest='tarfilename',help='Name of the tarfile you want to use',default="NA")
    parser.add_option('--notimestamp',dest='notimestamp',help='Flag to TURN OFF time stamp in the tag',default=False,action="store_true")
    parser.add_option('--filename',dest='filename',help='Name we want to use for output files',default="NONE")
    parser.add_option('--pdg',dest='pdg',help='The mode that we want to run the analysis in',default="NONE")

# Valid stages for the neutrinos (you can add more for your own analysis)
valid_stages=["eventLoop", "migration", "efficiency", "transwarp", "all", "closure"]
avail_playlists=["minervame6A", "minervame6B", "minervame6C", "minervame6D", "minervame6E", "minervame6F", "minervame6G", "minervame6H", "minervame6I", "minervame6J", "minervame5A", "minervame1A", "minervame1B", "minervame1C", "minervame1D", "minervame1E", "minervame1F", "minervame1G", "minervame1L", "minervame1M", "minervame1N", "minervame1O", "minervame1P"]
# Get user's name
_user_ = os.getenv("USER")
usage = "usage: %prog[opts]" 
parser = OptionParser(usage=usage)
writeOptions(parser)
(opts,args) = parser.parse_args()

# Now print some information for user's selected options
if opts.mail:
    print("############################################") 
    print("You have opted to Get mails once the jobs are finished of failed.") 
    print("#############################################") 

##############################################
#  TODO: Here I want to put the tag scheme...
##############################################
tag_name = str(opts.tag)
time_stamp = int(time.time())
if tag_name=="tag_":
    print("YOU DIDNT SPECIFY ANY ANY TAG...SO I WILL USE MY OWN TAGGING SCHEME") 
    # You can add in your own tagging scheme
    tag_name += "default_" 

# Add the time stamp to tag
if not opts.notimestamp:
    tag_name += str(time_stamp)

print("Tag for this version is ",tag_name)
print("******************************************************") 

# This is the output directory after the job is finished
output_dir = "$CONDOR_DIR_HISTS/" 

# Make outdir if not exist
if(not os.path.isdir(opts.outdir)):
    os.makedirs(opts.outdir)

#memory = 40000
memory = 10000
#memory = 3500
#disk = 10 #added in May20 b/c jobs going to 'held'

# Check the stage is valid or not
if opts.stage not in valid_stages:
    print(opts.stage,"Selected stage is not valid. Here are the valid options",valid_stages)
    sys.exit()

##############################################
#  Create wrapper
##############################################

wrapper_name = "%s_%s_wrapper_%s.sh"%(opts.stage,opts.playlist,tag_name)

mywrapper = open(wrapper_name,"w")
mywrapper.write("#!/bin/sh\n")
# Now create tarball
if (opts.sametar==False):
    createTarball(opts.outdir,tag_name,opts.basedir)
else:
    tarname = str(opts.tarfilename)
    if not os.path.exists(opts.outdir+opts.tarfilename):
        print("Tar File "+opts.outdir+opts.tarfilename+" doesn't Exist!") 
        sys.exit()
    #change the tag to the current one...
    cmd="cp "+opts.outdir+opts.tarfilename+" "+opts.outdir+"myareatar_"+tag_name+".tar.gz" 
    os.system(cmd)

# This will unpack the tarball we just made above
writeTarballProceedure(mywrapper,tag_name,opts.basedir)

# Now the add the command to run event loop 
if(opts.stage=="eventLoop"):
    writeEventLoop(mywrapper,opts.outdir,opts.playlist,opts.pdg,opts.filename)
elif (opts.stage=="migration"):
    writeMigration(mywrapper,opts.outdir,opts.playlist,opts.pdg,opts.filename)
elif (opts.stage=="efficiency"): 
    writeEfficiency(mywrapper,opts.outdir,opts.playlist,opts.pdg,opts.filename)
elif (opts.stage=="transwarp"):
    writeTransWarp(mywrapper,opts.outdir)
elif (opts.stage=="all"):
    writeAll(mywrapper,opts.outdir)
elif (opts.stage=="closure"):
    closureTest(mywrapper,opts.outdir)
mywrapper.close()

# Making the wrapper we just created readable, writable and executable by everyone
os.system("chmod 777 %s"%(wrapper_name))

# TODO: not sure if this is needed
configstring = "" 
# Get gcc release
#gccstring = "x86_64-slc7-gcc49-opt" 
# Now add the execute command
cmd = "" 
cmd += "jobsub_submit --group minerva " #Group of experiment
#cmd += "--debug "
#cmd += "--cmtconfig "+gccstring+" " #Setup minerva soft release built with minerva configuration

#cmd += "--OS sl7 " #Operating system #Not needed in SL7
#cmd += "--singularity-image /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest "
cmd += "--singularity-image /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest "
cmd += "-c has_avx2==True " #required architecture specification for grid node

if opts.mail:
    cmd += "-M " #this option to make decide if you want the mail or not
#cmd += "--subgroup=Nightly " #This is only for high priority jobs
#cmd += "--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC " 
cmd += "--role=Analysis " 
cmd += "--expected-lifetime 24h "  #normally 6h 
#cmd += "--cpu 0.25"
#cmd += "--memory "+str(memory)+"MB " 
#cmd += "--disk "+str(disk)+"GB " #added in May20 b/c jobs going to 'held'
cmd += configstring+" " #the environments for the tunes to bee applied
cmd += "-f "+opts.outdir+"/myareatar_"+tag_name+".tar.gz " 
#cmd += "-i /cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1"+" " 
cmd += "file://"+os.environ["PWD"]+"/"+wrapper_name
#print("jobsub command, not actually running for now though: ")
print(cmd)

#run the jobsub command
os.system(cmd)

print("Deleting the app area tar..... ") 
os.system("rm -rf /exp/minerva/app/users/$USER/myareatar_"+tag_name+".tar.gz")
print("Sleeping") 

time.sleep(1)
