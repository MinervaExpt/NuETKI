#!/bin/bash
set -x
echo Start  `date`
echo Site:${GLIDEIN_ResourceName}
echo "the worker node is " `hostname` "OS: "  `uname -a`
echo "You are running as user `whoami`"

IFDH_OPTION=""

# set group based on the EXPERIMENT evnironment variable set by jobsub
GROUP=$EXPERIMENT

if [ -z $GROUP ]; then

# try to figure out what group the user is in
GROUP=`id -gn`
fi


case $GROUP in

e938)
SCRATCH_DIR=/pnfs/minerva/scratch/users
;;
minerva)
SCRATCH_DIR=/pnfs/minerva/scratch/users
;;
e875)
#need to check this
SCRATCH_DIR=/pnfs/minos/scratch/users
;;
minos)
#need to check this                                                                                                                                                                                                                           
SCRATCH_DIR=/pnfs/minos/scratch/users
;;
mars)
SCRATCH_DIR=""
;;
lbnemars)
SCRATCH_DIR="/lbne/data/lbnemars/users/"
;;
marslbne)
SCRATCH_DIR="/pnfs/dune/scratch/users/"
;;
marsmu2e)
SCRATCH_DIR="/pnfs/mu2e/scratch/users"
;;
marsgm2)
SCRATCH_DIR=""
;;
marsaccel)
SCRATCH_DIR=""
;;
larrand)
#pnfs/scene? probably not....
SCRATCH_DIR=""
;;
nova)
SCRATCH_DIR=/pnfs/nova/scratch/users
;;
t-962)
SCRATCH_DIR="/pnfs/argoneut/scratch/users"
;;
argoneut)
SCRATCH_DIR="/pnfs/argoneut/scratch/users"
;;
mu2e)
SCRATCH_DIR=/pnfs/mu2e/scratch/users
;;
microboone)
SCRATCH_DIR=/pnfs/uboone/scratch/users
;;
uboone)
SCRATCH_DIR=/pnfs/uboone/scratch/users
;;
lbne)
SCRATCH_DIR=/pnfs/lbne/scratch/users
;;
dune)
SCRATCH_DIR=/pnfs/dune/scratch/users
;;
seaquest|e906)
SCRATCH_DIR="/pnfs/e906/scratch/users"
#IFDH_OPTION="--force=cpn"
;;
coupp)
SCRATCH_DIR=""
;;
gm2)
# g-2 does not allow the gm2ana user to write to pnfs so we have to use blueArc for now
#SCRATCH_DIR=/gm2/data/users
SCRATCH_DIR=/pnfs/GM2/scratch/users
;; 
t-1034)
# lariat... no pnfs yet
SCRATCH_DIR=/pnfs/lariat/scratch/users
;;
lariat)
SCRATCH_DIR=/pnfs/lariat/scratch/users
;;
darkside)
SCRATCH_DIR="/pnfs/darkside/scratch/users"
;;
lar1nd)
SCRATCH_DIR="/pnfs/lar1nd/scratch/users"
;;
lsst)
SCRATCH_DIR="/pnfs/lsst/scratch/users"
;;
annie)
SCRATCH_DIR=""
;;
numix)
SCRATCH_DIR=""
;;
fermilab)
SCRATCH_DIR="/pnfs/fermilab/volatile"
#SCRATCH_DIR="/grid/data"
export CPN_LOCK_GROUP=gpcf
;;
des)
SCRATCH_DIR="/pnfs/des/scratch"
;;
esac

if [ -z $X509_CERT_DIR ] && [ ! -d /etc/grid-security/certificates ]; then
    if [ -f /cvmfs/oasis.opensciencegrid.org/mis/osg-wn-client/current/el7-x86_64/setup.sh ]; then
	source /cvmfs/oasis.opensciencegrid.org/mis/osg-wn-client/current/el7-x86_64/setup.sh || echo "Failure to run OASIS software setup script!"
    else
	echo "X509_CERT_DIR is not set, and the /etc/grid-security/certificates directory is not present. No guarantees ifdh, etc. will work!"
    fi
fi
if [ -z "`which globus-url-copy`" ] || [ -z "`which uberftp`" ]; then
 if [ -f /cvmfs/oasis.opensciencegrid.org/mis/osg-wn-client/current/el7-x86_64/setup.sh ]; then
     source /cvmfs/oasis.opensciencegrid.org/mis/osg-wn-client/current/el7-x86_64/setup.sh || echo "Failure to run OASIS software setup script!"
 else
     echo "globus-url-copy or uberftp (or both) is not in PATH, and the oasis CVMFS software repo does not appear to be present. No guarantees ifdh, etc. will work!"
 fi
fi

# A few commands to lookk around to see what kind of environment we are in. Should match the requested container.
voms-proxy-info --all

lsb_release -a

cat /etc/redhat-release

ps -ef

source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup


setup ifdhc -z /cvmfs/fermilab.opensciencegrid.org/products/common/db || { echo "error setting up ifdhc!" ; sleep 300 ; exit 1; }

# run the basic setup_dune.sh script if we are a DUNE job
if [ $EXPERIMENT == "dune" ]; then
    source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh || { echo "error running setup_dune.sh!" ; sleep 300 ; exit 1; }
fi 


echo "Here is the your environment in this job: " > job_output_${CLUSTER}.${PROCESS}.log
env >> job_output_${CLUSTER}.${PROCESS}.log

echo "group = $GROUP"

# If GRID_USER is not set for some reason, try to get it from the proxy
if [ -z ${GRID_USER} ]; then
GRID_USER=`basename $X509_USER_PROXY | cut -d "_" -f 2`
fi

echo "GRID_USER = `echo $GRID_USER`"



export GLIDEIN_ToDie=`condor_config_val GLIDEIN_ToDie`
echo "GLIDEIN_ToDie = $GLIDEIN_ToDie"


# let's try an ldd on ifdh

ldd `which ifdh`
sleeptime=$RANDOM
if [ -z "${sleeptime}" ] ; then sleeptime=4 ; fi
sleeptime=$(( (($sleeptime % 10) + 1)*60 ))
sleep $sleeptime

umask 002

if [ -z "$SCRATCH_DIR" ]; then
    echo "Invalid scratch directory, not copying back."
    echo "I am going to dump the log file to the main job stdout in this case."
    cat job_output_${CLUSTER}.${PROCESS}.log
else

# Very useful for debugging problems with copies
export IFDH_DEBUG=1
export IFDH_CP_MAXRETRIES=2
export IFDH_GRIDFTP_EXTRA="-st 1000"

    # first do lfdh ls to check if directory exists. We put a zero on the end because we only want 
    # to check that the directory exists; we don't care what's in it (i.e recursion depth of 0).
    ifdh ls ${SCRATCH_DIR}/$GRID_USER 0

    # A non-zero exit value probably means it doesn't exist yet, or does not have group write permission, 
    # so send a useful message saying that is probably the issue
    if [ $? -ne 0 ] && [ -z "$IFDH_OPTION" ] ; then	
	echo "Unable to read ${SCRATCH_DIR}/$GRID_USER. Make sure that you have created this directory and given it group write permission (chmod g+w ${SCRATCH_DIR}/$GRID_USER)."
	exit 74
    else
        # directory already exists, so let's copy
	ifdh cp -D $IFDH_OPTION job_output_${CLUSTER}.${PROCESS}.log ${SCRATCH_DIR}/${GRID_USER}
	if [ $? -ne 0 ]; then
	    echo "Error $? when copying to dCache scratch area!"
	    echo "If you created ${SCRATCH_DIR}/${GRID_USER} yourself,"
	    echo " make sure that it has group write permission."
	    exit 73
	fi
    fi
fi

echo "End `date`"

exit 0
