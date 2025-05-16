#!/bin/sh
#export XRD_NETWORKSTACK=IPv4
#cd $CONDOR_DIR_INPUT
#tar -xvzf myareatar_tag_default_1745871086.tar.gz
#cd $CONDOR_DIR_INPUT
#export MINERVA_PREFIX=${CONDOR_DIR_INPUT}/exp/minerva/app/users/cpernas/MAT_AL9/opt
#cd $CONDOR_DIR_INPUT/exp/minerva/app/users/cpernas/MAT_AL9/
#cd $CONDOR_DIR_INPUT/exp/minerva/app/users/cpernas/MAT_AL9/
#source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
#spack load cmake@3.27.7
#spack load fife-utils@3.7.4
#spack load root@6.28.12
#export LD_LIBRARY_PATH=${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}
#cd $CONDOR_DIR_INPUT/exp/minerva/app/users/cpernas/MAT_AL9/opt/bin
#source setupGRID.sh
#export PYTHONPATH=/exp/minerva/app/users/cpernas/MAT_AL9/MAT-MINERvA/python:/exp/minerva/app/users/cpernas/MAT_AL9/MAT-MINERvA/python/PlotUtils
#cd $CONDOR_DIR_INPUT/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI
#pwd
#ls -a
#runEventLoop ${CONDOR_DIR_INPUT}/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/playlists/xrd_paths/p6/data_test_2.txt ${CONDOR_DIR_INPUT}/exp/minerva/app/users/cpernas/MAT_AL9/NuE_TKI/playlists/xrd_paths/p6/mc_test_2.txt
env
echo "Hello World" > testFile.txt
ls -al
ifdh cp ./testFile.txt /pnfs/minerva/scratch/users/cpernas
