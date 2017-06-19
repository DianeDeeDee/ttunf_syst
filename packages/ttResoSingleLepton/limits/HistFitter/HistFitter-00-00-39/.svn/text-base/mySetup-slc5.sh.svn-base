#. /afs/cern.ch/sw/lcg/external/gcc/4.6/x86_64-slc5/setup.sh
#export PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc5-gcc46-opt/bin:$PATH
#export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc5-gcc46-opt/lib:$LD_LIBRARY_PATH
#cd /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.08/x86_64-slc5-gcc46-opt/root
#source bin/thisroot.sh
#echo "ROOTSYS: $ROOTSYS"
#cd -
#export PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc6-gcc47-opt/bin:$PATH
#export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc6-gcc47-opt/lib:$LD_LIBRARY_PATH

###setupATLAS
### source ${ATLAS_LOCAL_ROOT_BASE}/packageSetups/atlasLocalROOTSetup.sh --rootVersion=5.34.08-x86_64-slc5-gcc4.3

setupATLAS
#export  ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
echo "sourcing atlasLocalSetup.sh"
#source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

#localSetupROOT -h
#localSetupROOT --rootVer=5.34.08-x86_64-slc5-gcc4.3

if [ ! $ROOTSYS ]; then
    echo "setting up root"
    ALRB_SKIP_XDR=1
    localSetupROOT --rootVer=5.34.14-x86_64-slc5-gcc4.3
    
    #pushd /a/data/xenia/users/altheim/root_slc5; . bin/thisroot.sh; popd

    
else
    which root
fi

#localSetupROOT --rootVer=5.34.11-x86_64-slc5-gcc4.7


echo "sourcing setup5.sh"
source setup5.sh

