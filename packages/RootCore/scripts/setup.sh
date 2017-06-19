if test "x$ROOTSYS" = "x"
then
    export ROOTSYS=`root-config --prefix`
    test $? -ne 0 && echo no valid root version found, please set up root && return 1
fi

test "x$ROOTSYS" != "x/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/5.34.19-x86_64-slc6-gcc4.7" && echo ROOTSYS=$ROOTSYS not /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/5.34.19-x86_64-slc6-gcc4.7 && echo please set up correct root version or follow the instructions at && echo \ \ https://twiki.cern.ch/twiki/bin/viewauth/Atlas/RootCore#If_Something_Goes_Wrong && return 1

export ROOTCOREDIR=/afs/cern.ch/work/d/dshoaleh/private/Analysis/ttunf_syst/packages/RootCore
test "x$ROOTCOREBIN" = "x" && export ROOTCOREBIN=/afs/cern.ch/work/d/dshoaleh/private/Analysis/ttunf_syst/packages/RootCore
test \! -d "$ROOTCOREBIN" && echo output directory not found, please create it && echo "  mkdir $ROOTCOREBIN" && return 1

export PATH=$ROOTCOREBIN/bin:$PATH
export LD_LIBRARY_PATH=$ROOTCOREBIN/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$ROOTCOREBIN/lib:$DYLD_LIBRARY_PATH
export PYTHONPATH=$ROOTCOREBIN/python:$PYTHONPATH

mkdir -p $ROOTCOREBIN/bin
test -f $ROOTCOREBIN/bin/rc || ln -sf $ROOTCOREDIR/scripts/rc.sh $ROOTCOREBIN/bin/rc
true
