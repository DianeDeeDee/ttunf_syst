#!/bin/bash

cd packages/RootCore
./configure
cd -

export ROOTCOREDIR=`pwd`/packages/RootCore
source ${ROOTCOREDIR}/scripts/setup.sh

export LD_LIBRARY_PATH=$ROOTCOREDIR/lib:$ROOTSYS/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$ROOTSYS/lib:$DYLD_LIBRARY_PATH

cd packages
RootCore/scripts/find_packages.sh
cd ..

export PYTHONPATH=.:$PYTHONPATH
export PATH=$ROOTSYS/bin:$PATH
