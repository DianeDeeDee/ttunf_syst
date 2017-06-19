#!/bin/bash

export ROOTCOREDIR=`pwd`/packages/RootCore

export PATH=$ROOTCOREDIR/bin:$PATH
export LD_LIBRARY_PATH=fastjet-install/lib:$ROOTCOREDIR/lib:$ROOTSYS/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$ROOTCOREDIR/lib:$ROOTSYS/lib:$DYLD_LIBRARY_PATH
export PYTHONPATH=$ROOTCOREDIR/python:$PYTHONPATH

export LHAPATH=.

./preselect $@

