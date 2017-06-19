#!/bin/bash
#

# Change directory to the package directory.
cd ../

pkgdir=$PWD
depfile="$pkgdir/cmt/tagdeps.RootCore"
if [[ ! -f $depfile ]]; then
  echo "Error: Tag dependency file $depfile not found"
  exit 1
fi

# check if the username should be overridden
user=$USER
if [[ -n $1 ]]; then
  user=$1
fi

prefix="svn+ssh://$user@svn.cern.ch/reps/"
cd ../

# Checkout all packages depedencies
while read line; do
  # Take the directory name from the end of the svn chechout string
  dir=$(echo $line | perl -ne '@substr=split(/\s+/, $_); print "$substr[$#substr]";')
  dir=$(basename $dir)

  # Check if the directory exists
  if [[ ! -d "$dir" ]]; then

    # Checkout the package if it does not exist
    cmd="svn co "
    if [[ "$line" == "atlasoff"* || "$line" == "atlasgrp"* || "$line" == "atlasusr"* ]]; then
      cmd="$cmd$prefix$line"
    else
      cmd="$cmd$line"
    fi

    echo ">> $cmd"
    $cmd
    if [[ $? != 0 ]]; then
      exit 1
    fi
  fi
done < $depfile

# Configure RootCore if needed (check Makefile.common)
if [[ ! -f "RootCore/Makefile-common" ]]; then
  echo ">> RootCore ./configure"
  cd RootCore
  ./configure
  if [[ $? != 0 ]]; then
    exit 2
  fi
  cd ../
fi

# Setup RootCore source setup script
source RootCore/scripts/setup.sh

# Find all packages and dependencies
$ROOTCOREDIR/scripts/find_packages.sh

# Build all packages, including this one
$ROOTCOREDIR/scripts/compile.sh

#-----------------------------------------------------
# Particular to this package
