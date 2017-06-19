#!/bin/bash
#
# W. H. Bell
# A script to build all TopWorkingGroup D3PD analysis packages.
#

base_dir="../../"
if [[ -n $1 ]]; then
  base_dir=$1
fi

svn_user=$USER
if [[ -n $CERN_USER ]]; then
  svn_user=$CERN_USER
fi

share_dir=$PWD
package_file="$share_dir/packages.txt"

which root-config &> /dev/null
if [[ $? != 0 ]]; then
  echo "Error: please setup ROOT before running this script."
  exit 1
fi

if [[ ! -f "$package_file" ]]; then
  echo "Error: file $package_file could not be found."
  exit 2
fi

cd $base_dir
if [[ $? != 0 ]]; then 
  exit 3
fi

# Check if RootCore is already present or not
if [[ ! -d "RootCore" ]]; then

  # Check if a RootCore version is provided or not
  rootCore="svn+ssh://"$svn_user"@svn.cern.ch/reps/"$(grep "PhysicsAnalysis/D3PDTools/RootCore" $package_file)
  if [[ $? != 0 ]]; then
    echo "Error: RootCore tag could not be found in packages.txt"
    exit 4
  fi

  # Checkout RootCore
  svn co $rootCore RootCore
  if [[ $? != 0 ]]; then
    echo "Error: failed to checkout RootCore from $rootCore"
    exit 5
  fi
fi

echo $PWD

# If the RootCore environmental variables are set in this shell
# the build script will fail.  Therefore, unset them if they are set
if [[ -n $ROOTCOREDIR ]]; then
  unset ROOTCOREDIR
fi
if [[ -n $ROOTCOREBIN ]]; then
  unset ROOTCOREBIN
fi

# Checkout and build all of the packages in the package list.
#RootCore/scripts/build.sh $package_file
