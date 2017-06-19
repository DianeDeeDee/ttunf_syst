#!/bin/bash

source setup_rootcore.sh

cd packages
RootCore/scripts/clean.sh
cd -
packages/RootCore/scripts/compile.sh
pkgfile=packages/RootCore/packages
for pkg in `cat $pkgfile`
do
    name=`basename $pkg`
    packages/link_again.sh $pkg
done

