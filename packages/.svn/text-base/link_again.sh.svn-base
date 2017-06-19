#!/bin/sh

# put this file in the packages directory

export ROOTCOREBIN=packages/RootCore

export pkg=$1
export name=`basename $pkg`

test -d packages/$name/$name && rm -f $ROOTCOREBIN/include/$name && ln -sf ../../$name/$name $ROOTCOREBIN/include/
test -d packages/$name/python && rm -f $ROOTCOREBIN/python/$name && ln -sf ../../$name/python $ROOTCOREBIN/python/$name
test -d packages/$name/scripts && rm -f $ROOTCOREBIN/user_scripts/$name && ln -sf ../../$name/scripts $ROOTCOREBIN/user_scripts/$name
(test -d packages/$name/data && rm -f $ROOTCOREBIN/data/$name && ln -sf ../../$name/data $ROOTCOREBIN/data/$name) || (test -d packages/$name/share && rm -f $ROOTCOREBIN/data/$name && ln -sf ../../$name/share $ROOTCOREBIN/data/$name)
test -f packages/$name/StandAlone/lib${name}.so && rm -f $ROOTCOREBIN/lib/lib${name}.so && ln -sf ../../$name/StandAlone/lib${name}.so $ROOTCOREBIN/lib/

cd $ROOTCOREBIN/bin
if test -d ../../$name/bin
then
    for file in ../../$name/bin/*
    do
	test -f $file && rm -f $ROOTCOREBIN/bin/`basename $file` && ln -sf $file $ROOTCOREBIN/bin/
    done
fi
cd -

