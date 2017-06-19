#!/bin/bash

#WEB=http://ppewww.physics.gla.ac.uk/~dferreira/submission/
#wget $WEB/lastsubmit.tar.gz
DIR=/afs/cern.ch/user/$L/$CERN_USER/public
cp $DIR/lastsubmit.tar.gz .

tar xvfz lastsubmit.tar.gz

./grid_exec.sh $@

