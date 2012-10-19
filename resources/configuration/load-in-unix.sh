#!/bin/bash
# Open the FCSys package in Dymola on a Unix/Linux platform.
#
# Created by Kevin Davies, 26-Jun-10

# Change to the FCSys base directory.
baseDir=`dirname $0`
cd $baseDir/../../

# Start Dymola and run a Modelica script to initialize the package.
dymola ./load.mos &
