#!/bin/bash
# Open the FCSys package in Dymola on a Unix/Linux platform.

# Change to the base directory.
baseDir=`dirname $0`
cd $baseDir/../../../

# Start Dymola and run a Modelica script to initialize the package.
dymola ./load.mos &
