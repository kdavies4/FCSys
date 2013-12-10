#!/bin/bash
# Count the number of classes in the Modelica files.
#
# Kevin Davies, 7/26/2013

f=`find . -type f \( -iname "*.mo" ! -iname "WorkInProgress.mo" ! -iname "Figures.mo" \)`

echo There are `cat $f | wc -l` lines of code.
wc -l $f
echo
echo There are `cat $f | grep -c "model [A-Za-z0-9_'+-]* "` models.
grep -c "model [A-Za-z0-9_'+-]* " $f
echo
echo There are `cat $f | grep -c "block [A-Za-z0-9_'+-]* "` blocks.
grep -c "block [A-Za-z0-9_'+-]* " $f
echo
echo There are `cat $f | grep -c "function [A-Za-z0-9_'+-]* "` functions.
grep -c "function [A-Za-z0-9_'+-]* " $f
echo
echo There are `cat $f | grep -c "record [A-Za-z0-9_'+-]* "` records.
grep -c "record [A-Za-z0-9_'+-]* " $f
echo
echo There are `cat $f | grep -c "package [A-Za-z0-9_'+-]* "` packages.
grep -c "package [A-Za-z0-9_'+-]* " $f
echo
