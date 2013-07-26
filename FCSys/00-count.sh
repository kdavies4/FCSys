#!/bin/bash
# Count the number of classes in the Modelica files.
#
# Kevin Davies, 7/26/2013

echo There are `cat *.mo | wc -l` lines of code.
wc -l *.mo
echo
echo There are `cat *.mo | grep -c "model [A-Za-z0-9_'+-]* "` models.
grep -c "model [A-Za-z0-9_'+-]* " *.mo
echo
echo There are `cat *.mo | grep -c "block [A-Za-z0-9_'+-]* "` blocks.
grep -c "block [A-Za-z0-9_'+-]* " *.mo
echo
echo There are `cat *.mo | grep -c "function [A-Za-z0-9_'+-]* "` functions.
grep -c "function [A-Za-z0-9_'+-]* " *.mo
echo
echo There are `cat *.mo | grep -c "record [A-Za-z0-9_'+-]* "` records.
grep -c "record [A-Za-z0-9_'+-]* " *.mo
echo
echo There are `cat *.mo | grep -c "package [A-Za-z0-9_'+-]* "` packages.
grep -c "package [A-Za-z0-9_'+-]* " *.mo
echo
