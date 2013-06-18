#!/bin/bash
# List the work in progress in the *.mo files (besides WorkInProgress.mo).
#
# Kevin Davies, 7/13/2012

function find_item {
   #grep -c $item --color --exclude WorkInProgress.mo *.mo
   grep -c $1 --color --exclude WorkInProgress.mo --exclude Systems.mo --exclude Blocks.mo *.mo
   echo
   #grep $item -n --color --exclude WorkInProgress.mo *.mo
   grep $1 -n --color --exclude WorkInProgress.mo --exclude Systems.mo --exclude Blocks.mo *.mo
}

echo Starred items:
echo --------------
find_item "\*\*"

echo
echo TODO items:
echo -----------
find_item TODO

echo
echo UnderConstruction items:
echo ------------------------
find_item UnderConstruction

echo
read -p "Press [Enter] to exit."
