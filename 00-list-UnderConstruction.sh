#!/bin/bash
# List the occurences of the "UnderConstruction" string in the *.mo files
# (besides WorkInProgress.mo).
#
# Kevin Davies, 7/13/12

#grep -c "\*\*" --color --exclude WorkInProgress.mo *.mo
grep -c "UnderConstruction;" --color --exclude WorkInProgress.mo --exclude Systems.mo --exclude Blocks.mo *.mo
echo
#grep "\*\*" --color --exclude WorkInProgress.mo *.mo
grep "UnderConstruction;" --color --exclude WorkInProgress.mo --exclude Systems.mo --exclude Blocks.mo *.mo
echo -n "Press enter to exit."
read answer
