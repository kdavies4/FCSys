#!/bin/bash
# List the work in progress in the *.mo files (besides WorkInProgress.mo).
#
# Kevin Davies, 7/13/12

echo Starred items:
echo --------------
#grep -c "\*\*" --color --exclude WorkInProgress.mo *.mo
grep -c "\*\*" --color --exclude WorkInProgress.mo --exclude Systems.mo --exclude Blocks.mo *.mo
echo
#grep "\*\*" --color --exclude WorkInProgress.mo *.mo
grep "\*\*" --color --exclude WorkInProgress.mo --exclude Systems.mo --exclude Blocks.mo *.mo

echo
echo TODO items:
echo -----------
#grep -c TODO --color --exclude WorkInProgress.mo *.mo
grep -c TODO --color --exclude WorkInProgress.mo --exclude Systems.mo --exclude Blocks.mo *.mo
echo
#grep TODO --color --exclude WorkInProgress.mo *.mo
grep TODO --color --exclude WorkInProgress.mo --exclude Systems.mo --exclude Blocks.mo *.mo

echo
echo UnderConstruction items:
echo ------------------------
#grep -c "UnderConstruction;" --color --exclude WorkInProgress.mo *.mo
grep -c "UnderConstruction;" --color --exclude WorkInProgress.mo --exclude Systems.mo --exclude Blocks.mo *.mo
echo
grep "UnderConstruction;" --color --exclude WorkInProgress.mo *.mo
grep "UnderConstruction;" --color --exclude WorkInProgress.mo --exclude Systems.mo --exclude Blocks.mo *.mo

echo -n "Press enter to exit."
read answer
