#!/bin/bash
# List the ** tags in the *.mo files (besides WorkInProgress.mo).
#
# Kevin Davies, 7/13/12

#grep "\*\*"  --color --exclude WorkInProgress.mo *.mo
grep "\*\*"  --color --exclude WorkInProgress.mo --exclude Systems.mo --exclude Blocks.mo *.mo
echo -n "Press enter to exit."
read answer
