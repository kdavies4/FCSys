#!/bin/bash
# List the tracked files in this branch.
#
# Created by Kevin Davies, 10/21/2012

branch=`git symbolic-ref HEAD 2>/dev/null | cut -d"/" -f 3`
echo These files are tracked on the $branch branch of the repository:
git ls-tree -r --name-only $branch
echo
echo -n Press enter to exit.
read answer
