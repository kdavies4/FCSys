#!/bin/bash
# Release the current branch.
#
# This does the following:
# 1. Pushes the web documentation (gh-pages) to origin.
# 2. Merges the branch into master, tags it, and pushes master to origin.
# 3. Merges the branch into development after removing the version information.
#
# Assumptions:
# 1. The repository has the same name as the Modelica package.
# 2. The branch is named with the version number (e.g., v1.0.0).
#
# See http://semver.org/ and
# http://nvie.com/posts/a-successful-git-branching-model/ .
#
# Kevin Davies, 7/8/13


# Get the name of the package (see assumption 1).
package="$(basename "$( pwd )" )"
#echo $package

# Get the branch and version number (see assumption 2).
branch=`git rev-parse --abbrev-ref HEAD`
version=$branch
#echo $version

# Go to the root of the repo.
cd "$(git rev-parse --show-toplevel)"


## Push the web documentation (action #1).
git push origin gh-pages


## Handle the master branch (action #2).
git checkout master
git merge --no-ff $branch
git tag -a $version
git push --tags origin master


## Handle the development branch (action #3).
git checkout $branch
git checkout -b $branch-temp # Temporary branch to reset the version info.

# Reset the version number in various items.
# Package folder
git mv "`echo $package*.*`" $package
# Load script
sed -i "s/$package[ 0-9.]*\/package.mo/$package\/package.mo" load.mos
# Python module
sed -i s/version='"'[0-9A-Za-z.]*'"',/version='""',/ $package*/Resources/Source/Python/setup.py
# Modelica release information
./00-reset.sh


## Finish.
git commit "Reset version info after merging $version"
git checkout development
git merge --no-ff $branch-temp
echo "Released $version.  Now on the development branch."
