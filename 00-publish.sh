#!/bin/bash
# Publish the current branch as a release version.
#
# This does the following:
# 1. Merges the branch into master, tags it, and pushes master to origin.
# 2. Merges the branch into development, removes the version information, and
#    pushes development to origin.
# 3. Pushes the web documentation (gh-pages) to origin.
#
# Assumptions:
# 1. The repository has the same name as the Modelica package.
# 2. The branch is named with the version number (e.g., v1.0.0).
#
# Kevin Davies, 7/8/13


# Get the name of the package (see assumption 1).
package="$(basename "$( pwd )" )"

# Get the branch and version number (see assumption 2).
branch=`git rev-parse --abbrev-ref HEAD`
version=$branch


## Handle the master branch.
git checkout master
git merge --no-ff $branch
git tag $version
git push --tags origin master


## Handle the development branch.
git checkout $branch
git checkout -b $branch-temp # Temporary branch to reset the version info.

# Reset the version information.
cd $name*
# Package folder
git mv $package* $package
# Load script
sed -i "s/$package[ 0-9.]*\/package.mo/$package" load.mos
# Modelica version string
cd $package
sed -i s/version='"'[0-9A-Za-z.]*'"',/version='""',/ package.mo
# Date modified
sed -i s/dateModified='"'[0-9\ Z:-]*'"'/dateModified='""'/ package.mo
# Abbreviated SHA of the last git commit
sed -i s/revisionID='"'[:\ 0-9A-Za-z]*'"'/revisionID='""'/ package.mo

git commit "Reset version info after merging $version"
git checkout development
git merge --no-ff $branch-temp
git push origin development


## Finish.
git push origin gh-pages
echo "Published release $version."
