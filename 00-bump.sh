#!/bin/bash
# Checkout a release branch from master and apply a version number to the
# Modelica files.
#
# The branch is named with the version number.
#
# Assumptions:
# 1. The repository has the same name as the Modelica package.
#
# See http://semver.org/ and
# http://nvie.com/posts/a-successful-git-branching-model/ .
#
# Kevin Davies, 7/7/13

# Get the version.
# TODO: Split the version number automatically (https://gist.github.com/pete-otaqui/4188238).
lasttag=`git describe --abbrev=0 --tags`
versiona=${lasttag%-*}
versiona=${versiona%.*}
versiona=${versiona#v*}
versionb=${lasttag#v*.*.}
read -p "Enter the major and minor version (last was $versiona): " versiona
read -p "Enter the patch and optional pre-release string (last was $versionb): " versionb
#echo v$versiona.$versionb

# Get the name of the package (see assumption 1).
package="$(basename "$( pwd )" )"
#echo $package

# Create the branch.
git checkout -b "v$versiona.$versionb" master

# Go to the root of the repo.
cd "$(git rev-parse --show-toplevel)"

# Update the version number in various items.
# Package folder
git mv "`echo $package*.*`" "$package $versiona"
# Load script
sed -i "s/$package[ 0-9.]*\/package.mo/$package $versiona\/package.mo/" load.mos
# Link to documentation
sed -i "s/$package[ 0-9.]*\/help/$package $versiona\/help/" Documentation.html
# Readme file
sed -i "s/($package[ 0-9.]*\//($package $versiona\//" README.md
# Python module
sed -i s/version='"'[0-9A-Za-z.]*'"',/version='"'$versiona.$versionb'"',/ "$package $versiona/Resources/Source/Python/setup.py"
# Modelica release information
./00-stamp.sh

# Finish.
git commit -am "Added release information for v$version"
echo "Now on branch v$versiona.$versionb from master."

