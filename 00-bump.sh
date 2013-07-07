#!/bin/bash
# Checkout a release branch from master and apply a version number to the
# Modelica files.

# Get the version.
tag=`git describe --tags release`
versiona=${tag%-*}
versiona=${versiona%.*}
versiona=${versiona#v*}
versionb=${tag#v*.*.}
read -p "Enter the major and minor version (last was $versiona): " versiona
read -p "Enter the patch and optional pre-release string (last was $versionb): " versionb
#echo v$versiona.$versionb

# Get the name of the repo (assume that the package has the same name as the
# repo).
name="$(basename "$( pwd )" )"

# Create the branch.
git checkout -b "v$versiona.$versionb" master

# Go to the root of the repo.
cd "$(git rev-parse --show-toplevel)"

# Update the version number in various items.
# Library folder name
git mv $name* "$name $versiona"
# Load script
sed -i "s/$name[ 0-9.]*\/package.mo/$name $versiona\/package.mo/" load.mos
# Modelica version string
cd $name*
sed -i s/version='"'[0-9A-Za-z.]*'"',/version='"'$versiona.$versionb'"',/ package.mo

# Commit the updates.
git commit -am "Initial commit for $versiona.$versionb"

echo "Now on branch v$versiona.$versionb from master; files updated and commited."
