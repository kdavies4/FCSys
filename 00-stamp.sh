#!/bin/bash
# Add the time and last git commit (SHA) to a Modelica package.

# Get the name of the repo (assume that the package has the same name as the
# repo).
name="$(basename "$( pwd )" )"

# Record the release information in the Modelica package.
cd $name*
timestamp=`date -u +"%Y-%m-%d %H:%M:%SZ"`
hash=`git log --pretty=format:'%h' -n 1`
# Date modified
sed -i s/dateModified='"'[0-9\ Z:-]*'"'/dateModified='"'"$timestamp"'"'/ package.mo
# Abbreviated SHA of the last git commit
sed -i s/revisionID='"'[:\ 0-9A-Za-z]*'"'/revisionID='"'"SHA: $hash"'"'/ package.mo

# Tag and commit (assume the name of the branch is the version number).
version=`git rev-parse --abbrev-ref HEAD`
git commit -am "Release stamp for $version"
git tag $version
#git push --tags origin $version

echo "Tagged release $version."
echo "Do 'git push --tags origin $version' once all ready."
