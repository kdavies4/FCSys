#!/bin/bash
# Add the version, time, and last git commit (SHA) to a Modelica package.
#
# Assumptions:
# 1. The repository has the same name as the Modelica package.
# 2. The branch is named with the version number (e.g., v1.0.0).
#
# See http://semver.org/ and
# http://nvie.com/posts/a-successful-git-branching-model/ .
#
# Kevin Davies, 7/7/13

# Get the name of the package.
package="$(basename "$( pwd )" )"

# Get the version.
version=`git rev-parse --abbrev-ref HEAD`
version=${version##v}

# Get the time stamp and SHA hash.
timestamp=`date -u +"%Y-%m-%d %H:%M:%SZ"`
hash=`git log --pretty=format:'%h' -n 1`

# Record the release information in the Modelica package.
# Modelica version string
cd "`echo $package*.*`"
sed -i s/version='"'[0-9A-Za-z.-]*'"',/version='"'"$version"'"',/ package.mo
# Date modified
sed -i s/dateModified='"'[0-9\ Z:-]*'"'/dateModified='"'"$timestamp"'"'/ package.mo
# Abbreviated SHA of the last git commit
sed -i s/revisionID='"'[:\ 0-9A-Za-z]*'"'/revisionID='"'"SHA: $hash"'"'/ package.mo

# Finish.
echo "Added release information (version, dateModified, and revisionID) for v$version."
