#!/bin/bash
# Reset the version, time, and last git commit (SHA) of a Modelica package.
#
# Assumptions:
# 1. The repository has the same name as the Modelica package.
#
# See http://semver.org/ and
# http://nvie.com/posts/a-successful-git-branching-model/ .
#
# Kevin Davies, 12/10/13

# Get the name of the package.
package="$(basename "$( pwd )" )"

# Reset the release information in the Modelica package.
cd $package*
# Version
sed -i s/version='"'[0-9A-Za-z.-]*'"',/version='""',/ package.mo
# Date modified
sed -i s/dateModified='"'[0-9\ Z:-]*'"'/dateModified='""'/ package.mo
# Abbreviated SHA of the last git commit
sed -i s/revisionID='"'[:\ 0-9A-Za-z]*'"'/revisionID='""'/ package.mo
