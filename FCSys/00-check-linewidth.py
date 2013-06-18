#!/usr/bin/env python
# Check that Modelica comments will fit without wrapping in the LaTeX listing.
#
# The first argument is the directory, which defaults to the current directory.
#
# Kevin Davies, 3/26/2013

import re, glob, sys, os

# Maximum width of commented lines allowed without warning
width = 76

# Directory specification
if (len(sys.argv) > 1):
    path = sys.argv[1]
else:
    path = '.'

# Replace strings.
for fname in glob.glob(os.path.join(path, '*.mo')):

    if fname[len(path)+1:] not in ['WorkInProgress.mo']:
        # Read the source file.
        print "Checking " + fname + "..."
        src = open(fname, 'r')
        lines = src.readlines()
        src.close()

        for i, line in enumerate(lines):
            match = re.match(r' *//(.*)', line)
            if match:
                if len(match.group(1)) > width - 3:
                    print '%i: //%s'%(i, match.group(1))
        print ''
