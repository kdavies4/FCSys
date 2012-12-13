#!/usr/bin/env python
# Make custom replacements in files.
#
# The first argument is the directory.
#
# Created by Kevin Davies, 5/30/12

import re
import glob
import sys
import os

## Settings

# Replacement pairs
rpls = [# Remove extra spacing.
        (r' +\n', r'\n'),
        (r'\n\n\n+', r'\n\n'),
        (r' +<br>', r'<br>'),
        (r'<br><br>(<br>)+', r'<br><br>'),
        # Use shortcuts for Units and Quantities.
        (r'FCSys\.Quantities\.', r'Q.'),
        (r'FCSys\.Units\.', r'U.'),
        # Sometimes Dymola adds a useless import.
        (r'import FCSys;\n', ''),
        # No empty line before "end x;".
        (r'\n(\n +end )([^; ]+);', r'\1\2;'),
        # Two spaces after "Note:" and "TODO:".
        (r'Note: +([^ ])', r'Note:  \1'),
        (r'TODO: +([^ ])', r'TODO:  \1'),
        # Use lowercase "e" for engineering notation.
        (r'([0-9]+)E(-?)([0-9]+)', r'\1e\2\3'),
        # Don't use a "+" for positive powers of 10 in engineering notation.
        (r'([0-9]+)e+?([0-9]+)', r'\1e\2'),
        # Use some contractions in comments,
        (r'(// .*)do not', r"\1don't"),
        (r'(// .*)does not', r"\1doesn't"),
        (r'(// .*)is not', r"\1isn't"),
        (r'(// .*)will not', r"\1won't"),
        (r'(// .*)cannot', r"\1can't"),
        (r'(// .*)there is', r"\1there's"),
        # but not others.
        (r"(// .*)it's", r'\1it is'),
       ]

# Directory specification
if (len(sys.argv) > 1):
    directory = sys.argv[1]
else:
    directory = ''

## Procedure
# Method to remove non-ASCII characters.
# See
# http://stackoverflow.com/questions/1342000/how-to-replace-non-ascii-characters-in-string.
def removeNonAscii(s): return "".join(i for i in s if ord(i)<128)

# Compile the regular expressions.
for i, rpl in enumerate(rpls):
    rpls[i] = (re.compile(rpl[0]), rpl[1])

# Replace strings.
for fname in glob.glob(os.path.join(directory, '*.mo')):
    # Read the source file.
    print("Processing %s ... " % fname)
    src = open(fname, 'r')
    text = src.read()
    src.close()
    #text = removeNonAscii(text)
    # **Add a method to warn about non-ASCII characters before removing them.

    # Make the replacements.
    for rpl in rpls:
        text = rpl[0].sub(rpl[1], text)

    nonASCII = "".join(i for i in text if ord(i)>=128)
    if len(nonASCII) > 0:
        print("Warning:  It contains these non-ASCII characters: " + nonASCII)

    # Re-write the file.
    src = open(fname, 'w')
    src.write(text)
    src.close()
