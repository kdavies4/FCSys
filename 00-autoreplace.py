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
rpls = [# Remove empty annotation tags.
        (r'\n? *Diagram\(graphics\), *', ' '),
        (r',\n? *Diagram\(graphics\)', ''),
        (r'\n? *Icon\(graphics\), *', ' '),
        (r',\n? *Icon\(graphics\)', ''),
        (r'\n? *Icon\(graphics\), *', ' '),
        (r'\n? *experimentSetupOutput, *', ' '),
        (r',\n? *experimentSetupOutput\)', ')'),
        # Remove spaces on the outside of bold and underline tags.
        ('<b> ', ' <b>'),
        (' </b>', '</b> '),
        ('<u> ', ' <u>'),
        (' </u>', '</u> '),
        # Remove unless open/close tags.
        ('</b>( *<u>)<b>', r'\1'),
        ('</b>(</u> *)<b>', r'\1'),
        # Remove extra spacing.
        (r' +\n', r'\n'),
        (r'\n\n\n+', r'\n\n'),
        (' +<br>', '<br>'),
        ('<br><br>(<br>)+', '<br><br>'),
        # Remove empty lines above annotations.
        (r'\n+(\n *annotation\()', r'\1'),
        # Use shortcuts for Units and Quantities.
        (r' FCSys\.Quantities\.', ' Q.'),
        (r' FCSys\.Units\.', ' U.'),
        # Don't use 1 when unnecessary.
        (r'=1*U\.', '=U.'),
        # Use absolute referencing for Connectors.
        (r' Connectors\.', r' FCSys.Connectors.'),
        # Sometimes Dymola adds a useless import.
        (r'import FCSys;\n', ''),
        # No empty line before "end x;".
        (r'\n(\n +end )([^; ]+);', r'\1\2;'),
        # Two spaces after "Note:" and "TODO:".
        ('Note: +([^ ])', r'Note:  \1'),
        ('TODO: +([^ ])', r'TODO:  \1'),
        # Use lowercase "e" for engineering notation.
        ('([0-9]+)E(-?)([0-9]+)', r'\1e\2\3'),
        # Don't use a "+" for positive powers of 10 in engineering notation.
        ('([0-9]+)e+?([0-9]+)', r'\1e\2'),
        # Use some contractions in comments.
        ('(// .*)do not', r"\1don't"),
        ('(// .*)does not', r"\1doesn't"),
        ('(// .*)is not', r"\1isn't"),
        ('(// .*)will not', r"\1won't"),
        ('(// .*)cannot', r"\1can't"),
        ('(// .*)there is', r"\1there's"),
        (r"(// .*)it is ", r"\1it's "), # Trailing space prevents mixup with "isn't"
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
