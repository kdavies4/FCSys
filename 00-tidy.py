#!/usr/bin/env python
# Make custom replacements in files.
#
# The first argument is the directory.
#
# Kevin Davies, 5/30/2012

import re
import glob
import sys
import os

## Settings

# Replacement pairs
rpls = [# Remove tabs.
        ('\t', ' '),
        # No space before newline
        (' +[\n\r]', '\n'),
        # Remove empty annotation tags.
        (r'\n? *Diagram\(graphics\), *', ' '),
        (r',\n? *Diagram\(graphics\)', ''),
        (r'annotation\(Diagram\(graphics\)\);', ''),
        (r'\n? *Icon\(graphics\), *', ' '),
        (r',\n? *Icon\(graphics\)', ''),
        (r'annotation\(Icon\(graphics\)\);', ''),
        (r'\n? *experimentSetupOutput, *', ' '),
        (r',\n? *experimentSetupOutput\)', ')'),
        # Remove spaces on the outside of bold and underline tags.
        ('<b> ', ' <b>'),
        (' </b>', '</b> '),
        ('<u> ', ' <u>'),
        (' </u>', '</u> '),
        # Insert two newlines between paragraphs.
        (r' *(</p>) *\n? *( *<p>)', r'\1\n\n\2'),
        # No newline after start of paragraph or before end
        (r'(<p>) *\n? *', r'\1'),
        (r' *\n? *(</p>)', r'\1'),
        # Remove unless open/close tags.
        ('</b>( *<u>)<b>', r'\1'),
        ('</b>(</u> *)<b>', r'\1'),
        # Start single-line comments with a space.
        (r'([ \n]//)([^ ])', r'\1 \2'),
        # Remove extra spacing.
        (r' +\n', r'\n'),
        (r'\n\n\n+', r'\n\n'),
        (' +<br>', '<br>'),
        ('<br><br>(<br>)+', '<br><br>'),
        # Remove empty lines above annotations.
        (r'\n+(\n *annotation\()', r'\1'),
        # Use shortcuts for Units and Quantities.
        (r' FCSys\.Quantities\.', ' Q.'), # Leading spaces distinguish these from hyperlinks.
        (r' FCSys\.Units\.([^*])', r' U.\1'),
        # Don't use factor of 1 if unnecessary.
        (r'=1*U\.', '=U.'),
        # Use relative references where possible.
        (r'  FCSys.Conditions\.', r'  Conditions.'), # Two spaces to prevent change of import statements
        (r'  FCSys.Assemblies\.', r'  Assemblies.'),
        #(r'  FCSys.Regions\.', r'  Regions.'),
        #(r'  FCSys.Subregions\.', r'  Subregions.'),
        (r'  FCSys.Connectors\.', r'  Connectors.'),
        (r'  FCSys.Characteristics\.', r'  Characteristics.'),
        # Remove useless import statements.
        (r'import FCSys;\n', ''),
        # One empty line before "end x;"
        (r'\n*(\n *end) +(.*;)', r'\n\1 \2'),
        # unless end if, for, loop, when, or while.
        (r'\n(\n *end if;)', r'\1'),
        (r'\n(\n *end for;)', r'\1'),
        (r'\n(\n *end loop;)', r'\1'),
        (r'\n(\n *end when;)', r'\1'),
        (r'\n(\n *end while;)', r'\1'),
        # No empty line before annotation
        (r'\n+(\n *annotation)', r'\1'),
        # One empty line before algorithm, equation, initial equation,
        # protected, and public but not after.
        (r'\n*(\n *algorithm\n)\n*', r'\n\1'),
        (r'\n*(\n *equation\n)\n*', r'\n\1'),
        (r'\n*(\n *initial equation\n)\n*', r'\n\1'),
        (r'\n*(\n *protected\n)\n*', r'\n\1'),
        (r'\n*(\n *public\n)\n*', r'\n\1'),
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
        # Indent the HTML for paragraphs by 2 spaces past the annotation. **doesn't work yet
        #(r'( *)(annotation *\(.*(?:\n.*).*?\n *Documentation\(info="<html>.*(?:\n.*)*?\n) *(<p>.*(?:\n.*).*?\n.*</p>)', r'\1\2\1  \3'),
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
    rpls[i] = (re.compile(rpl[0], re.MULTILINE), rpl[1])

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
