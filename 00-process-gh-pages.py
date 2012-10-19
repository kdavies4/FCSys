#!/usr/bin/env python
# Clean up the HTML help files by making custom replacements.
#
# The first argument is the directory, which defaults to ./.
#
# Created by Kevin Davies, 5/30/2012

import re, glob, sys, os

## Settings

# Replacement pairs
rpls = [
    # Update the image links.
    ('<img src="/[A-z.]+/', '<img src="'),
    ('<img src="[A-z.]+/', '<img src="'), # Repeated due to subdirectories
    ('<img src="[A-z.]+/', '<img src="'),
    ('<img src="[A-z.]+/', '<img src="'),
    ('<img src="[A-z.]+/', '<img src="'),
    ('<img src="[A-z.]+/', '<img src="'),
    ('<img src="[A-z.]+/', '<img src="'),
    ('<img src="[A-z.]+/', '<img src="'),
    ('([A-z.]+).png', r'images/\1.png'),
    # Use the SVG version of some images.
    ('"images/cell\.png"', '"images/cell.svg" width=450'),
    ('"images/FCSys\.Subassemblies\.Cells\.CellD\.png"', '"images/FCSys.Subassemblies.Cells.CellD.svg" width=600'),
    ('"images/FCSys\.Subassemblies\.Cells\.Examples\.CellProfileD\.png"', '"images/FCSys.Subassemblies.Cells.Examples.CellProfileD.svg" width=250'),
    # FCSys.html will be index.html.
    ('FCSys.html', r'index.html'),
    # Move the style sheet and icon.
    ('"../resources/www/modelicaDoc.css"', '"stylesheets/modelicaDoc.css"'),
    ('"../resources/images/favicon.ico"', '"images/favicon.ico"'),
    ]

# Directory specification
if (len(sys.argv) > 1):
    directory = sys.argv[1]
else:
    directory = '.'

## Procedure
# Compile the regular expressions.
for i, rpl in enumerate(rpls):
    rpls[i] = (re.compile(rpl[0]), rpl[1])

# Replace strings.
for fname in glob.glob(os.path.join(directory, '*.html')):
    # Read the source file.
    print "Processing " + fname + "..."
    src = open(fname, 'r')
    text = src.read()
    src.close()

    # Make the replacements.
    for rpl in rpls:
        text = rpl[0].sub(rpl[1], text)

    # Re-write the file.
    src = open(fname, 'w')
    src.write(text)
    src.close()
