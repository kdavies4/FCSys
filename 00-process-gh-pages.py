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
    ('([A-z.]+).ico', r'images/\1.ico'),
    ('([A-z.]+).svg', r'images/\1.svg'),
    # FCSys.html will be index.html.
    ('"FCSys.html"', '"index.html"'),
    # Add the download link.
    ('BaseClasses</a></li>\n  </ul>\n  </div>', 'BaseClasses</a></li>\n  </ul><h3>Download</h3>\n  <ul>\n    <li>Latest: <a href="release/FCSys-2.0.zip"\n           rel="nofollow">FCSys-2.0.zip</a> (**Check back soon)</li>\n  </ul>\n\n  </div>'),    ('FCSys.html', r'index.html'),
    # Move the style sheet.
    ('"../resources/www/modelicaDoc.css"', '"stylesheets/modelicaDoc.css"'),
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
