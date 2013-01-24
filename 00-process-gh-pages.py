#!/usr/bin/env python
# Clean up the HTML help files by making custom replacements.
#
# The first argument is the directory, which defaults to the present directory.
#
# Created by Kevin Davies, 5/30/2012

import re, glob, sys, os, datetime

## Settings

# Replacement pairs
now = datetime.datetime.now()
rpls = [
    # Update the image links.
    ('img src=".*?([^/]+\.png)', r'img src="images/\1'),
    ('img src=".*?([^/]+\.svg)', r'img src="images/\1'),
    ('img src=".*?([^/]+\.gif)', r'img src="images/\1'),
    ('img src=".*?([^/]+\.pdf)', r'img src="images/\1'),
    ('img src=".*?([^/]+\.ico)', r'img src="images/\1'),
    # FCSys.html will be index.html.
    ('FCSys\.html', 'index.html'),
    # Change the title of the main page.
    ('<title>FCSys</title>', '<title>FCSys &mdash; Modelica library of fuel cell models</title>'),
    # Add keywords.
    ('<meta name="HTML-Generator" content="Dymola">', r"""<meta name="keywords" content="fuel cell library, FCSys, fuel cell, PEM, proton exchange membrane, polymer exchange membrance, PEMFC, Modelica, Dymola, electrochemistry">
<meta name="date" content="%d-%d-%d">""" % (now.year, now.month, now.day)),
    # Add the download link.
    ('(BaseClasses</a></li>\n *</ul>\n)( *</div>)', r"""\1
  <h3>Download</h3>
    <ul>
      <li>Latest: <a href="release/FCSys-2.0.zip" rel="nofollow">FCSys-2.0.zip</a> (**Please check back soon or contact kdavies4 at gmail.com.)</li>
    </ul>
\2"""),
    # Move the style sheet.
    ('"\.\./resources/documentation/ModelicaDoc\.css"', '"stylesheets/ModelicaDoc.css"'),
    # Add the Google Analytics script.
    ('(<link rel="shortcut icon" href=".*\.ico">\n)(</head>)', r"""\1<script type="text/javascript" src="javascripts/analytics.js"></script>
\2"""),
    # Remove the self-reference.
    ('Updates to this package may be available at.*\n *<a href="http://kdavies4\.github\.com/FCSys/">.*\n *Development is being carried out at', 'The development site is'),
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
