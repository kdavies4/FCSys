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
    ('img src=".*?([^/]+\.png)', r'img src="images/\1'),
    ('img src=".*?([^/]+\.svg)', r'img src="images/\1'),
    ('img src=".*?([^/]+\.gif)', r'img src="images/\1'),
    ('img src=".*?([^/]+\.pdf)', r'img src="images/\1'),
    ('img src=".*?([^/]+\.ico)', r'img src="images/\1'),
    #('<img src="[A-Za-z0-9._]+/[A-Za-z0-9._]+/[A-Za-z0-9._]+/[A-Za-z0-9._]+/[A-Za-z0-9._]+/', '<img src="'),
    #('<img src="[A-Za-z0-9._]+/[A-Za-z0-9._]+/', '<img src="'), # Repeated for subdirectories
    #('<img src="[A-Za-z0-9._]+/', '<img src="'),
    #('([A-Za-z0-9._]+)\.png', r'images/\1.png'),
    #('([A-Za-z0-9._]+)\.svg', r'images/\1.svg'),
    #('([A-Za-z0-9._]+)\.gif', r'images/\1.gif'),


#src="../resources/documentation/Connectors/hierarchy.png">
#src="Connectors/images/hierarchy.png">


#    ('\.\./resources/documentation/(favicon\.ico)', r'images/\1'),
#    ('\.\./resources/documentation/(.+\.pdf)', r'images/\1'),
    # FCSys.html will be index.html.
    ('FCSys\.html', 'index.html'),
    # Change the title of the main page.
    ('<title>FCSys</title>', '<title>FCSys &mdash; Modelica library of fuel cell models</title>'),
    # Add keywords.
    ('(<meta name="description" content=".*">)', r"""\1
<meta name="keywords" content="Modelica, fuel cell, PEMFC, electrochemistry, model">"""),
    # Add the download link.
    ('(BaseClasses</a></li>\n *</ul>\n)( *</div>)', r"""\1
  <h3>Download</h3>
    <ul>
      <li>Latest: <a href="release/FCSys-2.0.zip" rel="nofollow">FCSys-2.0.zip</a> (**Check back soon.)</li>
    </ul>
\2"""),
    # Move the style sheet.
    ('"\.\./resources/documentation/ModelicaDoc\.css"', '"stylesheets/ModelicaDoc.css"'),
    # Add the Google Analytics script.
    ('(<link rel="shortcut icon" href=".*\.ico">\n)(</head>)', r"""\1<script type="text/javascript" src="javascripts/analytics.js"></script>
\2"""),
    # Remove the self-reference.
    ('Updates to this package may be available at the *\n *<a href="http://kdavies4\.github\.com/FCSys/">main project site</a>\. *\n *Development is being carried out at', 'The development site is'),
    # **Temporary:  Add note to the sidebar.
    ("""<h3>Table of Contents</h3>
    <ul>
    <li><a href="FCSys_UsersGuide\.html">User's Guide</a></li>
    <li><a href="FCSys_BCs\.html">BCs</a></li>
    <li><a href="FCSys_Sensors\.html">Sensors</a></li>
    <li><a href="FCSys_Assemblies\.html">Assemblies</a></li>
    <li><a href="FCSys_Regions\.html">Regions</a></li>
    <li><a href="FCSys_Subregions\.html">Subregions</a></li>
    <li><a href="FCSys_Connectors\.html">Connectors</a></li>
    <li><a href="FCSys_Characteristics\.html">Characteristics</a></li>
    <li><a href="FCSys_Units\.html">Units</a></li>
    <li><a href="FCSys_Quantities\.html">Quantities</a></li>
    <li><a href="FCSys_BaseClasses\.html">BaseClasses</a></li>
    </ul>""",
    """<h3>Table of Contents</h3>
    <ul>
    <li><a href="FCSys_UsersGuide.html">User's Guide</a></li>
    <li><a href="FCSys_BCs.html">BCs</a>**</li>
    <li><a href="FCSys_Sensors.html">Sensors</a>**</li>
    <li><a href="FCSys_Assemblies.html">Assemblies</a>**</li>
    <li><a href="FCSys_Regions.html">Regions</a>**</li>
    <li><a href="FCSys_Subregions.html">Subregions</a></li>
    <li><a href="FCSys_Connectors.html">Connectors</a></li>
    <li><a href="FCSys_Characteristics.html">Characteristics</a></li>
    <li><a href="FCSys_Units.html">Units</a></li>
    <li><a href="FCSys_Quantities.html">Quantities</a></li>
    <li><a href="FCSys_BaseClasses.html">BaseClasses</a></li>
    <li>**Please check back soon or contact kdavies4 at gmail.com.</li>
    </ul>"""),
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
