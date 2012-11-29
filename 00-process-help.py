#!/usr/bin/env python
# Clean up the HTML help files by making custom replacements.
#
# The first argument is the directory, which defaults to ./.
#
# Created by Kevin Davies, 5/30/2012

import re, glob, sys, os

## Settings
stylesheet = '../resources/documentation/ModelicaDoc.css'
favicon = '../resources/documentation/favicon.ico'

# Replacement pairs
rpls = [
    # Use lowercase tags (for consistency).
    ('<(/?)HTML>',  r'<\1html>'),
    ('<(/?)TITLE>',  r'<\1title>'),
    ('<(/?)HEAD>',  r'<\1head>'),
    ('<(/?)BODY>',  r'<\1body>'),
    ('<(/?)B>',  r'<\1b>'),
    ('<(/?)H([1-9]+)>', r'<\1h\2>'),
    ('<A HREF=', '<a href='),
    ('<A NAME=', '<a name='),
    ('</A>', '</a>'),
    ('<HR>', '<hr>'),
    ('<(/?)TABLE>', r'<\1table>'),
    ('<(/?)TD>', r'<\1td>'),
    ('<(/?)TR>', r'<\1tr>'),
    ('<(/?)P>', r'<\1p>'),
    ('<(/?)PRE>', r'<\1pre>'),
    ('<IMG SRC=',  '<img src='),
    ('</IMG>',  '</img>'),
    ('<META name', '<meta name'),
    ('<(/?)TH *>', r'<\1th>'), # Dymola adds a space.
    ('<TABLE ', '<table '),
    (' BORDER *= *', ' border='),
    (' CELLSPACING *= *', ' cellspacing='),
    (' CELLPADDING *= *', ' cellpadding='),
    (' ALT *= *', ' alt='),
    (' WIDTH *= *', ' width='),
    (' HEIGHT *= *', ' height='),
    (' ALIGN *= *TOP *', ' align=top'), # Dymola adds spaces.
    (' ALIGN *= *', ' align='), # Dymola adds spaces.
    # Remove strange line breaks.
    ('"\n>', '">'),
    # Remove empty groups.
    (r'&quot;</pre>', r'&quot;'),
    (r'<pre>&quot;', r'&quot;'),
    ('<pre>\n*</pre>', r''),
    # Remove extra line spacing.
    (r' *<br> +(<br>)* *', r'<br><br>'),
    (r'<br><br>(<br>)+', r'<br><br>'),
    # Make the local links local again.
    ('<a href="FCSys\.html#Fig([1-3])">', r'<a href="#Fig\1">'),
    # Add the sidebar.
    ('<body><p>', """<body>
<div class="sidebar">
  <div class="sidebarwrapper">
  <a href="FCSys.html"><p class="sidebar-title">FCSys</p></a>
  <p class="logo"><a href="FCSys.html">
    <img src="../resources/documentation/icon.gif" class="logo" alt="Logo" width=150>
  </a></p>

  <h3>Table of Contents</h3>
    <ul>
    <li><a href="FCSys_UsersGuide.html">User's Guide</a></li>
    <li><a href="FCSys_BCs.html">BCs</a></li>
    <li><a href="FCSys_Sensors.html">Sensors</a></li>
    <li><a href="FCSys_Assemblies.html">Assemblies</a></li>
    <li><a href="FCSys_Regions.html">Regions</a></li>
    <li><a href="FCSys_Subregions.html">Subregions</a></li>
    <li><a href="FCSys_Connectors.html">Connectors</a></li>
    <li><a href="FCSys_Characteristics.html">Characteristics</a></li>
    <li><a href="FCSys_Units.html">Units</a></li>
    <li><a href="FCSys_Quantities.html">Quantities</a></li>
    <li><a href="FCSys_BaseClasses.html">BaseClasses</a></li>
    </ul>
  </div>
</div>

<div class="document">
  <div class="documentwrapper">
    <div class="bodywrapper">
      <div class="body">"""),
    # Adjust the footer.
    ('<hr>\n<address>HTML-documentation generated by [^>]+>[^>]+> ([^.]+).\n</address>', """      </div>
    </div>
  </div>
</div>
<div class="footer">
  &copy; Copyright 2012, Kevin Davies. Last updated %s.
</div>
""" % r'\1'),
    # Remove horizontal lines.
    #(r'<hr>', ''),
    # Remove the internal links.
    #(r'<a href="FCSys[^>]+\n*[^>]*>(.*)</a>', r'<code>\1</code>'),
    #(r'<a href="file://[^>]+>(.*)</A>', r'<code>\1</code>'),
    # Remove Dymola tags for the Microsoft office template.
    ("<!--\[if supportFields\]><span style='mso-element:field-begin'></span>\n<span style='mso-spacerun:yes'></span>XE [A-Za-z.]+<!\[endif\]-->\n<!--\[if supportFields\]><span style='mso-element:field-end'></span><!\[endif\]-->", ''),
    # Make the local links local again.
    ('<a href="FCSys\.html#Fig([1-3])">', r'<a href="#Fig\1">'),
    # Use the better versions of some images.
    ('/cell\.png"', '/cell.svg" width=450'),
    ('/FCSys\.Assemblies\.Cells\.CellD\.png"', '/FCSys.Assemblies.Cells.CellD.png" width=600'),
    # Remove some classes from the main page.
    ('<tr><td><img src="FCSys.BlocksS.png" alt="FCSys.Blocks" width=20  height=20 align=top>&nbsp;<a href="FCSys_Blocks.html#FCSys.Blocks">Blocks</a>\n</td><td>[^<]+</td></tr>\n', ''),
    ('<tr><td><img src="FCSys.BlocksS.png" alt="FCSys.Figures" width=20  height=20 align=top>&nbsp;<a href="FCSys_Figures.html#FCSys.Figures">Figures</a>\n</td><td>[^<]+</td></tr>\n', ''),
    ('<tr><td><img src="FCSys.SystemsS.png" alt="FCSys.Systems" width=20  height=20 align=top>&nbsp;<a href="FCSys_Systems.html#FCSys.Systems">Systems</a>\n</td><td>[^<]+</td></tr>\n', ''),
    ('<tr><td><img src="FCSys.SystemsS.png" alt="FCSys.WorkInProgress" width=20  height=20 align=top>&nbsp;<a href="FCSys_WorkInProgress.html#FCSys.WorkInProgress">WorkInProgress</a>\n</td><td>[^<]+</td></tr>\n', ''),
    # Use relative links again.
    ('/media/Storage/Documents/Dymola/FCSys/', '../'),
    # Remove nested quotes from meta description.
    ('(<meta name="description" content=")&quot;(.*)&quot;(">)', r'\1\2\3'),
    # Change the style.
    ('<style type="text/css"> *\n\*       \{ font-size: 10pt; font-family: Arial,sans-serif; \} *\npre     \{ font-size:  9pt; font-family: Courier,monospace;\} *\nh4      \{ font-size: 10pt; font-weight: bold; color: green; \} *\nh3      \{ font-size: 11pt; font-weight: bold; color: green; \} *\nh2      \{ font-size: 13pt; font-weight: bold; color: green; \} *\naddress \{                  font-weight: normal} *\ntd      \{ solid \#000; vertical-align:top; \} *\nth      \{ solid \#000; vertical-align:top; font-weight: bold; \} *\ntable   \{ solid \#000; border-collapse: collapse;\}\n</style>', '<link rel="stylesheet" type="text/css" charset="utf-8" media="all" href="%s">' % stylesheet + '\n<link rel="shortcut icon" href="%s">' % favicon),
    # Remove the custom style for the Modelica license.
    #('<style type=\"text/css\">\n\*       \{ font-size: 10pt; font-family: Arial,sans-serif; \}\ncode    \{ font-size:  9pt; font-family: Courier,monospace;\}\nh6      \{ font-size: 10pt; font-weight: bold; color: green; \}\nh5      { font-size: 11pt; font-weight: bold; color: green; \}\nh4      \{ font-size: 13pt; font-weight: bold; color: green; \}\naddress \{                  font-weight: normal\}\ntd      \{ solid #000; vertical-align:top; \}\nth      \{ solid #000; vertical-align:top; font-weight: bold; \}\ntable   \{ solid #000; border-collapse: collapse;\}\n</style>', ''),
    # Try to replace the internal Modelica references with the proper HTML
    # page.
    ('"modelica://([^.]+)"', r'"\1.html"'),
    ('"modelica://([^.]+)\.([^.]+)"', r'"\1.html#\1.\2"'),
    ('"modelica://([^.]+)\.([^.]+)\.([^.]+)"', r'"\1_\2.html#\1.\2.\3"'),
    ('"modelica://([^.]+)\.([^.]+)\.([^.]+)\.([^.]+)"', r'"\1_\2_\3.html#\1.\2.\3.\4"'),
    ('"modelica://([^.]+)\.([^.]+)\.([^.]+)\.([^.]+)\.([^.]+)"', r'"\1_\2_\3_\4.html#\1.\2.\3.\4.\5"'),
    ('"modelica://(v)\.([^.]+)\.([^.]+)\.([^.]+)\.([^.]+)\.([^.]+)"', r'"\1_\2_\3_\4_\5.html#\1.\2.\3.\4.\5.\6"'),
    ('"modelica://([^.]+)\.([^.]+)\.([^.]+)\.([^.]+)\.([^.]+)\.([^.]+)\.([^.]+)"', r'"\1_\2_\3_\4_\5_\6.html#\1.\2.\3.\4.\5.\6.\7"'),
    ]

# Directory specification
if (len(sys.argv) > 1):
    directory = sys.argv[1]
else:
    directory = './help'

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
