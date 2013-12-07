#!/usr/bin/env python
# Clean up the HTML help files by making custom replacements.
#
# The first argument is the directory, which defaults to the help directory.
#
# Kevin Davies, 5/30/2012

import re, glob, sys, os

## Settings
stylesheet = '../Resources/Documentation/ModelicaDoc.css'
favicon = '../Resources/Documentation/favicon.ico'

# Replacement pairs
rpls = [
    # Use simple line endings.
    (r'\r\n', r'\n'),
    # Use lowercase tags (for consistency).
    ('<(/?)HTML>',  r'<\1html>'),
    ('<(/?)TITLE>',  r'<\1title>'),
    ('<(/?)HEAD>',  r'<\1head>'),
    ('<(/?)BODY>',  r'<\1body>'),
    ('<(/?)B>',  r'<\1b>'),
    ('<(/?)H([1-9]+)>', r'<\1h\2>'),
    ('<A +HREF *= *', '<a href='),
    ('<A +NAME *= *', '<a name='),
    ('</A>', '</a>'),
    ('<HR>', '<hr>'),
    ('<(/?)TABLE>', r'<\1table>'),
    ('<(/?)TD>', r'<\1td>'),
    ('<(/?)TR>', r'<\1tr>'),
    ('<(/?)P>', r'<\1p>'),
    ('<(/?)PRE>', r'<\1pre>'),
    ('<IMG +SRC *= *',  '<img src='),
    ('</IMG>',  '</img>'),
    ('<META +name', '<meta name'),
    ('<(/?)TH *>', r'<\1th>'), # Dymola adds a space.
    ('<TABLE +', '<table '),
    (' BORDER *= *', ' border='),
    (' CELLSPACING *= *', ' cellspacing='),
    (' CELLPADDING *= *', ' cellpadding='),
    (' ALT *= *', ' alt='),
    (' WIDTH *= *', ' width='),
    (r'( width=\d) >', r'\1>'),
    (' HEIGHT *= *', ' height='),
    (' ALIGN *= *TOP *', ' align=top'), # Dymola adds spaces.
    (' ALIGN *= *', ' align='), # Dymola adds spaces.
    # Remove strange line breaks.
    ('"\n>', '">'),
    # No space before newline
    (' *\n', '\n'),
    # Remove extra line spacing.
    (r' *<br> *(<br>)+ *', r'<br><br>'),
    # Remove empty groups.
    ('<p>\n*</p>', ''),    
    ('<pre>\n*</pre>', ''),
    (r'<pre>&quot;', r'&quot;'),
    (r'&quot;</pre>', r'&quot;'),
    # Remove the per-file license summary
    ("""<p><b>Licensed by the Georgia Tech Research Corporation.*
.*
.*
.*
.*
.*
.*
.*ModelicaLicense2.*""", ""),
    # Add the sidebar.
    ('<div class="body">>', """<div class="body">"""),
    ("""<body>

""", """<body>
<div class="sidebar">
  <div class="sidebarwrapper">
  <a href="FCSys.html"><p class="sidebar-title">FCSys</p></a>
  <p class="logo"><a href="FCSys.html">
    <img src="../Resources/Documentation/icon.gif" class="logo" alt="Logo" width=150>
  </a></p>

  <h3>Table of Contents</h3>
    <ul>
    <li><a href="FCSys_UsersGuide.html">User's Guide</a></li>
    <li><a href="FCSys_Blocks.html">Blocks</a></li>
    <li><a href="FCSys_Conditions.html">Conditions</a></li>
    <li><a href="FCSys_Assemblies.html">Assemblies</a></li>
    <li><a href="FCSys_Regions.html">Regions</a></li>
    <li><a href="FCSys_Subregions.html">Subregions</a></li>
    <li><a href="FCSys_Phases.html">Phases</a></li>
    <li><a href="FCSys_Species.html">Species</a></li>
    <li><a href="FCSys_Chemistry.html">Chemistry</a></li>
    <li><a href="FCSys_Connectors.html">Connectors</a></li>
    <li><a href="FCSys_Characteristics.html">Characteristics</a></li>
    <li><a href="FCSys_Units.html">Units</a></li>
    <li><a href="FCSys_Quantities.html">Quantities</a></li>
    <li><a href="FCSys_Utilities.html">Utilities</a></li>
    <li><a href="FCSys_Icons.html">Icons</a></li>
    </ul>
  </div>
</div>

<div class="document">
  <div class="documentwrapper">
    <div class="bodywrapper">
      <div class="body">"""),
    # Adjust the footer.
    ("""<hr>
 *<address><a href.*>Automatically generated</a> *([^.]+)\.
 *</address>""", """      </div>
    </div>
  </div>
</div>
<div class="footer">
  &copy; Copyright 2012&ndash;2013, Kevin Davies, Georgia Tech Research Corporation. Last updated %s.
</div>
""" % r'\1'),
    # Remove Dymola's tags for the Microsoft Office template.
    ("<!--.+'mso-element:field-begin'.+\n.+\n.+'mso-element:field-end'.+-->", ''),
    # Remove the (useless?) textblock tags.
    ('<textblock type="annotcomp" expanded="false">', ''),
    ('<textblock type="annotconnect" expanded="false">', ''),
    ('<textblock type="model" expanded="false" path=".*">', ''),
    ('</textblock>', ''),
    # Make the local links local again.
    ('<a href="FCSys.*\.html(#\w*\d+)">', r'<a href="\1">'),
    # Use relative links again.
    ('/media/Storage/Documents/Modelica/FCSys/FCSys[^/]*/', '../'), # Linux
    ('file:///D:/Documents/Modelica/FCSys/FCSys[^/]*/', '../'), # Windows
    # Remove the return from footnote links (they don't work).
    ('<a href="#ref\d+">&#8629;</a>', ''),
    # Use online links to the Modelica Standard Library.
    ('"file:////opt/dymola/Modelica/Library/Modelica.*/help/(.*)"', r'"http://build.openmodelica.org/Documentation/\1"'),
    ('"(http://build\.openmodelica\.org/Documentation/[^"]+)_([^"]+)_([^"]+)_([^"]+)_([^"]+)_([^"]+)_([^"]+)_([^"]+)"', r'"\1.\2.\3.\4.\5.\6.\7.\8"'),
    ('"(http://build\.openmodelica\.org/Documentation/[^"]+)_([^"]+)_([^"]+)_([^"]+)"', r'"\1.\2.\3.\4"'),
    ('"(http://build\.openmodelica\.org/Documentation/[^"]+)_([^"]+)"', r'"\1.\2"'),
    ('"(http://build\.openmodelica\.org/Documentation/[^"]+)_([^"]+)"', r'"\1.\2"'),
    # Use the better versions of some images.
    ('/cell\.png"', '/cell.svg" width=450'),
    ('/FCSys\.Assemblies\.Cells\.CellD\.png"', '/FCSys.Assemblies.Cells.CellD.png" width=600'),
    # Remove some classes from the main page.
    ('<tr><td><img .+>.+<a .+>Figures</a>\n?</td><td>.+</td></tr>\n?', ''),
    ('<tr><td><img .+>.+<a .+>WorkInProgress</a>\n?</td><td>.+</td></tr>\n?', ''),
    # Remove nested quotes from meta description.
    ('(<meta name="description" content=")&quot;(.*)&quot;(">)', r'\1\2\3'),
    # Don't escape the quotation marks.
    (r'\\"', '"'),
    # Change the style.
    ("""<style type="text/css"> *
 *\*.*
 *pre.*
 *h4.*
 *h3.*
 *h2.*
 *address.*
 *td.*
 *th.*
 *table.*
 *</style>""", '<link rel="stylesheet" type="text/css" charset="utf-8" media="all" href="%s">' % stylesheet + '\n<link rel="shortcut icon" href="%s">' % favicon),
    # Try to replace the internal Modelica references with the proper HTML
    # page.
    (r'\\?"modelica://([^"\\]+)\.([^"\\]+)\.([^"\\]+)\.([^"\\]+)\.([^"\\]+)\.([^"\\]+)\.([^"\\]+)\\?"', r'"\1_\2_\3_\4_\5_\6.html#\1.\2.\3.\4.\5.\6.\7"'),
    (r'\\?"modelica://([^"\\]+)\.([^"\\]+)\.([^"\\]+)\.([^"\\]+)\.([^"\\]+)\.([^"\\]+)\\?"', r'"\1_\2_\3_\4_\5.html#\1.\2.\3.\4.\5.\6"'),
    (r'\\?"modelica://([^"\\]+)\.([^"\\]+)\.([^"\\]+)\.([^"\\]+)\.([^"\\]+)\\?"', r'"\1_\2_\3_\4.html#\1.\2.\3.\4.\5"'),
    (r'\\?"modelica://([^"\\]+)\.([^"\\]+)\.([^"\\]+)\.([^"\\]+)\\?"', r'"\1_\2_\3.html#\1.\2.\3.\4"'),
    (r'\\?"modelica://([^"\\]+)\.([^"\\]+)\.([^"\\]+)\\?"', r'"\1_\2.html#\1.\2.\3"'),
    (r'\\?"modelica://([^"\\]+)\.([^"\\]+)\\?"', r'"\1.html#\1.\2"'),
    (r'\\?"modelica://([^"\\]+)\\?"', r'"\1.html#\1"'),
    ('href="Modelica', 'href="http://build.openmodelica.org/Documentation/Modelica'),
    # Use the ellipsis code.
    ('\.\.\.', '&hellip;'),
    # Eliminate back-to-back paragraph starts.
    ('<p> *\n? *(<p>)+', '<p>'),
    # Remove paragraph marks around headings.
    ('<p>(<h\d+>.+</h\d+>)</p>', r'\1'),
    ('<p>(<h3>Package Content</h3>)<p>', r'\1'),
    # Newline before and after headings:
    ('(<h\d+>.+</h\d+>)', r'\n\1\n'),
    # but no double blank lines:
    ('(\n *){3,}', '\n\n'),
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
