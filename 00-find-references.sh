#!/bin/bash
# Check if the references are used in the Modelica code.
#
# Kevin Davies, 7/18/12

for s in Ashcroft1976 Avogadro1.03 Bejan2006 Bernardi1992 Cellier1991 DuPont2004N DuPont2005 DuPont2004NRE Dymond2002 Entegris2012 Gurau1998 Hogan2006 Incropera2002 BIPM2006 Lin2006 Kandlikar2009 Mark1999  McBride1996 McBride2002 NIST2010 Nitta2008 Present1958 Reichert2010 Salzman2004 Schetz1996 SGL2007 SGL2004 Springer1991 Spry2009 Svehla1995 Takenaka1990 Tissandier1998 Toray2010 Weber2004 Woo1995
do
    #grep --color $s --exclude package.mo --exclude WorkInProgress.mo *.mo
    find \( -name "*.mo" -o -name "*.cdf" \) | xargs grep -niP --color --exclude package.mo --exclude WorkInProgress.mo $s
    read dummy
    done
