#!/bin/bash
# Check that all of the references are used in the Modelica code.
#
# Kevin Davies, 7/18/12

for s in Ashcroft1976 Avogadro1.03 Bejan2006 BIPM2006 Bernardi1992 Cellier1991 Dassault2010 DuPont2004N DuPont2004NRE DuPont2005 Dymond2002 Entegris2012 Gurau1998 Hogan2006 Incropera2002 Kandlikar2009 Lin2006 Mark1999 McBride1996 McBride2002 Moran2004 NIST2010 Nitta2008 Present1958 Rao1997 Reichert2010 Salzman2004 Schetz1996 SGL2004 SGL2007 Shah2009 Springer1991 Spry2009 Svehla1995 Takenaka1990 Tissandier1998 Toray2010 Weber2004 Woo1995
do
    hits=`find \( -name "*.mo" -o -name "*.cdf" \) | xargs grep --files-with-matches --max-count=1 --exclude package.mo --exclude WorkInProgress.mo $s`
    if [ -z "$hits" ]; then
        echo \"$s\" was not found.
    fi
done
read -p "Done.  Press [Enter] to exit."
