#!/bin/bash
# Check for quantities and references that aren't used in the Modelica code.
#
# Kevin Davies, 7/18/12

# Quantities
for s in Acceleration Amount AmountReciprocal Angle Angle2 Area AreaSpecific Capacitance ConductanceElectrical Current CurrentAreic CurrentAreicAbsolute CurrentRate Density DensityRate Diffusivity Energy Fluidity Force Frequency Inductance Length LengthSpecific MagneticFlux MagneticFluxAreic MagneticFluxReciprocal Mass MassRate MassSpecific Mobility MomentumRotational Number NumberAbsolute Permeability Permittivity PermittivityReciprocal Potential PotentialAbsolute PotentialLineic PotentialPerWavenumber PotentialRate Power PowerArea PowerAreic PowerAreicPerPotential4 PowerRadiant Pressure PressureAbsolute PressureRate PressureReciprocal ResistanceElectrical Resistivity ResistivityMaterial Time TimeAbsolute TimeLineic Velocity Velocity2 Volume VolumeRate VolumeSpecific
do
    hits=`grep --color "Q.$s " --files-with-matches --max-count=1 --exclude WorkInProgress.mo --exclude Quantities.mo *.mo`
    if [ -z "$hits" ]; then
        echo Quantity \"$s\" was not found.
    fi
done

# References
for s in Aronsson2009 Avogadro1.03 Bejan2006 Broman2008 Brown2011 BIPM2006 Bernardi1992 Dassault2010 DuPont2004N DuPont2004NRE DuPont2005 Dymond2002 Entegris2012 Fritzson2004 Greiner1995 Gurau1998 Hess2008 Incropera2002 Kandlikar2009 Larminie2003 Lin2006 Mark1999 Mattsson2008 Mattsson1993B McBride1996 McBride2002 Modelica3.2 Moran2004 NIST2010 Nitta2008 Present1958 Rao1997 Rapaport2004 Reichert2010 Salzman2004 Schetz1996 SGL2004 SGL2007 Shah2009 Springer1991 Spry2009 Svehla1995 Takenaka1990 Tissandier1998 Toray2010 Weber2004 Woo1995
do
    hits=`find \( -iname "*.mo" -o -iname "*.cdf" \) | xargs grep --files-with-matches --max-count=1 --exclude package.mo --exclude WorkInProgress.mo $s`
    if [ -z "$hits" ]; then
        echo Reference \"$s\" was not found.
    fi
done
read -p "Done.  Press [Enter] to exit."
