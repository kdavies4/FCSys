#!/bin/bash
# Check for quantities and references that aren't used in the Modelica code.
#
# Kevin Davies, 7/18/12

# Quantities
for s in Acceleration Amount AmountReciprocal AmountVolumic Angle Angle2 Area Capacitance CompressibilityDynamic ConductanceElectrical Current CurrentAreic CurrentRate Energy FluidityDynamic Force Frequency Inductance Length LengthSpecific MagneticFlux MagneticFluxAreic MagneticFluxReciprocal Mass MassSpecific MomentumAngular Number NumberAbsolute Permeability Permittivity PermittivityReciprocal Potential PotentialAbsolute PotentialPerWavenumber PotentialRate Power PowerArea PowerAreic PowerAreicPerPotential4 PowerRadiant Pressure PressureAbsolute PressureRate PressureReciprocal ResistanceElectrical Resistivity Time Velocity Velocity2 Volume VolumeRate VolumeSpecific VolumeSpecificAbsolute VolumeSpecificRate Wavenumber
do
    hits=`grep --color "Q.$s " --files-with-matches --max-count=1 --exclude WorkInProgress.mo --exclude Quantities.mo *.mo`
    if [ -z "$hits" ]; then
        echo Quantity \"$s\" was not found.
    fi
done

# References
for s in Avogadro1.03 Bejan2006 BIPM2006 Bernardi1992 Dassault2010 DuPont2004N DuPont2004NRE DuPont2005 Dymond2002 Entegris2012 Gurau1998 Hogan2006 Incropera2002 Kandlikar2009 Lin2006 Mark1999 McBride1996 McBride2002 Moran2004 NIST2010 Nitta2008 Present1958 Rao1997 Reichert2010 Salzman2004 Schetz1996 SGL2004 SGL2007 Shah2009 Springer1991 Spry2009 Svehla1995 Takenaka1990 Tissandier1998 Toray2010 Weber2004 Woo1995
do
    hits=`find \( -iname "*.mo" -o -iname "*.cdf" \) | xargs grep --files-with-matches --max-count=1 --exclude package.mo --exclude WorkInProgress.mo $s`
    if [ -z "$hits" ]; then
        echo Reference \"$s\" was not found.
    fi
done
read -p "Done.  Press [Enter] to exit."
