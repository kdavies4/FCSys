#!/bin/bash
# Check that all of the quantities are used in the Modelica code.
#
# Kevin Davies, 7/18/12

for s in Number NumberAbsolute AmountReciprocal NumberRate Angle MagneticFluxReciprocal Wavenumber Frequency Angle2 Length Permeability Force Fluidity Resistivity Acceleration LengthSpecific Velocity Acceleration Area MagneticFlux MomentumAngular PowerRadiant PotentialRate Potential PotentialAbsolute PotentialRate ResistanceElectrical Inductance Energy Power Velocity2 Volume PotentialPerWavenumber PermittivityReciprocal VolumeSpecific VolumeSpecificAbsolute VolumeSpecific VolumeRate PowerArea Mass PowerAreicPerPotential4 MagneticFluxAreic PressureRate Pressure PressureAbsolute PressureRate MassSpecific PowerAreic Amount CurrentAreic AmountVolumicRate CurrentRate AmountVolumic Current CurrentRate ConductanceElectrical Capacitance Permittivity NumberRate Time Fusivity
do
    hits=`grep --color "Q.$s " --files-with-matches --max-count=1 --exclude WorkInProgress.mo --exclude Quantities.mo *.mo`
    if [ -z "$hits" ]; then
        echo \"$s\" was not found.
    fi
done
echo Done.  Press enter to exit.
read answer
