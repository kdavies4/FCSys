#!/bin/bash
# Check that all of the quantities are used in the Modelica code.
#
# Kevin Davies, 7/18/12

for s in Acceleration Amount AmountReciprocal AmountVolumic AmountVolumicRate Angle Angle2 Area Capacitance CompressibilityDynamic ConductanceElectrical Current CurrentAreic CurrentRate Energy FluidityDynamic Force Frequency Inductance Length LengthSpecific MagneticFlux MagneticFluxAreic MagneticFluxReciprocal Mass MassSpecific MomentumAngular Number NumberAbsolute Permeability Permittivity PermittivityReciprocal Potential PotentialAbsolute PotentialPerWavenumber PotentialRate Power PowerArea PowerAreic PowerAreicPerPotential4 PowerRadiant Pressure PressureAbsolute PressureRate ResistanceElectrical Resistivity Time Velocity Velocity2 Volume VolumeRate VolumeSpecific VolumeSpecificAbsolute Wavenumber
do
    hits=`grep --color "Q.$s " --files-with-matches --max-count=1 --exclude WorkInProgress.mo --exclude Quantities.mo *.mo`
    if [ -z "$hits" ]; then
        echo \"$s\" was not found.
    fi
done
read -p "Done.  Press [Enter] to exit."
