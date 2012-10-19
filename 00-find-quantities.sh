#!/bin/bash
# Check if the quantities are used in the Modelica code.
#
# Kevin Davies, 7/18/12

for s in Acceleration Angle Angle2 Area Capacitance ConductanceElectrical Current CurrentAreic CurrentRate Energy EnergyMassic Force Frequency Inductance Length LengthSpecific MagneticFlux MagneticFluxAreic MagneticFluxReciprocal Mass MassSpecific MomentumAngular Number NumberAbsolute NumberRate ParticleNumber ParticleNumberReciprocal ParticleNumberVolumic Permeability Permittivity PermittivityReciprocal Potential PotentialAbsolute PotentialPerWavenumber PotentialRate Power PowerArea PowerAreic PowerAreicPerPotential4 PowerRadiant Pressure PressureAbsolute PressureRate Resistance ResistanceElectrical Resistivity Time Velocity VelocityMassSpecific Volume VolumeRate VolumeSpecific VolumeSpecificRate Wavenumber
do
    grep --color Q.$s --exclude WorkInProgress.mo --exclude Quantities.mo *.mo
    read dummy
    done
