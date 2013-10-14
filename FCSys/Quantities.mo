within FCSys;
package Quantities "Quantities to represent physical properties"
  package Examples "Examples"
    extends Modelica.Icons.ExamplesPackage;

    model Display "Demonstrate the display units for the quantities"
      extends Modelica.Icons.Example;

      ExampleModel doubleClickMe annotation (Placement(transformation(
            extent={{-20,-10},{20,10}},
            rotation=0,
            origin={0,0})));
      annotation (Diagram(graphics));

    end Display;

    model ExampleModel "Model that uses all of the quantities"
      extends FCSys.Icons.Blocks.Continuous;

      // Generated from FCSys/Resources/quantities.xls, 2013-10-11
      parameter Q.Acceleration Acceleration=1*U.m/U.s^2 "Acceleration";
      parameter Q.Amount Amount=1*U.C "Amount";
      parameter Q.AmountReciprocal AmountReciprocal=1/U.C
        "Reciprocal of amount";
      parameter Q.Angle Angle=1*U.rad "Angle";
      parameter Q.Angle2 Angle2=1*U.sr "Solid angle";
      parameter Q.Area Area=1*U.m^2 "Area";
      parameter Q.AreaSpecific AreaSpecific=1*U.m^2/U.mol "Specific area";
      parameter Q.Capacitance Capacitance=1*U.F "Capacitance";
      parameter Q.Concentration Concentration=1*U.C/U.m^3 "Concentration";
      parameter Q.ConcentrationRate ConcentrationRate=1*U.C/(U.m^3*U.s)
        "Rate of concentration";
      parameter Q.ConductanceElectrical ConductanceElectrical=1*U.S
        "Electrical conductance";
      parameter Q.ConductivityElectrical ConductivityElectrical=1*U.S
        "Electrical conductivity";
      parameter Q.Current Current=1*U.A "Current";
      parameter Q.CurrentAreic CurrentAreic=1*U.A/U.m^2 "Areic current";
      parameter Q.CurrentAreicAbsolute CurrentAreicAbsolute=1*U.A/U.m^2
        "Absolute areic current";
      parameter Q.CurrentRate CurrentRate=1*U.A/U.s "Rate of current";
      parameter Q.Diffusivity Diffusivity=1*U.m^2/U.s "Diffusivity";
      parameter Q.Energy Energy=1*U.J "Energy";
      parameter Q.Fluidity Fluidity=1/(U.Pa*U.s) "Fluidity";
      parameter Q.Force Force=1*U.N "Force";
      parameter Q.ForceSpecific ForceSpecific=1*U.V/U.m "Specific force";
      parameter Q.Frequency Frequency=1*U.rad/U.s "Frequency";
      parameter Q.Inductance Inductance=1*U.H "Inductance";
      parameter Q.Length Length=1*U.m "Length";
      parameter Q.LengthSpecific LengthSpecific=1*U.m/U.C "Specific length";
      parameter Q.MagneticFlux MagneticFlux=1*U.Wb "Magnetic flux";
      parameter Q.MagneticFluxAreic MagneticFluxAreic=1*U.T
        "Areic magnetic flux";
      parameter Q.MagneticFluxReciprocal MagneticFluxReciprocal=1/U.Wb
        "Reciprocal of magnetic flux";
      parameter Q.Mass Mass=1*U.kg "Mass";
      parameter Q.MassSpecific MassSpecific=1*U.micro*U.g/U.C "Specific mass";
      parameter Q.MassVolumic MassVolumic=1*U.kg/U.m^3 "Volumic mass";
      parameter Q.Mobility Mobility=1*U.C*U.s/U.g "Mobility";
      parameter Q.MomentumRotational MomentumRotational=1*U.J*U.s/U.rad
        "Rotational momentum";
      parameter Q.Number Number=1 "Number";
      parameter Q.NumberAbsolute NumberAbsolute=1*U.J/(U.mol*U.K)
        "Absolute number";
      parameter Q.Permeability Permeability=1*U.H/U.m "Permeability";
      parameter Q.Permittivity Permittivity=1*U.F/U.m "Permittivity";
      parameter Q.PermittivityReciprocal PermittivityReciprocal=1*U.m/U.H
        "Reciprocal of permittivity";
      parameter Q.Potential Potential=1*U.V "Potential";
      parameter Q.PotentialAbsolute PotentialAbsolute=1*U.K
        "Absolute potential";
      parameter Q.PotentialPerWavenumber PotentialPerWavenumber=1*U.V*U.m/U.rad
        "Potential per wavenumber";
      parameter Q.PotentialRate PotentialRate=1*U.V/U.s "Rate of potential";
      parameter Q.Power Power=1*U.W "Power";
      parameter Q.PowerArea PowerArea=1*U.W*U.m^2 "Power times area";
      parameter Q.PowerAreic PowerAreic=1*U.W/U.m^2 "Areic power";
      parameter Q.PowerAreicPerPotential4 PowerAreicPerPotential4=1*U.W/(U.m^2*
          U.K^4) "Areic power per 4th power of potential";
      parameter Q.PowerRadiant PowerRadiant=1*U.'cd' "Radiant power";
      parameter Q.Pressure Pressure=1*U.Pa "Pressure";
      parameter Q.PressureAbsolute PressureAbsolute=1*U.Pa "Absolute pressure";
      parameter Q.PressureRate PressureRate=1*U.Pa/U.s "Rate of pressure";
      parameter Q.PressureReciprocal PressureReciprocal=1/U.Pa
        "Reciprocal of pressure";
      parameter Q.ResistanceElectrical ResistanceElectrical=1*U.ohm
        "Electrical resistance";
      parameter Q.Resistivity Resistivity=1*U.m/U.A "Resistivity";
      parameter Q.ResistivityMaterial ResistivityMaterial=1*U.s/U.m^2
        "Material resistivity";
      parameter Q.Time Time=1*U.s "Time";
      parameter Q.TimeAbsolute TimeAbsolute=1*U.s "Absolute time";
      parameter Q.TimeLineic TimeLineic=1*U.s/U.m "Lineic time";
      parameter Q.Velocity Velocity=1*U.m/U.s "Velocity";
      parameter Q.Velocity2 Velocity2=1*U.Sv "Squared velocity";
      parameter Q.Volume Volume=1*U.m^3 "Volume";
      parameter Q.VolumeRate VolumeRate=1*U.m^3/U.s "Rate of volume";
      parameter Q.VolumeSpecific VolumeSpecific=1*U.m^3/U.C "Specific volume";
      parameter Q.VolumeSpecificAbsolute VolumeSpecificAbsolute=1*U.m^3/U.C
        "Absolute specific volume";
      parameter Q.VolumeSpecificRate VolumeSpecificRate=1*U.m^3/(U.C*U.s)
        "Rate of specific volume";
      parameter Q.Wavenumber Wavenumber=1*U.rad/U.m "Wavenumber";
      // -------- end from FCSys/Resources/quantities.xls
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{
                -100,-100},{100,100}})), Icon(coordinateSystem(
              preserveAspectRatio=false, extent={{-100,-100},{100,100}})));

    end ExampleModel;

  end Examples;

  extends Modelica.Icons.TypesPackage;
  import Modelica.Icons.TypeReal;

  // Generated from FCSys/Resources/quantities.xls, 2013-10-11
  type Acceleration = TypeReal (final unit="L/T2");
  type Amount = TypeReal (final unit="N", min=0);
  type AmountReciprocal = TypeReal (final unit="1/N", min=0)
    "Reciprocal of amount";
  type Angle = TypeReal (final unit="A");
  type Angle2 = TypeReal (final unit="A2") "Solid angle";
  type Area = TypeReal (final unit="L2", min=0);
  type AreaSpecific = TypeReal (final unit="L2/N", min=0) "Specific area";
  type Capacitance = TypeReal (final unit="N2.T2/(L2.M)", min=0);
  type Concentration = TypeReal (final unit="N/L3", min=0);
  type ConcentrationRate = TypeReal (final unit="N/(L3.T)")
    "Rate of concentration";
  type ConductanceElectrical = TypeReal (final unit="N2.T/(L2.M)", min=0)
    "Electrical conductance";
  type ConductivityElectrical = TypeReal (final unit="N2.T/(L3.M)", min=0)
    "Electrical conductivity";
  type Current = TypeReal (final unit="N/T");
  type CurrentAreic = TypeReal (final unit="N/(L2.T)") "Areic current";
  type CurrentAreicAbsolute = TypeReal (final unit="N/(L2.T)", min=0)
    "Absolute areic current";
  type CurrentRate = TypeReal (final unit="N/T2") "Rate of current";
  type Diffusivity = TypeReal (final unit="L2/T", min=0);
  type Energy = TypeReal (final unit="L2.M/T2");
  type Fluidity = TypeReal (final unit="L.T/M", min=0);
  type Force = TypeReal (final unit="L.M/T2");
  type ForceSpecific = TypeReal (final unit="L.M/(N.T2)") "Specific force";
  type Frequency = TypeReal (final unit="A/T");
  type Inductance = TypeReal (final unit="L2.M/N2", min=0);
  type Length = TypeReal (final unit="L", min=0);
  type LengthSpecific = TypeReal (final unit="L/N", min=0) "Specific length";
  type MagneticFlux = TypeReal (final unit="L2.M/(A.N.T)") "Magnetic flux";
  type MagneticFluxAreic = TypeReal (final unit="M/(A.N.T)")
    "Areic magnetic flux";
  type MagneticFluxReciprocal = TypeReal (final unit="A.N.T/(L2.M)")
    "Reciprocal of magnetic flux";
  type Mass = TypeReal (final unit="M", min=0);
  type MassSpecific = TypeReal (final unit="M/N", min=0) "Specific mass";
  type MassVolumic = TypeReal (final unit="M/L3", min=0) "Volumic mass";
  type Mobility = TypeReal (final unit="N.T/M", min=0);
  type MomentumRotational = TypeReal (final unit="L2.M/(A.T)")
    "Rotational momentum";
  type Number = TypeReal (final unit="1");
  type NumberAbsolute = TypeReal (final unit="1", min=0) "Absolute number";
  type Permeability = TypeReal (final unit="L.M/N2", min=0);
  type Permittivity = TypeReal (final unit="N2.T2/(L3.M)", min=0);
  type PermittivityReciprocal = TypeReal (final unit="L3.M/(N2.T2)", min=0)
    "Reciprocal of permittivity";
  type Potential = TypeReal (final unit="L2.M/(N.T2)");
  type PotentialAbsolute = TypeReal (final unit="L2.M/(N.T2)", min=0)
    "Absolute potential";
  type PotentialPerWavenumber = TypeReal (final unit="L3.M/(A.N.T2)")
    "Potential per wavenumber";
  type PotentialRate = TypeReal (final unit="L2.M/(N.T3)") "Rate of potential";
  type Power = TypeReal (final unit="L2.M/T3");
  type PowerArea = TypeReal (final unit="L4.M/T3") "Power times area";
  type PowerAreic = TypeReal (final unit="M/T3") "Areic power";
  type PowerAreicPerPotential4 = TypeReal (final unit="M.T5/L8")
    "Areic power per 4th power of potential";
  type PowerRadiant = TypeReal (final unit="L2.M/(A2.T3)") "Radiant power";
  type Pressure = TypeReal (final unit="M/(L.T2)");
  type PressureAbsolute = TypeReal (final unit="M/(L.T2)", min=0)
    "Absolute pressure";
  type PressureRate = TypeReal (final unit="M/(L.T3)") "Rate of pressure";
  type PressureReciprocal = TypeReal (final unit="L.T2/M", min=0)
    "Reciprocal of pressure";
  type ResistanceElectrical = TypeReal (final unit="L2.M/(N2.T)", min=0)
    "Electrical resistance";
  type Resistivity = TypeReal (final unit="L.T/N", min=0);
  type ResistivityMaterial = TypeReal (final unit="M/(N.T)", min=0)
    "Material resistivity";
  type Time = TypeReal (final unit="T");
  type TimeAbsolute = TypeReal (final unit="T", min=0) "Absolute time";
  type TimeLineic = TypeReal (final unit="T/L") "Lineic time";
  type Velocity = TypeReal (final unit="L/T");
  type Velocity2 = TypeReal (final unit="L2/T2") "Squared velocity";
  type Volume = TypeReal (final unit="L3", min=0);
  type VolumeRate = TypeReal (final unit="L3/T") "Rate of volume";
  type VolumeSpecific = TypeReal (final unit="L3/N") "Specific volume";
  type VolumeSpecificAbsolute = TypeReal (final unit="L3/N", min=0)
    "Absolute specific volume";
  type VolumeSpecificRate = TypeReal (final unit="L3/(N.T)")
    "Rate of specific volume";
  type Wavenumber = TypeReal (final unit="A/L");

  // -------- end from FCSys/Resources/quantities.xls

  // Aliases that imply special display units
  type CapacityThermal = Amount (displayUnit="J/K") "Thermal capacity";
  type CapacityThermalSpecific = NumberAbsolute (displayUnit="J/(mol.K)")
    "Specific thermal capacity";
  type PotentialChemical = Potential (displayUnit="J/mol") "Chemical potential";
  type Temperature = Potential (displayUnit="K");
  type TemperatureAbsolute = PotentialAbsolute (displayUnit="degC")
    "Absolute temperature";
  type TemperatureRate = PotentialRate (displayUnit="K/s")
    "Rate of temperature";
  type ResistivityThermal = Resistivity (displayUnit="m.K/W")
    "Thermal resistivity";
  type Conductance = Current (displayUnit="W/K") "Conductance";
  annotation (Documentation(info="<html><p>In this package the <code>unit</code> attribute of each <code>Real</code> variable actually denotes the
  dimension.<sup><a href=\"#fn1\" id=\"ref1\">1</a></sup>  The dimensions are
  angle (A), length (L), mass (M), particle number (N), and time (T).  These
  are combined according to the rules established for unit strings
  [<a href=\"modelica://FCSys.UsersGuide.References\">Modelica3.2</a>, p. 210].  In
  <a href=\"modelica://FCSys.FCSys\">FCSys</a>, temperature and charge are considered to be derived dimensions
  (see the <a href=\"modelica://FCSys.Units\">Units</a> package).</p>

  <p>The <code>quantity</code> attribute is not used since the type <i>is</i> the quantity.  The <code>displayUnit</code> attribute is
  only used for quantities that imply a certain display unit.</p>

  <p>It is helpful to check that the terms of each model equation have the same dimension.  Fortunately, methods for unit checking
  have been established [<a href=\"modelica://FCSys.UsersGuide.References\">Mattsson2008</a>,
  <a href=\"modelica://FCSys.UsersGuide.References\">Broman2008</a>,
  <a href=\"modelica://FCSys.UsersGuide.References\">Aronsson2009</a>] and can, in theory, be applied to
  dimension checking instead.</p>

  <p>The quantities are generally named with adjectives following the noun so that the
  quantities are grouped when alphabetized.</p>

  <p>The <a href=\"modelica://FCSys.Quantities\">Quantities</a> package is abbreviated as <code>Q</code> throughout
  the rest of <a href=\"modelica://FCSys.FCSys\">FCSys</a>.</p>

    <hr>

    <small>
    <p id=\"fn1\">1. This misnomer is necessary because <code>Real</code> variables do not have a <code>dimension</code>
    attribute.<a href=\"#ref1\" title=\"Jump back to footnote 1 in the text.\">&#8629;</a></p>
    </small>

  <p><b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
Copyright 2007&ndash;2013, <a href=\"http://www.gtrc.gatech.edu/\">Georgia Tech Research Corporation</a>.</p>

<p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.License\">
FCSys.UsersGuide.License</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p>
</html>"));

end Quantities;
