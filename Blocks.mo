within FCSys;
package Blocks "Imperative models (inputs and/or outputs only)"
  extends Modelica.Icons.Package;

  package UnitConversions
    "Blocks to convert to or from quantities expressed in units"
    extends Modelica.Icons.Package;
    block From_degC
      "Convert from temperature in degree Celsius to thermodynamic temperature"
      extends BaseClasses.PartialUnitConversion(y(final unit="l2.m/(N.T2)",
            displayUnit="U.K"));

    equation
      y = U.from_degC(u);
      annotation (Documentation(info="<html>
    <p>See the documentation in the
    <a href=\"modelica://FCSys.Blocks.UnitConversions\">UnitConversions</a> package.</p>
    </html>"));
    end From_degC;

    block To_degC
      "Convert from thermodynamic temperature to temperature in degree Celsius"
      extends BaseClasses.PartialUnitConversion(u(final unit="l2.m/(N.T2)",
            displayUnit="U.K"));

    equation
      y = U.to_degC(u);
      annotation (Documentation(info="<html>
    <p>See the documentation in the
    <a href=\"modelica://FCSys.Blocks.UnitConversions\">UnitConversions</a> package.</p>
    </html>"));
    end To_degC;

    block From_kPag
      "Convert from gauge pressure in kilopascals to absolute pressure"
      extends BaseClasses.PartialUnitConversion(y(final unit="m/(l.T2)"));

    equation
      y = U.from_kPag(u);
      annotation (Documentation(info="<html>
    <p>See the documentation in the
    <a href=\"modelica://FCSys.Blocks.UnitConversions\">UnitConversions</a> package.</p>
    </html>"));
    end From_kPag;

    block To_kPag "Convert absolute pressure to gauge pressure in kilopascals"
      extends BaseClasses.PartialUnitConversion(u(final unit="m/(l.T2)"));

    equation
      y = U.to_kPag(u);
      annotation (Documentation(info="<html>
    <p>See the documentation in the
    <a href=\"modelica://FCSys.Blocks.UnitConversions\">UnitConversions</a> package.</p>
    </html>"));
    end To_kPag;

    block UnitConversion
      "Generic block to convert to or from quantities expressed in units"
      extends BaseClasses.PartialUnitConversion;

      replaceable function f = U.from_degC "Conversion function" annotation (
          __Dymola_choicesFromPackage=true);

    equation
      y = f(u);
      annotation (Documentation(info="<html>
    <p>See the documentation in the
    <a href=\"modelica://FCSys.Blocks.UnitConversions\">UnitConversions</a> package.</p>
    </html>"));
    end UnitConversion;

    package BaseClasses "Base classes (not generally for direct use)"
      extends Modelica.Icons.BasesPackage;

      partial block PartialUnitConversion
        "Partial block to convert to or from quantities expressed in units"
        //extends FCSys.BaseClasses.Icons.Blocks.ContinuousShort;
        extends FCSys.BaseClasses.Icons.Names.Top2;

        Connectors.RealInput u "Real input for value in source representation"
          annotation (Placement(transformation(extent={{-120,-10},{-100,10}},
                rotation=0), iconTransformation(extent={{-120,-10},{-100,10}})));
        Connectors.RealOutput y
          "Real output for quantity in desired representation" annotation (
            Placement(transformation(extent={{100,-10},{120,10}}, rotation=0),
              iconTransformation(extent={{100,-10},{120,10}})));
        annotation (Documentation(info="<html>
<p>See the documentation in the
  <a href=\"modelica://FCSys.Blocks.UnitConversions\">UnitConversions</a> package.</p></html>"),
            Icon(graphics={Rectangle(
                      extent={{-100,40},{100,-40}},
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Solid,
                      lineColor={0,0,0}),Line(points={{-70,0},{10,0}}, color={
                191,0,0}),Polygon(
                      points={{70,0},{10,20},{10,-20},{70,0}},
                      lineColor={191,0,0},
                      fillColor={191,0,0},
                      fillPattern=FillPattern.Solid)}));

      end PartialUnitConversion;

    end BaseClasses;
    annotation (Documentation(info="<html>
  <p>This package contains blocks to convert quantities from the unit system in
  <a href=\"modelica://FCSys\">FCSys</a> (see <a href=\"modelica://FCSys.Units\">FCSys.Units</a>)
  to quantities expressed in units or vice versa.
  The <a href=\"modelica://FCSys.Blocks.UnitConversions.UnitConversion\">UnitConversion</a> block
  has a parameter that configures it
  to convert to or from any of the supported units.
  All other blocks convert to or from a predefined unit.
  </p>

<p>Blocks are only included for
  units that involve offsets or other functions besides simple scaling.
  For conversions that require just a scaling factor, it is best to use
the <a href=\"modelica://Modelica.Blocks.Math.Gain\">Modelica.Blocks.Math.Gain</a>
block with the appropriate factors from <a href=\"modelica://FCSys.Units\">Units</a> package
(<code>U</code>).  For example, to convert from potential in volts use:
<br><code>Modelica.Blocks.Math.Gain from_V(k=U.V);</code>
<br>To convert to current in amperes use:
<br><code>Modelica.Blocks.Math.Gain to_A(k=1/U.A);</code>
</p>

</html>"));

  end UnitConversions;

end Blocks;
