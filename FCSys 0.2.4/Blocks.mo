within FCSys;
package Blocks "Imperative models (inputs and outputs only)"
  extends Modelica.Icons.Package;
  package UnitConversions
    "Blocks to convert to or from quantities expressed in units"
    extends Modelica.Icons.Package;
    block From_degC
      "Convert from temperature in degree Celsius to thermodynamic temperature"
      extends Partial(y(final unit="l2.m/(N.T2)", displayUnit="U.K"));

    equation
      y = U.from_degC(u);
      annotation (Documentation(info="<html>
    <p>Please see the documentation of the
    <a href=\"modelica://FCSys.Blocks.UnitConversions\">UnitConversions</a> package.</p>
    </html>"));
    end From_degC;

    block To_degC
      "Convert from thermodynamic temperature to temperature in degree Celsius"
      extends Partial(u(final unit="l2.m/(N.T2)", displayUnit="U.K"));

    equation
      y = U.to_degC(u);
      annotation (Documentation(info="<html>
    <p>Please see the documentation of the
    <a href=\"modelica://FCSys.Blocks.UnitConversions\">UnitConversions</a> package.</p>
    </html>"));
    end To_degC;

    block From_kPag
      "Convert from gauge pressure in kilopascals to absolute pressure"
      extends Partial(y(final unit="m/(l.T2)"));

    equation
      y = U.from_kPag(u);
      annotation (Documentation(info="<html>
    <p>Please see the documentation of the
    <a href=\"modelica://FCSys.Blocks.UnitConversions\">UnitConversions</a> package.</p>
    </html>"));
    end From_kPag;

    block To_kPag "Convert absolute pressure to gauge pressure in kilopascals"
      extends Partial(u(final unit="m/(l.T2)"));

    equation
      y = U.to_kPag(u);
      annotation (Documentation(info="<html>
    <p>Please see the documentation of the
    <a href=\"modelica://FCSys.Blocks.UnitConversions\">UnitConversions</a> package.</p>
    </html>"));
    end To_kPag;

    block UnitConversion
      "Generic block to convert to or from quantities expressed in units"
      extends Partial;

      replaceable function f = U.from_degC constrainedby
        Modelica.SIunits.Conversions.ConversionIcon "Conversion function"
        annotation (__Dymola_choicesFromPackage=true);

    equation
      y = f(u);
      annotation (Documentation(info="<html>
    <p>Please see the documentation of the
    <a href=\"modelica://FCSys.Blocks.UnitConversions\">UnitConversions</a> package.</p>
    </html>"));
    end UnitConversion;

  protected
    partial block Partial
      "Base block to convert to or from quantities expressed in units"
      // extends FCSys.Icons.Blocks.ContinuousShort;
      extends FCSys.Icons.Names.Top2;

      Connectors.RealInput u "Real input for value in source representation"
        annotation (Placement(transformation(extent={{-120,-10},{-100,10}},
              rotation=0), iconTransformation(extent={{-120,-10},{-100,10}})));
      Connectors.RealOutput y
        "Real output for quantity in desired representation" annotation (
          Placement(transformation(extent={{100,-10},{120,10}}, rotation=0),
            iconTransformation(extent={{100,-10},{120,10}})));
      annotation (Documentation(info="<html>
<p>Please see the documentation of the
  <a href=\"modelica://FCSys.Blocks.UnitConversions\">UnitConversions</a> package.</p></html>"),
          Icon(graphics={
            Rectangle(
              extent={{-100,40},{100,-40}},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              lineColor={0,0,0}),
            Line(points={{-70,0},{10,0}}, color={191,0,0}),
            Polygon(
              points={{70,0},{10,20},{10,-20},{70,0}},
              lineColor={191,0,0},
              fillColor={191,0,0},
              fillPattern=FillPattern.Solid)}));

    end Partial;
    annotation (Documentation(info="<html>
  <p>This package contains blocks to convert quantities from the unit system in
  <a href=\"modelica://FCSys\">FCSys</a> (see <a href=\"modelica://FCSys.Units\">FCSys.Units</a>)
  to quantities expressed in units or vice versa.
  The <a href=\"modelica://FCSys.Blocks.UnitConversions.UnitConversion\">UnitConversion</a> block
  has a parameter that configures it
  to convert to or from any of the supported units.
  All other blocks convert to or from a predefined unit.</p>

<p>Blocks are only included for
  units that involve offsets or other functions besides simple scaling.
  For conversions that require just a scaling factor, it is best to use
the <a href=\"modelica://Modelica.Blocks.Math.Gain\">Modelica.Blocks.Math.Gain</a>
block with the appropriate factors from <a href=\"modelica://FCSys.Units\">Units</a> package
(<code>U</code>).  For example, to convert from potential in volts use:<br>
<code>Modelica.Blocks.Math.Gain from_V(k=U.V);</code><br>
To convert to current in amperes use:<br>
<code>Modelica.Blocks.Math.Gain to_A(k=1/U.A);</code></p>
</html>"));

  end UnitConversions;

  annotation (Documentation(info="
<html>
  <p><b>Licensed by the Hawaii Natural Energy Institute under the Modelica License 2</b><br>
Copyright 2007&ndash;2014, <a href=\"http://www.hnei.hawaii.edu/\">Hawaii Natural Energy Institute</a> and <a href=\"http://www.gtrc.gatech.edu/\">Georgia Tech Research Corporation</a>.</p>

<p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.License\">
FCSys.UsersGuide.License</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p>
</html>"), Icon(graphics={
        Rectangle(
          origin={0,35.1488},
          fillColor={255,255,255},
          extent={{-30.0,-20.1488},{30.0,20.1488}}),
        Line(origin={51.25,0}, points={{-21.25,35.0},{13.75,35.0},{13.75,-35.0},
              {-6.25,-35.0}}),
        Rectangle(
          origin={0,-34.8512},
          fillColor={255,255,255},
          extent={{-30.0,-20.1488},{30.0,20.1488}}),
        Line(origin={-51.25,0}, points={{21.25,-35.0},{-13.75,-35.0},{-13.75,
              35.0},{6.25,35.0}}),
        Polygon(
          origin={-40,35},
          pattern=LinePattern.None,
          fillPattern=FillPattern.Solid,
          points={{10.0,0.0},{-5.0,5.0},{-5.0,-5.0}}),
        Polygon(
          origin={40,-35},
          pattern=LinePattern.None,
          fillPattern=FillPattern.Solid,
          points={{-10.0,0.0},{5.0,5.0},{5.0,-5.0}})}));

end Blocks;
