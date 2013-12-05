within FCSys;
package Icons "Icons to annotate and represent classes"
  extends Modelica.Icons.IconsPackage;

  package Blocks "Icons for blocks (imperative or causal models)"
    extends Modelica.Icons.Package;
    partial class Continuous "Icon for a continuous-time block"
      // extends Names.Middle;
      // This has been modified from Modelica.Blocks.Interfaces.BlockIcon.
      annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}), graphics={Rectangle(
                  extent={{-100,-100},{100,100}},
                  lineColor={0,0,127},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),Text(
                  extent={{-100,-20},{100,20}},
                  textString="%name",
                  lineColor={0,0,0})}));

    end Continuous;

    partial class ContinuousShort "Short icon for a continuous block"
      extends Names.Middle;
      annotation (Icon(graphics={Rectangle(
                  extent={{-100,40},{100,-40}},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  lineColor={0,0,0}),Text(
                  extent={{-100,-20},{100,20}},
                  textString="%name",
                  lineColor={0,0,0})}));

    end ContinuousShort;

    partial class ContinuousShortWide
      "Short and wide icon for a continuous block"
      extends Names.Middle;
      annotation (Icon(graphics={Rectangle(
                  extent={{-120,40},{120,-40}},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  lineColor={0,0,0}),Text(
                  extent={{-120,-20},{120,20}},
                  textString="%name",
                  lineColor={0,0,0})}));

    end ContinuousShortWide;

    partial class Discrete "Icon for a discrete-time block"
      extends Names.Top5;
      // This has been modified from
      // Modelica.Blocks.Interfaces.DiscreteBlockIcon.
      annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}), graphics={Rectangle(
                  extent={{-100,-100},{100,100}},
                  lineColor={0,0,127},
                  fillColor={223,223,159},
                  fillPattern=FillPattern.Solid)}));

    end Discrete;

  end Blocks;

  package Conditions "Icons for conditions"
    extends Modelica.Icons.Package;
    partial class Pair "Icon for a two-connector boundary condition"
      // extends Names.Middle;
      annotation (Icon(graphics={Rectangle(
                  extent={{-100,100},{100,-100}},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Line(
                  points={{-100,100},{100,100}},
                  pattern=LinePattern.None,
                  smooth=Smooth.None),Line(
                  points={{-100,-100},{-100,100}},
                  color={0,0,0},
                  smooth=Smooth.None,
                  pattern=LinePattern.Dash),Text(
                  extent={{-150,-20},{150,20}},
                  textString="%name",
                  lineColor={0,0,0}),Line(
                  points={{-100,-100},{100,-100}},
                  pattern=LinePattern.None,
                  smooth=Smooth.None),Line(
                  points={{100,-100},{100,100}},
                  color={0,0,0},
                  smooth=Smooth.None,
                  pattern=LinePattern.Dash)}));

    end Pair;

    partial class PairShort "Short icon for a two-connector boundary condition"
      // extends Names.Middle;
      annotation (Icon(graphics={Rectangle(
                  extent={{-100,40},{100,-40}},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Line(
                  points={{-100,40},{100,40}},
                  pattern=LinePattern.None,
                  smooth=Smooth.None),Line(
                  points={{-100,-40},{-100,40}},
                  color={0,0,0},
                  smooth=Smooth.None,
                  pattern=LinePattern.Dash),Text(
                  extent={{-150,-20},{150,20}},
                  textString="%name",
                  lineColor={0,0,0}),Line(
                  points={{-100,-40},{100,-40}},
                  pattern=LinePattern.None,
                  smooth=Smooth.None),Line(
                  points={{100,-40},{100,40}},
                  color={0,0,0},
                  smooth=Smooth.None,
                  pattern=LinePattern.Dash)}));

    end PairShort;

    partial class Single "Icon for a single-connector boundary condition"
      // extends Names.Middle;
      annotation (Icon(graphics={Rectangle(
                  extent={{-100,100},{100,-100}},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Line(
                  points={{-100,-100},{-100,100},{100,100},{100,-100}},
                  pattern=LinePattern.None,
                  smooth=Smooth.None),Line(
                  points={{-100,-100},{100,-100}},
                  color={0,0,0},
                  smooth=Smooth.None,
                  pattern=LinePattern.Dash),Text(
                  extent={{-100,-20},{100,20}},
                  textString="%name",
                  lineColor={0,0,0})}));

    end Single;

    partial class SingleShort
      "Short icon for a single-connector boundary condition"
      // extends Names.Middle;
      annotation (Icon(graphics={
            Rectangle(
              extent={{-100,40},{100,-40}},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Line(
              points={{-100,-40},{-100,40},{100,40},{100,-40}},
              pattern=LinePattern.None,
              smooth=Smooth.None),
            Line(
              points={{-100,-40},{100,-40}},
              color={0,0,0},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),
            Text(
              extent={{-100,-20},{100,20}},
              textString="%name",
              lineColor={0,0,0})}));

    end SingleShort;

  end Conditions;

  package ChemistryPackage "Icon for packages containing chemical models"
    extends Modelica.Icons.Package;
    annotation (Icon(graphics={
          Polygon(
            points={{-20,70},{-20,-20},{-20,-20},{-70,-60},{-70,-80},{70,-80},{
                70,-60},{20,-20},{20,-20},{20,70},{20,70},{30,80},{-30,80},{-20,
                70},{-20,70}},
            smooth=Smooth.Bezier,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Line(
            points={{20,56},{0,56}},
            color={175,175,175},
            smooth=Smooth.None),
          Line(
            points={{20,36},{0,36}},
            color={175,175,175},
            smooth=Smooth.None),
          Line(
            points={{20,16},{0,16}},
            color={175,175,175},
            smooth=Smooth.None),
          Line(
            points={{20,-4},{0,-4}},
            color={175,175,175},
            smooth=Smooth.None),
          Polygon(
            points={{-20,70},{-20,-20},{-20,-20},{-70,-60},{-70,-80},{70,-80},{
                70,-60},{20,-20},{20,-20},{20,70},{20,70},{30,80},{-30,80},{-20,
                70},{-20,70}},
            lineColor={64,64,64},
            smooth=Smooth.Bezier)}));

  end ChemistryPackage;

  package Names "Icons labeled with the name of the class at various positions"
    extends Modelica.Icons.Package;

    partial class Top10

      annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}), graphics={Text(
                  extent={{-100,200},{100,240}},
                  textString="%name",
                  lineColor={0,0,0})}));

    end Top10;

    partial class Top9

      annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}), graphics={Text(
                  extent={{-100,180},{100,220}},
                  textString="%name",
                  lineColor={0,0,0})}));

    end Top9;

    partial class Top8

      annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}), graphics={Text(
                  extent={{-100,160},{100,200}},
                  textString="%name",
                  lineColor={0,0,0})}));

    end Top8;

    partial class Top7

      annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}), graphics={Text(
                  extent={{-100,140},{100,180}},
                  textString="%name",
                  lineColor={0,0,0})}));

    end Top7;

    partial class Top6

      annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}), graphics={Text(
                  extent={{-100,120},{100,160}},
                  textString="%name",
                  lineColor={0,0,0})}));

    end Top6;

    partial class Top5

      annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}), graphics={Text(
                  extent={{-100,100},{100,140}},
                  textString="%name",
                  lineColor={0,0,0})}));

    end Top5;

    partial class Top4

      annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}), graphics={Text(
                  extent={{-100,80},{100,120}},
                  textString="%name",
                  lineColor={0,0,0})}));

    end Top4;

    partial class Top3

      annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}), graphics={Text(
              extent={{-100,60},{100,100}},
              textString="%name",
              lineColor={0,0,0})}));

    end Top3;

    partial class Top2

      annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}), graphics={Text(
                  extent={{-100,40},{100,80}},
                  textString="%name",
                  lineColor={0,0,0})}));

    end Top2;

    partial class Top1

      annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}), graphics={Text(
              extent={{-100,20},{100,60}},
              textString="%name",
              lineColor={0,0,0})}));

    end Top1;

    partial class Middle

      annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}), graphics={Text(
                  extent={{-100,-20},{100,20}},
                  textString="%name",
                  lineColor={0,0,0})}));

    end Middle;

  end Names;

  partial class Cell "Icon for a cell"
    // extends Names.Top3;
    annotation (Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          initialScale=0.1), graphics={
          Polygon(
            points={{-4,52},{-14,42},{6,42},{16,52},{-4,52}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-10,52},{-20,42},{-14,42},{-4,52},{-10,52}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,0,0}),
          Rectangle(
            extent={{6,42},{12,-52}},
            fillPattern=FillPattern.Solid,
            fillColor={0,0,0},
            pattern=LinePattern.None),
          Polygon(
            points={{16,52},{6,42},{12,42},{22,52},{16,52}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,0,0}),
          Line(
            points={{-40,42},{-40,-52}},
            pattern=LinePattern.None,
            smooth=Smooth.None),
          Polygon(
            points={{-46,64},{-66,44},{-46,44},{-26,64},{-46,64}},
            lineColor={0,0,0},
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-39.6277,31.7996},{-67.912,17.6573}},
            lineColor={0,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            rotation=45,
            fillColor={255,255,255},
            origin={56.5067,67.5353}),
          Rectangle(
            extent={{-14,42},{6,-52}},
            lineColor={0,0,0},
            fillPattern=FillPattern.VerticalCylinder,
            fillColor={255,255,255}),
          Line(points={{-30,52},{32,52}}, color={0,0,0}),
          Rectangle(
            extent={{-5.21738,-5.21961},{-33.5017,-33.5041}},
            lineColor={0,0,170},
            fillPattern=FillPattern.HorizontalCylinder,
            rotation=45,
            fillColor={0,0,240},
            origin={31.9983,69.3803}),
          Rectangle(
            extent={{12,42},{52,-52}},
            lineColor={0,0,170},
            fillPattern=FillPattern.VerticalCylinder,
            fillColor={0,0,240}),
          Polygon(
            points={{-26,64},{-46,44},{-46,-64},{-26,-44},{-26,64}},
            lineColor={0,0,0},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-5.21774,-5.2196},{-33.502,-33.5042}},
            lineColor={170,0,0},
            fillPattern=FillPattern.HorizontalCylinder,
            rotation=45,
            fillColor={240,0,0},
            origin={-30.001,79.3803}),
          Rectangle(
            extent={{-60,42},{-20,-52}},
            lineColor={179,0,0},
            fillPattern=FillPattern.VerticalCylinder,
            fillColor={240,0,0}),
          Rectangle(
            extent={{-76.648,66.211},{-119.073,52.0689}},
            lineColor={95,95,95},
            fillPattern=FillPattern.HorizontalCylinder,
            rotation=45,
            fillColor={135,135,135},
            origin={65.0166,81.3801}),
          Rectangle(
            extent={{-76.648,66.211},{-119.073,52.0689}},
            lineColor={95,95,95},
            fillPattern=FillPattern.HorizontalCylinder,
            rotation=45,
            fillColor={135,135,135},
            origin={157,81}),
          Rectangle(
            extent={{26,44},{46,-64}},
            lineColor={95,95,95},
            fillPattern=FillPattern.VerticalCylinder,
            fillColor={135,135,135}),
          Polygon(
            points={{-26,64},{-26,52},{-30,52},{-30,60},{-26,64}},
            smooth=Smooth.None,
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None,
            lineColor={0,0,0}),
          Ellipse(
            extent={{-44,62},{-36,58}},
            lineColor={135,135,135},
            fillColor={250,0,0},
            fillPattern=FillPattern.Sphere),
          Ellipse(
            extent={{36,50},{44,46}},
            lineColor={135,135,135},
            fillColor={0,0,240},
            fillPattern=FillPattern.Sphere),
          Polygon(
            points={{-26,64},{-26,52},{-30,52},{-40,42},{-40,-52},{-34,-52},{-46,
                -64},{-46,44},{-26,64}},
            smooth=Smooth.None,
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None,
            lineColor={0,0,0}),
          Polygon(
            points={{66,64},{46,44},{46,-64},{66,-44},{66,64}},
            lineColor={0,0,0},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-20,42},{-14,-52}},
            fillPattern=FillPattern.Solid,
            fillColor={0,0,0},
            pattern=LinePattern.None),
          Rectangle(extent={{26,44},{46,-64}}, lineColor={0,0,0}),
          Rectangle(
            extent={{-66,44},{-46,-64}},
            lineColor={95,95,95},
            fillPattern=FillPattern.VerticalCylinder,
            fillColor={135,135,135}),
          Rectangle(extent={{-66,44},{-46,-64}}, lineColor={0,0,0}),
          Line(
            points={{-34,-52},{-46,-64}},
            color={0,0,0},
            smooth=Smooth.None),
          Polygon(
            points={{66,74},{66,64},{46,64},{34,52},{-26,52},{-26,64},{-46,64},
                {-46,74},{66,74}},
            smooth=Smooth.None,
            fillPattern=FillPattern.Solid,
            fillColor={255,255,255},
            pattern=LinePattern.None),
          Polygon(
            points={{46,64},{66,64},{46,44},{26,44},{46,64}},
            lineColor={0,0,0},
            smooth=Smooth.None),
          Polygon(
            points={{-46,64},{-26,64},{-46,44},{-66,44},{-46,64}},
            lineColor={0,0,0},
            smooth=Smooth.None),
          Line(points={{26,42},{-40,42},{-30,52},{34,52}}, color={0,0,0}),
          Line(
            points={{-26,64},{-26,52}},
            color={0,0,0},
            smooth=Smooth.None),
          Line(
            points={{-40,42},{-40,-52},{26,-52}},
            color={0,0,0},
            smooth=Smooth.None),
          Text(
            extent={{-100,60},{100,100}},
            textString="%name",
            lineColor={0,0,0})}), Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          initialScale=0.1)));

  end Cell;

  partial record Record "Icon for records"
    extends Names.Top3;
    annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
              {100,100}}), graphics={
          Rectangle(
            origin={0.0,-25.0},
            lineColor={64,64,64},
            fillColor={255,215,136},
            fillPattern=FillPattern.Solid,
            extent={{-100.0,-75.0},{100.0,75.0}},
            radius=25.0),
          Line(points={{-100.0,0.0},{100.0,0.0}}, color={64,64,64}),
          Line(
            origin={0.0,-50.0},
            points={{-100.0,0.0},{100.0,0.0}},
            color={64,64,64}),
          Line(
            origin={0.0,-25.0},
            points={{0.0,75.0},{0.0,-75.0}},
            color={64,64,64})}), Documentation(info="<html>
<p>
This icon is indicates a record.</p>
</html>"));

  end Record;
  annotation (Documentation(info="
<html>
  <p><b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
Copyright 2007&ndash;2013, <a href=\"http://www.gtrc.gatech.edu/\">Georgia Tech Research Corporation</a>.</p>

<p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.License\">
FCSys.UsersGuide.License</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p>
</html>"));

end Icons;
