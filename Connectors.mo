within FCSys;
package Connectors "Declarative and imperative connectors"
  extends Modelica.Icons.InterfacesPackage;

  expandable connector ChemicalBus
    "<html>Bus of <a href=\"modelica://FCSys.Connectors.ChemicalInput\">ChemicalInput</a> and <a href=\"modelica://FCSys.Connectors.ChemicalOutput\">ChemicalOutput</a> connectors (for multiple species)</html>"

    annotation (
      defaultComponentName="chemical",
      Documentation(info="<html><p>There is no minimal set of variables.  Species are included by connecting instances
    of a <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialChemical\">Chemical</a> connector
    (<a href=\"modelica://FCSys.Connectors.ChemicalInput\">ChemicalInput</a> or
    <a href=\"modelica://FCSys.Connectors.ChemicalOutput\">ChemicalOutput</a>).
    </p></html>"),
      Icon(graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={208,104,0},
            fillPattern=FillPattern.Solid,
            fillColor={255,128,0},
            lineThickness=0.5)}),
      Diagram(graphics={Text(
              extent={{-100,36},{100,76}},
              textString="%name",
              lineColor={0,0,0}),Ellipse(
              extent={{-30,30},{30,-30}},
              lineColor={208,104,0},
              fillPattern=FillPattern.Solid,
              fillColor={255,128,0},
              lineThickness=0.5)}));
  end ChemicalBus;

  expandable connector ChemicalBusInternal
    "<html>Internal bus of <a href=\"modelica://FCSys.Connectors.ChemicalInput\">ChemicalInput</a> and <a href=\"modelica://FCSys.Connectors.ChemicalOutput\">ChemicalOutput</a> connectors (for multiple species)</html>"

    annotation (
      defaultComponentPrefixes="protected",
      defaultComponentName="chemical",
      Documentation(info="<html><p>
    This is copy of the <a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a> connector, except that it
    has a smaller icon and a default <code>protected</code> prefix.  For more information, see that connector.</p></html>"),

      Icon(graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={208,104,0},
            fillPattern=FillPattern.Solid,
            fillColor={255,128,0},
            lineThickness=0.5)}),
      Diagram(graphics={Ellipse(
              extent={{-10,10},{10,-10}},
              lineColor={208,104,0},
              fillPattern=FillPattern.Solid,
              fillColor={255,128,0},
              lineThickness=0.5),Text(
              extent={{-100,20},{100,60}},
              textString="%name",
              lineColor={0,0,0})}));

  end ChemicalBusInternal;

  connector ChemicalInput
    "Connector to exchange material while advecting linear momentum and enthalpy, with characteristic data as input"

    extends FCSys.Connectors.BaseClasses.PartialChemical;

    input Q.MassSpecific m "Specific mass";
    input String formula(start="") "Chemical formula";
    // Note:  The start value prevents a warning when checked in Dymola 7.4.
    annotation (
      defaultComponentName="chemical",
      Documentation(info="<html><p>For more information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(graphics={Ellipse(
            extent={{-60,60},{60,-60}},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            lineColor={208,104,0}), Ellipse(extent={{-100,100},{100,-100}},
              lineColor={208,104,0})}),
      Diagram(graphics={Ellipse(
            extent={{-18,18},{18,-18}},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            lineColor={208,104,0})}));

  end ChemicalInput;

  connector ChemicalOutput
    "Connector to exchange material while advecting linear momentum and enthalpy, with characteristic data as output"

    extends FCSys.Connectors.BaseClasses.PartialChemical;

    output Q.MassSpecific m "Specific mass";
    output String formula "Chemical formula";

    annotation (
      defaultComponentName="chemical",
      Documentation(info="<html><p>For more information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(graphics={Ellipse(extent={{-100,100},{100,-100}}, lineColor={208,104,
                0})}),
      Diagram(graphics));

  end ChemicalOutput;

  expandable connector FaceBus
    "<html>Bus of <a href=\"modelica://FCSys.Connectors.Face\">Face</a>, <a href=\"modelica://FCSys.Connectors.FaceX\">FaceX</a>, <a href=\"modelica://FCSys.Connectors.FaceY\">FaceY</a>, <a href=\"modelica://FCSys.Connectors.FaceZ\">FaceZ</a> connectors (for multiple species)</html>"

    annotation (
      defaultComponentName="face",
      Documentation(info="<html><p>There is no minimal set of variables.  Species are included by connecting instances
    of a <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">Face</a> connector
    (<a href=\"modelica://FCSys.Connectors.Face\">Face</a>,
    <a href=\"modelica://FCSys.Connectors.FaceX\">FaceX</a>,
    <a href=\"modelica://FCSys.Connectors.FaceY\">FaceY</a>, or
    <a href=\"modelica://FCSys.Connectors.FaceZ\">FaceZ</a>).  In order to allow
    the subconnectors of the <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">Face</a> connectors
    (instances of <a href=\"modelica://FCSys.Connectors.Material\">Material</a>,
    <a href=\"modelica://FCSys.Connectors.Mechanical\">Mechanical</a>,
    and <a href=\"modelica://FCSys.Connectors.Thermal\">Thermal</a>) to be included independently, those
    subconnectors are connected explicitly.  For example,
    <blockquote>
        <code>
        connect(species.xNegative.material, xNegative.species.material);<br>
        connect(species.xNegative.mechanicalY, xNegative.species.mechanicalY);<br>
        connect(species.xNegative.mechanicalZ, xNegative.species.mechanicalZ);<br>
        connect(species.xNegative.thermal, xNegative.species.thermal);
        </code>
    </blockquote>
    where <code>species</code> is an instance of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a>
    model in which <code>xNegative</code> is a <a href=\"modelica://FCSys.Connectors.FaceX\">FaceX</a> connector.
    Meanwhile, <code>xNegative</code> (not <code>species.xNegative</code>) is an instance of this
    connector (<a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a>).
    </p></html>"),
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            fillColor={191,191,191},
            fillPattern=FillPattern.Solid,
            lineColor={127,127,127},
            lineThickness=0.5)}),
      Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
              100,100}}), graphics={Text(
            extent={{-100,36},{100,76}},
            textString="%name",
            lineColor={0,0,0}), Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={127,127,127},
            fillColor={191,191,191},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5)}));
  end FaceBus;

  expandable connector FaceBusInternal
    "<html>Internal bus of <a href=\"modelica://FCSys.Connectors.Face\">Face</a>, <a href=\"modelica://FCSys.Connectors.FaceX\">FaceX</a>, <a href=\"modelica://FCSys.Connectors.FaceY\">FaceY</a>, or <a href=\"modelica://FCSys.Connectors.FaceZ\">FaceZ</a> connectors (for multiple species)</html>"

    annotation (
      defaultComponentPrefixes="protected",
      defaultComponentName="face",
      Documentation(info="<html><p>
    This is copy of the <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, except that it
    has a smaller icon and a default <code>protected</code> prefix.  For more information, see that connector.</p></html>"),

      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            fillColor={191,191,191},
            fillPattern=FillPattern.Solid,
            lineColor={127,127,127},
            lineThickness=0.5)}),
      Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
              100,100}}), graphics={Ellipse(
              extent={{-10,10},{10,-10}},
              lineColor={127,127,127},
              fillColor={191,191,191},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),Text(
              extent={{-100,20},{100,60}},
              textString="%name",
              lineColor={0,0,0})}));

  end FaceBusInternal;

  connector Face
    "Connector to transport material, displacement, and heat of a single species"
    extends BaseClasses.PartialFace;

    parameter Axis axis=Axis.x "Axis normal to the face";
    parameter Boolean slipX=false if axis <> Axis.x
      "Allow slip along the x axis" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        compact=true,
        group="Boundary conditions (will remove subconnectors)",
        enable=axis <> 1));
    parameter Boolean slipY=false if axis <> Axis.y
      "Allow slip along the y axis" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        compact=true,
        group="Boundary conditions (will remove subconnectors)",
        enable=axis <> 2));
    parameter Boolean slipZ=false if axis <> Axis.z
      "Allow slip along the z axis" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        compact=true,
        group="Boundary conditions (will remove subconnectors)",
        enable=axis <> 3));
    // Note:  Dymola 7.4 doesn't recognize enumerations in the dialog enable
    // option, e.g.,
    //     enable=axis <> Axis.x.
    // Therefore, the values of the enumerations are specified numerically.

    MechanicalTransport mechanicalX if axis <> Axis.x and slipY
      "Subconnector for X-axis displacement";
    MechanicalTransport mechanicalY if axis <> Axis.y and slipY
      "Subconnector for Y-axis displacement";
    MechanicalTransport mechanicalZ if axis <> Axis.z and slipZ
      "Subconnector for Z-axis displacement";

    annotation (
      Documentation(info="<html>
  <p>To specify a no-slip condition at the boundary, set the
  corresponding <code>slipX</code>, <code>slipY</code> or <code>slipZ</code> parameter to <code>false</code>.</p>

  <p>For more information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package and the
    <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">PartialFace</a> connector.</p></html>"),

      Diagram(coordinateSystem(preserveAspectRatio=true,extent={{-100,-100},{
              100,100}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Ellipse(extent={{-100,100},{100,-100}},
              lineColor={127,127,127})}));

  end Face;

  connector FaceX
    "X-axis connector to transport material, displacement, and heat of a single species"
    extends BaseClasses.PartialFace;

    parameter Boolean slipY=false "Allow slip along the y axis" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true,group=
            "Boundary conditions (will remove subconnectors)"));
    parameter Boolean slipZ=false "Allow slip along the z axis" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true, group=
            "Boundary conditions (will remove subconnectors)"));

    MechanicalTransport mechanicalY if slipY
      "Subconnector for Y-axis displacement";
    MechanicalTransport mechanicalZ if slipZ
      "Subconnector for Z-axis displacement";

    annotation (
      Documentation(info="<html>
  <p>To specify a no-slip condition at the boundary, set the
  corresponding <code>slipY</code> or <code>slipZ</code> parameter to <code>false</code>.</p>

  <p>For more information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package and the
    <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">PartialFace</a> connector.</p></html>"),

      Diagram(coordinateSystem(preserveAspectRatio=true,extent={{-100,-100},{
              100,100}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Ellipse(extent={{-100,100},{100,-100}},
              lineColor={127,127,127})}));

  end FaceX;

  connector FaceY
    "Y-axis connector to transport material, displacement, and heat of a single species"
    extends BaseClasses.PartialFace;

    parameter Boolean slipZ=false "Allow slip along the z axis" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true,group=
            "Boundary conditions (will remove subconnectors)"));
    parameter Boolean slipX=false "Allow slip along the x axis" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true, group=
            "Boundary conditions (will remove subconnectors)"));

    MechanicalTransport mechanicalZ if slipZ
      "Subconnector for Z-axis displacement";
    MechanicalTransport mechanicalX if slipX
      "Subconnector for X-axis displacement";

    annotation (
      Documentation(info="<html>
  <p>To specify a no-slip condition at the boundary, set the
  corresponding <code>slipZ</code> or <code>slipX</code> parameter to <code>false</code>.</p>

  <p>For more information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package and the
    <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">PartialFace</a> connector.</p></html>"),

      Diagram(coordinateSystem(preserveAspectRatio=true,extent={{-100,-100},{
              100,100}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Ellipse(extent={{-100,100},{100,-100}},
              lineColor={127,127,127})}));

  end FaceY;

  connector FaceZ
    "Z-axis connector to transport material, displacement, and heat of a single species"
    extends BaseClasses.PartialFace;

    parameter Boolean slipX=false "Allow slip along the x axis" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true,group=
            "Boundary conditions (will remove subconnectors)"));
    parameter Boolean slipY=false "Allow slip along the y axis" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true, group=
            "Boundary conditions (will remove subconnectors)"));

    MechanicalTransport mechanicalX if slipX
      "Subconnector for X-axis displacement";
    MechanicalTransport mechanicalY if slipY
      "Subconnector for Y-axis displacement";

    annotation (
      Documentation(info="<html>
  <p>To specify a no-slip condition at the boundary, set the
  corresponding <code>slipX</code> or <code>slipY</code> parameter to <code>false</code>.</p>

  <p>For more information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package and the
    <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">PartialFace</a> connector.</p></html>"),

      Diagram(coordinateSystem(preserveAspectRatio=true,extent={{-100,-100},{
              100,100}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Ellipse(extent={{-100,100},{100,-100}},
              lineColor={127,127,127})}));

  end FaceZ;

  connector Inert "Connector to exchange linear momentum and heat by diffusion"

    parameter Integer n_lin(
      final min=0,
      final max=3) = 0
      "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>"
      annotation (HideResult=true);

    BaseClasses.MechanicalExchange mechanical(final n_lin=n_lin)
      "Subconnector for linear momentum";
    BaseClasses.Thermal thermal "Subconnector for heat";

    annotation (
      Documentation(info="<html>
    <p>For more information, see the <a href=\"modelica://FCSys.Connectors.BaseClasses.MechanicalExchange\">MechanicalExchange</a>
    subconnector or the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(graphics={
          Text(
            extent={{-100,36},{100,76}},
            textString="%name",
            lineColor={0,0,0}),
          Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={72,90,180},
            fillPattern=FillPattern.Solid,
            fillColor={102,128,255}),
          Ellipse(
            extent={{-18,18},{18,-18}},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            lineColor={72,90,180})}),
      Icon(graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={72,90,180},
            fillPattern=FillPattern.Solid,
            fillColor={102,128,255}), Ellipse(
            extent={{-60,60},{60,-60}},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            lineColor={72,90,180})}));

  end Inert;

  connector InertAmagat
    "<html><a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> connector with additivity of volume</html>"

    extends BaseClasses.MechanicalExchange;
    extends BaseClasses.Thermal;

    // Additivity of volume
    Q.PressureAbsolute p(nominal=1*U.atm) "Pressure";
    flow Q.Volume V(min=-Modelica.Constants.inf, nominal=1*U.cm^3) "Volume";

    annotation (
      defaultComponentName="inert",
      Documentation(info="<html><p>
    The concept of \"additivity of volume\" is defined by
    <a href=\"http://en.wikipedia.org/wiki/Amagat's_law\">Amagat's law</a> or the Law of Partial Volumes, which
    states that the partial extensive volumes of the components of a mixture sum to the total
    extensive volume of the mixture [<a href=\"modelica://FCSys.UsersGuide.References\">Bejan2006</a>, p. 194].
    The specific or molar volumes of the species are each evaluated at the temperature and the total pressure of the
    mixture.</p>

    <p>This concept loses its physical meaning once the species are mixed [<a href=\"modelica://FCSys.UsersGuide.References\">Woo1995</a>].
    If the species are completely mixed, then it is impossible to distinguish their particles and thus
    determine their partial volumes.
    Therefore, the concept is only used to allow distinct phases to exist within the same subregion&mdash;not
    to species within a phase.
    If a system contains only a solid phase and a gas phase, it is assumed that the
    partial volumes of the mixtures are additive and the mixtures exist at the same pressure.  Within
    a phase, the species are mixed according to Dalton's law (see the <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector).</p>

    <p>In order to implement Amagat's law, this connector includes volume (not rate of volume) as a flow variable.
    The effort variable is pressure such that the effort and flow variables are conjugates of
    energy (not power).
    </p>

    <p>See also the <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> and
    <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connectors and
    the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(graphics={
          Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={72,90,180},
            fillPattern=FillPattern.Solid,
            fillColor={102,128,255}),
          Text(
            extent={{-30,20},{30,-20}},
            lineColor={0,0,0},
            textString="A"),
          Text(
            extent={{-100,36},{100,76}},
            textString="%name",
            lineColor={0,0,0})}),
      Icon(graphics={
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={72,90,180},
            fillPattern=FillPattern.Solid,
            fillColor={102,128,255}),
          Ellipse(extent={{-100,100},{100,-100}}, lineColor={72,90,180}),
          Text(
            extent={{-100,80},{100,-90}},
            lineColor={0,0,0},
            textString="A")}));

  end InertAmagat;

  connector InertDalton
    "<html><a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> connector with additivity of pressure</html>"

    extends BaseClasses.MechanicalExchange;
    extends BaseClasses.Thermal;

    // Additivity of pressure
    Q.Volume V(nominal=1*U.cm^3) "Volume";
    flow Q.Pressure p(nominal=1*U.atm) "Pressure";

    annotation (
      defaultComponentName="inert",
      Documentation(info="<html><p>
    The concept of \"additivity of pressure\" is defined by
    <a href=\"http://en.wikipedia.org/wiki/Dalton's_law\">Dalton's law</a> or the Law of Partial Pressures,
    which states that the partial pressures of the components of a mixture sum to the total
    pressure of the mixture [<a href=\"modelica://FCSys.UsersGuide.References\">Bejan2006</a>, p. 192].
    The partial pressures of the species are evaluated at the temperature and the total volume of the
    mixture.</p>

    <p>In order to implement Dalton's law, this connector includes pressure as a flow variable.
    The effort variable is volume such that the effort and flow variables are conjugates of
    energy (not power).
    </p>

    <p>See also the <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> and
    <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> connectors and
    the documentation in the
  <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(graphics={
          Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={72,90,180},
            fillPattern=FillPattern.Solid,
            fillColor={102,128,255}),
          Text(
            extent={{-30,20},{30,-20}},
            lineColor={0,0,0},
            textString="D"),
          Text(
            extent={{-100,36},{100,76}},
            textString="%name",
            lineColor={0,0,0})}),
      Icon(graphics={
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={72,90,180},
            fillPattern=FillPattern.Solid,
            fillColor={102,128,255}),
          Ellipse(extent={{-100,100},{100,-100}}, lineColor={72,90,180}),
          Text(
            extent={{-100,80},{100,-90}},
            lineColor={0,0,0},
            textString="D")}));

  end InertDalton;

  connector InertInternal
    "<html>Internal <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> connector</html>"

    parameter Integer n_lin(
      final min=0,
      final max=3) = 0
      "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>"
      annotation (HideResult=true);
    parameter Boolean uniformVelocity=false "Uniform velocity" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true,group=
            "Boundary conditions (will remove subconnectors)"));
    parameter Boolean isothermal=false "Isothermal" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true,group=
            "Boundary conditions (will remove subconnectors)"));

    BaseClasses.MechanicalExchange mechanical(final n_lin=n_lin) if
      uniformVelocity "Subconnector for linear momentum";
    BaseClasses.Thermal thermal if isothermal "Subconnector for heat";

    annotation (
      defaultComponentPrefixes="protected",
      defaultComponentName="inert",
      Documentation(info="<html><p>
    This is copy of the <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> connector, except that it
    has a smaller icon, a default <code>protected</code> prefix, and the subconnectors are conditional.
    For more information, see that connector.</p></html>"),
      Icon(graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={72,90,180},
            fillPattern=FillPattern.Solid,
            fillColor={102,128,255}),Ellipse(
            extent={{-60,60},{60,-60}},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            lineColor={72,90,180})}),
      Diagram(graphics={Ellipse(
              extent={{-10,10},{10,-10}},
              lineColor={72,90,180},
              fillPattern=FillPattern.Solid,
              fillColor={102,128,255}),Text(
              extent={{-100,20},{100,60}},
              textString="%name",
              lineColor={0,0,0}),Ellipse(
              extent={{-4,4},{4,-4}},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              lineColor={72,90,180})}));
  end InertInternal;

  connector MaterialTransport "Connector to transport material"

    Q.PressureAbsolute p(nominal=1*U.Pa) "Pressure";
    flow Q.Current Ndot(nominal=1*U.A) "Current";

    annotation (
      defaultComponentName="material",
      Documentation(info="<html><p>For more information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(graphics={Text(
              extent={{-100,36},{100,76}},
              textString="%name",
              lineColor={0,0,0}),Rectangle(
              extent={{-30,30},{30,-30}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),Line(
              points={{-30,0},{30,0}},
              color={0,0,0},
              smooth=Smooth.None)}),
      Icon(graphics={Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid), Line(
            points={{-100,0},{100,0}},
            color={0,0,0},
            smooth=Smooth.None)}));

  end MaterialTransport;

  connector MechanicalTransport "Connector for transverse displacement"

    Q.Pressure tau(nominal=1*U.Pa) "Shear stress";
    flow Q.VolumeRate APdot(nominal=1*U.cm^3/U.s)
      "Shear velocity times surface area";
    annotation (
      defaultComponentName="mechanical",
      Documentation(info="<html>

  <p>Note that the geometric orientations of shear stress and velocity
  are referenced locally.  As defined by
  <a href=\"modelica://FCSys.BaseClasses.Utilities.inSign\">inSign</a>(),
  force and velocity in the positive direction on the negative side or face
  of a region or subregion is in the globally-referenced positive direction.  The
  positive direction on the positive side is globally negative.
  </p>
    <p>For more information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(graphics={Text(
              extent={{-100,36},{100,76}},
              textString="%name",
              lineColor={0,0,0}),Rectangle(
              extent={{-30,30},{30,-30}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),Line(
              points={{0,30},{0,-30}},
              color={0,0,0},
              smooth=Smooth.None)}),
      Icon(graphics={Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid), Line(
            points={{0,100},{0,-100}},
            color={0,0,0},
            smooth=Smooth.None)}));

  end MechanicalTransport;

  connector Thermal "Connector to transport heat"

    extends BaseClasses.Thermal;

    annotation (
      Documentation(info="<html>For information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(graphics={Text(
              extent={{-100,36},{100,76}},
              textString="%name",
              lineColor={0,0,0}),Rectangle(
              extent={{-30,30},{30,-30}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),Line(
              points={{-30,30},{30,-30}},
              color={0,0,0},
              smooth=Smooth.None),Line(
              points={{30,30},{-30,-30}},
              color={0,0,0},
              smooth=Smooth.None)}),
      Icon(graphics={
          Rectangle(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Line(
            points={{-100,100},{100,-100}},
            color={0,0,0},
            smooth=Smooth.None),
          Line(
            points={{-100,-100},{100,100}},
            color={0,0,0},
            smooth=Smooth.None)}));

  end Thermal;

  connector RealInput = input Real
    "<html>\"<code>input Real</code>\" as a connector</html>" annotation (
    defaultComponentName="u",
    Icon(graphics={Polygon(
          points={{-100,100},{100,0},{-100,-100},{-100,100}},
          lineColor={0,0,127},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid)}, coordinateSystem(
        extent={{-100,-100},{100,100}},
        preserveAspectRatio=true,
        initialScale=0.1,
        grid={2,2})),
    Diagram(coordinateSystem(
        preserveAspectRatio=true,
        initialScale=0.1,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={Polygon(
          points={{0,50},{100,0},{0,-50},{0,50}},
          lineColor={0,0,127},
          fillColor={0,0,127},
          fillPattern=FillPattern.Solid), Text(
          extent={{-200,50},{200,90}},
          textString="%name",
          lineColor={0,0,0})}),
    Documentation(info="<html>
<p>
Connector with one input signal of type <code>Real</code>.</p>
</html>"));
  connector RealInputInternal = input Real
    "<html>Internal \"<code>input Real</code>\" as a connector</html>"
    annotation (
    defaultComponentPrefixes="protected",
    defaultComponentName="u",
    Icon(graphics={Polygon(
          points={{-100,100},{100,0},{-100,-100},{-100,100}},
          lineColor={0,0,127},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid)}, coordinateSystem(
        extent={{-100,-100},{100,100}},
        preserveAspectRatio=true,
        initialScale=0.1,
        grid={2,2})),
    Diagram(coordinateSystem(
        preserveAspectRatio=true,
        initialScale=0.1,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={Polygon(
            points={{0,20},{40,0},{0,-20},{0,20}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid),Text(
            extent={{-200,24},{200,64}},
            textString="%name",
            lineColor={0,0,0})}),
    Documentation(info="<html>
<p>
Protected connector with one input signal of type <code>Real</code>.</p>
</html>"));
  expandable connector RealInputBus
    "<html>Bus of <a href=\"modelica://FCSys.Connectors.RealInput\">RealInput</a> connectors</html>"

    annotation (
      defaultComponentName="u",
      Documentation(info="<html><p>There is no minimal set of variables.
    Signals are included by connecting instances
   of the <a href=\"modelica://FCSys.Connectors.RealInput\">RealInput</a> connector.<p></html>"),

      Icon(graphics={Polygon(
            points={{-100,100},{100,0},{-100,-100},{-100,100}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5)}, coordinateSystem(
          extent={{-100,-100},{100,100}},
          preserveAspectRatio=true,
          initialScale=0.1,
          grid={2,2})),
      Diagram(coordinateSystem(
          preserveAspectRatio=true,
          initialScale=0.1,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics={Polygon(
            points={{0,50},{100,0},{0,-50},{0,50}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5), Text(
            extent={{-200,50},{200,90}},
            textString="%name",
            lineColor={0,0,0})}));

  end RealInputBus;

  expandable connector RealInputBusInternal
    "<html>Internal bus of <a href=\"modelica://FCSys.Connectors.RealInput\">RealInput</a> connectors</html>"

    annotation (
      defaultComponentPrefixes="protected",
      defaultComponentName="u",
      Documentation(info="<html><p>There is no minimal set of variables.
    Signals are included by connecting instances
   of the <a href=\"modelica://FCSys.Connectors.RealInput\">RealInput</a> connector.<p></html>"),

      Icon(graphics={Polygon(
            points={{-100,100},{100,0},{-100,-100},{-100,100}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5)}, coordinateSystem(
          extent={{-100,-100},{100,100}},
          preserveAspectRatio=true,
          initialScale=0.1,
          grid={2,2})),
      Diagram(coordinateSystem(
          preserveAspectRatio=true,
          initialScale=0.1,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics={Text(
              extent={{-200,24},{200,64}},
              textString="%name",
              lineColor={0,0,0}),Polygon(
              points={{0,20},{40,0},{0,-20},{0,20}},
              lineColor={0,0,127},
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5)}));

  end RealInputBusInternal;

  connector RealOutput = output Real
    "<html>\"<code>output Real</code>\" as a connector</html>" annotation (
    defaultComponentName="y",
    Icon(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={Polygon(
          points={{-100,100},{100,0},{-100,-100},{-100,100}},
          lineColor={0,0,127},
          fillColor={0,0,127},
          fillPattern=FillPattern.Solid)}),
    Diagram(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={Polygon(
            points={{-100,50},{0,0},{-100,-50},{-100,50}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),Text(
            extent={{-200,50},{200,90}},
            textString="%name",
            lineColor={0,0,0})}),
    Documentation(info="<html>
<p>
Connector with one output signal of type <code>Real</code>.</p>
</html>"));
  connector RealOutputInternal = output Real
    "<html>Internal \"<code>output Real</code>\" as a connector</html>"
    annotation (
    defaultComponentPrefixes="protected",
    defaultComponentName="y",
    Icon(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={Polygon(
          points={{-100,100},{100,0},{-100,-100},{-100,100}},
          lineColor={0,0,127},
          fillColor={0,0,127},
          fillPattern=FillPattern.Solid)}),
    Diagram(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={Polygon(
            points={{-40,20},{0,0},{-40,-20},{-40,20}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),Text(
            extent={{-200,24},{200,64}},
            textString="%name",
            lineColor={0,0,0})}),
    Documentation(info="<html>
<p>
Protected connector with one output signal of type <code>Real</code>.</p>
</html>"));
  expandable connector RealOutputBus
    "<html>Bus of <a href=\"modelica://FCSys.Connectors.RealOutput\">RealOutput</a> connectors</html>"

    annotation (
      defaultComponentName="y",
      Documentation(info="<html><p>There is no minimal set of variables.
    Signals are included by connecting instances
   of the <a href=\"modelica://FCSys.Connectors.RealOutput\">RealOutput</a> connector.<p></html>"),

      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics={Polygon(
            points={{-100,100},{100,0},{-100,-100},{-100,100}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5)}),
      Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics={Polygon(
              points={{-100,50},{0,0},{-100,-50},{-100,50}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),Text(
              extent={{-200,50},{200,90}},
              textString="%name",
              lineColor={0,0,0})}));

  end RealOutputBus;

  expandable connector RealOutputBusInternal
    "<html>Internal bus of <a href=\"modelica://FCSys.Connectors.RealOutput\">RealOutput</a> connectors</html>"

    annotation (
      defaultComponentPrefixes="protected",
      defaultComponentName="y",
      Documentation(info="<html><p>There is no minimal set of variables.
    Signals are included by connecting instances
   of the <a href=\"modelica://FCSys.Connectors.RealOutput\">RealOutput</a> connector.<p></html>"),

      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics={Polygon(
            points={{-100,100},{100,0},{-100,-100},{-100,100}},
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5)}),
      Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics={Polygon(
              points={{-40,20},{0,0},{-40,-20},{-40,20}},
              lineColor={0,0,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),Text(
              extent={{-200,24},{200,64}},
              textString="%name",
              lineColor={0,0,0})}));

  end RealOutputBusInternal;

  package BaseClasses "Base classes (not for direct use)"
    extends Modelica.Icons.BasesPackage;

    connector Thermal "Connector to transport heat"

      Q.TemperatureAbsolute T(nominal=298.15*U.K) "Temperature";
      flow Q.Power Qdot(nominal=1*U.W) "Heat flow rate";

      annotation (Documentation(info="<html>For information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"));
    end Thermal;

    connector MechanicalExchange "Connector for linear momentum"

      parameter Integer n_lin(
        final min=0,
        final max=3) = 1
        "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>"
        annotation (HideResult=true);

      Q.Velocity phi[n_lin](each nominal=1*U.cm/U.s) "Velocity";
      flow Q.VolumeRate mPhidot[n_lin](each nominal=1*U.N) "Force";
      annotation (defaultComponentName="mechanical",Documentation(info="<html>
  <p>Note that the geometric orientation is globally referenced.  Thus,
  force is the rate of globally-referenced linear momentum into the component.
  </p>

    <p>For more information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"));
    end MechanicalExchange;

    partial connector PartialChemical
      "Partial connector to exchange material while advecting linear momentum and enthalpy"

      // Material
      Q.Number muPerT(nominal=1*U.V)
        "Quotient of electrochemical potential and temperature";
      flow Q.Current Ndot(nominal=1*U.A) "Current";

      // For linear momentum
      extends MechanicalExchange;

      // For enthalpy
      Q.Velocity2 hbar(nominal=1*U.V) "Massic enthalpy";
      flow Q.Power Hdot(nominal=10*U.W) "Enthalpy flow rate";
      annotation (
        Documentation(info="<html><p>For more information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

        Icon(graphics={Ellipse(
              extent={{-100,100},{100,-100}},
              lineColor={208,104,0},
              fillPattern=FillPattern.Solid,
              fillColor={255,128,0})}),
        Diagram(graphics={Text(
              extent={{-100,36},{100,76}},
              textString="%name",
              lineColor={0,0,0}), Ellipse(
              extent={{-30,30},{30,-30}},
              lineColor={208,104,0},
              fillPattern=FillPattern.Solid,
              fillColor={255,128,0})}));

    end PartialChemical;

    connector PartialFace
      "Partial connector to transport material, displacement, and heat of a single species"

      parameter ThermoOpt thermoOpt=ThermoOpt.OpenDiabatic
        "Thermodynamic options" annotation (HideResult=true);

      FCSys.Connectors.MaterialTransport material if thermoOpt == ThermoOpt.OpenDiabatic
        "Subconnector for material transport";
      FCSys.Connectors.Thermal thermal if not thermoOpt == ThermoOpt.ClosedAdiabatic
        "Subconnector for thermal transport";

      annotation (
        Documentation(info="<html>
    <p>The face may be closed and adiabatic, closed and diabatic, or open and diabatic (using
    <code>thermoOpt</code>).  The
    closed and diabatic combination is not allowed because it prevents the energy balance
    from being determined correctly.</p>

    <p>For more information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package and the
    <a href=\"modelica://FCSys.Connectors.Material\">Material</a> and
    <a href=\"modelica://FCSys.Connectors.Thermal\">Thermal</a> connectors.</p></html>"),

        Diagram(coordinateSystem(preserveAspectRatio=true,extent={{-100,-100},{
                100,100}}), graphics={Text(
              extent={{-100,36},{100,76}},
              textString="%name",
              lineColor={0,0,0}), Ellipse(
              extent={{-30,30},{30,-30}},
              lineColor={127,127,127},
              fillColor={191,191,191},
              fillPattern=FillPattern.Solid)}),
        Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                100,100}}), graphics={Ellipse(
              extent={{-100,100},{100,-100}},
              lineColor={127,127,127},
              fillColor={191,191,191},
              fillPattern=FillPattern.Solid)}));

    end PartialFace;

    type ThermoOpt = enumeration(
        ClosedAdiabatic "Closed and adiabatic",
        ClosedDiabatic "Closed and diabatic",
        OpenDiabatic "Open and diabatic") "Thermodynamic options for a face";
  end BaseClasses;

  annotation (Documentation(info="<html>
  <p>Three types of physical connectors are used in <a href=\"modelica://FCSys\">FCSys</a>.
  The chemical connectors (<a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a>,
  <a href=\"modelica://FCSys.Connectors.ChemicalBusInternal\">ChemicalBusInternal</a>,
  <a href=\"modelica://FCSys.Connectors.ChemicalInput\">ChemicalInput</a>, and
  <a href=\"modelica://FCSys.Connectors.ChemicalOutput\">ChemicalOutput</a>)
  represent advective exchange among species that react chemically
  within a subregion.
  The inert connectors
  (<a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a>,
  <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a>,
  <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a>, and
  <a href=\"modelica://FCSys.Connectors.InertInternal\">InertInternal</a>)
  describe diffusive exchange among species within a subregion.
  The face connectors (<a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a>,
  <a href=\"modelica://FCSys.Connectors.FaceBusInternal\">FaceBusInternal</a>,
  <a href=\"modelica://FCSys.Connectors.Face\">Face</a>,
  <a href=\"modelica://FCSys.Connectors.FaceX\">FaceX</a>,
  <a href=\"modelica://FCSys.Connectors.FaceY\">FaceY</a>, and
  <a href=\"modelica://FCSys.Connectors.FaceZ\">FaceZ</a>)
  describe advective and diffusive transport between instances of a single
  species in
  neighboring regions or subregions.</p>

  <p>The effort/flow pairs of the connectors are listed in <a href=\"#Tab1\">Table 1</a>.
  The dimensions of the efforts and flows are noted in terms of mass (M),
  amount or particle number (N),
  length (L), and time (T).  The pairs are used to describe the
  transfer of material, linear momentum, and energy associated with material,
  mechanical, and thermal interactions.  Both advection (e.g., dynamic force and
  thermal convection) and diffusion (e.g., friction and thermal conduction)
  are resolved (see the
  <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model).</p>

    <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"bottom\" id=\"Tab1\"><b>Table 1:</b> Summary of the effort/flow pairs.</caption>
      <tr>
        <th>Within icon(s)</th>
        <th>Name or quantity</th>
        <th>Effort</th>
        <th>Flow</th>
      </tr>
      <tr>
        <td>
          <a href=\"modelica://FCSys.Connectors.ChemicalInput\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalInputI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.ChemicalOutput\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalOutputI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.ChemicalBus\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalBusI.png\"></a></td>
          <!--<a href=\"modelica://FCSys.Connectors.ChemicalBusInternal\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalBusInternalI.png\"></a></td>-->
        <td valign=middle>Material exchange</td>
        <td valign=middle>Electrochemical potential divided by temperature<br>&mu;/<i>T</i> [1]</td>
        <td valign=middle>Current<br><i>N&#775;</i> [N T<sup>-1</sup>]</td>
      </tr>
      <tr>
        <td>
          <a href=\"modelica://FCSys.Connectors.ChemicalInput\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalInputI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.ChemicalOutput\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalOutputI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.ChemicalBus\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalBusI.png\"></a>
          <!--<a href=\"modelica://FCSys.Connectors.ChemicalBusInternal\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalBusInternalI.png\"></a>-->
<br>
          <a href=\"modelica://FCSys.Connectors.Inert\"><img src=\"modelica://FCSys/help/FCSys.Connectors.InertI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.InertAmagat\"><img src=\"modelica://FCSys/help/FCSys.Connectors.InertAmagatI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.InertDalton\"><img src=\"modelica://FCSys/help/FCSys.Connectors.InertDaltonI.png\"></a></td>
        <td valign=middle>Mechanical exchange</td>
        <td valign=middle>Velocity<br>&phi; [L T<sup>-1</sup>]</td>
        <td valign=middle>Force<br><i>m</i>&Phi;dot [L M T<sup>-2</sup>]</td>
      </tr>
      <tr>
        <td>
          <a href=\"modelica://FCSys.Connectors.ChemicalInput\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalInputI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.ChemicalOutput\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalOutputI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.ChemicalBus\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalBusI.png\"></a></td>
          <!--<a href=\"modelica://FCSys.Connectors.ChemicalBusInternal\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalBusInternalI.png\"></a></td>-->
        <td valign=middle>Fluid</td>
        <td valign=middle>Massic enthalpy<br><i>h&#772;</i> [L<sup>2</sup> M N<sup>-1</sup> T<sup>-2</sup>]</td>
        <td valign=middle>Enthalpy flow rate<br><i>H&#775;</i> [L<sup>2</sup>  M T<sup>-3</sup>]</td>
      </tr>
      <tr>
        <td>
          <a href=\"modelica://FCSys.Connectors.Material\"><img src=\"modelica://FCSys/help/FCSys.Connectors.MaterialI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.Face\"><img src=\"modelica://FCSys/help/FCSys.Connectors.FaceI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.FaceBus\"><img src=\"modelica://FCSys/help/FCSys.Connectors.FaceBusI.png\"></a></td>
        <td valign=middle>Material transport</td>
        <td valign=middle>Pressure<br><i>p</i> [M L<sup>-1</sup> T<sup>-2</sup>]</td>
        <td valign=middle>Current<br><i>N&#775;</i> [N T<sup>-1</sup>]</td>
      </tr>
      <tr>
        <td>
          <a href=\"modelica://FCSys.Connectors.Mechanical\"><img src=\"modelica://FCSys/help/FCSys.Connectors.MechanicalI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.Face\"><img src=\"modelica://FCSys/help/FCSys.Connectors.FaceI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.FaceBus\"><img src=\"modelica://FCSys/help/FCSys.Connectors.FaceBusI.png\"></a></td>
        <td valign=middle>Mechanical transport</td>
        <td valign=middle>Shear stress<br>&tau; [M L<sup>-1</sup> T<sup>-2</sup>]</td>
        <td valign=middle>Surface area times shear velocity<br><i>AP&#775;</i> [L<sup>3</sup> T<sup>-1</sup>]</td>
      </tr>
      <tr>
        <td>
          <a href=\"modelica://FCSys.Connectors.Thermal\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ThermalI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.Face\"><img src=\"modelica://FCSys/help/FCSys.Connectors.FaceI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.FaceBus\"><img src=\"modelica://FCSys/help/FCSys.Connectors.FaceBusI.png\"></a>
<br>
          <a href=\"modelica://FCSys.Connectors.Inert\"><img src=\"modelica://FCSys/help/FCSys.Connectors.InertI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.InertAmagat\"><img src=\"modelica://FCSys/help/FCSys.Connectors.InertAmagatI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.InertDalton\"><img src=\"modelica://FCSys/help/FCSys.Connectors.InertDaltonI.png\"></a></td>
        <td valign=middle>Thermal</td>
        <td valign=middle>Temperature<br><i>T</i> [L<sup>2</sup> M N<sup>-1</sup> T<sup>-2</sup>]</td>
        <td valign=middle>Heat flow rate<sup></sup><br><i>Q&#775;</i> [L<sup>2</sup> M T<sup>-3</sup>]</td>
      </tr>
      <tr>
        <td>
          <a href=\"modelica://FCSys.Connectors.InertAmagat\"><img src=\"modelica://FCSys/help/FCSys.Connectors.InertAmagatI.png\"></a></td>
        <td valign=middle>Additivity of volume<br>(Amagat's law)</td>
        <td valign=middle>Pressure<br><i>p</i> [M L<sup>-1</sup> T<sup>-2</sup>]</td>
        <td valign=middle>Volume<br><i>V</i> [L<sup>3</sup>]</td>
      </tr>
      <tr>
        <td>
        <a href=\"modelica://FCSys.Connectors.InertDalton\"><img src=\"modelica://FCSys/help/FCSys.Connectors.InertDaltonI.png\"></a></td>
        <td valign=middle>Additivity of pressure<br>(Dalton's law)</td>
        <td valign=middle>Volume<br><i>V</i> [L<sup>3</sup>]</td>
        <td valign=middle>Pressure<br><i>p</i> [M L<sup>-1</sup> T<sup>-2</sup>]</td>
      </tr>
    </table>

  <p>In addition to the material exchange pair, the chemical
  connectors have mechanical exchange and fluid pairs that are used to describe the purely advective
  (non-diffusive) flow associated with the material exchange (see the
  <a href=\"modelica://FCSys.Subregions.Reaction\">Reaction</a> model for details).  Although
  the physical variables are acausal, the
  <a href=\"modelica://FCSys.Connectors.ChemicalInput\">ChemicalInput</a> and
  <a href=\"modelica://FCSys.Connectors.ChemicalOutput\">ChemicalOutput</a> connectors have inputs and
  outputs to pass characteristic data of the species&mdash;the chemical formula and specific mass.
  This information is used to determine the appropriate
  stoichiometric and advective equations in the
  <a href=\"modelica://FCSys.Subregions.Reaction\">Reaction</a>
  model.</p>

  <p>There are two specialized types of inert
  connectors.  The <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector
  (with an \"A\" in the icon)
  imposes additivity of volume and is used to combine material phases within a subregion.
  The
  <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector (with a \"D\" in the icon)
  applies additivity of pressure to mix species within a material phase (e.g., N<sub>2</sub> and O<sub>2</sub>
  within a gas).
  The two cannot be directly connected because the effort/flow designations
  are opposite.  An adapter must be used
  (e.g., <a href=\"modelica://FCSys.Subregions.PhaseBoundary\">FCSys.Subregions.PhaseBoundary</a>).</p>

  <p>The face connectors contain <a href=\"modelica://FCSys.Connectors.Material\">Material</a>,
  <a href=\"modelica://FCSys.Connectors.Mechanical\">Mechanical</a>, and
  <a href=\"modelica://FCSys.Connectors.Thermal\">Thermal</a> subconnectors.
  <a href=\"modelica://FCSys.Connectors.Face\">Face</a>,
  <a href=\"modelica://FCSys.Connectors.FaceX\">FaceX</a>,
  <a href=\"modelica://FCSys.Connectors.FaceY\">FaceY</a>, and
  <a href=\"modelica://FCSys.Connectors.FaceZ\">FaceZ</a>
  differ only in the name of the mechanical subconnectors, which describe transverse displacement.
  For example, <a href=\"modelica://FCSys.Connectors.FaceX\">FaceX</a> has
  <a href=\"modelica://FCSys.Connectors.Material\">Mechanical</a> subconnectors named
  <code>mechanicalY</code> and <code>mechanicalZ</code>.</p>

  <p>The effort-flow products of the mechanical exchange and transport
  pairs are the rates of energy associated with the interactions.
  However, many of the pairs are different.
  The effort of the material exchange pair is electrochemical potential divided by
  temperature because it is used in the chemical equilibrium
  (see the
  <a href=\"modelica://FCSys.Subregions.Reaction\">Reaction</a> model).  The
  fluid pair's effort is massic enthalpy (enthalpy divided by mass) since
  it is used for advection only.  The effort of the material transport pair
  is pressure rather than electrochemical potential so that the dynamic normal
  force can be resolved without nonlinear systems of equations (see the
  <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model).
  The flow variable of the thermal pair is heat flow rate instead of entropy flow
  rate so that the transport equations are linear and follow the traditional
  form (e.g., Fourier's law; see
  <a href=\"modelica://Modelica.Thermal.HeatTransfer\">Modelica.Thermal.HeatTransfer</a>).
  The additivity of volume and
  pressure pairs describe static balances; therefore, their flows are the quantities
  rather than the rates of the quantities.
  </p>

  <p><a href=\"#Fig1\">Figure 1</a> shows the instantiation hierarchy of the connectors.
  <a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a>
  is an expandable connector that groups
  <a href=\"modelica://FCSys.Connectors.ChemicalInput\">ChemicalInput</a> and
  <a href=\"modelica://FCSys.Connectors.ChemicalInput\">ChemicalOutput</a> connectors
  of multiple species.  Likewise,
  <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> groups
  <a href=\"modelica://FCSys.Connectors.Face\">Face</a>,
  <a href=\"modelica://FCSys.Connectors.FaceX\">FaceX</a>,
  <a href=\"modelica://FCSys.Connectors.FaceY\">FaceY</a>, or
  <a href=\"modelica://FCSys.Connectors.FaceZ\">FaceZ</a>
  connectors of multiple species.
  The
  <a href=\"modelica://FCSys.Connectors.ChemicalBusInternal\">ChemicalBusInternal</a> and
  <a href=\"modelica://FCSys.Connectors.FaceBusInternal\">FaceBusInternal</a> connectors (not shown)
  are versions of
  <a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a> and
  <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a>
  that merely have a default <code>protected</code> prefix and smaller icons to
  indicate that the connectors are internal to a model.  The inert connectors are
  flat; they do not have subconnectors.</p>

  <p align=center id=\"Fig1\"><img src=\"modelica://FCSys/resources/documentation/Connectors/hierarchy.png\">
<br><b>Figure 1:</b> Instantiation hierarchy of the connectors.</p>

    <p><b>Relation to thermodynamics and justification of the Amagat and Dalton connectors:</b><p>

    <p>In order to describe the dynamic behavior of a physical system, a model must include conservation laws or rate balances.  These
    equations involve the storage and flow of extensive quantities within (among species) and into the system.  In chemical/thermal systems, the extensive
    quantities of interest are particle number (or mass) and entropy or energy.
    For the sake of simplicity, momentum will be excluded from the present discussion; assume that the fluid is macroscopically stagnant.
    Also assume that there is only one inlet or outlet to the system.
    In terms of mathematics, we have introduced four variables (2 flows and 2 quantities) but only two equations (material and energy
    conservation).
    </p>

    <p>Two additional equations involve flow rates; these are transport equations with spatial nature&mdash;separate from the temporal
    conservation equations.  An
    extensive body of empirical evidence indicates that the
    the flows are related to differences in efforts
    or generalized \"driving forces.\"  The most appropriate efforts are conjugate to the quantities with respect to energy or entropy.
    For the chemical/thermal system, the efforts are then electrochemical
    potential and temperature.
    Yet these are intensive
    properties&mdash;distinct from the quantities, which are extensive.
    So far, there are two rate balances to relate extensive quantities to flows and two transport equations
    (4 equations in all) and six variables (2 quantities, 2 flows, and 2 efforts or intensive properties).
    </p>

    <p>One extensive quantity can be divided by the other to yield an intensive property.
    For example, internal energy can be divided by particle number
    to give internal potential (the relationship is not as direct for electrochemical potential, but the concept still holds).
    The other equation involves the spatial extent of the system, for example, the extensive
    volume of the system divided by particle number to give specific volume.
    This introduces another variable (extensive volume);
    now there are six equations and seven variables.</p>

    <p>Fortunately, we may assume that the extensive volume of the system is fixed (i.e., that the system
    is a \"control volume\").  If there is only one species in the system, we assume that it fills the entire volume
    (e.g., no macroscopically observable regions of vacuum).  If another species is included in the system, the number of variables is
    doubled.  All of the equations may be repeated except that the specific volume of each species is its own extensive volume or
    \"partial volume\" divided by its particle number to give \"partial specific volume.\"
    It is reasonable to assume that the sum of the partial volumes is equal to the total volume of the system
    (again, no voids).
    This is a generalization
    of the previous equation that set the volume of the single species equal to the volume of the system or control volume.  However, now there are three volumes
    (of each species and of the system) instead of two (of the one species and of the system) but no additional equations.</p>

    <p>In general,
    an additional equation may be added to exchange volume between the two species such that they reach equilibrium.
    This could be modeled by another transport-like equation, with a slight modification.
    However, in the
    <a href=\"modelica://FCSys\">FCSys</a> package, it is assumed that this equilibrium already/always exists.
    Since we wish to impose that the sum of the two
    partial volumes is equal to the total volume, it is appropriate to set the flow variable to be the
    quantity itself (volume) rather than
    the rate of the quantity.  Then, there is no need for another rate balance to relate the quantity to the flow; the quantity <i>is</i>
    the flow.  In this case,
    the most
    appropriate effort variable is pressure (or ideal density).  The relationship among pressure,
    specific volume, and temperature is given by an equation of state.
    This \"additivity of volume\"
    interaction occurs through the
    <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector.
    </p>

    <p>If the species are mixed, it may be more appropriate to assume that the pressures
    of the components of a mixture sum to the total pressure of the mixture.  This \"additivity of pressure\"
    is described by connections of the
    <a href=\"modelica://FCSys.Connectors.InertDelton\">InertDalton</a> connector (instead of
    instead of <a href=\"modelica://FCSys.Connectors.InertDelton\">InertAmagat</a>).
    </p>

  <p>
  <b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
  Copyright 2007&ndash;2013, Georgia Tech Research Corporation.
  </p>
  <p>
  <i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
  it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
  disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.ModelicaLicense2\">
  FCSys.UsersGuide.ModelicaLicense2</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
  http://www.modelica.org/licenses/ModelicaLicense2</a>.</i>
  </p>
  </html>"));
end Connectors;
