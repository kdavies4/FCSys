within FCSys;
package Connectors "Declarative and imperative connectors"
  extends Modelica.Icons.InterfacesPackage;

  expandable connector ChemicalBus
    "<html>Bus of <a href=\"modelica://FCSys.Connectors.BaseClasses.Chemical\">Chemical</a> connectors (for multiple species)</html>"

    annotation (
      defaultComponentName="chemical",
      Documentation(info="<html><p>There is no minimal set of variables.  Species are included by connecting instances
    of a <a href=\"modelica://FCSys.Connectors.BaseClasses.Chemical\">Chemical</a> connector
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
            lineColor={0,0,0}), Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={208,104,0},
            fillPattern=FillPattern.Solid,
            fillColor={255,128,0},
            lineThickness=0.5)}));
  end ChemicalBus;

  expandable connector ChemicalBusInternal
    "<html>Internal bus of <a href=\"modelica://FCSys.Connectors.BaseClasses.Chemical\">Chemical</a> connectors (for multiple species)</html>"

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
            lineThickness=0.5), Text(
            extent={{-100,20},{100,60}},
            textString="%name",
            lineColor={0,0,0})}));

  end ChemicalBusInternal;

  connector ChemicalInput
    "Connector to exchange material with advection of linear momentum and energy, with chemical formula as input"

    extends FCSys.Connectors.BaseClasses.Chemical;
    input String formula(start="") "Chemical formula of the species";
    // Note:  The start value prevents a warning when checked in Dymola
    // 7.4.
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
    "Connector to exchange material with advection of linear momentum and energy, with chemical formula as output"

    extends FCSys.Connectors.BaseClasses.Chemical;
    output String formula "Chemical formula of the species";

    annotation (
      defaultComponentName="chemical",
      Documentation(info="<html><p>For more information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(graphics={Ellipse(extent={{-100,100},{100,-100}}, lineColor={208,104,
                0})}),
      Diagram(graphics));

  end ChemicalOutput;

  expandable connector FaceBus
    "<html>Bus of <a href=\"modelica://FCSys.Connectors.FaceX\">FaceX</a>, <a href=\"modelica://FCSys.Connectors.FaceY\">FaceY</a>, <a href=\"modelica://FCSys.Connectors.FaceZ\">FaceZ</a>, or <a href=\"modelica://FCSys.Connectors.FaceGeneric\">FaceGeneric</a> connectors (for multiple species)</html>"

    annotation (
      defaultComponentName="face",
      Documentation(info="<html><p>There is no minimal set of variables.  Species are included by connecting instances
    of a <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">Face</a> connector
    (<a href=\"modelica://FCSys.Connectors.FaceX\">FaceX</a>,
    <a href=\"modelica://FCSys.Connectors.FaceY\">FaceY</a>,
    <a href=\"modelica://FCSys.Connectors.FaceZ\">FaceZ</a>, or
    <a href=\"modelica://FCSys.Connectors.FaceGeneric\">FaceGeneric</a>).  In order to allow
    the subconnectors of the <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">Face</a> connectors
    (instances of <a href=\"modelica://FCSys.Connectors.Material\">Material</a>,
    <a href=\"modelica://FCSys.Connectors.MomentumLineic\">MomentumLineic</a>,
    and <a href=\"modelica://FCSys.Connectors.Entropy\">Entropy</a>) to be included independently, those
    subconnectors are connected explicitly.  For example,
    <blockquote>
        <code>
        connect(species.xNegative.material, xNegative.species.material);<br>
        connect(species.xNegative.momentumY, xNegative.species.momentumY);<br>
        connect(species.xNegative.momentumZ, xNegative.species.momentumZ);<br>
        connect(species.xNegative.entropy, xNegative.species.entropy);
        </code>
    </blockquote>
    where <code>species</code> is an instance of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a>
    model in which <code>xNegative</code> is an instance of the <a href=\"modelica://FCSys.Connectors.FaceX\">FaceX</a> connector.
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
    "<html>Internal bus of <a href=\"modelica://FCSys.Connectors.FaceX\">FaceX</a>, <a href=\"modelica://FCSys.Connectors.FaceY\">FaceY</a>, <a href=\"modelica://FCSys.Connectors.FaceZ\">FaceZ</a>, or <a href=\"modelica://FCSys.Connectors.FaceGeneric\">FaceGeneric</a> connectors (for multiple species)</html>"

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

  connector FaceX
    "X-axis connector to transport material, linear momentum, and entropy of a single species"

    extends BaseClasses.PartialFace;

    parameter Boolean viscousY=false "Viscous along the y axis" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true));
    parameter Boolean viscousZ=false "Viscous along the z axis" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true));

    FCSys.Connectors.MomentumLineic momentumY if viscousY "Y-axis momentum";
    FCSys.Connectors.MomentumLineic momentumZ if viscousZ "Z-axis momentum";

    annotation (
      Documentation(info="<html>
  <p>The options <code>viscousY</code> and <code>viscousZ</code> are slight misnomers.
  For example, even if <code>viscousY == true</code>, it is possible that the flow rate of y-axis momentum
  is zero, which is an inviscous condition.  When these options are enabled, it merely indicates that there is
  a <i>possibility</i> of nonzero shear force.</p>

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
    "Y-axis connector to transport material, linear momentum, and entropy of a single species"
    extends BaseClasses.PartialFace;

    parameter Boolean viscousZ=false "Viscous along the z axis" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true));
    parameter Boolean viscousX=false "Viscous along the x axis" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true));

    FCSys.Connectors.MomentumLineic momentumZ if viscousZ "Z-axis momentum";
    FCSys.Connectors.MomentumLineic momentumX if viscousX "X-axis momentum";

    annotation (
      Documentation(info="<html>
  <p>The options <code>viscousZ</code> and <code>viscousX</code> are slight misnomers.
  For example, even if <code>viscousZ == true</code>, it is possible that the flow rate of z-axis momentum
  is zero, which is an inviscous condition.  When these options are enabled, it merely indicates that there is
  a <i>possibility</i> of nonzero shear force.</p>

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
    "Z-axis connector to transport material, linear momentum, and entropy of a single species"
    extends BaseClasses.PartialFace;

    parameter Boolean viscousX=false "Viscous along the x axis" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true));
    parameter Boolean viscousY=false "Viscous along the y axis" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true));

    FCSys.Connectors.MomentumLineic momentumX if viscousX "X-axis momentum";
    FCSys.Connectors.MomentumLineic momentumY if viscousY "Y-axis momentum";

    annotation (
      Documentation(info="<html>
  <p>The options <code>viscousX</code> and <code>viscousY</code> are slight misnomers.
  For example, even if <code>viscousX == true</code>, it is possible that the flow rate of x-axis momentum
  is zero, which is an inviscous condition.  When these options are enabled, it merely indicates that there is
  a <i>possibility</i> of nonzero shear force.</p>

  <p>For more information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package and the
    <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">PartialFace</a> connector.</p></html>"),

      Diagram(coordinateSystem(preserveAspectRatio=true,extent={{-100,-100},{
              100,100}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Ellipse(extent={{-100,100},{100,-100}},
              lineColor={127,127,127})}));

  end FaceZ;

  connector FaceGeneric
    "Connector to transport material, linear momentum, and entropy of a single species along a generic axis"
    extends BaseClasses.PartialFace;

    parameter Boolean viscous1=false
      "<html>Viscous along the 1<sup>st</sup> transverse axis</html>"
      annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true));
    parameter Boolean viscous2=false
      "<html>Viscous along the 2<sup>nd</sup> transverse axis</html>"
      annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true));

    FCSys.Connectors.MomentumLineic momentum1 if viscous1
      "1st transverse momentum";
    FCSys.Connectors.MomentumLineic momentum2 if viscous2
      "2nd transverse momentum";

    annotation (
      defaultComponentName="face",
      Documentation(info="<html>
  <p>The options <code>viscous1</code> and <code>viscous2</code> are slight misnomers.
  For example, even if <code>viscous1 == true</code>, it is possible that the flow rate of the first transverse momentum
  is zero, which is an inviscous condition.  When these options are enabled, it merely indicates that there is
  a <i>possibility</i> of nonzero shear force.</p>

  <p>For more information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package and the
    <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">PartialFace</a> connector.</p></html>"),

      Diagram(coordinateSystem(preserveAspectRatio=true,extent={{-100,-100},{
              100,100}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Ellipse(extent={{-100,100},{100,-100}},
              lineColor={127,127,127})}));

  end FaceGeneric;

  connector InertAmagat
    "Connector to exchange linear momentum and entropy by diffusion, with additivity of volume"

    extends FCSys.Connectors.BaseClasses.PartialInert;

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
    partial volumes are additive and that the mixtures exist at the same pressure.  Within
    a phase, the species are mixed according to Dalton's law (see the <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector).</p>

    <p>In order to implement Amagat's law, this connector includes volume (not rate of volume) as a flow variable.
    The effort variable is pressure such that the effort and flow variables are conjugates of
    energy (not power).
    </p>

    <p>See also the <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialInert\">PartialInert</a> and
    <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connectors and
    the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(graphics={Text(
            extent={{-30,20},{30,-20}},
            lineColor={0,0,0},
            textString="A")}),
      Icon(graphics={Ellipse(extent={{-100,100},{100,-100}}, lineColor={72,90,
                180}), Text(
            extent={{-100,80},{100,-90}},
            lineColor={0,0,0},
            textString="A")}));

  end InertAmagat;

  connector InertDalton
    "Connector to exchange linear momentum and entropy by diffusion, with additivity of pressure"

    extends FCSys.Connectors.BaseClasses.PartialInert;

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
    <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialInert\">PartialInert</a> connectors and
    the documentation in the
  <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(graphics={Text(
            extent={{-30,20},{30,-20}},
            lineColor={0,0,0},
            textString="D")}),
      Icon(graphics={Ellipse(extent={{-100,100},{100,-100}}, lineColor={72,90,
                180}), Text(
            extent={{-100,80},{100,-90}},
            lineColor={0,0,0},
            textString="D")}));

  end InertDalton;

  connector Material "Connector for material"

    Q.Potential mu(nominal=1*U.V) "Electrochemical potential";
    flow Q.Current Ndot(nominal=1*U.A) "Current";

    annotation (
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

  end Material;

  connector MomentumLineic
    "Connector for linear momentum in a single direction"

    Q.Velocity phi(nominal=1*U.cm/U.s) "Velocity" annotation (HideResult=true);
    flow Q.Force mPhidot(nominal=1*U.N) "Force" annotation (HideResult=true);
    annotation (
      defaultComponentName="momentum",
      Documentation(info="<html>
  <p>Note that the orientation is referenced globally.  For example, force may be in the positive-x
  direction&mdash;not from the interface into component as is pressure.  Thus, force is the rate of
  globally-referenced linear momentum into the component.
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

  end MomentumLineic;

  connector Entropy "Connector for entropy"
    Q.TemperatureAbsolute T(nominal=298.15*U.K) "Temperature";
    flow Q.Current Sdot(nominal=1*U.W/(298.15*U.K)) "Entropy flow rate";
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

  end Entropy;

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
            fillPattern=FillPattern.Solid),Text(
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

    connector Chemical
      "Connector to exchange material with advection of linear momentum and energy"

      parameter Integer n_lin(
        final min=0,
        final max=3) = 1
        "<html>Number of axes of linear momentum (<i>n</i><sub>lin</sub>)</html>"
        annotation (HideResult=true);

      // Material
      Q.Number muPerT(nominal=1*U.V)
        "Electrochemical potential divided by temperature";
      flow Q.Current Ndot(nominal=1*U.A) "Current";

      // For linear momentum
      Q.VelocityMassSpecific mphi[n_lin](each nominal=1*U.g*U.cm/(U.mol*U.s))
        "Specific mass times velocity";
      flow Q.Force mPhidot[n_lin](each nominal=1*U.N) "Force due to advection";

      // For energy
      Q.Potential h(nominal=1*U.V) "Specific enthalpy";
      flow Q.Power Hdot(nominal=10*U.W) "Rate of advection of enthalpy";
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

    end Chemical;

    connector PartialFace
      "Partial connector to transport material, linear momentum, and entropy of a single species"
      parameter MaterialEntropyOpt matEntOpt=MaterialEntropyOpt.ClosedAdiabatic
        "Options for material and thermal transport"
        annotation (HideResult=true,Dialog(compact=true));

      FCSys.Connectors.Material material if matEntOpt == MaterialEntropyOpt.OpenDiabatic;
      FCSys.Connectors.Entropy entropy if matEntOpt == MaterialEntropyOpt.ClosedDiabatic
         or matEntOpt == MaterialEntropyOpt.OpenDiabatic;

      annotation (
        Documentation(info="<html>
  <p>For more information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package and the
    <a href=\"modelica://FCSys.Connectors.MomentumLineic\">MomentumLineic</a> and
    <a href=\"modelica://FCSys.Connectors.Entropy\">Entropy</a> connectors.</p></html>"),

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

    partial connector PartialInert
      "Partial connector to exchange linear momentum and entropy by diffusion"

      parameter Integer n_lin(
        final min=0,
        final max=3) = 0
        "<html>Number of axes of linear momentum (<i>n</i><sub>lin</sub>)</html>"
        annotation (HideResult=true);

      // Linear momentum
      Q.Velocity phi[n_lin](each nominal=1*U.cm/U.s) "Velocity";
      flow Q.Force mPhidot[n_lin](each nominal=1*U.N) "Force";

      // Entropy
      Q.TemperatureAbsolute T(nominal=298.15*U.K,displayUnit="K") "Temperature";
      flow Q.Current Sdot(nominal=1*U.W/(298.15*U.K)) "Entropy flow rate";

      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="inert",
        Documentation(info="<html>
  <p>Note that the geometric orientation is globally referenced.  For example, force may be in the positive-x
  direction&mdash;not from the interface into component as is pressure.  Thus, force is the rate of
  globally-referenced linear momentum into the component.
  </p>
    <p>For more information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

        Diagram(graphics={Text(
              extent={{-100,36},{100,76}},
              textString="%name",
              lineColor={0,0,0}), Ellipse(
              extent={{-30,30},{30,-30}},
              lineColor={72,90,180},
              fillPattern=FillPattern.Solid,
              fillColor={102,128,255})}),
        Icon(graphics={Ellipse(
              extent={{-100,100},{100,-100}},
              lineColor={72,90,180},
              fillPattern=FillPattern.Solid,
              fillColor={102,128,255})}));

    end PartialInert;

    type MaterialEntropyOpt = enumeration(
        ClosedAdiabatic "Closed and adiabatic",
        ClosedDiabatic "Closed and diabatic",
        OpenDiabatic "Open and diabatic")
      "Choices for material and thermal transport through a face" annotation (
        Documentation(info="<html>
  <p>This enumeration is used instead of independent options for open/closed and diabatic/adiabatic
  because the open and adiabatic combination is invalid.  Since the driving force for material
  transport is electrochemical potential, yet the contribution of material to the energy balance
  is specific enthalpy, entropy must be advected if material is transported
  (<i>h</i> <i>N&#775;</i> = (g + <i>T</i> <i>s</i>) <i>N&#775;</i> ).
  </p>
  </html>"));
  end BaseClasses;

  annotation (Documentation(info="<html>
  <p>The connectors of <a href=\"modelica://FCSys\">FCSys</a> are nested according to the hierarchy shown in Figure 1.
  The icons on the bottom row represent the fundamental connectors, which each have only one
  effort/flow pair.
  The fundamental connectors are base classes for
  the <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialInert\">Inert</a>,
  <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">Face</a>, and
  <a href=\"modelica://FCSys.Connectors.BaseClasses.Chemical\">Chemical</a>
  connectors. The <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialInert\">Inert</a>
  connector describes diffusive exchange among <a href=\"modelica://FCSys.Subregions.Species\">Species</a>
  instances within a
  <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a>.  The
  <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">Face</a> connector
  describes advective and diffusive transport between instances of a single
  <a href=\"modelica://FCSys.Subregions.Species\">Species</a> in
  neighboring subregions.  The <a href=\"modelica://FCSys.Connectors.BaseClasses.Chemical\">Chemical</a>
  connector represents the advective exchange among <a href=\"modelica://FCSys.Subregions.Species\">Species</a> that react chemically
  within a
  <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a>.  <a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a>
  and <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> are expandable connectors
  that group
  the <a href=\"modelica://FCSys.Connectors.BaseClasses.Chemical\">Chemical</a> and
  <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">Face</a> connectors
  of multiple species.  The
  <a href=\"modelica://FCSys.Connectors.ChemicalBusInternal\">ChemicalBusInternal</a> and
  <a href=\"modelica://FCSys.Connectors.FaceBusInternal\">FaceBusInternal</a> connectors (not shown)
  are versions of the
  <a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a> and
  <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a>
  connectors that merely have a default <code>protected</code> prefix and smaller icons to
  indicate that the connectors are internal to a model.</p>

    <p align=center><img src=\"modelica://FCSys/resources/images/Connectors/hierarchy.png\">
<br><b>Figure 1:</b> Hierarchy of the connectors.</p>

  <p>The contents of the fundamental connectors are summarized in Table 1.
  The dimensions are noted in terms of mass (M), length (L), time (T), and particle number (N).
  Each effort variable is chosen such that the product of effort and the rate of the quantity is an energy
  rate<a href=\"#footnote-1\"><sup>1</sup></a>.  In most connectors, the flow variable is the
  rate of the quantity, but the additivity of volume and
  pressure connectors use the quantity itself as the flow variable.
  </p>

  <p>**Update this: In addition to the material effort and flow pair, the <a href=\"modelica://FCSys.Connectors.BaseClasses.Chemical\">Chemical</a>
  connector has <code>stream</code> variables&mdash;specific mass times velocity (m&phi;)
  and specific entropy times temperature (<i>sT</i> ).  These variables describe the purely advective
  (non-diffusive) flow associated with the chemical exchange.  The rate of advection of linear momentum is the
  product of specific mass of the source, the velocity of the source, and the current
  (<i>m</i> &phi;<i>N&#775;</i>).  The rate of thermal advection is the
  product of specific entropy of the source, the temperature of the source, and the current
  (<i>sT</i><i>N&#775;</i>). Since the chemical potential is given by the Gibbs potential
  (<i>g</i> = <i>h</i> - <i>sT</i> )
  and its stoichiometrically weighted sum is zero, the rate of thermal advection from reactants
  to products (or vice versa) is also the rate of advection of enthalpy.</p>

  <p>There are two types of chemical
  connectors.  The <a href=\"modelica://FCSys.Connectors.ChemicalInput\">ChemicalInput</a> connector
  has inputs for the chemical formula (<code>formula</code> string) of the associated species.
  The <a href=\"modelica://FCSys.Connectors.ChemicalInput\">ChemicalInput</a> connector
  has outputs for the same variables.  This information is used to determine the appropriate reaction
  stoichiometry and conservation equations in the
  <a href=\"modelica://FCSys.Subregions.Reaction\">Reaction</a>
  model.</p>

  <p>There are also two types of inert
  connectors.  The <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector
  (with an \"A\" in the icon)
  uses additivity of volume and the
  <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector (with a \"D\" in the icon)
  uses
  additivity of pressure.  The two cannot be directly connected; a failure will occur during translation
  because the effort/flow designations
  are opposite.  An adapter must be used (e.g., <a href=\"modelica://FCSys.Subregions.PhaseBoundary\">FCSys.Subregions.PhaseBoundary</a>).</p>

  <p>There are multiple <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">Face</a>
  connectors&mdash;<a href=\"modelica://FCSys.Connectors.FaceX\">FaceX</a>,
  <a href=\"modelica://FCSys.Connectors.FaceY\">FaceY</a>,
  <a href=\"modelica://FCSys.Connectors.FaceZ\">FaceZ</a>, and
  <a href=\"modelica://FCSys.Connectors.FaceGeneric\">FaceGeneric</a>.
  They differ only in the way that linear momentum is labeled.  The linear momentum is oriented
  transverse to the face.
  For example, in the <a href=\"modelica://FCSys.Connectors.FaceX\">FaceX</a> connector it is labeled
  <code>momentumY</code> and <code>momentumZ</code>.
  </p>

    <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"bottom\"><b>Table 1:</b> Summary of the fundamental connectors.</caption>
      <tr>
        <th>Within icon(s)</th>
        <th>Name or quantity</th>
        <th>Effort</th>
        <th>Flow</th>
      </tr>
      <tr>
        <td>
          <a href=\"modelica://FCSys.Connectors.Material\"><img src=\"modelica://FCSys/help/FCSys.Connectors.MaterialI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\"><img src=\"modelica://FCSys/help/FCSys.Connectors.BaseClasses.PartialFaceI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.BaseClasses.Chemical\"><img src=\"modelica://FCSys/help/FCSys.Connectors.BaseClasses.ChemicalI.png\"></a></td>
        <td valign=middle>Material</td>
        <td valign=middle>Electrochemical potential<br>&mu; [L<sup>2</sup> M N<sup>-1</sup> T<sup>-2</sup>]</td>
        <td valign=middle>Current<br><i>N&#775;</i> [N T<sup>-1</sup>]</td>
      </tr>
      <tr>
        <td>
          <a href=\"modelica://FCSys.Connectors.MomentumLineic\"><img src=\"modelica://FCSys/help/FCSys.Connectors.MomentumLineicI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\"><img src=\"modelica://FCSys/help/FCSys.Connectors.BaseClasses.PartialFaceI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.InertAmagat\"><img src=\"modelica://FCSys/help/FCSys.Connectors.InertAmagatI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.InertDalton\"><img src=\"modelica://FCSys/help/FCSys.Connectors.InertDaltonI.png\"></a>
        <td valign=middle>Linear momentum</td>
        <td valign=middle>Velocity<br>&phi; [L T<sup>-1</sup>]</td>
        <td valign=middle>Force<br><i>m</i>&Phi;dot [L M T<sup>-2</sup>]</td>
      </tr>
      <tr>
        <td>
          <a href=\"modelica://FCSys.Connectors.Entropy\"><img src=\"modelica://FCSys/help/FCSys.Connectors.EntropyI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\"><img src=\"modelica://FCSys/help/FCSys.Connectors.BaseClasses.PartialFaceI.png\"></a>
        <a href=\"modelica://FCSys.Connectors.InertAmagat\"><img src=\"modelica://FCSys/help/FCSys.Connectors.InertAmagatI.png\"></a>
        <a href=\"modelica://FCSys.Connectors.InertDalton\"><img src=\"modelica://FCSys/help/FCSys.Connectors.InertDaltonI.png\"></a></td>
        <td valign=middle>Entropy</td>
        <td valign=middle>Temperature<br><i>T</i> [L<sup>2</sup> M N<sup>-1</sup> T<sup>-2</sup>]</td>
        <td valign=middle>Entropy flow rate<br><i>S&#775;</i> [N T<sup>-1</sup>]</td>
      </tr>
      <tr>
        <td><a href=\"modelica://FCSys.Connectors.BaseClasses.PartialInert\"><img src=\"modelica://FCSys/resources/images/Connectors/VolumeOrPressureI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.InertAmagat\"><img src=\"modelica://FCSys/help/FCSys.Connectors.InertAmagatI.png\"></a></td>
        <td valign=middle>Additivity of volume<br>(Amagat's law)</td>
        <td valign=middle>Pressure<br><i>p</i> [M L<sup>-1</sup> T<sup>-2</sup>]</td>
        <td valign=middle>Volume<br><i>V</i> [L<sup>3</sup>]</td>
      </tr>
      <tr>
        <td><a href=\"modelica://FCSys.Connectors.BaseClasses.PartialInert\"><img src=\"modelica://FCSys/resources/images/Connectors/VolumeOrPressureI.png\"></a>
        <a href=\"modelica://FCSys.Connectors.InertDalton\"><img src=\"modelica://FCSys/help/FCSys.Connectors.InertDaltonI.png\"></a></td>
        <td valign=middle>Additivity of pressure<br>(Dalton's law)</td>
        <td valign=middle>Volume<br><i>V</i> [L<sup>3</sup>]</td>
        <td valign=middle>Pressure<br><i>p</i> [M L<sup>-1</sup> T<sup>-2</sup>]</td>
      </tr>
    </table>

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

    <p>One extensive
    quantity can be divided by the other to yield an intensive property.  For example, internal energy can be divided by particle number
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

<p id=\"footnote-1\"><sup>1</sup> The <a href=\"modelica://Modelica.Thermal.HeatTransfer\">Modelica.Thermal.HeatTransfer</a>
package uses heat rate as the flow variable rather than entropy rate.  The reason is that it is awkward to
apply Fourier's Law (thermal conduction) in terms of entropy flow rate to a component that does not store heat,
since entropy is generated [<a href=\"modelica://FCSys.UsersGuide.References\">Hogan2006</a>],
[<a href=\"modelica://FCSys.UsersGuide.References\">Cellier1991</a>, p. 304].  However, the
<a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">Species</a> model is transient and does not directly
impose that the heat transfer rates sum to zero.  Therefore, it is just as convenient to write Fourier's law in terms of
entropy flow rate.  This approach is chosen in <a href=\"modelica://FCSys\">FCSys</a> so that all transient effort/flow variables
are conjugates of an energy rate.</p>

<p>
<b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
Copyright 2007&ndash;2012, Georgia Tech Research Corporation.
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
