within FCSys;
package Connectors "Declarative and imperative connectors"
  extends Modelica.Icons.InterfacesPackage;

  expandable connector ChemicalIOBus
    "<html>Bus of <a href=\"modelica://FCSys.Connectors.ChemicalInput\">ChemicalInput</a> and <a href=\"modelica://FCSys.Connectors.ChemicalOutput\">ChemicalOutput</a> connectors (for multiple species)</html>"

    annotation (
      defaultComponentName="chemical",
      Documentation(info="<html><p>There is no minimal set of variables.  Species are included by connecting instances
    of a <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialChemical\">Chemical</a> connector
    (<a href=\"modelica://FCSys.Connectors.ChemicalInput\">ChemicalInput</a> or
    <a href=\"modelica://FCSys.Connectors.ChemicalOutput\">ChemicalOutput</a>).</p></html>"),

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

  end ChemicalIOBus;

  connector ChemicalInput
    "**Connector to exchange material while advecting linear momentum and enthalpy, with characteristic data as input"

    input String formula "Chemical formula";
    input Integer nu "Stoichiometric coefficient";
    input Integer n_spec(min=0) "Number of species in the reaction";
    input Q.CurrentAbsolute Io "Exchange current";
    // Note:  The start value prevents a warning when checked in Dymola 7.4.
    annotation (
      defaultComponentName="chemical",
      Documentation(info="<html><p>See the information in the
    <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialChemical\">PartialChemical</a>
    connector and the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={208,104,0},
            fillColor={255,128,0},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.Dash)}),
      Diagram(graphics={Text(
              extent={{-100,36},{100,76}},
              textString="%name",
              lineColor={0,0,0}),Ellipse(
              extent={{-30,30},{30,-30}},
              lineColor={208,104,0},
              fillPattern=FillPattern.Solid,
              fillColor={255,128,0},
              pattern=LinePattern.Dash)}));

  end ChemicalInput;

  connector ChemicalOutput
    "**Connector to exchange material while advecting linear momentum and enthalpy, with characteristic data as input"

    output String formula "Chemical formula";
    output Integer nu "Stoichiometric coefficient";
    output Integer n_spec(min=0) "Number of species in the reaction";
    output Q.CurrentAbsolute Io "Exchange current";
    // Note:  The start value prevents a warning when checked in Dymola 7.4.
    annotation (
      defaultComponentName="chemical",
      Documentation(info="<html><p>See the information in the
    <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialChemical\">PartialChemical</a>
    connector and the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={208,104,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.Dash)}),
      Diagram(graphics={Text(
            extent={{-100,36},{100,76}},
            textString="%name",
            lineColor={0,0,0}), Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={208,104,0},
            fillPattern=FillPattern.Solid,
            fillColor={255,255,255},
            pattern=LinePattern.Dash)}));

  end ChemicalOutput;

  partial connector Chemical
    "Partial connector to exchange material while advecting linear momentum and enthalpy"

    // Material
    Q.Number sprime(nominal=10) "Modified specific entropy";
    flow Q.Number sprime_eq(nominal=10)
      "Modified specific entropy at equilibrium";

    // Mechanical
    extends FCSys.Connectors.Mechanical;

    // Fluid
    extends FCSys.Connectors.Thermal;
    annotation (
      Documentation(info="<html><p>For more information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={208,104,0},
            fillPattern=FillPattern.Solid,
            fillColor={255,128,0})}),
      Diagram(graphics={Ellipse(
              extent={{-30,30},{30,-30}},
              lineColor={208,104,0},
              fillPattern=FillPattern.Solid,
              fillColor={255,128,0})}));

  end Chemical;

  expandable connector FaceBus
    "<html>Bus of <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connectors (for multiple species)</html>"

    annotation (
      defaultComponentName="face",
      Documentation(info="<html><p>There is no minimal set of variables.  Species are included by connecting instances
    of <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connectors.</p></html>"),

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
    "<html>Internal bus of <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connectors (for multiple species)</html>"

    annotation (
      defaultComponentPrefixes="protected",
      defaultComponentName="face",
      Documentation(info="<html><p>This is copy of the <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, except that it
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
    "Connector to transport linear momentum and heat of a single species"

    // Material
    Q.AmountVolumic rho(nominal=298.15*U.K/U.atm) "Density";
    flow Q.Current Ndot(nominal=U.A) "Current";

    // Normal
    extends Mechanical(final n_lin=1);

    // Transverse
    Q.Pressure tau[2] "Shear stress";
    Q.VolumeRate PAdot[2] "**Perimeter rate times surface area";

    // Thermal
    extends Thermal;
    annotation (
      Documentation(info="<html>
    <p>Note that the geometric orientation is global.
    Force is the rate of globally-referenced linear momentum into the component.
    Areic current is positive in the globally positive direction (not inward).
    This is different than the
    <a href=\"modelica://FCSys.Connectors.ChemicalInput\">ChemicalInput</a> and
    <a href=\"modelica://FCSys.Connectors.ChemicalOutput\">ChemicalOutput</a> connectors,
    where current is positive into the component.</p>

    <p>For more information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors.BaseClasses.Mechanical\">Mechanical</a> base class and the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(coordinateSystem(preserveAspectRatio=true,extent={{-100,-100},{
              100,100}}), graphics={Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={127,127,127},
            fillColor={191,191,191},
            fillPattern=FillPattern.Solid)}),
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={127,127,127},
            fillColor={191,191,191},
            fillPattern=FillPattern.Solid)}));

  end Face;

  connector Inert "Connector to exchange linear momentum and heat by diffusion"

    parameter Integer n_lin(
      final min=0,
      final max=3) = 0
      "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>"
      annotation (HideResult=true);

    Connectors.Mechanical mechanical(final n_lin=n_lin)
      "Subconnector for linear momentum";
    Connectors.Thermal thermal "Subconnector for heat";
    annotation (
      Documentation(info="<html>
    <p>For more information, see the documentation in the
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

    extends FCSys.Connectors.Mechanical;
    extends FCSys.Connectors.Thermal;

    // Additivity of volume
    Q.PressureAbsolute p(nominal=U.atm) "Pressure";
    flow Q.Volume V(min=-Modelica.Constants.inf, nominal=U.cc) "Volume";
    annotation (
      defaultComponentName="inert",
      Documentation(info="<html><p>The concept of \"additivity of volume\" is defined by
    <a href=\"http://en.wikipedia.org/wiki/Amagat's_law\">Amagat's law</a> or the Law of Partial Volumes, which
    states that the partial extensive volumes of the components of a mixture sum to the total
    extensive volume of the mixture [<a href=\"modelica://FCSys.UsersGuide.References\">Bejan2006</a>, p. 194].
    The specific or molar volumes of the species are each evaluated at the temperature and the total pressure of the
    mixture.</p>

    <p>This concept loses its physical meaning once the species are mixed [<a href=\"modelica://FCSys.UsersGuide.References\">Woo1995</a>].
    If the species are truly mixed, then it is impossible to distinguish their particles and thus
    determine their partial volumes.
    Therefore, the concept is only used for distinct phases within the same subregion&mdash;not
    for species within a phase.
    For example, if a system contains a solid phase and a gas phase, then it is assumed that the
    partial volumes of the mixtures are additive and the mixtures exist at the same pressure.  Within
    a phase, the species are mixed according to Dalton's law (see the <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector).</p>

    <p>In order to implement Amagat's law, this connector includes volume (not rate of volume) as a flow variable.
    The effort variable is pressure such that the effort and flow variables are conjugates of
    energy (not power).</p>

    <p>Note that this connector extends from
    the <a href=\"modelica://FCSys.Connectors.BaseClasses.Mechanical\">Mechanical</a> and
    <a href=\"modelica://FCSys.Connectors.BaseClasses.Thermal\">Thermal</a> connectors,
    whereas the <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> and
    <a href=\"modelica://FCSys.Connectors.InertInternal\">InertInternal</a> connectors includes
    them as subconnectors.</p>

    <p>See also the <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> and
    <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connectors and
    the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(graphics={Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={72,90,180},
            fillPattern=FillPattern.Solid,
            fillColor={102,128,255}), Text(
            extent={{-30,20},{30,-20}},
            lineColor={0,0,0},
            textString="A")}),
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

    extends FCSys.Connectors.Mechanical;
    extends FCSys.Connectors.Thermal;

    // Additivity of pressure
    Q.Volume V(nominal=U.cc) "Volume";
    flow Q.Pressure p(nominal=U.atm) "Pressure";
    annotation (
      defaultComponentName="inert",
      Documentation(info="<html><p>The concept of \"additivity of pressure\" is defined by
    <a href=\"http://en.wikipedia.org/wiki/Dalton's_law\">Dalton's law</a> or the Law of Partial Pressures,
    which states that the partial pressures of the components of a mixture sum to the total
    pressure of the mixture [<a href=\"modelica://FCSys.UsersGuide.References\">Bejan2006</a>, p. 192].
    The partial pressures of the species are evaluated at the temperature and the total volume of the
    mixture.</p>

    <p>In order to implement Dalton's law, this connector includes pressure as a flow variable.
    The effort variable is volume such that the effort and flow variables are conjugates of
    energy (not power).</p>

    <p>Note that this connector extends from
    the <a href=\"modelica://FCSys.Connectors.BaseClasses.Mechanical\">Mechanical</a> and
    <a href=\"modelica://FCSys.Connectors.BaseClasses.Thermal\">Thermal</a> connectors,
    whereas the <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> and
    <a href=\"modelica://FCSys.Connectors.InertInternal\">InertInternal</a> connectors includes
    them as subconnectors.</p>

    <p>See also the <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> and
    <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> connectors and
    the documentation in the
  <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(graphics={Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={72,90,180},
            fillPattern=FillPattern.Solid,
            fillColor={102,128,255}), Text(
            extent={{-30,20},{30,-20}},
            lineColor={0,0,0},
            textString="D")}),
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
    parameter Boolean inclMechanical=true "Include the mechanical subconnector"
      annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true));
    parameter Boolean inclThermal=true "Include the thermal subconnector"
      annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true));

    Connectors.Mechanical mechanical(final n_lin=n_lin) if inclMechanical
      "Subconnector for linear momentum";
    Connectors.Thermal thermal if inclThermal "Subconnector for heat";
    annotation (
      defaultComponentPrefixes="protected",
      defaultComponentName="inert",
      Documentation(info="<html><p>This is copy of the <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> connector, except that it
    has a smaller icon, a default <code>protected</code> prefix, and the subconnectors are conditional.
    For more information, see that connector.</p></html>"),
      Icon(graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={72,90,180},
            fillPattern=FillPattern.Solid,
            fillColor={102,128,255}), Ellipse(
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

  connector Mechanical "Connector to exchange or transport linear momentum"

    parameter Integer n_lin(
      final min=0,
      final max=3) = 1
      "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>"
      annotation (HideResult=true);

    Q.Velocity phi[n_lin](each nominal=U.cm/U.s) "Velocity";
    flow Q.Force mPhidot[n_lin](each nominal=U.N) "Force";
    annotation (
      Documentation(info="<html><p>This connector is used to transport the transverse
  orientation of linear momentum only.  The normal orientation is transported by the Normal
  effort-flow pair (see the top-level documentation in the
  <a href=\"modelica://FCSys.Connectors\">Connectors</a> package).</p>
  </html>"),
      Icon(graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={127,127,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}),
      Diagram(graphics={Text(
            extent={{-100,36},{100,76}},
            textString="%name",
            lineColor={0,0,0}), Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={127,127,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));

  end Mechanical;

  connector Thermal "Connector to exchange or transport heat"

    Q.TemperatureAbsolute T(nominal=298.15*U.K) "Temperature";
    flow Q.Power Qdot(nominal=U.W) "Heat flow rate";
    annotation (
      Documentation(info="<html>
    <p>See the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={127,127,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}),
      Diagram(graphics={Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={127,127,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));

    // Note:  The icon doesn't contain a "%name" tag because it overlaps
    // with the tag in the Mechanical connector when both are included in
    // another connector by extension.  In Dymola 7.4, the overlap isn't
    // perfect, so the text appears bolder than it should.

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
            fillPattern=FillPattern.Solid),Text(
            extent={{-200,50},{200,90}},
            textString="%name",
            lineColor={0,0,0})}),
    Documentation(info="<html>
<p>Connector with one input signal of type <code>Real</code>.</p>
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
<p>Protected connector with one input signal of type <code>Real</code>.</p>
</html>"));
  expandable connector RealInputBus
    "<html>Bus of <a href=\"modelica://FCSys.Connectors.RealInput\">RealInput</a> connectors</html>"

    annotation (
      defaultComponentName="u",
      Documentation(info="<html><p>There is no minimal set of variables.
    Signals are included by connecting instances
   of the <a href=\"modelica://FCSys.Connectors.RealInput\">RealInput</a> connector.<p></html>"),

      Icon(graphics={Polygon(
            points={{-102,100},{98,0},{-102,-100},{-102,100}},
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
              lineThickness=0.5),Text(
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
<p>Connector with one output signal of type <code>Real</code>.</p>
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
<p>Protected connector with one output signal of type <code>Real</code>.</p>
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

  package BaseClasses "Base classes (not generally for direct use)"
    extends Modelica.Icons.BasesPackage;

  end BaseClasses;
  annotation (Documentation(info="<html>
  **review and update all doc within the Connectors package.

  <p>Three types of physical connectors are used in <a href=\"modelica://FCSys\">FCSys</a>.
  The chemical connectors (<a href=\"modelica://FCSys.Connectors.ChemicalInput\">ChemicalInput</a>,
  <a href=\"modelica://FCSys.Connectors.ChemicalOutput\">ChemicalOutput</a>,
  <a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a>, and
  <a href=\"modelica://FCSys.Connectors.ChemicalBusInternal\">ChemicalBusInternal</a>)
  represent advective exchange among species that react chemically
  within a subregion.
  The inert connectors
  (<a href=\"modelica://FCSys.Connectors.Inert\">Inert</a>,
  <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a>,
  <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a>, and
  <a href=\"modelica://FCSys.Connectors.InertInternal\">InertInternal</a>)
  describe diffusive (non-chemical) exchange among species within a subregion.
  The face connectors (<a href=\"modelica://FCSys.Connectors.Face\">Face</a>,
  <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a>, and
  <a href=\"modelica://FCSys.Connectors.FaceBusInternal\">FaceBusInternal</a>)
  describe advective and diffusive transport between instances of a single
  species in
  neighboring regions or subregions.</p>

  <p>The effort-flow pairs of the connectors are listed in <a href=\"#Tab1\">Table 1</a>.
  The dimensions of the efforts and flows are noted in terms of mass (M),
  amount or particle number (N),
  length (L), and time (T).  The pairs are used to describe the
  transfer of material, linear momentum, and energy associated with material,
  mechanical, and thermal interactions.  Both advection (e.g., dynamic force and
  thermal convection) and diffusion (e.g., friction and thermal conduction)
  are resolved (see the
  <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model).</p>

  <p>The effort-flow product of the mechanical or transverse
  pair is the rate of energy associated with the interaction.
  However, the other pairs are different.
  The effort of the material pair is electrochemical potential divided by
  temperature because it is used in the chemical equilibrium
  (see the
  <a href=\"modelica://FCSys.Subregions.Reaction\">Reaction</a> model).  The
  fluid pair's effort is massic enthalpy (enthalpy divided by mass) since
  it is used for advection only.  The normal pair has areic current as the
  effort rather than current as the flow because it is more straightforward
  for the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.
  The corresponding flow is force rather than electrochemical potential (as the effort)
  so that the other thermodynamic properties can be calculated without nonlinear systems
  of equations
  (see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model).
  The flow variable of the thermal pair is heat flow rate instead of entropy flow
  rate so that the transport equations are linear and follow the traditional
  form (e.g., Fourier's law; see [<a href=\"modelica://FCSys.UsersGuide.References\">Hogan2006</a>],
  [<a href=\"modelica://FCSys.UsersGuide.References\">Cellier1991</a>], and
  <a href=\"modelica://Modelica.Thermal.HeatTransfer\">Modelica.Thermal.HeatTransfer</a>).
  The additivity of volume and
  pressure pairs describe static balances; therefore, their flows are the quantities
  rather than the rates of the quantities.</p>

    <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"bottom\" id=\"Tab1\">Table 1: Summary of the effort-flow pairs.</caption>
      <tr>
        <th>Within icon(s)</th>
        <th>Name</th>
        <th>Effort</th>
        <th>Flow</th>
      </tr>
      <tr>
        <td>
          <a href=\"modelica://FCSys.Connectors.ChemicalInput\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalInputI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.ChemicalOutput\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalOutputI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.ChemicalBus\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalBusI.png\"></a></td>
          <!--<a href=\"modelica://FCSys.Connectors.ChemicalBusInternal\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalBusInternalI.png\"></a></td>-->
        <td valign=middle>Material</td>
        <td valign=middle>Electrochemical potential<br>&mu; [L<sup>2</sup> M N<sup>-1</sup> T<sup>-2</sup>]</td>
        <td valign=middle>Current<br><i>N&#775;</i> [N T<sup>-1</sup>]</td>
      </tr>
      <tr>
        <td>
          <a href=\"modelica://FCSys.Connectors.Face\"><img src=\"modelica://FCSys/help/FCSys.Connectors.FaceI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.FaceBus\"><img src=\"modelica://FCSys/help/FCSys.Connectors.FaceBusI.png\"></a></td>
        <td valign=middle>Material</td>
        <td valign=middle>Density<br><i>&rho;</i> [N L<sup>-3</sup>]</td>
        <td valign=middle>Current<br><i>N&#775;</i> [N T<sup>-1</sup>]</td>
      </tr>
      <tr>
        <td>
          <a href=\"modelica://FCSys.Connectors.ChemicalInput\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalInputI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.ChemicalOutput\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalOutputI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.ChemicalBus\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalBusI.png\"></a>
          <!--<a href=\"modelica://FCSys.Connectors.ChemicalBusInternal\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalBusInternalI.png\"></a>-->
<br>
          <a href=\"modelica://FCSys.Connectors.Face\"><img src=\"modelica://FCSys/help/FCSys.Connectors.FaceI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.FaceBus\"><img src=\"modelica://FCSys/help/FCSys.Connectors.FaceBusI.png\"></a>
<br>
          <a href=\"modelica://FCSys.Connectors.Inert\"><img src=\"modelica://FCSys/help/FCSys.Connectors.InertI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.InertAmagat\"><img src=\"modelica://FCSys/help/FCSys.Connectors.InertAmagatI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.InertDalton\"><img src=\"modelica://FCSys/help/FCSys.Connectors.InertDaltonI.png\"></a></td>
        <td valign=middle>Mechanical</td>
        <td valign=middle>Velocity<br>&phi; [L T<sup>-1</sup>]</td>
        <td valign=middle>Force<br><i>M</i>&phi;dot [L M T<sup>-2</sup>]</td>
      </tr>
      <tr>
        <td>
          <a href=\"modelica://FCSys.Connectors.ChemicalInput\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalInputI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.ChemicalOutput\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalOutputI.png\"></a>
          <a href=\"modelica://FCSys.Connectors.ChemicalBus\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalBusI.png\"></a></td>
          <!--<a href=\"modelica://FCSys.Connectors.ChemicalBusInternal\"><img src=\"modelica://FCSys/help/FCSys.Connectors.ChemicalBusInternalI.png\"></a></td>-->
        <td valign=middle>Fluid</td>
        <td valign=middle>Temperature times specific entropy<br><i>Ts</i> [L<sup>2</sup> M N<sup>-1</sup> T<sup>-2</sup>]</td>
        <td valign=middle>Heat flow rate<br><i>Q&#775;</i> [L<sup>2</sup> M T<sup>-3</sup>]</td>
      </tr>
      <tr>
        <td>
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
        <td valign=middle>Amagat<br>(additive volume)</td>
        <td valign=middle>Pressure<br><i>p</i> [M L<sup>-1</sup> T<sup>-2</sup>]</td>
        <td valign=middle>Volume<br><i>V</i> [L<sup>3</sup>]</td>
      </tr>
      <tr>
        <td>
        <a href=\"modelica://FCSys.Connectors.InertDalton\"><img src=\"modelica://FCSys/help/FCSys.Connectors.InertDaltonI.png\"></a></td>
        <td valign=middle>Dalton<br>(additive pressure)</td>
        <td valign=middle>Volume<br><i>V</i> [L<sup>3</sup>]</td>
        <td valign=middle>Pressure<br><i>p</i> [M L<sup>-1</sup> T<sup>-2</sup>]</td>
      </tr>
    </table>

  <p><a href=\"#Fig1\">Figure 1</a> depicts the hierarchy of the connectors.
  The bottom row contains icons that represent the effort-flow pairs in
  <a href=\"#Tab1\">Table 1</a>.  Most of the connectors on the middle row are
  flat; they build on the row below by extension.
  The exceptions, which use instantiation, are the
  <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> and
  <a href=\"modelica://FCSys.Connectors.InertInternal\">InertInternal</a> connectors.
  The top row contains expandable connectors
  (<a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> and
  <a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a>) which group
  the <a href=\"modelica://FCSys.Connectors.Face\">Face</a> and chemical (<a href=\"modelica://FCSys.Connectors.ChemicalInput\">ChemicalInput</a> or
  <a href=\"modelica://FCSys.Connectors.ChemicalInput\">ChemicalOutput</a>) connectors of multiple species.
  The
  <a href=\"modelica://FCSys.Connectors.FaceBusInternal\">FaceBusInternal</a> and
  <a href=\"modelica://FCSys.Connectors.ChemicalBusInternal\">ChemicalBusInternal</a>
  connectors (not shown) are versions of
  <a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a> and
  <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a>
  that merely have a default <code>protected</code> prefix and smaller icons to
  indicate that they are internal to a model.</p>

  <p>Although the physical variables are acausal, the
  <a href=\"modelica://FCSys.Connectors.ChemicalInput\">ChemicalInput</a> and
  <a href=\"modelica://FCSys.Connectors.ChemicalOutput\">ChemicalOutput</a> connectors have inputs and
  outputs to pass characteristic data of the species&mdash;the chemical formula and the specific mass.
  That information is used to determine the appropriate
  stoichiometric and advective equations in the
  <a href=\"modelica://FCSys.Subregions.Reaction\">Reaction</a>
  model.</p>

  <p>There are two specialized types of inert
  connectors.  The <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector
  (with an \"A\" in the icon)
  imposes Amagat's law or additivity of volume and is used to combine material phases within a subregion.
  The
  <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector (with a \"D\" in the icon)
  applies Dalton's law or additivity of pressure to mix species within a material phase (e.g.,
  N<sub>2</sub> and O<sub>2</sub> within a gas).
  The two cannot be directly connected because the effort-flow designations
  are opposite.  An adapter must be used
  (e.g., <a href=\"modelica://FCSys.Subregions.PhaseBoundary\">PhaseBoundary</a>).</p>

  <p align=center id=\"Fig1\"><img src=\"modelica://FCSys/resources/documentation/Connectors/hierarchy.png\">
<br>Figure 1: Hierarchy of the connectors.</p>

    <p><b>Relation to thermodynamics and justification of the Amagat and Dalton connectors:</b><p><p>In order to describe the dynamic behavior of a physical system, a model must include conservation laws or rate balances.  These
    equations involve the storage and flow of extensive quantities within (among species) and into the system.  In chemical/thermal systems, the extensive
    quantities of interest are particle number (or mass) and entropy or energy.
    For the sake of simplicity, momentum will be excluded from the present discussion; assume that the fluid is macroscopically stagnant.
    Also assume that there is only one inlet or outlet to the system.
    In terms of mathematics, we have introduced four variables (2 flows and 2 quantities) but only two equations (material and energy
    conservation).</p>

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
    (4 equations in all) and six variables (2 quantities, 2 flows, and 2 efforts or intensive properties).</p>

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
    <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector.</p>

    <p>If the species are mixed, it may be more appropriate to assume that the pressures
    of the components of a mixture sum to the total pressure of the mixture.  This \"additivity of pressure\"
    is described by connections of the
    <a href=\"modelica://FCSys.Connectors.InertDelton\">InertDalton</a> connector (instead of
    instead of <a href=\"modelica://FCSys.Connectors.InertDelton\">InertAmagat</a>).</p>

  <p><b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
  Copyright 2007&ndash;2013, Georgia Tech Research Corporation.</p>

<p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
  it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
  disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.ModelicaLicense2\">
  FCSys.UsersGuide.ModelicaLicense2</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
  http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p>
  </html>"));

end Connectors;
