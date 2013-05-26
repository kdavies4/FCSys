within FCSys;
package Connectors "Declarative and imperative connectors"
  extends Modelica.Icons.InterfacesPackage;

  expandable connector PhysicalBus
    "<html>Bus of <a href=\"modelica://FCSys.Connectors.Chemical\">Chemical</a> connectors</html>**Activity and Potential"
    annotation (
      defaultComponentName="physical",
      Documentation(info="<html><p>See the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Ellipse(
            extent={{-80,80},{80,-80}},
            lineColor={38,196,52},
            lineThickness=0.5), Ellipse(
            extent={{-100,100},{100,-100}},
            fillColor={38,196,52},
            fillPattern=FillPattern.Solid,
            lineColor={2,157,21},
            lineThickness=0.5)}),
      Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
              100,100}}), graphics={Ellipse(
              extent={{-30,30},{30,-30}},
              lineColor={2,157,21},
              fillColor={38,196,52},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.Solid,
              lineThickness=0.5),Text(
              extent={{-120,36},{120,76}},
              textString="%name",
              lineColor={0,0,0})}));

  end PhysicalBus;

  expandable connector PhysicalBusInternal
    "<html>Internal bus of <a href=\"modelica://FCSys.Connectors.Chemical\">Chemical</a> connectors</html>**Activity and Potential"
    annotation (
      defaultComponentName="physical",
      Documentation(info="<html><p>This is copy of the <a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a> connector, except that it
    has a smaller icon and a default <code>protected</code> prefix.  For more information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Ellipse(
            extent={{-80,80},{80,-80}},
            lineColor={38,196,52},
            lineThickness=0.5), Ellipse(
            extent={{-100,100},{100,-100}},
            fillColor={38,196,52},
            fillPattern=FillPattern.Solid,
            lineColor={2,157,21},
            lineThickness=0.5)}),
      Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
              100,100}}), graphics={Ellipse(
              extent={{-10,10},{10,-10}},
              lineColor={2,157,21},
              fillColor={38,196,52},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),Text(
              extent={{-120,20},{120,60}},
              textString="%name",
              lineColor={0,0,0})}));

  end PhysicalBusInternal;

  connector Physical
    "Connector to exchange material while advecting translational momentum and thermal energy**"

    parameter Integer n_trans(
      final min=0,
      final max=3) = 1
      "<html>Number of components of translational momentum (<i>n</i><sub>trans</sub>)</html>"
      annotation (HideResult=true);
    parameter String formula(start="") "Chemical formula of the species";
    // The start value prevents a warning in Dymola 7.4.

    // Material exchange
    Q.Potential a(nominal=1) "Activity";
    flow Q.Current Ndot(nominal=U.A) "Diffusion current";

    // For translational advection
    stream Q.Velocity phi[n_trans](each nominal=U.cm/U.s) "Velocity";

    // For thermal advection
    stream Q.PotentialAbsolute sT(nominal=U.V)
      "Specific entropy-temperature product";
    annotation (
      Documentation(info="<html>
<p>See the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(graphics={Ellipse(extent={{-80,80},{80,-80}}, lineColor={38,196,52}),
            Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={2,157,21},
            fillPattern=FillPattern.Solid,
            fillColor={38,196,52})}),
      Diagram(graphics={Ellipse(
              extent={{-30,30},{30,-30}},
              lineColor={2,157,21},
              fillPattern=FillPattern.Solid,
              fillColor={38,196,52}),Text(
              extent={{-120,36},{120,76}},
              textString="%name",
              lineColor={0,0,0})}));

  end Physical;

  expandable connector ChemicalBus
    "<html>Bus of <a href=\"modelica://FCSys.Connectors.Chemical\">Chemical</a> connectors</html>**Activity and Potential"
    annotation (
      defaultComponentName="chemical",
      Documentation(info="<html><p>See the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Ellipse(
            extent={{-80,80},{80,-80}},
            lineColor={255,195,38},
            lineThickness=0.5), Ellipse(
            extent={{-100,100},{100,-100}},
            fillColor={255,195,38},
            fillPattern=FillPattern.Solid,
            lineColor={239,142,1},
            lineThickness=0.5,
            pattern=LinePattern.Solid)}),
      Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
              100,100}}), graphics={Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={239,142,1},
            fillColor={255,195,38},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.Solid,
            lineThickness=0.5), Text(
            extent={{-120,36},{120,76}},
            textString="%name",
            lineColor={0,0,0})}));

  end ChemicalBus;

  expandable connector ChemicalBusInternal
    "<html>Internal bus of <a href=\"modelica://FCSys.Connectors.Chemical\">Chemical</a> connectors</html>**Activity and Potential"
    annotation (
      defaultComponentName="chemical",
      Documentation(info="<html><p>This is copy of the <a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a> connector, except that it
    has a smaller icon and a default <code>protected</code> prefix.  For more information, see the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Ellipse(
            extent={{-80,80},{80,-80}},
            lineColor={255,195,38},
            lineThickness=0.5), Ellipse(
            extent={{-100,100},{100,-100}},
            fillColor={255,195,38},
            fillPattern=FillPattern.Solid,
            lineColor={239,142,1},
            lineThickness=0.5,
            pattern=LinePattern.Solid)}),
      Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
              100,100}}), graphics={Ellipse(
              extent={{-10,10},{10,-10}},
              lineColor={239,142,1},
              fillColor={255,195,38},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),Text(
              extent={{-120,20},{120,60}},
              textString="%name",
              lineColor={0,0,0})}));

  end ChemicalBusInternal;

  connector ChemicalA
    "Connector to exchange material while advecting translational momentum and thermal energy**"

    parameter Integer n_trans(
      final min=0,
      final max=3) = 1
      "<html>Number of components of translational momentum (<i>n</i><sub>trans</sub>)</html>"
      annotation (HideResult=true);

    output String formula "Chemical formula of the species";
    output Q.MassSpecific m "Specific mass";

    // Material exchange
    Q.Potential mu(nominal=U.V) "Chemical potential";
    flow Q.Current Ndot(nominal=U.A) "Diffusion current";

    // For translational advection
    stream Q.Velocity phi[n_trans](each nominal=U.cm/U.s) "Velocity";

    // For thermal advection
    stream Q.PotentialAbsolute sT(nominal=U.V)
      "Specific entropy-temperature product";
    annotation (
      defaultComponentName="chemical",
      Documentation(info="<html>
<p>See the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(graphics={Ellipse(extent={{-80,80},{80,-80}}, lineColor={255,195,38}),
            Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={239,142,1},
            fillPattern=FillPattern.Solid,
            fillColor={255,195,38})}),
      Diagram(graphics={Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={239,142,1},
            fillPattern=FillPattern.Solid,
            fillColor={255,195,38}), Text(
            extent={{-100,36},{100,76}},
            textString="%name",
            lineColor={0,0,0})}));

  end ChemicalA;

  connector ChemicalB
    "Connector to exchange material while advecting translational momentum and thermal energy**"

    parameter Integer n_trans(
      final min=0,
      final max=3) = 1
      "<html>Number of components of translational momentum (<i>n</i><sub>trans</sub>)</html>"
      annotation (HideResult=true);

    input String formula(start="") "Chemical formula of the species";
    // The start value prevents a warning in Dymola 7.4.
    input Q.MassSpecific m "Specific mass";

    // Material exchange
    Q.Potential mu(nominal=U.V) "Chemical potential";
    flow Q.Current Ndot(nominal=U.A) "Diffusion current";

    // For translational advection
    stream Q.Velocity phi[n_trans](each nominal=U.cm/U.s) "Velocity";

    // For thermal advection
    stream Q.PotentialAbsolute sT(nominal=U.V)
      "Specific entropy-temperature product";
    annotation (
      defaultComponentName="chemical",
      Documentation(info="<html>
<p>See the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(graphics={Ellipse(extent={{-80,80},{80,-80}}, lineColor={255,195,38}),
            Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={239,142,1},
            fillPattern=FillPattern.Solid,
            fillColor={255,255,255})}),
      Diagram(graphics={Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={239,142,1},
            fillPattern=FillPattern.Solid,
            fillColor={255,255,255}), Text(
            extent={{-100,36},{100,76}},
            textString="%name",
            lineColor={0,0,0})}));

  end ChemicalB;

  connector ChemicalReaction
    "Connector to add chemical potentials at a reaction rate"

    // Material exchange
    Q.Current Ndot(nominal=U.A) "Rate of reaction";
    flow Q.Potential mu(nominal=U.V) "Chemical potential";
    annotation (
      defaultComponentName="reaction",
      Documentation(info="<html>
<p>See the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(graphics={
          Ellipse(extent={{-80,80},{80,-80}}, lineColor={255,195,38}),
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={239,142,1},
            fillPattern=FillPattern.Solid,
            fillColor={255,255,255}),
          Ellipse(
            extent={{-60,60},{60,-60}},
            fillColor={255,195,38},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None)}),
      Diagram(graphics={
          Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={239,142,1},
            fillPattern=FillPattern.Solid,
            fillColor={255,255,255}),
          Ellipse(
            extent={{-16,16},{16,-16}},
            fillColor={255,195,38},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Text(
            extent={{-100,36},{100,76}},
            textString="%name",
            lineColor={0,0,0})}));

  end ChemicalReaction;

  expandable connector FaceBus
    "<html>Bus of <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connectors (for multiple configurations)</html>"
    annotation (
      defaultComponentName="face",
      Documentation(info="<html><p>There is no minimal set of variables.  Species are included by connecting instances
    of the <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connector.</p></html>"),

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

  connector Face
    "Connector to transport material, translational momentum, and thermal energy by diffusion"

    // Material
    Q.Density rho(nominal=300*U.K/U.atm) "Density";
    flow Q.Current Ndot(nominal=U.A) "Diffusion current";

    // Translational
    Q.Velocity phi[3](each nominal=U.cm/U.s) "Velocity";
    flow Q.Force mPhidot[3](each nominal=U.N) "Force";

    // Thermal
    extends ThermalDiffusion;
    annotation (
      Documentation(info="<html><p>This connector applies to a single species in a single phase.
    For multiple species or phases, use the <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a>
    connector.  For more information, see the documentation in the
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

  connector Inert
    "Connector to exchange translational momentum and thermal energy by diffusion"

    parameter Integer n_trans(
      final min=0,
      final max=3) = 0
      "<html>Number of components of translational momentum (<i>n</i><sub>trans</sub>)</html>"
      annotation (HideResult=true);

    Translational translational(final n_trans=n_trans)
      "Subconnector for translational diffusion";
    ThermalDiffusion thermal "Subconnector for thermal diffusion";
    annotation (
      Documentation(info="<html>
    <p>See the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(graphics={Text(
              extent={{-100,36},{100,76}},
              textString="%name",
              lineColor={0,0,0}),Ellipse(
              extent={{-30,30},{30,-30}},
              lineColor={11,43,197},
              fillPattern=FillPattern.Solid,
              fillColor={47,107,251})}),
      Icon(graphics={Ellipse(extent={{-80,80},{80,-80}}, lineColor={47,107,251}),
            Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={134,134,201},
            fillPattern=FillPattern.Solid,
            fillColor={47,107,251})}));

  end Inert;

  connector Inert2
    "Connector to exchange translational momentum and thermal energy by diffusion"

    parameter Integer n_trans(
      final min=0,
      final max=3) = 0
      "<html>Number of components of translational momentum (<i>n</i><sub>trans</sub>)</html>"
      annotation (HideResult=true);

    TranslationalSingle translational[n_trans]
      "Subconnector for translational diffusion";
    ThermalDiffusion thermal "Subconnector for thermal diffusion";
    annotation (
      Documentation(info="<html>
    <p>See the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(graphics={Text(
            extent={{-100,36},{100,76}},
            textString="%name",
            lineColor={0,0,0}), Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={11,43,197},
            fillPattern=FillPattern.Solid,
            fillColor={47,107,251})}),
      Icon(graphics={Ellipse(extent={{-80,80},{80,-80}}, lineColor={47,107,251}),
            Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={134,134,201},
            fillPattern=FillPattern.Solid,
            fillColor={47,107,251})}));

  end Inert2;

  connector InertAmagat
    "<html><a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> connector with additivity of volume (direct, without subconnectors)</html>"

    // Additivity of volume
    Q.PressureAbsolute p(nominal=U.atm) "Pressure";
    flow Q.Volume V(min=-Modelica.Constants.inf, nominal=U.cc) "Volume";

    // Translational diffusion
    extends Translational;

    // Thermal diffusion
    //extends ThermalDiffusion;
    Q.TemperatureAbsolute T(nominal=300*U.K) "Temperature";
    flow Q.Power Qdot(nominal=U.W) "Rate of thermal diffusion";
    // Note:  The Thermal connector is copied here rather than included
    // by extension because its "%name" tag overlaps with the tag from
    // the Translational connector.  In Dymola 7.4, the overlap isn't
    // perfect, so the text appears bolder than it should.
    annotation (
      Documentation(info="<html><p>The concept of \"additivity of volume\" is defined by
    <a href=\"http://en.wikipedia.org/wiki/Amagat's_law\">Amagat's law of partial volumes</a>, which
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
    The effort variable is pressure.  This means that the effort and flow variables are conjugates of
    energy (not power).</p>

    <p>See also the
    <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector and
    the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(graphics={
          Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={11,43,197},
            fillPattern=FillPattern.Solid,
            fillColor={47,107,251}),
          Text(
            extent={{-28,26},{24,-26}},
            lineColor={255,255,255},
            textString="A"),
          Text(
            extent={{-24,26},{28,-26}},
            lineColor={255,255,255},
            textString="A"),
          Text(
            extent={{-26,26},{26,-26}},
            lineColor={255,255,255},
            textString="A")}),
      Icon(graphics={
          Ellipse(extent={{-80,80},{80,-80}}, lineColor={47,107,251}),
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={11,43,197},
            fillPattern=FillPattern.Solid,
            fillColor={47,107,251}),
          Ellipse(extent={{-100,100},{100,-100}}, lineColor={11,43,197}),
          Text(
            extent={{-100,96},{100,-96}},
            lineColor={255,255,255},
            textString="A"),
          Text(
            extent={{-96,96},{104,-96}},
            lineColor={255,255,255},
            textString="A"),
          Text(
            extent={{-92,96},{108,-96}},
            lineColor={255,255,255},
            textString="A")}));

  end InertAmagat;

  connector InertDalton
    "<html><a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> connector with additivity of pressure (direct, without subconnectors)</html>"

    // Additivity of pressure
    Q.Volume V(nominal=U.cc) "Volume";
    flow Q.Pressure p(nominal=U.atm) "Pressure";

    // Translational diffusion
    extends Translational;

    // Thermal diffusion
    //extends ThermalDiffusion;
    Q.TemperatureAbsolute T(nominal=300*U.K) "Temperature";
    flow Q.Power Qdot(nominal=U.W) "Rate of thermal diffusion";
    // Note:  The Thermal connector is copied here rather than included
    // by extension because its "%name" tag overlaps with the tag from
    // the Translational connector.  In Dymola 7.4, the overlap isn't
    // perfect, so the text appears bolder than it should.
    annotation (
      Documentation(info="<html><p>The concept of \"additivity of pressure\" is defined by
    <a href=\"http://en.wikipedia.org/wiki/Dalton's_law\">Dalton's law of partial pressures</a>,
    which states that the partial pressures of the components of a mixture sum to the total
    pressure of the mixture [<a href=\"modelica://FCSys.UsersGuide.References\">Bejan2006</a>, p. 192].
    The partial pressures of the species are evaluated at the temperature and the total volume of the
    mixture.</p>

    <p>In order to implement Dalton's law, this connector includes pressure as a flow variable.
    The effort variable is volume.  This means that the effort and flow variables are conjugates of
    energy (not power).</p>

    <p>See also the <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a>
    connector and
    the documentation in the
  <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(graphics={
          Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={11,43,197},
            fillPattern=FillPattern.Solid,
            fillColor={47,107,251}),
          Text(
            extent={{-24,26},{28,-26}},
            lineColor={255,255,255},
            textString="D"),
          Text(
            extent={{-22,26},{30,-26}},
            lineColor={255,255,255},
            textString="D"),
          Text(
            extent={{-24,24},{28,-28}},
            lineColor={255,255,255},
            textString="D"),
          Text(
            extent={{-22,24},{30,-28}},
            lineColor={255,255,255},
            textString="D")}),
      Icon(graphics={
          Ellipse(extent={{-80,80},{80,-80}}, lineColor={47,107,251}),
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={11,43,197},
            fillPattern=FillPattern.Solid,
            fillColor={47,107,251}),
          Text(
            extent={{-94,92},{106,-88}},
            lineColor={255,255,255},
            textString="D"),
          Ellipse(extent={{-100,100},{100,-100}}, lineColor={11,43,197}),
          Text(
            extent={{-90,92},{110,-88}},
            lineColor={255,255,255},
            textString="D"),
          Text(
            extent={{-84,92},{116,-88}},
            lineColor={255,255,255},
            textString="D"),
          Text(
            extent={{-84,86},{116,-94}},
            lineColor={255,255,255},
            textString="D"),
          Text(
            extent={{-88,86},{112,-94}},
            lineColor={255,255,255},
            textString="D"),
          Text(
            extent={{-94,86},{106,-94}},
            lineColor={255,255,255},
            textString="D")}));

  end InertDalton;

  connector InertCond
    "Connector to exchange translational momentum and thermal energy by diffusion**"

    parameter Integer n_trans(
      final min=0,
      final max=3) = 0
      "<html>Number of components of translational momentum (<i>n</i><sub>trans</sub>)</html>"
      annotation (HideResult=true);

    parameter Boolean inclTranslational=true
      "Include the translational subconnector" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true));
    parameter Boolean inclThermal=true "Include the thermal subconnector"
      annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true));

    Translational translational(final n_trans=n_trans) if inclTranslational
      "Subconnector for translational diffusion";
    ThermalDiffusion thermal if inclThermal
      "Subconnector for thermal diffusion";
    annotation (
      Documentation(info="<html>
    <p>See the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(graphics={Text(
            extent={{-100,36},{100,76}},
            textString="%name",
            lineColor={0,0,0}), Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={11,43,197},
            fillPattern=FillPattern.Solid,
            fillColor={47,107,251})}),
      Icon(graphics={Ellipse(extent={{-80,80},{80,-80}}, lineColor={47,107,251}),
            Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={134,134,201},
            fillPattern=FillPattern.Solid,
            fillColor={47,107,251})}));

  end InertCond;

  connector InertInternal
    "<html>Internal <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> connector with conditional subconnectors</html>"

    parameter Integer n_trans(
      final min=0,
      final max=3) = 0
      "<html>Number of components of translational momentum (<i>n</i><sub>trans</sub>)</html>"
      annotation (HideResult=true);
    parameter Boolean inclTranslational=true
      "Include the translational subconnector" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true));
    parameter Boolean inclThermal=true "Include the thermal subconnector"
      annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(compact=true));

    Translational translational(final n_trans=n_trans) if inclTranslational
      "Subconnector for translational diffusion";
    ThermalDiffusion thermal if inclThermal
      "Subconnector for thermal diffusion";
    annotation (
      defaultComponentPrefixes="protected",
      defaultComponentName="inert",
      Documentation(info="<html><p>This is copy of the <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> connector, except that it
    has a smaller icon, a default <code>protected</code> prefix, and the subconnectors are conditional.
    Please see that connector for more information.</p></html>"),
      Icon(graphics={Ellipse(extent={{-80,80},{80,-80}}, lineColor={47,107,251}),
            Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={11,43,197},
            fillPattern=FillPattern.Solid,
            fillColor={47,107,251})}),
      Diagram(graphics={Ellipse(
              extent={{-10,10},{10,-10}},
              lineColor={11,43,197},
              fillPattern=FillPattern.Solid,
              fillColor={47,107,251}),Text(
              extent={{-100,20},{100,60}},
              textString="%name",
              lineColor={0,0,0})}));

  end InertInternal;

  connector Translational
    "Connector for advection or diffusion of translational momentum"

    parameter Integer n_trans(
      final min=0,
      final max=3) = 1
      "<html>Number of components of translational momentum (<i>n</i><sub>trans</sub>)</html>"
      annotation (HideResult=true);

    Q.Velocity phi[n_trans](each nominal=U.cm/U.s) "Velocity";
    flow Q.Force mPhidot[n_trans](each nominal=U.N) "Force";
    annotation (
      Documentation(info="<html><p>See the documentation in the
  <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p>
  </html>"),
      Icon(graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={127,127,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}),
      Diagram(graphics={Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={127,127,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid), Text(
            extent={{-100,36},{100,76}},
            textString="%name",
            lineColor={0,0,0})}));

  end Translational;

  connector TranslationalSingle
    "<html><a href=\"modelica://FCSys.Connectors.Translational\">Translational</a> connector with only one component of translational momentum</html>"

    Q.Velocity phi(nominal=U.cm/U.s) "Velocity";
    flow Q.Force mPhidot(nominal=U.N) "Force";
    annotation (
      Documentation(info="<html><p>See the documentation in the
  <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p>
  </html>"),
      Icon(graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={127,127,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}),
      Diagram(graphics={Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={127,127,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid), Text(
            extent={{-100,36},{100,76}},
            textString="%name",
            lineColor={0,0,0})}));

  end TranslationalSingle;

  connector ThermalDiffusion "Connector for diffusion of thermal energy"

    Q.TemperatureAbsolute T(nominal=300*U.K) "Temperature";
    flow Q.Power Qdot(nominal=U.W) "Rate of thermal diffusion";
    annotation (
      Documentation(info="<html>
    <p>See the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(graphics={Ellipse(extent={{-80,80},{80,-80}}, lineColor={221,23,47}),
            Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={196,11,40},
            fillColor={221,23,47},
            fillPattern=FillPattern.Solid)}),
      Diagram(graphics={Text(
            extent={{-100,36},{100,76}},
            textString="%name",
            lineColor={0,0,0}), Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={170,0,0},
            fillColor={221,23,47},
            fillPattern=FillPattern.Solid)}));

  end ThermalDiffusion;

  connector ThermalAdvection "Connector for advection of thermal energy"

    Q.PotentialAbsolute sT(nominal=3000*U.K)
      "Specific entropy-temperature product";
    flow Q.Power Qdot(nominal=U.W) "Rate of thermal diffusion";
    annotation (
      defaultComponentName="thermal",
      Documentation(info="<html>
    <p>See the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(graphics={
          Ellipse(extent={{-80,80},{80,-80}}, lineColor={221,23,47}),
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={170,0,0},
            fillPattern=FillPattern.Solid,
            fillColor={255,255,255}),
          Ellipse(
            extent={{-60,60},{60,-60}},
            fillColor={221,23,47},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None,
            lineColor={0,0,0})}),
      Diagram(graphics={
          Text(
            extent={{-100,36},{100,76}},
            textString="%name",
            lineColor={0,0,0}),
          Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={170,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-16,16},{16,-16}},
            fillColor={221,23,47},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None,
            lineColor={0,0,0})}));

  end ThermalAdvection;

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
   of the <a href=\"modelica://FCSys.Connectors.RealInput\">RealInput</a> connector.</html>"),

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
   of the <a href=\"modelica://FCSys.Connectors.RealInput\">RealInput</a> connector.</p></html>"),

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
   of the <a href=\"modelica://FCSys.Connectors.RealOutput\">RealOutput</a> connector.</html>"),

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
   of the <a href=\"modelica://FCSys.Connectors.RealOutput\">RealOutput</a> connector.</p></html>"),

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


  annotation (Documentation(info="<html>**Update (chemical and electrical).
  <p>Three types of physical connectors are used in <a href=\"modelica://FCSys\">FCSys</a>.
  The chemical connectors (<a href=\"modelica://FCSys.Connectors.Chemical\">Chemical</a>,
  <a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a>, and
  <a href=\"modelica://FCSys.Connectors.ChemicalBusInternal\">ChemicalBusInternal</a>)
  represent advective exchange among configurations (i.e., species in a particular phase) that
  react chemically or by phase change
  within a subregion.  The inert connectors
  (<a href=\"modelica://FCSys.Connectors.Inert\">Inert</a>,
  <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a>,
  <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a>, and
  <a href=\"modelica://FCSys.Connectors.InertInternal\">InertInternal</a>)
  describe diffusive (non-chemical) exchange among configuration within a subregion.
  The face connectors (<a href=\"modelica://FCSys.Connectors.Face\">Face</a> and
  <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a>)
  describe diffusive transport between instances of a single configuration in
  neighboring regions or subregions. Their variables are sufficient to also
  resolve the rates of advection.</p>

  <p><a href=\"#Fig1\">Figure 1</a> shows the hierarchy of the physical connectors.
  The top row contains <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a>, which expands to group
  the <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connectors of multiple species.
  Most of the connectors on the middle row are
  flat; they build on the connectors of the bottom row by extension.
  The exceptions, which use instantiation, are the
  <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> and
  <a href=\"modelica://FCSys.Connectors.InertInternal\">InertInternal</a> connectors.  The
  <a href=\"modelica://FCSys.Connectors.InertInternal\">InertInternal</a> connector has conditionally
  instantiated subconnectors to directly couple the velocity and temperature of configurations within a subregion.</p>

  The <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> and
  <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connectors have one more effort/flow
  pair than the <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> and
  <a href=\"modelica://FCSys.Connectors.InertInternal\">InertInternal</a> connectors.
  The <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector
  (with an \"A\" in the icon)
  imposes Amagat's law of partial volumes and is used to combine material phases within a subregion.
  The
  <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector (with a \"D\" in the icon)
  applies Dalton's law of partial pressures to mix species within a phase (e.g.,
  N<sub>2</sub> and O<sub>2</sub> within a gas).
  The two cannot be directly connected because the effort/flow designations
  are opposite.  An adapter must be used
  (e.g., <a href=\"modelica://FCSys.Subregions.PhaseBoundary\">PhaseBoundary</a>).</p>

  <p align=center id=\"Fig1\"><img src=\"modelica://FCSys/resources/documentation/Connectors/ConnectorHierarchy.png\">
<br>Figure 1: Hierarchy of the connectors.**Update once final.</p>

  <p>In addition to the physical connectors, there are connectors with inputs and outputs.
  The <a href=\"modelica://FCSys.Connectors.RealInput\">RealInput</a>,
  <a href=\"modelica://FCSys.Connectors.RealInputInternal\">RealInputInternal</a>,
  <a href=\"modelica://FCSys.Connectors.RealInputBus\">RealInputBus</a>, and
  <a href=\"modelica://FCSys.Connectors.RealInputBusInternal\">RealInputBusInternal</a> connectors
  contain only <code>Real input</code> variables.
  The <a href=\"modelica://FCSys.Connectors.RealOutput\">RealOutput</a>,
  <a href=\"modelica://FCSys.Connectors.RealOutputInternal\">RealOutputInternal</a>,
  <a href=\"modelica://FCSys.Connectors.RealOutputBus\">RealOutputBus</a>, and
  <a href=\"modelica://FCSys.Connectors.RealOutputBusInternal\">RealOutputBusInternal</a> connectors
  contain only <code>Real output</code> variables.</p>

    <p><b>Relation to Thermodynamics:</b></p>

    <p>In order to describe the dynamic behavior of a physical system, a model must include conservation
    laws or rate balances.  These
    equations involve the storage and flow of extensive quantities within (among configurations) and into the system.
    In chemical/thermal systems, the extensive
    quantities of interest are particle number (or mass) and energy.
    For the sake of simplicity, momentum will be excluded from the present discussion; assume that the fluid is
    macroscopically stagnant.
    Also assume that there is only one inlet or outlet to the system.
    In terms of mathematics, we have introduced four variables (2 flows and 2 quantities) but only two equations
    (material and energy conservation).</p>

    <p>Two additional equations involve flow rates; these are transport equations with spatial nature&mdash;separate from the temporal
    conservation equations.  Empirical evidence indicates that the
    the flows are related to differences in efforts
    or generalized \"driving forces.\"  The efforts are usually conjugate to the quantities with respect to energy.
    For the chemical/thermal system, the efforts are then electrochemical potential and temperature.  Yet these are intensive
    properties&mdash;distinct from the quantities, which are extensive.
    So far, there are two rate balances to relate extensive quantities to flows and two transport equations
    (4 equations in all) and six variables (2 quantities, 2 flows, and 2 efforts or intensive properties).</p>

    <p>One extensive quantity can be divided by the other to yield an intensive property.
    For example, internal energy can be divided by particle number
    to give internal potential (the relationship is not as direct for electrochemical potential, but the concept is the same).
    The other equation involves the spatial extent of the system, for example, the extensive
    volume of the system divided by particle number to give specific volume.
    This introduces another variable (extensive volume);
    now there are six equations and seven variables.</p>

    <p>In a Eulerian frame of reference, we assume that the extensive volume of the system is fixed (i.e., that the system
    is a \"control volume\").<sup><a href=\"#fn1\" id=\"ref1\">1</a></sup>  If there is only one species in the system,
    then we can assume that it fills the entire volume (e.g., no macroscopically observable regions of vacuum).
    If another species is included in the system, the number of variables is
    doubled.  All of the equations may be repeated except that the specific volume of each species is its own extensive volume or
    \"partial volume\" divided by its particle number to give \"partial specific volume.\"
    It is reasonable to assume that the sum of the partial volumes is equal to the total volume of the system
    (again, no voids).  This is a generalization
    of the previous equation that set the volume of the single species equal to the volume of the system or control volume.
    However, now there are three volumes (of each configuration and of the system) instead of two (of the one configuration and of
    the system) but no additional equations.</p>

    <p>In general,
    an additional equation may be added to exchange volume between the two species such that they reach equilibrium.
    This could be modeled by another transport-like equation.
    However, in the
    <a href=\"modelica://FCSys\">FCSys</a> package, it is assumed that this equilibrium already, always exists.
    Since we wish to impose that the sum of the two
    partial volumes is equal to the total volume, it is appropriate to set the flow variable to be the
    quantity itself (volume) rather than
    the rate of the quantity.  Then, there is no need for another rate balance to relate the quantity to the flow; the quantity <i>is</i>
    the flow.  In this case,
    the most appropriate effort variable is pressure.  The relationship among pressure,
    specific volume, and temperature is given by an equation of state.
    This \"additivity of volume\" interaction occurs through the
    <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector.</p>

    <p>If the species are mixed, it may be more appropriate to assume that the pressures
    of the components of a mixture sum to the total pressure of the mixture.  This \"additivity of pressure\"
    is described by connections of the
    <a href=\"modelica://FCSys.Connectors.InertDelton\">InertDalton</a> connector (instead of
    instead of <a href=\"modelica://FCSys.Connectors.InertDelton\">InertAmagat</a>).</p>

    <hr>

    <small>
    <p id=\"fn1\">1. In a Lagrangian frame of reference, the amount of material is fixed and thermal
    energy is reduced to random motion since particles are tracked directly.  There are only the momentum
    conservation equations.<a href=\"#ref1\" title=\"Jump back to footnote 1 in the text.\">&#8629;</a></p>
    </small>

  <p><b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
  Copyright 2007&ndash;2013, Georgia Tech Research Corporation.</p>

  <p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
  it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
  disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.ModelicaLicense2\">
  FCSys.UsersGuide.ModelicaLicense2</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
  http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p>
  </html>"));
end Connectors;
