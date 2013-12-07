within FCSys;
package Connectors "Declarative and imperative interfaces"
  extends Modelica.Icons.InterfacesPackage;

  connector Chemical
    "Connector for a species in a reaction or phase change process"

    parameter Integer n_trans(min=1,max=3)
      "Number of components of translational momentum" annotation (Evaluate=
          true, Dialog(__Dymola_label="<html><i>n</i><sub>trans</sub></html>"));

    // Material diffusion
    Q.Potential g(nominal=U.V) "Chemical potential";
    flow Q.Current Ndot(nominal=U.A) "Current";

    // For translational advection
    stream Q.Velocity phi[n_trans](each nominal=U.cm/U.s,each start=0)
      "Velocity upon outflow";

    // For thermal advection
    stream Q.PotentialAbsolute sT(nominal=3000*U.K)
      "Specific entropy-temperature product upon outflow";
    annotation (
      Documentation(info="<html>

        <p>Note that the product of the translational <code>stream</code> variable (velocity, &phi;) and the <code>flow</code>
    variable (current, <i>N&#775;</i>) is only conserved among species with the same specific mass.
    This connector is not directly connected among different species.  It is first converted to the

    <a href=\"modelica://FCSys.Connectors.Reaction\">Reaction</a> connector using the

    <a href=\"modelica://FCSys.Conditions.Adapters.ChemicalReaction\">ChemicalReaction</a> adapter,

    which reconciles the different specific masses.</p>

<p>For more information, please see the documentation of the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(graphics={Ellipse(extent={{-80,80},{80,-80}}, lineColor={255,195,38}),
            Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={239,142,1},
            fillPattern=FillPattern.Solid,
            fillColor={255,195,38})}),
      Diagram(graphics={Text(
            extent={{-100,32},{100,72}},
            textString="%name",
            lineColor={0,0,0}), Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={239,142,1},
            fillPattern=FillPattern.Solid,
            fillColor={255,195,38})}));

  end Chemical;

  connector Reaction
    "Connector for the combination of species in a chemical reaction"

    // Material diffusion
    Q.Current Ndot(nominal=U.A,start=0) "Rate of reaction";
    flow Q.Potential g(nominal=U.V) "Chemical potential";

    // Translational advection
    extends Translational;

    // Thermal advection
    Q.PotentialAbsolute sT(nominal=3000*U.K)
      "Product of specific entropy and temperature";
    flow Q.Power Qdot(nominal=U.W) "Rate of thermal advection";

    annotation (
      Documentation(info="<html>

<p>Please see the documentation of the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics={
          Ellipse(extent={{-80,80},{80,-80}}, lineColor={255,195,38}),
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={239,142,1},
            fillPattern=FillPattern.Solid,
            fillColor={255,255,255}),
          Ellipse(
            extent={{-50,50},{50,-50}},
            fillColor={255,195,38},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None,
            lineColor={0,0,0})}),
      Diagram(graphics={Ellipse(
              extent={{-30,30},{30,-30}},
              lineColor={239,142,1},
              fillPattern=FillPattern.Solid,
              fillColor={255,255,255}),Ellipse(
              extent={{-15,15},{15,-15}},
              fillColor={255,195,38},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None,
              lineColor={0,0,0})}));

  end Reaction;

  expandable connector BoundaryBus
    "<html>Bus of <a href=\"modelica://FCSys.Connectors.Boundary\">Boundary</a> connectors (for multiple configurations)</html>"

    annotation (
      defaultComponentName="boundary",
      Documentation(info="<html><p>There is no minimal set of variables.  Species are included by connecting instances
    of the <a href=\"modelica://FCSys.Connectors.Boundary\">Boundary</a> connector.</p>

    <p>For more information, please see the documentation of the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            fillColor={191,191,191},
            fillPattern=FillPattern.Solid,
            lineColor={127,127,127},
            lineThickness=0.5)}),
      Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
              100,100}}), graphics={Text(
            extent={{-100,32},{100,72}},
            textString="%name",
            lineColor={0,0,0}), Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={127,127,127},
            fillColor={191,191,191},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5)}));

  end BoundaryBus;

  connector Boundary
    "Connector to transport material, translational momentum, and thermal energy"

    // Material
    Q.PressureAbsolute p(nominal=U.atm) "Pressure";
    flow Q.Current Ndot(nominal=U.A) "Current";

    // Translational
    Q.Velocity phi[Orient](each nominal=U.cm/U.s,each start=0) "Velocity";
    flow Q.Force mPhidot[Orient](each nominal=U.N) "Force";

    // Thermal
    extends ThermalDiffusive;
    annotation (
      Documentation(info="<html><p>This connector applies to a single species in a single phase.
    For multiple species or phases, use the <a href=\"modelica://FCSys.Connectors.BoundaryBus\">BoundaryBus</a>
    connector.</p>

    <p>For more information, please see the documentation of the
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

  end Boundary;

  connector Amagat "Connector for additivity of volume"

    Q.PressureAbsolute p(nominal=U.atm) "Pressure";
    flow Q.Volume V(min=-Modelica.Constants.inf, nominal=U.cc) "Volume";

    annotation (
      Documentation(info="<html><p>The concept of \"additivity of volume\" is defined by
    <a href=\"http://en.wikipedia.org/wiki/Amagat's_law\">Amagat's law of partial volumes</a>, which
    states that the partial extensive volumes of the components of a mixture sum to the total
    extensive volume of the mixture [<a href=\"modelica://FCSys.UsersGuide.References.Bejan2006\">Bejan2006</a>, p. 194].
    The components are assumed to each exist at the total pressure of the
    mixture.</p>

    <p>This concept loses its physical meaning once the species are mixed [<a href=\"modelica://FCSys.UsersGuide.References.Woo1995\">Woo1995</a>].
    If the species are truly mixed, then it is impossible to distinguish their particles and thus
    determine their partial volumes.
    Therefore, additivity of volume is only used for distinct phases within the same subregion&mdash;not
    for species within a phase.
    For example, if a system contains a solid phase and a gas phase, then it is assumed that the
    volumes of the phases are additive.  Within each phase, the pressures of the species are added according to

    Dalton's law (see the <a href=\"modelica://FCSys.Connectors.Dalton\">Dalton</a> connector).</p>

    <p>In order to implement Amagat's law, this connector includes volume (not rate of volume) as a flow variable.
    The effort variable is pressure.  This implies that the effort and flow variables are conjugates of
    energy (not power).</p>

    <p>See also the
    <a href=\"modelica://FCSys.Connectors.Dalton\">Dalton</a> connector and
    the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(graphics={
          Text(
            extent={{-100,32},{100,72}},
            textString="%name",
            lineColor={0,0,0}),
          Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={11,43,197},
            fillPattern=FillPattern.Solid,
            fillColor={255,255,255}),
          Ellipse(
            extent={{-15,15},{15,-15}},
            fillColor={47,107,251},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None,
            lineColor={0,0,0})}),
      Icon(graphics={
          Ellipse(extent={{-80,80},{80,-80}}, lineColor={47,107,251}),
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={11,43,197},
            fillPattern=FillPattern.Solid,
            fillColor={255,255,255}),
          Ellipse(
            extent={{-50,50},{50,-50}},
            fillColor={47,107,251},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None,
            lineColor={0,0,0})}));

  end Amagat;

  connector Dalton "Connector for additivity of pressure"

    Q.Volume V(nominal=U.cc) "Volume";
    flow Q.Pressure p(nominal=U.atm) "Pressure";

    annotation (
      Documentation(info="<html><p>The concept of \"additivity of pressure\" is defined by
    <a href=\"http://en.wikipedia.org/wiki/Dalton's_law\">Dalton's law of partial pressures</a>,
    which states that the partial pressures of the components of a mixture sum to the total
    pressure of the mixture [<a href=\"modelica://FCSys.UsersGuide.References.Bejan2006\">Bejan2006</a>, p. 192].
    The components are assumed to exist at the total volume of the
    mixture.</p>

    <p>In order to implement Dalton's law, this connector includes pressure as a flow variable.
    The effort variable is volume.  This implies that the effort and flow variables are conjugates of
    energy (not power).</p>

    <p>See also the <a href=\"modelica://FCSys.Connectors.Amagat\">Amagat</a>
    connector and the documentation of the
  <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(graphics={Text(
            extent={{-100,32},{100,72}},
            textString="%name",
            lineColor={0,0,0}), Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={11,43,197},
            fillPattern=FillPattern.Solid,
            fillColor={47,107,251})}),
      Icon(graphics={Ellipse(extent={{-76,76},{84,-84}}, lineColor={47,107,251}),
            Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={11,43,197},
            fillPattern=FillPattern.Solid,
            fillColor={47,107,251})}));

  end Dalton;

  connector DaltonNode "Internal node for additivity of pressure"

    Q.Volume V(nominal=U.cc) "Volume";
    flow Q.Pressure p(nominal=U.atm) "Pressure";

    annotation (
      defaultComponentPrefixes="protected",
      defaultComponentName="dalton",
      Documentation(info="<html><p>This connector is identical to <a href=\"modelica://FCSys.Connectors.Dalton\">Dalton</a> except 
    that it has a smaller icon and a default <code>protected</code> prefix.  For more information, please see that connector.</p></html>"),

      Diagram(graphics={Ellipse(
            extent={{-10,10},{10,-10}},
            lineColor={11,43,197},
            fillColor={47,107,251},
            fillPattern=FillPattern.Solid), Text(
            extent={{-100,10},{100,50}},
            textString="%name",
            lineColor={0,0,0})}),
      Icon(graphics={Ellipse(extent={{-76,76},{84,-84}}, lineColor={47,107,251}),
            Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={11,43,197},
            fillPattern=FillPattern.Solid,
            fillColor={47,107,251})}));

  end DaltonNode;

  connector Intra
    "Connector for translational and thermal exchange among species within a phase"

    extends Translational;

    // extends ThermalDiffusive;
    // Note:  Inheritance isn't used to prevent overlap of the labels in Dymola
    // 2014.
    Q.TemperatureAbsolute T(nominal=300*U.K) "Temperature";
    flow Q.Power Qdot(nominal=U.W) "Rate of thermal conduction";

    annotation (
      Documentation(info="<html><p>This connector is identical to <a href=\"modelica://FCSys.Connectors.Inter\">Inter</a> except for the icon.
    It should not be connected directly among species.  Rather, it should be connected to

    an <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> or

    <a href=\"modelica://FCSys.Connectors.InertNode\">InertNode</a> connector so that <i>k</i><sub>&Phi;</sub>

    and <i>k<sub>Q</sub></i> (essentially parameters) are set.</p>

    <p>Please see <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a>, <a href=\"modelica://FCSys.Connectors.InertNode\">InertNode</a>,

    and the documentation of the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(graphics={Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={170,0,0},
            fillPattern=FillPattern.Solid,
            fillColor={221,23,47}), Text(
            extent={{-30,30},{30,-30}},
            lineColor={255,255,255},
            textStyle={TextStyle.Bold},
            textString="1")}),
      Icon(graphics={
          Ellipse(extent={{-80,80},{80,-80}}, lineColor={221,23,47}),
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={170,0,0},
            fillPattern=FillPattern.Solid,
            fillColor={221,23,47}),
          Text(
            extent={{-100,100},{100,-100}},
            lineColor={255,255,255},
            textString="1",
            textStyle={TextStyle.Bold})}));

  end Intra;

  connector Inter
    "Connector for translational and thermal exchange among species, regardless of phase"

    extends Translational;

    // extends ThermalDiffusive;
    // Note:  Inheritance isn't used to prevent overlap of the labels in
    // Dymola 2014.
    Q.TemperatureAbsolute T(nominal=300*U.K) "Temperature";
    flow Q.Power Qdot(nominal=U.W) "Rate of thermal conduction";

    annotation (
      Documentation(info="<html><p>This connector is identical to <a href=\"modelica://FCSys.Connectors.Intra\">Intra</a> except for the icon.
    It should not be connected directly among species.  Rather, it should be connected to

    an <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> or

    <a href=\"modelica://FCSys.Connectors.InertNode\">InertNode</a> connector so that <i>k</i><sub>&Phi;</sub>

    and <i>k<sub>Q</sub></i> (essentially parameters) are set.</p>

    <p>Please see <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a>, <a href=\"modelica://FCSys.Connectors.InertNode\">InertNode</a>,

    and the documentation of the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(graphics={Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={170,0,0},
            fillPattern=FillPattern.Solid,
            fillColor={221,23,47}), Text(
            extent={{-30,30},{30,-30}},
            lineColor={255,255,255},
            textStyle={TextStyle.Bold},
            textString="2")}),
      Icon(graphics={
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={127,127,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(extent={{-80,80},{80,-80}}, lineColor={221,23,47}),
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={170,0,0},
            fillPattern=FillPattern.Solid,
            fillColor={221,23,47}),
          Text(
            extent={{-100,100},{100,-100}},
            lineColor={255,255,255},
            textStyle={TextStyle.Bold},
            textString="2")}));

  end Inter;

  connector Inert "Connector for translational and thermal exchange"

    extends Translational;

    // extends ThermalDiffusive;
    // Note:  Inheritance isn't used to prevent overlap of the labels in
    // Dymola 2014.
    Q.TemperatureAbsolute T(nominal=300*U.K) "Temperature";
    flow Q.Power Qdot(nominal=U.W) "Rate of thermal conduction";

    annotation (
      Documentation(info="<html><p>This connector is identical to <a href=\"modelica://FCSys.Connectors.InertNode\">InertNode</a> except for the icon.
    It is used as a node to connect

    <a href=\"modelica://FCSys.Connectors.Intra\">Intra</a> or <a href=\"modelica://FCSys.Connectors.Inter\">Inter</a>
    connectors.  It contains outputs (<i>k</i><sub>&Phi;</sub> and <i>k<sub>Q</sub></i>) that are essentially
    parameters to scale the translational and thermal independences of the connection.</p>

    <p>Please see the documentation of the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics={Ellipse(
              extent={{-30,30},{30,-30}},
              lineColor={170,0,0},
              fillPattern=FillPattern.Solid,
              fillColor={221,23,47})}),
      Icon(graphics={Ellipse(extent={{-80,80},{80,-80}}, lineColor={221,23,47}),
            Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={170,0,0},
            fillPattern=FillPattern.Solid,
            fillColor={221,23,47})}));

  end Inert;

  expandable connector InertNode
    "<html>Internal node for <a href=\"modelica://FCSys.Connectors.Intra\">Intra</a> and <a href=\"modelica://FCSys.Connectors.Intra\">Inter</a></html>"

    annotation (
      Documentation(info="<html><p>This connector is used as an internal node to connect

    <a href=\"modelica://FCSys.Connectors.Intra\">Intra</a> or <a href=\"modelica://FCSys.Connectors.Inter\">Inter</a>
    connectors.  It is an expandable connector with no intrinsic variables.</p>

    <p>For more information, please see the documentation of the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics={Ellipse(extent={{-80,80},{80,-80}}, lineColor={
                221,23,47}), Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={170,0,0},
            fillPattern=FillPattern.Solid,
            fillColor={221,23,47})}),
      Diagram(coordinateSystem(preserveAspectRatio=true,extent={{-100,-100},{
              100,100}}), graphics={Ellipse(
            extent={{-10,10},{10,-10}},
            lineColor={170,0,0},
            fillColor={221,23,47},
            fillPattern=FillPattern.Solid,
            lineThickness=0.5), Text(
            extent={{-100,10},{100,50}},
            textString="%name",
            lineColor={0,0,0})}));

  end InertNode;

  connector Translational
    "Connector for the advection or diffusion of translational momentum"

    parameter Integer n_trans(min=1,max=3)
      "Number of components of translational momentum" annotation (HideResult=
          true,Dialog(__Dymola_label="<html><i>n</i><sub>trans</sub></html>"));

    Q.Velocity phi[n_trans](each nominal=U.cm/U.s,each start=0) "Velocity";
    flow Q.Force mPhidot[n_trans](each nominal=U.N) "Force";

    annotation (
      defaultComponentName="trans",
      Documentation(info="<html><p>Please see the documentation of the
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
            extent={{-100,32},{100,72}},
            textString="%name",
            lineColor={0,0,0})}));

  end Translational;

  connector ThermalDiffusive "Connector for the diffusion of thermal energy"

    Q.TemperatureAbsolute T(nominal=300*U.K) "Temperature";
    flow Q.Power Qdot(nominal=U.W) "Rate of thermal conduction";
    annotation (
      defaultComponentName="therm",
      Documentation(info="<html>
    <p>Please see the documentation of the
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
            fillPattern=FillPattern.Solid), Text(
            extent={{-100,32},{100,72}},
            textString="%name",
            lineColor={0,0,0})}));

  end ThermalDiffusive;

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
          fillPattern=FillPattern.Solid), Text(
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
          fillPattern=FillPattern.Solid), Text(
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
            lineThickness=0.5), Text(
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

  annotation (Documentation(info="<html>
  <p><a href=\"modelica://FCSys\">FCSys</a> uses four types of declarative connectors.
  The chemical connectors (<a href=\"modelica://FCSys.Connectors.Chemical\">Chemical</a> and
  <a href=\"modelica://FCSys.Connectors.Reaction\">Reaction</a>)
  represent the diffusion of material and the advection of other quantities among
  configurations (i.e., species in particular phases) that
  react chemically within a subregion.
  The inert connectors
  (<a href=\"modelica://FCSys.Connectors.Intra\">Intra</a>,
  <a href=\"modelica://FCSys.Connectors.Inter\">Inter</a>,

  <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a>, and
  <a href=\"modelica://FCSys.Connectors.InertNode\">InertNode</a>)
  describe the diffusive exchange of momentum and energy among configurations within a subregion.
  The mixing connectors (<a href=\"modelica://FCSys.Connectors.Amagat\">Amagat</a> and
  <a href=\"modelica://FCSys.Connectors.Dalton\">Dalton</a>) describe

  how species are combined within a phase and how phases are combined within a subregion.
  The boundary connectors (<a href=\"modelica://FCSys.Connectors.Boundary\">Boundary</a> and
  <a href=\"modelica://FCSys.Connectors.BoundaryBus\">BoundaryBus</a>)
  describe the transport between neighboring regions or subregions.</p>

  <p><a href=\"#Fig1\">Figure 1</a> shows the hierarchy of the declarative connectors.
  The top row contains a bus connector
  (<a href=\"modelica://FCSys.Connectors.BoundaryBus\">BoundaryBus</a>), which
  expands to group the

  <a href=\"modelica://FCSys.Connectors.Boundary\">Boundary</a> connectors of multiple species.
  The connectors in the middle row are
  flat; they build on the connectors of the bottom row by extension.
  Each icon on the bottom row represents one effort/flow
  pair, which may or may not be implemented as a separate connector.

  The

  <a href=\"modelica://FCSys.Connectors.Chemical\">Chemical</a> connector also
  has stream variables to represent the advection of translational momentum and
  thermal energy.</p>

  <p align=center id=\"Fig1\"><img src=\"modelica://FCSys/Resources/Documentation/Connectors.png\">
<br>Figure 1: Hierarchy of the connectors.</p>

  <p>The <a href=\"modelica://FCSys.Connectors.Chemical\">Chemical</a>
  connector is used for a single species in a chemical reaction or phase change process.

  It expresses the rate of
  consumption or generation of a species at a chemical potential.  The
  <a href=\"modelica://FCSys.Connectors.Reaction\">Reaction</a> connector is
  used for the chemical reaction as a whole.  It has chemical potential as a
  flow and current as an effort (opposite designations of the
  <a href=\"modelica://FCSys.Connectors.Chemical\">Chemical</a>
  connector).
  It sums the stoichiometrically weighted chemical potentials of the species
  participating in a reaction.  Its effort variable is the
  rate of the reaction.  <a href=\"modelica://FCSys.Connectors.Chemical\">Chemical</a> and
  <a href=\"modelica://FCSys.Connectors.Reaction\">Reaction</a> cannot be directly connected because they have opposite efforts and flows.

  An adapter must be used
  (e.g., <a href=\"modelica://FCSys.Conditions.Adapters.ChemicalReaction\">ChemicalReaction</a>),
  which applies the stoichiometry.</p>

  The <a href=\"modelica://FCSys.Connectors.Amagat\">Amagat</a> connector
  imposes Amagat's law of partial volumes and is used to combine material phases within a subregion.
  The
  <a href=\"modelica://FCSys.Connectors.Dalton\">Dalton</a> connector
  applies Dalton's law of partial pressures to mix species within a phase (e.g.,
  N<sub>2</sub> and O<sub>2</sub> within a gas).
  The two cannot be directly connected because they have opposite efforts and flows.  An adapter must be used
  (e.g., <a href=\"modelica://FCSys.Conditions.Adapters.AmagatDalton\">AmagatDalton</a>).</p>

  <p>In addition to the declarative connectors, there are connectors with inputs and outputs.
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
    For the chemical/thermal system, the efforts are then chemical potential and temperature.  Yet these are intensive
    properties&mdash;distinct from the quantities, which are extensive.
    So far, there are two rate balances to relate extensive quantities to flows and two transport equations
    (4 equations in all) and six variables (2 quantities, 2 flows, and 2 efforts or intensive properties).</p>

    <p>One extensive quantity can be divided by the other to yield an intensive property.
    For example, internal energy can be divided by particle number
    to give internal potential (the relationship is not as direct for chemical potential, but the concept is the same).
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
    This additivity-of-volume interaction occurs in the
    <a href=\"modelica://FCSys.Conditions.ByConnector.Amagat.VolumeFixed\">VolumeFixed</a> model.</p>

    <p>If the species are mixed, it may be more appropriate to assume that the pressures
    of the components of a mixture sum to the total pressure of the mixture.  This additivity of pressure
    is described by connections of the
    <a href=\"modelica://FCSys.Connectors.Dalton\">Dalton</a> connector.</p>

    <hr>

    <small>
    <p id=\"fn1\">1. In a Lagrangian frame of reference, the amount of material is fixed and thermal
    energy is reduced to random motion since particles are tracked directly.  There are only the momentum
    conservation equations.<a href=\"#ref1\" title=\"Jump back to footnote 1 in the text.\">&#8629;</a></p>
    </small>

  <p><b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
  Copyright 2007&ndash;2013, <a href=\"http://www.gtrc.gatech.edu/\">Georgia Tech Research Corporation</a>.</p>

  <p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
  it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.License\">
FCSys.UsersGuide.License</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
  http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p>
  </html>"));

end Connectors;
