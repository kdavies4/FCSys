within FCSys;
package Phases "Mixtures of species"
  extends Modelica.Icons.Package;

  package Examples "Examples"
    extends Modelica.Icons.ExamplesPackage;

    model TestCapacitor
      "<html>Test the <a href=\"modelica://FCSys.Phases.Examples.Capacitor\">Capacitor</a> model</html>"
      extends Modelica.Icons.Example;
      Modelica.Electrical.Analog.Basic.Capacitor capacitorMSL(C=1, v(fixed=true,
            start=1)) "Modelica electrical capacitor" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-60,0})));
      Modelica.Electrical.Analog.Basic.Resistor resistor1(R=1) annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-20,0})));
      Modelica.Electrical.Analog.Basic.Ground ground1
        annotation (Placement(transformation(extent={{-50,-50},{-30,-30}})));
      Capacitor capacitorFCSys(C=1, v(fixed=true, start=1))
        "Capacitor using FCSys" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={20,0})));
      Modelica.Electrical.Analog.Basic.Resistor resistor2(R=1) annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={60,0})));
      Modelica.Electrical.Analog.Basic.Ground ground2
        annotation (Placement(transformation(extent={{30,-50},{50,-30}})));
    equation

      connect(capacitorMSL.p, resistor1.p) annotation (Line(
          points={{-60,10},{-60,20},{-20,20},{-20,10}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(capacitorMSL.n, ground1.p) annotation (Line(
          points={{-60,-10},{-60,-20},{-40,-20},{-40,-30}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(resistor1.n, ground1.p) annotation (Line(
          points={{-20,-10},{-20,-20},{-40,-20},{-40,-30}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(capacitorFCSys.p, resistor2.p) annotation (Line(
          points={{20,10},{20,20},{60,20},{60,10}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(capacitorFCSys.n, ground2.p) annotation (Line(
          points={{20,-10},{20,-20},{40,-20},{40,-30}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(resistor2.n, ground2.p) annotation (Line(
          points={{60,-10},{60,-20},{40,-20},{40,-30}},
          color={0,0,255},
          smooth=Smooth.None));
      annotation (
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics),
        experiment(StopTime=5),
        __Dymola_experimentSetupOutput,
        Commands(file=
              "Resources/Scripts/Dymola/Phases.Examples.TestCapacitor.mos"
            "Phases.Examples.TestCapacitor.mos"));
    end TestCapacitor;

    model Capacitor
      "<html><a href=\"modelica://Modelica\">Modelica</a>-equivalent capacitor built using the <a href=\"modelica://FCSys.Phases.Dielectric\">Dielectric</a> model</html>"
      import SI = Modelica.SIunits;
      extends Modelica.Electrical.Analog.Interfaces.TwoPin(v(start=0));

      parameter SI.Capacitance C(start=1) "Capacitance in farads";
      output SI.Current i=p.i "Current in amperes";

      Dielectric dielectric(final epsilon=(C*U.F)*dielectric.L/dielectric.A)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Conditions.ByConnector.Amagat.Pressure pressure
        "Required to terminate the amagat connector of the dielectric"
        annotation (Placement(transformation(extent={{-10,-6},{10,-26}})));

      Plate positivePlate
        annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
      Plate negativePlate(N(stateSelect=StateSelect.never))
        annotation (Placement(transformation(extent={{60,-10},{40,10}})));
      // Note:  This avoids dynamic state selection in Dymola 2014.
    equation
      connect(pressure.amagat, dielectric.amagat) annotation (Line(
          points={{0,-12},{0,0}},
          color={47,107,251},
          smooth=Smooth.None));
      connect(dielectric.positive, negativePlate.electrical) annotation (Line(
          points={{6,0},{46,0}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(dielectric.negative, positivePlate.electrical) annotation (Line(
          points={{-6,0},{-46,0}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(positivePlate.pin, p) annotation (Line(
          points={{-54,0},{-100,0}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(negativePlate.pin, n) annotation (Line(
          points={{54,0},{100,0}},
          color={0,0,255},
          smooth=Smooth.None));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{
                -100,-100},{100,100}}), graphics), Icon(graphics={Line(points={
              {-6,28},{-6,-28}}, color={0,0,255}),Line(points={{-90,0},{-6,0}},
              color={0,0,255}),Line(points={{6,0},{90,0}}, color={0,0,255}),
              Line(points={{6,28},{6,-28}}, color={0,0,255}),Text(
                  extent={{-136,-60},{136,-92}},
                  lineColor={0,0,0},
                  textString="C=%C_SI F"),Text(
                  extent={{-150,85},{150,45}},
                  textString="%name",
                  lineColor={0,0,255})}));
    end Capacitor;

    model Plate "Plate of a capacitor (for charge storage)"

      extends FCSys.Icons.Names.Top1;

      parameter Q.Number z=1 "Charge number";
      Q.Amount N(start=0) "Amount of material";

      Modelica.Electrical.Analog.Interfaces.Pin pin "Modelica electrical pin"
        annotation (Placement(transformation(extent={{-50,-10},{-30,10}}),
            iconTransformation(extent={{-50,-10},{-30,10}})));
      Connectors.Electrostatic electrical "FCSys electrical connector"
        annotation (Placement(transformation(extent={{30,-10},{50,10}}),
            iconTransformation(extent={{30,-10},{50,10}})));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{
                -100,-100},{100,100}}), graphics), Icon(coordinateSystem(
              preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
            graphics={Polygon(
                  points={{-20,0},{0,20},{20,0},{0,-20},{-20,0}},
                  lineColor={127,127,127},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),Line(
                  points={{20,0},{30,0}},
                  color={255,195,38},
                  smooth=Smooth.None),Line(
                  points={{-30,0},{-20,0}},
                  color={0,0,255},
                  smooth=Smooth.None)}));

    equation
      // Aliases
      electrical.Z = z*N;

      // Equal potential
      electrical.w = pin.v*U.V;

      // Conservation
      der(N)/U.s = (pin.i*U.A)/z "Material";

    end Plate;
  end Examples;

  model Gas "Gas phase"
    import Modelica.Math.BooleanVectors.countTrue;

    extends Phase(final n_spec=countTrue({inclH2,inclH2O,inclN2,inclO2}));

    // Conditionally include species.
    parameter Boolean inclH2=false "Include H2" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Hydrogen (H<sub>2</sub>)</html>",
        __Dymola_joinNext=true));
    replaceable FCSys.Species.H2.Gas.Fixed H2(final n_trans) if inclH2
      constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "H2 model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>H<sub>2</sub> model</html>",
        enable=inclH2),
      Placement(transformation(extent={{-70,-10},{-50,10}})));

    parameter Boolean inclH2O=false "Include H2O" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Water (H<sub>2</sub>O)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.H2O.Gas.Fixed H2O(final n_trans) if inclH2O
      constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "H2O model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>H<sub>2</sub>O model</html>",
        enable=inclH2O),
      Placement(transformation(extent={{-30,-10},{-10,10}})));

    parameter Boolean inclN2=false "Include N2" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Nitrogen (N<sub>2</sub>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.N2.Gas.Fixed N2(final n_trans) if inclN2
      constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "N2 model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>N<sub>2</sub> model</html>",
        enable=inclN2),
      Placement(transformation(extent={{10,-10},{30,10}})));

    parameter Boolean inclO2=false "Include O2" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Oxygen (O<sub>2</sub>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.O2.Gas.Fixed O2(final n_trans) if inclO2
      constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "O2 model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>O<sub>2</sub> model</html>",
        enable=inclO2),
      Placement(transformation(extent={{50,-10},{70,10}})));

    Connectors.FaceBus xNegative if inclTrans[Axis.x]
      "Negative face along the x axis" annotation (Placement(transformation(
            extent={{-120,-10},{-100,10}}), iconTransformation(extent={{-90,-10},
              {-70,10}})));
    Connectors.FaceBus yNegative if inclTrans[Axis.y]
      "Negative face along the y axis" annotation (Placement(transformation(
            extent={{-96,-34},{-76,-14}}), iconTransformation(extent={{-10,-94},
              {10,-74}})));
    Connectors.FaceBus zNegative if inclTrans[Axis.z]
      "Negative face along the z axis" annotation (Placement(transformation(
            extent={{88,2},{108,22}}),iconTransformation(extent={{40,40},{60,60}})));
    Connectors.FaceBus xPositive if inclTrans[Axis.x]
      "Positive face along the x axis" annotation (Placement(transformation(
            extent={{100,-10},{120,10}}), iconTransformation(extent={{70,-10},{
              90,10}})));
    Connectors.FaceBus yPositive if inclTrans[Axis.y]
      "Positive face along the y axis" annotation (Placement(transformation(
            extent={{76,14},{96,34}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.FaceBus zPositive if inclTrans[Axis.z]
      "Positive face along the z axis" annotation (Placement(transformation(
            extent={{-108,-22},{-88,-2}}), iconTransformation(extent={{-90,-90},
              {-70,-70}})));
    Connectors.Amagat amagat(final V=-V) if n_spec > 0
      "Connector for additivity of volume" annotation (Placement(transformation(
            extent={{-120,26},{-100,46}}), iconTransformation(extent={{70,-90},
              {90,-70}})));
    Connectors.Inter inter[n_inter](each final n_trans=n_trans) if n_spec > 0
      "Connector to exchange momentum and energy with other phases" annotation
      (Placement(transformation(extent={{100,-46},{120,-26}}),
          iconTransformation(extent={{-60,60},{-40,40}})));

    // Auxiliary variables (for analysis)
    output Q.PressureAbsolute p(stateSelect=StateSelect.never) = amagat.p if
      n_spec > 0 and environment.analysis "Total thermodynamic pressure";

    // Aliases
    Q.Velocity phi[n_trans](each stateSelect=StateSelect.prefer) = direct.translational.phi
      if n_spec > 1 and reduceTrans "Velocity";
    Q.Temperature T(each stateSelect=StateSelect.prefer) = direct.thermal.T if
      n_spec > 1 and reduceThermal "Temperature";
    // These make the selected states more readable.

    Connectors.Chemical connH2(final n_trans=n_trans) if inclH2
      "Chemical connector for H2" annotation (Placement(transformation(extent={
              {-74,40},{-54,60}}), iconTransformation(extent={{-50,-50},{-30,-30}})));
    Connectors.Chemical connH2O(final n_trans=n_trans) if inclH2O
      "Chemical connector for H2O" annotation (Placement(transformation(extent=
              {{-34,40},{-14,60}}), iconTransformation(extent={{-10,-50},{10,-30}})));
    Connectors.Chemical connO2(final n_trans=n_trans) if inclO2
      "Chemical connector for O2" annotation (Placement(transformation(extent={
              {46,40},{66,60}}), iconTransformation(extent={{30,-50},{50,-30}})));

  protected
    Conditions.Adapters.AmagatDalton amagatDalton if n_spec > 0
      "Adapter between additivity of volume and additivity of pressure"
      annotation (Placement(transformation(extent={{-100,26},{-80,46}})));
    FCSys.Connectors.DirectNode direct(
      final n_trans=n_trans,
      final inclTrans=reduceTrans,
      final inclThermal=reduceThermal) if n_spec > 0 and (reduceTrans or
      reduceThermal)
      "Connector to directly couple velocities and temperatures within the phase"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={
              4,-48}), iconTransformation(extent={{-10,-10},{10,10}}, origin={
              26,-52})));

  equation
    // Chemical exchange
    connect(O2.chemical, connO2) annotation (Line(
        points={{56,9},{56,50}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(H2.chemical, connH2) annotation (Line(
        points={{-64,9},{-64,50}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(H2O.chemical, connH2O) annotation (Line(
        points={{-24,9},{-24,50}},
        color={221,23,47},
        smooth=Smooth.None));

    // Inert exchange
    connect(H2.inter, inter) annotation (Line(
        points={{-51,-3.8},{-49,-6},{-49,-36},{110,-36}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2O.inter, inter) annotation (Line(
        points={{-11,-3.8},{-9,-6},{-9,-36},{110,-36}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(N2.inter, inter) annotation (Line(
        points={{29,-3.8},{31,-6},{31,-36},{110,-36}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(O2.inter, inter) annotation (Line(
        points={{69,-3.8},{71,-6},{71,-36},{110,-36}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2.direct.translational, direct.translational) annotation (Line(
        points={{-56.2,-9},{-56.2,-48},{4,-48}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2.direct.thermal, direct.thermal) annotation (Line(
        points={{-56.2,-9},{-56.2,-48},{4,-48}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2O.direct.translational, direct.translational) annotation (Line(
        points={{-16.2,-9},{-16.2,-48},{4,-48}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2O.direct.thermal, direct.thermal) annotation (Line(
        points={{-16.2,-9},{-16.2,-48},{4,-48}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(N2.direct.translational, direct.translational) annotation (Line(
        points={{23.8,-9},{23.8,-48},{4,-48}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(N2.direct.thermal, direct.thermal) annotation (Line(
        points={{23.8,-9},{23.8,-48},{4,-48}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(O2.direct.translational, direct.translational) annotation (Line(
        points={{63.8,-9},{63.8,-48},{4,-48}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(O2.direct.thermal, direct.thermal) annotation (Line(
        points={{63.8,-9},{63.8,-48},{4,-48}},
        color={38,196,52},
        smooth=Smooth.None));

    // Mixing
    connect(H2.dalton, amagatDalton.dalton) annotation (Line(
        points={{-69,4},{-69,36},{-86,36}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(H2O.dalton, amagatDalton.dalton) annotation (Line(
        points={{-29,4},{-29,36},{-86,36}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(N2.dalton, amagatDalton.dalton) annotation (Line(
        points={{11,4},{11,36},{-86,36}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(O2.dalton, amagatDalton.dalton) annotation (Line(
        points={{51,4},{51,36},{-86,36}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(amagatDalton.amagat, amagat) annotation (Line(
        points={{-94,36},{-110,36}},
        color={47,107,251},
        smooth=Smooth.None));

    // Transport
    // ---------
    // H2
    connect(H2.faces[transCart[Axis.x], Side.n], xNegative.H2) annotation (Line(
        points={{-60,6.10623e-016},{-110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2.faces[transCart[Axis.x], Side.p], xPositive.H2) annotation (Line(
        points={{-60,6.10623e-016},{110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2.faces[transCart[Axis.y], Side.n], yNegative.H2) annotation (Line(
        points={{-60,6.10623e-016},{-60,-24},{-86,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2.faces[transCart[Axis.y], Side.p], yPositive.H2) annotation (Line(
        points={{-60,6.10623e-016},{-60,24},{86,24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2.faces[transCart[Axis.z], Side.n], zNegative.H2) annotation (Line(
        points={{-60,6.10623e-016},{-48,12},{98,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(H2.faces[transCart[Axis.z], Side.p], zPositive.H2) annotation (Line(
        points={{-60,6.10623e-016},{-72,-12},{-98,-12}},
        color={127,127,127},
        smooth=Smooth.None));
    // H2O
    connect(H2O.faces[transCart[Axis.x], Side.n], xNegative.H2O) annotation (
        Line(
        points={{-20,6.10623e-016},{-110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[transCart[Axis.x], Side.p], xPositive.H2O) annotation (
        Line(
        points={{-20,6.10623e-016},{110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[transCart[Axis.y], Side.n], yNegative.H2O) annotation (
        Line(
        points={{-20,6.10623e-016},{-20,-24},{-86,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[transCart[Axis.y], Side.p], yPositive.H2O) annotation (
        Line(
        points={{-20,6.10623e-016},{-20,24},{86,24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[transCart[Axis.z], Side.n], zNegative.H2O) annotation (
        Line(
        points={{-20,6.10623e-016},{-8,12},{98,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(H2O.faces[transCart[Axis.z], Side.p], zPositive.H2O) annotation (
        Line(
        points={{-20,6.10623e-016},{-32,-12},{-98,-12}},
        color={127,127,127},
        smooth=Smooth.None));
    // N2
    connect(N2.faces[transCart[Axis.x], Side.n], xNegative.N2) annotation (Line(
        points={{20,6.10623e-016},{-110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(N2.faces[transCart[Axis.x], Side.p], xPositive.N2) annotation (Line(
        points={{20,6.10623e-016},{110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(N2.faces[transCart[Axis.y], Side.n], yNegative.N2) annotation (Line(
        points={{20,6.10623e-016},{20,-24},{-86,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(N2.faces[transCart[Axis.y], Side.p], yPositive.N2) annotation (Line(
        points={{20,6.10623e-016},{20,24},{86,24}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(N2.faces[transCart[Axis.z], Side.n], zNegative.N2) annotation (Line(
        points={{20,6.10623e-016},{20,0},{32,12},{98,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(N2.faces[transCart[Axis.z], Side.p], zPositive.N2) annotation (Line(
        points={{20,6.10623e-016},{8,-12},{-98,-12}},
        color={127,127,127},
        smooth=Smooth.None));
    // O2
    connect(O2.faces[transCart[Axis.x], Side.n], xNegative.O2) annotation (Line(
        points={{60,6.10623e-016},{-110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(O2.faces[transCart[Axis.x], Side.p], xPositive.O2) annotation (Line(
        points={{60,6.10623e-016},{110,5.55112e-016}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(O2.faces[transCart[Axis.y], Side.n], yNegative.O2) annotation (Line(
        points={{60,6.10623e-016},{60,-24},{-86,-24}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(O2.faces[transCart[Axis.y], Side.p], yPositive.O2) annotation (Line(
        points={{60,6.10623e-016},{60,24},{86,24}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(O2.faces[transCart[Axis.z], Side.n], zNegative.O2) annotation (Line(
        points={{60,6.10623e-016},{72,12},{98,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(O2.faces[transCart[Axis.z], Side.p], zPositive.O2) annotation (Line(
        points={{60,6.10623e-016},{48,-12},{-98,-12}},
        color={127,127,127},
        smooth=Smooth.None));
    annotation (
      Documentation(info="<html>
<p>Please see the documentation of the <a href=\"modelica://FCSys.Phases.Partial\">Partial</a> model.</p></html>"),

      Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-120,-60},{
              120,60}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
              100}}), graphics));

  end Gas;

  model Graphite "Graphite phase"
    import assert = FCSys.Utilities.assertEval;
    import Modelica.Math.BooleanVectors.countTrue;

    extends Phase(
      final reduceTrans=false,
      final reduceThermal=true,
      final n_spec=countTrue({'inclC+','incle-'}));

    // Conditionally include species.
    parameter Boolean 'inclC+'=false "Include C+" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Carbon plus (C<sup>+</sup>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.'C+'.Graphite.Fixed 'C+'(final n_trans) if
      'inclC+' constrainedby FCSys.Species.Solid(
      n_trans=n_trans,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "C+ model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>C<sup>+</sup> model</html>",
        enable='inclC+'),
      Placement(transformation(extent={{-30,-30},{-10,-10}})));

    parameter Boolean 'incle-'=false "Include e-" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Electrons (e<sup>-</sup>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.'e-'.Graphite.Fixed 'e-'(final n_trans) if
      'incle-' constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "e- model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>'e-' model</html>",
        enable='incle-'),
      Placement(transformation(extent={{10,-30},{30,-10}})));

    Connectors.FaceBus xNegative if inclTrans[Axis.x]
      "Negative face along the x axis" annotation (Placement(transformation(
            extent={{-80,-30},{-60,-10}}),iconTransformation(extent={{-90,-10},
              {-70,10}})));
    Connectors.FaceBus yNegative if inclTrans[Axis.y]
      "Negative face along the y axis" annotation (Placement(transformation(
            extent={{-56,-54},{-36,-34}}),iconTransformation(extent={{-10,-94},
              {10,-74}})));
    Connectors.FaceBus zNegative if inclTrans[Axis.z]
      "Negative face along the z axis" annotation (Placement(transformation(
            extent={{48,-18},{68,2}}),iconTransformation(extent={{40,40},{60,60}})));
    Connectors.FaceBus xPositive if inclTrans[Axis.x]
      "Positive face along the x axis" annotation (Placement(transformation(
            extent={{60,-30},{80,-10}}),iconTransformation(extent={{70,-10},{90,
              10}})));
    Connectors.FaceBus yPositive if inclTrans[Axis.y]
      "Positive face along the y axis" annotation (Placement(transformation(
            extent={{36,-6},{56,14}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.FaceBus zPositive if inclTrans[Axis.z]
      "Positive face along the z axis" annotation (Placement(transformation(
            extent={{-68,-42},{-48,-22}}),iconTransformation(extent={{-90,-90},
              {-70,-70}})));
    Connectors.Inter inter[n_inter](each final n_trans=n_trans) if n_spec > 0
      "Connector to exchange momentum and energy with other phases" annotation
      (Placement(transformation(extent={{60,-66},{80,-46}}), iconTransformation(
            extent={{-60,60},{-40,40}})));

    Connectors.Amagat amagat(final V=-V) if n_spec > 0
      "Connector for additivity of volume" annotation (Placement(transformation(
            extent={{-80,4},{-60,24}}), iconTransformation(extent={{70,-90},{90,
              -70}})));

    // Aliases
    Q.Velocity phi[n_trans](each stateSelect=StateSelect.prefer) = direct.translational.phi
      if n_spec > 1 and reduceTrans "Velocity";
    Q.Temperature T(each stateSelect=StateSelect.prefer) = direct.thermal.T if
      n_spec > 1 and reduceThermal "Temperature";
    // These make the selected states more readable.

    Connectors.Chemical 'conne-'(final n_trans=n_trans) if 'incle-'
      "Chemical connector for e-" annotation (Placement(transformation(extent={
              {6,20},{26,40}}), iconTransformation(extent={{-10,-50},{10,-30}})));
    Connectors.Electrostatic electrostatic if 'incle-' and n_spec > 0
      "Interface with the dielectric" annotation (Placement(transformation(
            extent={{60,4},{80,24}}), iconTransformation(extent={{90,-50},{110,
              -30}})));
    // Note:  ('incle-' and n_spec > 0) can be logically reduced to 'incle-',
    // but n_spec is included so that the connector always appears in the icon
    // layer in Dymola 2014.

  protected
    Connectors.DirectNode direct(
      final n_trans=n_trans,
      final inclTrans=reduceTrans,
      final inclThermal=reduceThermal) if n_spec > 0 and (reduceTrans or
      reduceThermal)
      "Connector to directly couple velocities and temperatures within the phase"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={
              0,-68}), iconTransformation(extent={{-10,-10},{10,10}}, origin={
              26,-52})));

    Conditions.Adapters.AmagatDalton amagatDalton if n_spec > 0
      "Adapter between additivity of volume and additivity of pressure"
      annotation (Placement(transformation(extent={{-60,4},{-40,24}})));

  equation
    // Chemical exchange
    connect('e-'.chemical, 'conne-') annotation (Line(
        points={{16,-11},{16,30}},
        color={221,23,47},
        smooth=Smooth.None));

    // Electrical storage
    connect(electrostatic, 'e-'.electrostatic) annotation (Line(
        points={{70,14},{13,14},{13,-13}},
        color={255,195,38},
        smooth=Smooth.None));

    // Inert exchange
    connect('C+'.inter, inter) annotation (Line(
        points={{-11,-23.8},{-9,-26},{-9,-56},{70,-56}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('e-'.inter, inter) annotation (Line(
        points={{29,-23.8},{31,-26},{31,-56},{70,-56}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('C+'.direct.translational, direct.translational) annotation (Line(
        points={{-16.2,-29},{-16.2,-68},{0,-68}},
        color={47,107,251},
        smooth=Smooth.None));
    connect('C+'.direct.thermal, direct.thermal) annotation (Line(
        points={{-16.2,-29},{-16.2,-68},{0,-68}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('e-'.direct.translational, direct.translational) annotation (Line(
        points={{23.8,-29},{23.8,-68},{0,-68}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('e-'.direct.thermal, direct.thermal) annotation (Line(
        points={{23.8,-29},{23.8,-68},{0,-68}},
        color={38,196,52},
        smooth=Smooth.None));

    // Mixing
    connect(amagatDalton.amagat, amagat) annotation (Line(
        points={{-54,14},{-70,14}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(amagatDalton.dalton, 'C+'.dalton) annotation (Line(
        points={{-46,14},{-29,14},{-29,-16}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(amagatDalton.dalton, 'e-'.dalton) annotation (Line(
        points={{-46,14},{11,14},{11,-16}},
        color={47,107,251},
        smooth=Smooth.None));

    // Transport
    // ---------
    // C+
    connect('C+'.faces[transCart[Axis.x], Side.n], xNegative.'C+') annotation (
        Line(
        points={{-20,-20},{-70,-20}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('C+'.faces[transCart[Axis.x], Side.p], xPositive.'C+') annotation (
        Line(
        points={{-20,-20},{70,-20}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('C+'.faces[transCart[Axis.y], Side.n], yNegative.'C+') annotation (
        Line(
        points={{-20,-20},{-20,-44},{-46,-44}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('C+'.faces[transCart[Axis.y], Side.p], yPositive.'C+') annotation (
        Line(
        points={{-20,-20},{-20,4},{46,4}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('C+'.faces[transCart[Axis.z], Side.n], zNegative.'C+') annotation (
        Line(
        points={{-20,-20},{-8,-8},{58,-8}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('C+'.faces[transCart[Axis.z], Side.p], zPositive.'C+') annotation (
        Line(
        points={{-20,-20},{-32,-32},{-58,-32}},
        color={127,127,127},
        smooth=Smooth.None));

    // e-
    connect('e-'.faces[transCart[Axis.x], Side.n], xNegative.'e-') annotation (
        Line(
        points={{20,-20},{-70,-20}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.faces[transCart[Axis.x], Side.p], xPositive.'e-') annotation (
        Line(
        points={{20,-20},{70,-20}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.faces[transCart[Axis.y], Side.n], yNegative.'e-') annotation (
        Line(
        points={{20,-20},{20,-44},{-46,-44}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.faces[transCart[Axis.y], Side.p], yPositive.'e-') annotation (
        Line(
        points={{20,-20},{20,4},{46,4}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.faces[transCart[Axis.z], Side.n], zNegative.'e-') annotation (
        Line(
        points={{20,-20},{32,-8},{58,-8}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('e-'.faces[transCart[Axis.z], Side.p], zPositive.'e-') annotation (
        Line(
        points={{20,-20},{8,-32},{-58,-32}},
        color={127,127,127},
        smooth=Smooth.None));

    annotation (
      Documentation(info="<html>
    <p>See <a href=\"modelica://FCSys.Species.'e-'.Graphite.Fixed\">Species.'e-'.Graphite.Fixed</a> for assumptions.
    For more information, see the
 <a href=\"modelica://FCSys.Phases.Partial\">Partial</a> model.</p></html>"),
      Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-80,-80},{80,
              40}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
              100}}), graphics));
  end Graphite;

  model Ionomer "Ionomer phase"
    import Modelica.Math.BooleanVectors.countTrue;

    extends Phase(final n_spec=countTrue({'inclSO3-','inclH+',inclH2O}));

    parameter Q.NumberAbsolute k_common=1
      "Coupling factor for exchange among all species within the phase"
      annotation (Dialog(group="Geometry", __Dymola_label=
            "<html><i>k</i><sub>common</sub></html>"));

    // Conditionally include species.
    parameter Boolean 'inclSO3-'=false
      "Include C19HF37O5S- (abbreviated as SO3-)" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label=
            "<html>Nafion sulfonate minus (abbreviated as SO<sub>3</sub><sup>-</sup>)</html>",

        __Dymola_joinNext=true));

    replaceable FCSys.Species.'SO3-'.Ionomer.Fixed 'SO3-'(final n_trans) if
      'inclSO3-' constrainedby FCSys.Species.Solid(
      n_trans=n_trans,
      n_intra=1,
      k_intra={k_common},
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "SO3- model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>SO<sub>3</sub><sup>-</sup> model</html>",
        enable='inclSO3-'),
      Placement(transformation(extent={{30,-20},{50,0}})));

    parameter Boolean 'inclH+'=false "Include H+" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Protons (H<sup>+</sup>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.'H+'.Ionomer.Fixed 'H+'(final n_trans) if
      'inclH+' constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      n_intra=1,
      k_intra={k_common},
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "H+ model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>H<sup>+</sup> model</html>",
        enable='inclH+'),
      Placement(transformation(extent={{-50,-20},{-30,0}})));

    parameter Boolean inclH2O=false "Include H2O" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Water (H<sub>2</sub>O)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.H2O.Ionomer.Fixed H2O(final n_trans, n_intra=1)
      if inclH2O constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      n_intra=1,
      k_intra={k_common},
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "H2O model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>H<sub>2</sub>O model</html>",
        enable=inclH2O),
      Placement(transformation(extent={{-10,-20},{10,0}})));

    Connectors.FaceBus xNegative if inclTrans[Axis.x]
      "Negative face along the x axis" annotation (Placement(transformation(
            extent={{-100,-20},{-80,0}}), iconTransformation(extent={{-90,-10},
              {-70,10}})));
    Connectors.FaceBus yNegative if inclTrans[Axis.y]
      "Negative face along the y axis" annotation (Placement(transformation(
            extent={{-76,-44},{-56,-24}}),iconTransformation(extent={{-10,-94},
              {10,-74}})));
    Connectors.FaceBus zNegative if inclTrans[Axis.z]
      "Negative face along the z axis" annotation (Placement(transformation(
            extent={{68,-8},{88,12}}),iconTransformation(extent={{40,40},{60,60}})));
    Connectors.FaceBus xPositive if inclTrans[Axis.x]
      "Positive face along the x axis" annotation (Placement(transformation(
            extent={{80,-20},{100,0}}), iconTransformation(extent={{70,-10},{90,
              10}})));
    Connectors.FaceBus yPositive if inclTrans[Axis.y]
      "Positive face along the y axis" annotation (Placement(transformation(
            extent={{56,4},{76,24}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.FaceBus zPositive if inclTrans[Axis.z]
      "Positive face along the z axis" annotation (Placement(transformation(
            extent={{-88,-32},{-68,-12}}),iconTransformation(extent={{-90,-90},
              {-70,-70}})));
    Connectors.Amagat amagat(final V=-V) if n_spec > 0
      "Connector for additivity of volume" annotation (Placement(transformation(
            extent={{-100,16},{-80,36}}), iconTransformation(extent={{70,-90},{
              90,-70}})));
    Connectors.Inter inter[n_inter](each final n_trans=n_trans) if n_spec > 0
      "Connector to exchange momentum and energy with other phases" annotation
      (Placement(transformation(extent={{80,-56},{100,-36}}),
          iconTransformation(extent={{-60,60},{-40,40}})));
    Connectors.Chemical connH2O(final n_trans=n_trans) if inclH2O
      "Chemical connector for H2O" annotation (Placement(transformation(extent=
              {{-14,30},{6,50}}), iconTransformation(extent={{10,-50},{30,-30}})));
    Connectors.Chemical 'connH+'(final n_trans=n_trans) if 'inclH+'
      "Chemical connector for H+" annotation (Placement(transformation(extent={
              {-54,30},{-34,50}}), iconTransformation(extent={{-30,-50},{-10,-30}})));
    Connectors.Electrostatic electrostatic if 'inclH+' and n_spec > 0
      "Interface with the dielectric" annotation (Placement(transformation(
            extent={{-100,4},{-80,24}}), iconTransformation(extent={{-110,-50},
              {-90,-30}})));
    // Note:  ('inclH+' and n_spec > 0) can be logically reduced to 'inclH+',
    // but n_spec is included so that the connector always appears in the icon
    // layer in Dymola 2014.

  protected
    Conditions.Adapters.AmagatDalton amagatDalton if n_spec > 0
      "Adapter between additivity of volume and additivity of pressure"
      annotation (Placement(transformation(extent={{-80,16},{-60,36}})));
    FCSys.Connectors.DirectNode direct(
      final n_trans=n_trans,
      final inclTrans=reduceTrans,
      final inclThermal=reduceThermal) if n_spec > 0 and (reduceTrans or
      reduceThermal)
      "Connector to directly couple velocities and temperatures within the phase"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={
              24,-70}), iconTransformation(extent={{-10,-10},{10,10}}, origin={
              26,-52})));

    Connectors.InertNode intra "Connection node for exchange within the phase"
      annotation (Placement(transformation(extent={{18,-68},{38,-48}}),
          iconTransformation(extent={{48,-76},{68,-56}})));

  equation
    // Chemical exchange
    connect(H2O.chemical, connH2O) annotation (Line(
        points={{-4,-1},{-4,40}},
        color={221,23,47},
        smooth=Smooth.None));
    connect('H+'.chemical, 'connH+') annotation (Line(
        points={{-44,-1},{-44,40}},
        color={221,23,47},
        smooth=Smooth.None));

    // Electrical storage
    connect(electrostatic, 'H+'.electrostatic) annotation (Line(
        points={{-90,14},{-47,14},{-47,-3}},
        color={255,195,38},
        smooth=Smooth.None));

    // Inert exchange
    connect('SO3-'.inter, inter) annotation (Line(
        points={{49,-13.8},{52,-13.8},{52,-46},{90,-46}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('SO3-'.direct.translational, direct.translational) annotation (Line(
        points={{43.8,-19},{43.8,-70},{24,-70}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('SO3-'.direct.thermal, direct.thermal) annotation (Line(
        points={{43.8,-19},{43.8,-70},{24,-70}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('H+'.inter, inter) annotation (Line(
        points={{-31,-13.8},{-28,-13.8},{-28,-46},{90,-46}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('H+'.direct.translational, direct.translational) annotation (Line(
        points={{-36.2,-19},{-36.2,-70},{24,-70}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('H+'.direct.thermal, direct.thermal) annotation (Line(
        points={{-36.2,-19},{-36.2,-70},{24,-70}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2O.inter, inter) annotation (Line(
        points={{9,-13.8},{12,-13.8},{12,-46},{90,-46}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2O.direct.translational, direct.translational) annotation (Line(
        points={{3.8,-19},{3.8,-70},{24,-70}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2O.direct.thermal, direct.thermal) annotation (Line(
        points={{3.8,-19},{3.8,-70},{24,-70}},
        color={38,196,52},
        smooth=Smooth.None));

    connect('SO3-'.intra[1], intra.exchange) annotation (Line(
        points={{47,-17},{47,-58},{28,-58}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('H+'.intra[1], intra.exchange) annotation (Line(
        points={{-33,-17},{-33,-58},{28,-58}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2O.intra[1], intra.exchange) annotation (Line(
        points={{7,-17},{7,-58},{28,-58}},
        color={38,196,52},
        smooth=Smooth.None));

    // Mixing
    connect('SO3-'.dalton, amagatDalton.dalton) annotation (Line(
        points={{31,-6},{31,26},{-66,26}},
        color={47,107,251},
        smooth=Smooth.None));
    connect('H+'.dalton, amagatDalton.dalton) annotation (Line(
        points={{-49,-6},{-49,26},{-66,26}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(H2O.dalton, amagatDalton.dalton) annotation (Line(
        points={{-9,-6},{-9,26},{-66,26}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(amagatDalton.amagat, amagat) annotation (Line(
        points={{-74,26},{-90,26}},
        color={47,107,251},
        smooth=Smooth.None));

    // Transport
    // ---------
    // SO3-
    connect('SO3-'.faces[transCart[Axis.x], Side.n], xNegative.'SO3-')
      annotation (Line(
        points={{40,-10},{-90,-10}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('SO3-'.faces[transCart[Axis.x], Side.p], xPositive.'SO3-')
      annotation (Line(
        points={{40,-10},{90,-10}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('SO3-'.faces[transCart[Axis.y], Side.n], yNegative.'SO3-')
      annotation (Line(
        points={{40,-10},{40,-34},{-66,-34}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('SO3-'.faces[transCart[Axis.y], Side.p], yPositive.'SO3-')
      annotation (Line(
        points={{40,-10},{40,14},{66,14}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('SO3-'.faces[transCart[Axis.z], Side.n], zNegative.'SO3-')
      annotation (Line(
        points={{40,-10},{52,2},{78,2}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('SO3-'.faces[transCart[Axis.z], Side.p], zPositive.'SO3-')
      annotation (Line(
        points={{40,-10},{28,-22},{-78,-22}},
        color={127,127,127},
        smooth=Smooth.None));
    // 'H+'
    connect('H+'.faces[transCart[Axis.x], Side.n], xNegative.'H+') annotation (
        Line(
        points={{-40,-10},{-90,-10}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.faces[transCart[Axis.x], Side.p], xPositive.'H+') annotation (
        Line(
        points={{-40,-10},{90,-10}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.faces[transCart[Axis.y], Side.n], yNegative.'H+') annotation (
        Line(
        points={{-40,-10},{-40,-34},{-66,-34}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.faces[transCart[Axis.y], Side.p], yPositive.'H+') annotation (
        Line(
        points={{-40,-10},{-40,14},{66,14}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.faces[transCart[Axis.z], Side.n], zNegative.'H+') annotation (
        Line(
        points={{-40,-10},{-28,2},{78,2}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('H+'.faces[transCart[Axis.z], Side.p], zPositive.'H+') annotation (
        Line(
        points={{-40,-10},{-52,-22},{-78,-22}},
        color={127,127,127},
        smooth=Smooth.None));
    // H2O
    connect(H2O.faces[transCart[Axis.x], Side.n], xNegative.H2O) annotation (
        Line(
        points={{0,-10},{-90,-10}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[transCart[Axis.x], Side.p], xPositive.H2O) annotation (
        Line(
        points={{0,-10},{90,-10}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[transCart[Axis.y], Side.n], yNegative.H2O) annotation (
        Line(
        points={{0,-10},{0,-34},{-66,-34}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[transCart[Axis.y], Side.p], yPositive.H2O) annotation (
        Line(
        points={{0,-10},{0,14},{66,14}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[transCart[Axis.z], Side.n], zNegative.H2O) annotation (
        Line(
        points={{0,-10},{12,2},{78,2}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(H2O.faces[transCart[Axis.z], Side.p], zPositive.H2O) annotation (
        Line(
        points={{0,-10},{0,-22},{-78,-22}},
        color={127,127,127},
        smooth=Smooth.None));

    annotation (
      Documentation(info="<html><p>Assumptions:<ol>
    <li>The water in the ionomer does not participate in the reaction (only the water vapor does).</li>
    </ol</p>
    
    <p>See <a href=\"modelica://FCSys.Species.'H+'.Ionomer.Fixed\">Species.'H+'.Ionomer.Fixed</a> 
    for additional assumptions.
    For more information, see the
 <a href=\"modelica://FCSys.Phases.Partial\">Partial</a> model.</p></html>"),
      Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-80},{
              100,60}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
              100}}), graphics));
  end Ionomer;

  model Liquid "Liquid phase"
    extends Phase(final n_spec=if inclH2O then 1 else 0);

    // Conditionally include species.
    parameter Boolean inclH2O=false "Include H2O" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Water (H<sub>2</sub>O)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.H2O.Liquid.Fixed H2O(final n_trans) if inclH2O
      constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "H2O model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>H<sub>2</sub>O model</html>",
        enable=inclH2O),
      Placement(transformation(extent={{-20,-20},{0,0}})));

    Connectors.FaceBus xNegative if inclTrans[Axis.x]
      "Negative face along the x axis" annotation (Placement(transformation(
            extent={{-60,-20},{-40,0}}), iconTransformation(extent={{-90,-10},{
              -70,10}})));
    Connectors.FaceBus yNegative if inclTrans[Axis.y]
      "Negative face along the y axis" annotation (Placement(transformation(
            extent={{-20,-60},{0,-40}}), iconTransformation(extent={{-10,-94},{
              10,-74}})));
    Connectors.FaceBus zNegative if inclTrans[Axis.z]
      "Negative face along the z axis" annotation (Placement(transformation(
            extent={{0,0},{20,20}}), iconTransformation(extent={{40,40},{60,60}})));
    Connectors.FaceBus xPositive if inclTrans[Axis.x]
      "Positive face along the x axis" annotation (Placement(transformation(
            extent={{20,-20},{40,0}}), iconTransformation(extent={{70,-10},{90,
              10}})));
    Connectors.FaceBus yPositive if inclTrans[Axis.y]
      "Positive face along the y axis" annotation (Placement(transformation(
            extent={{-20,20},{0,40}}),iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.FaceBus zPositive if inclTrans[Axis.z]
      "Positive face along the z axis" annotation (Placement(transformation(
            extent={{-40,-40},{-20,-20}}), iconTransformation(extent={{-90,-90},
              {-70,-70}})));
    Connectors.Inter inter[n_inter](each final n_trans=n_trans) if n_spec > 0
      "Connector to exchange momentum and energy with other phases" annotation
      (Placement(transformation(extent={{20,-40},{40,-20}}), iconTransformation(
            extent={{-60,60},{-40,40}})));
    Connectors.Amagat amagat(final V=-V) if n_spec > 0
      "Connector for additivity of volume" annotation (Placement(transformation(
            extent={{-60,0},{-40,20}}), iconTransformation(extent={{70,-90},{90,
              -70}})));
    Connectors.Chemical connH2O(final n_trans=n_trans) if inclH2O
      "Chemical connector for H2O" annotation (Placement(transformation(extent=
              {{-24,40},{-4,60}}), iconTransformation(extent={{-10,-50},{10,-30}})));

  protected
    Conditions.Adapters.AmagatDalton amagatDalton if n_spec > 0
      "Adapter between additivity of volume and additivity of pressure"
      annotation (Placement(transformation(extent={{-40,0},{-20,20}})));

  equation
    // Chemical exchange
    connect(H2O.chemical, connH2O) annotation (Line(
        points={{-14,-1},{-14,50}},
        color={221,23,47},
        smooth=Smooth.None));

    // Inert exchange
    connect(H2O.inter, inter) annotation (Line(
        points={{-1,-13.8},{1,-16},{1,-30},{30,-30}},
        color={38,196,52},
        smooth=Smooth.None));

    // Mixing
    connect(amagatDalton.amagat, amagat) annotation (Line(
        points={{-34,10},{-50,10}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(amagatDalton.dalton, H2O.dalton) annotation (Line(
        points={{-26,10},{-19,10},{-19,-6}},
        color={47,107,251},
        smooth=Smooth.None));

    // Transport
    // ---------
    // H2O
    connect(H2O.faces[transCart[Axis.x], Side.n], xNegative.H2O) annotation (
        Line(
        points={{-10,-10},{-50,-10}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[transCart[Axis.x], Side.p], xPositive.H2O) annotation (
        Line(
        points={{-10,-10},{30,-10}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[transCart[Axis.y], Side.n], yNegative.H2O) annotation (
        Line(
        points={{-10,-10},{-10,-50}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[transCart[Axis.y], Side.p], yPositive.H2O) annotation (
        Line(
        points={{-10,-10},{-10,30}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[transCart[Axis.z], Side.n], zNegative.H2O) annotation (
        Line(
        points={{-10,-10},{10,10}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(H2O.faces[transCart[Axis.z], Side.p], zPositive.H2O) annotation (
        Line(
        points={{-10,-10},{-30,-30}},
        color={127,127,127},
        smooth=Smooth.None));

    annotation (
      Documentation(info="<html>
    <p>Please see the documentation of the
 <a href=\"modelica://FCSys.Phases.Partial\">Partial</a> model.</p></html>"),
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-60,-60},{40,
              60}}), graphics));
  end Liquid;

protected
  partial model Phase "Base model for a phase"
    import Modelica.Math.BooleanVectors.index;
    // extends FCSys.Icons.Names.Middle;

    parameter Integer n_spec(start=0) "Number of species"
      annotation (HideResult=true);
    inner parameter Integer n_inter=0
      "Number of independent couplings with other phases"
      annotation (Dialog(connectorSizing=true),HideResult=n_spec == 0);

    // Geometric parameters
    inner parameter Q.NumberAbsolute k[Axis](
      each min=Modelica.Constants.small,
      each final nominal=1) = {1,1,1} if n_spec > 0
      "Scaling factor for diffusive transport" annotation (Dialog(group=
            "Geometry", __Dymola_label="<html><b><i>k</i></b></html>"));
    inner parameter Q.NumberAbsolute k_inter[n_inter]=ones(n_inter) if n_spec
       > 0 "Coupling factors with other phases" annotation (Dialog(group=
            "Geometry", __Dymola_label="<html><i>k</i><sub>inter</sub></html>"));
    parameter Integer n_trans=1 "Number of transport axes"
      annotation (HideResult=true);
    // This cannot be an inner/outer parameter in Dymola 2014.

    // Assumptions
    inner parameter Boolean reduceTrans=false "Same velocity for all species"
      annotation (
      HideResult=true,
      Dialog(tab="Assumptions", enable=n_spec > 1),
      choices(__Dymola_checkBox=true));
    inner parameter Boolean reduceThermal=false
      "Same temperature for all species" annotation (
      HideResult=true,
      Dialog(tab="Assumptions", enable=n_spec > 1),
      choices(__Dymola_checkBox=true));

    inner Q.Volume V if n_spec > 0 "Volume";

    // Auxiliary variables (for analysis)
    output Q.Number epsilon=V/product(L) if n_spec > 0 and environment.analysis
      "Volumetric fill fraction";

  protected
    outer parameter Q.Length L[Axis] if n_spec > 0 "Length of the subregion"
      annotation (missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Q.Area A[Axis] if n_spec > 0
      "Cross-sectional area of the subregion" annotation (missingInnerMessage="This model should be used within a subregion model.
");
    final inner Q.Length Lprime[n_trans]=k[cartTrans]*V ./ L[cartTrans] .^ 2
      if n_spec > 0 "Effective area divided by transport length";
    outer parameter Integer cartTrans[:]
      "Cartesian-axis indices of the components of translational momentum"
      annotation (missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Integer transCart[:]
      "Face-pair indices of the Cartesian axes" annotation (missingInnerMessage
        ="This model should be used within a subregion model.
");
    outer parameter Boolean inclTrans[Axis]
      "true, if each pairs of faces is included" annotation (
        missingInnerMessage="This model should be used within a subregion model.
");

    outer Conditions.Environment environment "Environmental conditions";
    // This component is conditional to prevent a mathematical singularity
    // when two or more empty phases (without any species included) are
    // connected.

    annotation (
      defaultComponentPrefixes="replaceable",
      defaultComponentName="phase",
      Documentation(info="<html><p>The scaling factor for diffusive transport (<b><i>k</i></b>) is a vector which directly affects 
    the rates of diffusion of material, transverse translational momentum, and energy of all of the species
    within the phase.  It can be used to introduce minor head loss or the effects of
    porosity or tortousity.  These effects may be anisotropic. Using the
    Bruggeman correction factor [<a href=\"modelica://FCSys.UsersGuide.References\">Weber2004</a>, p. 4696],
    the scaling factor for diffusive transport (<b><i>k</i></b>) within a phase should be set to &epsilon;<sup>3/2</sup>
    along each axis, where &epsilon; is the volumetric filling ratio, or the ratio of the volume of the phase to the total volume of the subregion.</p>


</html>"),
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics={
          Ellipse(
            extent={{-40,100},{40,20}},
            lineColor={127,127,127},
            startAngle=30,
            endAngle=149,
            pattern=LinePattern.Dash,
            fillPattern=FillPattern.Solid,
            fillColor={225,225,225}),
          Ellipse(
            extent={{20,-4},{100,-84}},
            lineColor={127,127,127},
            startAngle=270,
            endAngle=390,
            pattern=LinePattern.Dash,
            fillPattern=FillPattern.Solid,
            fillColor={225,225,225}),
          Ellipse(
            extent={{-100,-4},{-20,-84}},
            lineColor={127,127,127},
            startAngle=149,
            endAngle=270,
            pattern=LinePattern.Dash,
            fillPattern=FillPattern.Solid,
            fillColor={225,225,225}),
          Polygon(
            points={{60,-84},{-60,-84},{-94.5,-24},{-34.5,80},{34.5,80},{94.5,-24},
                {60,-84}},
            pattern=LinePattern.None,
            fillPattern=FillPattern.Sphere,
            smooth=Smooth.None,
            fillColor={225,225,225},
            lineColor={0,0,0}),
          Line(
            points={{-60,-84.1},{60,-84.1}},
            color={127,127,127},
            pattern=LinePattern.Dash,
            smooth=Smooth.None),
          Line(
            points={{34.5,80},{94.5,-24}},
            color={127,127,127},
            pattern=LinePattern.Dash,
            smooth=Smooth.None),
          Line(
            points={{-34.5,80},{-94.5,-24}},
            color={127,127,127},
            pattern=LinePattern.Dash,
            smooth=Smooth.None),
          Text(
            extent={{-100,-20},{100,20}},
            textString="%name",
            lineColor={0,0,0})}),
      Diagram(graphics));
  end Phase;

public
  model Dielectric "Dielectric gap"
    extends FCSys.Icons.Names.Top2;

    parameter Q.Length L=U.um "Length of the dielectric (not of the region)"
      annotation (Dialog(group="Geometry", __Dymola_label=
            "<html><i>L</i></html>"));
    parameter Q.Area A=U.cm^2
      "Cross-sectional area of the dielectric (not of the region)" annotation (
        Dialog(group="Geometry", __Dymola_label="<html><i>A</i></html>"));
    parameter Q.Permittivity epsilon=U.epsilon_0 "Permittivity"
      annotation (Dialog(__Dymola_label="<html>&epsilon;</html>"));
    final parameter Q.Capacitance C=epsilon*A/L "Capacitance";

    Q.Amount Z "Amount of charge shifted in the positive direction";
    Q.Potential w "Electrical potential";

    Connectors.Electrostatic negative "Reference electrical connector"
      annotation (Placement(transformation(extent={{-30,-10},{-10,10}}),
          iconTransformation(extent={{-70,-10},{-50,10}})));
    Connectors.Electrostatic positive "Positive electrical connector"
      annotation (Placement(transformation(extent={{10,-10},{30,10}}),
          iconTransformation(extent={{50,-10},{70,10}})));
    Connectors.Amagat amagat
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

  equation
    // Aliases
    w = positive.w - negative.w;
    Z = negative.Z;

    // Connector equations
    Z = C*w "Electrical capacitance";
    0 = negative.Z + positive.Z "Charge neutrality";
    0 = amagat.V + A*L "Conservation of volume";

    annotation (Documentation(info="<html><p>Assumptions:<ol>
<li>No material
<li>No momentum (follows from #1)</li>
<li>No heat capacity (follows from #1)</li>
<li>The charges exist on parallel planes (used to calculate capacitance).</li> 
</ol></p></html>"), Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics={
          Line(
            points={{-20,30},{-20,-30}},
            color={255,195,38},
            smooth=Smooth.None),
          Line(
            points={{20,30},{20,-30}},
            color={255,195,38},
            smooth=Smooth.None),
          Line(
            points={{-20,0},{-50,0}},
            color={255,195,38},
            smooth=Smooth.None),
          Line(
            points={{50,0},{20,0}},
            color={255,195,38},
            smooth=Smooth.None)}));
  end Dielectric;
  annotation (Documentation(info="
<html><p>The graphite, ionomer, and
liquid phases can only be used with a compressible phase (gas).</p>

  <p><b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
Copyright 2007&ndash;2013, <a href=\"http://www.gtrc.gatech.edu/\">Georgia Tech Research Corporation</a>.</p>

<p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.License\">
FCSys.UsersGuide.License</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p>
</html>"));
end Phases;
