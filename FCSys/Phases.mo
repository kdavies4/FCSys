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

    equation
      // Aliases
      electrical.Z = z*N;

      // Equal potential
      electrical.w = pin.v*U.V;

      // Conservation
      der(N)/U.s = (pin.i*U.A)/z "Material";

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
    end Plate;
  end Examples;

  model Gas "Gas phase"
    import Modelica.Math.BooleanVectors.countTrue;

    extends Partial(final n_spec=countTrue({inclH2,inclH2O,inclN2,inclO2}));

    // Conditionally include species.
    parameter Boolean inclH2=false "Include H2" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Hydrogen (H<sub>2</sub>)</html>",
        __Dymola_joinNext=true));
    replaceable FCSys.Species.H2.Gas.Fixed H2(final n_trans,final n_inter) if
      inclH2 constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      n_inter=n_inter,
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

    replaceable FCSys.Species.H2O.Gas.Fixed H2O(final n_trans,final n_inter)
      if inclH2O constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      n_inter=n_inter,
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

    replaceable FCSys.Species.N2.Gas.Fixed N2(final n_trans,final n_inter) if
      inclN2 constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      n_inter=n_inter,
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

    replaceable FCSys.Species.O2.Gas.Fixed O2(final n_trans,final n_inter) if
      inclO2 constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      n_inter=n_inter,
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

    Connectors.Electrochemical connH2(final n_trans=n_trans) if inclH2
      "Chemical connector for H2" annotation (Placement(transformation(extent={
              {-74,40},{-54,60}}), iconTransformation(extent={{-50,-50},{-30,-30}})));
    Connectors.Electrochemical connH2O(final n_trans=n_trans) if inclH2O
      "Chemical connector for H2O" annotation (Placement(transformation(extent=
              {{-34,40},{-14,60}}), iconTransformation(extent={{-10,-50},{10,-30}})));
    Connectors.Electrochemical connO2(final n_trans=n_trans) if inclO2
      "Chemical connector for O2" annotation (Placement(transformation(extent={
              {46,40},{66,60}}), iconTransformation(extent={{30,-50},{50,-30}})));

    // Aliases
    Q.Velocity phi[n_trans](each stateSelect=StateSelect.prefer) = direct.translational.phi
      if n_spec > 0 and reduceTrans "Velocity";
    Q.Temperature T(stateSelect=StateSelect.prefer) = direct.thermal.T if
      n_spec > 0 and reduceThermal "Temperature";
    // These make the selected states more readable.

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
    connect(O2.electrochemical, connO2) annotation (Line(
        points={{56,9},{56,50}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(H2.electrochemical, connH2) annotation (Line(
        points={{-64,9},{-64,50}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(H2O.electrochemical, connH2O) annotation (Line(
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

    extends Partial(
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

    replaceable FCSys.Species.'C+'.Graphite.Fixed 'C+'(final n_trans, final
        n_inter) if 'inclC+' constrainedby FCSys.Species.Solid(
      n_trans=n_trans,
      n_intra=1,
      n_inter=n_inter,
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
      Placement(transformation(extent={{-30,-10},{-10,10}})));

    parameter Boolean 'incle-'=true "Include e-" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Electrons (e<sup>-</sup>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.'e-'.Graphite.Fixed 'e-'(final n_trans,final
        n_inter) if 'incle-' constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      n_intra=1,
      n_inter=0,
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
      Placement(transformation(extent={{10,-10},{30,10}})));

    Connectors.FaceBus xNegative if inclTrans[Axis.x]
      "Negative face along the x axis" annotation (Placement(transformation(
            extent={{-80,-10},{-60,10}}), iconTransformation(extent={{-90,-10},
              {-70,10}})));
    Connectors.FaceBus yNegative if inclTrans[Axis.y]
      "Negative face along the y axis" annotation (Placement(transformation(
            extent={{-56,-34},{-36,-14}}),iconTransformation(extent={{-10,-94},
              {10,-74}})));
    Connectors.FaceBus zNegative if inclTrans[Axis.z]
      "Negative face along the z axis" annotation (Placement(transformation(
            extent={{48,2},{68,22}}), iconTransformation(extent={{40,40},{60,60}})));
    Connectors.FaceBus xPositive if inclTrans[Axis.x]
      "Positive face along the x axis" annotation (Placement(transformation(
            extent={{60,-10},{80,10}}), iconTransformation(extent={{70,-10},{90,
              10}})));
    Connectors.FaceBus yPositive if inclTrans[Axis.y]
      "Positive face along the y axis" annotation (Placement(transformation(
            extent={{36,14},{56,34}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.FaceBus zPositive if inclTrans[Axis.z]
      "Positive face along the z axis" annotation (Placement(transformation(
            extent={{-68,-22},{-48,-2}}), iconTransformation(extent={{-90,-90},
              {-70,-70}})));
    Connectors.Inter inter[n_inter](each final n_trans=n_trans) if n_spec > 0
      "Connector to exchange momentum and energy with other phases" annotation
      (Placement(transformation(extent={{60,-46},{80,-26}}), iconTransformation(
            extent={{-60,60},{-40,40}})));

    Connectors.Amagat amagat(final V=-V) if n_spec > 0
      "Connector for additivity of volume" annotation (Placement(transformation(
            extent={{-80,24},{-60,44}}), iconTransformation(extent={{70,-90},{
              90,-70}})));

    Connectors.Electrochemical 'conne-'(final n_trans=n_trans) if 'incle-'
      "Electrochemical connector for e-" annotation (Placement(transformation(
            extent={{6,40},{26,60}}), iconTransformation(extent={{-10,-50},{10,
              -30}})));
    Connectors.Electrostatic electrostatic if 'incle-' and n_spec > 0
      "Interface with the dielectric" annotation (Placement(transformation(
            extent={{60,24},{80,44}}), iconTransformation(extent={{90,-50},{110,
              -30}})));
    // Note:  ('incle-' and n_spec > 0) is logically equivalent to 'incle-',
    // but n_spec is included so that the connector always appears in the icon
    // layer in Dymola 2014.

    // Aliases
    Q.Velocity phi[n_trans](each stateSelect=StateSelect.prefer) = direct.translational.phi
      if n_spec > 0 and reduceTrans "Velocity";
    Q.Temperature T(stateSelect=StateSelect.prefer) = direct.thermal.T if
      n_spec > 0 and reduceThermal "Temperature";
    // These make the selected states more readable.

  protected
    Connectors.DirectNode direct(
      final n_trans=n_trans,
      final inclTrans=reduceTrans,
      final inclThermal=reduceThermal) if n_spec > 0 and (reduceTrans or
      reduceThermal)
      "Connector to directly couple velocities and temperatures within the phase"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={
              0,-60}), iconTransformation(extent={{-10,-10},{10,10}}, origin={
              26,-52})));

    Conditions.Adapters.AmagatDalton amagatDalton if n_spec > 0
      "Adapter between additivity of volume and additivity of pressure"
      annotation (Placement(transformation(extent={{-60,24},{-40,44}})));

    Connectors.InertNode electronsSolid
      "Connection node for drag betweene- and C+ (electrical resistance)"
      annotation (Placement(transformation(extent={{36,-58},{56,-38}}),
          iconTransformation(extent={{48,-76},{68,-56}})));
  equation
    // Chemical exchange
    connect('e-'.electrochemical, 'conne-') annotation (Line(
        points={{16,9},{16,50}},
        color={221,23,47},
        smooth=Smooth.None));

    // Electrical storage
    connect(electrostatic, 'e-'.electrostatic) annotation (Line(
        points={{70,34},{13,34},{13,7}},
        color={255,195,38},
        smooth=Smooth.None));

    // Inert exchange
    connect('C+'.inter, inter) annotation (Line(
        points={{-11,-3.8},{-9,-6},{-9,-36},{70,-36}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('C+'.direct.translational, direct.translational) annotation (Line(
        points={{-16.2,-9},{-16.2,-60},{0,-60}},
        color={47,107,251},
        smooth=Smooth.None));
    connect('C+'.direct.thermal, direct.thermal) annotation (Line(
        points={{-16.2,-9},{-16.2,-60},{0,-60}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('e-'.direct.translational, direct.translational) annotation (Line(
        points={{23.8,-9},{23.8,-60},{0,-60}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('e-'.direct.thermal, direct.thermal) annotation (Line(
        points={{23.8,-9},{23.8,-60},{0,-60}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('e-'.intra[1], electronsSolid.exchange) annotation (Line(
        points={{27,-7},{27,-48},{46,-48}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('C+'.intra[1], electronsSolid.exchange) annotation (Line(
        points={{-13,-7},{-13,-48},{46,-48}},
        color={38,196,52},
        smooth=Smooth.None));

    // Mixing
    connect(amagatDalton.amagat, amagat) annotation (Line(
        points={{-54,34},{-70,34}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(amagatDalton.dalton, 'C+'.dalton) annotation (Line(
        points={{-46,34},{-29,34},{-29,4}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(amagatDalton.dalton, 'e-'.dalton) annotation (Line(
        points={{-46,34},{11,34},{11,4}},
        color={47,107,251},
        smooth=Smooth.None));

    // Transport
    // ---------
    // C+
    connect('C+'.faces[transCart[Axis.x], Side.n], xNegative.'C+') annotation (
        Line(
        points={{-20,0},{-70,0}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('C+'.faces[transCart[Axis.x], Side.p], xPositive.'C+') annotation (
        Line(
        points={{-20,0},{70,0}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('C+'.faces[transCart[Axis.y], Side.n], yNegative.'C+') annotation (
        Line(
        points={{-20,0},{-20,-24},{-46,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('C+'.faces[transCart[Axis.y], Side.p], yPositive.'C+') annotation (
        Line(
        points={{-20,0},{-20,24},{46,24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('C+'.faces[transCart[Axis.z], Side.n], zNegative.'C+') annotation (
        Line(
        points={{-20,0},{-8,12},{58,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('C+'.faces[transCart[Axis.z], Side.p], zPositive.'C+') annotation (
        Line(
        points={{-20,0},{-32,-12},{-58,-12}},
        color={127,127,127},
        smooth=Smooth.None));

    // e-
    connect('e-'.faces[transCart[Axis.x], Side.n], xNegative.'e-') annotation (
        Line(
        points={{20,0},{-70,0}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.faces[transCart[Axis.x], Side.p], xPositive.'e-') annotation (
        Line(
        points={{20,0},{70,0}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.faces[transCart[Axis.y], Side.n], yNegative.'e-') annotation (
        Line(
        points={{20,0},{20,-24},{-46,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.faces[transCart[Axis.y], Side.p], yPositive.'e-') annotation (
        Line(
        points={{20,0},{20,24},{46,24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.faces[transCart[Axis.z], Side.n], zNegative.'e-') annotation (
        Line(
        points={{20,0},{32,12},{58,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('e-'.faces[transCart[Axis.z], Side.p], zPositive.'e-') annotation (
        Line(
        points={{20,0},{8,-12},{-58,-12}},
        color={127,127,127},
        smooth=Smooth.None));

    annotation (
      Documentation(info="<html>
    <p>See <a href=\"modelica://FCSys.Species.'e-'.Graphite.Fixed\">Species.'e-'.Graphite.Fixed</a> for assumptions.
    For more information, please see the
 <a href=\"modelica://FCSys.Phases.Partial\">Partial</a> model.</p></html>"),
      Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-80,-80},{80,
              60}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
              100}}), graphics));
  end Graphite;

  model Ionomer "Ionomer phase"
    import Modelica.Math.BooleanVectors.countTrue;

    extends Partial(final n_spec=countTrue({'inclSO3-','inclH+',inclH2O}));

    parameter Q.NumberAbsolute k_EOD=0.5
      "Coupling factor for electro-osmotic drag" annotation (Dialog(group=
            "Geometry", __Dymola_label="<html><i>k</i><sub>EOD</sub></html>"));
    // TODO:  Set an appropriate value (possibly use a new parameter for the current ratio).

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

    replaceable FCSys.Species.'SO3-'.Ionomer.Fixed 'SO3-'(final n_trans, final
        n_inter) if 'inclSO3-' constrainedby FCSys.Species.Solid(
      n_trans=n_trans,
      n_intra=2,
      n_inter=n_inter,
      k_intra={Modelica.Constants.inf,Modelica.Constants.inf},
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
      Placement(transformation(extent={{30,10},{50,30}})));

    parameter Boolean 'inclH+'=false "Include H+" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Protons (H<sup>+</sup>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.'H+'.Ionomer.Fixed 'H+'(final n_trans,final
        n_inter) if 'inclH+' constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      n_intra=2,
      n_inter=0,
      k_intra={1,k_EOD},
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
      Placement(transformation(extent={{-50,10},{-30,30}})));

    parameter Boolean inclH2O=false "Include H2O" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Water (H<sub>2</sub>O)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.H2O.Ionomer.Fixed H2O(final n_trans,final n_inter)
      if inclH2O constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      n_intra=2,
      n_inter=0,
      k_intra={1,k_EOD},
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
      Placement(transformation(extent={{-10,10},{10,30}})));

    Connectors.FaceBus xNegative if inclTrans[Axis.x]
      "Negative face along the x axis" annotation (Placement(transformation(
            extent={{-100,10},{-80,30}}), iconTransformation(extent={{-90,-10},
              {-70,10}})));
    Connectors.FaceBus yNegative if inclTrans[Axis.y]
      "Negative face along the y axis" annotation (Placement(transformation(
            extent={{-76,-14},{-56,6}}), iconTransformation(extent={{-10,-94},{
              10,-74}})));
    Connectors.FaceBus zNegative if inclTrans[Axis.z]
      "Negative face along the z axis" annotation (Placement(transformation(
            extent={{68,22},{88,42}}),iconTransformation(extent={{40,40},{60,60}})));
    Connectors.FaceBus xPositive if inclTrans[Axis.x]
      "Positive face along the x axis" annotation (Placement(transformation(
            extent={{80,10},{100,30}}), iconTransformation(extent={{70,-10},{90,
              10}})));
    Connectors.FaceBus yPositive if inclTrans[Axis.y]
      "Positive face along the y axis" annotation (Placement(transformation(
            extent={{56,34},{76,54}}),iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.FaceBus zPositive if inclTrans[Axis.z]
      "Positive face along the z axis" annotation (Placement(transformation(
            extent={{-88,-2},{-68,18}}), iconTransformation(extent={{-90,-90},{
              -70,-70}})));
    Connectors.Amagat amagat(final V=-V) if n_spec > 0
      "Connector for additivity of volume" annotation (Placement(transformation(
            extent={{-100,46},{-80,66}}), iconTransformation(extent={{70,-90},{
              90,-70}})));
    Connectors.Inter inter[n_inter](each final n_trans=n_trans) if n_spec > 0
      "Connector to exchange momentum and energy with other phases" annotation
      (Placement(transformation(extent={{80,-26},{100,-6}}), iconTransformation(
            extent={{-60,60},{-40,40}})));
    Connectors.Electrochemical connH2O(final n_trans=n_trans) if inclH2O
      "Chemical connector for H2O" annotation (Placement(transformation(extent=
              {{-14,60},{6,80}}), iconTransformation(extent={{10,-50},{30,-30}})));
    Connectors.Electrochemical 'connH+'(final n_trans=n_trans) if 'inclH+'
      "Electrochemical connector for H+" annotation (Placement(transformation(
            extent={{-54,60},{-34,80}}), iconTransformation(extent={{-30,-50},{
              -10,-30}})));
    Connectors.Electrostatic electrostatic if 'inclH+' and n_spec > 0
      "Interface with the dielectric" annotation (Placement(transformation(
            extent={{-100,34},{-80,54}}), iconTransformation(extent={{-110,-50},
              {-90,-30}})));
    // Note:  ('inclH+' and n_spec > 0) is logically equivalent to 'inclH+',
    // but n_spec is included so that the connector always appears in the icon
    // layer in Dymola 2014.

    // Aliases
    Q.Velocity phi[n_trans](each stateSelect=StateSelect.prefer) = direct.translational.phi
      if n_spec > 0 and reduceTrans "Velocity";
    Q.Temperature T(stateSelect=StateSelect.prefer) = direct.thermal.T if
      n_spec > 0 and reduceThermal "Temperature";
    // These make the selected states more readable.

  protected
    Conditions.Adapters.AmagatDalton amagatDalton if n_spec > 0
      "Adapter between additivity of volume and additivity of pressure"
      annotation (Placement(transformation(extent={{-80,46},{-60,66}})));
    Connectors.DirectNode direct(
      final n_trans=n_trans,
      final inclTrans=reduceTrans,
      final inclThermal=reduceThermal) if n_spec > 0 and (reduceTrans or
      reduceThermal)
      "Connector to directly couple velocities and temperatures within the phase"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={
              24,-64}), iconTransformation(extent={{-10,-10},{10,10}}, origin={
              26,-52})));

    Connectors.InertNode protonsSolid
      "Connection node for drag between H+ and SO3- (electrical resistance)"
      annotation (Placement(transformation(extent={{56,-38},{76,-18}}),
          iconTransformation(extent={{48,-76},{68,-56}})));
    Connectors.InertNode waterSolid
      "Connection node for drag between H2O and SO3-" annotation (Placement(
          transformation(extent={{56,-50},{76,-30}}), iconTransformation(extent
            ={{48,-76},{68,-56}})));
    Connectors.InertNode EOD
      "Connection node for electro-osmotic drag (between H+ and H2O)"
      annotation (Placement(transformation(extent={{56,-62},{76,-42}}),
          iconTransformation(extent={{48,-76},{68,-56}})));
  equation
    // Chemical exchange
    connect(H2O.electrochemical, connH2O) annotation (Line(
        points={{-4,29},{-4,70}},
        color={221,23,47},
        smooth=Smooth.None));
    connect('H+'.electrochemical, 'connH+') annotation (Line(
        points={{-44,29},{-44,70}},
        color={221,23,47},
        smooth=Smooth.None));

    // Electrical storage
    connect(electrostatic, 'H+'.electrostatic) annotation (Line(
        points={{-90,44},{-47,44},{-47,27}},
        color={255,195,38},
        smooth=Smooth.None));

    // Inert exchange
    connect('SO3-'.inter, inter) annotation (Line(
        points={{49,16.2},{52,16.2},{52,-16},{90,-16}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('SO3-'.direct.translational, direct.translational) annotation (Line(
        points={{43.8,11},{43.8,-64},{24,-64}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('SO3-'.direct.thermal, direct.thermal) annotation (Line(
        points={{43.8,11},{43.8,-64},{24,-64}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('H+'.direct.translational, direct.translational) annotation (Line(
        points={{-36.2,11},{-36.2,-64},{24,-64}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('H+'.direct.thermal, direct.thermal) annotation (Line(
        points={{-36.2,11},{-36.2,-64},{24,-64}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2O.direct.translational, direct.translational) annotation (Line(
        points={{3.8,11},{3.8,-64},{24,-64}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2O.direct.thermal, direct.thermal) annotation (Line(
        points={{3.8,11},{3.8,-64},{24,-64}},
        color={38,196,52},
        smooth=Smooth.None));

    connect('SO3-'.intra[1], protonsSolid.exchange) annotation (Line(
        points={{47,13},{47,-28},{66,-28}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('SO3-'.intra[2], waterSolid.exchange) annotation (Line(
        points={{47,13},{47,-40},{66,-40}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('H+'.intra[1], protonsSolid.exchange) annotation (Line(
        points={{-33,13},{-33,-28},{66,-28}},
        color={38,196,52},
        smooth=Smooth.None));
    connect('H+'.intra[2], EOD.exchange) annotation (Line(
        points={{-33,13},{-33,-52},{66,-52}},
        color={38,196,52},
        smooth=Smooth.None));

    connect(H2O.intra[1], waterSolid.exchange) annotation (Line(
        points={{7,13},{7,-40},{66,-40}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(H2O.intra[2], EOD.exchange) annotation (Line(
        points={{7,13},{7,-52},{66,-52}},
        color={38,196,52},
        smooth=Smooth.None));

    // Mixing
    connect('SO3-'.dalton, amagatDalton.dalton) annotation (Line(
        points={{31,24},{31,56},{-66,56}},
        color={47,107,251},
        smooth=Smooth.None));
    connect('H+'.dalton, amagatDalton.dalton) annotation (Line(
        points={{-49,24},{-49,56},{-66,56}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(H2O.dalton, amagatDalton.dalton) annotation (Line(
        points={{-9,24},{-9,56},{-66,56}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(amagatDalton.amagat, amagat) annotation (Line(
        points={{-74,56},{-90,56}},
        color={47,107,251},
        smooth=Smooth.None));

    // Transport
    // ---------
    // SO3-
    connect('SO3-'.faces[transCart[Axis.x], Side.n], xNegative.'SO3-')
      annotation (Line(
        points={{40,20},{-90,20}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('SO3-'.faces[transCart[Axis.x], Side.p], xPositive.'SO3-')
      annotation (Line(
        points={{40,20},{90,20}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('SO3-'.faces[transCart[Axis.y], Side.n], yNegative.'SO3-')
      annotation (Line(
        points={{40,20},{40,-4},{-66,-4}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('SO3-'.faces[transCart[Axis.y], Side.p], yPositive.'SO3-')
      annotation (Line(
        points={{40,20},{40,44},{66,44}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('SO3-'.faces[transCart[Axis.z], Side.n], zNegative.'SO3-')
      annotation (Line(
        points={{40,20},{52,32},{78,32}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('SO3-'.faces[transCart[Axis.z], Side.p], zPositive.'SO3-')
      annotation (Line(
        points={{40,20},{28,8},{-78,8}},
        color={127,127,127},
        smooth=Smooth.None));
    // 'H+'
    connect('H+'.faces[transCart[Axis.x], Side.n], xNegative.'H+') annotation (
        Line(
        points={{-40,20},{-90,20}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.faces[transCart[Axis.x], Side.p], xPositive.'H+') annotation (
        Line(
        points={{-40,20},{90,20}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.faces[transCart[Axis.y], Side.n], yNegative.'H+') annotation (
        Line(
        points={{-40,20},{-40,-4},{-66,-4}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.faces[transCart[Axis.y], Side.p], yPositive.'H+') annotation (
        Line(
        points={{-40,20},{-40,44},{66,44}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.faces[transCart[Axis.z], Side.n], zNegative.'H+') annotation (
        Line(
        points={{-40,20},{-28,32},{78,32}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('H+'.faces[transCart[Axis.z], Side.p], zPositive.'H+') annotation (
        Line(
        points={{-40,20},{-52,8},{-78,8}},
        color={127,127,127},
        smooth=Smooth.None));
    // H2O
    connect(H2O.faces[transCart[Axis.x], Side.n], xNegative.H2O) annotation (
        Line(
        points={{0,20},{-90,20}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[transCart[Axis.x], Side.p], xPositive.H2O) annotation (
        Line(
        points={{0,20},{90,20}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[transCart[Axis.y], Side.n], yNegative.H2O) annotation (
        Line(
        points={{0,20},{0,-4},{-66,-4}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[transCart[Axis.y], Side.p], yPositive.H2O) annotation (
        Line(
        points={{0,20},{0,44},{66,44}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[transCart[Axis.z], Side.n], zNegative.H2O) annotation (
        Line(
        points={{0,20},{12,32},{78,32}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(H2O.faces[transCart[Axis.z], Side.p], zPositive.H2O) annotation (
        Line(
        points={{0,20},{0,8},{-78,8}},
        color={127,127,127},
        smooth=Smooth.None));

    annotation (
      Documentation(info="<html><p>Assumptions:<ol>
    <li>The water in the ionomer does not participate in the reaction (only the water vapor does).</li>
    </ol</p>
    
    <p>See <a href=\"modelica://FCSys.Species.'H+'.Ionomer.Fixed\">Species.'H+'.Ionomer.Fixed</a> 
    for additional assumptions.
    For more information, please see the
 <a href=\"modelica://FCSys.Phases.Partial\">Partial</a> model.</p></html>"),
      Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-80},{
              100,80}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
              100}}), graphics));
  end Ionomer;

  model Liquid "Liquid phase"
    extends Partial(final n_spec=if inclH2O then 1 else 0,reduceThermal=true);

    // Conditionally include species.
    parameter Boolean inclH2O=false "Include H2O" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Water (H<sub>2</sub>O)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.H2O.Liquid.Fixed H2O(final n_trans,final n_inter)
      if inclH2O constrainedby FCSys.Species.Fluid(
      n_trans=n_trans,
      n_inter=n_inter,
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
    Connectors.Electrochemical connH2O(final n_trans=n_trans) if inclH2O
      "Chemical connector for H2O" annotation (Placement(transformation(extent=
              {{-24,40},{-4,60}}), iconTransformation(extent={{-10,-50},{10,-30}})));

  protected
    Conditions.Adapters.AmagatDalton amagatDalton if n_spec > 0
      "Adapter between additivity of volume and additivity of pressure"
      annotation (Placement(transformation(extent={{-40,0},{-20,20}})));

  equation
    // Chemical exchange
    connect(H2O.electrochemical, connH2O) annotation (Line(
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
  partial model Partial "Base model for a phase"
    import Modelica.Math.BooleanVectors.index;
    // extends FCSys.Icons.Names.Middle;

    parameter Integer n_spec(start=0) "Number of species"
      annotation (HideResult=true);
    parameter Integer n_inter=0
      "Number of exchange connections with other phases"
      annotation (Dialog(connectorSizing=true),HideResult=true);

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
      annotation (Dialog(tab="Assumptions", enable=n_spec > 1), choices(
          __Dymola_checkBox=true));

    inner parameter Boolean reduceThermal=false
      "Same temperature for all species" annotation (Dialog(tab="Assumptions",
          enable=n_spec > 1), choices(__Dymola_checkBox=true));

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
    the conductance of material, transverse translational momentum, and energy of all of the species
    within the phase.  It can be used to introduce minor head loss or the effects of
    porosity or tortousity.  These effects may be anisotropic.  Using
    Bruggeman correction [<a href=\"modelica://FCSys.UsersGuide.References\">Weber2004</a>, p. 4696],
    the factor (<b><i>k</i></b>) within a phase should be set to &epsilon;<sup>1/2</sup>
    along each axis, where &epsilon; is the volumetric filling ratio, or the ratio of the volume of the phase to the total volume of the subregion.
    The Bruggeman factor itself is &epsilon;<sup>3/2</sup>, but a factor of &epsilon; is included inherently.</p>
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
  end Partial;

public
  model Dielectric "Dielectric gap"
    extends FCSys.Icons.Names.Top2;

    parameter Q.Length L=1e-9*U.m
      "Length of the dielectric (not of the region)" annotation (Dialog(group=
            "Geometry",__Dymola_label="<html><i>L</i></html>"));
    parameter Q.Area A=100*U.m^2
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
