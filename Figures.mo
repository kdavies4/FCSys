within FCSys;
package Figures "Graphical layouts for documentation"
  extends Modelica.Icons.Package;
  package DeclarativeVsImperative
    import SI = Modelica.SIunits;
    extends Modelica.Icons.Package;
    // See notes from 12/7/10 for the derivation of the equations in this
    // package.

    model ResistorDeclarative "Electrical resistor in declarative formalism"

      parameter SI.Resistance r=1 "Resistance";

      SI.Voltage Deltav "Voltage difference";
      SI.Current i "Current";

      Pin n "Electrical pin on the negative side";
      Pin p "Electrical pin on the positive side";

    protected
      connector Pin
        SI.Voltage v;
        flow SI.Current i;

      end Pin;

    equation
      Deltav = p.v - n.v "Voltage across the resistor";
      i = (n.i - p.i)/2 "Current through the resistor";
      Deltav = i*r "Transport (Ohm's law)";
      0 = n.i + p.i "Conservation of charge (no storage)";

    end ResistorDeclarative;

    function ResistorImperative "Electrical resistor in imperative formalism"
      extends Modelica.Icons.Function;

      parameter SI.Resistance r=1 "Resistance";

      input SI.Current i "Current";
      output SI.Voltage Deltav "Voltage difference";

    algorithm
      Deltav := i*r "Transport (Ohm's law)";

    end ResistorImperative;

    package SwitchCausality
      "Examples showing opposite causalities of electrical circuit"
      package Examples "Examples"
        extends Modelica.Icons.ExamplesPackage;

        model Test_vi
          extends Modelica.Icons.Example;
          Modelica.Blocks.Sources.Pulse pulse(
            period=1,
            offset=-0.5,
            startTime=0)
            annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
          Declarative_vi declarative_vi
            annotation (Placement(transformation(extent={{0,30},{20,50}})));
          Imperative_vi imperative_vi
            annotation (Placement(transformation(extent={{0,-10},{20,10}})));
          ImperativeTF_vi imperativeTF_vi
            annotation (Placement(transformation(extent={{0,-50},{20,-30}})));

        equation
          connect(pulse.y, declarative_vi.v) annotation (Line(
              points={{-19,6.10623e-16},{-10,6.10623e-16},{-10,40},{-1,40}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(imperativeTF_vi.v, pulse.y) annotation (Line(
              points={{-1,-40},{-10,-40},{-10,6.10623e-16},{-19,6.10623e-16}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(imperative_vi.v, pulse.y) annotation (Line(
              points={{-1,6.10623e-16},{-5.5,6.10623e-16},{-5.5,6.10623e-16},{-10,
                  6.10623e-16},{-10,6.10623e-16},{-19,6.10623e-16}},
              color={0,0,127},
              smooth=Smooth.None));
          annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent
                  ={{-100,-100},{100,100}}), graphics), experiment(StopTime=2));
        end Test_vi;

        model Test_iv
          extends Modelica.Icons.Example;
          Modelica.Blocks.Sources.Pulse pulse(
            period=1,
            offset=-0.5,
            startTime=0)
            annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
          Declarative_iv declarative_iv
            annotation (Placement(transformation(extent={{0,30},{20,50}})));
          FCSys.Figures.DeclarativeVsImperative.SwitchCausality.Imperative_iv
            imperative_iv
            annotation (Placement(transformation(extent={{0,-10},{20,10}})));
          ImperativeTF_iv imperativeTF_iv
            annotation (Placement(transformation(extent={{0,-50},{20,-30}})));

        equation
          connect(pulse.y, declarative_iv.i) annotation (Line(
              points={{-19,6.10623e-16},{-10,6.10623e-16},{-10,40},{-1,40}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(imperative_iv.i, pulse.y) annotation (Line(
              points={{-1,6.10623e-16},{-5.5,6.10623e-16},{-5.5,6.10623e-16},{-10,
                  6.10623e-16},{-10,6.10623e-16},{-19,6.10623e-16}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(imperativeTF_iv.i, pulse.y) annotation (Line(
              points={{-1,-40},{-10,-40},{-10,6.10623e-16},{-19,6.10623e-16}},
              color={0,0,127},
              smooth=Smooth.None));
          annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent
                  ={{-100,-100},{100,100}}), graphics), experiment(StopTime=2));
        end Test_iv;

        model Declarative
          extends
            FCSys.Figures.DeclarativeVsImperative.SwitchCausality.BaseClasses.Parameters;
          extends Modelica.Icons.Example;
          Modelica.Electrical.Analog.Basic.Resistor res1(R=R1) annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={30,10})));
          Modelica.Electrical.Analog.Basic.Resistor res2(R=R2) annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={60,10})));
          Modelica.Electrical.Analog.Basic.Inductor ind(L=L) annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={60,-20})));
          Modelica.Electrical.Analog.Basic.Capacitor cap(C=C) annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={30,-20})));
          Modelica.Electrical.Analog.Basic.Ground ground
            annotation (Placement(transformation(extent={{20,-60},{40,-40}})));

        equation
          connect(res2.p, res1.p) annotation (Line(
              points={{60,20},{60,30},{30,30},{30,20}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cap.p, res1.n) annotation (Line(
              points={{30,-10},{30,-7.5},{30,-7.5},{30,-5},{30,1.22125e-15},{30,
                  1.22125e-15}},
              color={0,0,255},
              smooth=Smooth.None));

          connect(ind.p, res2.n) annotation (Line(
              points={{60,-10},{60,1.22125e-15}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(ind.n, ground.p) annotation (Line(
              points={{60,-30},{60,-40},{30,-40}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cap.n, ground.p) annotation (Line(
              points={{30,-30},{30,-35},{30,-35},{30,-40}},
              color={0,0,255},
              smooth=Smooth.None));
          annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent
                  ={{-80,-60},{80,60}}), graphics));
        end Declarative;

      end Examples;
      extends Modelica.Icons.Package;

      model Declarative_vi
        "Declarative-based circuit with voltage in, current out"
        extends BaseClasses.Parameters;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
        Connectors.RealInput v annotation (Placement(transformation(extent={{-60,
                  -30},{-40,-10}}), iconTransformation(extent={{-120,-10},{-100,
                  10}})));
        Connectors.RealOutput i annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={10,10}), iconTransformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={110,0})));
        Modelica.Electrical.Analog.Basic.Resistor res1(R=R1) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={30,10})));
        Modelica.Electrical.Analog.Basic.Capacitor cap(C=C) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={30,-20})));
        Modelica.Electrical.Analog.Basic.Inductor ind(L=L) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={60,-20})));
        Modelica.Electrical.Analog.Basic.Ground ground
          annotation (Placement(transformation(extent={{20,-60},{40,-40}})));
        Modelica.Electrical.Analog.Basic.Resistor res2(R=R2) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={60,10})));
        Modelica.Electrical.Analog.Sources.SignalVoltage voltageSource
          annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=270,
              origin={-20,-20})));
        Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor
          annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=270,
              origin={-20,10})));

      equation
        connect(res2.p, res1.p) annotation (Line(
            points={{60,20},{60,30},{30,30},{30,20}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(cap.p, res1.n) annotation (Line(
            points={{30,-10},{30,-7.5},{30,-7.5},{30,-5},{30,1.22125e-15},{30,
                1.22125e-15}},
            color={0,0,255},
            smooth=Smooth.None));

        connect(ind.p, res2.n) annotation (Line(
            points={{60,-10},{60,1.22125e-15}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(ind.n, ground.p) annotation (Line(
            points={{60,-30},{60,-40},{30,-40}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(voltageSource.n, ground.p) annotation (Line(
            points={{-20,-30},{-20,-40},{30,-40}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(v, voltageSource.v) annotation (Line(
            points={{-50,-20},{-27,-20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(i, currentSensor.i) annotation (Line(
            points={{10,10},{0,10},{0,10},{-10,10}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(cap.n, ground.p) annotation (Line(
            points={{30,-30},{30,-35},{30,-35},{30,-40}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(i, i) annotation (Line(
            points={{10,10},{10,10}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(currentSensor.p, res1.p) annotation (Line(
            points={{-20,20},{-20,30},{30,30},{30,20}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(currentSensor.n, voltageSource.p) annotation (Line(
            points={{-20,-5.55112e-16},{-20,-2.5},{-20,-2.5},{-20,-5},{-20,-10},
                {-20,-10}},
            color={0,0,255},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-80,-60},{80,60}}), graphics));
      end Declarative_vi;

      model Declarative_iv
        "Declarative-based circuit with current in, voltage out"
        extends BaseClasses.Parameters;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
        Connectors.RealInput i annotation (Placement(transformation(extent={{-80,
                  -30},{-60,-10}}), iconTransformation(extent={{-120,-10},{-100,
                  10}})));
        Connectors.RealOutput v annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={10,10}), iconTransformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={110,0})));
        Modelica.Electrical.Analog.Basic.Resistor res1(R=R1) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={30,10})));
        Modelica.Electrical.Analog.Basic.Capacitor cap(C=C) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={30,-20})));
        Modelica.Electrical.Analog.Basic.Inductor ind(L=L) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={60,-20})));
        Modelica.Electrical.Analog.Basic.Ground ground
          annotation (Placement(transformation(extent={{20,-60},{40,-40}})));
        Modelica.Electrical.Analog.Basic.Resistor res2(R=R2) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={60,10})));
        Modelica.Electrical.Analog.Sources.SignalCurrent currentSource
          annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=270,
              origin={-40,-20})));
        Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
          annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=270,
              origin={-20,10})));

      equation
        connect(currentSource.i, i) annotation (Line(
            points={{-47,-20},{-70,-20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(voltageSensor.v, v) annotation (Line(
            points={{-10,10},{-4,10},{-4,10},{10,10}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(res1.p, res2.p) annotation (Line(
            points={{30,20},{30,30},{60,30},{60,20}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(res1.n, cap.p) annotation (Line(
            points={{30,1.22125e-15},{30,-2.5},{30,-2.5},{30,-5},{30,-10},{30,-10}},

            color={0,0,255},
            smooth=Smooth.None));

        connect(res2.n, ind.p) annotation (Line(
            points={{60,1.22125e-15},{60,-10}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(currentSource.n, ground.p) annotation (Line(
            points={{-40,-30},{-40,-40},{30,-40}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(cap.n, ground.p) annotation (Line(
            points={{30,-30},{30,-35},{30,-35},{30,-40}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(ind.n, ground.p) annotation (Line(
            points={{60,-30},{60,-40},{30,-40}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(currentSource.p, res1.p) annotation (Line(
            points={{-40,-10},{-40,30},{30,30},{30,20}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(voltageSensor.p, currentSource.p) annotation (Line(
            points={{-20,20},{-20,30},{-40,30},{-40,-10}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(voltageSensor.n, ground.p) annotation (Line(
            points={{-20,-5.55112e-16},{-20,-40},{30,-40}},
            color={0,0,255},
            smooth=Smooth.None));
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics), Diagram(coordinateSystem(
                preserveAspectRatio=true, extent={{-80,-60},{80,60}}), graphics));
      end Declarative_iv;

      model Imperative_vi "Imperative circuit with voltage in, current out"
        extends BaseClasses.Parameters;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
        Connectors.RealInput v annotation (Placement(transformation(extent={{-100,
                  20},{-80,40}}), iconTransformation(extent={{-120,-10},{-100,
                  10}})));
        Connectors.RealOutput i annotation (Placement(transformation(extent={{
                  80,20},{100,40}}), iconTransformation(extent={{100,-10},{120,
                  10}})));
        Modelica.Blocks.Continuous.Integrator ind(k=1/L)
          annotation (Placement(transformation(extent={{-20,-10},{0,10}})));
        Modelica.Blocks.Math.Gain res1(k=1/R1)
          annotation (Placement(transformation(extent={{-20,50},{0,70}})));
        Modelica.Blocks.Math.Add diff1(k1=-1, k2=1)
          annotation (Placement(transformation(extent={{-60,50},{-40,70}})));
        Modelica.Blocks.Math.Add sum(k1=-1, k2=-1)
          annotation (Placement(transformation(extent={{50,20},{70,40}})));
        Modelica.Blocks.Math.Add diff2(k2=-1)
          annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
        Modelica.Blocks.Continuous.Integrator cap(k=1/C)
          annotation (Placement(transformation(extent={{20,50},{40,70}})));
        Modelica.Blocks.Math.Gain res2(k=R2)
          annotation (Placement(transformation(extent={{20,-10},{40,10}})));

      equation
        connect(diff1.y, res1.u) annotation (Line(
            points={{-39,60},{-22,60}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(sum.y, i) annotation (Line(
            points={{71,30},{90,30}},
            color={0,200,0},
            smooth=Smooth.None));
        connect(diff2.y, ind.u) annotation (Line(
            points={{-39,6.10623e-16},{-34,0},{-30,1.27676e-15},{-30,
                6.66134e-16},{-22,6.66134e-16}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(cap.u, res1.y) annotation (Line(
            points={{18,60},{1,60}},
            color={0,200,0},
            smooth=Smooth.None));
        connect(res2.u, ind.y) annotation (Line(
            points={{18,6.66134e-16},{14,0},{10,1.27676e-15},{10,6.10623e-16},{
                1,6.10623e-16}},
            color={0,200,0},
            smooth=Smooth.None));

        connect(diff1.u2, v) annotation (Line(
            points={{-62,54},{-70,54},{-70,30},{-90,30}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(diff2.u1, v) annotation (Line(
            points={{-62,6},{-70,6},{-70,30},{-90,30}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(ind.y, sum.u2) annotation (Line(
            points={{1,6.10623e-16},{10,6.10623e-16},{10,24},{48,24}},
            color={0,200,0},
            smooth=Smooth.None));
        connect(sum.u1, res1.y) annotation (Line(
            points={{48,36},{10,36},{10,60},{1,60}},
            color={0,200,0},
            smooth=Smooth.None));
        connect(cap.y, diff1.u1) annotation (Line(
            points={{41,60},{50,60},{50,80},{-70,80},{-70,66},{-62,66}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(res2.y, diff2.u2) annotation (Line(
            points={{41,6.10623e-16},{50,6.10623e-16},{50,-20},{-70,-20},{-70,-6},
                {-62,-6}},
            color={255,195,38},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-100},{100,100}}), graphics={Line(
                      points={{62,76},{70,76}},
                      color={255,195,38},
                      smooth=Smooth.None),Rectangle(
                      extent={{-16,8},{16,-8}},
                      lineColor={0,0,0},
                      origin={74,72},
                      rotation=180),Text(
                      extent={{70,78},{88,74}},
                      lineColor={0,0,0},
                      textString="Voltage"),Text(
                      extent={{70,70},{88,66}},
                      lineColor={0,0,0},
                      textString="Current"),Line(
                      points={{62,68},{70,68}},
                      color={0,200,0},
                      smooth=Smooth.None)}));
      end Imperative_vi;

      model Imperative_iv "Imperative circuit with current in, voltage out"
        extends BaseClasses.Parameters;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;

      public
        Connectors.RealInput i annotation (Placement(transformation(extent={{-100,
                  20},{-80,40}}), iconTransformation(extent={{-120,-10},{-100,
                  10}})));
        Connectors.RealOutput v annotation (Placement(transformation(extent={{
                  80,20},{100,40}}), iconTransformation(extent={{100,-10},{120,
                  10}})));
        Modelica.Blocks.Math.Gain res1(k=R1)
          annotation (Placement(transformation(extent={{-10,0},{10,20}})));
        Modelica.Blocks.Math.Add sum
          annotation (Placement(transformation(extent={{30,20},{50,40}})));
        Modelica.Blocks.Math.Add diff1(k2=-1, k1=-1)
          annotation (Placement(transformation(extent={{-50,20},{-30,40}})));
        Modelica.Blocks.Continuous.Integrator cap(k=1/C)
          annotation (Placement(transformation(extent={{-10,40},{10,60}})));
        Modelica.Blocks.Continuous.Integrator ind(k=1/L)
          annotation (Placement(transformation(extent={{10,-40},{-10,-20}})));
        Modelica.Blocks.Math.Add diff2(k2=-1)
          annotation (Placement(transformation(extent={{50,-40},{30,-20}})));
        Modelica.Blocks.Math.Gain res2(k=R2)
          annotation (Placement(transformation(extent={{-10,-80},{10,-60}})));

      equation
        connect(res2.u, ind.y) annotation (Line(
            points={{-12,-70},{-60,-70},{-60,-30},{-11,-30}},
            color={0,200,0},
            smooth=Smooth.None));
        connect(diff2.y, ind.u) annotation (Line(
            points={{29,-30},{12,-30}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(i, diff1.u1) annotation (Line(
            points={{-90,30},{-70,30},{-70,36},{-52,36}},
            color={0,200,0},
            smooth=Smooth.None));
        connect(diff1.u2, ind.y) annotation (Line(
            points={{-52,24},{-60,24},{-60,-30},{-11,-30}},
            color={0,200,0},
            smooth=Smooth.None));
        connect(cap.u, diff1.y) annotation (Line(
            points={{-12,50},{-20,50},{-20,30},{-29,30}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(diff1.y, res1.u) annotation (Line(
            points={{-29,30},{-20,30},{-20,10},{-12,10}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(cap.y, sum.u1) annotation (Line(
            points={{11,50},{20,50},{20,36},{28,36}},
            color={0,200,0},
            smooth=Smooth.None));
        connect(sum.u2, res1.y) annotation (Line(
            points={{28,24},{20,24},{20,10},{11,10}},
            color={0,200,0},
            smooth=Smooth.None));
        connect(diff2.u2, res2.y) annotation (Line(
            points={{52,-36},{60,-36},{60,-70},{11,-70}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(sum.y, diff2.u1) annotation (Line(
            points={{51,30},{60,30},{60,-24},{52,-24}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(sum.y, v) annotation (Line(
            points={{51,30},{90,30}},
            color={255,195,38},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-100},{100,100}}), graphics={Line(
                      points={{60,58},{68,58}},
                      color={255,195,38},
                      smooth=Smooth.None),Rectangle(
                      extent={{-16,8},{16,-8}},
                      lineColor={0,0,0},
                      origin={72,54},
                      rotation=180),Text(
                      extent={{68,60},{86,56}},
                      lineColor={0,0,0},
                      textString="Voltage"),Text(
                      extent={{68,52},{86,48}},
                      lineColor={0,0,0},
                      textString="Current"),Line(
                      points={{60,50},{68,50}},
                      color={0,200,0},
                      smooth=Smooth.None)}));
      end Imperative_iv;

      model ImperativeTF_vi
        "Equivalent transfer function for voltage in, current out"
        extends BaseClasses.Parameters;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
        Connectors.RealInput v annotation (Placement(transformation(extent={{-100,
                  20},{-80,40}}), iconTransformation(extent={{-120,-10},{-100,
                  10}})));
        Connectors.RealOutput i annotation (Placement(transformation(extent={{
                  80,20},{100,40}}), iconTransformation(extent={{100,-10},{120,
                  10}})));
        Modelica.Blocks.Continuous.TransferFunction transferFunction(a={R1*L*C,
              R1*R2*C + L,R2}, b={-L*C,-(R1 + R2)*C,-1})
          annotation (Placement(transformation(extent={{-10,20},{10,40}})));

      equation
        connect(v, transferFunction.u) annotation (Line(
            points={{-90,30},{-12,30}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(transferFunction.y, i) annotation (Line(
            points={{11,30},{90,30}},
            color={0,200,0},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-100},{100,100}}), graphics));
      end ImperativeTF_vi;

      model ImperativeTF_iv
        "Equivalent transfer function for voltage in, current out"
        extends BaseClasses.Parameters;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;

      public
        Connectors.RealInput i annotation (Placement(transformation(extent={{-100,
                  20},{-80,40}}), iconTransformation(extent={{-120,-10},{-100,
                  10}})));
        Connectors.RealOutput v annotation (Placement(transformation(extent={{
                  80,20},{100,40}}), iconTransformation(extent={{100,-10},{120,
                  10}})));
        Modelica.Blocks.Continuous.TransferFunction transferFunction(b={R1*L*C,
              R1*R2*C + L,R2}, a={-L*C,-(R1 + R2)*C,-1})
          annotation (Placement(transformation(extent={{-10,20},{10,40}})));

      equation
        connect(transferFunction.y, v) annotation (Line(
            points={{11,30},{90,30}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(transferFunction.u, i) annotation (Line(
            points={{-12,30},{-90,30}},
            color={0,200,0},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-100},{100,100}}), graphics));
      end ImperativeTF_iv;

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;

        partial model Parameters
          "Base model with parameters for electrical circuit"
          parameter SI.Resistance R1=10 "Resistance at temperature T_ref";
          parameter SI.Resistance R2=100 "Resistance at temperature T_ref";
          parameter SI.Inductance L=0.1 "Inductance";
          parameter SI.Capacitance C=0.01 "Capacitance";

        end Parameters;

      end BaseClasses;

    end SwitchCausality;

    package Cascade "Examples showing cascaded electrical circuits"
      package Examples "Examples"
        extends Modelica.Icons.ExamplesPackage;

        model Test
          extends Modelica.Icons.Example;
          Modelica.Blocks.Sources.Pulse pulse(
            period=1,
            offset=-0.5,
            startTime=0)
            annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
          DeclarativeAB declarativeAB
            annotation (Placement(transformation(extent={{0,50},{20,70}})));
          ImperativeAB imperativeAB
            annotation (Placement(transformation(extent={{0,10},{20,30}})));
          ImperativeABTF imperativeABTF
            annotation (Placement(transformation(extent={{0,-30},{20,-10}})));
          ImperativeABIncorrect imperativeABIncorrect
            annotation (Placement(transformation(extent={{0,-70},{20,-50}})));

        equation
          connect(declarativeAB.vIn, pulse.y) annotation (Line(
              points={{-1,60},{-10,60},{-10,6.10623e-16},{-19,6.10623e-16}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(pulse.y, imperativeAB.vIn) annotation (Line(
              points={{-19,6.10623e-16},{-10,6.10623e-16},{-10,20},{-1,20}},
              color={0,0,127},
              smooth=Smooth.None));
          connect(pulse.y, imperativeABIncorrect.vIn) annotation (Line(
              points={{-19,6.10623e-16},{-10,6.10623e-16},{-10,-60},{-1,-60}},
              color={0,0,127},
              smooth=Smooth.None));

          connect(imperativeABTF.vIn, pulse.y) annotation (Line(
              points={{-1,-20},{-10,-20},{-10,6.10623e-16},{-19,6.10623e-16}},
              color={0,0,127},
              smooth=Smooth.None));
          annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent
                  ={{-100,-100},{100,100}}), graphics), experiment(StopTime=2));
        end Test;

      end Examples;
      extends Modelica.Icons.Package;

      model DeclarativeA "First circuit in declarative formalism"
        extends BaseClasses.Parameters;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
        Connectors.RealInput vIn annotation (Placement(transformation(extent={{
                  -100,-10},{-80,10}}), iconTransformation(extent={{-120,-10},{
                  -100,10}})));
        Connectors.RealOutput vOut annotation (Placement(transformation(extent=
                  {{40,-10},{60,10}}), iconTransformation(extent={{100,-10},{
                  120,10}})));
        Modelica.Electrical.Analog.Sources.SignalVoltage voltageSource
          annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=270,
              origin={-50,-20})));
        Modelica.Electrical.Analog.Basic.Resistor res1(R=R1) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-20,20})));
        Modelica.Electrical.Analog.Basic.Capacitor cap(C=C) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-20,-20})));
        Modelica.Electrical.Analog.Basic.Ground ground
          annotation (Placement(transformation(extent={{-30,-60},{-10,-40}})));
        Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=90,
              origin={10,-20})));

      equation
        connect(cap.n, ground.p) annotation (Line(
            points={{-20,-30},{-20,-35},{-20,-35},{-20,-40}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(cap.p, res1.n) annotation (Line(
            points={{-20,-10},{-20,-5},{-20,-5},{-20,0},{-20,10},{-20,10}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(voltageSource.p, res1.p) annotation (Line(
            points={{-50,-10},{-50,40},{-20,40},{-20,30}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(vIn, voltageSource.v) annotation (Line(
            points={{-90,5.55112e-16},{-82,5.55112e-16},{-82,0},{-70,0},{-70,-20},
                {-57,-20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(voltageSource.n, ground.p) annotation (Line(
            points={{-50,-30},{-50,-40},{-20,-40}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(vOut, voltageSensor.v) annotation (Line(
            points={{50,5.55112e-16},{30,5.55112e-16},{30,-20},{20,-20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(voltageSensor.n, ground.p) annotation (Line(
            points={{10,-30},{10,-40},{-20,-40}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(voltageSensor.p, res1.n) annotation (Line(
            points={{10,-10},{10,0},{-20,0},{-20,10}},
            color={0,0,255},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-60},{100,60}}), graphics));
      end DeclarativeA;

      model DeclarativeB "First circuit in declarative formalism"
        extends BaseClasses.Parameters;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
        Connectors.RealInput vIn annotation (Placement(transformation(extent={{
                  -60,-10},{-40,10}}), iconTransformation(extent={{-120,-10},{-100,
                  10}})));
        Connectors.RealOutput vOut annotation (Placement(transformation(extent=
                  {{80,-10},{100,10}}), iconTransformation(extent={{100,-10},{
                  120,10}})));
        Modelica.Electrical.Analog.Sources.SignalVoltage voltageSource
          annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=270,
              origin={-10,-20})));
        Modelica.Electrical.Analog.Basic.Resistor res2(R=R2) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={20,20})));
        Modelica.Electrical.Analog.Basic.Resistor res3(R=R3) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={20,-20})));
        Modelica.Electrical.Analog.Basic.Ground ground
          annotation (Placement(transformation(extent={{10,-60},{30,-40}})));
        Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=90,
              origin={50,-20})));

      equation
        connect(voltageSource.p, res2.p) annotation (Line(
            points={{-10,-10},{-10,40},{20,40},{20,30}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(voltageSource.n, ground.p) annotation (Line(
            points={{-10,-30},{-10,-40},{20,-40}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(res3.p, res2.n) annotation (Line(
            points={{20,-10},{20,-5},{20,-5},{20,0},{20,10},{20,10}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(res3.n, ground.p) annotation (Line(
            points={{20,-30},{20,-35},{20,-35},{20,-40}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(voltageSensor.v, vOut) annotation (Line(
            points={{60,-20},{70,-20},{70,5.55112e-16},{90,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(vIn, voltageSource.v) annotation (Line(
            points={{-50,5.55112e-16},{-30,5.55112e-16},{-30,-20},{-17,-20}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(voltageSensor.p, res2.n) annotation (Line(
            points={{50,-10},{50,0},{20,0},{20,10}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(voltageSensor.n, ground.p) annotation (Line(
            points={{50,-30},{50,-40},{20,-40}},
            color={0,0,255},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-60},{100,60}}), graphics));
      end DeclarativeB;

      model DeclarativeAB
        "Cascaded first and second circuits in declarative formalism"
        extends BaseClasses.Parameters;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
        Connectors.RealInput vIn annotation (Placement(transformation(extent={{
                  -100,-10},{-80,10}}), iconTransformation(extent={{-120,-10},{
                  -100,10}})));
        Connectors.RealOutput vOut annotation (Placement(transformation(extent=
                  {{80,-10},{100,10}}), iconTransformation(extent={{100,-10},{
                  120,10}})));
        Modelica.Electrical.Analog.Sources.SignalVoltage voltageSource
          annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=270,
              origin={-50,-20})));
        Modelica.Electrical.Analog.Basic.Resistor res1(R=R1) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-20,20})));
        Modelica.Electrical.Analog.Basic.Capacitor cap(C=C) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-20,-20})));
        Modelica.Electrical.Analog.Basic.Resistor res2(R=R2) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={20,20})));
        Modelica.Electrical.Analog.Basic.Resistor res3(R=R3) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={20,-20})));

      public
        Modelica.Electrical.Analog.Basic.Ground ground
          annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
        Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=90,
              origin={50,-20})));

      equation
        connect(voltageSource.p, res1.p) annotation (Line(
            points={{-50,-10},{-50,40},{-20,40},{-20,30}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(voltageSource.n, ground.p) annotation (Line(
            points={{-50,-30},{-50,-40},{4.996e-16,-40}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(vIn, voltageSource.v) annotation (Line(
            points={{-90,5.55112e-16},{-82,5.55112e-16},{-82,0},{-70,0},{-70,-20},
                {-57,-20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(res3.p, res2.n) annotation (Line(
            points={{20,-10},{20,-5},{20,-5},{20,0},{20,10},{20,10}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(res3.n, ground.p) annotation (Line(
            points={{20,-30},{20,-40},{4.996e-16,-40}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(cap.n, ground.p) annotation (Line(
            points={{-20,-30},{-20,-40},{4.996e-16,-40}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(res1.n, cap.p) annotation (Line(
            points={{-20,10},{-20,5},{-20,5},{-20,0},{-20,-10},{-20,-10}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(res1.n, res2.p) annotation (Line(
            points={{-20,10},{-20,0},{0,0},{0,40},{20,40},{20,30}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(voltageSensor.v, vOut) annotation (Line(
            points={{60,-20},{70,-20},{70,5.55112e-16},{90,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(voltageSensor.p, res2.n) annotation (Line(
            points={{50,-10},{50,0},{20,0},{20,10}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(voltageSensor.n, ground.p) annotation (Line(
            points={{50,-30},{50,-40},{4.996e-16,-40}},
            color={0,0,255},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-60},{100,60}}), graphics));
      end DeclarativeAB;

      model ImperativeA "First circuit in imperative formalism"
        extends BaseClasses.Parameters;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
        Connectors.RealInput vIn annotation (Placement(transformation(extent={{
                  -160,-10},{-140,10}}), iconTransformation(extent={{-120,-10},
                  {-100,10}})));
        Connectors.RealOutput vOut annotation (Placement(transformation(extent=
                  {{0,-10},{20,10}}), iconTransformation(extent={{100,-10},{120,
                  10}})));
        Modelica.Blocks.Math.Gain res1(k=R1)
          annotation (Placement(transformation(extent={{-120,-10},{-100,10}})));
        Modelica.Blocks.Continuous.Derivative ind(k=C)
          annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
        Modelica.Blocks.Math.Add KVL1
          annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));

      equation
        connect(KVL1.y, vOut) annotation (Line(
            points={{-9,6.10623e-16},{-10,6.10623e-16},{-10,5.55112e-16},{10,
                5.55112e-16}},
            color={255,195,38},
            smooth=Smooth.None));

        connect(ind.y, KVL1.u1) annotation (Line(
            points={{-59,6.10623e-16},{-50,6.10623e-16},{-50,6},{-32,6}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(res1.y, ind.u) annotation (Line(
            points={{-99,6.10623e-16},{-94,6.10623e-16},{-94,0},{-90,0},{-90,
                6.66134e-16},{-82,6.66134e-16}},
            color={0,200,0},
            smooth=Smooth.None));
        connect(res1.u, vIn) annotation (Line(
            points={{-122,6.66134e-16},{-134,6.66134e-16},{-134,5.55112e-16},{-150,
                5.55112e-16}},
            color={255,195,38},
            smooth=Smooth.None));

        connect(KVL1.u2, vIn) annotation (Line(
            points={{-32,-6},{-40,-6},{-40,-20},{-130,-20},{-130,5.55112e-16},{
                -150,5.55112e-16}},
            color={255,195,38},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-160,-40},{160,40}}), graphics={Rectangle(
                      extent={{-15,7},{15,-7}},
                      lineColor={0,0,0},
                      origin={-5,27},
                      rotation=180,
                      fillPattern=FillPattern.Solid,
                      fillColor={255,255,255}),Line(
                      points={{-16,30},{-10,30}},
                      color={255,195,38},
                      smooth=Smooth.None),Text(
                      extent={{-10,32},{10,28}},
                      lineColor={0,0,0},
                      textString="Voltage"),Text(
                      extent={{-10,26},{10,22}},
                      lineColor={0,0,0},
                      textString="Current"),Line(
                      points={{-16,24},{-10,24}},
                      color={0,200,0},
                      smooth=Smooth.None)}), Icon(coordinateSystem(
                preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
              graphics));
      end ImperativeA;

      model ImperativeB "Second circuit in imperative formalism"
        extends BaseClasses.Parameters;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
        Connectors.RealInput vIn annotation (Placement(transformation(extent={{
                  -30,-10},{-10,10}}), iconTransformation(extent={{-120,-10},{-100,
                  10}})));
        Connectors.RealOutput vOut annotation (Placement(transformation(extent=
                  {{140,-10},{160,10}}), iconTransformation(extent={{100,-10},{
                  120,10}})));

      public
        Modelica.Blocks.Math.Gain res2(k=1/R2)
          annotation (Placement(transformation(extent={{60,-10},{80,10}})));

      public
        Modelica.Blocks.Math.Gain res3(k=R3)
          annotation (Placement(transformation(extent={{100,-10},{120,10}})));
        Modelica.Blocks.Math.Add KVL2(k2=-1, k1=1)
          annotation (Placement(transformation(extent={{20,-10},{40,10}})));

      equation
        connect(res3.y, vOut) annotation (Line(
            points={{121,6.10623e-16},{136,6.10623e-16},{136,5.55112e-16},{150,
                5.55112e-16}},
            color={255,195,38},
            smooth=Smooth.None));

        connect(res2.y, res3.u) annotation (Line(
            points={{81,6.10623e-16},{86,0},{90,1.27676e-15},{90,6.66134e-16},{
                98,6.66134e-16}},
            color={0,200,0},
            smooth=Smooth.None));

        connect(KVL2.y, res2.u) annotation (Line(
            points={{41,6.10623e-16},{46,0},{50,1.27676e-15},{50,6.66134e-16},{
                58,6.66134e-16}},
            color={255,195,38},
            smooth=Smooth.None));

        connect(KVL2.u1, vIn) annotation (Line(
            points={{18,6},{0,6},{0,5.55112e-16},{-20,5.55112e-16}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(KVL2.u2, res3.y) annotation (Line(
            points={{18,-6},{10,-6},{10,-20},{130,-20},{130,6.10623e-16},{121,
                6.10623e-16}},
            color={255,195,38},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-160,-40},{160,40}}), graphics={Rectangle(
                      extent={{-15,7},{15,-7}},
                      lineColor={0,0,0},
                      origin={135,27},
                      rotation=180,
                      fillPattern=FillPattern.Solid,
                      fillColor={255,255,255}),Line(
                      points={{124,30},{130,30}},
                      color={255,195,38},
                      smooth=Smooth.None),Text(
                      extent={{130,32},{150,28}},
                      lineColor={0,0,0},
                      textString="Voltage"),Text(
                      extent={{130,26},{150,22}},
                      lineColor={0,0,0},
                      textString="Current"),Line(
                      points={{124,24},{130,24}},
                      color={0,200,0},
                      smooth=Smooth.None)}), Icon(coordinateSystem(
                preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
              graphics));
      end ImperativeB;

      model ImperativeAB
        "Cascaded first and second circuits in imperative formalism"
        extends BaseClasses.Parameters;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
        Connectors.RealInput vIn annotation (Placement(transformation(extent={{
                  -110,30},{-90,50}}), iconTransformation(extent={{-120,-10},{-100,
                  10}})));
        Connectors.RealOutput vOut annotation (Placement(transformation(extent=
                  {{80,30},{100,50}}), iconTransformation(extent={{100,-10},{
                  120,10}})));
        Modelica.Blocks.Continuous.Integrator cap(k=1/C)
          annotation (Placement(transformation(extent={{-40,-10},{-60,10}})));
        Modelica.Blocks.Math.Gain res1(k=1/R1)
          annotation (Placement(transformation(extent={{-20,30},{0,50}})));
        Modelica.Blocks.Math.Gain res2(k=1/R2)
          annotation (Placement(transformation(extent={{-20,-50},{0,-30}})));
        Modelica.Blocks.Math.Gain res3(k=R3)
          annotation (Placement(transformation(extent={{40,-50},{60,-30}})));
        Modelica.Blocks.Math.Add KVL1(k2=-1, k1=1)
          annotation (Placement(transformation(extent={{-60,30},{-40,50}})));
        Modelica.Blocks.Math.Add KCL(k2=-1, k1=1)
          annotation (Placement(transformation(extent={{0,-10},{-20,10}})));
        Modelica.Blocks.Math.Add KVL2(k2=-1)
          annotation (Placement(transformation(extent={{-60,-50},{-40,-30}})));

      equation
        connect(res3.y, vOut) annotation (Line(
            points={{61,-40},{70,-40},{70,40},{90,40}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(KVL1.u1, vIn) annotation (Line(
            points={{-62,46},{-80,46},{-80,40},{-100,40}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(cap.y, KVL1.u2) annotation (Line(
            points={{-61,6.10623e-16},{-70,6.10623e-16},{-70,34},{-62,34}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(KCL.u1, res1.y) annotation (Line(
            points={{2,6},{10,6},{10,40},{1,40}},
            color={0,200,0},
            smooth=Smooth.None));
        connect(KCL.y, cap.u) annotation (Line(
            points={{-21,6.10623e-16},{-26,0},{-30,1.27676e-15},{-30,
                6.66134e-16},{-38,6.66134e-16}},
            color={0,200,0},
            smooth=Smooth.None));
        connect(KVL1.y, res1.u) annotation (Line(
            points={{-39,40},{-22,40}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(KVL2.y, res2.u) annotation (Line(
            points={{-39,-40},{-22,-40}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(res2.y, res3.u) annotation (Line(
            points={{1,-40},{38,-40}},
            color={0,200,0},
            smooth=Smooth.None));
        connect(KCL.u2, res2.y) annotation (Line(
            points={{2,-6},{30,-6},{30,-40},{1,-40}},
            color={0,200,0},
            smooth=Smooth.None));
        connect(KVL2.u1, cap.y) annotation (Line(
            points={{-62,-34},{-70,-34},{-70,6.10623e-16},{-61,6.10623e-16}},
            color={255,195,38},
            smooth=Smooth.None));

        connect(KVL2.u2, res3.y) annotation (Line(
            points={{-62,-46},{-70,-46},{-70,-60},{70,-60},{70,-40},{61,-40}},
            color={255,195,38},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-120,-80},{120,80}}), graphics={Rectangle(
                      extent={{-15,7},{15,-7}},
                      lineColor={0,0,0},
                      origin={95,-67},
                      rotation=180,
                      fillPattern=FillPattern.Solid,
                      fillColor={255,255,255}),Line(
                      points={{84,-64},{90,-64}},
                      color={255,195,38},
                      smooth=Smooth.None),Text(
                      extent={{90,-62},{110,-66}},
                      lineColor={0,0,0},
                      textString="Voltage"),Text(
                      extent={{90,-68},{110,-72}},
                      lineColor={0,0,0},
                      textString="Current"),Line(
                      points={{84,-70},{90,-70}},
                      color={0,200,0},
                      smooth=Smooth.None)}));
      end ImperativeAB;

      model ImperativeABIncorrect
        "Incorrectly cascaded first and second circuits in imperative formalism"
        extends BaseClasses.Parameters;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
        Connectors.RealInput vIn annotation (Placement(transformation(extent={{
                  -160,-10},{-140,10}}), iconTransformation(extent={{-120,-10},
                  {-100,10}})));
        Connectors.RealOutput vOut annotation (Placement(transformation(extent=
                  {{140,-10},{160,10}}), iconTransformation(extent={{100,-10},{
                  120,10}})));
        Modelica.Blocks.Math.Gain res1(k=R1)
          annotation (Placement(transformation(extent={{-120,-10},{-100,10}})));
        Modelica.Blocks.Continuous.Derivative ind(k=C)
          annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
        Modelica.Blocks.Math.Add KVL1
          annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));

      public
        Modelica.Blocks.Math.Gain res2(k=1/R2)
          annotation (Placement(transformation(extent={{60,-10},{80,10}})));

      public
        Modelica.Blocks.Math.Gain res3(k=R3)
          annotation (Placement(transformation(extent={{100,-10},{120,10}})));
        Modelica.Blocks.Math.Add KVL2(k2=-1, k1=1)
          annotation (Placement(transformation(extent={{20,-10},{40,10}})));

      equation
        connect(ind.y, KVL1.u1) annotation (Line(
            points={{-59,6.10623e-16},{-50,6.10623e-16},{-50,6},{-32,6}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(res1.y, ind.u) annotation (Line(
            points={{-99,6.10623e-16},{-94,0},{-90,1.27676e-15},{-90,
                6.66134e-16},{-82,6.66134e-16}},
            color={0,200,0},
            smooth=Smooth.None));
        connect(vIn, res1.u) annotation (Line(
            points={{-150,5.55112e-16},{-144,0},{-136,1.22125e-15},{-136,
                6.66134e-16},{-122,6.66134e-16}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(vIn, KVL1.u2) annotation (Line(
            points={{-150,5.55112e-16},{-130,5.55112e-16},{-130,-20},{-40,-20},
                {-40,-6},{-32,-6}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(res2.y, res3.u) annotation (Line(
            points={{81,6.10623e-16},{86,0},{90,1.27676e-15},{90,6.66134e-16},{
                98,6.66134e-16}},
            color={0,200,0},
            smooth=Smooth.None));

        connect(KVL2.y, res2.u) annotation (Line(
            points={{41,6.10623e-16},{46,0},{50,1.27676e-15},{50,6.66134e-16},{
                58,6.66134e-16}},
            color={0,200,0},
            smooth=Smooth.None));

        connect(KVL2.u2, res3.y) annotation (Line(
            points={{18,-6},{10,-6},{10,-20},{130,-20},{130,6.10623e-16},{121,
                6.10623e-16}},
            color={255,195,38},
            smooth=Smooth.None));

        connect(KVL1.y, KVL2.u1) annotation (Line(
            points={{-9,6.10623e-16},{0,6.10623e-16},{0,6},{18,6}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(res3.y, vOut) annotation (Line(
            points={{121,6.10623e-16},{132.5,6.10623e-16},{132.5,5.55112e-16},{
                150,5.55112e-16}},
            color={255,195,38},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-160,-40},{160,40}}), graphics={Rectangle(
                      extent={{-15,7},{15,-7}},
                      lineColor={0,0,0},
                      origin={135,27},
                      rotation=180,
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Solid),Text(
                      extent={{130,32},{150,28}},
                      lineColor={0,0,0},
                      textString="Voltage"),Line(
                      points={{124,30},{130,30}},
                      color={255,195,38},
                      smooth=Smooth.None),Text(
                      extent={{130,26},{150,22}},
                      lineColor={0,0,0},
                      textString="Current"),Line(
                      points={{124,24},{130,24}},
                      color={0,200,0},
                      smooth=Smooth.None)}));
      end ImperativeABIncorrect;

      model ImperativeABTF
        "Cascaded first and second circuits as a transfer function"
        extends BaseClasses.Parameters;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
        Connectors.RealInput vIn annotation (Placement(transformation(extent={{
                  -110,30},{-90,50}}), iconTransformation(extent={{-120,-10},{-100,
                  10}})));
        Connectors.RealOutput vOut annotation (Placement(transformation(extent=
                  {{80,30},{100,50}}), iconTransformation(extent={{100,-10},{
                  120,10}})));
        Modelica.Blocks.Continuous.TransferFunction transferFunction(b={R3}, a=
              {C*R1*(R2 + R3),R1 + R2 + R3})
          annotation (Placement(transformation(extent={{-20,30},{0,50}})));

      equation
        connect(transferFunction.y, vOut) annotation (Line(
            points={{1,40},{90,40}},
            color={255,195,38},
            smooth=Smooth.None));
        connect(transferFunction.u, vIn) annotation (Line(
            points={{-22,40},{-100,40}},
            color={255,195,38},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-120,-80},{120,80}}), graphics));
      end ImperativeABTF;

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;

        partial model Parameters "Base model with parameters"
          parameter SI.Resistance R1=1 "Resistance at temperature T_ref";
          parameter SI.Resistance R2=1 "Resistance at temperature T_ref";
          parameter SI.Resistance R3=1 "Resistance at temperature T_ref";
          parameter SI.Capacitance C=100e-3 "Capacitance";

        end Parameters;

      end BaseClasses;

    end Cascade;

  end DeclarativeVsImperative;

  model EnvironmentIcon

    Conditions.Environment Environment
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
    annotation (Diagram(graphics));

  end EnvironmentIcon;

  model RouterCrossOver

    Conditions.Router router(crossOver=true)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

  end RouterCrossOver;

  model RouterPassThrough

    Conditions.Router router
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

  end RouterPassThrough;

  model CellIcon

    Assemblies.Cells.Cell Cell
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

  end CellIcon;

  model Logo
    extends FCSys.BaseClasses.Icons.Cell;
    annotation (Icon(graphics={Rectangle(
              extent={{-100,100},{100,65}},
              fillPattern=FillPattern.Solid,
              fillColor={255,255,255},
              pattern=LinePattern.None,
              lineColor={0,0,0})}));

  end Logo;

  model AnFPIcon "Anode flow plate"

    Regions.AnFPs.AnFP AnFP
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

  end AnFPIcon;

  model AnGDLIcon "Anode gas diffusion layer"

    Regions.AnGDLs.AnGDL AnGDL
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

  end AnGDLIcon;

  model AnCLIcon "Anode catalyst layer"

    Regions.AnCLs.AnCL AnCL
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

  end AnCLIcon;

  model PEMIcon "Proton exchange membrane"

    Regions.PEMs.PEM PEM
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

  end PEMIcon;

  model CaCLIcon "Cathode catalyst layer"

    Regions.CaCLs.CaCL CaCL
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

  end CaCLIcon;

  model CaGDLIcon "Cathode gas diffusion layer"

    Regions.CaGDLs.CaGDL CaGDL
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

  end CaGDLIcon;

  model CaFPIcon "Cathode flow plate"

    Regions.CaFPs.CaFP CaFP
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

  end CaFPIcon;

  model RegionIcon "Region"

    Regions.AnFPs.AnFP Region
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
    annotation (Diagram(graphics));

  end RegionIcon;

  model Matrix3D

    FCSys.Subregions.BaseClasses.EmptySubregion subreg000(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{-30,-30},{-10,-10}})));
    FCSys.Subregions.BaseClasses.EmptySubregion subreg001(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{-50,-50},{-30,-30}})));
    FCSys.Subregions.BaseClasses.EmptySubregion subreg010(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{-30,20},{-10,40}})));
    FCSys.Subregions.BaseClasses.EmptySubregion subreg011(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{-50,0},{-30,20}})));
    FCSys.Subregions.BaseClasses.EmptySubregion subreg100(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{30,-30},{50,-10}})));
    FCSys.Subregions.BaseClasses.EmptySubregion subreg101(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{10,-50},{30,-30}})));
    FCSys.Subregions.BaseClasses.EmptySubregion subreg110(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{30,20},{50,40}})));
    FCSys.Subregions.BaseClasses.EmptySubregion subreg111(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{10,0},{30,20}})));

  equation
    connect(subreg010.yNegative, subreg000.yPositive) annotation (Line(
        points={{-20,20},{-20,-10}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg011.yNegative, subreg001.yPositive) annotation (Line(
        points={{-40,-5.55112e-16},{-40,-30}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg011.zNegative, subreg010.zPositive) annotation (Line(
        points={{-35,15},{-25,25}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg111.zNegative, subreg110.zPositive) annotation (Line(
        points={{25,15},{35,25}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg110.yNegative, subreg100.yPositive) annotation (Line(
        points={{40,20},{40,-10}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg110.xNegative, subreg010.xPositive) annotation (Line(
        points={{30,30},{-10,30}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg000.xPositive, subreg100.xNegative) annotation (Line(
        points={{-10,-20},{30,-20}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg001.xPositive, subreg101.xNegative) annotation (Line(
        points={{-30,-40},{10,-40}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg001.zNegative, subreg000.zPositive) annotation (Line(
        points={{-35,-35},{-25,-25}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg101.zNegative, subreg100.zPositive) annotation (Line(
        points={{25,-35},{35,-25}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg111.xNegative, subreg011.xPositive) annotation (Line(
        points={{10,10},{-30,10}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg111.yNegative, subreg101.yPositive) annotation (Line(
        points={{20,-5.55112e-16},{20,-30}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    annotation (Diagram(graphics));
  end Matrix3D;

  model SubregionIcon

    FCSys.Subregions.BaseClasses.EmptySubregion Subregion(
      inclX=true,
      inclZ=true,
      inclY=true)
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

  end SubregionIcon;

  partial model PhaseIcon

    FCSys.Subregions.Phases.BaseClasses.EmptyPhase Phase
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
    annotation (structurallyIncomplete=true);

  end PhaseIcon;

  partial model SpeciesIcon

    Subregions.Species.Species Species
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
    annotation (structurallyIncomplete=true);

  end SpeciesIcon;

  model PhaseBoundaryIcon

    Subregions.PhaseBoundary PhaseBoundary
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

  end PhaseBoundaryIcon;

  model ReactionIcon

    FCSys.Subregions.Reaction2 Reaction
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

  end ReactionIcon;

  model VolumeIcon

    Subregions.Volume Volume
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));

  end VolumeIcon;

  model ConnectorHierarchy
    "Extension and instantiation hierarchy of the connectors"

  protected
    Connectors.PhysicalBusInternal Internal
      annotation (Placement(transformation(extent={{54,8},{74,28}})));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-120,
              -40},{100,60}}), graphics={Line(
              points={{16,-8.5},{16,-30}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Line(
              points={{-32,-8.5},{-32,-30}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Line(
              points={{13,-6},{-32,-30}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Line(
              points={{-8,-6},{-32,-30}},
              color={0,0,0},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Line(
              points={{-58,-8},{-80,-30}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Line(
              points={{-30,-8},{-8,-30}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Line(
              points={{64,-8.5},{64,-30}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Line(
              points={{88,-8.5},{88,-30}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Line(
              points={{88,18},{88,-6}},
              color={191,191,191},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Line(
              points={{64,18},{88,-6}},
              color={191,191,191},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Line(
              points={{-8,-6},{-8,-30}},
              color={0,0,0},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Line(
              points={{14,-8},{-8,-30}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Line(
              points={{-54,-8},{-32,-30}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Line(
              points={{18,-8},{40,-30}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Line(
              points={{-53,-6},{-8,-30}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Rectangle(
              extent={{60,56},{96,36}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),Line(
              points={{-32,18},{-32,-6}},
              color={191,191,191},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Line(
              points={{-34,-8},{-56,-30}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Line(
              points={{-100,24},{-100,-26}},
              color={0,0,0},
              smooth=Smooth.None),Line(
              points={{-102,20},{-100,24},{-98,20}},
              color={0,0,0},
              smooth=Smooth.None),Text(
              extent={{-112,30},{-88,25.8}},
              lineColor={0,0,0},
              textString="Composite"),Text(
              extent={{-114,-28},{-86,-32}},
              lineColor={0,0,0},
              textString="Basic"),Text(
              extent={{-40,30},{-24,26}},
              lineColor={0,0,0},
              textString="Face"),Text(
              extent={{-40,26},{-24,22}},
              lineColor={0,0,0},
              textString="Bus"),Ellipse(
              extent={{-35,21},{-29,15}},
              lineColor={127,127,127},
              fillColor={191,191,191},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),Text(
              extent={{65.5,42},{100,38}},
              lineColor={0,0,0},
              textString="    Instantiation"),Text(
              extent={{66,43.8},{100,48}},
              lineColor={0,0,0},
              textString=" Expansion"),Text(
              extent={{66,54},{100,50}},
              lineColor={0,0,0},
              textString="Extension"),Ellipse(
              extent={{-59,-27},{-53,-33}},
              lineColor={127,127,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),Ellipse(
              extent={{-35,-27},{-29,-33}},
              lineColor={170,0,0},
              fillColor={221,23,47},
              fillPattern=FillPattern.Solid),Ellipse(
              extent={{85,-27},{91,-33}},
              lineColor={127,127,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),Text(
              extent={{-20,6},{4,2}},
              lineColor={0,0,0},
              textString="Inert and"),Text(
              extent={{-20,1.8},{4,-2.2}},
              lineColor={0,0,0},
              textString="Inert Internal"),Ellipse(
              extent={{-11,-3},{-5,-9}},
              lineColor={11,43,197},
              fillColor={47,107,251},
              fillPattern=FillPattern.Solid),Ellipse(
              extent={{-10,-4},{-6,-8}},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              lineColor={11,43,197}),Text(
              extent={{-64,6},{-48,2}},
              lineColor={0,0,0},
              textString="Inert"),Ellipse(
              extent={{-83,-27},{-77,-33}},
              lineColor={127,127,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),Text(
              extent={{-20,-22},{4,-26}},
              lineColor={0,0,0},
              textString="Translational"),Text(
              extent={{-92,-22},{-68,-26}},
              lineColor={0,0,0},
              textString="Dalton"),Ellipse(
              extent={{37,-27},{43,-33}},
              lineColor={127,127,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),Line(
              points={{62,46},{72,46}},
              color={191,191,191},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Line(
              points={{72,40},{62,40}},
              color={0,0,0},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Line(
              points={{3,0},{-7,0}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None},
              origin={69,52},
              rotation=360),Text(
              extent={{-44,-22},{-20,-26}},
              lineColor={0,0,0},
              textString="Diffusion"),Text(
              extent={{-68,-22},{-44,-26}},
              lineColor={0,0,0},
              textString="Material"),Ellipse(
              extent={{-11,-27},{-5,-33}},
              lineColor={127,127,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),Text(
              extent={{-44,-18},{-20,-22}},
              lineColor={0,0,0},
              textString="Thermal"),Text(
              extent={{4,-22},{28,-26}},
              lineColor={0,0,0},
              textString="Advection"),Ellipse(
              extent={{13,-27},{19,-33}},
              lineColor={127,127,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),Ellipse(
              extent={{61,-27},{67,-33}},
              lineColor={127,127,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),Ellipse(
              extent={{85,21},{91,15}},
              lineColor={2,157,21},
              fillColor={38,196,52},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),Ellipse(
              extent={{61,-3},{67,-9}},
              lineColor={239,142,1},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),Ellipse(
              extent={{62.5,-4.5},{65.5,-7.5}},
              fillColor={255,195,38},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),Text(
              extent={{52,2},{76,-2.2}},
              lineColor={0,0,0},
              textString="Species"),Ellipse(
              extent={{85,-3},{91,-9}},
              lineColor={2,157,21},
              fillColor={38,196,52},
              fillPattern=FillPattern.Solid),Text(
              extent={{76,30},{100,25.8}},
              lineColor={0,0,0},
              textString="Physical"),Text(
              extent={{76,2},{100,-2.2}},
              lineColor={0,0,0},
              textString="Physical"),Text(
              extent={{52,28},{76,23.8}},
              lineColor={0,0,0},
              textString="Physical Bus"),Text(
              extent={{52,-22},{76,-26}},
              lineColor={0,0,0},
              textString="Potential"),Text(
              extent={{28,-18},{52,-22}},
              lineColor={0,0,0},
              textString="Net"),Text(
              extent={{76,-22},{100,-26.2}},
              lineColor={0,0,0},
              textString="Activity"),Text(
              extent={{28,-22},{52,-26}},
              lineColor={0,0,0},
              textString="Potential"),Text(
              extent={{80,26},{96,22}},
              lineColor={0,0,0},
              textString="Bus"),Text(
              extent={{48,-18},{80,-22}},
              lineColor={0,0,0},
              textString="Electrochemical"),Text(
              extent={{4,6},{28,2}},
              lineColor={0,0,0},
              textString="Chemical"),Text(
              extent={{4,-18},{28,-22}},
              lineColor={0,0,0},
              textString="Thermal"),Text(
              extent={{52,6},{76,2}},
              lineColor={0,0,0},
              textString="Chemical"),Ellipse(
              extent={{-35,-3},{-29,-9}},
              lineColor={127,127,127},
              fillColor={191,191,191},
              fillPattern=FillPattern.Solid),Ellipse(
              extent={{-59,-3},{-53,-9}},
              lineColor={11,43,197},
              fillColor={47,107,251},
              fillPattern=FillPattern.Solid),Text(
              extent={{-64,2},{-48,-2}},
              lineColor={0,0,0},
              textString="Dalton"),Text(
              extent={{-40,2},{-24,-2}},
              lineColor={0,0,0},
              textString="Face"),Ellipse(
              extent={{13,-3},{19,-9}},
              lineColor={239,142,1},
              fillColor={255,195,38},
              fillPattern=FillPattern.Solid),Text(
              extent={{8,2},{24,-2}},
              lineColor={0,0,0},
              textString="Reaction")}));

  end ConnectorHierarchy;

end Figures;
