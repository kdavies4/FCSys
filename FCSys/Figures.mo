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
            FCSys.Figures.DeclarativeVsImperative.SwitchCausality.Parameters;
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
        extends Parameters;
        extends FCSys.Icons.Blocks.Continuous;
        FCSys.Connectors.RealInput v annotation (Placement(transformation(
                extent={{-60,-30},{-40,-10}}), iconTransformation(extent={{-120,
                  -10},{-100,10}})));
        FCSys.Connectors.RealInput i annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={0,10}), iconTransformation(
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
            points={{-5.55112e-16,10},{-5,10},{-5,10},{-10,10}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(cap.n, ground.p) annotation (Line(
            points={{30,-30},{30,-35},{30,-35},{30,-40}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(i, i) annotation (Line(
            points={{-5.55112e-16,10},{-5.55112e-16,10}},
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
        extends Parameters;
        extends FCSys.Icons.Blocks.Continuous;
        FCSys.Connectors.RealInput i annotation (Placement(transformation(
                extent={{-80,-30},{-60,-10}}), iconTransformation(extent={{-120,
                  -10},{-100,10}})));
        FCSys.Connectors.RealOutput v annotation (Placement(transformation(
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
        extends Parameters;
        extends FCSys.Icons.Blocks.Continuous;
        FCSys.Connectors.RealInput v annotation (Placement(transformation(
                extent={{-100,20},{-80,40}}), iconTransformation(extent={{-120,
                  -10},{-100,10}})));
        FCSys.Connectors.RealOutput i annotation (Placement(transformation(
                extent={{80,20},{100,40}}), iconTransformation(extent={{100,-10},
                  {120,10}})));
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
        extends Parameters;
        extends FCSys.Icons.Blocks.Continuous;

      public
        FCSys.Connectors.RealInput i annotation (Placement(transformation(
                extent={{-90,26},{-70,46}}), iconTransformation(extent={{-120,-10},
                  {-100,10}})));
        FCSys.Connectors.RealOutput v annotation (Placement(transformation(
                extent={{80,20},{100,40}}), iconTransformation(extent={{100,-10},
                  {120,10}})));
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
            points={{-80,36},{-52,36}},
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
        extends Parameters;
        extends FCSys.Icons.Blocks.Continuous;
        FCSys.Connectors.RealInput v annotation (Placement(transformation(
                extent={{-100,20},{-80,40}}), iconTransformation(extent={{-120,
                  -10},{-100,10}})));
        FCSys.Connectors.RealOutput i annotation (Placement(transformation(
                extent={{80,20},{100,40}}), iconTransformation(extent={{100,-10},
                  {120,10}})));
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
        extends Parameters;
        extends FCSys.Icons.Blocks.Continuous;

      public
        FCSys.Connectors.RealInput i annotation (Placement(transformation(
                extent={{-100,20},{-80,40}}), iconTransformation(extent={{-120,
                  -10},{-100,10}})));
        FCSys.Connectors.RealOutput v annotation (Placement(transformation(
                extent={{80,20},{100,40}}), iconTransformation(extent={{100,-10},
                  {120,10}})));
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

    protected
      partial model Parameters
        "Base model with parameters for electrical circuit"
        parameter SI.Resistance R1=10 "Resistance at temperature T_ref";
        parameter SI.Resistance R2=100 "Resistance at temperature T_ref";
        parameter SI.Inductance L=0.1 "Inductance";
        parameter SI.Capacitance C=0.01 "Capacitance";

      end Parameters;
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
        extends Parameters;
        extends FCSys.Icons.Blocks.Continuous;
        FCSys.Connectors.RealInput vIn annotation (Placement(transformation(
                extent={{-100,-10},{-80,10}}), iconTransformation(extent={{-120,
                  -10},{-100,10}})));
        FCSys.Connectors.RealOutput vOut annotation (Placement(transformation(
                extent={{40,-10},{60,10}}), iconTransformation(extent={{100,-10},
                  {120,10}})));
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
        extends Parameters;
        extends FCSys.Icons.Blocks.Continuous;
        FCSys.Connectors.RealInput vIn annotation (Placement(transformation(
                extent={{-60,-10},{-40,10}}), iconTransformation(extent={{-120,
                  -10},{-100,10}})));
        FCSys.Connectors.RealOutput vOut annotation (Placement(transformation(
                extent={{80,-10},{100,10}}), iconTransformation(extent={{100,-10},
                  {120,10}})));
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
        extends Parameters;
        extends FCSys.Icons.Blocks.Continuous;
        FCSys.Connectors.RealInput vIn annotation (Placement(transformation(
                extent={{-100,-10},{-80,10}}), iconTransformation(extent={{-120,
                  -10},{-100,10}})));
        FCSys.Connectors.RealOutput vOut annotation (Placement(transformation(
                extent={{80,-10},{100,10}}), iconTransformation(extent={{100,-10},
                  {120,10}})));
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
        extends Parameters;
        extends FCSys.Icons.Blocks.Continuous;
        FCSys.Connectors.RealInput vIn annotation (Placement(transformation(
                extent={{-160,-10},{-140,10}}), iconTransformation(extent={{-120,
                  -10},{-100,10}})));
        FCSys.Connectors.RealOutput vOut annotation (Placement(transformation(
                extent={{0,-10},{20,10}}), iconTransformation(extent={{100,-10},
                  {120,10}})));
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
        extends Parameters;
        extends FCSys.Icons.Blocks.Continuous;
        FCSys.Connectors.RealInput vIn annotation (Placement(transformation(
                extent={{-30,-10},{-10,10}}), iconTransformation(extent={{-120,
                  -10},{-100,10}})));
        FCSys.Connectors.RealOutput vOut annotation (Placement(transformation(
                extent={{140,-10},{160,10}}), iconTransformation(extent={{100,-10},
                  {120,10}})));

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
        extends Parameters;
        extends FCSys.Icons.Blocks.Continuous;
        FCSys.Connectors.RealInput vIn annotation (Placement(transformation(
                extent={{-110,30},{-90,50}}), iconTransformation(extent={{-120,
                  -10},{-100,10}})));
        FCSys.Connectors.RealOutput vOut annotation (Placement(transformation(
                extent={{80,30},{100,50}}), iconTransformation(extent={{100,-10},
                  {120,10}})));
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
        extends Parameters;
        extends FCSys.Icons.Blocks.Continuous;
        FCSys.Connectors.RealInput vIn annotation (Placement(transformation(
                extent={{-160,-10},{-140,10}}), iconTransformation(extent={{-120,
                  -10},{-100,10}})));
        FCSys.Connectors.RealOutput vOut annotation (Placement(transformation(
                extent={{140,-10},{160,10}}), iconTransformation(extent={{100,-10},
                  {120,10}})));
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
        extends Parameters;
        extends FCSys.Icons.Blocks.Continuous;
        FCSys.Connectors.RealInput vIn annotation (Placement(transformation(
                extent={{-110,30},{-90,50}}), iconTransformation(extent={{-120,
                  -10},{-100,10}})));
        FCSys.Connectors.RealOutput vOut annotation (Placement(transformation(
                extent={{80,30},{100,50}}), iconTransformation(extent={{100,-10},
                  {120,10}})));
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

    protected
      partial model Parameters "Base model with parameters"
        parameter SI.Resistance R1=1 "Resistance at temperature T_ref";
        parameter SI.Resistance R2=1 "Resistance at temperature T_ref";
        parameter SI.Resistance R3=1 "Resistance at temperature T_ref";
        parameter SI.Capacitance C=100e-3 "Capacitance";

      end Parameters;
    end Cascade;

  end DeclarativeVsImperative;

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
    annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-120,
              -120},{120,120}})));
  end CellIcon;

  model Logo
    extends FCSys.Icons.Cell;
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
    annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-120,
              -120},{120,120}})));

  end RegionIcon;

  model Matrix

    FCSys.Subregions.Subregion subreg000(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{-30,-30},{-10,-10}})));
    FCSys.Subregions.Subregion subreg001(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{-50,-50},{-30,-30}})));
    FCSys.Subregions.Subregion subreg010(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{-30,20},{-10,40}})));
    FCSys.Subregions.Subregion subreg011(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{-50,0},{-30,20}})));
    FCSys.Subregions.Subregion subreg100(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{30,-30},{50,-10}})));
    FCSys.Subregions.Subregion subreg101(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{10,-50},{30,-30}})));
    FCSys.Subregions.Subregion subreg110(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{30,20},{50,40}})));
    FCSys.Subregions.Subregion subreg111(
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
  end Matrix;

  model SubregionIcon

    Subregions.Subregion Subregion
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-120,
              -120},{120,120}})));
  end SubregionIcon;

  partial model PhaseIcon

    Phases.Gas Phase
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
    annotation (structurallyIncomplete=true, Diagram(coordinateSystem(
            preserveAspectRatio=true, extent={{-120,-120},{120,120}})));

  end PhaseIcon;

  partial model SpeciesIcon

    FCSys.Species.Fluid Species
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
    annotation (structurallyIncomplete=true, Diagram(coordinateSystem(
            preserveAspectRatio=true, extent={{-120,-120},{120,120}})));

  end SpeciesIcon;

  model Connectors "Extension and instantiation hierarchy of the connectors"

    import FCSys.Connectors;

  protected
    Connectors.Amagat Amagat
      annotation (Placement(transformation(extent={{-88,-38},{-68,-18}})));

    Connectors.Dalton Dalton
      annotation (Placement(transformation(extent={{-68,-38},{-48,-18}})));
    Connectors.Reaction Reaction
      annotation (Placement(transformation(extent={{16,-14},{36,6}})));
    Connectors.Direct Direct
      annotation (Placement(transformation(extent={{40,-14},{60,6}})));
    Connectors.Translational Material
      annotation (Placement(transformation(extent={{88,-38},{108,-18}})));

    Connectors.Face Face
      annotation (Placement(transformation(extent={{88,-14},{108,6}})));

    Connectors.Translational Translational
      annotation (Placement(transformation(extent={{40,-38},{60,-18}})));

    annotation (Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-120,
              -40},{120,40}}), graphics={Line(
              points={{74,-6.5},{74,-29.9771},{75.815,-26.369},{72,-28}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Line(
              points={{96,-6},{74,-28}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Line(
              points={{26,-6.5},{26,-29.9771},{27.8145,-26.369},{24,-28}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Line(
              points={{50,-4},{74,-28}},
              color={0,0,0},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Line(
              points={{95,-5},{50,-28}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Line(
              points={{50,-4},{50,-28}},
              color={0,0,0},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Line(
              points={{28,-6},{50,-28}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Line(
              points={{24,-6},{2,-28}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Rectangle(
              extent={{-81,32},{-45,12}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),Line(
              points={{98,20},{98,-4}},
              color={191,191,191},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Line(
              points={{-100,-24},{-100,26}},
              color={0,0,0},
              smooth=Smooth.None,
              arrow={Arrow.None,Arrow.Filled}),Text(
              extent={{-112,32},{-88,27.8}},
              lineColor={0,0,0},
              textString="Composite"),Text(
              extent={{-114,-26},{-86,-30}},
              lineColor={0,0,0},
              textString="Basic"),Ellipse(
              extent={{95,23},{101,17}},
              lineColor={127,127,127},
              fillColor={191,191,191},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),Text(
              extent={{-75.5,18},{-41,14}},
              lineColor={0,0,0},
              textString="    Instantiation"),Text(
              extent={{-75,19.8},{-41,24}},
              lineColor={0,0,0},
              textString=" Expansion"),Text(
              extent={{-75,30},{-41,26}},
              lineColor={0,0,0},
              textString="Extension"),Ellipse(
              extent={{-1,-25},{5,-31}},
              lineColor={127,127,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),Line(
              points={{-79,22},{-69,22}},
              color={191,191,191},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Line(
              points={{-69,16},{-79,16}},
              color={0,0,0},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Ellipse(
              extent={{23,-25},{29,-31}},
              lineColor={127,127,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),Ellipse(
              extent={{71,-25},{77,-31}},
              lineColor={127,127,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),Ellipse(
              extent={{71,-1},{77,-7}},
              lineColor={2,157,21},
              fillColor={38,196,52},
              fillPattern=FillPattern.Solid),Line(
              points={{98,-6.5},{98,-29.9771},{99.815,-26.369},{96,-28}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Line(
              points={{72,-6},{50,-28}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Text(
              extent={{88,30},{108,26}},
              lineColor={0,0,0},
              textString="Face
Bus"),Text(   extent={{62,-18},{86,-22}},
              lineColor={0,0,0},
              textString="Thermal
Diffusive"),Text(
              extent={{14,-18},{38,-22}},
              lineColor={0,0,0},
              textString="Thermal
Advective"),Text(
              extent={{-10,-18},{14,-22}},
              lineColor={0,0,0},
              textString="Stoichio-
metric"),Text(extent={{62,6},{86,1.8}},
              lineColor={0,0,0},
              textString="Inter, Intra, 
and Inert"),Line(
              points={{-69,28},{-79,28},{-69,28},{-79,28}},
              color={191,191,191},
              smooth=Smooth.None,
              arrow={Arrow.Filled,Arrow.None}),Text(
              extent={{-30,-18},{-6,-22}},
              lineColor={0,0,0},
              textString="Electro-
chemical"),Text(
              extent={{-50,-18},{-26,-22}},
              lineColor={0,0,0},
              textString="Electro-
static"),Ellipse(
              extent={{-41,-25},{-35,-31}},
              lineColor={239,142,1},
              fillColor={255,195,38},
              fillPattern=FillPattern.Solid),Ellipse(
              extent={{-21,-25},{-15,-31}},
              lineColor={170,0,0},
              fillColor={221,23,47},
              fillPattern=FillPattern.Solid)}));

  end Connectors;

  model FPtoFP "Test one flow plate to the other"

    extends Modelica.Icons.Example;
    parameter Q.NumberAbsolute n_O2=0.21;
    parameter Q.NumberAbsolute anStoich=1.5 "Anode stoichiometric ratio";
    parameter Q.NumberAbsolute caStoich=2.0 "Cathode stoichiometric ratio";
    parameter Q.NumberAbsolute anInletRH=0.8;
    parameter Q.NumberAbsolute caInletRH=0.5;
    parameter Real T_degC=60;
    parameter Real p_kPag=48.3;
    final parameter Q.TemperatureAbsolute T=U.from_degC(T_degC);
    final parameter Q.PressureAbsolute p=U.from_kPag(p_kPag);
    final parameter Q.PressureAbsolute p_H2O_an_in=anInletRH*
        Modelica.Media.Air.MoistAir.saturationPressureLiquid(T/U.K)*U.Pa
      "Pressure of H2O vapor";
    final parameter Q.PressureAbsolute p_H2O_ca_in=caInletRH*
        Modelica.Media.Air.MoistAir.saturationPressureLiquid(T/U.K)*U.Pa
      "Pressure of H2O vapor";

    output Q.Potential w[n_y, n_z]=(anBC.graphite.'e-'.face.mPhidot[1] ./ anBC.graphite.
        'e-'.face.rho - caBC.graphite.'e-'.face.mPhidot[1] ./ caBC.graphite.
        'e-'.face.rho)/PEM.A[Axis.x] "Potential";
    output Q.CurrentAreic zJ[n_y, n_z]=-anBC.graphite.'e-'.face.phi[1] .* anBC.graphite.
        'e-'.face.rho "Current density";
    output Q.Number zJ_Apercm2[n_y, n_z]=zJ*U.cm^2/U.A
      "Electrical current density, in A/cm2";
    output Q.Current zI=sum(zJ .* outerProduct(L_y, L_z))
      "Total electrical current";
    Q.Current I_H2_an_in=max(0*50*0.18*U.A, zI*anStoich);
    Q.Current I_O2_ca_in=max(0*50*0.18*U.A, zI*caStoich);

    parameter Q.Length L_y[:]={U.m} "Lengths along the channel" annotation (
        Dialog(group="Geometry", __Dymola_label=
            "<html><i>L</i><sub>y</sub></html>"));
    parameter Q.Length L_z[:]={5*U.mm} "Lengths across the channel" annotation
      (Dialog(group="Geometry", __Dymola_label=
            "<html><i>L</i><sub>z</sub></html>"));
    final parameter Integer n_y=size(L_y, 1)
      "Number of regions along the channel" annotation (HideResult=true);
    final parameter Integer n_z=size(L_z, 1)
      "Number of regions across the channel" annotation (HideResult=true);
    final parameter Q.Area A_an_seg[anFP.n_x, n_z]=outerProduct(anFP.L_x, L_z)
      "Areas of the segments over the xz plane of the anode flow plate";
    final parameter Q.Area A_ca_seg[caFP.n_x, n_z]=outerProduct(caFP.L_x, L_z)
      "Areas of the segments over the xz plane of the cathode flow plate";
    final parameter Q.Area A_an=sum(anFP.L_x)*sum(L_z)
      "Total cross-sectional area of the anode flow plate in the xz plane";
    final parameter Q.Area A_ca=sum(caFP.L_x)*sum(L_z)
      "Total cross-sectional area of the cathode flow plate in the xz plane";

    Regions.AnFPs.AnFP anFP(
      final L_y=L_y,
      final L_z=L_z,
      Subregion(
        gas(H2(p_IC=p - p_H2O_an_in), H2O(T_IC=T, p_IC=p_H2O_an_in)),
        liquid(H2O(T_IC=T)),
        graphite('C+'(T_IC=T))))
      annotation (Placement(transformation(extent={{-70,30},{-50,50}})));

    Regions.AnGDLs.AnGDL anGDL(
      final L_y=L_y,
      final L_z=L_z,
      Subregion(
        gas(H2(p_IC=p - p_H2O_an_in), H2O(T_IC=T, p_IC=p_H2O_an_in)),
        liquid(H2O(T_IC=T)),
        graphite('C+'(T_IC=T))))
      annotation (Placement(transformation(extent={{-50,30},{-30,50}})));

    Regions.AnCLs.AnCL anCL(
      final L_y=L_y,
      final L_z=L_z,
      Subregion(
        ionomer('SO3-'(T_IC=T)),
        gas(H2(p_IC=p - p_H2O_an_in), H2O(T_IC=T, p_IC=p_H2O_an_in)),
        liquid(H2O(T_IC=T)),
        graphite('C+'(T_IC=T))))
      annotation (Placement(transformation(extent={{-30,30},{-10,50}})));

    Regions.PEMs.PEM PEM(
      final L_y=L_y,
      final L_z=L_z,
      Subregion(ionomer('SO3-'(T_IC=T))))
      annotation (Placement(transformation(extent={{-10,30},{10,50}})));

    Regions.CaCLs.CaCL caCL(
      final L_y=L_y,
      final L_z=L_z,
      Subregion(
        ionomer('SO3-'(T_IC=T)),
        gas(
          H2O(T_IC=T, p_IC=p_H2O_ca_in),
          N2(p_IC=(p - p_H2O_ca_in)*(1 - environment.n_O2)),
          O2(p_IC=(p - p_H2O_ca_in)*environment.n_O2)),
        liquid(H2O(T_IC=T)),
        graphite('C+'(T_IC=T))))
      annotation (Placement(transformation(extent={{10,30},{30,50}})));

    Regions.CaGDLs.CaGDL caGDL(
      final L_y=L_y,
      final L_z=L_z,
      Subregion(
        gas(
          H2O(T_IC=T, p_IC=p_H2O_ca_in),
          N2(p_IC=(p - p_H2O_ca_in)*(1 - environment.n_O2)),
          O2(p_IC=(p - p_H2O_ca_in)*environment.n_O2)),
        liquid(H2O(T_IC=T)),
        graphite('C+'(T_IC=T))))
      annotation (Placement(transformation(extent={{30,30},{50,50}})));

    Regions.CaFPs.CaFP caFP(
      final L_y=L_y,
      final L_z=L_z,
      Subregion(
        gas(
          H2O(T_IC=T, p_IC=p_H2O_ca_in),
          N2(p_IC=(p - p_H2O_ca_in)*(1 - environment.n_O2)),
          O2(p_IC=(p - p_H2O_ca_in)*environment.n_O2)),
        liquid(H2O(T_IC=T)),
        graphite('C+'(T_IC=T))))
      annotation (Placement(transformation(extent={{50,30},{70,50}})));

    Conditions.ByConnector.FaceBus.Single.Flows anBC[n_y, n_z](each graphite(
        'inclC+'=true,
        'incle-'=true,
        'C+'(redeclare function thermalSpec =
              FCSys.Conditions.ByConnector.Face.Single.Thermal.temperature,
            thermalSet(y=T)),
        'e-'(redeclare function thermalSpec =
              FCSys.Conditions.ByConnector.Face.Single.Thermal.temperature,
            thermalSet(y=T)))) annotation (Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=270,
          origin={-84,40})));

    Conditions.ByConnector.FaceBus.Single.Flows caBC[n_y, n_z](each graphite(
        'inclC+'=true,
        'incle-'=true,
        'C+'(redeclare function thermalSpec =
              FCSys.Conditions.ByConnector.Face.Single.Thermal.temperature,
            thermalSet(y=T)),
        'e-'(
          redeclare function normalSpec =
              Conditions.ByConnector.Face.Single.TranslationalNormal.currentDensity,

          redeclare Modelica.Blocks.Sources.Ramp normalSet(
            height=-2.25*U.A/U.cm^2,
            duration=225.1,
            startTime=10.1),
          redeclare function thermalSpec =
              FCSys.Conditions.ByConnector.Face.Single.Thermal.temperature,
          thermalSet(y=T)))) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={84,40})));

    Conditions.ByConnector.FaceBus.Single.Efforts anSource[anFP.n_x, n_z](gas(
        each inclH2=true,
        each inclH2O=true,
        H2(
          redeclare each function materialSpec =
              FCSys.Conditions.ByConnector.Face.Single.Material.current,
          materialSet(y=(anSource.gas.H2.normalOut.y - fill(
                  I_H2_an_in/(2*A_an),
                  anFP.n_x,
                  n_z)) .* A_an_seg),
          redeclare each function normalMeas =
              FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.currentDensity,

          redeclare each function normalSpec =
              FCSys.Conditions.ByConnector.Face.Single.Translational.force,
          each thermalSet(y=T)),
        H2O(
          redeclare each function materialSpec =
              FCSys.Conditions.ByConnector.Face.Single.Material.current,
          materialSet(y=(anSource.gas.H2O.normalOut.y - fill(
                  I_H2_an_in*p_H2O_an_in/(2*(p - p_H2O_an_in)*A_an),
                  anFP.n_x,
                  n_z)) .* A_an_seg),
          redeclare each function normalMeas =
              FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.currentDensity,

          redeclare each function normalSpec =
              FCSys.Conditions.ByConnector.Face.Single.Translational.force,
          each thermalSet(y=T)))) annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-60,16})));

    Conditions.ByConnector.FaceBus.Single.Flows anSink[anFP.n_x, n_z](gas(
        each inclH2=true,
        each inclH2O=true,
        H2(redeclare each function materialMeas =
              FCSys.Conditions.ByConnector.Face.Single.Material.pressure (
                redeclare package Data = FCSys.Characteristics.IdealGas),
            normalSet(y=p_an_out_noneq*A_an_seg .* anSink.gas.H2.face.rho ./ (
                anSink.gas.H2.face.rho + anSink.gas.H2O.face.rho))),
        H2O(redeclare each function materialMeas =
              FCSys.Conditions.ByConnector.Face.Single.Material.pressure (
                redeclare package Data = FCSys.Characteristics.IdealGas),
            normalSet(y=p_an_out_noneq*A_an_seg .* anSink.gas.H2O.face.rho ./ (
                anSink.gas.H2.face.rho + anSink.gas.H2O.face.rho)))))
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-60,64})));

    Conditions.ByConnector.FaceBus.Single.Efforts caSource[caFP.n_x, n_z](gas(
        each inclH2O=true,
        each inclN2=true,
        each inclO2=true,
        H2O(
          redeclare each function materialSpec =
              FCSys.Conditions.ByConnector.Face.Single.Material.current,
          materialSet(y=(caSource.gas.H2O.normalOut.y - fill(
                  I_O2_ca_in*p_H2O_ca_in/(4*(p - p_H2O_ca_in)*n_O2*A_ca),
                  caFP.n_x,
                  n_z)) .* A_ca_seg),
          redeclare each function normalMeas =
              FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.currentDensity,

          redeclare each function normalSpec =
              FCSys.Conditions.ByConnector.Face.Single.Translational.force,
          each thermalSet(y=T)),
        N2(
          redeclare each function materialSpec =
              FCSys.Conditions.ByConnector.Face.Single.Material.current,
          materialSet(y=(caSource.gas.N2.normalOut.y - fill(
                  I_O2_ca_in*(1 - n_O2)/(4*n_O2*A_ca),
                  caFP.n_x,
                  n_z)) .* A_ca_seg),
          redeclare each function normalMeas =
              FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.currentDensity,

          redeclare each function normalSpec =
              FCSys.Conditions.ByConnector.Face.Single.Translational.force,
          each thermalSet(y=T)),
        O2(
          redeclare each function materialSpec =
              FCSys.Conditions.ByConnector.Face.Single.Material.current,
          materialSet(y=(caSource.gas.O2.normalOut.y - fill(
                  I_O2_ca_in/(4*A_ca),
                  caFP.n_x,
                  n_z)) .* A_ca_seg),
          redeclare each function normalMeas =
              FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.currentDensity,

          redeclare each function normalSpec =
              FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.force,

          each thermalSet(y=T)))) annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={60,16})));

    Conditions.ByConnector.FaceBus.Single.Flows caSink[caFP.n_x, n_z](gas(
        each inclH2O=true,
        each inclN2=true,
        each inclO2=true,
        H2O(redeclare each function materialMeas =
              FCSys.Conditions.ByConnector.Face.Single.Material.pressure (
                redeclare package Data = FCSys.Characteristics.IdealGas),
            normalSet(y=p_ca_out_noneq*A_ca_seg .* caSink.gas.H2O.face.rho ./ (
                caSink.gas.H2O.face.rho + caSink.gas.N2.face.rho + caSink.gas.O2.face.rho))),

        N2(redeclare each function materialMeas =
              FCSys.Conditions.ByConnector.Face.Single.Material.pressure (
                redeclare package Data = FCSys.Characteristics.IdealGas),
            normalSet(y=p_ca_out_noneq*A_ca_seg .* caSink.gas.N2.face.rho ./ (
                caSink.gas.H2O.face.rho + caSink.gas.N2.face.rho + caSink.gas.O2.face.rho))),

        O2(redeclare each function materialMeas =
              FCSys.Conditions.ByConnector.Face.Single.Material.pressure (
                redeclare package Data = FCSys.Characteristics.IdealGas),
            normalSet(y=p_ca_out_noneq*A_ca_seg .* caSink.gas.O2.face.rho ./ (
                caSink.gas.H2O.face.rho + caSink.gas.N2.face.rho + caSink.gas.O2.face.rho)))))
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={60,64})));

    inner Conditions.Environment environment(analysis=true, p=48.3*U.kPa + U.atm)
      "Environmental conditions"
      annotation (Placement(transformation(extent={{-10,70},{10,90}})));

  equation
    connect(anFP.xPositive, anGDL.xNegative) annotation (Line(
        points={{-50,40},{-50,40}},
        color={240,0,0},
        thickness=0.5,
        smooth=Smooth.None));
    connect(anGDL.xPositive, anCL.xNegative) annotation (Line(
        points={{-30,40},{-30,40}},
        color={240,0,0},
        smooth=Smooth.None,
        thickness=0.5));

    connect(anCL.xPositive, PEM.xNegative) annotation (Line(
        points={{-10,40},{-10,40}},
        color={240,0,0},
        smooth=Smooth.None,
        thickness=0.5));

    connect(PEM.xPositive, caCL.xNegative) annotation (Line(
        points={{10,40},{10,40}},
        color={0,0,240},
        smooth=Smooth.None,
        thickness=0.5));

    connect(caCL.xPositive, caGDL.xNegative) annotation (Line(
        points={{30,40},{30,40}},
        color={0,0,240},
        smooth=Smooth.None,
        thickness=0.5));

    connect(caGDL.xPositive, caFP.xNegative) annotation (Line(
        points={{50,40},{50,40}},
        color={0,0,240},
        thickness=0.5,
        smooth=Smooth.None));

    connect(anSource.face, anFP.yNegative) annotation (Line(
        points={{-60,20},{-60,30}},
        color={240,0,0},
        thickness=0.5,
        smooth=Smooth.None));
    connect(anSink.face, anFP.yPositive) annotation (Line(
        points={{-60,60},{-60,50}},
        color={240,0,0},
        thickness=0.5,
        smooth=Smooth.None));
    connect(caSource.face, caFP.yNegative) annotation (Line(
        points={{60,20},{60,30}},
        color={0,0,240},
        thickness=0.5,
        smooth=Smooth.None));
    connect(caSink.face, caFP.yPositive) annotation (Line(
        points={{60,60},{60,50}},
        color={0,0,240},
        thickness=0.5,
        smooth=Smooth.None));
    connect(anBC.face, anFP.xNegative) annotation (Line(
        points={{-80,40},{-70,40}},
        color={127,127,127},
        thickness=0.5,
        smooth=Smooth.None));

    connect(caFP.xPositive, caBC.face) annotation (Line(
        points={{70,40},{80,40}},
        color={127,127,127},
        thickness=0.5,
        smooth=Smooth.None));

    annotation (
      Commands(file="Resources/Scripts/Dymola/Regions.Examples.FPtoFP.mos"
          "Regions.Examples.FPtoFP.mos"),
      experiment(
        StopTime=230,
        Tolerance=1e-05,
        Algorithm="Dassl"),
      Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
              100,100}}), graphics),
      __Dymola_experimentSetupOutput);
  end FPtoFP;

  model InternalFlow
    extends Subregions.Examples.InternalFlow(final Vdot_large=Vdot_large_SI*U.m
          ^3/U.s);

    // Conditions
    parameter Q.VolumeRate Vdot_large_SI=1e-5
      "Prescribed large-signal volumetric flow rate in m3/s"
      annotation (Dialog(__Dymola_label="<html><i>V&#775;</i></html>"));

  end InternalFlow;
end Figures;
