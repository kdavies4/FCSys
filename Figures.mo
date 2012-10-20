within FCSys;
package Figures "Layouts for documentation"
  extends Modelica.Icons.Package;
  package DeclarativeVsImperative
    extends Modelica.Icons.Package;
    // See notes from 12/7/10 for the derivation of the equations in this
    // package.
    model DeclarativeResistor
      parameter SI.Resistance r=1;
      SI.Voltage Deltav;
      SI.Current i;
      ElecPin n;
      ElecPin p;
    protected
      connector ElecPin
        SI.Voltage v;
        flow SI.Current i;
      end ElecPin;
    equation
      // Voltage across the resistor
      Deltav = p.v - n.v;
      // Current through the resistor
      i = (n.i - p.i)/2;
      // No charge storage
      0 = n.i + p.i;
      // Ohm's law
      Deltav = i*r;
    end DeclarativeResistor;

    function ImperativeResistor
      extends Modelica.Icons.Function;
      parameter SI.Resistance r=1;
      input SI.Current i;
      output SI.Voltage Deltav;
    algorithm
      // Ohm's law
      Deltav := i*r;
    end ImperativeResistor;

    package SwitchCausality
      extends Modelica.Icons.Package;
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
        annotation (
          Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                  {100,100}}), graphics),
          experiment(StopTime=2),
          experimentSetupOutput);
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
        annotation (
          Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                  {100,100}}), graphics),
          experiment(StopTime=2),
          experimentSetupOutput);
      end Test_iv;

      model Declarative
        extends Params;
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
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-80,-60},{80,60}}), graphics));
      end Declarative;

      model Declarative_vi
        extends Params;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
        FCSys.Connectors.RealInput v annotation (Placement(transformation(
                extent={{-60,-30},{-40,-10}}), iconTransformation(extent={{-120,
                  -10},{-100,10}})));
        FCSys.Connectors.RealOutput i annotation (Placement(transformation(
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
                  {-80,-60},{80,60}}), graphics), Icon(graphics));
      end Declarative_vi;

      model Declarative_iv
        extends Params;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
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

      model Imperative_vi
        extends Params;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
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
            color={200,0,0},
            smooth=Smooth.None));
        connect(sum.y, i) annotation (Line(
            points={{71,30},{90,30}},
            color={0,200,0},
            smooth=Smooth.None));
        connect(diff2.y, ind.u) annotation (Line(
            points={{-39,6.10623e-16},{-34,0},{-30,1.27676e-15},{-30,
                6.66134e-16},{-22,6.66134e-16}},
            color={200,0,0},
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
            color={200,0,0},
            smooth=Smooth.None));
        connect(diff2.u1, v) annotation (Line(
            points={{-62,6},{-70,6},{-70,30},{-90,30}},
            color={200,0,0},
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
            color={200,0,0},
            smooth=Smooth.None));
        connect(res2.y, diff2.u2) annotation (Line(
            points={{41,6.10623e-16},{50,6.10623e-16},{50,-20},{-70,-20},{-70,-6},
                {-62,-6}},
            color={200,0,0},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-100},{100,100}}), graphics={Line(
                      points={{62,76},{70,76}},
                      color={200,0,0},
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
                      smooth=Smooth.None)}), Icon(graphics));
      end Imperative_vi;

      model Imperative_iv
        extends Params;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
      public
        FCSys.Connectors.RealInput i annotation (Placement(transformation(
                extent={{-100,20},{-80,40}}), iconTransformation(extent={{-120,
                  -10},{-100,10}})));
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
            points={{-12,-70},{-20,-70},{-20,-30},{-11,-30}},
            color={0,200,0},
            smooth=Smooth.None));
        connect(diff2.y, ind.u) annotation (Line(
            points={{29,-30},{12,-30}},
            color={200,0,0},
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
            color={200,0,0},
            smooth=Smooth.None));
        connect(diff1.y, res1.u) annotation (Line(
            points={{-29,30},{-20,30},{-20,10},{-12,10}},
            color={200,0,0},
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
            color={200,0,0},
            smooth=Smooth.None));
        connect(sum.y, diff2.u1) annotation (Line(
            points={{51,30},{60,30},{60,-24},{52,-24}},
            color={200,0,0},
            smooth=Smooth.None));
        connect(sum.y, v) annotation (Line(
            points={{51,30},{90,30}},
            color={200,0,0},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-100},{100,100}}), graphics={Line(
                      points={{60,58},{68,58}},
                      color={200,0,0},
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
                      smooth=Smooth.None)}), Icon(graphics));
      end Imperative_iv;

      model ImperativeTF_vi
        extends Params;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
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
            color={200,0,0},
            smooth=Smooth.None));
        connect(transferFunction.y, i) annotation (Line(
            points={{11,30},{90,30}},
            color={0,200,0},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-100},{100,100}}), graphics), Icon(graphics));
      end ImperativeTF_vi;

      model ImperativeTF_iv
        extends Params;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
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
            color={200,0,0},
            smooth=Smooth.None));
        connect(transferFunction.u, i) annotation (Line(
            points={{-12,30},{-90,30}},
            color={0,200,0},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-100},{100,100}}), graphics), Icon(graphics));
      end ImperativeTF_iv;

      partial model Params
        parameter SI.Resistance R1=10 "Resistance at temperature T_ref";
        parameter SI.Resistance R2=100 "Resistance at temperature T_ref";
        parameter SI.Inductance L=0.1 "Inductance";
        parameter SI.Capacitance C=0.01 "Capacitance";
      end Params;
    end SwitchCausality;

    package Cascade
      extends Modelica.Icons.Package;
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
        annotation (
          Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                  {100,100}}), graphics),
          experiment(StopTime=2),
          experimentSetupOutput);
      end Test;

      model DeclarativeA
        extends ABParams;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
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
                  {-100,-60},{100,60}}), graphics), Icon(graphics));
      end DeclarativeA;

      model DeclarativeB
        extends ABParams;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
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
                  {-100,-60},{100,60}}), graphics), Icon(graphics));
      end DeclarativeB;

      model DeclarativeAB
        extends ABParams;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
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
                  {-100,-60},{100,60}}), graphics), Icon(graphics));
      end DeclarativeAB;

      model ImperativeA
        extends ABParams;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
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
            color={200,0,0},
            smooth=Smooth.None));

        connect(ind.y, KVL1.u1) annotation (Line(
            points={{-59,6.10623e-16},{-50,6.10623e-16},{-50,6},{-32,6}},
            color={200,0,0},
            smooth=Smooth.None));
        connect(res1.y, ind.u) annotation (Line(
            points={{-99,6.10623e-16},{-94,6.10623e-16},{-94,0},{-90,0},{-90,
                6.66134e-16},{-82,6.66134e-16}},
            color={0,200,0},
            smooth=Smooth.None));
        connect(res1.u, vIn) annotation (Line(
            points={{-122,6.66134e-16},{-134,6.66134e-16},{-134,5.55112e-16},{-150,
                5.55112e-16}},
            color={200,0,0},
            smooth=Smooth.None));
        connect(KVL1.u2, vIn) annotation (Line(
            points={{-32,-6},{-40,-6},{-40,-20},{-130,-20},{-130,5.55112e-16},{
                -150,5.55112e-16}},
            color={200,0,0},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-160,-40},{160,40}}), graphics={Line(
                      points={{-16,30},{-10,30}},
                      color={200,0,0},
                      smooth=Smooth.None),Rectangle(
                      extent={{-15,7},{15,-7}},
                      lineColor={0,0,0},
                      origin={-5,27},
                      rotation=180),Text(
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

      model ImperativeB
        extends ABParams;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
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
            color={200,0,0},
            smooth=Smooth.None));
        connect(res2.y, res3.u) annotation (Line(
            points={{81,6.10623e-16},{86,0},{90,1.27676e-15},{90,6.66134e-16},{
                98,6.66134e-16}},
            color={0,200,0},
            smooth=Smooth.None));
        connect(KVL2.y, res2.u) annotation (Line(
            points={{41,6.10623e-16},{46,0},{50,1.27676e-15},{50,6.66134e-16},{
                58,6.66134e-16}},
            color={200,0,0},
            smooth=Smooth.None));
        connect(KVL2.u1, vIn) annotation (Line(
            points={{18,6},{0,6},{0,5.55112e-16},{-20,5.55112e-16}},
            color={200,0,0},
            smooth=Smooth.None));
        connect(KVL2.u2, res3.y) annotation (Line(
            points={{18,-6},{10,-6},{10,-20},{130,-20},{130,6.10623e-16},{121,
                6.10623e-16}},
            color={200,0,0},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-160,-40},{160,40}}), graphics={Line(
                      points={{124,30},{130,30}},
                      color={200,0,0},
                      smooth=Smooth.None),Rectangle(
                      extent={{-15,7},{15,-7}},
                      lineColor={0,0,0},
                      origin={135,27},
                      rotation=180),Text(
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
        extends ABParams;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
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
            color={200,0,0},
            smooth=Smooth.None));
        connect(KVL1.u1, vIn) annotation (Line(
            points={{-62,46},{-80,46},{-80,40},{-100,40}},
            color={200,0,0},
            smooth=Smooth.None));
        connect(cap.y, KVL1.u2) annotation (Line(
            points={{-61,6.10623e-16},{-70,6.10623e-16},{-70,34},{-62,34}},
            color={200,0,0},
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
            color={200,0,0},
            smooth=Smooth.None));
        connect(KVL2.y, res2.u) annotation (Line(
            points={{-39,-40},{-22,-40}},
            color={200,0,0},
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
            color={200,0,0},
            smooth=Smooth.None));

        connect(KVL2.u2, res3.y) annotation (Line(
            points={{-62,-46},{-70,-46},{-70,-60},{70,-60},{70,-40},{61,-40}},
            color={200,0,0},
            smooth=Smooth.None));

        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-120,-80},{120,80}}), graphics={Line(
                      points={{84,-64},{90,-64}},
                      color={200,0,0},
                      smooth=Smooth.None),Rectangle(
                      extent={{-15,7},{15,-7}},
                      lineColor={0,0,0},
                      origin={95,-67},
                      rotation=180),Text(
                      extent={{90,-62},{110,-66}},
                      lineColor={0,0,0},
                      textString="Voltage"),Text(
                      extent={{90,-68},{110,-72}},
                      lineColor={0,0,0},
                      textString="Current"),Line(
                      points={{84,-70},{90,-70}},
                      color={0,200,0},
                      smooth=Smooth.None)}), Icon(graphics));
      end ImperativeAB;

      model ImperativeABTF
        extends ABParams;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
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
            color={200,0,0},
            smooth=Smooth.None));
        connect(transferFunction.u, vIn) annotation (Line(
            points={{-22,40},{-100,40}},
            color={200,0,0},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-120,-80},{120,80}}), graphics), Icon(graphics));
      end ImperativeABTF;

      model ImperativeABIncorrect
        extends ABParams;
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
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
            color={200,0,0},
            smooth=Smooth.None));
        connect(res1.y, ind.u) annotation (Line(
            points={{-99,6.10623e-16},{-94,0},{-90,1.27676e-15},{-90,
                6.66134e-16},{-82,6.66134e-16}},
            color={0,200,0},
            smooth=Smooth.None));
        connect(vIn, res1.u) annotation (Line(
            points={{-150,5.55112e-16},{-144,0},{-136,1.22125e-15},{-136,
                6.66134e-16},{-122,6.66134e-16}},
            color={200,0,0},
            smooth=Smooth.None));
        connect(vIn, KVL1.u2) annotation (Line(
            points={{-150,5.55112e-16},{-130,5.55112e-16},{-130,-20},{-40,-20},
                {-40,-6},{-32,-6}},
            color={200,0,0},
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
            color={200,0,0},
            smooth=Smooth.None));
        connect(KVL1.y, KVL2.u1) annotation (Line(
            points={{-9,6.10623e-16},{0,6.10623e-16},{0,6},{18,6}},
            color={200,0,0},
            smooth=Smooth.None));
        connect(res3.y, vOut) annotation (Line(
            points={{121,6.10623e-16},{132.5,6.10623e-16},{132.5,5.55112e-16},{
                150,5.55112e-16}},
            color={200,0,0},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-160,-40},{160,40}}), graphics={Line(
                      points={{124,30},{130,30}},
                      color={200,0,0},
                      smooth=Smooth.None),Rectangle(
                      extent={{-15,7},{15,-7}},
                      lineColor={0,0,0},
                      origin={135,27},
                      rotation=180),Text(
                      extent={{130,32},{150,28}},
                      lineColor={0,0,0},
                      textString="Voltage"),Text(
                      extent={{130,26},{150,22}},
                      lineColor={0,0,0},
                      textString="Current"),Line(
                      points={{124,24},{130,24}},
                      color={0,200,0},
                      smooth=Smooth.None)}), Icon(graphics));
      end ImperativeABIncorrect;

      partial model ABParams
        parameter SI.Resistance R1=1 "Resistance at temperature T_ref";
        parameter SI.Resistance R2=1 "Resistance at temperature T_ref";
        parameter SI.Resistance R3=1 "Resistance at temperature T_ref";
        parameter SI.Capacitance C=100e-3 "Capacitance";
      end ABParams;
    end Cascade;
  end DeclarativeVsImperative;

  model DefaultsIcon

    FCSys.BCs.Defaults Defaults
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
    annotation (Diagram(graphics));
  end DefaultsIcon;

  model CellIcon

    FCSys.Assemblies.Cells.Cell Cell
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
  end CellIcon;

  model AnFPIcon "Anode flow plate"

    FCSys.Regions.AnFPs.AnFP AnFP
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
  end AnFPIcon;

  model AnGDLIcon "Anode gas diffusion layer"

    FCSys.Regions.AnGDLs.AnGDL AnGDL
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
  end AnGDLIcon;

  model AnCLIcon "Anode catalyst layer"

    FCSys.Regions.AnCLs.AnCL AnCL
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
  end AnCLIcon;

  model PEMIcon "Proton exchange membrane"

    FCSys.Regions.PEMs.PEM PEM
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
  end PEMIcon;

  model CaCLIcon "Cathode catalyst layer"

    FCSys.Regions.CaCLs.CaCL CaCL
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
  end CaCLIcon;

  model CaGDLIcon "Cathode gas diffusion layer"

    FCSys.Regions.CaGDLs.CaGDL CaGDL
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
  end CaGDLIcon;

  model CaFPIcon "Cathode flow plate"

    FCSys.Regions.CaFPs.CaFP CaFP
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
  end CaFPIcon;

  model SubregionIcon

    FCSys.Subregions.BaseClasses.PartialSubregion Subregion(
      inclX=true,
      inclZ=true,
      inclY=true)
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
  end SubregionIcon;

  model Logo
    extends FCSys.BaseClasses.Icons.Cell;

    annotation (Icon(graphics={Rectangle(
              extent={{-100,100},{100,65}},
              fillPattern=FillPattern.Solid,
              fillColor={255,255,255},
              pattern=LinePattern.None,
              lineColor={0,0,0})}));
  end Logo;

  partial model PhaseIcon

    FCSys.Subregions.Phases.Phase Phase
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
    annotation (structurallyIncomplete=true);
  end PhaseIcon;

  partial model SpeciesIcon

    FCSys.Subregions.Species.Species Species
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
    annotation (structurallyIncomplete=true);
  end SpeciesIcon;

  model PhaseBoundaryIcon

    FCSys.Subregions.PhaseBoundary PhaseBoundary
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
  end PhaseBoundaryIcon;

  model ReactionIcon

    FCSys.Subregions.Reaction Reaction
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
  end ReactionIcon;

  model VolumeIcon

    FCSys.Subregions.Volume Volume
      annotation (Placement(transformation(extent={{-100,-100},{100,100}})));
  end VolumeIcon;

  model Matrix3D

    FCSys.Subregions.BaseClasses.PartialSubregion subreg000(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{-30,-30},{-10,-10}})));
    FCSys.Subregions.BaseClasses.PartialSubregion subreg001(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{-50,-50},{-30,-30}})));
    FCSys.Subregions.BaseClasses.PartialSubregion subreg010(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{-30,20},{-10,40}})));
    FCSys.Subregions.BaseClasses.PartialSubregion subreg011(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{-50,0},{-30,20}})));
    FCSys.Subregions.BaseClasses.PartialSubregion subreg100(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{30,-30},{50,-10}})));
    FCSys.Subregions.BaseClasses.PartialSubregion subreg101(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{10,-50},{30,-30}})));
    FCSys.Subregions.BaseClasses.PartialSubregion subreg110(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{30,20},{50,40}})));
    FCSys.Subregions.BaseClasses.PartialSubregion subreg111(
      inclX=true,
      inclY=true,
      inclZ=true)
      annotation (Placement(transformation(extent={{10,0},{30,20}})));
  equation
    connect(subreg010.negativeY, subreg000.positiveY) annotation (Line(
        points={{-20,20},{-20,-10}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg011.negativeY, subreg001.positiveY) annotation (Line(
        points={{-40,-5.55112e-16},{-40,-30}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg011.negativeZ, subreg010.positiveZ) annotation (Line(
        points={{-35,15},{-25,25}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg111.negativeZ, subreg110.positiveZ) annotation (Line(
        points={{25,15},{35,25}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg110.negativeY, subreg100.positiveY) annotation (Line(
        points={{40,20},{40,-10}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg110.negativeX, subreg010.positiveX) annotation (Line(
        points={{30,30},{-10,30}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg000.positiveX, subreg100.negativeX) annotation (Line(
        points={{-10,-20},{30,-20}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg001.positiveX, subreg101.negativeX) annotation (Line(
        points={{-30,-40},{10,-40}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg001.negativeZ, subreg000.positiveZ) annotation (Line(
        points={{-35,-35},{-25,-25}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg101.negativeZ, subreg100.positiveZ) annotation (Line(
        points={{25,-35},{35,-25}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg111.negativeX, subreg011.positiveX) annotation (Line(
        points={{10,10},{-30,10}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(subreg111.negativeY, subreg101.positiveY) annotation (Line(
        points={{20,-5.55112e-16},{20,-30}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    annotation (Diagram(graphics));
  end Matrix3D;

  model ConnectorHierarchy

    Connectors.BaseClasses.PartialInert Inert annotation (Placement(
          transformation(extent={{-13,-3},{-7,3}}), iconTransformation(extent={
              {-30,-20},{-10,0}})));
    Connectors.ChemicalOutput Chemical annotation (Placement(transformation(
            extent={{20,-10},{40,10}}), iconTransformation(extent={{10,-20},{30,
              0}})));
    Connectors.FaceGeneric Face annotation (Placement(transformation(extent={{0,
              -10},{20,10}}), iconTransformation(extent={{-10,-20},{10,0}})));
    //protected
    Connectors.Material Material annotation (Placement(transformation(extent={{
              20,-30},{40,-10}}), iconTransformation(extent={{10,-40},{30,-20}})));
    Connectors.Entropy Entropy annotation (Placement(transformation(extent={{0,
              -30},{20,-10}}), iconTransformation(extent={{-10,-40},{10,-20}})));
    Connectors.LinearMomentum momentum annotation (Placement(transformation(
            extent={{-20,-30},{0,-10}}), iconTransformation(extent={{-30,-40},{
              -10,-20}})));
    VolumeOrPressure pressure annotation (Placement(transformation(extent={{-40,
              -30},{-20,-10}}), iconTransformation(extent={{-50,-40},{-30,-20}})));
    annotation (Diagram(graphics={Line(
              points={{10,20},{10,0}},
              color={127,127,127},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Line(
              points={{30,20},{30,0}},
              color={127,127,127},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Line(
              points={{-10,0},{-30,-20}},
              color={127,127,127},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Line(
              points={{-10,0},{-10,-20}},
              color={127,127,127},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Line(
              points={{30,0},{30,-20}},
              color={127,127,127},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Line(
              points={{10,0},{10,-20}},
              color={127,127,127},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Line(
              points={{10,0},{-10,-20}},
              color={127,127,127},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Line(
              points={{10,0},{30,-20}},
              color={127,127,127},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Line(
              points={{-10,0},{10,-20}},
              color={127,127,127},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),Text(
              extent={{-40,-8},{-20,-12}},
              lineColor={0,0,0},
              textString="Volume or"),Line(
              points={{56,16},{56,-16}},
              color={0,0,0},
              smooth=Smooth.None),Line(
              points={{54,12},{56,16},{58,12}},
              color={0,0,0},
              smooth=Smooth.None),Text(
              extent={{48,22},{64,18}},
              lineColor={0,0,0},
              textString="High-level"),Text(
              extent={{44,-18},{68,-22}},
              lineColor={0,0,0},
              textString="Fundamental"),Text(
              extent={{22,32},{38,28}},
              lineColor={0,0,0},
              textString="Chemical"),Text(
              extent={{22,28},{38,24}},
              lineColor={0,0,0},
              textString="bus"),Ellipse(
              extent={{27,23},{33,17}},
              lineColor={208,104,0},
              fillColor={255,128,0},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),Text(
              extent={{2,32},{18,28}},
              lineColor={0,0,0},
              textString="Face"),Text(
              extent={{2,28},{18,24}},
              lineColor={0,0,0},
              textString="bus"),Ellipse(
              extent={{7,23},{13,17}},
              lineColor={127,127,127},
              fillColor={191,191,191},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5),Text(
              extent={{-18,8},{-2,4}},
              lineColor={0,0,0},
              textString="Inert"),Text(
              extent={{-18,-8},{-2,-12}},
              lineColor={0,0,0},
              textString="Linear")}), Icon(graphics));
  end ConnectorHierarchy;

  connector VolumeOrPressure "Icon for additivity of volume or pressure"

    annotation (Diagram(graphics={Text(
              extent={{-100,36},{100,76}},
              textString="%name",
              lineColor={0,0,0}),Ellipse(
              extent={{-30,30},{30,-30}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}), Icon(graphics={Ellipse(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
  end VolumeOrPressure;

  model RouterPassThrough

    BCs.Router router
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
  end RouterPassThrough;

  model RouterCrossOver

    BCs.Router router(crossOver=true)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
  end RouterCrossOver;

  package ReactionComparison
    "Comparison of traditional reaction rate to linear alternative"
    extends Modelica.Icons.Package;
    model Test
      extends Modelica.Icons.Example;

      // Nonessential variables for analysis
      output Q.Current err(start=0,fixed=true) "Integrated squared error";

      FCSys.Figures.ReactionComparison.TraditionalReaction traditionalReaction(
          k=0.125*U.A*{(U.cm^3/U.C)^sum(max(traditionalReaction.nu[i], 0) for i
             in 1:size(traditionalReaction.nu, 1)),1.333*(U.cm^3/U.C)^sum(-min(
            traditionalReaction.nu[i], 0) for i in 1:size(traditionalReaction.nu,
            1))})
        annotation (Placement(transformation(extent={{-10,30},{10,50}})));
      Species species1(final V=traditionalReaction.V,final N_IC=
            traditionalReaction.N_IC[1])
        annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));
      Species species2(final V=traditionalReaction.V,final N_IC=
            traditionalReaction.N_IC[2])
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Species species3(final V=traditionalReaction.V,final N_IC=
            traditionalReaction.N_IC[3])
        annotation (Placement(transformation(extent={{30,-10},{50,10}})));
      FCSys.Figures.ReactionComparison.Reaction reaction(final nu=
            traditionalReaction.nu)
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

    equation
      der(err) = (reaction.Xidot - traditionalReaction.Xidot)^2;

      connect(species1.material, reaction.material[1]) annotation (Line(
          points={{-40,6.10623e-16},{-40,-20},{5.55112e-16,-20},{5.55112e-16,-40}},

          color={0,0,0},
          smooth=Smooth.None));

      connect(species2.material, reaction.material[2]) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{6.10623e-16,-20},{0,-40},{
              5.55112e-16,-40}},
          color={0,0,0},
          smooth=Smooth.None));

      connect(species3.material, reaction.material[3]) annotation (Line(
          points={{40,6.10623e-16},{40,-20},{5.55112e-16,-20},{5.55112e-16,-40}},

          color={0,0,0},
          smooth=Smooth.None));

      annotation (
        Diagram(graphics),
        experiment(StopTime=10),
        experimentSetupOutput,
        Commands);
    end Test;

    model TraditionalReaction

      //extends FCSys.BaseClasses.Icons.Names.Top2;

      parameter Integer nu[:]={1,1,-1}
        "<html>stoichiometric coefficients (&nu;)</html>"
        annotation (Evaluate=true);

      parameter Real k[2]=1*U.A*{(U.cm^3/U.C)^sum(max(nu[i], 0) for i in 1:size(
          nu, 1)),(U.cm^3/U.C)^sum(-min(nu[i], 0) for i in 1:size(nu, 1))}
        "Forward and reverse reaction coefficients";
      // Note:  The units are not consistent, but the form is correct.

      parameter Q.Volume V=1*U.cm^3 "Volume";
      parameter Q.ParticleNumber N_IC[n_spec]=fill(4*U.C, n_spec)
        "<html>Initial particle number (<i>N</i><sub>IC</sub></html>";

      Q.Current Xidot(nominal=1*U.A) "Reaction rate";
      Q.ParticleNumber N[n_spec](
        each nominal=1*U.C,
        final start=N_IC,
        each fixed=true) "Particle number";

      // Alias variable
      Q.ParticleNumberVolumic rho[n_spec](each nominal=1*U.C/U.cm^3)
        "Volumic particle number";

    protected
      final parameter Integer n_spec=size(nu, 1) "Number of species";

    equation
      // Alias
      V*rho = N;

      // Reaction rate
      Xidot = k[1]*product(if nu[i] > 0 then rho[i]^nu[i] else 1 for i in 1:
        n_spec) - k[2]*product(if nu[i] < 0 then rho[i]^(-nu[i]) else 1 for i
         in 1:n_spec) "Net reaction rate";

      // Material conservation
      der(N)/U.s = -nu*Xidot "stoichiometry";

      annotation (Icon(graphics={Rectangle(
                  extent={{-100,100},{100,-100}},
                  lineColor={0,0,0},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),Text(
                  extent={{-100,40},{100,80}},
                  textString="%name",
                  lineColor={0,0,0})}));
    end TraditionalReaction;

    model Reaction
      "Model for a reaction characterized by equilibrium of volumetric density"
      //extends FCSys.BaseClasses.Icons.Names.Top2;

      parameter Integer nu[:]={1,1,-1}
        "<html>stoichiometric coefficients (&nu;)</html>";

      Q.Current Xidot(nominal=1*U.A) "Reaction rate";

      Material material[n_spec] "Connection for material exchange"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

    protected
      final parameter Integer n_spec=size(nu, 1) "Number of species";

    equation
      nu*material.rho = 0 - 1 "Equilibrium";
      material.Ndot = nu*Xidot "stoichiometry";

      annotation (Diagram(graphics), Icon(graphics={Rectangle(
                  extent={{-100,100},{100,-100}},
                  lineColor={0,0,0},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),Text(
                  extent={{-100,40},{100,80}},
                  textString="%name",
                  lineColor={0,0,0})}));
    end Reaction;

    model Species
      "Model for a single species with linear density-driven reaction rate"

      // extends FCSys.BaseClasses.Icons.Names.Top2;

      parameter Real R=1*U.s/U.cm^3 "Material resistance";
      // Note:  The dimensions are not consistent, but the form is correct.
      parameter Q.Volume V=1*U.cm^3 "Volume";
      parameter Q.ParticleNumber N_IC=4*U.C
        "<html>Initial particle number (<i>N</i><sub>IC</sub></html>";

      Q.ParticleNumber N(
        nominal=1*U.C,
        start=N_IC,
        fixed=true) "Particle number";

      // Alias variable
      Q.ParticleNumberVolumic rho(nominal=1*U.C/U.cm^3)
        "Volumic particle number";

      Material material "Connections for material exchange"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

    equation
      // Alias
      V*rho = N;

      // Reaction rate
      R*material.Ndot = material.rho - rho;

      // Material conservation
      der(N)/U.s = material.Ndot;
      annotation (Icon(graphics={Rectangle(
                  extent={{-100,100},{100,-100}},
                  lineColor={0,0,0},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),Text(
                  extent={{-100,40},{100,80}},
                  textString="%name",
                  lineColor={0,0,0})}));
    end Species;

    connector Material "Connector for density-driven reaction"
      Q.ParticleNumberVolumic rho(nominal=1*U.C/U.cm^3) "Volumetric density";
      flow Q.Current Ndot(nominal=1*U.A) "Current";
    end Material;
  end ReactionComparison;

end Figures;
