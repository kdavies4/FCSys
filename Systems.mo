within FCSys;
package Systems
  "Models of systems containing both hardware devices and control blocks"

  extends Modelica.Icons.Package;
  extends FCSys.BaseClasses.Icons.PackageUnderConstruction;
  package Examples "Examples and tests"
    extends Modelica.Icons.ExamplesPackage;
    model FCPlantNoRecirc
      extends Modelica.Icons.Example;

      package Medium_an = FCSys.WorkInProgress.MolarBasis.H2O_H2_CO_Molar
        "Anode medium model";
      package Medium_ca = FCSys.WorkInProgress.MolarBasis.H2OAndO2_N2_Molar
        "Cathode medium model";

      replaceable Modelica.Blocks.Sources.Constant valveSetpoint(k=0.5)
        constrainedby Modelica.Blocks.Interfaces.SO annotation (
          choicesAllMatching=true, Placement(transformation(
            origin={-90,-60},
            extent={{10,-10},{-10,10}},
            rotation=180)));
      FCSys.Systems.FC.FCNoRecirc fCPlant annotation (Placement(transformation(
              extent={{0,-20},{40,20}}, rotation=0)));
      Modelica.Fluid.Sources.FixedBoundary_pTX ambientCa(
        redeclare replaceable package Medium = Medium_ca,
        T=ambient.default_T_ambient,
        flowDirection=Modelica.Fluid.Types.SourceFlowDirection.InToPort,
        p=ambient.p_default_ambient,
        X=Medium_ca.X_default) annotation (Placement(transformation(extent={{
                100,20},{80,40}}, rotation=0)));
      Modelica.Fluid.Sources.FixedBoundary_pTX ambientAn(
        redeclare replaceable package Medium = Medium_an,
        T=ambient.default_T_ambient,
        flowDirection=Modelica.Fluid.Types.SourceFlowDirection.InToPort,
        p=ambient.p_default_ambient,
        X=Medium_an.X_default) annotation (Placement(transformation(extent={{
                100,-40},{80,-20}}, rotation=0)));
      replaceable Modelica.Blocks.Sources.Constant fanSpeed(k=0.5)
        constrainedby Modelica.Blocks.Interfaces.SO annotation (
          choicesAllMatching=true, Placement(transformation(
            origin={-90,-30},
            extent={{10,-10},{-10,10}},
            rotation=180)));
      inner Modelica.Fluid.Ambient ambient(p_default_ambient=
            SI.Conversions.from_bar(1.01325), default_T_ambient=
            SI.Conversions.from_degC(20)) annotation (Placement(transformation(
              extent={{80,80},{100,100}}, rotation=0)));
      Modelica.Thermal.HeatTransport.Sources.FixedTemperature fixedTemperature(
          T=ambient.default_T_ambient) annotation (Placement(transformation(
              extent={{-100,60},{-80,80}}, rotation=0)));
      Modelica.Electrical.Analog.Basic.Resistor resistor(R=1) annotation (
          Placement(transformation(
            origin={-20,-80},
            extent={{-10,-10},{10,10}},
            rotation=270)));
      Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(
            transformation(extent={{0,-110},{20,-90}}, rotation=0)));
      Modelica.Fluid.Pipes.BaseClasses.PortVolume H2tank(
        redeclare package Medium =
            FCSys.WorkInProgress.MolarBasis.H2O_H2_CO_Molar,
        V=1e3,
        p_start=10e5,
        use_T_start=true,
        T_start=ambient.default_T_ambient,
        X_start={0,1,0},
        initType=Modelica.Fluid.Types.Init.NoInit) annotation (Placement(
            transformation(extent={{-70,-20},{-50,0}}, rotation=0)));

      Modelica.Thermal.HeatTransport.Species.ThermalConductor thermalConductor2
        (G=100) annotation (Placement(transformation(
            origin={-60,30},
            extent={{-10,-10},{10,10}},
            rotation=270)));
      Modelica.Fluid.Pipes.BaseClasses.PortVolume O2tank(
        V=1e3,
        p_start=10e5,
        use_T_start=true,
        T_start=ambient.default_T_ambient,
        X_start={0,1,0},
        redeclare package Medium =
            FCSys.WorkInProgress.MolarBasis.H2OAndO2_N2_Molar,
        initType=Modelica.Fluid.Types.Init.NoInit) annotation (Placement(
            transformation(extent={{-40,0},{-20,20}}, rotation=0)));

      Modelica.Thermal.HeatTransport.Species.ThermalConductor thermalConductor1
        (G=100) annotation (Placement(transformation(
            origin={-30,50},
            extent={{-10,-10},{10,10}},
            rotation=270)));
    equation
      connect(valveSetpoint.y, fCPlant.exitValvePosAn)
        annotation (Line(points={{-79,-60},{42,-60},{42,2}}, color={0,0,127}));
      connect(fCPlant.exitValvePosCa, valveSetpoint.y) annotation (Line(points=
              {{42,18},{31.7,18},{31.7,-60},{-79,-60}}, color={0,0,127}));
      connect(fanSpeed.y, fCPlant.coolingFanSpeed) annotation (Line(points={{-79,
              -30},{10,-30},{10,-20.8}}, color={0,0,127}));
      connect(fCPlant.fluidPort_ca2, ambientCa.port) annotation (Line(points={{
              40,10},{60,10},{60,30},{80,30}}, color={0,127,255}));
      connect(fCPlant.fluidPort_an2, ambientAn.port) annotation (Line(points={{
              40,-10},{60,-10},{60,-30},{80,-30}}, color={0,127,255}));
      connect(fixedTemperature.port, fCPlant.port_a)
        annotation (Line(points={{-80,70},{20,70},{20,20}}, color={191,0,0}));
      connect(resistor.n, fCPlant.n) annotation (Line(points={{-20,-90},{48,-90},
              {48,0},{40,0}},color={0,0,255}));
      connect(resistor.p, fCPlant.p)
        annotation (Line(points={{-20,-70},{-20,0},{0,0}}, color={0,0,255}));
      connect(ground.p, resistor.n)
        annotation (Line(points={{10,-90},{-20,-90}}, color={0,0,255}));
      connect(H2tank.port, fCPlant.fluidPort_an1)
        annotation (Line(points={{-60,-10},{0,-10}}, color={0,127,255}));
      connect(thermalConductor2.port_b, H2tank.thermalPort)
        annotation (Line(points={{-60,20},{-60,0}}, color={191,0,0}));
      connect(thermalConductor2.port_a, fixedTemperature.port) annotation (Line(
            points={{-60,40},{-60,70},{-80,70}}, color={191,0,0}));
      connect(thermalConductor1.port_a, fixedTemperature.port) annotation (Line(
            points={{-30,60},{-30,70},{-80,70}}, color={191,0,0}));
      connect(thermalConductor1.port_b, O2tank.thermalPort) annotation (Line(
            points={{-30,40},{-30,30},{-30,30},{-30,20}}, color={191,0,0}));
      connect(O2tank.port, fCPlant.fluidPort_ca1)
        annotation (Line(points={{-30,10},{0,10}}, color={0,127,255}));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}), graphics));
    end FCPlantNoRecirc;
  end Examples;

  package FC
    extends Modelica.Icons.Package;
    model FCNoRecirc
      "Fuel cell system with fuel cell stack, H2 and H2O tanks, DC/DC converter, valves, pumps, heat exchangers, actuators, and sensors"

      extends FCSys.WorkInProgress.FCSysPlant;
      SI.Voltage v "Voltage drop between the two pins (= pinP.v - pinP.v)";

      FCSys.Connectors.FaceBus wireP
        "Positive pin Positive pin (potential pinP.v > pinP.v for positive voltage drop v)"
        annotation (Placement(transformation(extent={{170,110},{190,130}},
              rotation=0), iconTransformation(extent={{90,-50},{110,-70}})));
      FCSys.Connectors.FaceBus wireN "Negative pin" annotation (Placement(
            transformation(
            extent={{-10,10},{10,-10}},
            rotation=180,
            origin={180,-120}), iconTransformation(
            extent={{-10,10},{10,-10}},
            rotation=180,
            origin={-100,-60})));

      FCSys.WorkInProgress.AssembliesStacksStack stack(
        T_start=ambient.default_T_ambient,
        p_start=ambient.p_default_ambient,
        redeclare replaceable package Medium_an = Medium_an,
        redeclare replaceable package Medium_ca = Medium_ca) "Fuel cell model"
        annotation (Placement(transformation(
            origin={20,-30},
            extent={{10,10},{-10,-10}},
            rotation=90)));

      // TODO: Use inner/outer expandable connectors for control buses.
      FCSys.Systems.Humidifier.Humidifier anHumidifier annotation (Placement(
            transformation(extent={{-50,-14},{-30,6}}, rotation=0)));
      FCSys.Assemblies.ClosedTank tankH2 annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-190,-10})));
      FCSys.Systems.Pump.Pump coolingPump annotation (Placement(transformation(
            extent={{10,10},{-10,-10}},
            rotation=180,
            origin={-100,62})));
      FCSys.Systems.DCDC.DCDC dCDC annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=270,
            origin={174,-12})));
      FCSys.Assemblies.OpenTank H2OTank annotation (Placement(transformation(
              extent={{-30,10},{-10,30}}, rotation=0)));
      FCSys.Systems.Pump.Pump compressor annotation (Placement(transformation(
              extent={{-110,-60},{-90,-40}}, rotation=0)));
      FCSys.Systems.Pump.Pump condensatePump annotation (Placement(
            transformation(extent={{-110,22},{-90,42}}, rotation=0)));
      FCSys.Connectors.FaceBus ambientN annotation (Placement(transformation(
              extent={{-210,70},{-190,90}}), iconTransformation(extent={{-110,
                50},{-90,70}})));
      outer FCSys.Systems.FC.Interfaces.ActBusIn actBusIn annotation (Placement(
            transformation(extent={{-210,-90},{-190,-70}}), iconTransformation(
              extent={{-92,-10},{-72,10}})));
      outer FCSys.Systems.FC.Interfaces.SenBusOut senBusOut annotation (
          Placement(transformation(extent={{190,-90},{210,-70}}),
            iconTransformation(extent={{72,-10},{92,10}})));
      FCSys.Systems.Valve.Valve caExitValve annotation (Placement(
            transformation(extent={{110,-62},{130,-42}}, rotation=0)));
      FCSys.Systems.Valve.Valve anExitValve annotation (Placement(
            transformation(extent={{110,-20},{130,0}}, rotation=0)));
      FCSys.Systems.Valve.Valve anInletValve annotation (Placement(
            transformation(extent={{-180,-22},{-160,-2}}, rotation=0)));
      FCSys.Systems.FluidHeater.FluidHeater anPreheater
        annotation (Placement(transformation(extent={{-110,-24},{-90,-4}})));
      FCSys.Assemblies.Condenser condenser
        annotation (Placement(transformation(extent={{70,-54},{90,-34}})));
      FCSys.Systems.Humidifier.Humidifier caHumidifier
        annotation (Placement(transformation(extent={{-50,-54},{-30,-34}})));
      FCSys.WorkInProgress.AssembliesHeatExchanger radiator
        annotation (Placement(transformation(extent={{-60,66},{-40,86}})));
      FCSys.Connectors.FaceBus ambientP annotation (Placement(transformation(
              extent={{190,70},{210,90}}), iconTransformation(extent={{-110,50},
                {-90,70}})));
    equation
      v = pinP.v - pinP.v;

      connect(condensatePump.actBusIn, actBusIn);
      connect(actBusIn, compressor.actBusIn);
      connect(compressor.senBusOut, senBusOut);
      connect(condensatePump.senBusOut, senBusOut);
      connect(anInletValve.actBusIn, actBusIn);
      connect(anPreheater.actBusIn, actBusIn);
      connect(anInletValve.senBusOut, senBusOut);

      connect(anPreheater.senBusOut, senBusOut);
      connect(caExitValve.senBusOut, senBusOut);
      connect(anExitValve.senBusOut, senBusOut);
      connect(anExitValve.actBusIn, actBusIn);
      connect(caExitValve.actBusIn, actBusIn);
      connect(caHumidifier.actBusIn, actBusIn);
      connect(caHumidifier.actBusOut, senBusOut);
      connect(anHumidifier.actBusIn, actBusIn);
      connect(anHumidifier.actBusOut, senBusOut);
      connect(dCDC.actBusIn, actBusIn);

      connect(senBusOut, dCDC.senBusOut);
      connect(coolingPump.actBusIn, actBusIn);
      connect(coolingPump.senBusOut, senBusOut);

      connect(wireN, dCDC.matN2) annotation (Line(
          points={{180,-120},{180,-22}},
          color={0,200,0},
          smooth=Smooth.None));
      connect(dCDC.matP2, wireP) annotation (Line(
          points={{180,-2},{180,120}},
          color={0,200,0},
          smooth=Smooth.None));
      connect(radiator.pipeP2, ambientP) annotation (Line(
          points={{-40,80},{200,80}},
          color={0,128,255},
          smooth=Smooth.None));
      connect(radiator.pipe_N2, ambientN) annotation (Line(
          points={{-60,80},{-200,80}},
          color={0,128,255},
          smooth=Smooth.None));
      connect(compressor.pipeN, ambientN) annotation (Line(
          points={{-110,-50},{-150,-50},{-150,80},{-200,80}},
          color={0,128,255},
          smooth=Smooth.None));
      connect(caHumidifier.mixturePipeP, compressor.pipeP) annotation (Line(
          points={{-50,-50},{-90,-50}},
          color={0,128,255},
          smooth=Smooth.None));
      connect(coolingPump.pipeP, radiator.pipe_N1) annotation (Line(
          points={{-90,62},{-70,62},{-70,72},{-60,72}},
          color={192,0,192},
          smooth=Smooth.None));
      connect(condenser.mixturePipeP, caExitValve.pipeN) annotation (Line(
          points={{90,-50},{110,-50}},
          color={0,128,255},
          smooth=Smooth.None));
      connect(anExitValve.pipeP, ambientP) annotation (Line(
          points={{130,-8},{150,-8},{150,80},{200,80}},
          color={255,128,0},
          smooth=Smooth.None));

      connect(condensatePump.pipeP, H2OTank.matP) annotation (Line(
          points={{-90,32},{-20,32},{-20,24}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(H2OTank.matN, anHumidifier.H2OPipe) annotation (Line(
          points={{-20,16},{-20,2},{-40,2}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(caHumidifier.H2OPipe, H2OTank.matN) annotation (Line(
          points={{-40,-38},{-20,-38},{-20,16}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(condensatePump.wireP, coolingPump.wireP) annotation (Line(
          points={{-90,28},{-80,28},{-80,58},{-90,58}},
          color={0,200,0},
          smooth=Smooth.None));
      connect(compressor.wireP, condensatePump.wireP) annotation (Line(
          points={{-90,-54},{-80,-54},{-80,28},{-90,28}},
          color={0,200,0},
          smooth=Smooth.None));
      connect(condensatePump.wireN, compressor.wireN) annotation (Line(
          points={{-110,28},{-120,28},{-120,-54},{-110,-54}},
          color={0,200,0},
          smooth=Smooth.None));
      connect(coolingPump.wireN, compressor.wireN) annotation (Line(
          points={{-110,58},{-120,58},{-120,-54},{-110,-54}},
          color={0,200,0},
          smooth=Smooth.None));
      connect(anPreheater.wireP, compressor.wireP) annotation (Line(
          points={{-90,-18},{-80,-18},{-80,-54},{-90,-54}},
          color={0,200,0},
          smooth=Smooth.None));
      connect(anPreheater.wireN, compressor.wireN) annotation (Line(
          points={{-110,-18},{-120,-18},{-120,-54},{-110,-54}},
          color={0,200,0},
          smooth=Smooth.None));
      connect(tankH2.normal, anInletValve.pipeN) annotation (Line(
          points={{-190,-10},{-180,-10}},
          color={255,128,0},
          smooth=Smooth.None));
      connect(anHumidifier.mixturePipeN, stack.anPipeN) annotation (Line(
          points={{-30,-10},{0,-10},{0,-24},{10,-24}},
          color={255,128,0},
          smooth=Smooth.None));
      connect(stack.anPipeP, anExitValve.pipeN) annotation (Line(
          points={{30,-24},{40,-24},{40,-8},{110,-8}},
          color={255,128,0},
          smooth=Smooth.None));
      connect(stack.caPipeP, condenser.mixturePipeN) annotation (Line(
          points={{30,-36},{40,-36},{40,-50},{70,-50}},
          color={0,128,255},
          smooth=Smooth.None));
      connect(radiator.pipeP1, condenser.coolantPipeP) annotation (Line(
          points={{-40,72},{100,72},{100,-40},{90,-40}},
          color={192,0,192},
          smooth=Smooth.None));
      connect(condenser.coolantPipeN, stack.coolantPipeP) annotation (Line(
          points={{70,-40},{50,-40},{50,-30},{30,-30}},
          color={192,0,192},
          smooth=Smooth.None));
      connect(stack.coolantPipeN, coolingPump.pipeN) annotation (Line(
          points={{10,-30},{-10,-30},{-10,-70},{-140,-70},{-140,62},{-110,62}},

          color={192,0,192},
          smooth=Smooth.None));

      connect(caExitValve.pipeP, ambientP) annotation (Line(
          points={{130,-50},{150,-50},{150,80},{200,80}},
          color={0,128,255},
          smooth=Smooth.None));
      connect(stack.anWire, dCDC.wireP1) annotation (Line(
          points={{18,-20},{18,100},{168,100},{168,-2}},
          color={0,200,0},
          smooth=Smooth.None));
      connect(coolingPump.wireN, stack.anWire) annotation (Line(
          points={{-110,58},{-120,58},{-120,100},{18,100},{18,-20}},
          color={0,200,0},
          smooth=Smooth.None));
      connect(dCDC.wireN1, stack.caWire) annotation (Line(
          points={{168,-22},{168,-100},{18,-100},{18,-40}},
          color={0,200,0},
          smooth=Smooth.None));
      connect(condensatePump.pipeN, condenser.H2OPipe) annotation (Line(
          points={{-110,32},{-130,32},{-130,-84},{80,-84},{80,-54}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(compressor.wireN, stack.caWire) annotation (Line(
          points={{-110,-54},{-120,-54},{-120,-100},{18,-100},{18,-40}},
          color={0,200,0},
          smooth=Smooth.None));
      connect(anPreheater.pipeN, anInletValve.pipeP) annotation (Line(
          points={{-110,-10},{-160,-10}},
          color={255,128,0},
          smooth=Smooth.None));
      connect(caHumidifier.mixturePipeN, stack.caPipeN) annotation (Line(
          points={{-30,-50},{0,-50},{0,-36},{10,-36}},
          color={0,128,255},
          smooth=Smooth.None));
      connect(anHumidifier.mixturePipeP, anPreheater.pipeP) annotation (Line(
          points={{-50,-10},{-90,-10}},
          color={255,128,0},
          smooth=Smooth.None));
      annotation (
        Diagram(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-200,-120},{200,120}},
            initialScale=0.1), graphics),
        experiment(
          StopTime=2000,
          Tolerance=1e-005,
          Algorithm="Dassl"),
        experimentSetupOutput,
        Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                100,100}}), graphics));
    end FCNoRecirc;

    package Interfaces
      extends Modelica.Icons.InterfacesPackage;

      expandable connector ActBusIn
        "Data bus for inputs to a plant's actuators"
        extends FCSys.BaseClasses.Icons.SignalBuses.In;
        input SI.Voltage v_ref;
        input Real valvePos_an_in;
        input Real valvePos_an_out;
        input Real valvePos_ca_out;
        input Real valvePos_humid_an;
        input Real valvePos_humid_ca;
        input Real qdot_compressor;
        input Real qdot_H2Opump;
        input Real qdot_coolingPump;

        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics));
      end ActBusIn;

      expandable connector ActBusOut
        "Data bus for inputs to a plant's actuators"
        extends FCSys.BaseClasses.Icons.SignalBuses.Out;
        output SI.Current current_ref;
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics));
      end ActBusOut;

      expandable connector SenBusIn
        "Data bus for outputs from a plant's sensors"
        extends FCSys.BaseClasses.Icons.SignalBuses.In;
        input Q.Voltage current;
        input Q.Rotation angle;
        input Q.RotationalVelocity speed;
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics));
      end SenBusIn;

      expandable connector SenBusOut
        "Data bus for outputs from a plant's sensors"
        extends FCSys.BaseClasses.Icons.SignalBuses.Out;
        output Q.Voltage current;
        output Q.Rotation angle;
        output Q.RotationalVelocity speed;
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics));
      end SenBusOut;

      model ActArrayToBus

        FCSys.Systems.FC.Interfaces.ActBusOut actBusOut annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={100,0}),iconTransformation(
              extent={{-80,-20},{-40,20}},
              rotation=0,
              origin={100,0})));
        FCSys.Connectors.RealInput actIN[1] annotation (Placement(
              transformation(extent={{-110,-10},{-90,10}}), iconTransformation(
                extent={{-60,-20},{-20,20}})));
      equation

        connect(actIN[1], actBusOut.current_ref) annotation (Line(
            points={{-100,0},{-4,0},{-4,5.55112e-16},{100,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{6,3},{6,3}}));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-100},{100,100}}), graphics), Icon(coordinateSystem(
                preserveAspectRatio=true, extent={{-60,-40},{60,40}}), graphics
              ={Rectangle(
                      extent={{-40,40},{40,-40}},
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Solid,
                      pattern=LinePattern.Dash,
                      lineColor={0,0,0}),Line(
                      points={{-28,4},{28,4}},
                      color={0,0,127},
                      smooth=Smooth.None),Line(
                      points={{-28,-4},{34,-4}},
                      color={0,0,127},
                      smooth=Smooth.None)}));
      end ActArrayToBus;

      model SenBusToArray

        SenBusIn senBusIn annotation (Placement(transformation(extent={{-110,-10},
                  {-90,10}}), iconTransformation(extent={{-60,-20},{-20,20}})));
        FCSys.Connectors.RealOutput senOut[3] annotation (Placement(
              transformation(extent={{90,-10},{110,10}}), iconTransformation(
                extent={{20,-20},{60,20}})));
      equation

        connect(senBusIn.angle, senOut[3]) annotation (Line(
            points={{-100,5.55112e-16},{0,5.55112e-16},{0,6.66667},{100,6.66667}},

            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        connect(senBusIn.speed, senOut[2]) annotation (Line(
            points={{-100,5.55112e-16},{0,5.55112e-16},{0,0},{100,0}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));
        connect(senBusIn.current, senOut[1]) annotation (Line(
            points={{-100,5.55112e-16},{0,5.55112e-16},{0,-6.66667},{100,-6.66667}},

            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-100},{100,100}}), graphics), Icon(coordinateSystem(
                preserveAspectRatio=true, extent={{-60,-40},{60,40}}), graphics
              ={Rectangle(
                      extent={{-40,40},{40,-40}},
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Solid,
                      pattern=LinePattern.Dash,
                      lineColor={0,0,0}),Line(
                      points={{-34,-4},{20,-4}},
                      color={0,0,127},
                      smooth=Smooth.None),Line(
                      points={{-28,4},{20,4}},
                      color={0,0,127},
                      smooth=Smooth.None)}));
      end SenBusToArray;
    end Interfaces;
  end FC;

  package Humidifier
    extends Modelica.Icons.Package;

    model Humidifier
      // TODO: Complete this.

      extends FCSys.BaseClasses.Icons.Name.Top4;

      FCSys.Connectors.FaceBus mixturePipeN annotation (Placement(
            transformation(extent={{90,-70},{110,-50}}), iconTransformation(
              extent={{90,-70},{110,-50}})));
      FCSys.Systems.Humidifier.Interfaces.ActBusIn actBusIn annotation (
          Placement(transformation(extent={{-32,60},{-12,80}}),
            iconTransformation(extent={{-32,10},{-12,30}})));
      FCSys.Systems.Humidifier.Interfaces.ActBusOut actBusOut annotation (
          Placement(transformation(extent={{14,60},{34,80}}),
            iconTransformation(extent={{14,10},{34,30}})));
      FCSys.Systems.Valve.Valve valve annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,50})));
      FCSys.Connectors.FaceBus H2OPipe annotation (Placement(transformation(
              extent={{-10,50},{10,70}}), iconTransformation(extent={{-10,50},{
                10,70}})));
      FCSys.Connectors.FaceBus mixturePipeP annotation (Placement(
            transformation(extent={{90,-10},{110,10}}), iconTransformation(
              extent={{-110,-70},{-90,-50}})));
      annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}), graphics={Line(
                  points={{10,-30},{0,-50}},
                  color={0,0,255},
                  thickness=0.5),Line(
                  points={{-10,-30},{0,-50}},
                  color={0,0,255},
                  thickness=0.5),Line(
                  points={{0,50},{0,-50}},
                  color={0,0,255},
                  thickness=0.5),Rectangle(
                  extent={{-10,0},{10,-20}},
                  lineColor={0,0,0},
                  lineThickness=0.5,
                  fillPattern=FillPattern.VerticalCylinder,
                  fillColor={135,135,135}),Rectangle(
                  extent={{-40,40},{40,0}},
                  lineThickness=0.5,
                  fillPattern=FillPattern.VerticalCylinder,
                  fillColor={0,0,0},
                  pattern=LinePattern.None),Rectangle(
                  extent={{-10,100},{10,-100}},
                  lineColor={0,0,0},
                  lineThickness=0.5,
                  fillPattern=FillPattern.VerticalCylinder,
                  fillColor={135,135,135},
                  origin={0,-60},
                  rotation=90)}), Diagram(coordinateSystem(preserveAspectRatio=
                true, extent={{-100,-100},{100,100}}), graphics));
    end Humidifier;

    package Interfaces
      extends Modelica.Icons.InterfacesPackage;

      expandable connector ActBusIn
        "Data bus for inputs to a plant's actuators"
        extends FCSys.BaseClasses.Icons.SignalBuses.In;
        input SI.Voltage v_ref;
        input Real valvePos_an_in;
        input Real valvePos_an_out;
        input Real valvePos_ca_out;
        input Real valvePos_humid_an;
        input Real valvePos_humid_ca;
        input Real qdot_compressor;
        input Real qdot_H2Opump;
        input Real qdot_coolingPump;

        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics));
      end ActBusIn;

      expandable connector ActBusOut
        "Data bus for inputs to a plant's actuators"
        extends FCSys.BaseClasses.Icons.SignalBuses.Out;
        output SI.Current current_ref;
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics));
      end ActBusOut;

      expandable connector SenBusIn
        "Data bus for outputs from a plant's sensors"
        extends FCSys.BaseClasses.Icons.SignalBuses.In;
        input Q.Voltage current;
        input Q.Rotation angle;
        input Q.RotationalVelocity speed;
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics));
      end SenBusIn;

      expandable connector SenBusOut
        "Data bus for outputs from a plant's sensors"
        extends FCSys.BaseClasses.Icons.SignalBuses.Out;
        output Q.Voltage current;
        output Q.Rotation angle;
        output Q.RotationalVelocity speed;
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics));
      end SenBusOut;

      model ActArrayToBus

        FCSys.Systems.Humidifier.Interfaces.ActBusOut actBusOut annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={100,0}),iconTransformation(
              extent={{-80,-20},{-40,20}},
              rotation=0,
              origin={100,0})));
        FCSys.Connectors.RealInput actIN[1] annotation (Placement(
              transformation(extent={{-110,-10},{-90,10}}), iconTransformation(
                extent={{-60,-20},{-20,20}})));
      equation

        connect(actIN[1], actBusOut.current_ref) annotation (Line(
            points={{-100,0},{-4,0},{-4,5.55112e-16},{100,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{6,3},{6,3}}));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-100},{100,100}}), graphics), Icon(coordinateSystem(
                preserveAspectRatio=true, extent={{-60,-40},{60,40}}), graphics
              ={Rectangle(
                      extent={{-40,40},{40,-40}},
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Solid,
                      pattern=LinePattern.Dash,
                      lineColor={0,0,0}),Line(
                      points={{-28,4},{28,4}},
                      color={0,0,127},
                      smooth=Smooth.None),Line(
                      points={{-28,-4},{34,-4}},
                      color={0,0,127},
                      smooth=Smooth.None)}));
      end ActArrayToBus;

      model SenBusToArray

        SenBusIn senBusIn annotation (Placement(transformation(extent={{-110,-10},
                  {-90,10}}), iconTransformation(extent={{-60,-20},{-20,20}})));
        FCSys.Connectors.RealOutput senOut[3] annotation (Placement(
              transformation(extent={{90,-10},{110,10}}), iconTransformation(
                extent={{20,-20},{60,20}})));
      equation

        connect(senBusIn.angle, senOut[3]) annotation (Line(
            points={{-100,5.55112e-16},{0,5.55112e-16},{0,6.66667},{100,6.66667}},

            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        connect(senBusIn.speed, senOut[2]) annotation (Line(
            points={{-100,5.55112e-16},{0,5.55112e-16},{0,0},{100,0}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));
        connect(senBusIn.current, senOut[1]) annotation (Line(
            points={{-100,5.55112e-16},{0,5.55112e-16},{0,-6.66667},{100,-6.66667}},

            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-100},{100,100}}), graphics), Icon(coordinateSystem(
                preserveAspectRatio=true, extent={{-60,-40},{60,40}}), graphics
              ={Rectangle(
                      extent={{-40,40},{40,-40}},
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Solid,
                      pattern=LinePattern.Dash,
                      lineColor={0,0,0}),Line(
                      points={{-34,-4},{20,-4}},
                      color={0,0,127},
                      smooth=Smooth.None),Line(
                      points={{-28,4},{20,4}},
                      color={0,0,127},
                      smooth=Smooth.None)}));
      end SenBusToArray;
    end Interfaces;
  end Humidifier;

  package DCDC
    extends Modelica.Icons.Package;

    model DCDC

      // TODO: Complete this.

      extends FCSys.WorkInProgress.BaseClassesIconsTransformer;
      import DomainQ = FCSys.Domains.MagneticElectroMechanical.Quantities;

      parameter FCSys.Systems.BaseClasses.LHSchoice LHS=BaseClasses.LHSchoice.Delta_v1
        "Choice of variable for the left hand side (LHS) of the governing equation"
        annotation (Evaluate=true);
      parameter FCSys.Systems.BaseClasses.RHSchoice RHS=BaseClasses.RHSchoice.Delta_v2
        "Choice of variable for the right hand side (RHS) of the governing equation"
        annotation (Evaluate=true);
      parameter Boolean specAsParam=true
        "true: use a parameter to specify the governing equation; false: use dynamic input";

      parameter Q.Unity k_ratio=1 "Voltage or current ratio of side 1:side 2"
        annotation (Dialog(enable=((LHS == RHS) and specAsParam)));
      parameter FCSys.WorkInProgress.Magnetic.Conductance k_R=1*U.ohm
        "Hall resistance across side 1 (Delta_v1/Delta_qdot2)"
        annotation (Dialog(enable=((LHS == 1) and (RHS == 2) and specAsParam)));
      parameter FCSys.WorkInProgress.Magnetic2.Conductance k_sigma=1*U.S
        "Hall conductance through side 1 (Delta_qdot1/Delta_v2)"
        annotation (Dialog(enable=((LHS == 2) and (RHS == 1) and specAsParam)));
      parameter FCSys.WorkInProgress.Magnetic2.Effort k_Delta_v1=1*U.V
        "Voltage across side 1"
        annotation (Dialog(enable=((LHS == 1) and (RHS == 3) and specAsParam)));
      parameter FCSys.WorkInProgress.Magnetic2.Flow k_Delta_qdot1=1*U.A
        "Total current through side 1"
        annotation (Dialog(enable=((LHS == 2) and (RHS == 3) and specAsParam)));

      FCSys.WorkInProgress.Magnetic2.Effort Delta_v1
        "Electrical driving force (for net current positive from chargeP to chargeN)";
      FCSys.WorkInProgress.Magnetic2.Flow Delta_qdot1
        "Total current of both pins (positive by right hand rule from chargeP to chargeN)";
      FCSys.WorkInProgress.Magnetic2.Effort Sigma_v1
        "Electrical sink force (for net current into each flange)";
      FCSys.WorkInProgress.Magnetic2.Flow Sigma_qdot1
        "Net current stored or destroyed (positive into each pin)";

      FCSys.WorkInProgress.Magnetic2.Effort Delta_v2
        "Electrical driving force (for net current positive from chargeP to chargeN)";
      FCSys.WorkInProgress.Magnetic2.Flow Delta_qdot2
        "Total current of both pins (positive by right hand rule from chargeP to chargeN)";
      FCSys.WorkInProgress.Magnetic2.Effort Sigma_v2
        "Electrical sink force (for net current into each flange)";
      FCSys.WorkInProgress.Magnetic2.Flow Sigma_qdot2
        "Net current stored or destroyed (positive into each pin)";

      FCSys.Systems.DCDC.Interfaces.ActBusIn actBusIn annotation (Placement(
            transformation(extent={{-62,-10},{-42,10}}), iconTransformation(
              extent={{-52,-10},{-32,10}})));
      FCSys.Systems.DCDC.Interfaces.SenBusOut senBusOut annotation (Placement(
            transformation(extent={{42,-10},{62,10}}), iconTransformation(
              extent={{32,-10},{52,10}})));
      FCSys.Connectors.FaceBus wireN2 annotation (Placement(transformation(
              extent={{-106,54},{-94,66}}), iconTransformation(extent={{-110,50},
                {-90,70}})));
      FCSys.Connectors.FaceBus wireP2 annotation (Placement(transformation(
              extent={{94,54},{106,66}}), iconTransformation(extent={{90,50},{
                110,70}})));
      FCSys.Connectors.FaceBus wireN1 annotation (Placement(transformation(
              extent={{-106,-66},{-94,-54}}), iconTransformation(extent={{-110,
                -70},{-90,-50}})));
      FCSys.Connectors.FaceBus wireP1 annotation (Placement(transformation(
              extent={{94,-66},{106,-54}}), iconTransformation(extent={{90,-70},
                {110,-50}})));
    equation
      Delta_v1 = chargeP1.v - chargeN1.v;
      Delta_qdot1 = chargeP1.qdot - chargeN1.qdot;
      Sigma_v1 = chargeP1.v + chargeN1.v;
      Sigma_qdot1 = chargeP1.qdot + chargeN1.qdot;

      Delta_v2 = chargeP2.v - chargeN2.v;
      Delta_qdot2 = chargeP2.qdot - chargeN2.qdot;
      Sigma_v2 = chargeP2.v + chargeN2.v;
      Sigma_qdot2 = chargeP2.qdot + chargeN2.qdot;

      // Steady state charge balance (no current destroyed or stored)
      0 = Sigma_qdot1;
      0 = Sigma_qdot2;

      // Steady state energy balance
      0 = Delta_v1*Delta_qdot1 + Delta_v2*Delta_qdot2 + Sigma_v1*Sigma_qdot1 +
        Sigma_v2*Sigma_qdot2;

      // Internal connector
      if specAsParam then
        if (LHS == RHS) then
          k_internal = k_ratio;
        elseif (LHS == 1 and RHS == 2) then
          k_internal = k_R;
        elseif (LHS == 2 and RHS == 1) then
          k_internal = k_sigma;
        elseif (LHS == 1 and RHS == 3) then
          k_internal = k_Delta_v1;
        elseif (LHS == 2 and RHS == 3) then
          k_internal = k_Delta_qdot1;
        end if;
      end if;

      // Governing linear equation
      if (LHS == BaseClasses.LHSchoice.Delta_v1 and RHS == BaseClasses.RHSchoice.Delta_v2)
           then
        Delta_v1 = k_internal*Delta_v2;
      elseif (LHS == BaseClasses.LHSchoice.Delta_v1 and RHS == BaseClasses.RHSchoice.Delta_qdot2)
           then
        Delta_v1 = k_internal*Delta_qdot2;
      elseif (LHS == BaseClasses.LHSchoice.Delta_qdot1 and RHS == BaseClasses.RHSchoice.Delta_v2)
           then
        Delta_qdot1 = k_internal*Delta_v2;
      elseif (LHS == BaseClasses.LHSchoice.Delta_qdot1 and RHS == BaseClasses.RHSchoice.Delta_qdot2)
           then
        Delta_qdot1 = k_internal*Delta_qdot2;
      elseif (LHS == BaseClasses.LHSchoice.Delta_v1 and RHS == BaseClasses.RHSchoice.unity)
           then
        Delta_v1 = k_internal;
      elseif (LHS == BaseClasses.LHSchoice.Delta_qdot1 and RHS == BaseClasses.RHSchoice.unity)
           then
        Delta_qdot1 = k_internal;
      end if;

      annotation (Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            grid={2,2}), graphics={Rectangle(
                  extent={{-60,20},{60,-20}},
                  fillColor={0,0,0},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None,
                  lineColor={0,0,0}),Polygon(
                  points={{-40,20},{80,20},{20,-40},{-40,20}},
                  lineColor={0,192,0},
                  smooth=Smooth.None,
                  origin={20,-40},
                  rotation=180),Polygon(
                  points={{-40,20},{80,20},{20,-40},{-40,20}},
                  lineColor={0,192,0},
                  smooth=Smooth.None,
                  origin={-20,40},
                  rotation=360),Line(
                  points={{-90,60},{90,60}},
                  color={0,192,0},
                  smooth=Smooth.None),Line(
                  points={{-90,-60},{90,-60}},
                  color={0,192,0},
                  smooth=Smooth.None),Polygon(
                  points={{-20,0},{100,0},{40,-60},{-20,0}},
                  lineColor={0,192,0},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  origin={-40,60},
                  rotation=360),Polygon(
                  points={{-40,-6.89683e-15},{80,7.13232e-15},{20,-60},{-40,-6.89683e-15}},
                  lineColor={0,192,0},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  origin={20,-60},
                  rotation=180)}), Diagram(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            grid={2,2}), graphics));
    end DCDC;

    model IdealDCDC
      extends FCSys.Systems.DCDC.ElecElec(
        final LHS=BaseClasses.LHSchoice.Delta_v1,
        final RHS=BaseClasses.RHSchoice.unity,
        final specAsParam=false,
        final k_ratio=1,
        final k_R=1,
        final k_sigma=1,
        final k_Delta_v1=1,
        final k_Delta_qdot1=1);

      annotation (Diagram(graphics), Icon(graphics));
    end IdealDCDC;

    model ElecElec
      extends FCSys.WorkInProgress.TwoElecStreams;
      extends FCSys.WorkInProgress.BaseClassesIconsTransformer;

      parameter FCSys.Systems.BaseClasses.LHSchoice LHS=BaseClasses.LHSchoice.Delta_v1
        "Choice of variable for the left hand side (LHS) of the governing equation"
        annotation (Evaluate=true);
      parameter FCSys.Systems.BaseClasses.RHSchoice RHS=BaseClasses.RHSchoice.Delta_v2
        "Choice of variable for the right hand side (RHS) of the governing equation"
        annotation (Evaluate=true);
      parameter Boolean specAsParam=true
        "true: use a parameter to specify the governing equation; false: use dynamic input";

      parameter Q.Unity k_ratio=1 "Voltage or current ratio of side 1:side 2"
        annotation (Dialog(enable=((LHS == RHS) and specAsParam)));
      parameter FCSys.WorkInProgress.Magnetic.Conductance k_R=1*U.ohm
        "Hall resistance across side 1 (Delta_v1/Delta_qdot2)"
        annotation (Dialog(enable=((LHS == 1) and (RHS == 2) and specAsParam)));
      parameter FCSys.WorkInProgress.Magnetic2.Conductance k_sigma=1*U.S
        "Hall conductance through side 1 (Delta_qdot1/Delta_v2)"
        annotation (Dialog(enable=((LHS == 2) and (RHS == 1) and specAsParam)));
      parameter FCSys.WorkInProgress.Magnetic2.Effort k_Delta_v1=1*U.V
        "Voltage across side 1"
        annotation (Dialog(enable=((LHS == 1) and (RHS == 3) and specAsParam)));
      parameter FCSys.WorkInProgress.Magnetic2.QuantumRate k_Delta_qdot1=1*U.A
        "Total current through side 1"
        annotation (Dialog(enable=((LHS == 2) and (RHS == 3) and specAsParam)));

      FCSys.WorkInProgress.Magnetic2.Effort Delta_v1
        "Electrical driving force (for net current positive from elecP to elecN)";
      FCSys.WorkInProgress.Magnetic2.QuantumRate Delta_qdot1
        "Total current of both interfaces (positive by right hand rule from elecP to elecN)";
      FCSys.WorkInProgress.Magnetic2.Effort Sigma_v1
        "Electrical sink force (for net current into each interface)";
      FCSys.WorkInProgress.Magnetic2.QuantumRate Sigma_qdot1
        "Net current stored or destroyed (positive into each interface)";

      FCSys.WorkInProgress.Magnetic2.Effort Delta_v2
        "Electrical driving force (for net current positive from elecP to elecN)";
      FCSys.WorkInProgress.Magnetic2.QuantumRate Delta_qdot2
        "Total current of both interfaces (positive by right hand rule from elecP to elecN)";
      FCSys.WorkInProgress.Magnetic2.Effort Sigma_v2
        "Electrical sink force (for net current into each interface)";
      FCSys.WorkInProgress.Magnetic2.QuantumRate Sigma_qdot2
        "Net current stored or destroyed (positive into each interface)";

      FCSys.Connectors.RealInput k_set(final unit=if (LHS == RHS) then "1"
             else if (LHS == 1 and RHS == 2) then "1/Q2" else if (LHS == 2 and
            RHS == 1) then "Q2" else if (LHS == 1 and RHS == 3) then "1/(Q.T)"
             else "Q/T", displayUnit=if (LHS == RHS) then "1" else if (LHS == 1
             and RHS == 2) then "ohm" else if (LHS == 2 and RHS == 1) then "S"
             else if (LHS == 1 and RHS == 3) then "V" else "A") if (not
        specAsParam)
        "Coefficient for the RHS of the governing equation as an input"
        annotation (Placement(transformation(
            origin={0,-120},
            extent={{-20,-20},{20,20}},
            rotation=90), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=90,
            origin={0,-80})));

    protected
      FCSys.Connectors.RealInput k_internal(final unit=if (LHS == RHS) then "1"
             else if (LHS == 1 and RHS == 2) then "1/Q2" else if (LHS == 2 and
            RHS == 1) then "Q2" else if (LHS == 1 and RHS == 3) then "1/(Q.T)"
             else "Q/T", displayUnit=if (LHS == RHS) then "1" else if (LHS == 1
             and RHS == 2) then "ohm" else if (LHS == 2 and RHS == 1) then "S"
             else if (LHS == 1 and RHS == 3) then "V" else "A")
        "Coefficient for the RHS of the governing equation as an internal connector"
        annotation (Placement(transformation(
            origin={0,-80},
            extent={{-20,-20},{20,20}},
            rotation=90)));

    equation
      Delta_v1 = elecP1.Phidot_d - elecN1.Phidot_d;
      Delta_qdot1 = elecP1.qdot_t - elecN1.qdot_t;
      Sigma_v1 = elecP1.Phidot_d + elecN1.Phidot_d;
      Sigma_qdot1 = elecP1.qdot_t + elecN1.qdot_t;

      Delta_v2 = elecP2.Phidot_d - elecN2.Phidot_d;
      Delta_qdot2 = elecP2.qdot_t - elecN2.qdot_t;
      Sigma_v2 = elecP2.Phidot_d + elecN2.Phidot_d;
      Sigma_qdot2 = elecP2.qdot_t + elecN2.qdot_t;

      // Steady state charge balance (no current destroyed or stored)
      0 = Sigma_qdot1;
      0 = Sigma_qdot2;

      // Steady state energy balance
      0 = Delta_v1*Delta_qdot1 + Delta_v2*Delta_qdot2 + Sigma_v1*Sigma_qdot1 +
        Sigma_v2*Sigma_qdot2;

      // Internal connector
      if specAsParam then
        if (LHS == RHS) then
          k_internal = k_ratio;
        elseif (LHS == 1 and RHS == 2) then
          k_internal = k_R;
        elseif (LHS == 2 and RHS == 1) then
          k_internal = k_sigma;
        elseif (LHS == 1 and RHS == 3) then
          k_internal = k_Delta_v1;
        elseif (LHS == 2 and RHS == 3) then
          k_internal = k_Delta_qdot1;
        end if;
      end if;

      // Governing linear equation
      if (LHS == BaseClasses.LHSchoice.Delta_v1 and RHS == BaseClasses.RHSchoice.Delta_v2)
           then
        Delta_v1 = k_internal*Delta_v2;
      elseif (LHS == BaseClasses.LHSchoice.Delta_v1 and RHS == BaseClasses.RHSchoice.Delta_qdot2)
           then
        Delta_v1 = k_internal*Delta_qdot2;
      elseif (LHS == BaseClasses.LHSchoice.Delta_qdot1 and RHS == BaseClasses.RHSchoice.Delta_v2)
           then
        Delta_qdot1 = k_internal*Delta_v2;
      elseif (LHS == BaseClasses.LHSchoice.Delta_qdot1 and RHS == BaseClasses.RHSchoice.Delta_qdot2)
           then
        Delta_qdot1 = k_internal*Delta_qdot2;
      elseif (LHS == BaseClasses.LHSchoice.Delta_v1 and RHS == BaseClasses.RHSchoice.unity)
           then
        Delta_v1 = k_internal;
      elseif (LHS == BaseClasses.LHSchoice.Delta_qdot1 and RHS == BaseClasses.RHSchoice.unity)
           then
        Delta_qdot1 = k_internal;
      end if;

      connect(k_set, k_internal) annotation (Line(
          points={{1.11022e-15,-120},{1.11022e-15,-110},{1.11022e-15,-110},{
              1.11022e-15,-100},{1.11022e-15,-80},{1.11022e-15,-80}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation (Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            grid={2,2}), graphics={Polygon(
                  points={{-40,20},{80,20},{20,-40},{-40,20}},
                  lineColor={0,192,0},
                  smooth=Smooth.None,
                  origin={20,-40},
                  rotation=180),Polygon(
                  points={{-40,20},{80,20},{20,-40},{-40,20}},
                  lineColor={0,192,0},
                  smooth=Smooth.None,
                  origin={-20,40},
                  rotation=360),Line(
                  points={{-90,60},{90,60}},
                  color={0,192,0},
                  smooth=Smooth.None),Line(
                  points={{-90,-60},{90,-60}},
                  color={0,192,0},
                  smooth=Smooth.None)}), Diagram(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            grid={2,2}), graphics));
    end ElecElec;

    package Interfaces
      extends Modelica.Icons.InterfacesPackage;

      expandable connector ActBusIn
        "Data bus for inputs to a plant's actuators"
        extends FCSys.BaseClasses.Icons.SignalBuses.In;
        input SI.Voltage v_ref;
        input Real valvePos_an_in;
        input Real valvePos_an_out;
        input Real valvePos_ca_out;
        input Real valvePos_humid_an;
        input Real valvePos_humid_ca;
        input Real qdot_compressor;
        input Real qdot_H2Opump;
        input Real qdot_coolingPump;

        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics));
      end ActBusIn;

      expandable connector ActBusOut
        "Data bus for inputs to a plant's actuators"
        extends FCSys.BaseClasses.Icons.SignalBuses.Out;
        output SI.Current current_ref;
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics));
      end ActBusOut;

      expandable connector SenBusIn
        "Data bus for outputs from a plant's sensors"
        extends FCSys.BaseClasses.Icons.SignalBuses.In;
        input Q.Voltage current;
        input Q.Rotation angle;
        input Q.RotationalVelocity speed;
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics));
      end SenBusIn;

      expandable connector SenBusOut
        "Data bus for outputs from a plant's sensors"
        extends FCSys.BaseClasses.Icons.SignalBuses.Out;
        output Q.Voltage current;
        output Q.Rotation angle;
        output Q.RotationalVelocity speed;
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics));
      end SenBusOut;

      model ActArrayToBus

        FCSys.Systems.DCDC.Interfaces.ActBusOut actBusOut annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={100,0}),iconTransformation(
              extent={{-80,-20},{-40,20}},
              rotation=0,
              origin={100,0})));
        FCSys.Connectors.RealInput actIN[1] annotation (Placement(
              transformation(extent={{-110,-10},{-90,10}}), iconTransformation(
                extent={{-60,-20},{-20,20}})));
      equation

        connect(actIN[1], actBusOut.current_ref) annotation (Line(
            points={{-100,0},{-4,0},{-4,5.55112e-16},{100,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{6,3},{6,3}}));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-100},{100,100}}), graphics), Icon(coordinateSystem(
                preserveAspectRatio=true, extent={{-60,-40},{60,40}}), graphics
              ={Rectangle(
                      extent={{-40,40},{40,-40}},
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Solid,
                      pattern=LinePattern.Dash,
                      lineColor={0,0,0}),Line(
                      points={{-28,4},{28,4}},
                      color={0,0,127},
                      smooth=Smooth.None),Line(
                      points={{-28,-4},{34,-4}},
                      color={0,0,127},
                      smooth=Smooth.None)}));
      end ActArrayToBus;

      model SenBusToArray

        SenBusIn senBusIn annotation (Placement(transformation(extent={{-110,-10},
                  {-90,10}}), iconTransformation(extent={{-60,-20},{-20,20}})));
        FCSys.Connectors.RealOutput senOut[3] annotation (Placement(
              transformation(extent={{90,-10},{110,10}}), iconTransformation(
                extent={{20,-20},{60,20}})));
      equation

        connect(senBusIn.angle, senOut[3]) annotation (Line(
            points={{-100,5.55112e-16},{0,5.55112e-16},{0,6.66667},{100,6.66667}},

            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        connect(senBusIn.speed, senOut[2]) annotation (Line(
            points={{-100,5.55112e-16},{0,5.55112e-16},{0,0},{100,0}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));
        connect(senBusIn.current, senOut[1]) annotation (Line(
            points={{-100,5.55112e-16},{0,5.55112e-16},{0,-6.66667},{100,-6.66667}},

            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-100},{100,100}}), graphics), Icon(coordinateSystem(
                preserveAspectRatio=true, extent={{-60,-40},{60,40}}), graphics
              ={Rectangle(
                      extent={{-40,40},{40,-40}},
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Solid,
                      pattern=LinePattern.Dash,
                      lineColor={0,0,0}),Line(
                      points={{-34,-4},{20,-4}},
                      color={0,0,127},
                      smooth=Smooth.None),Line(
                      points={{-28,4},{20,4}},
                      color={0,0,127},
                      smooth=Smooth.None)}));
      end SenBusToArray;
    end Interfaces;
  end DCDC;

  package Pump
    extends Modelica.Icons.Package;
    model Pump
      // TODO: Complete this.

      extends FCSys.WorkInProgress.BaseClassesIconsPump;
      // extends FCSys.Processes.BaseClasses.PartialTransport.Middle;
      // extends FCSys.Processes.BaseClasses.PartialTransport.Bottom2;

      FCSys.Connectors.FaceBus wireN annotation (Placement(transformation(
              extent={{-30,-10},{-10,10}}), iconTransformation(extent={{-110,-50},
                {-90,-30}})));
      FCSys.Connectors.FaceBus wireP annotation (Placement(transformation(
              extent={{10,-10},{30,10}}), iconTransformation(extent={{90,-50},{
                110,-30}})));
      FCSys.Systems.Pump.Interfaces.ActBusIn actBusIn annotation (Placement(
            transformation(extent={{-42,-70},{-22,-50}}), iconTransformation(
              extent={{-32,-90},{-12,-70}})));
      FCSys.Systems.Pump.Interfaces.SenBusOut senBusOut annotation (Placement(
            transformation(extent={{2,-70},{22,-50}}), iconTransformation(
              extent={{10,-90},{30,-70}})));
      FCSys.Connectors.FaceBus pipeP annotation (Placement(transformation(
              extent={{94,-6},{106,6}}), iconTransformation(extent={{90,-10},{
                110,10}})));
      FCSys.Connectors.FaceBus pipeN annotation (Placement(transformation(
              extent={{-106,-6},{-94,6}}), iconTransformation(extent={{-110,-10},
                {-90,10}})));
    equation

      annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}), graphics), Icon(coordinateSystem(
              preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
            graphics={Line(
                  points={{-90,-40},{-30,-40}},
                  color={0,200,0},
                  smooth=Smooth.None),Line(
                  points={{30,-40},{90,-40}},
                  color={0,200,0},
                  smooth=Smooth.None),Line(
                  points={{20,0},{90,0}},
                  color={0,0,0},
                  smooth=Smooth.None),Line(
                  points={{-90,0},{-20,0}},
                  color={0,0,0},
                  smooth=Smooth.None)}));
    end Pump;

    package Interfaces
      extends Modelica.Icons.InterfacesPackage;
      expandable connector ActBusIn
        "Data bus for inputs to a plant's actuators"
        extends FCSys.BaseClasses.Icons.SignalBuses.In;
        input SI.Current current_ref;
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics));
      end ActBusIn;

      expandable connector ActBusOut
        "Data bus for inputs to a plant's actuators"
        extends FCSys.BaseClasses.Icons.SignalBuses.Out;
        output SI.Current current_ref;
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics));
      end ActBusOut;

      expandable connector SenBusIn
        "Data bus for outputs from a plant's sensors"
        extends FCSys.BaseClasses.Icons.SignalBuses.In;
        input Q.Voltage current;
        input Q.Rotation angle;
        input Q.RotationalVelocity speed;
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics));
      end SenBusIn;

      expandable connector SenBusOut
        "Data bus for outputs from a plant's sensors"
        extends FCSys.BaseClasses.Icons.SignalBuses.Out;
        output Q.Voltage current;
        output Q.Rotation angle;
        output Q.RotationalVelocity speed;
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics));
      end SenBusOut;

      model ActArrayToBus

        FCSys.Systems.Pump.Interfaces.ActBusOut actBusOut annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={100,0}),iconTransformation(
              extent={{-80,-20},{-40,20}},
              rotation=0,
              origin={100,0})));
        FCSys.Connectors.RealInput actIN[1] annotation (Placement(
              transformation(extent={{-110,-10},{-90,10}}), iconTransformation(
                extent={{-60,-20},{-20,20}})));
      equation

        connect(actIN[1], actBusOut.current_ref) annotation (Line(
            points={{-100,0},{-4,0},{-4,5.55112e-16},{100,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{6,3},{6,3}}));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-100},{100,100}}),graphics), Icon(coordinateSystem(
                preserveAspectRatio=true, extent={{-60,-40},{60,40}}), graphics
              ={Rectangle(
                      extent={{-40,40},{40,-40}},
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Solid,
                      pattern=LinePattern.Dash,
                      lineColor={0,0,0}),Line(
                      points={{-28,4},{28,4}},
                      color={0,0,127},
                      smooth=Smooth.None),Line(
                      points={{-28,-4},{34,-4}},
                      color={0,0,127},
                      smooth=Smooth.None)}));
      end ActArrayToBus;

      model SenBusToArray

        SenBusIn senBusIn annotation (Placement(transformation(extent={{-110,-10},
                  {-90,10}}), iconTransformation(extent={{-60,-20},{-20,20}})));
        FCSys.Connectors.RealOutput senOut[3] annotation (Placement(
              transformation(extent={{90,-10},{110,10}}), iconTransformation(
                extent={{20,-20},{60,20}})));
      equation

        connect(senBusIn.angle, senOut[3]) annotation (Line(
            points={{-100,5.55112e-16},{0,5.55112e-16},{0,6.66667},{100,6.66667}},

            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        connect(senBusIn.speed, senOut[2]) annotation (Line(
            points={{-100,5.55112e-16},{0,5.55112e-16},{0,0},{100,0}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));
        connect(senBusIn.current, senOut[1]) annotation (Line(
            points={{-100,5.55112e-16},{0,5.55112e-16},{0,-6.66667},{100,-6.66667}},

            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-100},{100,100}}),graphics), Icon(coordinateSystem(
                preserveAspectRatio=true, extent={{-60,-40},{60,40}}), graphics
              ={Rectangle(
                      extent={{-40,40},{40,-40}},
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Solid,
                      pattern=LinePattern.Dash,
                      lineColor={0,0,0}),Line(
                      points={{-34,-4},{20,-4}},
                      color={0,0,127},
                      smooth=Smooth.None),Line(
                      points={{-28,4},{20,4}},
                      color={0,0,127},
                      smooth=Smooth.None)}));
      end SenBusToArray;
    end Interfaces;
  end Pump;

  package Valve
    extends Modelica.Icons.Package;
    model Valve "Valve for water/steam flows with linear pressure drop"
      // TODO: Complete this.

      //  extends FCSys.Processes.BaseClasses.PartialTransport.Top1;
      extends FCSys.BaseClasses.Icons.Name.Top3;

      FCSys.Systems.Valve.Interfaces.ActBusIn actBusIn annotation (Placement(
            transformation(extent={{-50,-70},{-30,-50}}), iconTransformation(
              extent={{-32,-30},{-12,-10}})));
      FCSys.Systems.Valve.Interfaces.SenBusOut senBusOut annotation (Placement(
            transformation(extent={{30,-70},{50,-50}}), iconTransformation(
              extent={{12,-30},{32,-10}})));
      FCSys.Connectors.FaceBus pipeN annotation (Placement(transformation(
              extent={{-106,14},{-94,26}}), iconTransformation(extent={{-110,10},
                {-90,30}})));
      FCSys.Connectors.FaceBus pipeP annotation (Placement(transformation(
              extent={{94,14},{106,26}}), iconTransformation(extent={{90,10},{
                110,30}})));
      annotation (
        Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            grid={2,2}), graphics={Line(points={{0,20},{0,-30}}, color={0,0,0}),
              Rectangle(
                  extent={{-40,0},{40,-40}},
                  lineColor={0,0,0},
                  fillColor={0,0,0},
                  fillPattern=FillPattern.Solid),Polygon(
                  points={{-40,40},{40,0},{40,40},{0,20},{-40,0},{-40,40}},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),Line(
                  points={{-90,20},{-40,20}},
                  color={0,0,0},
                  smooth=Smooth.None),Line(
                  points={{40,20},{90,20}},
                  color={0,0,0},
                  smooth=Smooth.None)}),
        Diagram(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            grid={2,2}), graphics),
        Documentation(info="<html>
<p>This very simple model provides a pressure drop which is proportional to the flow rate and to the <code>face</code> input, without computing any fluid property. It can be used for testing purposes, when
a simple model of a variable pressure loss is needed.</p>
<p>A medium model must be nevertheless be specified, so that the fluid ports can be connected to other species using the same medium model.</p>
<p>The model is adiabatic (no heat losses to the ambient) and neglects changes in kinetic energy from the inlet to the outlet. </p>
</html>", revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"));
    end Valve;

    package Interfaces
      // TODO: Complete this.

      extends Modelica.Icons.InterfacesPackage;
      expandable connector ActBusIn
        "Data bus for inputs to a plant's actuators"
        extends FCSys.BaseClasses.Icons.SignalBuses.In;
        input SI.Voltage v_ref;
        input Real valvePos_an_in;
        input Real valvePos_an_out;
        input Real valvePos_ca_out;
        input Real valvePos_humid_an;
        input Real valvePos_humid_ca;
        input Real qdot_compressor;
        input Real qdot_H2Opump;
        input Real qdot_coolingPump;

        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics));
      end ActBusIn;

      expandable connector ActBusOut
        "Data bus for inputs to a plant's actuators"
        extends FCSys.BaseClasses.Icons.SignalBuses.Out;
        output SI.Current current_ref;
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics));
      end ActBusOut;

      expandable connector SenBusIn
        "Data bus for outputs from a plant's sensors"
        extends FCSys.BaseClasses.Icons.SignalBuses.In;
        input Q.Voltage current;
        input Q.Rotation angle;
        input Q.RotationalVelocity speed;
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics));
      end SenBusIn;

      expandable connector SenBusOut
        "Data bus for outputs from a plant's sensors"
        extends FCSys.BaseClasses.Icons.SignalBuses.Out;
        output Q.Voltage current;
        output Q.Rotation angle;
        output Q.RotationalVelocity speed;
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics));
      end SenBusOut;

      model ActArrayToBus

        FCSys.Systems.Valve.Interfaces.ActBusOut actBusOut annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={100,0}),iconTransformation(
              extent={{-80,-20},{-40,20}},
              rotation=0,
              origin={100,0})));
        FCSys.Connectors.RealInput actIN[1] annotation (Placement(
              transformation(extent={{-110,-10},{-90,10}}), iconTransformation(
                extent={{-60,-20},{-20,20}})));
      equation

        connect(actIN[1], actBusOut.current_ref) annotation (Line(
            points={{-100,0},{-4,0},{-4,5.55112e-16},{100,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{6,3},{6,3}}));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-100},{100,100}}), graphics), Icon(coordinateSystem(
                preserveAspectRatio=true, extent={{-60,-40},{60,40}}), graphics
              ={Rectangle(
                      extent={{-40,40},{40,-40}},
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Solid,
                      pattern=LinePattern.Dash,
                      lineColor={0,0,0}),Line(
                      points={{-28,4},{28,4}},
                      color={0,0,127},
                      smooth=Smooth.None),Line(
                      points={{-28,-4},{34,-4}},
                      color={0,0,127},
                      smooth=Smooth.None)}));
      end ActArrayToBus;

      model SenBusToArray

        SenBusIn senBusIn annotation (Placement(transformation(extent={{-110,-10},
                  {-90,10}}), iconTransformation(extent={{-60,-20},{-20,20}})));
        FCSys.Connectors.RealOutput senOut[3] annotation (Placement(
              transformation(extent={{90,-10},{110,10}}), iconTransformation(
                extent={{20,-20},{60,20}})));
      equation

        connect(senBusIn.angle, senOut[3]) annotation (Line(
            points={{-100,5.55112e-16},{0,5.55112e-16},{0,6.66667},{100,6.66667}},

            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        connect(senBusIn.speed, senOut[2]) annotation (Line(
            points={{-100,5.55112e-16},{0,5.55112e-16},{0,0},{100,0}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));
        connect(senBusIn.current, senOut[1]) annotation (Line(
            points={{-100,5.55112e-16},{0,5.55112e-16},{0,-6.66667},{100,-6.66667}},

            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={
                  {-100,-100},{100,100}}), graphics), Icon(coordinateSystem(
                preserveAspectRatio=true, extent={{-60,-40},{60,40}}), graphics
              ={Rectangle(
                      extent={{-40,40},{40,-40}},
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Solid,
                      pattern=LinePattern.Dash,
                      lineColor={0,0,0}),Line(
                      points={{-34,-4},{20,-4}},
                      color={0,0,127},
                      smooth=Smooth.None),Line(
                      points={{-28,4},{20,4}},
                      color={0,0,127},
                      smooth=Smooth.None)}));
      end SenBusToArray;
    end Interfaces;
  end Valve;

  package FluidHeater
    extends Modelica.Icons.Package;
    model FluidHeater
      // TODO: Complete this.

      extends FCSys.BaseClasses.Icons.Name.Top4;

      FCSys.Subregions.HeatExchanger fCConvection annotation (Placement(
            transformation(extent={{90,116},{110,136}}, rotation=0)));
      FCSys.Connectors.FaceBus wireN annotation (Placement(transformation(
              extent={{-110,-50},{-90,-30}}), iconTransformation(extent={{-110,
                -50},{-90,-30}})));
      FCSys.Connectors.FaceBus wireP annotation (Placement(transformation(
              extent={{90,-50},{110,-30}}), iconTransformation(extent={{90,-50},
                {110,-30}})));
      FCSys.Connectors.FaceBus pipeP annotation (Placement(transformation(
              extent={{90,30},{110,50}}), iconTransformation(extent={{90,30},{
                110,50}})));
      FCSys.Connectors.FaceBus pipeN annotation (Placement(transformation(
              extent={{-110,30},{-90,50}}), iconTransformation(extent={{-110,30},
                {-90,50}})));
      FCSys.Systems.Valve.Interfaces.ActBusIn actBusIn annotation (Placement(
            transformation(extent={{-50,-70},{-30,-50}}), iconTransformation(
              extent={{-34,-90},{-14,-70}})));
      FCSys.Systems.Valve.Interfaces.SenBusOut senBusOut annotation (Placement(
            transformation(extent={{30,-70},{50,-50}}), iconTransformation(
              extent={{10,-90},{30,-70}})));
    equation
      connect(fCConvection.matN1, fCConvection.matP1) annotation (Line(
          points={{90,122},{110,122}},
          color={0,200,0},
          smooth=Smooth.None,
          thickness=0.5));
      annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                -100},{100,100}}), graphics={Rectangle(
                  extent={{-72,-20},{72,-60}},
                  lineColor={215,215,215},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),Rectangle(
                  extent={{-72,60},{72,20}},
                  lineColor={215,215,215},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),Line(
                  points={{-90,40},{-72,40},{-60,60},{-36,20},{-12,60},{12,20},
                {36,60},{60,20},{72,40},{90,40}},
                  color={0,0,0},
                  smooth=Smooth.None),Rectangle(
                  extent={{-40,-60},{40,-100}},
                  fillColor={0,0,0},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Line(
                  points={{-38,-24},{-38,-24},{-38,-24},{-12,20},{12,-20},{36,
                20},{60,-20},{72,0}},
                  color={191,0,0},
                  smooth=Smooth.Bezier,
                  origin={0,-2},
                  rotation=90),Line(
                  points={{0,70},{4,56}},
                  color={191,0,0},
                  smooth=Smooth.Bezier),Line(
                  points={{0,70},{16,70}},
                  color={191,0,0},
                  smooth=Smooth.Bezier),Line(
                  points={{-90,-40},{-72,-40},{-60,-20},{-36,-60},{-12,-20},{12,
                -60},{36,-20},{60,-60},{72,-40},{90,-40}},
                  color={0,200,0},
                  smooth=Smooth.None)}), Diagram(coordinateSystem(
              preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
            graphics));
    end FluidHeater;
  end FluidHeater;

  package BaseClasses "Base classes (not for direct use)"
    extends Modelica.Icons.BasesPackage;

    type LHSchoice = enumeration(
        Delta_v1 "Voltage across side 1",
        Delta_qdot1 "Total current through side 1")
      "Enumeration defining the left hand side (LHS) variable of an electrical-electrical transformer"
      annotation (Evaluate=true);
    type RHSchoice = enumeration(
        Delta_v2 "Voltage across side 2",
        Delta_qdot2 "Total current through side 2",
        unity "Unity (independent of side 2)")
      "Enumeration defining the right hand side (RHS) variable of an electrical-electrical transformer"
      annotation (Evaluate=true);
  end BaseClasses;

  annotation (Documentation(info="<html><p>
<b>Licensed by Kevin Davies under the Modelica License 2</b><br>
Copyright 2007&ndash;2012, Kevin Davies.
</p>
<p>
<i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.ModelicaLicense2\">
FCSys.UsersGuide.ModelicaLicense2</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
http://www.modelica.org/licenses/ModelicaLicense2</a>.</i>
</p></html>"));
end Systems;
