within FCSys;
package WorkInProgress "Incomplete classes under development"
  extends Modelica.Icons.Package;
  extends Modelica.Icons.UnderConstruction;

  model ConditionsAdaptersPhasesIonomer
    "<html>Adapter for ionomer between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
    extends Conditions.Adapters.Phases.BaseClasses.PartialPhase;
    extends Modelica.Icons.UnderConstruction;
    Conditions.Adapters.Species.Solid 'SO3-';
    FCSys.WorkInProgress.ConditionsAdaptersSpeciesFluid 'H+'(redeclare package
        Data = FCSys.Characteristics.'H+'.Ionomer, redeclare package Medium =
          Modelica.Media.IdealGases.SingleGases.H2)
      annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));
    // **Use model for H instead.
    Conditions.Adapters.Species.FluidNeutral H2O(redeclare package Data =
          Characteristics.H2O.Ionomer, redeclare package Medium =
          Modelica.Media.IdealGases.SingleGases.H2O)
      annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
    Conditions.Adapters.Junctions.Junction2 junction2
      annotation (Placement(transformation(extent={{60,-50},{40,-30}})));
    Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final package
        Medium = Medium) "Modelica fluid port" annotation (Placement(
          transformation(extent={{70,-50},{90,-30}}), iconTransformation(extent
            ={{70,-50},{90,-30}})));
    Modelica.Electrical.Analog.Interfaces.NegativePin pin
      "Modelica electrical pin" annotation (Placement(transformation(extent={{
              70,30},{90,50}}), iconTransformation(extent={{70,30},{90,50}})));

  equation
    connect('SO3-'.boundary.thermal, boundary.'SO3-'.thermal) annotation (Line(
        points={{-8,20},{-40,20},{-40,5.55112e-16},{-80,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('H+'.boundary.normal, boundary.'H+'.normal) annotation (Line(
        points={{-8,-20},{-40,-20},{-40,5.55112e-16},{-80,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('H+'.boundary.thermal, boundary.'H+'.thermal) annotation (Line(
        points={{-8,-20},{-40,-20},{-40,5.55112e-16},{-80,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('H+'.pin, pin) annotation (Line(
        points={{8,-16},{60,-16},{60,40},{80,40}},
        color={0,0,255},
        smooth=Smooth.None));
    connect('H+'.heatPort, heatPort) annotation (Line(
        points={{8,-20},{40,-20},{40,5.55112e-16},{80,5.55112e-16}},
        color={191,0,0},
        smooth=Smooth.None));
    connect('SO3-'.heatPort, heatPort) annotation (Line(
        points={{8,20},{40,20},{40,5.55112e-16},{80,5.55112e-16}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(junction2.mixturePort, fluidPort) annotation (Line(
        points={{58,-40},{80,-40}},
        color={0,127,255},
        smooth=Smooth.None));
    annotation (Placement(transformation(extent={{-10,10},{10,30}})), Icon(
          graphics={Line(
              points={{0,60},{0,-60}},
              color={0,0,0},
              smooth=Smooth.None,
              pattern=LinePattern.Dash,
              thickness=0.5),Line(
              points={{0,0},{-80,0}},
              color={127,127,127},
              smooth=Smooth.None,
              thickness=0.5),Line(
              points={{0,40},{80,40}},
              color={0,0,255},
              smooth=Smooth.None),Line(
              points={{0,0},{80,0}},
              color={191,0,0},
              smooth=Smooth.None),Line(
              points={{0,-40},{80,-40}},
              color={0,127,255},
              smooth=Smooth.None)}));
  end ConditionsAdaptersPhasesIonomer;

  model CellModelica
    "<html>Cell interfaced to components from the <a href=\"modelica://Modelica\">Modelica</a> package</html>"
    extends Modelica.Icons.Example;
    extends Modelica.Icons.UnderConstruction;
    Assemblies.Cells.Examples.Cell cell(anFP(redeclare
          FCSys.Subregions.Subregion subregions(
          each final inclX=true,
          each inclY=true,
          each graphite('incle-'=true, 'e-'(perfectMaterialDiff={{{{true,false}}}})),

          each gas(inclH2=true, inclH2O=true))))
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

    inner FCSys.Conditions.Environment environment(analysis=false)
      annotation (Placement(transformation(extent={{40,60},{60,80}})));
    Conditions.Adapters.Phases.Graphite caModelicaAdapt(A=cell.L_y[1]*cell.L_z[
          1]) annotation (Placement(transformation(extent={{20,-10},{40,10}})));
    Conditions.Adapters.Phases.Graphite anModelicaAdapt(A=cell.L_y[1]*cell.L_z[
          1])
      annotation (Placement(transformation(extent={{-20,-10},{-40,10}})));
    FCSys.WorkInProgress.TanConduct tanConduct
      annotation (Placement(transformation(extent={{10,40},{-10,60}})));
    Modelica.Blocks.Sources.Ramp loadSweep(duration=1000)
      "This is the arctangent of conductance."
      annotation (Placement(transformation(extent={{-40,60},{-20,80}})));

  equation
    connect(anModelicaAdapt.normal, cell.anFPX[1, 1]) annotation (Line(
        points={{-20,6.10623e-16},{-16,6.10623e-16},{-16,5.55112e-16},{-10,
            5.55112e-16}},
        color={0,0,0},
        thickness=0.5,
        smooth=Smooth.None));

    connect(caModelicaAdapt.normal, cell.caFPX[1, 1]) annotation (Line(
        points={{20,6.10623e-16},{16,6.10623e-16},{16,5.55112e-16},{10,
            5.55112e-16}},
        color={0,0,0},
        thickness=0.5,
        smooth=Smooth.None));

    connect(loadSweep.y, tanConduct.atanGstar) annotation (Line(
        points={{-19,70},{-6.66134e-16,70},{-6.66134e-16,61}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(tanConduct.p, caModelicaAdapt.pin) annotation (Line(
        points={{10,50},{60,50},{60,4},{40,4}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(tanConduct.n, anModelicaAdapt.pin) annotation (Line(
        points={{-10,50},{-60,50},{-60,4},{-40,4}},
        color={0,0,255},
        smooth=Smooth.None));
    annotation (experiment(StopTime=1000));
  end CellModelica;

  function plot "Create plots using FCRes"
    extends Modelica.Icons.Function;
    extends Modelica.Icons.UnderConstruction;

  algorithm
    Modelica.Utilities.System.command("loadres");

  end plot;

  model EISPlaceholder
    "Placeholder model for electrochemical-impedance spectroscopy"
    extends FCSys.Icons.Blocks.Continuous;

    parameter Modelica.SIunits.CurrentDensity zJ_large_SI=0.01
      "Large-signal current density in SI base units";
    Modelica.Electrical.Analog.Basic.Resistor resistor2(R=0.1) annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={10,-10})));
    Modelica.Electrical.Analog.Basic.Capacitor capacitor(C=1e-3) annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={40,-10})));
    Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V=1)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-80,-10})));
    Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor annotation (
        Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=270,
          origin={70,-10})));
    Modelica.Electrical.Analog.Sources.ConstantCurrent constantCurrent(I=
          zJ_large_SI)
      annotation (Placement(transformation(extent={{-60,20},{-40,0}})));
    Modelica.Electrical.Analog.Sources.SignalCurrent signalCurrent
      annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
    Connectors.RealInput zJ_small_SI
      "Small-signal cell current density in SI base units" annotation (
        Placement(transformation(extent={{-110,40},{-90,60}}),
          iconTransformation(extent={{-120,-10},{-100,10}})));
    Connectors.RealOutput w_V "Cell potential in volts" annotation (Placement(
          transformation(extent={{90,-20},{110,0}}), iconTransformation(extent=
              {{100,-10},{120,10}})));
    Modelica.Electrical.Analog.Basic.Resistor resistor1(R=0.1) annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-10,10})));
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{60,-60},{80,-40}})));
  equation
    connect(zJ_small_SI, signalCurrent.i) annotation (Line(
        points={{-100,50},{-50,50},{-50,37}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(signalCurrent.n, constantCurrent.n) annotation (Line(
        points={{-40,30},{-30,30},{-30,10},{-40,10}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(signalCurrent.p, constantCurrent.p) annotation (Line(
        points={{-60,30},{-70,30},{-70,10},{-60,10}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(constantVoltage.p, constantCurrent.p) annotation (Line(
        points={{-80,5.55112e-16},{-80,10},{-60,10}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(voltageSensor.v, w_V) annotation (Line(
        points={{80,-10},{90,-10},{90,-10},{100,-10}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(resistor2.p, capacitor.p) annotation (Line(
        points={{10,5.55112e-16},{10,10},{40,10},{40,5.55112e-16}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(constantVoltage.n, resistor2.n) annotation (Line(
        points={{-80,-20},{-80,-30},{10,-30},{10,-20}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(capacitor.n, resistor2.n) annotation (Line(
        points={{40,-20},{40,-30},{10,-30},{10,-20}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(capacitor.n, voltageSensor.n) annotation (Line(
        points={{40,-20},{40,-30},{70,-30},{70,-20}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(capacitor.p, voltageSensor.p) annotation (Line(
        points={{40,5.55112e-16},{40,10},{70,10},{70,5.55112e-16}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(constantCurrent.n, resistor1.p) annotation (Line(
        points={{-40,10},{-20,10}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(resistor1.n, resistor2.p) annotation (Line(
        points={{5.55112e-16,10},{10,10},{10,5.55112e-16}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(ground.p, voltageSensor.n) annotation (Line(
        points={{70,-40},{70,-20}},
        color={0,0,255},
        smooth=Smooth.None));
    annotation (Diagram(graphics), Icon(graphics));
  end EISPlaceholder;

  model ElectroOsmoticDrag
    "<html>Example to calibrate the coupling between H<sup>+</sup> and H<sub>2</sub>O in the PEM</html>"

    extends Regions.Examples.CLtoCL(anCL(redeclare model Subregion =
            Subregions.Subregion (ionomer(redeclare FCSys.Species.H2O.Gas.Fixed
                H2O(
                redeclare package Data = FCSys.Characteristics.'H+'.Ionomer,
                p_IC=65536*U.kPa,
                consMaterial=ConsThermo.IC)))), caCL(redeclare model Subregion
          = Subregions.Subregion (ionomer(redeclare FCSys.Species.H2O.Gas.Fixed
                H2O(
                redeclare package Data = FCSys.Characteristics.'H+'.Ionomer,
                p_IC=65536*U.kPa,
                consMaterial=ConsThermo.IC)))));

    extends Modelica.Icons.UnderConstruction;
  end ElectroOsmoticDrag;

  model ConditionsExamplesTestStandsTestStand
    "Fuel cell test stand (applies boundary conditions)"

    import FCSys.Utilities.average;
    import FCSys.Utilities.inSign;
    import FCSys.Characteristics.H2O.p_sat;
    import 'Datae-' = FCSys.Characteristics.'e-'.Graphite;
    extends FCSys.Icons.Names.Top9;

    // Geometry
    final parameter Integer n_x_an=1 "Along the through-cell axis in anode FP"
      annotation (Dialog(group="Number of subregions", __Dymola_label=
            "<html><i>n</i><sub>x an</sub></html>"));
    final parameter Integer n_x_ca=1
      "Along the through-cell axis in cathode FP" annotation (Dialog(group=
            "Number of subregions", __Dymola_label=
            "<html><i>n</i><sub>x ca</sub></html>"));
    final parameter Integer n_y=1 "Along the channel" annotation (Dialog(group=
            "Number of subregions", __Dymola_label=
            "<html><i>n</i><sub>y/sub></html>"));
    final parameter Integer n_z=1 "Across the channel" annotation (Dialog(group
          ="Number of subregions", __Dymola_label=
            "<html><i>n</i><sub>z</sub></html>"));

    // Operating conditions
    // --------------------
    // Electrical
    parameter Assemblies.Cells.Examples.Enumerations.ElectricalSpec
      electricalSpec=ElectricalSpec.current "Type of specification"
      annotation (Dialog(tab="Electrical"));
    Connectors.RealInput u_electrical=U.A/U.cm^2 "Value of the specification"
      annotation (Dialog(tab="Electrical", __Dymola_label=
            "<html><i>u</i><sub>electrical</sub></html>"));

    Q.Current zI "Current";
    Q.Potential w "Voltage";
    Q.ResistanceElectrical R "Resistance";
    Q.Power P "Power";
    //
    // General anode conditions
    parameter Side anInletSide=Side.p "Side of the inlet"
      annotation (Dialog(tab="Anode"));
    Q.TemperatureAbsolute T_an_in=333.15*U.K "Inlet temperature" annotation (
        Dialog(tab="Anode", __Dymola_label=
            "<html><i>T</i><sub>an in</sub></html>"));

    Q.PressureAbsolute p_an_out=U.from_kPag(48.3) "Outlet pressure" annotation
      (Dialog(tab="Anode", __Dymola_label=
            "<html><i>p</i><sub>an out</sub></html>"));

    //
    // General cathode conditions
    parameter Side caInletSide=Side.p "Side of the inlet"
      annotation (Dialog(tab="Cathode"));
    Q.TemperatureAbsolute T_ca_in=333.15*U.K "Inlet temperature" annotation (
        Dialog(tab="Cathode", __Dymola_label=
            "<html><i>T</i><sub>ca in</sub></html>"));

    Q.PressureAbsolute p_ca_out=U.from_kPag(48.3) "Outlet pressure" annotation
      (Dialog(tab="Cathode", __Dymola_label=
            "<html><i>p</i><sub>ca out</sub></html>"));

    Q.NumberAbsolute psi_dry_O2_in(
      final max=1,
      displayUnit="%") = 0.20946
      "<html>Dry-gas concentration of O<sub>2</sub> at the inlet</html>"
      annotation (Dialog(tab="Cathode", __Dymola_label=
            "<html>&psi;<sub>dry O2 in</sub></html>"));

    //
    // Anode flow rate
    parameter FCSys.Assemblies.Cells.Examples.Enumerations.FlowSpec anFlowSpec=
        FlowSpec.stoich "Type of specification"
      annotation (Dialog(tab="Anode",group="Supply"));
    Connectors.RealInput u_an_flow=1.5 "Value of the specification" annotation
      (Dialog(
        tab="Anode",
        group="Supply",
        __Dymola_label="<html><i>u</i><sub>an flow</sub></html>"));

    Q.NumberAbsolute anStoich "Anode stoichiometric flow rate";
    Q.Current I_an "Equivalent current of anode supply";
    Q.PressureAbsolute p_an_in "Anode inlet pressure";
    //
    // Cathode flow rate
    parameter FCSys.Assemblies.Cells.Examples.Enumerations.FlowSpec caFlowSpec=
        FlowSpec.stoich "Type of specification"
      annotation (Dialog(tab="Cathode",group="Supply"));
    Connectors.RealInput u_ca_flow=2.0 "Value of the specification" annotation
      (Dialog(
        tab="Cathode",
        group="Supply",
        __Dymola_label="<html><i>u</i><sub>ca flow</sub></html>"));

    Q.NumberAbsolute caStoich "Cathode stoichiometric flow rate";
    Q.Current I_ca "Equivalent current of cathode supply";
    Q.PressureAbsolute p_ca_in "Cathode inlet pressure";
    //
    // Anode humidity
    parameter FCSys.Assemblies.Cells.Examples.Enumerations.HumiditySpec
      anHumiditySpec=HumiditySpec.relative "Type of specification"
      annotation (Dialog(tab="Anode",group="Humidity"));
    Connectors.RealInput u_an_humidity=0.8 "Value of the specification"
      annotation (Dialog(
        tab="Anode",
        group="Humidity",
        __Dymola_label="<html><i>u</i><sub>an humidity</sub></html>"));

    Q.NumberAbsolute anInletRH(displayUnit="%")
      "Relative humidity at anode inlet";
    Q.PressureAbsolute p_H2O_an_in "H2O vapor pressure at anode inlet";
    Q.TemperatureAbsolute T_sat_an_in "Dew point at anode inlet";
    //
    // Cathode humidity
    parameter FCSys.Assemblies.Cells.Examples.Enumerations.HumiditySpec
      caHumiditySpec=HumiditySpec.relative "Type of specification"
      annotation (Dialog(tab="Cathode", group="Humidity"));
    Connectors.RealInput u_ca_humidity=0.5 "Value of the specification"
      annotation (Dialog(
        tab="Cathode",
        group="Humidity",
        __Dymola_label="<html><i>u</i><sub>ca humidity</sub></html>"));

    Q.NumberAbsolute caInletRH(displayUnit="%")
      "Relative humidity at cathode inlet";
    Q.PressureAbsolute p_H2O_ca_in "H2O vapor pressure at cathode inlet";
    Q.TemperatureAbsolute T_sat_ca_in "Dew point at cathode inlet";
    //
    // Anode thermal
    parameter FCSys.Assemblies.Cells.Examples.Enumerations.ThermalSpec
      anThermalSpec=ThermalSpec.temperature "Type of specification"
      annotation (Dialog(tab="Anode",group="End plate thermal"));
    Connectors.RealInput u_an_thermal=333.15*U.K "Value of the specification"
      annotation (Dialog(
        tab="Anode",
        group="End plate thermal",
        __Dymola_label="<html><i>u</i><sub>an end plate</sub></html>"));

    Q.TemperatureAbsolute T_an "Temperature of anode end plate";
    Q.ResistanceThermal Theta_an
      "Thermal resistance between the cathode end plate and the environment";
    Q.Power Qdot_an "Heat flow rate from the anode end plate";
    //
    // Cathode thermal
    parameter FCSys.Assemblies.Cells.Examples.Enumerations.ThermalSpec
      caThermalSpec=ThermalSpec.temperature "Type of specification"
      annotation (Dialog(tab="Cathode",group="End plate thermal"));
    Connectors.RealInput u_ca_thermal=333.15*U.K "Value of the specification"
      annotation (Dialog(
        tab="Cathode",
        group="End plate thermal",
        __Dymola_label="<html><i>u</i><sub>ca end plate</sub></html>"));

    Q.TemperatureAbsolute T_ca "Temperature of cathode end plate";
    Q.ResistanceThermal Theta_ca
      "Thermal resistance between the cathode end plate and the environment";
    Q.Power Qdot_ca "Heat flow rate from the cathode end plate";

    // Derived and measured conditions
    Q.PressureAbsolute p_sat_an_in "Saturation pressure at the anode inlet";
    Q.PressureAbsolute p_sat_ca_in "Saturation pressure at the cathode inlet";

    // Auxiliary measurements
    output Q.Current Ndot_H2(stateSelect=StateSelect.never) = sum(anSource.gas.H2.boundary.Ndot
       + anSink.gas.H2.boundary.Ndot) if environment.analysis
      "Net rate of hydrogen into the cell";
    output Q.Current Ndot_H2O(stateSelect=StateSelect.never) = sum(anSource.gas.H2O.boundary.Ndot
       + anSink.gas.H2O.boundary.Ndot) + sum(caSource.gas.H2O.boundary.Ndot +
      caSink.gas.H2O.boundary.Ndot) if environment.analysis
      "Net rate of water from the cell";
    output Q.Current Ndot_O2(stateSelect=StateSelect.never) = sum(caSource.gas.O2.boundary.Ndot
       + caSink.gas.O2.boundary.Ndot) if environment.analysis
      "Net rate of oxygen into the cell";
    output Q.NumberAbsolute anOutletRH(
      stateSelect=StateSelect.never,
      displayUnit="%") = average(average(anSink.gas.H2O.boundary.p ./ p_sat(
      anSink.gas.H2O.boundary.T))) if environment.analysis
      "Relative humidity at the anode outlet";
    output Q.NumberAbsolute caOutletRH(
      stateSelect=StateSelect.never,
      displayUnit="%") = average(average(caSink.gas.H2O.boundary.p ./ p_sat(
      caSink.gas.H2O.boundary.T))) if environment.analysis
      "Relative humidity at the cathode outlet";

    Connectors.BoundaryBus an[n_y, n_z] "Interface to the anode end plate"
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-70,-10}), iconTransformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-160,0})));
    Connectors.BoundaryBus anNegative[n_x_an, n_z]
      "Negative interface to the anode flow channel" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-40,-50}),iconTransformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-40,-160})));
    Connectors.BoundaryBus anPositive[n_x_an, n_z]
      "Positive interface to the anode flow channel" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-40,30}), iconTransformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-40,160})));
    Connectors.BoundaryBus ca[n_y, n_z] "Interface to the cathode end plate"
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={70,-10}),iconTransformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={162,0})));
    Connectors.BoundaryBus caNegative[n_x_ca, n_z]
      "Negative interface to the cathode flow channel" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={40,-50}),iconTransformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={40,-160})));
    Connectors.BoundaryBus caPositive[n_x_ca, n_z]
      "Positive interface to the cathode flow channel" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={40,30}), iconTransformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={40,160})));

    Conditions.ByConnector.BoundaryBus.Single.Sink anBC[n_y, n_z](each graphite(
        'incle-'=true,
        'e-'(materialSet(y=U.bar)),
        'inclC+'=true,
        redeclare
          FCSys.Conditions.ByConnector.ThermoDiffusive.Single.Temperature 'C+'(
            source(y=T_an)))) "Anode end plate" annotation (Placement(
          transformation(
          extent={{10,10},{-10,-10}},
          rotation=90,
          origin={-56,-10})));
    Conditions.ByConnector.BoundaryBus.Single.Source caBC[n_y, n_z](each
        graphite(
        'incle-'=true,
        'e-'(
          redeclare function materialSpec =
              FCSys.Conditions.ByConnector.Boundary.Single.Material.pressure,
          materialSet(y=p_ca_elec),
          redeclare function thermalSpec =
              FCSys.Conditions.ByConnector.Boundary.Single.ThermoDiffusive.temperature,

          thermalSet(y=T_ca)),
        'inclC+'=true,
        'C+'(source(y=T_ca)))) "Cathode end plate" annotation (Placement(
          transformation(
          extent={{-10,10},{10,-10}},
          rotation=270,
          origin={56,-10})));

    Conditions.ByConnector.BoundaryBus.Single.Source anSource[n_x_an, n_z](
        each gas(
        inclH2=true,
        inclH2O=true,
        H2(
          redeclare function materialSpec =
              FCSys.Conditions.ByConnector.Boundary.Single.Material.pressure,
          materialSet(y=p_an_in - p_H2O_an_in),
          thermalSet(y=T_an_in)),
        H2O(
          redeclare function materialSpec =
              FCSys.Conditions.ByConnector.Boundary.Single.Material.pressure,
          materialSet(y=p_H2O_an_in),
          thermalSet(y=T_an_in)))) "Anode source" annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-24,-30})));

    Conditions.ByConnector.BoundaryBus.Single.Source caSource[n_x_ca, n_z](
        each gas(
        inclH2O=true,
        inclN2=true,
        inclO2=true,
        H2O(
          redeclare function materialSpec =
              FCSys.Conditions.ByConnector.Boundary.Single.Material.pressure,
          materialSet(y=p_H2O_ca_in),
          thermalSet(y=T_ca_in)),
        N2(
          redeclare function materialSpec =
              FCSys.Conditions.ByConnector.Boundary.Single.Material.pressure,
          materialSet(y=(p_ca_in - p_H2O_ca_in)*(1 - psi_dry_O2_in)),
          thermalSet(y=T_ca_in)),
        O2(
          redeclare function materialSpec =
              FCSys.Conditions.ByConnector.Boundary.Single.Material.pressure,
          materialSet(y=(p_ca_in - p_H2O_ca_in)*psi_dry_O2_in),
          thermalSet(y=T_ca_in)))) "Cathode source" annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={24,-30})));

    Conditions.ByConnector.BoundaryBus.Single.Sink anSink[n_x_an, n_z](gas(
        each inclH2=true,
        each inclH2O=true,
        H2O(materialSet(y=fill(
                  p_an_out,
                  n_x_an,
                  n_z) - anSink.gas.H2.p), redeclare each function materialSpec
            = FCSys.Conditions.ByConnector.Boundary.Single.Material.volumeRate),

        H2(materialSet(y=anSink.gas.H2O.materialOut.y), redeclare each function
            materialSpec =
              FCSys.Conditions.ByConnector.Boundary.Single.Material.volumeRate)),
        each liquid(H2O(materialSet(y=p_an_out)))) "Anode sink" annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-24,10})));

    Conditions.ByConnector.BoundaryBus.Single.Sink caSink[n_x_ca, n_z](gas(
        each inclO2=true,
        each inclN2=true,
        each inclH2O=true,
        H2O(materialSet(y=fill(
                  p_ca_out,
                  n_x_ca,
                  n_z) - caSink.gas.N2.p - caSink.gas.O2.p)),
        N2(materialSet(y=caSink.gas.H2O.materialOut.y), redeclare each function
            materialSpec =
              FCSys.Conditions.ByConnector.Boundary.Single.Material.volumeRate),

        O2(materialSet(y=caSink.gas.H2O.materialOut.y), redeclare each function
            materialSpec =
              FCSys.Conditions.ByConnector.Boundary.Single.Material.volumeRate)),
        each liquid(H2O(materialSet(y=p_ca_out)))) "Cathode sink" annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={24,10})));

  protected
    Conditions.Router anRouter[n_x_an, n_z](each final crossOver=anInletSide
           == Side.p) annotation (Dialog, Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-32,-10})));
    Conditions.Router caRouter[n_x_an, n_z](each final crossOver=caInletSide
           == Side.p) annotation (Dialog, Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={32,-10})));
    outer Conditions.Environment environment "Environmental conditions";

    Q.Pressure p_ca_elec "Electronic pressure on the cathode";

  equation
    // Electrical
    // ----------
    // Measurements
    zI = -sum(caBC.graphite.'e-'.boundary.Ndot);
    w = average(average('Datae-'.g(caBC.graphite.'e-'.boundary.T, caBC.graphite.
      'e-'.boundary.p) - 'Datae-'.g(anBC.graphite.'e-'.boundary.T, anBC.graphite.
      'e-'.boundary.p)));
    R*zI = w;
    P = w*zI;
    //
    // Setpoint
    if electricalSpec == ElectricalSpec.current then
      zI = u_electrical;
    elseif electricalSpec == ElectricalSpec.voltage then
      w = u_electrical;
    elseif electricalSpec == ElectricalSpec.resistance then
      R = u_electrical;
    else
      // electricalSpec == ElectricalSpec.power
      P = u_electrical;
    end if;

    // Anode flow rate
    // ---------------
    // Measurements
    // p_an_in is set in anSource.
    anStoich*zI = I_an;
    I_an = -sum(anSource.gas.H2.boundary.Ndot)/2;
    //
    // Setpoint
    if anFlowSpec == FlowSpec.stoich then
      anStoich = u_an_flow;
    elseif anFlowSpec == FlowSpec.current then
      I_an = u_an_flow;
    else
      // anFlowSpec == FlowSpec.pressure
      p_an_in = u_an_flow;
    end if;

    // Cathode flow rate
    // -----------------
    // Measurements
    caStoich*zI = I_ca;
    I_ca = -sum(caSource.gas.O2.boundary.Ndot)/4;
    //
    // Setpoint
    if caFlowSpec == FlowSpec.stoich then
      caStoich = u_ca_flow;
    elseif caFlowSpec == FlowSpec.current then
      I_ca = u_ca_flow;
    else
      // caFlowSpec == FlowSpec.pressure
      p_ca_in = u_ca_flow;
    end if;

    // Anode humidity
    // --------------
    // Measurements
    p_sat_an_in = p_sat(T_an_in);
    p_H2O_an_in = p_sat(T_sat_an_in);
    p_H2O_an_in = min(anInletRH, 1)*p_sat_an_in;
    //
    // Setpoint
    if anHumiditySpec == HumiditySpec.relative then
      anInletRH = u_an_humidity;
    elseif anHumiditySpec == HumiditySpec.pressure then
      p_H2O_an_in = u_an_humidity;
    else
      // anHumiditySpec == HumiditySpec.dewPoint
      T_sat_an_in = u_an_humidity;
    end if;

    // Cathode humidity
    // ----------------
    // Measurements
    p_sat_ca_in = p_sat(T_ca_in);
    p_H2O_ca_in = p_sat(T_sat_ca_in);
    p_H2O_ca_in = min(caInletRH, 1)*p_sat_ca_in;
    //
    // Setpoint
    if caHumiditySpec == HumiditySpec.relative then
      caInletRH = u_ca_humidity;
    elseif caHumiditySpec == HumiditySpec.pressure then
      p_H2O_ca_in = u_ca_humidity;
    else
      // caHumiditySpec == HumiditySpec.dewPoint
      T_sat_ca_in = u_ca_humidity;
    end if;

    // Anode thermal
    // -------------
    // Measurements
    // T_an is set in anBC.
    Theta_an*Qdot_an = T_an - environment.T;
    Qdot_an = sum(anBC.graphite.'C+'.thermo.Qdot + anBC.graphite.'e-'.boundary.Qdot);
    //
    // Setpoint
    if anThermalSpec == ThermalSpec.temperature then
      T_an = u_an_thermal;
    elseif anThermalSpec == ThermalSpec.resistance then
      Theta_an = u_an_thermal;
    else
      // anThermalSpec == ThermalSpec.heatRate
      Qdot_an = u_an_thermal;
    end if;

    // Cathode thermal
    // ---------------
    // Measurements
    // T_ca is set in caBC.
    Theta_ca*Qdot_ca = T_ca - environment.T;
    Qdot_ca = sum(caBC.graphite.'C+'.thermo.Qdot + caBC.graphite.'e-'.boundary.Qdot);
    //
    // Setpoint
    if caThermalSpec == ThermalSpec.temperature then
      T_ca = u_ca_thermal;
    elseif caThermalSpec == ThermalSpec.resistance then
      Theta_ca = u_ca_thermal;
    else
      // caThermalSpec == ThermalSpec.heatRate
      Qdot_ca = u_ca_thermal;
    end if;

    connect(caBC.boundary, ca) annotation (Line(
        points={{60,-10},{70,-10}},
        color={127,127,127},
        thickness=0.5,
        smooth=Smooth.None));

    connect(anRouter.positive2, anSink.boundary) annotation (Line(
        points={{-24,-6},{-24,6}},
        color={127,127,127},
        thickness=0.5,
        smooth=Smooth.None));
    connect(anRouter.positive1, anSource.boundary) annotation (Line(
        points={{-24,-14},{-24,-26}},
        color={127,127,127},
        thickness=0.5,
        smooth=Smooth.None));
    connect(anPositive, anRouter.negative2) annotation (Line(
        points={{-40,30},{-40,-6}},
        color={127,127,127},
        thickness=0.5,
        smooth=Smooth.None));
    connect(anRouter.negative1, anNegative) annotation (Line(
        points={{-40,-14},{-40,-50}},
        color={127,127,127},
        thickness=0.5,
        smooth=Smooth.None));
    connect(caSink.boundary, caRouter.negative2) annotation (Line(
        points={{24,6},{24,-6}},
        color={127,127,127},
        thickness=0.5,
        smooth=Smooth.None));
    connect(caSource.boundary, caRouter.negative1) annotation (Line(
        points={{24,-26},{24,-14}},
        color={127,127,127},
        thickness=0.5,
        smooth=Smooth.None));
    connect(caRouter.positive2, caPositive) annotation (Line(
        points={{40,-6},{40,30}},
        color={127,127,127},
        thickness=0.5,
        smooth=Smooth.None));
    connect(caRouter.positive1, caNegative) annotation (Line(
        points={{40,-14},{40,-50}},
        color={127,127,127},
        thickness=0.5,
        smooth=Smooth.None));
    connect(anBC.boundary, an) annotation (Line(
        points={{-60,-10},{-70,-10}},
        color={127,127,127},
        thickness=0.5,
        smooth=Smooth.None));
    annotation (
      Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-80,-60},{80,
              40}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-160,-160},{160,
              160}}), graphics={Rectangle(
              extent={{-160,160},{160,-160}},
              lineColor={191,191,191},
              fillColor={255,255,255},
              fillPattern=FillPattern.Backward),Rectangle(extent={{-160,160},{
            160,-160}}, lineColor={0,0,0})}),
      Documentation(info="
    <html>
    <p>Any of the settings for the operating conditions can be time-varying expressions.</p>
    
    <p><i>Equivalent current</i> is the rate of supply of a reactant required to support the
    given current
    assuming the reactant is entirely consumed (complete utilization).</p>

    <p>Assumptions:
    <ol>
    <li>The outer yz surface of each end plate is each uniform in temperature.</li>
    <li>No heat is conducted from the rest of the cell hardware.</li>
    <li>The electronic pressure is uniform across each end plate.</li>
    <li>There is no shear force on the fluid at either outlet.</li>
    <li>The pressure of each gas is uniform over the inlets.</li>
    <li>The volumetric flow rate of the gases is equal at each outlet.</li>
    <li>The outlet pressure is applied to the gas mixture by Dalton's law (additivity of pressure).</li>
    <li>At the outlet, the liquid has the same pressure as the gas (Amagat's law).</li>
    <li>There is no thermal conduction across either outlet.</li>
    </ol></p>
    </html>"));
  end ConditionsExamplesTestStandsTestStand;

  model ConditionsTestStandsTestStandEIS
    "Test stand to perform electrochemical impedance spectroscopy"
    extends WorkInProgress.ConditionsExamplesTestStandsTestStand(final
        electricalSpec=ElectricalSpec.current, final u_electrical=zI_large +
          zI_small_A*U.A);

    parameter Q.Current zI_large=U.A "Large-signal current" annotation (Dialog(
          __Dymola_label="<html><i>zJ</i><sub>large</sub></html>"));
    Connectors.RealInput zI_small_A "Small-signal current in amperes"
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=-45,
          origin={-107,107}), iconTransformation(
          extent={{-10,-10},{10,10}},
          rotation=-45,
          origin={-167,167})));
    Connectors.RealOutput w_V=w/U.V "Cell potential in volts" annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=-45,
          origin={107,-107}), iconTransformation(
          extent={{-10,-10},{10,10}},
          rotation=-45,
          origin={167,-167})));

    annotation (
      Documentation(info="<html><p>This model modulates the electrical current applied to the cell 
    according to an input.
    The current is the sum of a steady-state large-signal current and a small-signal 
    current introduced via the input <i>zI</i><sub>small A</sub>.</p>
       
    <p>For more information, please see the documentation in the
    <a href=\"modelica://FCSys.Conditions.TestStands.TestStand\">test stand</a> model.</p></html>"),

      Diagram(coordinateSystem(preserveAspectRatio=true,extent={{-100,-100},{
              100,100}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-160,-160},{160,
              160}}), graphics));

  end ConditionsTestStandsTestStandEIS;

  model SubregionsExamplesCapacitor
    "Test the storage of charge on the electrolytic double layer"
    extends Subregions.Examples.Subregion(
      subregion(ionomer('H+'(N(stateSelect=StateSelect.default)))),
      inclH2=true,
      'inclC+'=true,
      'incle-'=true,
      'inclSO3-'=true,
      'inclH+'=true);

    output Q.Potential Sigmag=subregion.graphite.'e-'.g + subregion.ionomer.
        'H+'.g "Sum of the electronic and protonic potentials";

    Conditions.ByConnector.BoundaryBus.Single.Source electrons(graphite(
          'incle-'=true, 'e-'(redeclare Modelica.Blocks.Sources.Ramp
            materialSet(height=-U.A, duration=1)))) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-24,0})));
    Conditions.ByConnector.BoundaryBus.Single.Source protons(ionomer('inclH+'=
            true,'H+'(redeclare function materialSpec =
              FCSys.Conditions.ByConnector.Boundary.Single.Material.pressure,
            materialSet(y=U.bar)))) annotation (Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=90,
          origin={24,0})));

  equation
    connect(protons.boundary, subregion.xPositive) annotation (Line(
        points={{20,0},{10,0}},
        color={127,127,127},
        thickness=0.5,
        smooth=Smooth.None));
    connect(electrons.boundary, subregion.xNegative) annotation (Line(
        points={{-20,0},{-10,0}},
        color={127,127,127},
        thickness=0.5,
        smooth=Smooth.None));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics));
  end SubregionsExamplesCapacitor;

  model AssembliesCellsExamplesTestStandEIS
    "Test stand to perform electrochemical impedance spectroscopy"
    extends Assemblies.Cells.Examples.TestStand(redeclare
        Assemblies.Cells.Examples.TestConditions testConditions(final
          electricalSpec=ElectricalSpec.current, electricalSet(y=zI_large +
              firstOrder.y*U.A)));

    parameter Q.Current zI_large=50*U.A "Large-signal current" annotation (
        Dialog(__Dymola_label="<html><i>zJ</i><sub>large</sub></html>"));
    Connectors.RealInput zI_small_A "Small-signal current in amperes"
      annotation (Placement(transformation(extent={{-24,-66},{-4,-46}})));
    Connectors.RealOutput w_V=w/U.V "Cell potential in volts";

    Modelica.Blocks.Continuous.FirstOrder firstOrder(
      initType=Modelica.Blocks.Types.Init.InitialOutput,
      y_start=0,
      T=1e-8)
      annotation (Placement(transformation(extent={{-24,-66},{-4,-46}})));
  equation
    connect(firstOrder.u, zI_small_A) annotation (Line(
        points={{-26,-56},{-20,-56},{-20,-56},{-14,-56}},
        color={0,0,127},
        smooth=Smooth.None));
    annotation (Documentation(info="<html><p>This model modulates the electrical current applied to the cell 
    according to an input.
    The current is the sum of a steady-state large-signal current and a small-signal 
    current introduced via the input <i>zI</i><sub>small A</sub>.</p>
       
    <p>For more information, please see the documentation in the
    <a href=\"modelica://FCSys.Conditions.TestStands.TestStand\">test stand</a> model.</p></html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-80,-60},{
              80,40}}), graphics));
  end AssembliesCellsExamplesTestStandEIS;
  annotation (Commands(
      file="../../units.mos"
        "Establish the constants and units in the workspace (first translate a model besides Units.Evaluate).",

      file="test/check.mos" "Check all of FCSys using Dymola's check function.",

      file="../../../LaTeX/Dissertation/Results/Cell/Simulation/sim.mos"));

end WorkInProgress;
