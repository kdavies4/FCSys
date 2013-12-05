within FCSys;
package WorkInProgress "Incomplete classes under development"
  extends Modelica.Icons.Package;

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

  model CellMSL
    "<html>Single-cell PEMFC with interfaces from the <a href=\"modelica://Modelica\">Modelica Standard Library</a></html>"
    extends FCSys.Icons.Cell;

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
    annotation (defaultComponentName="cell",experiment(StopTime=1000));
  end CellMSL;

  function plot "Create plots using FCRes"
    extends Modelica.Icons.Function;
    extends Modelica.Icons.UnderConstruction;

  algorithm
    Modelica.Utilities.System.command("loadres");

  end plot;

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
        redeclare Conditions.ByConnector.ThermalDiffusive.Single.Temperature
          'C+'(set(y=T_an)))) "Anode end plate" annotation (Placement(
          transformation(
          extent={{10,10},{-10,-10}},
          rotation=90,
          origin={-56,-10})));
    Conditions.ByConnector.BoundaryBus.Single.Source caBC[n_y, n_z](each
        graphite(
        'incle-'=true,
        'e-'(
          redeclare function materialSpec =
              Conditions.ByConnector.Boundary.Single.Material.pressure,
          materialSet(y=p_ca_elec),
          redeclare function thermalSpec =
              Conditions.ByConnector.Boundary.Single.Thermal.temperature,
          thermalSet(y=T_ca)),
        'inclC+'=true,
        'C+'(set(y=T_ca)))) "Cathode end plate" annotation (Placement(
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
              Conditions.ByConnector.Boundary.Single.Material.pressure,
          materialSet(y=p_an_in - p_H2O_an_in),
          thermalSet(y=T_an_in)),
        H2O(
          redeclare function materialSpec =
              Conditions.ByConnector.Boundary.Single.Material.pressure,
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
              Conditions.ByConnector.Boundary.Single.Material.pressure,
          materialSet(y=p_H2O_ca_in),
          thermalSet(y=T_ca_in)),
        N2(
          redeclare function materialSpec =
              Conditions.ByConnector.Boundary.Single.Material.pressure,
          materialSet(y=(p_ca_in - p_H2O_ca_in)*(1 - psi_dry_O2_in)),
          thermalSet(y=T_ca_in)),
        O2(
          redeclare function materialSpec =
              Conditions.ByConnector.Boundary.Single.Material.pressure,
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
              Conditions.ByConnector.Boundary.Single.Material.volumeRate)),
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
              Conditions.ByConnector.Boundary.Single.Material.volumeRate),
        O2(materialSet(y=caSink.gas.H2O.materialOut.y), redeclare each function
            materialSpec =
              Conditions.ByConnector.Boundary.Single.Material.volumeRate)),
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
    Qdot_an = sum(anBC.graphite.'C+'.therm.Qdot + anBC.graphite.'e-'.boundary.Qdot);
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
    Qdot_ca = sum(caBC.graphite.'C+'.therm.Qdot + caBC.graphite.'e-'.boundary.Qdot);
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

  model AssembliesCellsExamplesTestStandEIS
    "Test stand to perform electrochemical impedance spectroscopy"
    extends Assemblies.Cells.Examples.TestStand2(redeclare
        Assemblies.Cells.Examples.TestConditions testConditions(final
          electricalSpec=ElectricalSpec.current, electricalSet(y=zI_large +
              firstOrder.y*U.A)));
    Connectors.RealOutput w_V=w/U.V "Cell potential in volts";
    parameter Q.Current zI_large=50*U.A "Large-signal current" annotation (
        Dialog(__Dymola_label="<html><i>zJ</i><sub>large</sub></html>"));

    Connectors.RealInput zI_small_A "Small-signal current in amperes"
      annotation (Placement(transformation(extent={{-24,-66},{-4,-46}})));

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

  model AssembliesCellsExamplesTestStandO2
    "Simulate the fuel cell with prescribed conditions, with pure oxygen"
    extends Assemblies.Cells.Examples.TestStand2(
      p_ca_in=sum(caSource.gas.H2O.boundary.p + caSource.gas.O2.boundary.p),
      cell(final inclN2=false),
      testConditions(final psi_O2_dry=1),
      caSource(gas(each final inclN2=false)),
      caSink(gas(H2O(materialSet(y=fill(
                    testConditions.p,
                    cell.anFP.n_x,
                    cell.n_z) - caSink.gas.O2.p)), each final inclN2=false)));

    extends Modelica.Icons.UnderConstruction;

    annotation (experiment(StopTime=700, Tolerance=1e-005),
        __Dymola_experimentSetupOutput);

  end AssembliesCellsExamplesTestStandO2;

  model RegionsExamplesCLtoCLVoltage
    "Test one catalyst layer to the other, with prescribed voltage"

    extends Modelica.Icons.Example;
    extends Modelica.Icons.UnderConstruction;
    output Q.Potential w=anCL.subregions[1, 1, 1].graphite.'e-'.g_boundaries[1,
        Side.n] - caCL.subregions[1, 1, 1].graphite.'e-'.g_boundaries[1, Side.p]
      if environment.analysis "Electrical potential";
    output Q.Current zI=-sum(anCL.subregions[1, :, :].graphite.'e-'.boundaries[
        1, Side.n].Ndot) if environment.analysis "Electrical current";
    output Q.Number J_Apercm2=zI*U.cm^2/(anCL.A[Axis.x]*U.A)
      "Electrical current density, in A/cm2";

    parameter Q.Length L_y[:]={8}*U.cm "**Lengths in the y direction";
    parameter Q.Length L_z[:]={6.25}*U.cm "**Lengths in the z direction";
    Regions.AnCLs.AnCL anCL(
      final L_y=L_y,
      final L_z=L_z,
      subregions(graphite(each inclDL=true, 'e-Transfer'(each fromI=false))))
      annotation (Placement(transformation(extent={{-30,30},{-10,50}})));

    Regions.PEMs.PEM PEM(
      final L_y=L_y,
      final L_z=L_z,
      subregions(ionomer('H+'(each consTransX=ConsTrans.dynamic))))
      annotation (Placement(transformation(extent={{-10,30},{10,50}})));
    Regions.CaCLs.CaCL caCL(
      final L_y=L_y,
      final L_z=L_z,
      subregions(graphite(each inclDL=true, 'e-Transfer'(each fromI=false)),
          each ORR('e-'(reaction(Ndot(stateSelect=StateSelect.always))))))
      annotation (Placement(transformation(extent={{10,30},{30,50}})));

    // Conditions
    Conditions.ByConnector.BoundaryBus.Single.Sink anBC[anCL.n_y, anCL.n_z](
        each gas(
        inclH2=true,
        inclH2O=true,
        H2(
          materialSet(y=environment.p_dry),
          redeclare function thermalSpec =
              Conditions.ByConnector.Boundary.Single.Thermal.temperature,
          thermalSet(y=environment.T)),
        H2O(
          redeclare function materialSpec =
              Conditions.ByConnector.Boundary.Single.Material.pressure,
          materialSet(y=environment.p_H2O),
          redeclare function thermalSpec =
              Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
          thermalSet(y=0))), each graphite(
        'inclC+'=true,
        'incle-'=true,
        redeclare Conditions.ByConnector.ThermalDiffusive.Single.Temperature
          'C+'(set(y=environment.T)),
        'e-'(
          redeclare function thermalSpec =
              Conditions.ByConnector.Boundary.Single.Thermal.temperature,
          thermalSet(y=environment.T),
          redeclare function materialMeas =
              Conditions.ByConnector.Boundary.Single.Material.potential (
                redeclare package Data = FCSys.Characteristics.'e-'.Graphite))))
      annotation (Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=270,
          origin={-44,40})));

    Conditions.ByConnector.BoundaryBus.Single.Source caBC[caCL.n_y, caCL.n_z](
        each gas(
        inclH2O=true,
        inclO2=true,
        H2O(
          redeclare function materialSpec =
              Conditions.ByConnector.Boundary.Single.Material.pressure,
          materialSet(y=environment.p_H2O),
          redeclare function thermalSpec =
              Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
          thermalSet(y=0)),
        O2(
          redeclare function materialSpec =
              Conditions.ByConnector.Boundary.Single.Material.pressure,
          materialSet(y=environment.p_O2),
          thermalSet(y=environment.T))), graphite(
        each 'inclC+'=true,
        each 'incle-'=true,
        each 'C+'(set(y=environment.T)),
        'e-'(
          redeclare each function materialSpec =
              Conditions.ByConnector.Boundary.Single.Material.potential (
                redeclare package Data = FCSys.Characteristics.'e-'.Graphite),
          materialSet(y=anBC.graphite.'e-'.materialOut.y + fill(
                  -voltageSet.y,
                  caCL.n_y,
                  caCL.n_z)),
          each thermalSet(y=environment.T)))) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={44,40})));

    Modelica.Blocks.Sources.Ramp voltageSet(
      duration=300,
      offset=1.19997*U.V,
      height=-0.8*U.V)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    inner Conditions.Environment environment(
      analysis=true,
      T=333.15*U.K,
      p=U.from_kPag(48.3),
      RH=0.7) "Environmental conditions"
      annotation (Placement(transformation(extent={{-10,70},{10,90}})));
  equation
    connect(anCL.xPositive, PEM.xNegative) annotation (Line(
        points={{-10,40},{-10,40}},
        color={240,0,0},
        thickness=0.5,
        smooth=Smooth.None));
    connect(PEM.xPositive, caCL.xNegative) annotation (Line(
        points={{10,40},{10,40}},
        color={0,0,240},
        thickness=0.5,
        smooth=Smooth.None));
    connect(anBC.boundary, anCL.xNegative) annotation (Line(
        points={{-40,40},{-30,40}},
        color={240,0,0},
        thickness=0.5,
        smooth=Smooth.None));
    connect(caCL.xPositive, caBC.boundary) annotation (Line(
        points={{30,40},{40,40}},
        color={0,0,240},
        thickness=0.5,
        smooth=Smooth.None));
    annotation (
      Commands(file=
            "Resources/Scripts/Dymola/Regions.Examples.CLtoCLVoltage.mos"
          "Regions.Examples.CLtoCLVoltage.mos"),
      experiment(
        StopTime=600,
        Tolerance=1e-007,
        __Dymola_Algorithm="Dassl"),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics));
  end RegionsExamplesCLtoCLVoltage;

  model RegionsExamplesFPtoFPVoltage
    "Test one flow plate to the other, with prescribed voltage"
    extends Modelica.Icons.Example;
    extends Modelica.Icons.UnderConstruction;
    parameter Boolean inclLiq=false "Include liquid H2O";
    parameter Q.NumberAbsolute psi_H2O=environment.psi_H2O
      "Mole fraction of H2O at the inlet";
    parameter Q.NumberAbsolute psi_H2=environment.psi_dry
      "Mole fraction of H2 at the inlet";
    parameter Q.NumberAbsolute psi_O2=environment.psi_O2_dry*environment.psi_dry
      "Mole fraction of O2 at the inlet";
    parameter Q.NumberAbsolute psi_N2=(1 - environment.psi_O2_dry)*environment.psi_dry
      "Mole fraction of N2 at the inlet";
    output Q.Number J_Apercm2=zI*U.cm^2/(caFP.A[Axis.x]*U.A)
      "Electrical current density, in A/cm2";
    output Q.Potential w=anFP.subregions[1, 1, 1].graphite.'e-'.g_boundaries[1,
        Side.n] - caFP.subregions[end, 1, 1].graphite.'e-'.g_boundaries[1, Side.p]
      if environment.analysis "Electrical potential";

    parameter Q.Length L_y[:]={4,4}*U.cm "**Lengths in the y direction";
    parameter Q.Length L_z[:]={6.25}*U.cm "**Lengths in the z direction";

    // Layers
    Regions.AnFPs.AnFP anFP(
      final L_y=L_y,
      final L_z=L_z,
      subregions(each liquid(inclH2O=inclLiq)))
      annotation (Placement(transformation(extent={{-70,30},{-50,50}})));
    Regions.AnGDLs.AnGDL anGDL(
      final L_y=L_y,
      final L_z=L_z,
      subregions(each liquid(inclH2O=inclLiq)))
      annotation (Placement(transformation(extent={{-50,30},{-30,50}})));
    Regions.AnCLs.AnCL anCL(
      final L_y=L_y,
      final L_z=L_z,
      subregions(graphite(each inclDL=false, 'e-Transfer'(each fromI=true)),
          each liquid(inclH2O=inclLiq)))
      annotation (Placement(transformation(extent={{-30,30},{-10,50}})));
    Regions.PEMs.PEM PEM(
      final L_y=L_y,
      final L_z=L_z,
      subregions(ionomer('H+'(each consTransX=ConsTrans.dynamic))))
      annotation (Placement(transformation(extent={{-10,30},{10,50}})));
    Regions.CaCLs.CaCL caCL(
      final L_y=L_y,
      final L_z=L_z,
      subregions(
        graphite(each inclDL=false, 'e-Transfer'(each fromI=true)),
        each liquid(inclH2O=inclLiq),
        each ORR('e-'(reaction(Ndot(stateSelect=StateSelect.default))))))
      annotation (Placement(transformation(extent={{10,30},{30,50}})));
    Regions.CaGDLs.CaGDL caGDL(
      final L_y=L_y,
      final L_z=L_z,
      subregions(each liquid(inclH2O=inclLiq)))
      annotation (Placement(transformation(extent={{30,30},{50,50}})));
    Regions.CaFPs.CaFP caFP(
      final L_y=L_y,
      final L_z=L_z,
      subregions(graphite('e-'(each consTransY=ConsTrans.steady)), each liquid(
            inclH2O=inclLiq)))
      annotation (Placement(transformation(extent={{50,30},{70,50}})));

    // Conditions
    Conditions.ByConnector.BoundaryBus.Single.Sink anBC1[1, anFP.n_z](each gas(
        inclH2=true,
        inclH2O=true,
        H2(
          redeclare function materialSpec =
              Conditions.ByConnector.Boundary.Single.Material.current,
          materialSet(y=0),
          redeclare function afterSpec =
              Conditions.ByConnector.Boundary.Single.Translational.velocity,
          redeclare function beforeSpec =
              Conditions.ByConnector.Boundary.Single.Translational.velocity),
        H2O(
          redeclare function materialSpec =
              Conditions.ByConnector.Boundary.Single.Material.current,
          materialSet(y=0),
          redeclare function afterSpec =
              Conditions.ByConnector.Boundary.Single.Translational.velocity,
          redeclare function beforeSpec =
              Conditions.ByConnector.Boundary.Single.Translational.velocity)),
        each graphite(
        'inclC+'=true,
        'incle-'=true,
        redeclare Conditions.ByConnector.ThermalDiffusive.Single.Temperature
          'C+'(set(y=environment.T)),
        'e-'(
          redeclare function thermalSpec =
              Conditions.ByConnector.Boundary.Single.Thermal.temperature,
          thermalSet(y=environment.T),
          redeclare function materialMeas =
              Conditions.ByConnector.Boundary.Single.Material.potential (
                redeclare package Data = FCSys.Characteristics.'e-'.Graphite))))
      annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=90,
          origin={-84,54})));

    Conditions.ByConnector.BoundaryBus.Single.Sink anBC2[anFP.n_y - 1, anFP.n_z]
      (each gas(
        inclH2=true,
        inclH2O=true,
        H2(
          redeclare function materialSpec =
              Conditions.ByConnector.Boundary.Single.Material.current,
          materialSet(y=0),
          redeclare function afterSpec =
              Conditions.ByConnector.Boundary.Single.Translational.velocity,
          redeclare function beforeSpec =
              Conditions.ByConnector.Boundary.Single.Translational.velocity),
        H2O(
          redeclare function materialSpec =
              Conditions.ByConnector.Boundary.Single.Material.current,
          materialSet(y=0),
          redeclare function afterSpec =
              Conditions.ByConnector.Boundary.Single.Translational.velocity,
          redeclare function beforeSpec =
              Conditions.ByConnector.Boundary.Single.Translational.velocity)),
        each graphite(
        'inclC+'=true,
        'incle-'=true,
        redeclare Conditions.ByConnector.ThermalDiffusive.Single.Temperature
          'C+'(set(y=environment.T)),
        'e-'(
          redeclare function materialSpec =
              Conditions.ByConnector.Boundary.Single.Material.current,
          materialSet(y=0),
          redeclare function thermalSpec =
              Conditions.ByConnector.Boundary.Single.Thermal.temperature,
          thermalSet(y=environment.T)))) if anFP.n_y > 1 annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=90,
          origin={-84,26})));

    Conditions.ByConnector.BoundaryBus.Single.Source anSource[anFP.n_x, anFP.n_z]
      (each gas(
        inclH2=true,
        inclH2O=true,
        H2(materialSet(y=-Ndot_H2), thermalSet(y=environment.T)),
        H2O(materialSet(y=-Ndot_H2O_an), thermalSet(y=environment.T))))
      annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-60,16})));
    Conditions.ByConnector.BoundaryBus.Single.Sink anSink[anFP.n_x, anFP.n_z](
        gas(
        each inclH2=true,
        each inclH2O=true,
        H2(materialSet(y=anSink.gas.H2O.boundary.Ndot .* anFP.subregions[:,
                anFP.n_y, :].gas.H2O.v ./ anFP.subregions[:, anFP.n_y, :].gas.H2.v),
            redeclare each function materialSpec =
              Conditions.ByConnector.Boundary.Single.Material.current),
        H2O(materialSet(y=fill(
                  environment.p,
                  anFP.n_x,
                  anFP.n_z) - anSink.gas.H2.p))), each liquid(H2O(materialSet(y
              =environment.p)), inclH2O=inclLiq)) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-60,64})));

    Conditions.ByConnector.BoundaryBus.Single.Source caBC1[1, caFP.n_z](each
        gas(
        inclH2O=true,
        inclN2=true,
        inclO2=true,
        H2O(redeclare function thermalSpec =
              Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
            thermalSet(y=0)),
        N2(redeclare function thermalSpec =
              Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
            thermalSet(y=0)),
        O2(redeclare function thermalSpec =
              Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
            thermalSet(y=0))), graphite(
        each 'inclC+'=true,
        each 'incle-'=true,
        each 'C+'(set(y=environment.T)),
        'e-'(
          redeclare each function materialSpec =
              Conditions.ByConnector.Boundary.Single.Material.potential (
                redeclare package Data = FCSys.Characteristics.'e-'.Graphite),
          materialSet(y=anBC1.graphite.'e-'.materialOut.y + fill(
                  -voltageSet.y,
                  1,
                  caCL.n_z)),
          each thermalSet(y=environment.T)))) annotation (Placement(
          transformation(
          extent={{10,10},{-10,-10}},
          rotation=90,
          origin={84,54})));

    Conditions.ByConnector.BoundaryBus.Single.Source caBC2[caFP.n_y - 1, caFP.n_z]
      (each gas(
        inclH2O=true,
        inclN2=true,
        inclO2=true,
        H2O(redeclare function thermalSpec =
              Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
            thermalSet(y=0)),
        N2(redeclare function thermalSpec =
              Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
            thermalSet(y=0)),
        O2(redeclare function thermalSpec =
              Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
            thermalSet(y=0))),each graphite(
        'inclC+'=true,
        'incle-'=true,
        'C+'(set(y=environment.T)),
        'e-'(
          redeclare function materialSpec =
              Conditions.ByConnector.Boundary.Single.Material.current,
          materialSet(y=0),
          thermalSet(y=environment.T)))) if caFP.n_y > 1 annotation (Placement(
          transformation(
          extent={{10,10},{-10,-10}},
          rotation=90,
          origin={84,26})));

    Conditions.ByConnector.BoundaryBus.Single.Source caSource[caFP.n_x, caFP.n_z]
      (each gas(
        inclH2O=true,
        inclN2=true,
        inclO2=true,
        H2O(materialSet(y=-Ndot_H2O_ca), thermalSet(y=environment.T)),
        N2(materialSet(y=-Ndot_N2), thermalSet(y=environment.T)),
        O2(materialSet(y=-Ndot_O2), thermalSet(y=environment.T)))) annotation (
        Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={60,16})));

    Conditions.ByConnector.BoundaryBus.Single.Sink caSink[caFP.n_x, caFP.n_z](
        gas(
        each inclH2O=true,
        each inclN2=true,
        each inclO2=true,
        H2O(materialSet(y=fill(
                  environment.p,
                  caFP.n_x,
                  caFP.n_z) - caSink.gas.N2.p - caSink.gas.O2.p)),
        N2(redeclare function materialSpec =
              Conditions.ByConnector.Boundary.Single.Material.current,
            materialSet(y=caSink.gas.H2O.boundary.Ndot .* caFP.subregions[:,
                caFP.n_y, :].gas.H2O.v ./ caFP.subregions[:, caFP.n_y, :].gas.N2.v)),

        O2(redeclare function materialSpec =
              Conditions.ByConnector.Boundary.Single.Material.current,
            materialSet(y=caSink.gas.H2O.boundary.Ndot .* caFP.subregions[:,
                caFP.n_y, :].gas.H2O.v ./ caFP.subregions[:, caFP.n_y, :].gas.O2.v))),
        each liquid(H2O(materialSet(y=environment.p)), inclH2O=inclLiq))
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={60,64})));

    Modelica.Blocks.Sources.Ramp currentSet(
      offset=U.mA,
      height=100*U.A,
      duration=600,
      startTime=50)
      annotation (Placement(transformation(extent={{-80,-30},{-60,-10}})));
    Modelica.Blocks.Math.Gain stoichH2(k=1.5/2)
      annotation (Placement(transformation(extent={{-30,-30},{-10,-10}})));
    Modelica.Blocks.Math.Gain anStoichH2O(k=psi_H2O/psi_H2)
      annotation (Placement(transformation(extent={{30,-30},{50,-10}})));
    Modelica.Blocks.Math.Gain stoichO2(k=2/4)
      annotation (Placement(transformation(extent={{-30,-90},{-10,-70}})));
    Modelica.Blocks.Math.Gain caStoichH2O(k=psi_H2O/psi_O2)
      annotation (Placement(transformation(extent={{30,-70},{50,-50}})));
    Modelica.Blocks.Math.Gain stoichN2(k=psi_N2/psi_O2)
      annotation (Placement(transformation(extent={{30,-110},{50,-90}})));

    Modelica.Blocks.Sources.Ramp voltageSet(
      duration=600,
      offset=1.19838*U.V,
      startTime=50,
      height=-0.87*U.V)
      annotation (Placement(transformation(extent={{-80,-70},{-60,-50}})));

  protected
    Connectors.RealOutputInternal Ndot_H2O_an(unit="N/T")
      "Rate of supply of H2O into the anode"
      annotation (Placement(transformation(extent={{54,-30},{74,-10}})));
    Connectors.RealOutputInternal Ndot_O2(unit="N/T") "Rate of supply of O2"
      annotation (Placement(transformation(extent={{-6,-90},{14,-70}})));
    Connectors.RealOutputInternal Ndot_H2O_ca(unit="N/T")
      "Rate of supply of H2O into the cathode"
      annotation (Placement(transformation(extent={{54,-70},{74,-50}})));
    Connectors.RealOutputInternal Ndot_N2(unit="N/T") "Rate of supply of N2"
      annotation (Placement(transformation(extent={{54,-110},{74,-90}})));
    Connectors.RealOutputInternal zI(unit="N/T") "Electrical current"
      annotation (Placement(transformation(extent={{-54,-30},{-34,-10}})));
    Connectors.RealOutputInternal Ndot_H2(unit="N/T") "Rate of supply of H2"
      annotation (Placement(transformation(extent={{-6,-30},{14,-10}})));

    inner Conditions.Environment environment(
      a={0,0,0},
      analysis=true,
      T=333.15*U.K,
      p=U.from_kPag(48.3),
      RH=0.7) "Environmental conditions"
      annotation (Placement(transformation(extent={{-10,70},{10,90}})));

    /*
  Real states[:](each stateSelect=StateSelect.always) = {caFP.subregions[1, 1, 1].gas.O2.boundaries[
    2, 2].Ndot,caFP.subregions[1, 1, 1].gas.H2O.boundaries[2, 2].Ndot,anFP.subregions[
    1, 1, 1].gas.H2.boundaries[2, 2].Ndot} if anFP.n_y > 1 "Force state select";
  */

    Real states[:](each stateSelect=StateSelect.always) = {caFP.subregions[1, 1,
      1].gas.O2.boundaries[2, 2].Ndot,caFP.subregions[1, 1, 1].gas.H2O.boundaries[
      2, 2].Ndot,anFP.subregions[1, 1, 1].gas.H2.boundaries[2, 2].Ndot,caCL.subregions[
      1, 2, 1].ORR.'e-'.reaction.Ndot,caCL.subregions[1, 1, 1].gas.H2O.phi[1],
      caCL.subregions[1, 2, 1].gas.H2O.phi[1],PEM.subregions[1, 1, 1].ionomer.
      'H+'.phi[1]} if anFP.n_y > 1 "Force state select";

  initial equation
    /* **
  der(caCL.subregions[1, 1, 1].graphite.doubleLayer.w) = 0;
  der(caCL.subregions[1, 2, 1].graphite.doubleLayer.w) = 0;
  anCL.subregions[1, 1, 1].graphite.'e-Transfer'.I = 0;
  anCL.subregions[1, 2, 1].graphite.'e-Transfer'.I = 0;
  der(anCL.subregions[1, 1, 1].graphite.'e-Transfer'.I) = 0;
  der(anCL.subregions[1, 2, 1].graphite.'e-Transfer'.I) = 0;
    caCL.subregions[1, 1, 1].graphite.'e-Transfer'.I = 0;
  caCL.subregions[1, 2, 1].graphite.'e-Transfer'.I = 0;
  der(caCL.subregions[1, 1, 1].graphite.'e-Transfer'.I) = 0;
  der(caCL.subregions[1, 2, 1].graphite.'e-Transfer'.I) = 0;
    PEM.subregions[1, 1, 1].ionomer.'H+'.I[1] = U.mA;
  PEM.subregions[1, 2, 1].ionomer.'H+'.I[1] = U.mA;
*/

  equation
    connect(anCL.xPositive, PEM.xNegative) annotation (Line(
        points={{-10,40},{-10,40}},
        color={240,0,0},
        thickness=0.5,
        smooth=Smooth.None));
    connect(PEM.xPositive, caCL.xNegative) annotation (Line(
        points={{10,40},{10,40}},
        color={0,0,240},
        thickness=0.5,
        smooth=Smooth.None));
    connect(caCL.xPositive, caGDL.xNegative) annotation (Line(
        points={{30,40},{30,40}},
        color={0,0,240},
        thickness=0.5,
        smooth=Smooth.None));
    connect(anCL.xNegative, anGDL.xPositive) annotation (Line(
        points={{-30,40},{-30,40}},
        color={240,0,0},
        thickness=0.5,
        smooth=Smooth.None));
    connect(anBC1[1, :].boundary, anFP.xNegative[1, :]) annotation (Line(
        points={{-80,54},{-76,54},{-76,40},{-70,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(anBC2.boundary, anFP.xNegative[2:end, :]) annotation (Line(
        points={{-80,26},{-76,26},{-76,40},{-70,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(anSource.boundary, anFP.yNegative) annotation (Line(
        points={{-60,20},{-60,30}},
        color={240,0,0},
        thickness=0.5,
        smooth=Smooth.None));
    connect(anSink.boundary, anFP.yPositive) annotation (Line(
        points={{-60,60},{-60,50}},
        color={240,0,0},
        thickness=0.5,
        smooth=Smooth.None));
    connect(currentSet.y, zI) annotation (Line(
        points={{-59,-20},{-44,-20}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(zI, stoichH2.u) annotation (Line(
        points={{-44,-20},{-32,-20}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(stoichH2.y, Ndot_H2) annotation (Line(
        points={{-9,-20},{4,-20}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(Ndot_H2, anStoichH2O.u) annotation (Line(
        points={{4,-20},{28,-20}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(anStoichH2O.y, Ndot_H2O_an) annotation (Line(
        points={{51,-20},{64,-20}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(caBC1[1, :].boundary, caFP.xPositive[1, :]) annotation (Line(
        points={{80,54},{76,54},{76,40},{70,40}},
        color={127,127,127},
        thickness=0.5,
        smooth=Smooth.None));
    connect(caBC2.boundary, caFP.xPositive[2:end, :]) annotation (Line(
        points={{80,26},{76,26},{76,40},{70,40}},
        color={127,127,127},
        thickness=0.5,
        smooth=Smooth.None));
    connect(stoichO2.y, Ndot_O2) annotation (Line(
        points={{-9,-80},{4,-80}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(Ndot_O2, caStoichH2O.u) annotation (Line(
        points={{4,-80},{20,-80},{20,-60},{28,-60}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(caStoichH2O.y, Ndot_H2O_ca) annotation (Line(
        points={{51,-60},{64,-60}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(stoichN2.y, Ndot_N2) annotation (Line(
        points={{51,-100},{64,-100}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(stoichN2.u, Ndot_O2) annotation (Line(
        points={{28,-100},{20,-100},{20,-80},{4,-80}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(anGDL.xNegative, anFP.xPositive) annotation (Line(
        points={{-50,40},{-50,40}},
        color={240,0,0},
        thickness=0.5,
        smooth=Smooth.None));
    connect(caGDL.xPositive, caFP.xNegative) annotation (Line(
        points={{50,40},{50,40}},
        color={0,0,240},
        thickness=0.5,
        smooth=Smooth.None));
    connect(caFP.yNegative, caSource.boundary) annotation (Line(
        points={{60,30},{60,20}},
        color={0,0,240},
        thickness=0.5,
        smooth=Smooth.None));
    connect(caSink.boundary, caFP.yPositive) annotation (Line(
        points={{60,60},{60,50}},
        color={0,0,240},
        thickness=0.5,
        smooth=Smooth.None));
    connect(stoichO2.u, zI) annotation (Line(
        points={{-32,-80},{-40,-80},{-40,-20},{-44,-20}},
        color={0,0,127},
        smooth=Smooth.None));
    annotation (
      Commands(file=
            "Resources/Scripts/Dymola/Regions.Examples.FPtoFPVoltage.mos"
          "Regions.Examples.FPtoFPVoltage.mos", file=
            "Resources/Scripts/Dymola/Regions.Examples.FPtoFPVoltage-states.mos"
          "Regions.Examples.FPtoFPVoltage-states.mos"),
      experiment(StopTime=650, __Dymola_Algorithm="Dassl"),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-120},{
              100,100}}), graphics),
      __Dymola_experimentSetupOutput);
  end RegionsExamplesFPtoFPVoltage;

  model AssembliesCellsExamplesTestStandFourier
    "Simulate the fuel cell under prescribed conditions, with cyclical load"
    extends Assemblies.Cells.Examples.TestStand(redeclare
        Modelica.Electrical.Analog.Sources.SineCurrent load(
        I=48,
        freqHz=0.5,
        offset=50,
        startTime=5));
    Fourier fourier(F0=10, n=500) "FourierAnalysis"
      annotation (Placement(transformation(extent={{40,-70},{60,-50}})));
    annotation (experiment(
        StopTime=8,
        __Dymola_NumberOfIntervals=5000,
        Tolerance=1e-006,
        __Dymola_Algorithm="Dassl"), Commands(file=
            "Resources/Scripts/Dymola/Assemblies.Cells.Examples.TestStandCycle.mos"
          "Assemblies.Cells.Examples.TestStandCycle.mos"));

  end AssembliesCellsExamplesTestStandFourier;

  model Fourier "Compute Fourier coefficients of an input signal"

    parameter Modelica.SIunits.Frequency F0 "Base frequency for analysis";
    parameter Integer n "Number of harmonics";
    Modelica.Blocks.Interfaces.RealInput u "Input signal"
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
    Modelica.Blocks.Interfaces.RealOutput a0 "Signal bias"
      annotation (Placement(transformation(extent={{100,70},{120,90}})));
    Modelica.Blocks.Interfaces.RealOutput a[n]
      "Fourier coefficients for cosine terms"
      annotation (Placement(transformation(extent={{100,30},{120,50}})));
    Modelica.Blocks.Interfaces.RealOutput b[n]
      "Fourier coefficients for sine terms"
      annotation (Placement(transformation(extent={{100,-10},{120,10}})));
    import Modelica.Constants.pi;
    import Modelica.Math.atan2;
  protected
    parameter Modelica.SIunits.Time dt=1.0/F0 "Period at base frequency";
    Real s[n]={sin(2*pi*F0*i*time) for i in 1:n}
      "Sine waves at various frequencies";
    Real c[n]={cos(2*pi*F0*i*time) for i in 1:n}
      "Cosine waves at various frequencies";
    Real a0i "Integral of bias term";
    Real ai[n] "Integral of cosine terms";
    Real bi[n] "Integral of sine terms";
    Real f "Reconstructed function";
  public
    Modelica.Blocks.Interfaces.RealOutput mag[n] "Magnitude for each frequency"
      annotation (Placement(transformation(extent={{100,-50},{120,-30}})));
    Modelica.Blocks.Interfaces.RealOutput phase[n] "Phase for each frequency"
      annotation (Placement(transformation(extent={{100,-90},{120,-70}})));
  initial equation
    a0i = 0;
    ai = zeros(n);
    bi = zeros(n);
  equation
    der(a0i) = 2*u;
    der(ai) = 2*u*c;
    der(bi) = 2*u*s;
    f = a0/2 + a*c + b*s;
    when sample(0, dt) then
      a0 = pre(a0i);
      a = pre(ai);
      b = pre(bi);
      reinit(a0i, 0);
      for i in 1:n loop
        reinit(ai[i], 0);
        reinit(bi[i], 0);
      end for;
    end when;
    mag = {sqrt(a[i]^2 + b[i]^2) for i in 1:n};
    phase = {atan2(b[i], a[i]) for i in 1:n};
    annotation (
      Icon(graphics={Bitmap(extent={{-100,120},{100,-80}}, fileName=
            "modelica://FCSys/Resources/Documentation/Conditions/Fourier.jpg"),
            Text(
              extent={{-100,-58},{100,-98}},
              lineColor={0,0,255},
              textString="%name")}),
      Documentation(info="<html><p><i>This is a copy of the model by Mike Tiller.  It was copied from

    <a href=\"https://github.com/xogeny/Sensors\">https://github.com/xogeny/Sensors</a>

    on Dec. 4, 2013.  See also <a href=\"http://blog.xogeny.com/blog/fourier-analysis/\">http://blog.xogeny.com/blog/fourier-analysis/</a>.
  </i></p>

<p>This model performs a Fourier analysis of the input signal. The user must select a base frequency and the number of harmonics. The component will then analyze the signal over each period of the base frequency and output the coefficients for each harmonic as well as the steady state component of the input signal.</p>
<p>The outputs of this model are a0, a, and b where a and b are both vectors whose size is equal to the number of harmonics involved.</p>
<p>The original waveform can be reconstructed from the equation:</p>
<p>f = a0/2 + a[i]*sin(2*pi*F*i) + b[i]*cos(2*pi*F*i)</p>
<p>This equation is in identical notation, so there is an implied sum over i (from 1 to n).  Also note the factor of 1/2 on the a0 coefficient.</p>
</html>"),
      Diagram(graphics));
  end Fourier;
  annotation (Commands(
      file="../../units.mos"
        "Establish the constants and units in the workspace (first translate a model besides Units.Evaluate).",

      file="test/check.mos" "Check all of FCSys using Dymola's check function.",

      file="../../../LaTeX/Dissertation/Results/Cell/Simulation/sim.mos"), Icon(
        graphics={Polygon(
          points={{-80,-72},{0,68},{80,-72},{-80,-72}},
          lineColor={255,0,0},
          lineThickness=0.5)}));

end WorkInProgress;
