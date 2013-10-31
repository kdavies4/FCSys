within FCSys;
package Assemblies "Combinations of regions (e.g., cells)"
  extends Modelica.Icons.Package;

  package Cells "Single-cell PEMFC models"
    extends Modelica.Icons.Package;

    package Examples "Examples"

      extends Modelica.Icons.ExamplesPackage;

      model Polarization "Run a cell polarization"
        extends Modelica.Icons.Example;
        extends Modelica.Icons.UnderConstruction;
        /* **
                 params=dict(comp=['"O2"'],
                             anStoich=[1.5, 1.1, 2],
                             caStoich=[9.5, 7.5, 12.5],
                             anRH=[0.8, 0.6, 1],
                             caRH=[0.5, 0.3, 0.7],
                             T_degC=[60, 40, 80],
                             p_kPag=[48.3, 0, 202.7]),
                             */

        inner FCSys.Conditions.Environment environment
          "Environmental conditions"
          annotation (Placement(transformation(extent={{20,20},{40,40}})));

        Conditions.TestStands.TestStand testStand(
          zJ=currentDensity.y,
          anInletRH=0.8,
          caInletRH=0.5,
          T_an=333.15*U.K,
          T_ca=333.15*U.K,
          anStoich=1.5,
          caStoich=2.0) annotation (__Dymola_choicesFromPackage=true, Placement(
              transformation(extent={{-16,-16},{16,16}})));

        Cells.Cell cell "Fuel cell"
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      equation

        connect(testStand.an, cell.an) annotation (Line(
            points={{-16,9.4369e-16},{-14,9.4369e-16},{-14,6.10623e-16},{-10,
                6.10623e-16}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));

        connect(cell.ca, testStand.ca) annotation (Line(
            points={{10,6.10623e-16},{14,6.10623e-16},{14,9.4369e-16},{16.2,
                9.4369e-16}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));

        connect(testStand.anPositive, cell.anPositive) annotation (Line(
            points={{-4,16},{-4,10}},
            color={240,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(testStand.caNegative, cell.caNegative) annotation (Line(
            points={{4,-16},{4,-10}},
            color={0,0,240},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.anNegative, testStand.anNegative) annotation (Line(
            points={{-4,-10},{-4,-16}},
            color={240,0,0},
            thickness=0.5,
            smooth=Smooth.None));

        connect(cell.caPositive, testStand.caPositive) annotation (Line(
            points={{4,10},{4,16}},
            color={0,0,240},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (Commands(file=
                "Resources/Scripts/Dymola/Assemblies.Cells.Examples.Polarization.mos"
              "Assemblies.Cells.Examples.Polarization.mos"), Diagram(graphics));
      end Polarization;

      model EIS "Model for electrochemical-impedance spectroscopy"

        extends FCSys.Icons.Blocks.Continuous;

        parameter Modelica.SIunits.Current zI_large_A=100
          "Large-signal current in amperes";
        Connectors.RealInput zJ_small_SI
          "Small-signal current density in SI base units" annotation (Placement(
              transformation(extent={{-110,-10},{-90,10}}), iconTransformation(
                extent={{-120,-10},{-100,10}})));
        Connectors.RealOutput w_V "Cell potential in volts" annotation (
            Placement(transformation(extent={{90,-10},{110,10}}),
              iconTransformation(extent={{100,-10},{120,10}})));
        Conditions.TestStands.TestStandEIS testStandEIS
          annotation (Placement(transformation(extent={{-16,-16},{16,16}})));
        replaceable FCSys.Assemblies.Cells.Cell cell "Fuel cell model"
          annotation (__Dymola_choicesFromPackage=true, Placement(
              transformation(extent={{-10,-10},{10,10}})));
        inner FCSys.Conditions.Environment environment(
          analysis=false,
          p=149.6*U.kPa,
          T=333.15*U.K) "Environmental conditions"
          annotation (Placement(transformation(extent={{30,32},{50,52}})));
      equation
        connect(zJ_small_SI, testStandEIS.zJ_small_SI) annotation (Line(
            points={{-100,5.55112e-16},{-40,5.55112e-16},{-40,20},{-20,20},{-16,
                16.4},{-16.7,16.7}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(testStandEIS.w_V, w_V) annotation (Line(
            points={{16.7,-16.7},{16,-16},{20,-20},{40,-20},{40,0},{70,0},{70,
                5.55112e-16},{100,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(testStandEIS.caNegative, cell.caNegative) annotation (Line(
            points={{4,-16},{4,-10}},
            color={0,0,240},
            thickness=0.5,
            smooth=Smooth.None));
        connect(testStandEIS.anPositive, cell.anPositive) annotation (Line(
            points={{-4,16},{-4,10}},
            color={240,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(testStandEIS.an, cell.an) annotation (Line(
            points={{-16,9.4369e-16},{-10,6.10623e-16}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.anNegative, testStandEIS.anNegative) annotation (Line(
            points={{-4,-10},{-4,-16}},
            color={240,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.caPositive, testStandEIS.caPositive) annotation (Line(
            points={{4,10},{4,16}},
            color={0,0,240},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.ca, testStandEIS.ca) annotation (Line(
            points={{10,6.10623e-16},{16.2,6.10623e-16},{16.2,9.4369e-16}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (Diagram(graphics), Icon(graphics));
      end EIS;

      model FPtoFPCell "**"
        extends Modelica.Icons.Example;

        // temperature
        // pressure
        // RH
        /*  parameter Q.NumberAbsolute psi_O2_dry=0.21;
  parameter Q.NumberAbsolute anStoich=1.5 "Anode stoichiometric ratio";
  parameter Q.NumberAbsolute caStoich=2.0 "Cathode stoichiometric ratio";
  parameter Q.NumberAbsolute anInletRH=0.8;
  parameter Q.NumberAbsolute caInletRH=0.5;
  parameter Real T_degC=60;
  parameter Real p_kPag=48.3;
  final parameter Q.TemperatureAbsolute T=U.from_degC(T_degC);
  final parameter Q.PressureAbsolute p=U.from_kPag(p_kPag);
  final parameter Q.PressureAbsolute p_H2O_an_in=anInletRH*
      Characteristics.H2O.p_sat(T) "Pressure of H2O vapor";
  final parameter Q.PressureAbsolute p_H2O_ca_in=caInletRH*
      Characteristics.H2O.p_sat(T) "Pressure of H2O vapor";

  output Q.CurrentAreic zJ[n_y, n_z]=-anBC.graphite.'e-'.boundary.phi[1] .*
      anBC.graphite.'e-'.boundary.rho "Current density";
  output Q.Number zJ_Apercm2[n_y, n_z]=zJ*U.cm^2/U.A 
    "Electrical current density, in A/cm2";
  Q.Current I_H2_an_in=zI*anStoich/2;
  Q.Current I_O2_ca_in=zI*caStoich/4;

  parameter Q.Length L_y[:]=fill(U.m/1, 1) "Lengths along the channel" 
    annotation (Dialog(group="Geometry", __Dymola_label="<html><i>L</i><sub>y</sub></html>"))
    ;
  parameter Q.Length L_z[:]={5*U.mm} "Lengths across the channel" annotation (
      Dialog(group="Geometry", __Dymola_label="<html><i>L</i><sub>z</sub></html>"));
  final parameter Integer n_y=size(L_y, 1) 
    "Number of regions along the channel" annotation (HideResult=true);
  final parameter Integer n_z=size(L_z, 1) 
    "Number of regions across the channel" annotation (HideResult=true);
  final parameter Q.Area A_an_seg[cell.anFP.n_x, n_z]=outerProduct( cell.anFP.L_x, L_z) 
    "Areas of the segments over the xz plane of the anode flow plate";
  final parameter Q.Area A_ca_seg[cell.caFP.n_x, n_z]=outerProduct( cell.caFP.L_x, L_z) 
    "Areas of the segments over the xz plane of the cathode flow plate";
  final parameter Q.Area A_an=sum(cell.anFP.L_x)*sum(L_z) 
    "Total cross-sectional area of the anode flow plate in the xz plane";
  final parameter Q.Area A_ca=sum(cell.caFP.L_x)*sum(L_z) 
    "Total cross-sectional area of the cathode flow plate in the xz plane";

*/
        parameter Boolean inclLiq=false "Include liquid H2O";
        parameter Q.NumberAbsolute psi_H2O=environment.psi_H2O
          "Mole fraction of H2O at the inlet";
        parameter Q.NumberAbsolute psi_H2=environment.psi_dry
          "Mole fraction of H2 at the inlet";
        parameter Q.NumberAbsolute psi_O2=environment.psi_O2_dry*environment.psi_dry
          "Mole fraction of O2 at the inlet";
        parameter Q.NumberAbsolute psi_N2=(1 - environment.psi_O2_dry)*
            environment.psi_dry "Mole fraction of N2 at the inlet";
        output Q.Number J_Apercm2=zI*U.cm^2/(cell.caFP.A[Axis.x]*U.A)
          "Electrical current density, in A/cm2";
        output Q.Potential w=cell.anFP.subregions[1, 1, 1].graphite.'e-'.g_boundaries[
            1, Side.n] - cell.caFP.subregions[end, 1, 1].graphite.'e-'.g_boundaries[
            1, Side.p] if environment.analysis "Electrical potential";

        inner Conditions.Environment environment(
          analysis=true,
          RH=0.8,
          a={0,0,0}) "Environmental conditions"
          annotation (Placement(transformation(extent={{-10,70},{10,90}})));

        // Conditions
        Conditions.ByConnector.BoundaryBus.Single.Sink anBC[cell.n_y, cell.n_z]
          (each graphite(
            'incle-'=true,
            'e-'(materialSet(y=U.bar)),
            'inclC+'=true,
            redeclare FCSys.Conditions.ByConnector.ThermoDiffusive.Temperature
              'C+'(source(y=environment.T)))) annotation (Placement(
              transformation(
              extent={{10,-10},{-10,10}},
              rotation=90,
              origin={-84,40})));
        Conditions.ByConnector.BoundaryBus.Single.Source anSource[cell.anFP.n_x,
          cell.n_z](each gas(
            inclH2=true,
            inclH2O=true,
            H2(materialSet(y=-Ndot_H2), thermalSet(y=environment.T)),
            H2O(materialSet(y=-Ndot_H2O_an), thermalSet(y=environment.T))))
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={-60,16})));
        Conditions.ByConnector.BoundaryBus.Single.Sink anSink[cell.anFP.n_x,
          cell.n_z](gas(
            each inclH2=true,
            each inclH2O=true,
            H2O(materialSet(y=fill(
                          environment.p,
                          cell.anFP.n_x,
                          cell.n_z) - anSink.gas.H2.p)),
            H2(materialSet(y=anSink.gas.H2O.boundary.Ndot .* cell.anFP.subregions[
                    :, cell.n_y, :].gas.H2O.v ./ cell.anFP.subregions[:, cell.n_y,
                    :].gas.H2.v), redeclare each function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current)),
            each liquid(H2O(materialSet(y=environment.p)), inclH2O=inclLiq))
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-20,68})));
        Conditions.ByConnector.BoundaryBus.Single.Source caBC[cell.n_y, cell.n_z]
          (each graphite(
            'incle-'=true,
            'e-'(materialSet(y=-zI), thermalSet(y=environment.T)),
            'inclC+'=true,
            redeclare FCSys.Conditions.ByConnector.ThermoDiffusive.Temperature
              'C+'(source(y=environment.T)))) annotation (Placement(
              transformation(
              extent={{10,10},{-10,-10}},
              rotation=90,
              origin={84,40})));
        Conditions.ByConnector.BoundaryBus.Single.Source caSource[cell.caFP.n_x,
          cell.n_z](each gas(
            inclO2=true,
            inclN2=true,
            inclH2O=true,
            O2(materialSet(y=-Ndot_O2), thermalSet(y=environment.T)),
            N2(materialSet(y=-Ndot_N2), thermalSet(y=environment.T)),
            H2O(materialSet(y=-Ndot_H2O_ca), thermalSet(y=environment.T))))
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={60,16})));

        Conditions.ByConnector.BoundaryBus.Single.Sink caSink[cell.caFP.n_x,
          cell.n_z](gas(
            each inclO2=true,
            each inclN2=true,
            each inclH2O=true,
            H2O(materialSet(y=fill(
                          environment.p,
                          cell.caFP.n_x,
                          cell.n_z) - caSink.gas.N2.p - caSink.gas.O2.p)),
            N2(materialSet(y=caSink.gas.H2O.boundary.Ndot .* cell.caFP.subregions[
                    :, cell.n_y, :].gas.H2O.v ./ cell.caFP.subregions[:, cell.n_y,
                    :].gas.N2.v), redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current),

            O2(materialSet(y=caSink.gas.H2O.boundary.Ndot .* cell.caFP.subregions[
                    :, cell.n_y, :].gas.H2O.v ./ cell.caFP.subregions[:, cell.n_y,
                    :].gas.O2.v), redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current)),
            each liquid(H2O(materialSet(y=environment.p)), inclH2O=inclLiq))
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={60,64})));

        Modelica.Blocks.Sources.Ramp currentSet(
          offset=U.mA,
          height=100*U.A,
          duration=600)
          annotation (Placement(transformation(extent={{-80,-30},{-60,-10}})));
        Modelica.Blocks.Math.Gain stoichH2(k=1.2/2)
          annotation (Placement(transformation(extent={{-30,-30},{-10,-10}})));
        Modelica.Blocks.Math.Gain anStoichH2O(k=psi_H2O/psi_H2)
          annotation (Placement(transformation(extent={{30,-30},{50,-10}})));
        Modelica.Blocks.Math.Gain stoichO2(k=1.6/4)
          annotation (Placement(transformation(extent={{-30,-90},{-10,-70}})));
        Modelica.Blocks.Math.Gain caStoichH2O(k=psi_H2O/psi_O2)
          annotation (Placement(transformation(extent={{30,-70},{50,-50}})));
        Modelica.Blocks.Math.Gain stoichN2(k=psi_N2/psi_O2)
          annotation (Placement(transformation(extent={{30,-110},{50,-90}})));
        Assemblies.Cells.Cell cell(final inclLiq=inclLiq) annotation (Dialog,
            Placement(transformation(extent={{-10,30},{10,50}})));

      protected
        Connectors.RealOutputInternal Ndot_H2O_an(unit="N/T")
          "Rate of supply of H2O into the anode"
          annotation (Placement(transformation(extent={{54,-30},{74,-10}})));
        Connectors.RealOutputInternal Ndot_O2(unit="N/T")
          "Rate of supply of O2"
          annotation (Placement(transformation(extent={{-6,-90},{14,-70}})));
        Connectors.RealOutputInternal Ndot_H2O_ca(unit="N/T")
          "Rate of supply of H2O into the cathode"
          annotation (Placement(transformation(extent={{54,-70},{74,-50}})));
        Connectors.RealOutputInternal Ndot_N2(unit="N/T")
          "Rate of supply of N2"
          annotation (Placement(transformation(extent={{54,-110},{74,-90}})));
        Connectors.RealOutputInternal zI(unit="N/T") "Electrical current"
          annotation (Placement(transformation(extent={{-56,-30},{-36,-10}})));
        Connectors.RealOutputInternal Ndot_H2(unit="N/T")
          "Rate of supply of H2"
          annotation (Placement(transformation(extent={{-6,-30},{14,-10}})));

        Q.Velocity phi_states[:](
          each stateSelect=StateSelect.always,
          each start=0,
          each fixed=true) = {cell.anFP.subregions[1, 1, 1].gas.H2.phi[1],cell.anFP.subregions[
          1, 1, 1].gas.H2O.phi[1],cell.caFP.subregions[1, 1, 1].gas.H2O.phi[1],
          cell.caFP.subregions[1, 1, 1].gas.N2.phi[1],cell.caFP.subregions[1, 1,
          1].gas.O2.phi[1],cell.PEM.subregions[1, 1, 1].ionomer.H2O.phi[1],cell.caGDL.subregions[
          1, 1, 1].gas.H2O.phi[1],cell.caGDL.subregions[1, 1, 1].gas.N2.phi[1],
          cell.anGDL.subregions[1, 1, 1].gas.H2O.phi[1],cell.anGDL.subregions[1,
          1, 1].gas.H2.phi[1],cell.anCL.subregions[1, 1, 1].ionomer.H2O.phi[1]}
          "Velocities for state selection";

        //cell.caGDL.subregions[1, 1, 1].gas.O2.phi[1],

        // Note:  These must be forced to avoid dynamic state selection in Dymola 2014.
        Q.Velocity phi_states_liq[:](
          each stateSelect=StateSelect.always,
          each start=0,
          each fixed=true) = {cell.anGDL.subregions[1, 1, 1].liquid.H2O.phi[1],
          cell.anFP.subregions[1, 1, 1].liquid.H2O.phi[1],cell.caFP.subregions[
          1, 1, 1].liquid.H2O.phi[1],cell.caGDL.subregions[1, 1, 1].liquid.H2O.phi[
          1]} if inclLiq "Liquid velocities for state selection";
        Q.Current Ndot_states[:](
          each stateSelect=StateSelect.always,
          each start=0,
          each fixed=true) = {cell.anFP.subregions[1, 1, 1].gas.H2O.boundaries[
          2, 2].Ndot,cell.caFP.subregions[1, 1, 1].gas.H2O.boundaries[2, 2].Ndot}
          "Currents for state selection";

      equation
        connect(currentSet.y, zI) annotation (Line(
            points={{-59,-20},{-46,-20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(zI, stoichH2.u) annotation (Line(
            points={{-46,-20},{-32,-20}},
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
        connect(stoichO2.u, zI) annotation (Line(
            points={{-32,-80},{-40,-80},{-40,-20},{-46,-20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(anBC.boundary, cell.an) annotation (Line(
            points={{-80,40},{-10,40}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.ca, caBC.boundary) annotation (Line(
            points={{10,40},{80,40}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(anSink.boundary, cell.anPositive) annotation (Line(
            points={{-20,64},{-20,50},{-4,50}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.caPositive, caSink.boundary) annotation (Line(
            points={{4,50},{32,50},{32,60},{60,60}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.caNegative, caSource.boundary) annotation (Line(
            points={{4,30},{32,30},{32,20},{60,20}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.anNegative, anSource.boundary) annotation (Line(
            points={{-4,30},{-34,30},{-34,20},{-60,20}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (
          Commands(file="Resources/Scripts/Dymola/Regions.Examples.FPtoFP.mos"
              "Regions.Examples.FPtoFP.mos", file=
                "Resources/Scripts/Dymola/Regions.Examples.FPtoFP-states.mos"
              "Regions.Examples.FPtoFP-states.mos"),
          experiment(
            StopTime=600,
            Tolerance=1e-005,
            __Dymola_Algorithm="Dassl"),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-120},
                  {100,100}}), graphics),
          __Dymola_experimentSetupOutput);
      end FPtoFPCell;
    end Examples;

    model Cell "Single-cell PEMFC"
      extends FCSys.Icons.Cell;

      // Geometric parameters
      parameter Q.Length L_y[:]={U.m} "Lengths along the channel" annotation (
          Dialog(group="Geometry",__Dymola_label=
              "<html><i>L</i><sub>y</sub></html>"));
      parameter Q.Length L_z[:]={5*U.mm} "Lengths across the channel"
        annotation (Dialog(group="Geometry", __Dymola_label=
              "<html><i>L</i><sub>z</sub></html>"));
      final parameter Integer n_y=size(L_y, 1)
        "Number of subregions along the channel";
      final parameter Integer n_z=size(L_z, 1)
        "Number of subregions across the channel";

      // Assumptions
      parameter Boolean inclLiq=false "Include liquid H2O" annotation (Dialog(
            tab="Assumptions", compact=true), choices(__Dymola_checkBox=true));

      Connectors.BoundaryBus an[n_y, n_z] "Interface with the anode end plate"
        annotation (Placement(transformation(extent={{-100,-20},{-80,0}},
              rotation=0), iconTransformation(extent={{-110,-10},{-90,10}})));
      Connectors.BoundaryBus ca[n_y, n_z]
        "Interface with the cathode end plate" annotation (Placement(
            transformation(extent={{60,-20},{80,0}}, rotation=0),
            iconTransformation(extent={{90,-10},{110,10}})));
      Connectors.BoundaryBus anNegative[anFP.n_x, n_z]
        "Negative anode fluid port" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-70,-30}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-40,-100})));
      Connectors.BoundaryBus caNegative[caFP.n_x, n_z]
        "Negative cathode fluid port" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={50,-30}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={40,-100})));
      Connectors.BoundaryBus anPositive[anFP.n_x, n_z]
        "Positive anode fluid port" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-70,10}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-40,100})));
      Connectors.BoundaryBus caPositive[caFP.n_x, n_z]
        "Positive cathode fluid port" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={50,10}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={40,100})));

      replaceable Regions.AnFPs.AnFP anFP(
        final L_y=L_y,
        final L_z=L_z,
        subregions(each liquid(inclH2O=inclLiq))) "Anode flow plate"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-80,-20},{-60,0}})));
      replaceable Regions.AnGDLs.AnGDL anGDL(
        final L_y=L_y,
        final L_z=L_z,
        subregions(each liquid(inclH2O=inclLiq))) "Anode gas diffusion layer"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-60,-20},{-40,0}})));
      replaceable Regions.AnCLs.AnCL anCL(
        final L_y=L_y,
        final L_z=L_z,
        subregions(each liquid(inclH2O=inclLiq))) "Anode catalyst layer"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-40,-20},{-20,0}})));
      replaceable Regions.PEMs.PEM PEM(final L_y=L_y, final L_z=L_z)
        "Proton exchange membrane" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-20,-20},{0,0}})));
      replaceable Regions.CaCLs.CaCL caCL(
        final L_y=L_y,
        final L_z=L_z,
        subregions(each liquid(inclH2O=inclLiq))) "Cathode catalyst layer"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{0,-20},{20,0}})));
      replaceable Regions.CaGDLs.CaGDL caGDL(
        final L_y=L_y,
        final L_z=L_z,
        subregions(each liquid(inclH2O=inclLiq))) "Cathode gas diffusion layer"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{20,-20},{40,0}})));
      replaceable Regions.CaFPs.CaFP caFP(
        final L_y=L_y,
        final L_z=L_z,
        subregions(each liquid(inclH2O=inclLiq))) "Cathode flow plate"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{40,-20},{60,0}})));

    protected
      outer Conditions.Environment environment "Environmental conditions";

    equation
      // Internal connections (between layers)
      connect(anFP.xPositive, anGDL.xNegative) annotation (Line(
          points={{-60,-10},{-60,-10}},
          color={240,0,0},
          smooth=Smooth.None,
          thickness=0.5));
      connect(anGDL.xPositive, anCL.xNegative) annotation (Line(
          points={{-40,-10},{-40,-10}},
          color={240,0,0},
          smooth=Smooth.None,
          thickness=0.5));
      connect(anCL.xPositive, PEM.xNegative) annotation (Line(
          points={{-20,-10},{-20,-10}},
          color={240,0,0},
          smooth=Smooth.None,
          thickness=0.5));
      connect(PEM.xPositive, caCL.xNegative) annotation (Line(
          points={{0,-10},{0,-10}},
          color={0,0,240},
          smooth=Smooth.None,
          thickness=0.5));
      connect(caCL.xPositive, caGDL.xNegative) annotation (Line(
          points={{20,-10},{20,-10}},
          color={0,0,240},
          smooth=Smooth.None,
          thickness=0.5));
      connect(caGDL.xPositive, caFP.xNegative) annotation (Line(
          points={{40,-10},{40,-10}},
          color={0,0,240},
          smooth=Smooth.None,
          thickness=0.5));

      // External connections
      connect(an, anFP.xNegative) annotation (Line(
          points={{-90,-10},{-80,-10}},
          color={127,127,127},
          smooth=Smooth.None,
          thickness=0.5));

      connect(anFP.yNegative, anNegative) annotation (Line(
          points={{-70,-20},{-70,-30}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(anFP.yPositive, anPositive) annotation (Line(
          points={{-70,0},{-70,10}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));

      connect(caFP.xPositive, ca) annotation (Line(
          points={{60,-10},{70,-10}},
          color={127,127,127},
          smooth=Smooth.None,
          thickness=0.5));
      connect(caFP.yNegative, caNegative) annotation (Line(
          points={{50,-20},{50,-30}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      connect(caFP.yPositive, caPositive) annotation (Line(
          points={{50,0},{50,10}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        Documentation(info="
    <html><p>This is a model of a single-cell proton exchange membrane fuel cell (PEMFC).  An overview
    of a PEMFC is given in the <a href=\"modelica://FCSys\">top-level documentation of FCSys</a>.</p>
    </html>"),
        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-40},{
                80,20}}), graphics),
        Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.1), graphics={
            Line(
              points={{-40,-58},{-40,-100}},
              color={240,0,0},
              visible=inclY,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{-8,-1},{28,-1}},
              color={0,0,240},
              visible=inclX,
              thickness=0.5,
              origin={39,-92},
              rotation=90),
            Line(
              points={{-40,100},{-40,60}},
              color={240,0,0},
              visible=inclY,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{-66,0},{-100,0}},
              color={127,127,127},
              visible=inclX,
              thickness=0.5),
            Line(
              points={{-8,-1},{44,-1}},
              color={0,0,240},
              visible=inclX,
              thickness=0.5,
              origin={39,56},
              rotation=90),
            Line(
              points={{100,0},{56,0}},
              color={127,127,127},
              visible=inclX,
              thickness=0.5)}),
        experiment(StopTime=120, Tolerance=1e-06));

    end Cell;

    model SimpleCell
      "Cell model with integrated catalyst and gas diffusion layers"
      extends FCSys.Icons.Cell;

      // Geometric parameters
      parameter Q.Length L_y[:]={U.m} "Lengths along the channel" annotation (
          Dialog(group="Geometry",__Dymola_label=
              "<html><i>L</i><sub>y</sub></html>"));
      parameter Q.Length L_z[:]={5*U.mm} "Lengths across the channel"
        annotation (Dialog(group="Geometry", __Dymola_label=
              "<html><i>L</i><sub>z</sub></html>"));
      final parameter Integer n_y=size(L_y, 1)
        "Number of subregions along the channel";
      final parameter Integer n_z=size(L_z, 1)
        "Number of subregions across the channel";

      // Assumptions
      parameter Boolean inclLiq=false "Include liquid H2O" annotation (Dialog(
            tab="Assumptions", compact=true), choices(__Dymola_checkBox=true));

      Connectors.BoundaryBus an[n_y, n_z] "Interface with the anode end plate"
        annotation (Placement(transformation(extent={{-80,-20},{-60,0}},
              rotation=0), iconTransformation(extent={{-110,-10},{-90,10}})));
      Connectors.BoundaryBus ca[n_y, n_z]
        "Interface with the cathode end plate" annotation (Placement(
            transformation(extent={{40,-20},{60,0}}, rotation=0),
            iconTransformation(extent={{90,-10},{110,10}})));
      Connectors.BoundaryBus anNegative[anFP.n_x, n_z]
        "Negative anode fluid port" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-50,-30}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-40,-100})));
      Connectors.BoundaryBus caNegative[caFP.n_x, n_z]
        "Negative cathode fluid port" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={30,-30}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={40,-100})));
      Connectors.BoundaryBus anPositive[anFP.n_x, n_z]
        "Positive anode fluid port" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-50,10}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-40,100})));
      Connectors.BoundaryBus caPositive[caFP.n_x, n_z]
        "Positive cathode fluid port" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={30,10}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={40,100})));

      replaceable Regions.AnFPs.AnFP anFP(
        final L_y=L_y,
        final L_z=L_z,
        subregions(each liquid(inclH2O=inclLiq))) "Anode flow plate"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-60,-20},{-40,0}})));
      FCSys.Regions.AnCLs.AnCGDL anCGDL(
        final L_y=L_y,
        final L_z=L_z,
        subregions(each liquid(inclH2O=inclLiq)))
        "Anode catalyst and gas diffusion layer" annotation (Dialog(group=
              "Layers"), Placement(transformation(extent={{-40,-20},{-20,0}})));

      replaceable Regions.PEMs.PEM PEM(
        final L_y=L_y,
        final L_z=L_z,
        subregions(each liquid(inclH2O=inclLiq))) "Proton exchange membrane"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-20,-20},{0,0}})));
      FCSys.Regions.CaCLs.CaCGDL caCGDL(
        final L_y=L_y,
        final L_z=L_z,
        subregions(each liquid(inclH2O=inclLiq)))
        "Cathode catalyst and gas diffusion layer" annotation (Dialog(group=
              "Layers"), Placement(transformation(extent={{0,-20},{20,0}})));

      replaceable Regions.CaFPs.CaFP caFP(
        final L_y=L_y,
        final L_z=L_z,
        subregions(each liquid(inclH2O=inclLiq))) "Cathode flow plate"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{20,-20},{40,0}})));

    protected
      outer Conditions.Environment environment "Environmental conditions";

    equation
      // Internal connections (between layers)
      connect(anFP.xPositive, anCGDL.xNegative) annotation (Line(
          points={{-40,-10},{-40,-10}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(anCGDL.xPositive, PEM.xNegative) annotation (Line(
          points={{-20,-10},{-20,-10}},
          color={240,0,0},
          smooth=Smooth.None,
          thickness=0.5));
      connect(PEM.xPositive, caCGDL.xNegative) annotation (Line(
          points={{0,-10},{0,-10}},
          color={0,0,240},
          smooth=Smooth.None,
          thickness=0.5));
      connect(caCGDL.xPositive, caFP.xNegative) annotation (Line(
          points={{20,-10},{20,-10}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));

      // External connections
      connect(an, anFP.xNegative) annotation (Line(
          points={{-70,-10},{-60,-10}},
          color={127,127,127},
          smooth=Smooth.None,
          thickness=0.5));
      connect(anFP.yNegative, anNegative) annotation (Line(
          points={{-50,-20},{-50,-30}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(anFP.yPositive, anPositive) annotation (Line(
          points={{-50,0},{-50,10}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));

      connect(caFP.xPositive, ca) annotation (Line(
          points={{40,-10},{50,-10}},
          color={127,127,127},
          smooth=Smooth.None,
          thickness=0.5));
      connect(caFP.yNegative, caNegative) annotation (Line(
          points={{30,-20},{30,-30}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      connect(caFP.yPositive, caPositive) annotation (Line(
          points={{30,0},{30,10}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        defaultComponentName="cell",
        Documentation(info="<html>
    <p>This is a model of a single-cell proton exchange membrane fuel cell (PEMFC).  The catalyst layers and gas diffusion
    layers are integrated on each side to reduce the complexity of the model.  An overview
    of a PEMFC is given in the <a href=\"modelica://FCSys\">top-level documentation of FCSys</a>.</p>
    </html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-80,-40},{
                60,20}}), graphics),
        Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.1), graphics={
            Line(
              points={{-40,100},{-40,60}},
              color={240,0,0},
              visible=inclY,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{-8,-1},{44,-1}},
              color={0,0,240},
              visible=inclX,
              thickness=0.5,
              origin={39,56},
              rotation=90),
            Line(
              points={{100,0},{56,0}},
              color={127,127,127},
              visible=inclX,
              thickness=0.5),
            Line(
              points={{-8,-1},{28,-1}},
              color={0,0,240},
              visible=inclX,
              thickness=0.5,
              origin={39,-92},
              rotation=90),
            Line(
              points={{-40,-58},{-40,-100}},
              color={240,0,0},
              visible=inclY,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{-66,0},{-100,0}},
              color={127,127,127},
              visible=inclX,
              thickness=0.5)}));

    end SimpleCell;
  end Cells;
  annotation (Documentation(info="
<html>
  <p><b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
Copyright 2007&ndash;2013, <a href=\"http://www.gtrc.gatech.edu/\">Georgia Tech Research Corporation</a>.</p>

<p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.License\">
FCSys.UsersGuide.License</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p>
</html>"));

end Assemblies;
