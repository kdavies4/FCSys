within FCSys;
package Assemblies "Combinations of regions (e.g., cells)"
  extends Modelica.Icons.Package;

  package Cells "Single-cell PEMFC models"
    package Examples "Examples"

      extends Modelica.Icons.ExamplesPackage;

      model TestConditions "Fuel cell test conditions"
        extends FCSys.Icons.Record;
        import FCSys.Characteristics.H2O.p_sat;

        // Note:  This isn't a record because it contains time-varying variables.

        final parameter Q.NumberAbsolute psi_sat=environment.p_sat/environment.p
          "Mole fraction of H2O at saturation";

        // Electrical
        parameter Enumerations.ElectricalSpec electricalSpec=ElectricalSpec.current
          "Method of specification" annotation (Dialog(tab="Electrical"));
        replaceable Modelica.Blocks.Sources.RealExpression electricalSet
          constrainedby Modelica.Blocks.Interfaces.SO "Setpoint" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(tab="Electrical"),
          Placement(transformation(extent={{-30,-10},{-10,10}})));
        Q.Current zI "Current";
        Q.Power P "Power";
        //
        // Anode composition
        parameter Q.NumberAbsolute anRH(
          displayUnit="%",
          max=1) = 0.8 "Relative humidity (at inlet)"
          annotation (Dialog(tab="Anode"));
        final parameter Q.NumberAbsolute psi_H2O_an=anRH*psi_sat
          "Mole fraction of H2O at the anode inlet";
        final parameter Q.NumberAbsolute psi_H2=1 - psi_H2O_an
          "Mole fraction of H2 at the anode inlet";
        //
        // Cathode composition
        parameter Q.NumberAbsolute caRH(
          displayUnit="%",
          max=1) = 0.5 "Relative humidity (at inlet)"
          annotation (Dialog(tab="Cathode"));
        final parameter Q.NumberAbsolute psi_H2O_ca=caRH*psi_sat
          "Mole fraction of H2O at the cathode inlet";
        final parameter Q.NumberAbsolute psi_O2=environment.psi_O2_dry*(1 -
            psi_H2O_ca) "Mole fraction of O2 at the cathode inlet";
        final parameter Q.NumberAbsolute psi_N2=(1 - environment.psi_O2_dry)*(1
             - psi_H2O_ca) "Mole fraction of N2 at the cathode inlet";
        //
        // Anode flow rate
        parameter FCSys.Assemblies.Cells.Examples.Enumerations.FlowSpec
          anFlowSpec=FlowSpec.stoich "Method of specification"
          annotation (Dialog(tab="Anode",group="Rate of reactant supply"));
        replaceable Modelica.Blocks.Sources.RealExpression anFlowSet(y=1.5)
          constrainedby Modelica.Blocks.Interfaces.SO "Setpoint" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(tab="Anode",group="Rate of reactant supply"),
          Placement(transformation(extent={{10,0},{30,20}})));
        Q.NumberAbsolute anStoich "Anode stoichiometric flow rate";
        Q.Current I_an "Equivalent current of anode supply";
        //
        // Cathode flow rate
        parameter FCSys.Assemblies.Cells.Examples.Enumerations.FlowSpec
          caFlowSpec=FlowSpec.stoich "Method of specification"
          annotation (Dialog(tab="Cathode",group="Rate of reactant supply"));
        replaceable Modelica.Blocks.Sources.RealExpression caFlowSet(y=2.0)
          constrainedby Modelica.Blocks.Interfaces.SO "Setpoint" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(tab="Cathode",group="Rate of reactant supply"),
          Placement(transformation(extent={{10,-20},{30,0}})));
        Q.NumberAbsolute caStoich "Cathode stoichiometric flow rate";
        Q.Current I_ca "Equivalent current of cathode supply";

        // Measured conditions
        input Q.Potential w "Voltage"
          annotation (Dialog(tab="Measured conditions"));

      protected
        outer Conditions.Environment environment
          annotation (Placement(transformation(extent={{-10,30},{10,50}})));

      equation
        // Electrical
        // ----------
        // Relations
        P = w*zI;
        //
        // Setpoint
        if electricalSpec == ElectricalSpec.current then
          zI = electricalSet.y;
        elseif electricalSpec == ElectricalSpec.voltage then
          w = electricalSet.y;
        else
          // electricalSpec == ElectricalSpec.power
          P = electricalSet.y;
        end if;

        // Anode flow rate
        // ---------------
        // Relations
        anStoich*zI = I_an;
        //
        // Setpoint
        if anFlowSpec == FlowSpec.stoich then
          anStoich = anFlowSet.y;
        else
          // anFlowSpec == FlowSpec.current
          I_an = anFlowSet.y;
        end if;

        // Cathode flow rate
        // -----------------
        // Relations
        caStoich*zI = I_ca;
        //
        // Setpoint
        if caFlowSpec == FlowSpec.stoich then
          caStoich = caFlowSet.y;
        else
          // caFlowSpec == FlowSpec.current
          I_ca = caFlowSet.y;
        end if;
        annotation (Documentation(info="<html><p>Some conditions are taken from the outer <a href=\"modelica://FCSys.Conditions.Environment\">environment</a> model.  In particular,<ol>
    
    <li><code>environment.T</code> is used as the initial temperature throughout the cell, the temperature at each inlet, and the exterior temperature of each endplate in the yz plane.</li>
  
    <li><code>environment.p</code> is used as the initial pressure throughout the cell and the pressure at each outlet.</li>
  
    <li><code>environment.RH</code> is used as the initial relative humidity throughout the cell.</li>
 
    <li><code>environment.psi_O2_dry</code> is used as the dry-gas concentration of O<sub>2</sub> at the cathode inlet.</li>

 </ol></p>
 </html>"));
      end TestConditions;

      model TestStand2 "Simulate the fuel cell with prescribed conditions"
        import 'Datae-' = FCSys.Characteristics.'e-'.Graphite;
        import DataH2 = FCSys.Characteristics.H2.Gas;
        import DataH2O = FCSys.Characteristics.H2O.Gas;
        import DataO2 = FCSys.Characteristics.O2.Gas;
        import FCSys.Utilities.average;
        extends Modelica.Icons.Example;

        output Q.Number J_Apercm2=testConditions.zI*U.cm^2/(cell.caFP.A[Axis.x]
            *U.A) "Electrical current density, in A/cm2";
        output Q.Potential w=cell.anFP.subregions[1, 1, 1].graphite.'e-'.g_boundaries[
            1, Side.n] - cell.caFP.subregions[end, 1, 1].graphite.'e-'.g_boundaries[
            1, Side.p] if environment.analysis "Electrical potential";
        output Q.PressureAbsolute p_an_in=average(average(anSource.gas.H2.boundary.p
             + anSource.gas.H2O.boundary.p))
          "Total pressure at the anode inlet";
        output Q.PressureAbsolute p_ca_in=average(average(caSource.gas.H2O.boundary.p
             + p_N2_in + caSource.gas.O2.boundary.p))
          "Total pressure at the cathode inlet";
        output Q.Pressure Deltap_an=environment.p - p_an_in
          "Pressure difference down the anode channel";
        output Q.Pressure Deltap_ca=environment.p - p_ca_in
          "Pressure difference down the cathode channel";
        output Q.Potential Deltaw_O2=(DataO2.g(caSource[1, 1].gas.O2.boundary.T,
            caSource[1, 1].gas.O2.boundary.p) - cell.caCL.subregions[1, 1, 1].gas.O2.g)
            /4 if environment.analysis "Voltage loss due to O2 supply";
        output Q.Potential Deltaw_H2O=(DataH2O.g(caSource[1, 1].gas.H2O.boundary.T,
            caSource[1, 1].gas.H2O.boundary.p) - cell.caCL.subregions[1, 1, 1].gas.H2O.g)
            /2 if environment.analysis "Voltage loss due to H2O removal";
        output Q.Potential Deltaw_H2=(DataH2.g(anSource[1, 1].gas.H2.boundary.T,
            anSource[1, 1].gas.H2.boundary.p) - cell.anCL.subregions[1, 1, 1].gas.H2.g)
            /2 if environment.analysis "Voltage loss due to H2 supply";
        output Q.Potential 'Deltaw_e-'=cell.caFP.subregions[1, 1, 1].graphite.
            'e-'.g_boundaries[1, Side.p] - cell.caCL.subregions[1, 1, 1].graphite.
            'e-'.g + cell.anCL.subregions[1, 1, 1].graphite.'e-'.g - cell.anFP.subregions[
            1, 1, 1].graphite.'e-'.g_boundaries[1, Side.n] if environment.analysis
          "Voltage loss due to e- transport";
        output Q.Potential 'Deltaw_H+'=cell.anCL.subregions[1, 1, 1].ionomer.
            'H+'.g - cell.caCL.subregions[1, 1, 1].ionomer.'H+'.g if
          environment.analysis "Voltage loss due to H+ transport";
        output Q.Potential Deltaw_an=cell.anCL.subregions[1, 1, 1].graphite.
            'e-Transfer'.Deltag if environment.analysis "Anode overpotential";
        output Q.Potential Deltaw_ca=-cell.caCL.subregions[1, 1, 1].graphite.
            'e-Transfer'.Deltag if environment.analysis "Cathode overpotential";

        inner Conditions.Environment environment(
          a={0,0,0},
          T=333.13*U.K,
          p=U.from_kPag(48.3),
          RH=0.6) "Environmental conditions"
          annotation (Placement(transformation(extent={{-30,10},{-10,30}})));

        replaceable Cell cell(inclN2=environment.psi_O2_dry < 1 - Modelica.Constants.eps)
          constrainedby FCSys.Icons.Cell "Fuel cell" annotation (
          __Dymola_choicesFromPackage=true,
          choicesAllMatching=true,
          Placement(transformation(extent={{-10,-30},{10,-10}})));

        // Conditions
        Conditions.ByConnector.BoundaryBus.Single.Sink anBC[cell.n_y, cell.n_z]
          (
          each gas(
            inclH2=true,
            inclH2O=true,
            H2(
              redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current,

              materialSet(y=0),
              redeclare function afterSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Translational.velocity,

              redeclare function beforeSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Translational.velocity),

            H2O(
              redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current,

              materialSet(y=0),
              redeclare function afterSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Translational.velocity,

              redeclare function beforeSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Translational.velocity)),

          each graphite(
            'inclC+'=true,
            'incle-'=true,
            redeclare
              FCSys.Conditions.ByConnector.ThermalDiffusive.Single.Temperature
              'C+'(set(y=environment.T)),
            'e-'(materialSet(y=0))),
          each liquid(inclH2O=cell.inclLiq, H2O(
              redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current,

              materialSet(y=0),
              redeclare function thermalSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Thermal.heatRate,

              thermalSet(y=0)))) "Boundary condition for the cathode endplate"
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=90,
              origin={-24,-20})));

        Conditions.ByConnector.BoundaryBus.Single.Source anSource[cell.anFP.n_x,
          cell.n_z](each gas(
            inclH2=true,
            inclH2O=true,
            H2(materialSet(y=-testConditions.I_an/2), thermalSet(y=environment.T)),

            H2O(materialSet(y=-(testConditions.I_an/2)*(testConditions.psi_H2O_an
                    /testConditions.psi_H2)),thermalSet(y=environment.T))))
          "Source for the anode channel" annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={-56,-44})));

        Conditions.ByConnector.BoundaryBus.Single.Sink anSink[cell.anFP.n_x,
          cell.n_z](gas(
            each inclH2=true,
            each inclH2O=true,
            H2(redeclare each function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current,
                materialSet(y=anSink.gas.H2O.boundary.Ndot .* cell.anFP.subregions[
                    :, cell.n_y, :].gas.H2O.v ./ cell.anFP.subregions[:, cell.n_y,
                    :].gas.H2.v)),
            H2O(materialSet(y=fill(
                          environment.p,
                          cell.anFP.n_x,
                          cell.n_z) - anSink.gas.H2.p))), each liquid(inclH2O=
                cell.inclLiq, H2O(materialSet(y=environment.p))))
          "Sink for the anode channel" annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-56,4})));

        Conditions.ByConnector.BoundaryBus.Single.Source caBC[cell.n_y, cell.n_z]
          (
          each gas(
            inclH2O=true,
            inclN2=cell.inclN2,
            inclO2=true,
            H2O(redeclare function thermalSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
                thermalSet(y=0)),
            N2(redeclare function thermalSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
                thermalSet(y=0)),
            O2(redeclare function thermalSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
                thermalSet(y=0))),
          each graphite(
            'inclC+'=true,
            'incle-'=true,
            'C+'(set(y=environment.T)),
            'e-'(
              redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.pressure,

              materialSet(y=p_ca_elec),
              thermalSet(y=environment.T))),
          each liquid(inclH2O=cell.inclLiq, H2O(redeclare function thermalSpec
                = FCSys.Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
                thermalSet(y=0))))
          "Boundary condition for the cathode endplate" annotation (Placement(
              transformation(
              extent={{10,10},{-10,-10}},
              rotation=90,
              origin={24,-20})));

        Conditions.ByConnector.BoundaryBus.Single.Source caSource[cell.caFP.n_x,
          cell.n_z](gas(
            each inclH2O=true,
            each inclN2=cell.inclN2,
            each inclO2=true,
            each H2O(materialSet(y=-(testConditions.I_ca/4)*testConditions.psi_H2O_ca
                    /testConditions.psi_O2),thermalSet(y=environment.T)),
            N2(
              each materialSet(y=-(testConditions.I_ca/4)*testConditions.psi_N2
                    /testConditions.psi_O2),
              each thermalSet(y=environment.T),
              redeclare function materialMeas =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.pressure),

            each O2(materialSet(y=-testConditions.I_ca/4),thermalSet(y=
                    environment.T)))) "Source for the cathode channel"
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={56,-44})));

        Conditions.ByConnector.BoundaryBus.Single.Sink caSink[cell.caFP.n_x,
          cell.n_z](gas(
            each inclO2=true,
            each inclN2=cell.inclN2,
            each inclH2O=true,
            H2O(materialSet(y=fill(
                          environment.p,
                          cell.caFP.n_x,
                          cell.n_z) - p_N2_out - caSink.gas.O2.p), redeclare
                function materialMeas =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.volumeRate),

            N2(
              redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current,

              materialSet(y=caSink.gas.H2O.boundary.Ndot .* cell.caFP.subregions[
                    :, cell.n_y, :].gas.H2O.v ./ cell.caFP.subregions[:, cell.n_y,
                    :].gas.N2.v),
              redeclare function materialMeas =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.pressure),

            O2(redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current,
                materialSet(y=caSink.gas.H2O.boundary.Ndot .* cell.caFP.subregions[
                    :, cell.n_y, :].gas.H2O.v ./ cell.caFP.subregions[:, cell.n_y,
                    :].gas.O2.v))), each liquid(inclH2O=cell.inclLiq, H2O(
                materialSet(y=environment.p)))) "Sink for the cathode channel"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={56,4})));

        TestConditions testConditions(
          final zI=-sum(caBC.graphite.'e-'.boundary.Ndot),
          electricalSpec=FCSys.Assemblies.Cells.Examples.Enumerations.ElectricalSpec.current,

          redeclare replaceable Modelica.Blocks.Sources.Ramp electricalSet(
            offset=0,
            height=150*U.A,
            duration=3600,
            startTime=60)) "Test conditions" annotation (Dialog, Placement(
              transformation(extent={{10,10},{30,30}})));

      protected
        Q.Pressure p_ca_elec "Electronic pressure on the cathode";
        Conditions.Router anRouter[cell.anFP.n_x, cell.n_z] annotation (Dialog,
            Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=0,
              origin={-48,-20})));
        Conditions.Router caRouter[cell.anFP.n_x, cell.n_z] annotation (Dialog,
            Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=0,
              origin={48,-20})));

        Connectors.RealInput p_N2_in[cell.caFP.n_x, cell.n_z](each unit=
              "M/(L.T2)") "Pressure of N2 at inlet";
        Connectors.RealInput p_N2_out[cell.caFP.n_x, cell.n_z](each unit=
              "M/(L.T2)") "Pressure of N2 at outlet";

      equation
        testConditions.w = 'Datae-'.g(environment.T, 0) - 'Datae-'.g(
          environment.T, p_ca_elec);

        // Nitrogen pressures (since N2 is conditionally included)
        connect(p_N2_in, caSource.gas.N2.materialOut.y);
        connect(p_N2_out, caSink.gas.N2.materialOut.y);
        if not cell.inclN2 then
          p_N2_in = zeros(cell.caFP.n_x, cell.n_y);
          p_N2_out = zeros(cell.caFP.n_x, cell.n_y);
        end if;

        connect(cell.ca, caBC.boundary) annotation (Line(
            points={{10,-20},{20,-20}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.an, anBC.boundary) annotation (Line(
            points={{-10,-20},{-20,-20}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(anRouter.positive2, anSink.boundary) annotation (Line(
            points={{-56,-16},{-56,0}},
            color={240,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(anRouter.positive1, anSource.boundary) annotation (Line(
            points={{-56,-24},{-56,-40}},
            color={240,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.anPositive, anRouter.negative2) annotation (Line(
            points={{-4,-10},{-4,0},{-40,0},{-40,-16}},
            color={240,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.anNegative, anRouter.negative1) annotation (Line(
            points={{-4,-30},{-4,-40},{-40,-40},{-40,-24}},
            color={240,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(caSink.boundary, caRouter.negative2) annotation (Line(
            points={{56,0},{56,-16}},
            color={0,0,240},
            thickness=0.5,
            smooth=Smooth.None));
        connect(caSource.boundary, caRouter.negative1) annotation (Line(
            points={{56,-40},{56,-24}},
            color={0,0,240},
            thickness=0.5,
            smooth=Smooth.None));
        connect(caRouter.positive2, cell.caPositive) annotation (Line(
            points={{40,-16},{40,0},{4,0},{4,-10}},
            color={0,0,240},
            thickness=0.5,
            smooth=Smooth.None));
        connect(caRouter.positive1, cell.caNegative) annotation (Line(
            points={{40,-24},{40,-40},{4,-40},{4,-30}},
            color={0,0,240},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (
          Documentation(info="<html>
    
    <p>Please see <a href=\"FCSys.Assemblies.Cells.Examples.TestConditions\">TestConditions</a> regarding 
    how the properties of the environment are used.</p>
    
    <p>To run the cell with pure O<sub>2</sub> in the cathode (no N<sub>2</sub>), set <code>environment.psi_O2_dry</code> to 100%.</p>
    
    <p>Assumptions:
    <ol>
    <li>The outer surface of each end plate has uniform temperature in the yz plane.</li>
    <li>No heat is conducted from the rest of the cell hardware.</li>
    <li>All electronic current passes through the first segment (index <code>[1, 1]</code>) of each end plate in the yz plane.</li>
    <li>There is no shear force on the fluid at either outlet.</li>
    <li>The pressure of each gas species is uniform over each inlet.</li>
    <li>The volumetric flow rates of the gas species are equal at each outlet, where the volumetric flow rate 
    is approximated using the current at the outlet and the density of the gas within the last segment of the cell.</li>
    <li>The outlet pressure is applied to the gas mixture by Dalton's law (additivity of pressure).</li>
    <li>At the outlet, the liquid has the same pressure as the gas (Amagat's law).</li>
    <li>There is no thermal conduction across either outlet.</li>
    </ol></p>
    </html>"),
          Commands(file=
                "Resources/Scripts/Dymola/Assemblies.Cells.Examples.TestStand.mos"
              "Assemblies.Cells.Examples.TestStand.mos", file=
                "Resources/Scripts/Dymola/Assemblies.Cells.Examples.TestStand-states.mos"
              "Assemblies.Cells.Examples.TestStand-states.mos"),
          experiment(StopTime=3650, __Dymola_Algorithm="Dassl"),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-80,-60},
                  {80,40}}), graphics),
          __Dymola_experimentSetupOutput);
      end TestStand2;

      package Enumerations "Choices of options"

        extends Modelica.Icons.BasesPackage;

        type ElectricalSpec = enumeration(
            current "Current",
            voltage "Voltage",
            power "Power") "Ways to specify the electrical load";
        type FlowSpec = enumeration(
            stoich "Stoichiometric rate",
            current "Equivalent current",
            pressure "Inlet pressure") "Ways to specify the anode flow rate";

      end Enumerations;

      model TestStandSimple
        "Simulate the simple fuel cell model with prescribed conditions"

        import DataH2 = FCSys.Characteristics.H2.Gas;
        import DataH2O = FCSys.Characteristics.H2O.Gas;
        import DataO2 = FCSys.Characteristics.O2.Gas;

        extends TestStand(
          redeclare SimpleCell cell,
          Deltaw_O2=(DataO2.g(caSource[1, 1].gas.O2.boundary.T, caSource[1, 1].gas.O2.boundary.p)
               - cell.caCGDL.subregions[1, 1, 1].gas.O2.g)/4,
          Deltaw_H2O=(DataH2O.g(caSource[1, 1].gas.H2O.boundary.T, caSource[1,
              1].gas.H2O.boundary.p) - cell.caCGDL.subregions[1, 1, 1].gas.H2O.g)
              /2,
          Deltaw_H2=(DataH2.g(anSource[1, 1].gas.H2.boundary.T, anSource[1, 1].gas.H2.boundary.p)
               - cell.anCGDL.subregions[1, 1, 1].gas.H2.g)/2,
          'Deltaw_e-'=cell.caFP.subregions[1, 1, 1].graphite.'e-'.g_boundaries[
              1, Side.p] - cell.caCGDL.subregions[1, 1, 1].graphite.'e-'.g +
              cell.anCGDL.subregions[1, 1, 1].graphite.'e-'.g - cell.anFP.subregions[
              1, 1, 1].graphite.'e-'.g_boundaries[1, Side.n],
          'Deltaw_H+'=cell.anCGDL.subregions[1, 1, 1].ionomer.'H+'.g - cell.caCGDL.subregions[
              1, 1, 1].ionomer.'H+'.g,
          Deltaw_an=cell.anCGDL.subregions[1, 1, 1].graphite.'e-Transfer'.Deltag,

          Deltaw_ca=-cell.caCGDL.subregions[1, 1, 1].graphite.'e-Transfer'.Deltag);

        annotation (
          Commands(file=
                "Resources/Scripts/Dymola/Assemblies.Cells.Examples.TestStandSimple.mos"
              "Assemblies.Cells.Examples.TestStandSimple.mos"),
          experiment(
            StopTime=3650,
            Tolerance=1e-005,
            __Dymola_Algorithm="Dassl"),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-80,-60},
                  {80,40}}), graphics),
          __Dymola_experimentSetupOutput);
      end TestStandSimple;

      model TestStand "Simulate the fuel cell with prescribed conditions"
        import 'Datae-' = FCSys.Characteristics.'e-'.Graphite;
        import DataH2 = FCSys.Characteristics.H2.Gas;
        import DataH2O = FCSys.Characteristics.H2O.Gas;
        import DataO2 = FCSys.Characteristics.O2.Gas;
        import FCSys.Utilities.average;
        extends Modelica.Icons.Example;

        Q.Potential w "Electrical potential";
        output Q.Current zI=testConditions.zI if environment.analysis
          "Electrical current";
        output Q.Number J_Apercm2=zI*U.cm^2/(cell.caFP.A[Axis.x]*U.A) if
          environment.analysis "Electrical current density, in A/cm2";
        output Q.PressureAbsolute p_an_in=average(average(anSource.gas.H2.boundary.p
             + anSource.gas.H2O.boundary.p)) if environment.analysis
          "Total pressure at the anode inlet";
        output Q.PressureAbsolute p_ca_in=average(average(caSource.gas.H2O.boundary.p
             + p_N2_in + caSource.gas.O2.boundary.p)) if environment.analysis
          "Total pressure at the cathode inlet";
        output Q.Pressure Deltap_an=environment.p - p_an_in if environment.analysis
          "Pressure difference down the anode channel";
        output Q.Pressure Deltap_ca=environment.p - p_ca_in if environment.analysis
          "Pressure difference down the cathode channel";
        output Q.Potential Deltaw_O2=(DataO2.g(caSource[1, 1].gas.O2.boundary.T,
            caSource[1, 1].gas.O2.boundary.p) - cell.caCL.subregions[1, 1, 1].gas.O2.g)
            /4 if environment.analysis "Voltage loss due to O2 supply";
        output Q.Potential Deltaw_H2O=(DataH2O.g(caSource[1, 1].gas.H2O.boundary.T,
            caSource[1, 1].gas.H2O.boundary.p) - cell.caCL.subregions[1, 1, 1].gas.H2O.g)
            /2 if environment.analysis "Voltage loss due to H2O removal";
        output Q.Potential Deltaw_H2=(DataH2.g(anSource[1, 1].gas.H2.boundary.T,
            anSource[1, 1].gas.H2.boundary.p) - cell.anCL.subregions[1, 1, 1].gas.H2.g)
            /2 if environment.analysis "Voltage loss due to H2 supply";
        output Q.Potential 'Deltaw_e-'=cell.caFP.subregions[1, 1, 1].graphite.
            'e-'.g_boundaries[1, Side.p] - cell.caCL.subregions[1, 1, 1].graphite.
            'e-'.g + cell.anCL.subregions[1, 1, 1].graphite.'e-'.g - cell.anFP.subregions[
            1, 1, 1].graphite.'e-'.g_boundaries[1, Side.n] if environment.analysis
          "Voltage loss due to e- transport";
        output Q.Potential 'Deltaw_H+'=cell.anCL.subregions[1, 1, 1].ionomer.
            'H+'.g - cell.caCL.subregions[1, 1, 1].ionomer.'H+'.g if
          environment.analysis "Voltage loss due to H+ transport";
        output Q.Potential Deltaw_an=cell.anCL.subregions[1, 1, 1].graphite.
            'e-Transfer'.Deltag if environment.analysis "Anode overpotential";
        output Q.Potential Deltaw_ca=-cell.caCL.subregions[1, 1, 1].graphite.
            'e-Transfer'.Deltag if environment.analysis "Cathode overpotential";

        inner Conditions.Environment environment(
          a={0,0,0},
          T=333.13*U.K,
          p=U.from_kPag(48.3),
          RH=0.6) "Environmental conditions"
          annotation (Placement(transformation(extent={{-30,30},{-10,50}})));

        replaceable Cell cell(inclN2=environment.psi_O2_dry < 1 - Modelica.Constants.eps,
            L_y={8}*U.cm) constrainedby FCSys.Icons.Cell "Fuel cell"
          annotation (
          __Dymola_choicesFromPackage=true,
          choicesAllMatching=true,
          Placement(transformation(extent={{-10,-10},{10,10}})));

        // Conditions
        Conditions.ByConnector.BoundaryBus.Single.Sink anBC[cell.n_y, cell.n_z]
          (
          each gas(
            inclH2=true,
            inclH2O=true,
            H2(
              redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current,

              materialSet(y=0),
              redeclare function afterSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Translational.velocity,

              redeclare function beforeSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Translational.velocity),

            H2O(
              redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current,

              materialSet(y=0),
              redeclare function afterSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Translational.velocity,

              redeclare function beforeSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Translational.velocity)),

          each graphite('inclC+'=true, redeclare
              FCSys.Conditions.ByConnector.ThermalDiffusive.Single.Temperature
              'C+'(set(y=environment.T))),
          each liquid(inclH2O=cell.inclLiq, H2O(
              redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current,

              materialSet(y=0),
              redeclare function thermalSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Thermal.heatRate,

              thermalSet(y=0))))
          "Boundary condition for the anode endplate, except electrical"
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=90,
              origin={-24,0})));

        Conditions.ByConnector.BoundaryBus.Single.Source anSource[cell.anFP.n_x,
          cell.n_z](each gas(
            inclH2=true,
            inclH2O=true,
            H2(materialSet(y=-testConditions.I_an/2), thermalSet(y=environment.T)),

            H2O(materialSet(y=-(testConditions.I_an/2)*(testConditions.psi_H2O_an
                    /testConditions.psi_H2)),thermalSet(y=environment.T))))
          "Source for the anode reactant stream" annotation (Placement(
              transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={-56,-24})));

        Conditions.ByConnector.BoundaryBus.Single.Sink anSink[cell.anFP.n_x,
          cell.n_z](gas(
            each inclH2=true,
            each inclH2O=true,
            H2(redeclare each function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current,
                materialSet(y=anSink.gas.H2O.boundary.Ndot .* cell.anFP.subregions[
                    :, cell.n_y, :].gas.H2O.v ./ cell.anFP.subregions[:, cell.n_y,
                    :].gas.H2.v)),
            H2O(materialSet(y=fill(
                          environment.p,
                          cell.anFP.n_x,
                          cell.n_z) - anSink.gas.H2.p))), each liquid(inclH2O=
                cell.inclLiq, H2O(materialSet(y=environment.p))))
          "Sink for the anode reactant stream" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-56,24})));

        Conditions.ByConnector.BoundaryBus.Single.Source caBC[cell.n_y, cell.n_z]
          (
          each gas(
            inclH2O=true,
            inclN2=cell.inclN2,
            inclO2=true,
            H2O(redeclare function thermalSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
                thermalSet(y=0)),
            N2(redeclare function thermalSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
                thermalSet(y=0)),
            O2(redeclare function thermalSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
                thermalSet(y=0))),
          each graphite('inclC+'=true, 'C+'(set(y=environment.T))),
          each liquid(inclH2O=cell.inclLiq, H2O(redeclare function thermalSpec
                = FCSys.Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
                thermalSet(y=0))))
          "Boundary condition for the cathode endplate, except electrical"
          annotation (Placement(transformation(
              extent={{10,10},{-10,-10}},
              rotation=90,
              origin={24,0})));

        Conditions.ByConnector.BoundaryBus.Single.Source caSource[cell.caFP.n_x,
          cell.n_z](gas(
            each inclH2O=true,
            each inclN2=cell.inclN2,
            each inclO2=true,
            each H2O(materialSet(y=-(testConditions.I_ca/4)*testConditions.psi_H2O_ca
                    /testConditions.psi_O2),thermalSet(y=environment.T)),
            N2(
              each materialSet(y=-(testConditions.I_ca/4)*testConditions.psi_N2
                    /testConditions.psi_O2),
              each thermalSet(y=environment.T),
              redeclare function materialMeas =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.pressure),

            each O2(materialSet(y=-testConditions.I_ca/4),thermalSet(y=
                    environment.T)))) "Source for the cathode reactant stream"
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={56,-24})));

        Conditions.ByConnector.BoundaryBus.Single.Sink caSink[cell.caFP.n_x,
          cell.n_z](gas(
            each inclO2=true,
            each inclN2=cell.inclN2,
            each inclH2O=true,
            H2O(materialSet(y=fill(
                          environment.p,
                          cell.caFP.n_x,
                          cell.n_z) - p_N2_out - caSink.gas.O2.p), redeclare
                function materialMeas =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.volumeRate),

            N2(
              redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current,

              materialSet(y=caSink.gas.H2O.boundary.Ndot .* cell.caFP.subregions[
                    :, cell.n_y, :].gas.H2O.v ./ cell.caFP.subregions[:, cell.n_y,
                    :].gas.N2.v),
              redeclare function materialMeas =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.pressure),

            O2(redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current,
                materialSet(y=caSink.gas.H2O.boundary.Ndot .* cell.caFP.subregions[
                    :, cell.n_y, :].gas.H2O.v ./ cell.caFP.subregions[:, cell.n_y,
                    :].gas.O2.v))), each liquid(inclH2O=cell.inclLiq, H2O(
                materialSet(y=environment.p))))
          "Sink for the cathode reactant stream" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={56,24})));

        TestConditions testConditions(
          final w=w,
          electricalSpec=FCSys.Assemblies.Cells.Examples.Enumerations.ElectricalSpec.current,

          redeclare replaceable Modelica.Blocks.Sources.Ramp electricalSet(
            offset=U.mA,
            height=150*U.A,
            duration=3600,
            startTime=60)) "Test conditions" annotation (Dialog, Placement(
              transformation(extent={{10,30},{30,50}})));

        Real states[:](each stateSelect=StateSelect.always) = {cell.PEM.subregions[
          1, 2, 1].ionomer.'H+'.I[1],cell.caFP.subregions[1, 2, 1].gas.O2.phi[2],
          cell.caFP.subregions[1, 2, 1].gas.H2O.phi[2],cell.anFP.subregions[1,
          2, 1].gas.H2.phi[2]} if cell.n_y > 1
          "**clean up/integrate into layers";

      protected
        Conditions.Router anRouter[cell.anFP.n_x, cell.n_z]
          "Switch to route the anode for reverse flow" annotation (Dialog,
            Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=0,
              origin={-48,0})));
        Conditions.Router caRouter[cell.anFP.n_x, cell.n_z]
          "Switch to route the cathode for reverse flow" annotation (Dialog,
            Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=0,
              origin={48,0})));

        Connectors.RealInput p_N2_in[cell.caFP.n_x, cell.n_z](each unit=
              "M/(L.T2)") "Pressure of N2 at inlet";
        Connectors.RealInput p_N2_out[cell.caFP.n_x, cell.n_z](each unit=
              "M/(L.T2)") "Pressure of N2 at outlet";

      public
        replaceable Modelica.Electrical.Analog.Sources.RampCurrent load(
          duration=3600,
          offset=0.001,
          startTime=60,
          I=150) constrainedby Modelica.Electrical.Analog.Interfaces.TwoPin
          "Electrical load"
          annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
        Modelica.Electrical.Analog.Basic.Ground ground
          annotation (Placement(transformation(extent={{-20,-80},{0,-60}})));
        Conditions.Adapters.MSL.ElectronicOnly anAdapt
          "Interface with Modelica.Electrical on the anode side" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-10,-30})));
        Conditions.Adapters.MSL.ElectronicOnly caAdapt
          "Interface with Modelica.Electrical on the cathode side" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={10,-30})));

      equation
        w = load.v*U.V;

        // Nitrogen pressures (since N2 is conditionally included)
        connect(p_N2_in, caSource.gas.N2.materialOut.y) "Not shown in diagram";
        connect(p_N2_out, caSink.gas.N2.materialOut.y) "Not shown in diagram";
        if not cell.inclN2 then
          p_N2_in = zeros(cell.caFP.n_x, cell.n_y);
          p_N2_out = zeros(cell.caFP.n_x, cell.n_y);
        end if;

        connect(cell.an[1, 1], anAdapt.boundary) annotation (Line(
            points={{-10,0},{-10,-26}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.ca[1, 1], caAdapt.boundary) annotation (Line(
            points={{10,0},{10,-26}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.an, anBC.boundary) annotation (Line(
            points={{-10,0},{-20,0}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.ca, caBC.boundary) annotation (Line(
            points={{10,0},{20,0}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(anRouter.positive2, anSink.boundary) annotation (Line(
            points={{-56,4},{-56,20}},
            color={240,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(anRouter.positive1, anSource.boundary) annotation (Line(
            points={{-56,-4},{-56,-20}},
            color={240,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.anPositive, anRouter.negative2) annotation (Line(
            points={{-4,10},{-4,20},{-40,20},{-40,4}},
            color={240,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.anNegative, anRouter.negative1) annotation (Line(
            points={{-4,-10},{-4,-20},{-40,-20},{-40,-4}},
            color={240,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(caSink.boundary, caRouter.negative2) annotation (Line(
            points={{56,20},{56,4}},
            color={0,0,240},
            thickness=0.5,
            smooth=Smooth.None));
        connect(caSource.boundary, caRouter.negative1) annotation (Line(
            points={{56,-20},{56,-4}},
            color={0,0,240},
            thickness=0.5,
            smooth=Smooth.None));
        connect(caRouter.positive2, cell.caPositive) annotation (Line(
            points={{40,4},{40,20},{4,20},{4,10}},
            color={0,0,240},
            thickness=0.5,
            smooth=Smooth.None));
        connect(caRouter.positive1, cell.caNegative) annotation (Line(
            points={{40,-4},{40,-20},{4,-20},{4,-10}},
            color={0,0,240},
            thickness=0.5,
            smooth=Smooth.None));
        connect(load.p, ground.p) annotation (Line(
            points={{-10,-50},{-10,-60}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(anAdapt.pin, load.p) annotation (Line(
            points={{-10,-34},{-10,-50}},
            color={0,0,255},
            smooth=Smooth.None));
        connect(caAdapt.pin, load.n) annotation (Line(
            points={{10,-34},{10,-50}},
            color={0,0,255},
            smooth=Smooth.None));
        annotation (
          Documentation(info="<html>
    
    <p>Please see <a href=\"FCSys.Assemblies.Cells.Examples.TestConditions\">TestConditions</a> regarding 
    how the properties of the environment are used.</p>
    
    <p>To run the cell with pure O<sub>2</sub> in the cathode (no N<sub>2</sub>), set <code>environment.psi_O2_dry</code> to 100%.</p>
    
    <p>Assumptions:
    <ol>
    <li>The outer surface of each end plate has uniform temperature in the yz plane.</li>
    <li>No heat is conducted from the rest of the cell hardware.</li>
    <li>All electronic current passes through the first segment (index <code>[1, 1]</code>) of each end plate in the yz plane.</li>
    <li>There is no shear force on the fluid at either outlet.</li>
    <li>The pressure of each gas species is uniform over each inlet.</li>
    <li>The volumetric flow rates of the gas species are equal at each outlet, where the volumetric flow rate 
    is approximated using the current at the outlet and the density of the gas within the last segment of the cell.</li>
    <li>The outlet pressure is applied to the gas mixture by Dalton's law (additivity of pressure).</li>
    <li>At the outlet, the liquid has the same pressure as the gas (Amagat's law).</li>
    <li>There is no thermal conduction across either outlet.</li>
    </ol></p>
    </html>"),
          Commands(file=
                "Resources/Scripts/Dymola/Assemblies.Cells.Examples.TestStand.mos"
              "Assemblies.Cells.Examples.TestStand.mos", file=
                "Resources/Scripts/Dymola/Assemblies.Cells.Examples.TestStand-states.mos"
              "Assemblies.Cells.Examples.TestStand-states.mos"),
          experiment(StopTime=3660, __Dymola_Algorithm="Dassl"),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-80,-80},
                  {80,60}}), graphics),
          __Dymola_experimentSetupOutput);
      end TestStand;

      model TestStandVoltage
        "Simulate the fuel cell with prescribed conditions, voltage"
        import 'Datae-' = FCSys.Characteristics.'e-'.Graphite;
        import DataH2 = FCSys.Characteristics.H2.Gas;
        import DataH2O = FCSys.Characteristics.H2O.Gas;
        import DataO2 = FCSys.Characteristics.O2.Gas;
        import FCSys.Utilities.average;
        extends Modelica.Icons.Example;

        parameter Boolean inclLiq=false "Include liquid H2O";
        parameter Q.NumberAbsolute psi_H2O=(testConditions.psi_H2O_an +
            testConditions.psi_H2O_ca)/2
          "Average molar concentration of H2O at the inlets";
        parameter Q.NumberAbsolute RH=(testConditions.anRH + testConditions.caRH)
            /2 "Average relative humidities at the inlets";

        output Q.Number J_Apercm2=testConditions.zI*U.cm^2/(cell.caFP.A[Axis.x]
            *U.A) "Electrical current density, in A/cm2";
        output Q.Potential w=cell.anFP.subregions[1, 1, 1].graphite.'e-'.g_boundaries[
            1, Side.n] - cell.caFP.subregions[end, 1, 1].graphite.'e-'.g_boundaries[
            1, Side.p] if environment.analysis "Electrical potential";
        output Q.PressureAbsolute p_an_in=average(average(anSource.gas.H2.boundary.p
             + anSource.gas.H2O.boundary.p))
          "Total pressure at the anode inlet";
        output Q.PressureAbsolute p_ca_in=average(average(caSource.gas.H2O.boundary.p
             + caSource.gas.N2.boundary.p + caSource.gas.O2.boundary.p))
          "Total pressure at the cathode inlet";
        output Q.Pressure Deltap_an=testConditions.p - p_an_in
          "Pressure difference down the anode channel";
        output Q.Pressure Deltap_ca=testConditions.p - p_ca_in
          "Pressure difference down the cathode channel";
        output Q.Potential Deltaw_O2=(DataO2.g(caSource[1, 1].gas.O2.boundary.T,
            caSource[1, 1].gas.O2.boundary.p) - cell.caCL.subregions[1, 1, 1].gas.O2.g)
            /4 if environment.analysis "Voltage loss due to O2 supply";
        output Q.Potential Deltaw_H2O=(DataH2O.g(caSource[1, 1].gas.H2O.boundary.T,
            caSource[1, 1].gas.H2O.boundary.p) - cell.caCL.subregions[1, 1, 1].gas.H2O.g)
            /2 if environment.analysis "Voltage loss due to H2O removal";
        output Q.Potential Deltaw_H2=(DataH2.g(anSource[1, 1].gas.H2.boundary.T,
            anSource[1, 1].gas.H2.boundary.p) - cell.anCL.subregions[1, 1, 1].gas.H2.g)
            /2 if environment.analysis "Voltage loss due to H2 supply";
        output Q.Potential 'Deltaw_e-'=cell.caFP.subregions[1, 1, 1].graphite.
            'e-'.g_boundaries[1, Side.p] - cell.caCL.subregions[1, 1, 1].graphite.
            'e-'.g + cell.anCL.subregions[1, 1, 1].graphite.'e-'.g - cell.anFP.subregions[
            1, 1, 1].graphite.'e-'.g_boundaries[1, Side.n] if environment.analysis
          "Voltage loss due to e- transport";
        output Q.Potential 'Deltaw_H+'=cell.anCL.subregions[1, 1, 1].ionomer.
            'H+'.g - cell.caCL.subregions[1, 1, 1].ionomer.'H+'.g if
          environment.analysis "Voltage loss due to H+ transport";
        output Q.Potential Deltaw_an=cell.anCL.subregions[1, 1, 1].graphite.
            'e-Transfer'.Deltag if environment.analysis "Anode overpotential";
        output Q.Potential Deltaw_ca=-cell.caCL.subregions[1, 1, 1].graphite.
            'e-Transfer'.Deltag if environment.analysis "Cathode overpotential";

        inner Conditions.Environment environment(a={0,0,0}, analysis=true)
          "Environmental conditions"
          annotation (Placement(transformation(extent={{-30,10},{-10,30}})));

        Cell cell(
          inclLiq=inclLiq,
          anFP(subregions(gas(H2(each T_IC=testConditions.T, each p_IC=
                      testConditions.p*(1 - psi_H2O)), H2O(each p_IC=
                      testConditions.p*psi_H2O)), graphite('C+'(each T_IC=
                      testConditions.T)))),
          anGDL(subregions(gas(H2(each T_IC=testConditions.T, each p_IC=
                      testConditions.p*(1 - psi_H2O)), H2O(each p_IC=
                      testConditions.p*psi_H2O)), graphite('C+'(each T_IC=
                      testConditions.T)))),
          anCL(subregions(
              gas(H2(each T_IC=testConditions.T, each p_IC=testConditions.p*(1
                       - psi_H2O)), H2O(each p_IC=testConditions.p*psi_H2O)),
              graphite('C+'(each T_IC=testConditions.T)),
              ionomer(H2O(each lambda_IC=14*RH), 'SO3-'(each T_IC=
                      testConditions.T)))),
          PEM(subregions(ionomer(
                H2O(each lambda_IC=14*RH),
                'SO3-'(each T_IC=testConditions.T),
                'H+'(each consTransX=ConsTrans.dynamic)))),
          caCL(subregions(
              gas(
                H2O(each T_IC=testConditions.T, each p_IC=testConditions.p*
                      psi_H2O),
                O2(each p_IC=testConditions.p*(1 - psi_H2O)*testConditions.psi_O2_dry),

                N2(each p_IC=testConditions.p*(1 - psi_H2O)*(1 - testConditions.psi_O2_dry))),

              graphite('C+'(each T_IC=testConditions.T)),
              ionomer(H2O(each lambda_IC=14*RH), 'SO3-'(each T_IC=
                      testConditions.T)))),
          caGDL(subregions(gas(
                H2O(each T_IC=testConditions.T, each p_IC=testConditions.p*
                      psi_H2O),
                O2(each p_IC=testConditions.p*(1 - psi_H2O)*testConditions.psi_O2_dry),

                N2(each p_IC=testConditions.p*(1 - psi_H2O)*(1 - testConditions.psi_O2_dry))),
                graphite('C+'(each T_IC=testConditions.T)))),
          caFP(subregions(gas(
                H2O(each T_IC=testConditions.T, each p_IC=testConditions.p*
                      psi_H2O),
                O2(each p_IC=testConditions.p*(1 - psi_H2O)*testConditions.psi_O2_dry),

                N2(each p_IC=testConditions.p*(1 - psi_H2O)*(1 - testConditions.psi_O2_dry))),
                graphite('C+'(each T_IC=testConditions.T)))),
          L_y={4,4}*U.cm) "Fuel cell" annotation (__Dymola_choicesFromPackage=
              true,Placement(transformation(extent={{-10,-30},{10,-10}})));

        Real states[:](each stateSelect=StateSelect.always) = {cell.anFP.subregions[
          1, 1, 1].graphite.T,cell.anFP.subregions[1, 2, 1].graphite.T,cell.anGDL.subregions[
          1, 1, 1].graphite.T,cell.anGDL.subregions[1, 2, 1].graphite.T,cell.caFP.subregions[
          1, 2, 1].gas.H2O.phi[2],cell.caFP.subregions[1, 2, 1].gas.O2.phi[2],
          cell.anFP.subregions[1, 2, 1].gas.H2.phi[2],cell.caCL.subregions[1, 1,
          1].ORR.'e-'.reaction.Ndot,cell.caCL.subregions[1, 2, 1].ORR.'e-'.reaction.Ndot,
          cell.caCL.subregions[1, 1, 1].gas.H2O.phi[1],cell.caCL.subregions[1,
          2, 1].gas.H2O.phi[1]} if cell.n_y == 2;

        // Conditions
        Conditions.ByConnector.BoundaryBus.Single.Sink anBC[1, cell.anFP.n_z](
            each gas(
            inclH2=true,
            inclH2O=true,
            H2(
              redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current,

              materialSet(y=0),
              redeclare function afterSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Translational.velocity,

              redeclare function beforeSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Translational.velocity),

            H2O(
              redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current,

              materialSet(y=0),
              redeclare function afterSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Translational.velocity,

              redeclare function beforeSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Translational.velocity)),
            each graphite(
            'incle-'=true,
            'inclC+'=true,
            'e-'(materialSet(y=0)),
            redeclare
              FCSys.Conditions.ByConnector.ThermalDiffusive.Single.Temperature
              'C+'(set(y=testConditions.T)))) annotation (Placement(
              transformation(
              extent={{10,-10},{-10,10}},
              rotation=90,
              origin={-24,-20})));

        Conditions.ByConnector.BoundaryBus.Single.Sink anBC2[cell.anFP.n_y - 1,
          cell.anFP.n_z](each gas(
            inclH2=true,
            inclH2O=true,
            H2(
              redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current,

              materialSet(y=0),
              redeclare function afterSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Translational.velocity,

              redeclare function beforeSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Translational.velocity),

            H2O(
              redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current,

              materialSet(y=0),
              redeclare function afterSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Translational.velocity,

              redeclare function beforeSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Translational.velocity)),
            each graphite(
            'incle-'=true,
            'inclC+'=true,
            'e-'(redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current,
                materialSet(y=0)),
            redeclare
              FCSys.Conditions.ByConnector.ThermalDiffusive.Single.Temperature
              'C+'(set(y=testConditions.T)))) if cell.anFP.n_y > 1 annotation (
            Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=90,
              origin={-24,-20})));

        Conditions.ByConnector.BoundaryBus.Single.Source anSource[cell.anFP.n_x,
          cell.anFP.n_z](each gas(
            inclH2=true,
            inclH2O=true,
            H2(materialSet(y=-testConditions.I_an/2), thermalSet(y=
                    testConditions.T)),
            H2O(materialSet(y=-(testConditions.I_an/2)*(testConditions.psi_H2O_an
                    /testConditions.psi_H2)),thermalSet(y=testConditions.T))))
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={-56,-44})));

        Conditions.ByConnector.BoundaryBus.Single.Sink anSink[cell.anFP.n_x,
          cell.anFP.n_z](gas(
            each inclH2=true,
            each inclH2O=true,
            H2O(materialSet(y=fill(
                          testConditions.p,
                          cell.anFP.n_x,
                          cell.n_z) - anSink.gas.H2.p)),
            H2(materialSet(y=anSink.gas.H2O.boundary.Ndot .* cell.anFP.subregions[
                    :, cell.anFP.n_y, :].gas.H2O.v ./ cell.anFP.subregions[:,
                    cell.anFP.n_y, :].gas.H2.v), redeclare each function
                materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current)),
            each liquid(H2O(materialSet(y=testConditions.p)), inclH2O=inclLiq))
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-56,4})));
        Conditions.ByConnector.BoundaryBus.Single.Source caBC[1, cell.caFP.n_z]
          (each gas(
            inclH2O=true,
            H2O(redeclare function thermalSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
                thermalSet(y=0)),
            inclN2=true,
            N2(redeclare function thermalSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
                thermalSet(y=0)),
            inclO2=true,
            O2(redeclare function thermalSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
                thermalSet(y=0))),each graphite(
            'incle-'=true,
            'inclC+'=true,
            redeclare
              FCSys.Conditions.ByConnector.ThermalDiffusive.Single.Temperature
              'C+'(set(y=testConditions.T)),
            'e-'(
              redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.pressure,

              materialSet(y=p_ca_elec),
              thermalSet(y=testConditions.T)))) annotation (Placement(
              transformation(
              extent={{10,10},{-10,-10}},
              rotation=90,
              origin={24,-20})));

        Conditions.ByConnector.BoundaryBus.Single.Source caBC2[cell.caFP.n_y -
          1, cell.caFP.n_z](each gas(
            inclH2O=true,
            H2O(redeclare function thermalSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
                thermalSet(y=0)),
            inclN2=true,
            N2(redeclare function thermalSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
                thermalSet(y=0)),
            inclO2=true,
            O2(redeclare function thermalSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
                thermalSet(y=0))),each graphite(
            'incle-'=true,
            'inclC+'=true,
            redeclare
              FCSys.Conditions.ByConnector.ThermalDiffusive.Single.Temperature
              'C+'(set(y=testConditions.T)),
            'e-'(thermalSet(y=testConditions.T)))) if cell.caFP.n_y > 1
          annotation (Placement(transformation(
              extent={{10,10},{-10,-10}},
              rotation=90,
              origin={24,-20})));

        Conditions.ByConnector.BoundaryBus.Single.Source caSource[cell.caFP.n_x,
          cell.caFP.n_z](each gas(
            inclO2=true,
            inclN2=true,
            inclH2O=true,
            O2(materialSet(y=-testConditions.I_ca/4), thermalSet(y=
                    testConditions.T)),
            N2(materialSet(y=-(testConditions.I_ca/4)*testConditions.psi_N2/
                    testConditions.psi_O2), thermalSet(y=testConditions.T)),
            H2O(materialSet(y=-(testConditions.I_ca/4)*testConditions.psi_H2O_ca
                    /testConditions.psi_O2),thermalSet(y=testConditions.T))))
          annotation (Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=180,
              origin={56,-44})));

        Conditions.ByConnector.BoundaryBus.Single.Sink caSink[cell.caFP.n_x,
          cell.caFP.n_z](gas(
            each inclO2=true,
            each inclN2=true,
            each inclH2O=true,
            H2O(materialSet(y=fill(
                          testConditions.p,
                          cell.caFP.n_x,
                          cell.n_z) - caSink.gas.N2.p - caSink.gas.O2.p)),
            N2(materialSet(y=caSink.gas.H2O.boundary.Ndot .* cell.caFP.subregions[
                    :, cell.caFP.n_y, :].gas.H2O.v ./ cell.caFP.subregions[:,
                    cell.caFP.n_y, :].gas.N2.v), redeclare function
                materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current),

            O2(materialSet(y=caSink.gas.H2O.boundary.Ndot .* cell.caFP.subregions[
                    :, cell.caFP.n_y, :].gas.H2O.v ./ cell.caFP.subregions[:,
                    cell.caFP.n_y, :].gas.O2.v), redeclare function
                materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current)),
            each liquid(H2O(materialSet(y=testConditions.p)), inclH2O=inclLiq))
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={56,4})));

        TestConditions testConditions(
          final zI=-sum(caBC.graphite.'e-'.boundary.Ndot),
          electricalSpec=FCSys.Assemblies.Cells.Examples.Enumerations.ElectricalSpec.voltage,

          redeclare replaceable Modelica.Blocks.Sources.Ramp electricalSet(
            duration=600,
            startTime=50,
            height=-0.6*U.V,
            offset=1.2*U.V)) annotation (Dialog, Placement(transformation(
                extent={{10,10},{30,30}})));

      protected
        Q.Pressure p_ca_elec "Electronic pressure on the cathode";
        Conditions.Router anRouter[cell.anFP.n_x, cell.n_z] annotation (Dialog,
            Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=0,
              origin={-48,-20})));
        Conditions.Router caRouter[cell.anFP.n_x, cell.n_z] annotation (Dialog,
            Placement(transformation(
              extent={{10,-10},{-10,10}},
              rotation=0,
              origin={48,-20})));

      equation
        testConditions.w = 'Datae-'.g(testConditions.T, 0) - 'Datae-'.g(
          testConditions.T, p_ca_elec);

        connect(cell.ca[1, :], caBC[1, :].boundary) annotation (Line(
            points={{10,-20},{20,-20}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.an[1, :], anBC[1, :].boundary) annotation (Line(
            points={{-10,-20},{-20,-20}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.ca[2:end, :], caBC2.boundary) annotation (Line(
            points={{10,-20},{20,-20}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.an[2:end, :], anBC2.boundary) annotation (Line(
            points={{-10,-20},{-20,-20}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(anRouter.positive2, anSink.boundary) annotation (Line(
            points={{-56,-16},{-56,0}},
            color={240,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(anRouter.positive1, anSource.boundary) annotation (Line(
            points={{-56,-24},{-56,-40}},
            color={240,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.anPositive, anRouter.negative2) annotation (Line(
            points={{-4,-10},{-4,0},{-40,0},{-40,-16}},
            color={240,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.anNegative, anRouter.negative1) annotation (Line(
            points={{-4,-30},{-4,-40},{-40,-40},{-40,-24}},
            color={240,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(caSink.boundary, caRouter.negative2) annotation (Line(
            points={{56,0},{56,-16}},
            color={0,0,240},
            thickness=0.5,
            smooth=Smooth.None));
        connect(caSource.boundary, caRouter.negative1) annotation (Line(
            points={{56,-40},{56,-24}},
            color={0,0,240},
            thickness=0.5,
            smooth=Smooth.None));
        connect(caRouter.positive2, cell.caPositive) annotation (Line(
            points={{40,-16},{40,0},{4,0},{4,-10}},
            color={0,0,240},
            thickness=0.5,
            smooth=Smooth.None));
        connect(caRouter.positive1, cell.caNegative) annotation (Line(
            points={{40,-24},{40,-40},{4,-40},{4,-30}},
            color={0,0,240},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (
          Documentation(info="<html>
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
    </html>"),
          Commands(file=
                "Resources/Scripts/Dymola/Assemblies.Cells.Examples.TestStand.mos"
              "Assemblies.Cells.Examples.TestStand.mos", file=
                "Resources/Scripts/Dymola/Assemblies.Cells.Examples.TestStand-states.mos"
              "Assemblies.Cells.Examples.TestStand-states.mos"),
          experiment(
            StopTime=650,
            Tolerance=1e-006,
            __Dymola_Algorithm="Dassl"),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-80,-60},
                  {80,40}}), graphics),
          __Dymola_experimentSetupOutput);
      end TestStandVoltage;


    end Examples;
    extends Modelica.Icons.Package;

    model Cell "Single-cell PEMFC"
      extends FCSys.Icons.Cell;

      // Geometric parameters
      parameter Q.Length L_y[:]={8}*U.cm "Lengths in the y direction"
        annotation (Dialog(group="Geometry", __Dymola_label=
              "<html><i>L</i><sub>y</sub></html>"));
      parameter Q.Length L_z[:]={6.25}*U.cm "Lengths in the z direction"
        annotation (Dialog(group="Geometry", __Dymola_label=
              "<html><i>L</i><sub>z</sub></html>"));
      final parameter Integer n_y=size(L_y, 1)
        "Number of subregions in the y direction";
      final parameter Integer n_z=size(L_z, 1)
        "Number of subregions in the z direction";

      // Assumptions
      parameter Boolean inclLiq=true
        "<html>Include liquid H<sub>2</sub>O</html>" annotation (Dialog(tab=
              "Assumptions", compact=true), choices(__Dymola_checkBox=true));
      parameter Boolean inclN2=true "<html>Include N<sub>2</sub></html>"
        annotation (Dialog(tab="Assumptions", compact=true), choices(
            __Dymola_checkBox=true));

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
        subregions(liquid(each inclH2O=inclLiq))) "Anode flow plate"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-80,-20},{-60,0}})));
      replaceable Regions.AnGDLs.AnGDL anGDL(
        final L_y=L_y,
        final L_z=L_z,
        subregions(liquid(each inclH2O=inclLiq))) "Anode gas diffusion layer"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-60,-20},{-40,0}})));
      replaceable Regions.AnCLs.AnCL anCL(
        final L_y=L_y,
        final L_z=L_z,
        subregions(liquid(each inclH2O=inclLiq))) "Anode catalyst layer"
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
        subregions(gas(each inclN2=inclN2), liquid(each inclH2O=inclLiq)))
        "Cathode catalyst layer" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{0,-20},{20,0}})));
      replaceable Regions.CaGDLs.CaGDL caGDL(
        final L_y=L_y,
        final L_z=L_z,
        subregions(gas(each inclN2=inclN2), liquid(each inclH2O=inclLiq)))
        "Cathode gas diffusion layer" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{20,-20},{40,0}})));
      replaceable Regions.CaFPs.CaFP caFP(
        final L_y=L_y,
        final L_z=L_z,
        subregions(gas(each inclN2=inclN2), liquid(each inclH2O=inclLiq)))
        "Cathode flow plate" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{40,-20},{60,0}})));

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
      // Note:  This model must be marked as structurally incomplete to pass the check
      // in Dymola 2014.
      annotation (
        structurallyIncomplete=true,
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
      import Modelica.Constants.inf;
      extends FCSys.Icons.Cell;

      // Geometric parameters
      parameter Q.Length L_y[:]={8}*U.cm "Lengths in the y direction"
        annotation (Dialog(group="Geometry", __Dymola_label=
              "<html><i>L</i><sub>y</sub></html>"));
      parameter Q.Length L_z[:]={6.25}*U.cm "Lengths in the z direction"
        annotation (Dialog(group="Geometry", __Dymola_label=
              "<html><i>L</i><sub>z</sub></html>"));
      final parameter Integer n_y=size(L_y, 1)
        "Number of subregions in the y direction";
      final parameter Integer n_z=size(L_z, 1)
        "Number of subregions in the z direction";

      // Assumptions
      parameter Boolean inclLiq=false
        "<html>Include liquid H<sub>2</sub>O</html>" annotation (Dialog(tab=
              "Assumptions", compact=true), choices(__Dymola_checkBox=true));
      parameter Boolean inclN2=true "<html>Include N<sub>2</sub></html>"
        annotation (Dialog(tab="Assumptions", compact=true), choices(
            __Dymola_checkBox=true));

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
        subregions(
          common(each k_Q=0),
          gas(common(each k_Q=inf),H2(each initEnergy=Init.none,T(each
                  stateSelect=StateSelect.default))),
          liquid(each inclH2O=inclLiq, H2O(each initEnergy=Init.none,T(each
                  stateSelect=StateSelect.default))))) "Anode flow plate"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-60,-20},{-40,0}})));

      FCSys.Regions.AnCLs.AnCGDL anCGDL(
        final L_y=L_y,
        final L_z=L_z,
        subregions(
          common(each k_Q=0),
          gas(common(each k_Q=inf),H2(each initEnergy=Init.none,T(each
                  stateSelect=StateSelect.default))),
          ionomer('SO3-'(T(each fixed=false,each stateSelect=StateSelect.default))),

          liquid(each inclH2O=inclLiq, H2O(each initEnergy=Init.none,T(each
                  stateSelect=StateSelect.default)))))
        "Anode catalyst and gas diffusion layer" annotation (Dialog(group=
              "Layers"), Placement(transformation(extent={{-40,-20},{-20,0}})));

      replaceable Regions.PEMs.PEM PEM(final L_y=L_y, final L_z=L_z)
        "Proton exchange membrane" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-20,-20},{0,0}})));
      FCSys.Regions.CaCLs.CaCGDL caCGDL(
        final L_y=L_y,
        final L_z=L_z,
        subregions(
          common(each k_Q=0),
          gas(
            common(each k_Q=inf),
            each inclN2=inclN2,
            H2O(each initEnergy=Init.none,T(each stateSelect=StateSelect.default))),

          ionomer('SO3-'(T(each fixed=false,each stateSelect=StateSelect.default))),

          liquid(each inclH2O=inclLiq,H2O(each initEnergy=Init.none,T(each
                  stateSelect=StateSelect.default)))))
        "Cathode catalyst and gas diffusion layer" annotation (Dialog(group=
              "Layers"), Placement(transformation(extent={{0,-20},{20,0}})));

      replaceable Regions.CaFPs.CaFP caFP(
        final L_y=L_y,
        final L_z=L_z,
        subregions(
          common(each k_Q=0),
          gas(
            common(each k_Q=inf),
            each inclN2=inclN2,
            H2O(each initEnergy=Init.none,T(each stateSelect=StateSelect.default))),

          liquid(each inclH2O=inclLiq,H2O(each initEnergy=Init.none,T(each
                  stateSelect=StateSelect.default))))) "Cathode flow plate"
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
