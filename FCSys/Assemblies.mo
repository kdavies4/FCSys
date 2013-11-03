within FCSys;
package Assemblies "Combinations of regions (e.g., cells)"
  extends Modelica.Icons.Package;

  package Cells "Single-cell PEMFC models"
    extends Modelica.Icons.Package;

    package Examples "Examples"

      extends Modelica.Icons.ExamplesPackage;

      model TestConditions "Fuel cell test conditions"
        extends FCSys.Icons.Record;
        import FCSys.Characteristics.H2O.p_sat;

        // Note:  This isn't a record because it contains time-varying variables.

        /* **
                 params=dict(comp=['"O2"'],
                             anStoich=[1.5, 1.1, 2],
                             caStoich=[9.5, 7.5, 12.5],
                             anRH=[0.8, 0.6, 1],
                             caRH=[0.5, 0.3, 0.7],
                             T_degC=[60, 40, 80],
                             p_kPag=[48.3, 0, 202.7]),
                             */
        // Basic conditions
        parameter Q.TemperatureAbsolute T=333.15*U.K
          "Temperature (IC, end plate BC, and inlets)"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        parameter Q.PressureAbsolute p(nominal=U.atm) = U.from_kPag(48.3)
          "Pressure (IC and outlets)"
          annotation (Dialog(__Dymola_label="<html><i>p</i></html>"));
        final parameter Q.NumberAbsolute psi_sat=Characteristics.H2O.p_sat(T)/p
          "Mole fraction of H2O at saturation";

        // Electrical
        parameter Enumerations.ElectricalSpec electricalSpec=ElectricalSpec.voltage
          "Method of specification" annotation (Dialog(tab="Electrical"));
        replaceable Modelica.Blocks.Sources.Ramp electricalSet(
          offset=U.mA,
          height=100*U.A,
          duration=600) constrainedby Modelica.Blocks.Interfaces.SO "Setpoint"
          annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(tab="Electrical"),
          Placement(transformation(extent={{-30,-10},{-10,10}})));
        Q.Potential w "Voltage";
        Q.ResistanceElectrical R "Resistance";
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
        parameter Q.NumberAbsolute psi_O2_dry(
          final max=1,
          displayUnit="%") = environment.psi_O2_dry
          "<html>Dry-gas concentration of O<sub>2</sub> (at inlet)</html>"
          annotation (Dialog(tab="Cathode", __Dymola_label=
                "<html>&psi;<sub>O2 dry</sub></html>"));
        final parameter Q.NumberAbsolute psi_H2O_ca=caRH*psi_sat
          "Mole fraction of H2O at the cathode inlet";
        final parameter Q.NumberAbsolute psi_O2=psi_O2_dry*(1 - psi_H2O_ca)
          "Mole fraction of O2 at the cathode inlet";
        final parameter Q.NumberAbsolute psi_N2=(1 - psi_O2_dry)*(1 -
            psi_H2O_ca) "Mole fraction of N2 at the cathode inlet";
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
        input Q.Current zI "Current"
          annotation (Dialog(tab="Measured conditions"));

      protected
        inner Conditions.Environment environment
          annotation (Placement(transformation(extent={{-10,30},{10,50}})));

      equation
        // Electrical
        // ----------
        // Relations
        zI = R*w "**";
        P = w*zI;
        //
        // Setpoint
        if electricalSpec == ElectricalSpec.current then
          zI = electricalSet.y;
        elseif electricalSpec == ElectricalSpec.voltage then
          w = electricalSet.y;
        elseif electricalSpec == ElectricalSpec.resistance then
          R = electricalSet.y;
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
        annotation (Documentation(info="<html>    <p><i>Equivalent current</i> is the rate of supply of a reactant required to support the
    given current
    assuming the reactant is entirely consumed (complete utilization).</p></html>"));
      end TestConditions;

      model TestStandEIS
        "Test stand to perform electrochemical impedance spectroscopy"
        extends TestStand(testConditions(final electricalSpec=ElectricalSpec.current,
              final u_electrical=zI_large + zI_small_A*U.A));

        parameter Q.Current zI_large=U.A "Large-signal current" annotation (
            Dialog(__Dymola_label="<html><i>zJ</i><sub>large</sub></html>"));
        Connectors.RealInput zI_small_A "Small-signal current in amperes";
        Connectors.RealOutput w_V=w/U.V "Cell potential in volts";

        annotation (Documentation(info="<html><p>This model modulates the electrical current applied to the cell 
    according to an input.
    The current is the sum of a steady-state large-signal current and a small-signal 
    current introduced via the input <i>zI</i><sub>small A</sub>.</p>
       
    <p>For more information, please see the documentation in the
    <a href=\"modelica://FCSys.Conditions.TestStands.TestStand\">test stand</a> model.</p></html>"));

      end TestStandEIS;

      package Enumerations "Choices of options"

        extends Modelica.Icons.BasesPackage;

        type ElectricalSpec = enumeration(
            current "Current",
            voltage "Voltage",
            resistance "Resistance",
            power "Power") "Ways to specify the electrical load";
        type FlowSpec = enumeration(
            stoich "Stoichiometric rate",
            current "Equivalent current",
            pressure "Inlet pressure") "Ways to specify the anode flow rate";

      end Enumerations;

      model TestStand "Simulate the fuel cell with prescribed conditions"
        import 'Datae-' = FCSys.Characteristics.'e-'.Graphite;
        extends Modelica.Icons.Example;

        parameter Boolean inclLiq=false "Include liquid H2O";
        output Q.Number J_Apercm2=testConditions.zI*U.cm^2/(cell.caFP.A[Axis.x]
            *U.A) "Electrical current density, in A/cm2";
        output Q.Potential w=cell.anFP.subregions[1, 1, 1].graphite.'e-'.g_boundaries[
            1, Side.n] - cell.caFP.subregions[end, 1, 1].graphite.'e-'.g_boundaries[
            1, Side.p] if environment.analysis "Electrical potential";

        inner Conditions.Environment environment(analysis=true, a={0,0,0})
          "Environmental conditions"
          annotation (Placement(transformation(extent={{-30,10},{-10,30}})));

        Cell cell(
          inclLiq=inclLiq,
          anFP(subregions(gas(H2O(each T_IC=testConditions.T,each p_IC=
                      testConditions.p*testConditions.psi_H2O_an), H2(each p_IC
                    =testConditions.p*testConditions.psi_H2)), graphite('C+'(
                    each T_IC=testConditions.T)))),
          anGDL(subregions(gas(H2O(each T_IC=testConditions.T,each p_IC=
                      testConditions.p*testConditions.psi_H2O_an), H2(each p_IC
                    =testConditions.p*testConditions.psi_H2)), graphite('C+'(
                    each T_IC=testConditions.T)))),
          anCL(subregions(
              gas(H2O(each T_IC=testConditions.T,each p_IC=testConditions.p*
                      testConditions.psi_H2O_an),H2(each p_IC=testConditions.p*
                      testConditions.psi_H2)),
              graphite('C+'(each T_IC=testConditions.T)),
              ionomer('SO3-'(each T_IC=testConditions.T)))),
          PEM(subregions(ionomer('SO3-'(each T_IC=testConditions.T)))),
          caCL(subregions(
              gas(
                H2O(each T_IC=testConditions.T, each p_IC=testConditions.p*
                      testConditions.psi_H2O_ca),
                O2(each p_IC=testConditions.p*testConditions.psi_O2),
                N2(each p_IC=testConditions.p*testConditions.psi_N2)),
              graphite('C+'(each T_IC=testConditions.T)),
              ionomer('SO3-'(each T_IC=testConditions.T)))),
          caGDL(subregions(gas(
                H2O(each T_IC=testConditions.T, each p_IC=testConditions.p*
                      testConditions.psi_H2O_ca),
                O2(each p_IC=testConditions.p*testConditions.psi_O2),
                N2(each p_IC=testConditions.p*testConditions.psi_N2)), graphite(
                  'C+'(each T_IC=testConditions.T)))),
          caFP(subregions(gas(
                H2O(each T_IC=testConditions.T, each p_IC=testConditions.p*
                      testConditions.psi_H2O_ca),
                O2(each p_IC=testConditions.p*testConditions.psi_O2),
                N2(each p_IC=testConditions.p*testConditions.psi_N2)), graphite(
                  'C+'(each T_IC=testConditions.T))))) "Fuel cell" annotation (
            __Dymola_choicesFromPackage=true, Placement(transformation(extent={
                  {-10,-30},{10,-10}})));

        // Conditions
        Conditions.ByConnector.BoundaryBus.Single.Sink anBC[cell.anFP.n_y, cell.anFP.n_z]
          (each graphite(
            'incle-'=true,
            'inclC+'=true,
            'e-'(materialSet(y=U.bar)),
            redeclare
              FCSys.Conditions.ByConnector.ThermoDiffusive.Single.Temperature
              'C+'(source(y=testConditions.T)))) annotation (Placement(
              transformation(
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
                          cell.anFP.n_z) - anSink.gas.H2.p)),
            H2(materialSet(y=anSink.gas.H2O.boundary.Ndot .* cell.anFP.subregions[
                    :, cell.anFP.n_y, :].gas.H2O.v ./ cell.anFP.subregions[:,
                    cell.anFP.n_y, :].gas.H2.v), redeclare each function
                materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.current)),
            each liquid(H2O(materialSet(y=testConditions.p_an_out)), inclH2O=
                inclLiq)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-56,4})));
        Conditions.ByConnector.BoundaryBus.Single.Source caBC[cell.caFP.n_y,
          cell.caFP.n_z](each graphite(
            'incle-'=true,
            'inclC+'=true,
            redeclare
              FCSys.Conditions.ByConnector.ThermoDiffusive.Single.Temperature
              'C+'(source(y=testConditions.T)),
            'e-'(
              redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Boundary.Single.Material.pressure,

              materialSet(y=p_ca_elec),
              thermalSet(y=testConditions.T)))) annotation (Placement(
              transformation(
              extent={{10,10},{-10,-10}},
              rotation=90,
              origin={24,-20})));

        Conditions.ByConnector.BoundaryBus.Single.Source caSource[cell.caFP.n_x,
          cell.caFP.n_z](each gas(
            inclO2=true,
            inclN2=true,
            inclH2O=true,
            O2(materialSet(y=-testConditions.I_an/4), thermalSet(y=
                    testConditions.T)),
            N2(materialSet(y=-(testConditions.I_an/4)*testConditions.psi_N2/
                    testConditions.psi_O2), thermalSet(y=testConditions.T)),
            H2O(materialSet(y=-(testConditions.I_an/4)*testConditions.psi_H2O_ca
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
                          cell.caFP.n_z) - caSink.gas.N2.p - caSink.gas.O2.p)),

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
            each liquid(H2O(materialSet(y=testConditions.p_ca_out)), inclH2O=
                inclLiq)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={56,4})));

        TestConditions testConditions(
          zI=-sum(caBC.graphite.'e-'.boundary.Ndot),
          electricalSpec=FCSys.Assemblies.Cells.Examples.Enumerations.ElectricalSpec.resistance,

          electricalSet(
            startTime=100,
            height=1e-4*U.S,
            offset=1e-10*U.S)) "Test conditions" annotation (Dialog, Placement(
              transformation(extent={{10,10},{30,30}})));

      protected
        Q.Velocity phi_states[:](
          each stateSelect=StateSelect.always,
          each start=0,
          each fixed=true) = {cell.anFP.subregions[1, 1, 1].gas.H2.phi[1],cell.anFP.subregions[
          1, 1, 1].gas.H2O.phi[1],cell.caFP.subregions[1, 1, 1].gas.H2O.phi[1],
          cell.caFP.subregions[1, 1, 1].gas.N2.phi[1],cell.caFP.subregions[1, 1,
          1].gas.O2.phi[1],cell.caGDL.subregions[1, 1, 1].gas.O2.phi[1],cell.PEM.subregions[
          1, 1, 1].ionomer.H2O.phi[1],cell.caGDL.subregions[1, 1, 1].gas.H2O.phi[
          1],cell.caGDL.subregions[1, 1, 1].gas.N2.phi[1],cell.anGDL.subregions[
          1, 1, 1].gas.H2O.phi[1],cell.anGDL.subregions[1, 1, 1].gas.H2.phi[1],
          cell.anCL.subregions[1, 1, 1].ionomer.H2O.phi[1]}
          "Velocities for state selection";

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
          2, 2].Ndot,cell.caFP.subregions[1, 1, 1].gas.H2O.boundaries[2, 2].Ndot,
          cell.caCL.subregions[1, 1, 1].ORR.'e-'.reaction.Ndot}
          "Currents for state selection";

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

        testConditions.w = 'Datae-'.g(testConditions.T, U.bar) - 'Datae-'.g(
          testConditions.T, p_ca_elec);

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
            StopTime=600,
            Tolerance=1e-005,
            __Dymola_Algorithm="Dassl"),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-80,-60},
                  {80,40}}), graphics),
          __Dymola_experimentSetupOutput);
      end TestStand;
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
            initialScale=0.1), graphics={Line(
                  points={{-40,-58},{-40,-100}},
                  color={240,0,0},
                  visible=inclY,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{-8,-1},{28,-1}},
                  color={0,0,240},
                  visible=inclX,
                  thickness=0.5,
                  origin={39,-92},
                  rotation=90),Line(
                  points={{-40,100},{-40,60}},
                  color={240,0,0},
                  visible=inclY,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{-66,0},{-100,0}},
                  color={127,127,127},
                  visible=inclX,
                  thickness=0.5),Line(
                  points={{-8,-1},{44,-1}},
                  color={0,0,240},
                  visible=inclX,
                  thickness=0.5,
                  origin={39,56},
                  rotation=90),Line(
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
            initialScale=0.1), graphics={Line(
                  points={{-40,100},{-40,60}},
                  color={240,0,0},
                  visible=inclY,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{-8,-1},{44,-1}},
                  color={0,0,240},
                  visible=inclX,
                  thickness=0.5,
                  origin={39,56},
                  rotation=90),Line(
                  points={{100,0},{56,0}},
                  color={127,127,127},
                  visible=inclX,
                  thickness=0.5),Line(
                  points={{-8,-1},{28,-1}},
                  color={0,0,240},
                  visible=inclX,
                  thickness=0.5,
                  origin={39,-92},
                  rotation=90),Line(
                  points={{-40,-58},{-40,-100}},
                  color={240,0,0},
                  visible=inclY,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
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
