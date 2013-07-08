within FCSys;
package Assemblies "Combinations of regions (e.g., cells)"
  extends Modelica.Icons.Package;

  package Cells "Single-cell PEMFC models"
    extends Modelica.Icons.Package;
    package Examples "Examples"

      model Cell "<html>Isolated <code>Cell</code> model</html>"
        extends Modelica.Icons.Example;
        extends Modelica.Icons.UnderConstruction;
        inner FCSys.Conditions.Environment environment(
          analysis=false,
          p=149.6*U.kPa,
          T=333.15*U.K)
          annotation (Placement(transformation(extent={{20,20},{40,40}})));
        replaceable Cells.Cell cell annotation (__Dymola_choicesFromPackage=
              true, Placement(transformation(extent={{-10,-10},{10,10}})));
        annotation (experiment(StopTime=1e-24, Tolerance=1e-06), Commands(file=
                "Resources/Scripts/Dymola/Assemblies.Cells.Examples.Cell.mos"
              "Assemblies.Cells.Examples.Cell.mos"));

      end Cell;
      extends Modelica.Icons.ExamplesPackage;
      model CellProfile
        "Apply boundary conditions to a cell according to a test profile"

        extends FCSys.Conditions.TestStands.TestProfile(anEnd(each graphite(
                'incle-'=true, 'e-'(redeclare Modelica.Blocks.Sources.Ramp
                  materialSpec(height=10000*U.A, duration=500)))), caEnd(each
              graphite('incle-'=true, 'e-'(redeclare
                  Modelica.Blocks.Sources.Ramp materialSpec(height=-10000*U.A,
                    duration=500)))));
        extends Modelica.Icons.UnderConstruction;
        replaceable Cell cell annotation (__Dymola_choicesFromPackage=true,
            Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        connect(cell.anFPX, anEnd.face) annotation (Line(
            points={{-10,6.10623e-16},{-14,6.10623e-16},{-14,3.65701e-16},{-26,
                3.65701e-16}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));

        connect(cell.caFPX, caEnd.face) annotation (Line(
            points={{10,6.10623e-16},{14,6.10623e-16},{14,1.23436e-15},{26,
                1.23436e-15}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));

        connect(anSource.face, cell.anFPNegY) annotation (Line(
            points={{-20,-26},{-20,-20.5},{-4,-20.5},{-4,-10}},
            color={253,52,56},
            thickness=0.5,
            smooth=Smooth.None));

        connect(cell.anFPPosY, anSink.face) annotation (Line(
            points={{-4,10},{-4,20},{-20,20},{-20,26}},
            color={253,52,56},
            thickness=0.5,
            smooth=Smooth.None));
        connect(caSink.face, cell.caFPPosY) annotation (Line(
            points={{20,26},{20,20},{4,20},{4,10}},
            color={0,0,240},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.caFPNegY, caSource.face) annotation (Line(
            points={{4,-10},{4,-20},{20,-20},{20,-26}},
            color={0,0,240},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (
          experiment(StopTime=100, Tolerance=1e-06),
          Commands(file=
                "Resources/Scripts/Dymola/Assemblies.Cells.Examples.CellProfile.mos"
              "Assemblies.Cells.Examples.CellProfile.mos"),
          experiment(StopTime=600, Tolerance=1e-08));
      end CellProfile;

      model Polarization "Run a cell polarization"
        extends CellProfile;
        extends Modelica.Icons.UnderConstruction;
        annotation (Commands(file=
                "Resources/Scripts/Dymola/Assemblies.Cells.Examples.Polarization.mos"
              "Assemblies.Cells.Examples.Polarization.mos"));

      end Polarization;

      model CellProfileIO
        "Apply Conditions to a cell according to a test profile, with inputs and outputs"
        extends CellProfile(testStand(final inclIO=true));
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;
        extends Modelica.Icons.UnderConstruction;
        Connectors.RealInputBus u "Input bus" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=315,
              origin={-30,30}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-2,100})));
        Connectors.RealOutputBus y "Output bus" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=315,
              origin={30,-30}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-100})));

      equation
        connect(y, testStand.y) annotation (Line(
            points={{30,-30},{16,-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        connect(u, testStand.u) annotation (Line(
            points={{-30,30},{-16,16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (experiment(StopTime=600, Tolerance=1e-08), Commands(file=
                "Resources/Scripts/Dymola/Assemblies.Examples.CellPolarizationstoich.mos"
              "Assemblies.Examples.CellPolarizationstoich.mos"));
      end CellProfileIO;

      model CellModelica
        "<html>Cell interfaced to components from the <a href=\"modelica://Modelica\">Modelica</a> package</html>"
        extends Modelica.Icons.Example;
        extends Modelica.Icons.UnderConstruction;
        Cell cell(anFP(redeclare FCSys.Subregions.Subregion subregions(
              each final inclX=true,
              each inclY=true,
              each graphite('incle-'=true, 'e-'(perfectMaterialDiff={{{{true,
                      false}}}})),
              each gas(inclH2=true, inclH2O=true))))
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

        inner FCSys.Conditions.Environment environment(analysis=false)
          annotation (Placement(transformation(extent={{40,60},{60,80}})));
        Conditions.Adapters.Phases.Graphite caModelicaAdapt(A=cell.L_y[1]*cell.L_z[
              1])
          annotation (Placement(transformation(extent={{20,-10},{40,10}})));
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

      model PolarizationPlaceholder "**temp"
        extends Modelica.Icons.Example;

        /*
                 params=dict(comp=['"O2"'],
                             anStoich=[1.5, 1.1, 2],
                             caStoich=[9.5, 7.5, 12.5],
                             anRH=[0.8, 0.6, 1],
                             caRH=[0.5, 0.3, 0.7],
                             T_degC=[60, 40, 80],
                             p_kPag=[48.3, 0, 202.7]),
                             */
      end PolarizationPlaceholder;

      model EISPlaceholder
        "Placeholder model for electro-impedance spectroscopy"
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;

        parameter Modelica.SIunits.Current zI_large_A=100
          "Large-signal current in amperes";
        Modelica.Electrical.Analog.Basic.Resistor resistor2(R=0.1) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={10,-10})));
        Modelica.Electrical.Analog.Basic.Capacitor capacitor(C=1e-3)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={40,-10})));
        Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V=1)
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-80,-10})));
        Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor
          annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=270,
              origin={70,-10})));
        Modelica.Electrical.Analog.Sources.ConstantCurrent constantCurrent(I=
              zI_large_A)
          annotation (Placement(transformation(extent={{-60,20},{-40,0}})));
        Modelica.Electrical.Analog.Sources.SignalCurrent signalCurrent
          annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
        Connectors.RealInput zI_small_A "Small-signal cell current in amperes"
          annotation (Placement(transformation(extent={{-110,40},{-90,60}}),
              iconTransformation(extent={{-120,-10},{-100,10}})));
        Connectors.RealOutput w_V "Cell potential in volts" annotation (
            Placement(transformation(extent={{90,-20},{110,0}}),
              iconTransformation(extent={{100,-10},{120,10}})));
        Modelica.Electrical.Analog.Basic.Resistor resistor1(R=0.1) annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-10,10})));
        Modelica.Electrical.Analog.Basic.Ground ground
          annotation (Placement(transformation(extent={{60,-60},{80,-40}})));
      equation
        connect(zI_small_A, signalCurrent.i) annotation (Line(
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


    end Examples;

    model Cell "Single-cell PEMFC"
      import FCSys.BaseClasses.Utilities.average;
      extends FCSys.BaseClasses.Icons.Cell;

      // Geometric parameters
      parameter Q.Length L_y[:]={U.m}
        "<html>Lengths along the channel (L<sub>y</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      parameter Q.Length L_z[:]={5*U.mm}
        "<html>Lengths across the channel (L<sub>z</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      final parameter Integer n_y=size(L_y, 1)
        "Number of subregions along the channel";
      final parameter Integer n_z=size(L_z, 1)
        "Number of subregions across the channel";

      // Essential analysis variables
      final parameter Q.Area A=anFP.A[Axis.x] "Cross-sectional area";
      output Q.Power Wdot(stateSelect=StateSelect.never) = -sum(anFP.subregions[
        1, :, :].graphite.'e-'.faces[1, Side.n].phi[1] .* anFP.subregions[1, :,
        :].graphite.'e-'.faces[1, Side.n].mPhidot[1] + caFP.subregions[caFP.n_x,
        :, :].graphite.'e-'.faces[1, Side.p].phi[1] .* caFP.subregions[caFP.n_x,
        :, :].graphite.'e-'.faces[1, Side.p].mPhidot[1])
        "Electrical power output";
      output Q.Potential w(stateSelect=StateSelect.never) = sum(anFP.subregions[
        1, :, :].graphite.'e-'.faces[1, Side.n].mPhidot[1] ./ anFP.subregions[1,
        :, :].graphite.'e-'.faces[1, Side.n].rho + caFP.subregions[caFP.n_x, :,
        :].graphite.'e-'.faces[1, Side.p].mPhidot[1] ./ caFP.subregions[caFP.n_x,
        :, :].graphite.'e-'.faces[1, Side.p].rho)/A "Electrical potential";
      output Q.Current zI(stateSelect=StateSelect.never) = Wdot/w
        "Electrical current";

      // Auxiliary variables (for analysis)
      final parameter Q.Length L[Axis]={anFP.L[Axis.x] + anGDL.L[Axis.x] + anCL.L[
          Axis.x] + PEM.L[Axis.x] + caCL.L[Axis.x] + caGDL.L[Axis.x] + caFP.L[
          Axis.x],anFP.L[Axis.y],anFP.L[Axis.z]} if environment.analysis
        "Total lengths along the x, y, and z axes";
      final parameter Q.Volume V=product(L) if environment.analysis "Volume";
      output Q.Power Wdot_yz[n_y, n_z](each stateSelect=StateSelect.never) = -
        anFP.subregions[1, :, :].graphite.'e-'.faces[1, Side.n].phi[1] .* anFP.subregions[
        1, :, :].graphite.'e-'.faces[1, Side.n].mPhidot[1] - caFP.subregions[
        caFP.n_x, :, :].graphite.'e-'.faces[1, Side.p].phi[1] .* caFP.subregions[
        caFP.n_x, :, :].graphite.'e-'.faces[1, Side.p].mPhidot[1] if
        environment.analysis "Electrical power of the segments (x axis)";
      output Q.CurrentAreic zJ_yz[n_y, n_z](each stateSelect=StateSelect.never)
         = -anFP.subregions[1, :, :].graphite.'e-'.faces[1, Side.n].phi[1] .*
        anFP.subregions[1, :, :].graphite.'e-'.faces[1, Side.n].rho if
        environment.analysis
        "Areic electrical current of the segments (x axis)";
      output Q.Current zI_yz[n_y, n_z](each stateSelect=StateSelect.never) =
        zJ_yz .* anFP.subregions[1, :, :].A[Axis.x] if environment.analysis
        "Electrical current of the segments (x axis)";
      output Q.Current Ndot_H2(stateSelect=StateSelect.never) = sum(anFP.subregions[
        :, 1, :].gas.H2.faces[2, Side.n].Ndot + anFP.subregions[:, n_y, :].gas.H2.faces[
        2, Side.p].Ndot) if environment.analysis
        "Rate of hydrogen intake (y axis)";
      output Q.Current Ndot_H2O(stateSelect=StateSelect.never) = sum(anFP.subregions[
        :, 1, :].gas.H2O.faces[2, Side.n].Ndot + anFP.subregions[:, n_y, :].gas.H2O.faces[
        2, Side.p].Ndot) + sum(caFP.subregions[:, 1, :].gas.H2O.faces[1, Side.n].Ndot
         + caFP.subregions[:, n_y, :].gas.H2O.faces[1, Side.p].Ndot) if
        environment.analysis "Rate of water intake (y axis)";
      output Q.Current Ndot_O2(stateSelect=StateSelect.never) = sum(caFP.subregions[
        :, 1, :].gas.O2.faces[1, Side.n].Ndot + caFP.subregions[:, n_y, :].gas.O2.faces[
        1, Side.p].Ndot) if environment.analysis
        "Rate of oxygen intake (y axis)";

      Connectors.FaceBus anFPX[n_y, n_z] "Interface with the anode end plate"
        annotation (Placement(transformation(extent={{-90,-10},{-70,10}},
              rotation=0), iconTransformation(extent={{-110,-10},{-90,10}})));
      Connectors.FaceBus caFPX[n_y, n_z] "Interface with the cathode end plate"
        annotation (Placement(transformation(extent={{70,-10},{90,10}},
              rotation=0), iconTransformation(extent={{90,-10},{110,10}})));
      Connectors.FaceBus anFPNegY[anFP.n_x, n_z] "Negative anode fluid port"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-60,-20}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-40,-100})));
      Connectors.FaceBus caFPNegY[caFP.n_x, n_z] "Negative cathode fluid port"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={60,-20}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={40,-100})));
      Connectors.FaceBus anFPPosY[anFP.n_x, n_z] "Positive anode fluid port"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-60,20}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-40,100})));
      Connectors.FaceBus caFPPosY[caFP.n_x, n_z] "Positive cathode fluid port"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={60,20}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={40,100})));

      replaceable Regions.AnFPs.AnFP anFP(final L_y=L_y, final L_z=L_z)
        "Anode flow plate" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-70,-10},{-50,10}})));
      replaceable Regions.AnGDLs.AnGDL anGDL(final L_y=L_y, final L_z=L_z)
        "Anode gas diffusion layer" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-50,-10},{-30,10}})));
      replaceable Regions.AnCLs.AnCL anCL(final L_y=L_y, final L_z=L_z)
        "Anode catalyst layer" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-30,-10},{-10,10}})));
      replaceable Regions.PEMs.PEM PEM(final L_y=L_y, final L_z=L_z)
        "Proton exchange membrane" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-10,-10},{10,10}})));
      replaceable Regions.CaCLs.CaCL caCL(final L_y=L_y, final L_z=L_z)
        "Cathode catalyst layer" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{10,-10},{30,10}})));
      replaceable Regions.CaGDLs.CaGDL caGDL(final L_y=L_y, final L_z=L_z)
        "Cathode gas diffusion layer" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{30,-10},{50,10}})));
      replaceable Regions.CaFPs.CaFP caFP(final L_y=L_y, final L_z=L_z)
        "Cathode flow plate" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{50,-10},{70,10}})));

    protected
      outer Conditions.Environment environment "Environmental conditions";

    equation
      // Internal connections (between layers)
      connect(anFP.xPositive, anGDL.xNegative) annotation (Line(
          points={{-50,6.10623e-16},{-50,6.10623e-16}},
          color={240,0,0},
          smooth=Smooth.None,
          thickness=0.5));
      connect(anGDL.xPositive, anCL.xNegative) annotation (Line(
          points={{-30,6.10623e-16},{-30,6.10623e-16}},
          color={240,0,0},
          smooth=Smooth.None,
          thickness=0.5));
      connect(anCL.xPositive, PEM.xNegative) annotation (Line(
          points={{-10,6.10623e-16},{-10,6.10623e-16}},
          color={240,0,0},
          smooth=Smooth.None,
          thickness=0.5));
      connect(PEM.xPositive, caCL.xNegative) annotation (Line(
          points={{10,6.10623e-16},{16,-3.36456e-22},{10,6.10623e-16}},
          color={0,0,240},
          smooth=Smooth.None,
          thickness=0.5));
      connect(caCL.xPositive, caGDL.xNegative) annotation (Line(
          points={{30,6.10623e-16},{36,-3.36456e-22},{30,6.10623e-16}},
          color={0,0,240},
          smooth=Smooth.None,
          thickness=0.5));
      connect(caGDL.xPositive, caFP.xNegative) annotation (Line(
          points={{50,6.10623e-16},{56,-3.36456e-22},{50,6.10623e-16}},
          color={0,0,240},
          smooth=Smooth.None,
          thickness=0.5));

      // External connections
      connect(anFPX, anFP.xNegative) annotation (Line(
          points={{-80,5.55112e-16},{-80,6.10623e-16},{-70,6.10623e-16}},
          color={127,127,127},
          smooth=Smooth.None,
          thickness=0.5));
      connect(anFP.yNegative, anFPNegY) annotation (Line(
          points={{-60,-10},{-60,-20}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(anFP.yPositive, anFPPosY) annotation (Line(
          points={{-60,10},{-60,20}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));

      connect(caFP.xPositive, caFPX) annotation (Line(
          points={{70,6.10623e-16},{80,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None,
          thickness=0.5));
      connect(caFP.yNegative, caFPNegY) annotation (Line(
          points={{60,-10},{60,-20}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      connect(caFP.yPositive, caFPPosY) annotation (Line(
          points={{60,10},{60,20}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="cell",
        Documentation(info="
    <html><p>This model presents a single-cell proton exchange membrane fuel cell (PEMFC).  An overview
    of a PEMFC is given in the <a href=\"modelica://FCSys\">top-level documentation of FCSys</a>.</p>

    <p>The output variable <code>Wdot</code> is the total electrical power output through the
    ends of the flowplates over the yz plane.
    The output variable <code>w</code> is the electrical potential difference along the x axis.
    It is averaged over the cell segments on an areal basis.
    The output variable <code>zI</code> is the power-effective electrical current, i.e.,
    the current that would yield the electrical power (<code>Wdot</code>)
    at the electrical potential (<code>w</code>).  The power-effective electrical current (<code>zI</code>)
    is not generally equal to the sum of the actual electrical
    currents of the segments (<code>zI_yz</code>).</p>
    </html>"),
        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-80,-20},{
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
      extends FCSys.BaseClasses.Icons.Cell;

      // Geometric parameters
      parameter Q.Length L_y[:]={U.m}
        "<html>Lengths along the channel (L<sub>y</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      parameter Q.Length L_z[:]={5*U.mm}
        "<html>Lengths across the channel (L<sub>z</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      final parameter Integer n_y=size(L_y, 1)
        "Number of subregions along the channel";
      final parameter Integer n_z=size(L_z, 1)
        "Number of subregions across the channel";

      // Essential analysis variables
      final parameter Q.Area A=anFP.A[Axis.x] "Cross-sectional area";
      output Q.Power Wdot(stateSelect=StateSelect.never) = -sum(anFP.subregions[
        1, :, :].graphite.'e-'.faces[1, Side.n].phi[1] .* anFP.subregions[1, :,
        :].graphite.'e-'.faces[1, Side.n].mPhidot[1] + caFP.subregions[caFP.n_x,
        :, :].graphite.'e-'.faces[1, Side.p].phi[1] .* caFP.subregions[caFP.n_x,
        :, :].graphite.'e-'.faces[1, Side.p].mPhidot[1])
        "Electrical power output";
      output Q.Potential w(stateSelect=StateSelect.never) = sum(anFP.subregions[
        1, :, :].graphite.'e-'.faces[1, Side.n].mPhidot[1] ./ anFP.subregions[1,
        :, :].graphite.'e-'.faces[1, Side.n].rho + caFP.subregions[caFP.n_x, :,
        :].graphite.'e-'.faces[1, Side.p].mPhidot[1] ./ caFP.subregions[caFP.n_x,
        :, :].graphite.'e-'.faces[1, Side.p].rho)/A "Electrical potential";
      output Q.Current zI(stateSelect=StateSelect.never) = Wdot/w
        "Electrical current";

      // Auxiliary variables (for analysis)
      final parameter Q.Length L[Axis]={anFP.L[Axis.x] + anCGDL.L[Axis.x] + PEM.L[
          Axis.x] + caCGDL.L[Axis.x] + caFP.L[Axis.x],anFP.L[Axis.y],anFP.L[
          Axis.z]} if environment.analysis
        "Total lengths along the x, y, and z axes";
      final parameter Q.Volume V=product(L) if environment.analysis "Volume";
      output Q.Power Wdot_yz[n_y, n_z](each stateSelect=StateSelect.never) = -
        anFP.subregions[1, :, :].graphite.'e-'.faces[1, Side.n].phi[1] .* anFP.subregions[
        1, :, :].graphite.'e-'.faces[1, Side.n].mPhidot[1] - caFP.subregions[
        caFP.n_x, :, :].graphite.'e-'.faces[1, Side.p].phi[1] .* caFP.subregions[
        caFP.n_x, :, :].graphite.'e-'.faces[1, Side.p].mPhidot[1] if
        environment.analysis "Electrical power of the segments (x axis)";
      output Q.CurrentAreic zJ_yz[n_y, n_z](each stateSelect=StateSelect.never)
         = -anFP.subregions[1, :, :].graphite.'e-'.faces[1, Side.n].phi[1] .*
        anFP.subregions[1, :, :].graphite.'e-'.faces[1, Side.n].rho if
        environment.analysis
        "Areic electrical current of the segments (x axis)";
      output Q.Current zI_yz[n_y, n_z](each stateSelect=StateSelect.never) =
        zJ_yz .* anFP.subregions[1, :, :].A[Axis.x] if environment.analysis
        "Electrical current of the segments (x axis)";
      output Q.Current Ndot_H2(stateSelect=StateSelect.never) = sum(anFP.subregions[
        :, 1, :].gas.H2.faces[2, Side.n].Ndot + anFP.subregions[:, n_y, :].gas.H2.faces[
        2, Side.p].Ndot) if environment.analysis
        "Rate of hydrogen intake (y axis)";
      output Q.Current Ndot_H2O(stateSelect=StateSelect.never) = sum(anFP.subregions[
        :, 1, :].gas.H2O.faces[2, Side.n].Ndot + anFP.subregions[:, n_y, :].gas.H2O.faces[
        2, Side.p].Ndot) + sum(caFP.subregions[:, 1, :].gas.H2O.faces[1, Side.n].Ndot
         + caFP.subregions[:, n_y, :].gas.H2O.faces[1, Side.p].Ndot) if
        environment.analysis "Rate of water intake (y axis)";
      output Q.Current Ndot_O2(stateSelect=StateSelect.never) = sum(caFP.subregions[
        :, 1, :].gas.O2.faces[1, Side.n].Ndot + caFP.subregions[:, n_y, :].gas.O2.faces[
        1, Side.p].Ndot) if environment.analysis
        "Rate of oxygen intake (y axis)";

      Connectors.FaceBus anFPX[n_y, n_z] "Interface with the anode end plate"
        annotation (Placement(transformation(extent={{-70,-10},{-50,10}},
              rotation=0), iconTransformation(extent={{-110,-10},{-90,10}})));
      Connectors.FaceBus caFPX[n_y, n_z] "Interface with the cathode end plate"
        annotation (Placement(transformation(extent={{50,-10},{70,10}},
              rotation=0), iconTransformation(extent={{90,-10},{110,10}})));
      Connectors.FaceBus anFPNegY[anFP.n_x, n_z] "Negative anode fluid port"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-40,-20}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-40,-100})));
      Connectors.FaceBus caFPNegY[caFP.n_x, n_z] "Negative cathode fluid port"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={40,-20}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={40,-100})));
      Connectors.FaceBus anFPPosY[anFP.n_x, n_z] "Positive anode fluid port"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-40,20}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-40,100})));
      Connectors.FaceBus caFPPosY[caFP.n_x, n_z] "Positive cathode fluid port"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={40,20}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={40,100})));

      replaceable Regions.AnFPs.AnFP anFP(final L_y=L_y, final L_z=L_z)
        "Anode flow plate" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-50,-10},{-30,10}})));
      replaceable Regions.AnCLs.AnCGDL anCGDL(final L_y=L_y, final L_z=L_z)
        "Anode catalyst and gas diffusion layer" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-30,-10},{-10,10}})));
      replaceable Regions.PEMs.PEM PEM(final L_y=L_y, final L_z=L_z)
        "Proton exchange membrane" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-10,-10},{10,10}})));
      replaceable Regions.CaCLs.CaCGDL caCGDL(final L_y=L_y, final L_z=L_z)
        "Cathode catalyst and gas diffusion layer" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{10,-10},{30,10}})));
      replaceable Regions.CaFPs.CaFP caFP(final L_y=L_y, final L_z=L_z)
        "Cathode flow plate" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{30,-10},{50,10}})));

    protected
      outer Conditions.Environment environment "Environmental conditions";

    equation
      // Internal connections (between layers)
      connect(anFP.xPositive, anCGDL.xNegative) annotation (Line(
          points={{-30,6.10623e-16},{-30,6.10623e-16}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(anCGDL.xPositive, PEM.xNegative) annotation (Line(
          points={{-10,6.10623e-16},{-10,6.10623e-16}},
          color={240,0,0},
          smooth=Smooth.None,
          thickness=0.5));
      connect(PEM.xPositive, caCGDL.xNegative) annotation (Line(
          points={{10,6.10623e-16},{16,-3.36456e-22},{10,6.10623e-16}},
          color={0,0,240},
          smooth=Smooth.None,
          thickness=0.5));
      connect(caCGDL.xPositive, caFP.xNegative) annotation (Line(
          points={{30,6.10623e-16},{26,-3.36456e-22},{36,0},{30,6.10623e-16}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));

      // External connections
      connect(anFPX, anFP.xNegative) annotation (Line(
          points={{-60,5.55112e-16},{-60,6.10623e-16},{-50,6.10623e-16}},
          color={127,127,127},
          smooth=Smooth.None,
          thickness=0.5));
      connect(anFP.yNegative, anFPNegY) annotation (Line(
          points={{-40,-10},{-40,-20}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(anFP.yPositive, anFPPosY) annotation (Line(
          points={{-40,10},{-40,20}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));

      connect(caFP.xPositive, caFPX) annotation (Line(
          points={{50,6.10623e-16},{60,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None,
          thickness=0.5));
      connect(caFP.yNegative, caFPNegY) annotation (Line(
          points={{40,-10},{40,-20}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      connect(caFP.yPositive, caFPPosY) annotation (Line(
          points={{40,10},{40,20}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="cell",
        Documentation(info="

<html><p>Please see the documentation of the
  <a href=\"modelica://FCSys.Assemblies.Cells.Cell\">Cell</a> model.</p></html>"),

        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-60,-20},{
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
Copyright 2007&ndash;2012, Georgia Tech Research Corporation.</p>

<p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.ModelicaLicense2\">
FCSys.UsersGuide.ModelicaLicense2</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p>
</html>"));

end Assemblies;
