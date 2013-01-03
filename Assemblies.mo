within FCSys;
package Assemblies "Combinations of regions (e.g., cells)"
  extends Modelica.Icons.Package;
  extends FCSys.BaseClasses.Icons.PackageUnderConstruction;
  package Cells "Single-cell PEMFC models"
    extends Modelica.Icons.Package;
    package Examples "Examples and tests"

      model Cell "<html>Isolated <code>Cell</code> model</html>"
        extends Modelica.Icons.Example;
        inner FCSys.BCs.Defaults defaults(
          analysis=false,
          p=149.6*U.kPa,
          T=333.15*U.K)
          annotation (Placement(transformation(extent={{24,-8},{44,12}})));
        replaceable FCSys.Assemblies.Cells.CellSSIC cell annotation (
            __Dymola_choicesFromPackage=true, Placement(transformation(extent={
                  {-10,-10},{10,10}})));
        annotation (
          experiment(StopTime=1e-24, Tolerance=1e-06),
          experimentSetupOutput,
          Commands(file=
                "resources/scripts/Dymola/Assemblies.Cells.Examples.Cell.mos"),

          Icon(graphics));

      end Cell;
      extends Modelica.Icons.ExamplesPackage;
      model CellProfile
        "Apply boundary conditions to a cell according to a test profile"

        extends FCSys.BCs.TestStands.TestProfile(anEnd(each graphite('incle-'=
                  true, 'e-'(redeclare Modelica.Blocks.Sources.Ramp
                  materialSpec(height=10000*U.A, duration=500)))), caEnd(each
              graphite('incle-'=true, 'e-'(redeclare
                  Modelica.Blocks.Sources.Ramp materialSpec(height=-10000*U.A,
                    duration=500)))));
        replaceable FCSys.Assemblies.Cells.Cell cell annotation (
            __Dymola_choicesFromPackage=true, Placement(transformation(extent={
                  {-10,-10},{10,10}})));

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
            color={240,0,0},
            thickness=0.5,
            smooth=Smooth.None));

        connect(cell.anFPPosY, anSink.face) annotation (Line(
            points={{-4,10},{-4,20},{-20,20},{-20,26}},
            color={240,0,0},
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
          experimentSetupOutput,
          Commands(file=
                "resources/scripts/Dymola/Assemblies.Cells.Examples.CellProfile.mos"),

          experiment(StopTime=600, Tolerance=1e-08),
          experimentSetupOutput,
          Icon(graphics));
      end CellProfile;

      model Polarization "Run a cell polarization"
        extends CellProfile;
        annotation (Commands(file=
                "resources/scripts/Dymola/Assemblies.Cells.Examples.Polarization.mos"));
      end Polarization;

      model CellProfileIO
        "Apply BCs to a cell according to a test profile, with inputs and outputs"
        extends FCSys.Assemblies.Cells.Examples.CellProfile(testStand(final
              inclIO=true));
        extends FCSys.BaseClasses.Icons.Blocks.Continuous;

        FCSys.Connectors.RealInputBus u "Input bus" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=315,
              origin={-30,30}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-2,100})));
        FCSys.Connectors.RealOutputBus y "Output bus" annotation (Placement(
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
        annotation (
          experiment(StopTime=600, Tolerance=1e-08),
          experimentSetupOutput,
          Commands(file=
                "resources/scripts/Dymola/Assemblies.Examples.CellPolarizationstoich.mos"),

          Icon(graphics));
      end CellProfileIO;

      model CellModelica
        "<html>Cell interfaced to components from the <a href=\"modelica://Modelica\">Modelica</a> package</html>"
        extends Modelica.Icons.Example;
        extends Modelica.Icons.UnderConstruction;
        FCSys.Assemblies.Cells.Cell cell(anFP(redeclare
              FCSys.Subregions.Subregion subregions(
              each final inclX=true,
              each inclY=true,
              each graphite('incle-'=true, 'e-'(perfectMaterialDiff={{{{true,
                      false}}}})),
              each gas(inclH2=true, inclH2O=true))))
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

        inner FCSys.BCs.Defaults defaults(analysis=false)
          annotation (Placement(transformation(extent={{40,60},{60,80}})));
        FCSys.BCs.Adapters.'AdaptBuse-' caModelicaAdapt(A=cell.L_y[1]*cell.L_z[
              1])
          annotation (Placement(transformation(extent={{20,-10},{40,10}})));
        FCSys.BCs.Adapters.'AdaptBuse-' anModelicaAdapt(A=cell.L_y[1]*cell.L_z[
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
        annotation (experiment(StopTime=1000), experimentSetupOutput);
      end CellModelica;
    end Examples;

    model Cell "Default single-cell PEMFC"

      extends FCSys.BaseClasses.Icons.Cell;

      // Geometric parameters
      parameter Q.Length L_y[:]=fill(1*U.m/1, 1)
        "<html>Lengths along the channel (L<sub>y</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      parameter Q.Length L_z[:]=fill(5*U.mm/1, 1)
        "<html>Lengths across the channel (L<sub>z</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      final parameter Integer n_y=size(L_y, 1)
        "Number of subregions along the channel";
      final parameter Integer n_z=size(L_z, 1)
        "Number of subregions across the channel";
      // TODO:  For GM cell, use dissimilar L_y and L_z for anode, cathode, and PEM.

      // Essential analysis variables (x-axis electrical voltage, power, and current)
      output Q.Potential v=average(average(anFP.subregions[1, :, :].graphite.
          'e-'.mu_face[1, 1] - caFP.subregions[caFP.n_x, :, :].graphite.'e-'.mu_face[
          1, 2])) "Average electrical potential (x axis)";
      output Q.Power Wdot=-sum(anFP.subregions[1, :, :].graphite.'e-'.mu_face[1,
          1] .* anFP.subregions[1, :, :].graphite.'e-'.Ndot_face[1, 1] + caFP.subregions[
          caFP.n_x, :, :].graphite.'e-'.mu_face[1, 2] .* caFP.subregions[caFP.n_x,
          :, :].graphite.'e-'.Ndot_face[1, 2])
        "Electrical power output (x axis)";
      output Q.Current I=Wdot/v "Current";
      // This is the power-effective current that would create the actual
      // electrical power at the average voltage (above).

      // Auxiliary variables (for analysis)
      final parameter Q.Length L[Axis](each stateSelect=StateSelect.never) = {
        sum(anFP.L_x) + sum(anGDL.L_x) + sum(anCL.L_x) + sum(pEM.L_x) + sum(
        caCL.L_x) + sum(caGDL.L_x) + sum(caFP.L_x),sum(L_y),sum(L_z)} if
        defaults.analysis "Total lengths along the x, y, and z axes";
      final parameter Q.Area A[Axis](each stateSelect=StateSelect.never) = {L[
        cartWrap(ax + 2)]*L[cartWrap(ax + 2)] for ax in 1:3} if defaults.analysis
        "Cross-sectional areas";
      final parameter Q.Volume V=product(L) if defaults.analysis "Volume";
      output Q.Potential Deltav_x_seg[n_y, n_z](each stateSelect=StateSelect.never)
         = caFP.subregions[caFP.n_x, :, :].graphite.'e-'.mu_face[1, 2] - anFP.subregions[
        1, :, :].graphite.'e-'.mu_face[1, 1] if defaults.analysis
        "Electrical potential differences of the segments (x axis)";
      output Q.Power 'Wdot_e-_x'[n_y, n_z](each stateSelect=StateSelect.never)
         = -(anFP.subregions[1, :, :].graphite.'e-'.mu_face[1, 1] .* anFP.subregions[
        1, :, :].graphite.'e-'.Ndot_face[1, 1] + caFP.subregions[caFP.n_x, :, :].graphite.
        'e-'.mu_face[1, 2] .* caFP.subregions[caFP.n_x, :, :].graphite.'e-'.Ndot_face[
        1, 2]) if defaults.analysis "Electrical power of the segments (x axis)";
      output Q.Current I_x_seg[n_y, n_z](each stateSelect=StateSelect.never) =
        (caFP.subregions[caFP.n_x, :, :].graphite.'e-'.Ndot_face[1, 2] - anFP.subregions[
        1, :, :].graphite.'e-'.Ndot_face[1, 1])/2 if defaults.analysis
        "Electrical currents of the segments (x axis)";
      output Q.CurrentAreic Iprimeprime_x[n_y, n_z](each stateSelect=
            StateSelect.never) = {I_x_seg[i_y, i_z]/(L_y[i_y]*L_z[i_z]) for i_z
         in 1:n_z, i_y in 1:n_y} if defaults.analysis
        "Areic electrical current of the segments (x axis)";
      output Q.CurrentAreic Iprimeprime_x_avg(stateSelect=StateSelect.never) =
        I/A[1] if defaults.analysis "Average areic electrical current (x axis)";
      output Q.Current Ndot_H2(stateSelect=StateSelect.never) = sum(anFP.subregions[
        1, :, :].gas.H2.Ndot_face[1, 1]) + sum(anFP.subregions[:, 1, :].gas.H2.Ndot_face[
        2, 1] + anFP.subregions[:, n_y, :].gas.H2.Ndot_face[2, 2]) if defaults.analysis
        "Rate of hydrogen intake";
      output Q.Current Ndot_H2O(stateSelect=StateSelect.never) = sum(anFP.subregions[
        1, :, :].gas.H2O.Ndot_face[1, 1] + caFP.subregions[caFP.n_x, :, :].gas.H2O.Ndot_face[
        1, 2]) + sum(anFP.subregions[:, 1, :].gas.H2O.Ndot_face[1, 1] + anFP.subregions[
        :, n_y, :].gas.H2O.Ndot_face[1, 2]) + sum(caFP.subregions[:, 1, :].gas.H2O.Ndot_face[
        1, 1] - caFP.subregions[:, n_y, :].gas.H2O.Ndot_face[1, 2]) if defaults.analysis
        "Rate of water intake";
      output Q.Current Ndot_O2(stateSelect=StateSelect.never) = sum(caFP.subregions[
        caFP.n_x, :, :].gas.O2.Ndot_face[1, 2]) + sum(caFP.subregions[:, 1, :].gas.O2.Ndot_face[
        1, 1] - caFP.subregions[:, n_y, :].gas.O2.Ndot_face[1, 2]) if defaults.analysis
        "Rate of oxygen intake";

      FCSys.Connectors.FaceBus anFPX[n_y, n_z] "Anode plate face" annotation (
          Placement(transformation(extent={{-90,-10},{-70,10}},rotation=0),
            iconTransformation(extent={{-110,-10},{-90,10}})));
      FCSys.Connectors.FaceBus caFPX[n_y, n_z] "Cathode plate face" annotation
        (Placement(transformation(extent={{70,-10},{90,10}}, rotation=0),
            iconTransformation(extent={{90,-10},{110,10}})));
      FCSys.Connectors.FaceBus anFPPosY[anFP.n_x, n_z]
        "Positive anode flow plate face along the y axis" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-60,20}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-40,100})));
      FCSys.Connectors.FaceBus caFPPosY[caFP.n_x, n_z]
        "Positive anode flow plate face along the y axis" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={60,20}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={40,100})));
      FCSys.Connectors.FaceBus caFPNegY[caFP.n_x, n_z]
        "Negative cathode flow plate face along the y axis" annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={60,-20}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={40,-100})));
      FCSys.Connectors.FaceBus anFPNegY[anFP.n_x, n_z]
        "Negative anode flow plate face along the y axis" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-60,-20}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-40,-100})));

      replaceable FCSys.Regions.AnFPs.AnFP anFP(final L_y=L_y, final L_z=L_z)
        "Anode flow plate" annotation (Dialog(group="Layers"), Placement(
            transformation(extent={{-70,-10},{-50,10}})));

      replaceable FCSys.Regions.AnGDLs.AnGDL anGDL(final L_y=L_y, final L_z=L_z)
        "Anode gas diffusion layer" annotation (Dialog(group="Layers"),
          Placement(transformation(extent={{-50,-10},{-30,10}})));

      replaceable FCSys.Regions.AnCLs.AnCL anCL(final L_y=L_y, final L_z=L_z)
        "Anode catalyst layer" annotation (Dialog(group="Layers"), Placement(
            transformation(extent={{-30,-10},{-10,10}})));

      replaceable FCSys.Regions.PEMs.PEM pEM(final L_y=L_y, final L_z=L_z)
        "Proton exchange membrane" annotation (Dialog(group="Layers"),
          Placement(transformation(extent={{-10,-10},{10,10}})));

      replaceable FCSys.Regions.CaCLs.CaCL caCL(final L_y=L_y, final L_z=L_z)
        "Cathode catalyst layer" annotation (Dialog(group="Layers"), Placement(
            transformation(extent={{10,-10},{30,10}})));

      replaceable FCSys.Regions.CaGDLs.CaGDL caGDL(final L_y=L_y, final L_z=L_z)
        "Cathode gas diffusion layer" annotation (Dialog(group="Layers"),
          Placement(transformation(extent={{30,-10},{50,10}})));

      replaceable FCSys.Regions.CaFPs.CaFP caFP(final L_y=L_y, final L_z=L_z)
        "Cathode flow plate" annotation (Dialog(group="Layers"), Placement(
            transformation(extent={{50,-10},{70,10}})));

    protected
      outer FCSys.BCs.Defaults defaults "Environmental properties and settings";

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
      connect(anCL.xPositive, pEM.xNegative) annotation (Line(
          points={{-10,6.10623e-16},{-10,6.10623e-16}},
          color={240,0,0},
          smooth=Smooth.None,
          thickness=0.5));
      connect(pEM.xPositive, caCL.xNegative) annotation (Line(
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
of a PEMFC is given in the top-level documentation of <a href=\"modelica://FCSys\">FCSys</a>.</html>"),

        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-80,-20},{
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
        experimentSetupOutput,
        experiment(StopTime=120, Tolerance=1e-06));
    end Cell;

    model CalibratedCell
      "Cell model with calibration parameters for exchange and transport"
      import FCSys.Subregions.Species;

      // Exchange of linear momentum
      parameter Q.NumberAbsolute k_alpha_tau_C(
        final min=0,
        final nominal=1) = 1
        "<html>For C (<i>k</i><sub>&alpha; &Phi; C</sub>)</html>" annotation (
          Dialog(tab="Calibration factors", group="Exchange of linear momentum"));
      parameter Q.NumberAbsolute k_alpha_tau_C19HF37O5S(
        final min=0,
        final nominal=1) = 1
        "<html>For C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S (<i>k</i><sub>&alpha; &Phi; C19HF37O5S</sub>)</html>"
        annotation (Dialog(tab="Calibration factors", group=
              "Exchange of linear momentum"));
      parameter Q.NumberAbsolute 'k_alpha_tau_e-'(
        final min=0,
        final nominal=1) = 1
        "<html>For e<sup>-</sup> (<i>k</i><sub>&alpha; &Phi; e-</sub>)</html>"
        annotation (Dialog(tab="Calibration factors", group=
              "Exchange of linear momentum"));
      parameter Q.NumberAbsolute k_alpha_tau_H2(
        final min=0,
        final nominal=1) = 1
        "<html>For H<sub>2</sub> (<i>k</i><sub>&alpha; &Phi; H2</sub>)</html>"
        annotation (Dialog(tab="Calibration factors", group=
              "Exchange of linear momentum"));
      parameter Q.NumberAbsolute k_alpha_tau_H2O(
        final min=0,
        final nominal=1) = 1
        "<html>For H<sub>2</sub>O (<i>k</i><sub>&alpha; &Phi; H2O</sub>)</html>"
        annotation (Dialog(tab="Calibration factors", group=
              "Exchange of linear momentum"));
      parameter Q.NumberAbsolute 'k_alpha_tau_H+'(
        final min=0,
        final nominal=1) = 1
        "<html>For H<sup>+</sup> (<i>k</i><sub>&alpha; &Phi; H+</sub>)</html>"
        annotation (Dialog(tab="Calibration factors", group=
              "Exchange of linear momentum"));
      parameter Q.NumberAbsolute k_alpha_tau_N2(
        final min=0,
        final nominal=1) = 1
        "<html>For N<sub>2</sub> (<i>k</i><sub>&alpha; &Phi; N2</sub>)</html>"
        annotation (Dialog(tab="Calibration factors", group=
              "Exchange of linear momentum"));
      parameter Q.NumberAbsolute k_alpha_tau_O2(
        final min=0,
        final nominal=1) = 1
        "<html>For O<sub>2</sub> (<i>k</i><sub>&alpha; &Phi; O2</sub>)</html>"
        annotation (Dialog(tab="Calibration factors", group=
              "Exchange of linear momentum"));

      // Thermal exchange
      parameter Q.NumberAbsolute k_alpha_Qdot_C(
        final min=0,
        final nominal=1) = 1
        "<html>For C (<i>k</i><sub><i>S</i> C</sub>)</html>" annotation (Dialog(
            tab="Calibration factors", group="Thermal exchange"));
      parameter Q.NumberAbsolute k_alpha_Qdot_C19HF37O5S(
        final min=0,
        final nominal=1) = 1
        "<html>For C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S (<i>k</i><sub><i>S</i> C19HF37O5S</sub>)</html>"
        annotation (Dialog(tab="Calibration factors", group="Thermal exchange"));
      parameter Q.NumberAbsolute 'k_alpha_Qdot_e-'(
        final min=0,
        final nominal=1) = 1
        "<html>For e<sup>-</sup> (<i>k</i><sub><i>S</i> e-</sub>)</html>"
        annotation (Dialog(tab="Calibration factors", group="Thermal exchange"));
      parameter Q.NumberAbsolute k_alpha_Qdot_H2(
        final min=0,
        final nominal=1) = 1
        "<html>For H<sub>2</sub> (<i>k</i><sub><i>S</i> H2</sub>)</html>"
        annotation (Dialog(tab="Calibration factors", group="Thermal exchange"));
      parameter Q.NumberAbsolute k_alpha_Qdot_H2O(
        final min=0,
        final nominal=1) = 1
        "<html>For H<sub>2</sub>O (<i>k</i><sub><i>S</i> H2O</sub>)</html>"
        annotation (Dialog(tab="Calibration factors", group="Thermal exchange"));
      parameter Q.NumberAbsolute 'k_alpha_Qdot_H+'(
        final min=0,
        final nominal=1) = 1
        "<html>For H<sup>+</sup> (<i>k</i><sub><i>S</i> H+</sub>)</html>"
        annotation (Dialog(tab="Calibration factors", group="Thermal exchange"));
      parameter Q.NumberAbsolute k_alpha_Qdot_N2(
        final min=0,
        final nominal=1) = 1
        "<html>For N<sub>2</sub> (<i>k</i><sub><i>S</i> N2</sub>)</html>"
        annotation (Dialog(tab="Calibration factors", group="Thermal exchange"));
      parameter Q.NumberAbsolute k_alpha_Qdot_O2(
        final min=0,
        final nominal=1) = 1
        "<html>For O<sub>2</sub> (<i>k</i><sub><i>S</i> O2</sub>)</html>"
        annotation (Dialog(tab="Calibration factors", group="Thermal exchange"));

      // Material transport
      /*
  parameter Q.NumberAbsolute k_alpha_Ndot_C(
    final min=0,
    final nominal=1) = 1 "<html>For C (<i>k</i><sub>&alpha; <i>N&#775;</i> &nbsp;C</sub>)</html>"
    annotation (Dialog(tab="Calibration factors", group="Material transport"));
  parameter Q.NumberAbsolute k_alpha_Ndot_C19HF37O5S(
    final min=0,
    final nominal=1) = 1
    "<html>For C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S (<i>k</i><sub>&alpha; <i>N&#775;</i> &nbsp;C19HF37O5S</sub>)</html>"
    annotation (Dialog(tab="Calibration factors", group="Material transport"));
  */
      parameter Q.NumberAbsolute 'k_alpha_Ndot_e-'(
        final min=0,
        final nominal=1) = 1
        "<html>For e<sup>-</sup> (<i>k</i><sub>&alpha; <i>N&#775;</i> &nbsp;e-</sub>)</html>"
        annotation (Dialog(tab="Calibration factors", group=
              "Material transport"));
      parameter Q.NumberAbsolute k_alpha_Ndot_H2(
        final min=0,
        final nominal=1) = 1
        "<html>For H<sub>2</sub> (<i>k</i><sub>&alpha; <i>N&#775;</i> &nbsp;H2</sub>)</html>"
        annotation (Dialog(tab="Calibration factors", group=
              "Material transport"));
      parameter Q.NumberAbsolute k_alpha_Ndot_H2O(
        final min=0,
        final nominal=1) = 1
        "<html>For H<sub>2</sub>O (<i>k</i><sub>&alpha; <i>N&#775;</i> &nbsp;H2O</sub>)</html>"
        annotation (Dialog(tab="Calibration factors", group=
              "Material transport"));
      parameter Q.NumberAbsolute 'k_alpha_Ndot_H+'(
        final min=0,
        final nominal=1) = 1
        "<html>For H<sup>+</sup> (<i>k</i><sub>&alpha; <i>N&#775;</i> &nbsp;H+</sub>)</html>"
        annotation (Dialog(tab="Calibration factors", group=
              "Material transport"));
      parameter Q.NumberAbsolute k_alpha_Ndot_N2(
        final min=0,
        final nominal=1) = 1
        "<html>For N<sub>2</sub> (<i>k</i><sub>&alpha; <i>N&#775;</i> &nbsp;N2</sub>)</html>"
        annotation (Dialog(tab="Calibration factors", group=
              "Material transport"));
      parameter Q.NumberAbsolute k_alpha_Ndot_O2(
        final min=0,
        final nominal=1) = 1
        "<html>For O<sub>2</sub> (<i>k</i><sub>&alpha; <i>N&#775;</i> &nbsp;O2</sub>)</html>"
        annotation (Dialog(tab="Calibration factors", group=
              "Material transport"));

      // Transport of linear momentum
      /*
  parameter Q.NumberAbsolute k_alpha_tau_C(
    final min=0,
    final nominal=1) = 1 "<html>For C (<i>k</i><sub>&Phi; C</sub>)</html>"
    annotation (Dialog(tab="Calibration factors", group="Transport of linear momentum"))
    ;
  parameter Q.NumberAbsolute k_alpha_tau_C19HF37O5S(
    final min=0,
    final nominal=1) = 1
    "<html>For C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S (<i>k</i><sub>&Phi; C19HF37O5S</sub>)</html>"
    annotation (Dialog(tab="Calibration factors", group="Transport of linear momentum"))
    ;
  */

      // Thermal transport

      extends Cell(
        anFP(subregions(each gas(
              redeclare Species.H2.Gas.Calibrated H2(
                final k_alpha_tau=k_alpha_tau_H2,
                final k_alpha_Qdot=k_alpha_Qdot_H2,
                final k_alpha_Ndot=k_alpha_Ndot_H2,
                final k_alpha_tau=k_alpha_tau_H2,
                final k_alpha_Qdot=k_alpha_Qdot_H2),
              redeclare Species.H2O.Gas.Calibrated H2O(
                final k_alpha_tau=k_alpha_tau_H2O,
                final k_alpha_Qdot=k_alpha_Qdot_H2O,
                final k_alpha_Ndot=k_alpha_Ndot_H2O,
                final k_alpha_tau=k_alpha_tau_H2O,
                final k_alpha_Qdot=k_alpha_Qdot_H2O),
              redeclare Species.N2.Gas.Calibrated N2(
                final k_alpha_tau=k_alpha_tau_N2,
                final k_alpha_Qdot=k_alpha_Qdot_N2,
                final k_alpha_Ndot=k_alpha_Ndot_N2,
                final k_alpha_tau=k_alpha_tau_N2,
                final k_alpha_Qdot=k_alpha_Qdot_N2),
              redeclare Species.O2.Gas.Calibrated O2(
                final k_alpha_tau=k_alpha_tau_O2,
                final k_alpha_Qdot=k_alpha_Qdot_O2,
                final k_alpha_Ndot=k_alpha_Ndot_O2,
                final k_alpha_tau=k_alpha_tau_O2,
                final k_alpha_Qdot=k_alpha_Qdot_O2)), each graphite(redeclare
                Species.C.Graphite.Calibrated C(
                final k_alpha_tau=k_alpha_tau_C,
                final k_alpha_Qdot=k_alpha_Qdot_C,
                final k_alpha_Qdot=k_alpha_Qdot_C), redeclare
                Species.'e-'.solid.Calibrated 'e-'(
                final k_alpha_tau='k_alpha_tau_e-',
                final k_alpha_Qdot='k_alpha_Qdot_e-',
                final k_alpha_Ndot='k_alpha_Ndot_e-',
                final k_alpha_tau='k_alpha_tau_e-',
                final k_alpha_Qdot='k_alpha_Qdot_e-')))),
        anGDL(subregions(each gas(
              redeclare Species.H2.Gas.Calibrated H2(
                final k_alpha_tau=k_alpha_tau_H2,
                final k_alpha_Qdot=k_alpha_Qdot_H2,
                final k_alpha_Ndot=k_alpha_Ndot_H2,
                final k_alpha_tau=k_alpha_tau_H2,
                final k_alpha_Qdot=k_alpha_Qdot_H2),
              redeclare Species.H2O.Gas.Calibrated H2O(
                final k_alpha_tau=k_alpha_tau_H2O,
                final k_alpha_Qdot=k_alpha_Qdot_H2O,
                final k_alpha_Ndot=k_alpha_Ndot_H2O,
                final k_alpha_tau=k_alpha_tau_H2O,
                final k_alpha_Qdot=k_alpha_Qdot_H2O),
              redeclare Species.N2.Gas.Calibrated N2(
                final k_alpha_tau=k_alpha_tau_N2,
                final k_alpha_Qdot=k_alpha_Qdot_N2,
                final k_alpha_Ndot=k_alpha_Ndot_N2,
                final k_alpha_tau=k_alpha_tau_N2,
                final k_alpha_Qdot=k_alpha_Qdot_N2),
              redeclare Species.O2.Gas.Calibrated O2(
                final k_alpha_tau=k_alpha_tau_O2,
                final k_alpha_Qdot=k_alpha_Qdot_O2,
                final k_alpha_Ndot=k_alpha_Ndot_O2,
                final k_alpha_tau=k_alpha_tau_O2,
                final k_alpha_Qdot=k_alpha_Qdot_O2)), each graphite(redeclare
                Species.C.Graphite.Calibrated C(
                final k_alpha_tau=k_alpha_tau_C,
                final k_alpha_Qdot=k_alpha_Qdot_C,
                final k_alpha_Qdot=k_alpha_Qdot_C), redeclare
                Species.'e-'.solid.Calibrated 'e-'(
                final k_alpha_tau='k_alpha_tau_e-',
                final k_alpha_Qdot='k_alpha_Qdot_e-',
                final k_alpha_Ndot='k_alpha_Ndot_e-',
                final k_alpha_tau='k_alpha_tau_e-',
                final k_alpha_Qdot='k_alpha_Qdot_e-')))),
        anCL(subregions(
            each gas(
              redeclare Species.H2.Gas.Calibrated H2(
                final k_alpha_tau=k_alpha_tau_H2,
                final k_alpha_Qdot=k_alpha_Qdot_H2,
                final k_alpha_Ndot=k_alpha_Ndot_H2,
                final k_alpha_tau=k_alpha_tau_H2,
                final k_alpha_Qdot=k_alpha_Qdot_H2),
              redeclare Species.H2O.Gas.Calibrated H2O(
                final k_alpha_tau=k_alpha_tau_H2O,
                final k_alpha_Qdot=k_alpha_Qdot_H2O,
                final k_alpha_Ndot=k_alpha_Ndot_H2O,
                final k_alpha_tau=k_alpha_tau_H2O,
                final k_alpha_Qdot=k_alpha_Qdot_H2O),
              redeclare Species.N2.Gas.Calibrated N2(
                final k_alpha_tau=k_alpha_tau_N2,
                final k_alpha_Qdot=k_alpha_Qdot_N2,
                final k_alpha_Ndot=k_alpha_Ndot_N2,
                final k_alpha_tau=k_alpha_tau_N2,
                final k_alpha_Qdot=k_alpha_Qdot_N2),
              redeclare Species.O2.Gas.Calibrated O2(
                final k_alpha_tau=k_alpha_tau_O2,
                final k_alpha_Qdot=k_alpha_Qdot_O2,
                final k_alpha_Ndot=k_alpha_Ndot_O2,
                final k_alpha_tau=k_alpha_tau_O2,
                final k_alpha_Qdot=k_alpha_Qdot_O2)),
            each graphite(redeclare Species.C.Graphite.Calibrated C(
                final k_alpha_tau=k_alpha_tau_C,
                final k_alpha_Qdot=k_alpha_Qdot_C,
                final k_alpha_Qdot=k_alpha_Qdot_C), redeclare
                Species.'e-'.solid.Calibrated 'e-'(
                final k_alpha_tau='k_alpha_tau_e-',
                final k_alpha_Qdot='k_alpha_Qdot_e-',
                final k_alpha_Ndot='k_alpha_Ndot_e-',
                final k_alpha_tau='k_alpha_tau_e-',
                final k_alpha_Qdot='k_alpha_Qdot_e-')),
            each ionomer(
              redeclare Species.C19HF37O5S.Solid.Calibrated C19HF37O5S(
                final k_alpha_tau=k_alpha_tau_C19HF37O5S,
                final k_alpha_Qdot=k_alpha_Qdot_C19HF37O5S,
                final k_alpha_Qdot=k_alpha_Qdot_C19HF37O5S),
              redeclare Species.H2O.Gas.Calibrated H2O(
                final k_alpha_tau=k_alpha_tau_H2O,
                final k_alpha_Qdot=k_alpha_Qdot_H2O,
                final k_alpha_Ndot=k_alpha_Ndot_H2O,
                final k_alpha_tau=k_alpha_tau_H2O,
                final k_alpha_Qdot=k_alpha_Qdot_H2O),
              redeclare Species.'H+'.Solid.Calibrated 'H+'(
                final k_alpha_tau='k_alpha_tau_H+',
                final k_alpha_Qdot='k_alpha_Qdot_H+',
                final k_alpha_Ndot='k_alpha_Ndot_H+',
                final k_alpha_tau='k_alpha_tau_H+',
                final k_alpha_Qdot='k_alpha_Qdot_H+')))),
        pEM(subregions(each gas(
              redeclare Species.H2.Gas.Calibrated H2(
                final k_alpha_tau=k_alpha_tau_H2,
                final k_alpha_Qdot=k_alpha_Qdot_H2,
                final k_alpha_Ndot=k_alpha_Ndot_H2,
                final k_alpha_tau=k_alpha_tau_H2,
                final k_alpha_Qdot=k_alpha_Qdot_H2),
              redeclare Species.H2O.Gas.Calibrated H2O(
                final k_alpha_tau=k_alpha_tau_H2O,
                final k_alpha_Qdot=k_alpha_Qdot_H2O,
                final k_alpha_Ndot=k_alpha_Ndot_H2O,
                final k_alpha_tau=k_alpha_tau_H2O,
                final k_alpha_Qdot=k_alpha_Qdot_H2O),
              redeclare Species.N2.Gas.Calibrated N2(
                final k_alpha_tau=k_alpha_tau_N2,
                final k_alpha_Qdot=k_alpha_Qdot_N2,
                final k_alpha_Ndot=k_alpha_Ndot_N2,
                final k_alpha_tau=k_alpha_tau_N2,
                final k_alpha_Qdot=k_alpha_Qdot_N2),
              redeclare Species.O2.Gas.Calibrated O2(
                final k_alpha_tau=k_alpha_tau_O2,
                final k_alpha_Qdot=k_alpha_Qdot_O2,
                final k_alpha_Ndot=k_alpha_Ndot_O2,
                final k_alpha_tau=k_alpha_tau_O2,
                final k_alpha_Qdot=k_alpha_Qdot_O2)), each ionomer(
              redeclare Species.C19HF37O5S.Solid.Calibrated C19HF37O5S(
                final k_alpha_tau=k_alpha_tau_C19HF37O5S,
                final k_alpha_Qdot=k_alpha_Qdot_C19HF37O5S,
                final k_alpha_Qdot=k_alpha_Qdot_C19HF37O5S),
              redeclare Species.H2O.Gas.Calibrated H2O(
                final k_alpha_tau=k_alpha_tau_H2O,
                final k_alpha_Qdot=k_alpha_Qdot_H2O,
                final k_alpha_Ndot=k_alpha_Ndot_H2O,
                final k_alpha_tau=k_alpha_tau_H2O,
                final k_alpha_Qdot=k_alpha_Qdot_H2O),
              redeclare Species.'H+'.Solid.Calibrated 'H+'(
                final k_alpha_tau='k_alpha_tau_H+',
                final k_alpha_Qdot='k_alpha_Qdot_H+',
                final k_alpha_Ndot='k_alpha_Ndot_H+',
                final k_alpha_tau='k_alpha_tau_H+',
                final k_alpha_Qdot='k_alpha_Qdot_H+')))),
        caCL(subregions(
            each gas(
              redeclare Species.H2.Gas.Calibrated H2(
                final k_alpha_tau=k_alpha_tau_H2,
                final k_alpha_Qdot=k_alpha_Qdot_H2,
                final k_alpha_Ndot=k_alpha_Ndot_H2,
                final k_alpha_tau=k_alpha_tau_H2,
                final k_alpha_Qdot=k_alpha_Qdot_H2),
              redeclare Species.H2O.Gas.Calibrated H2O(
                final k_alpha_tau=k_alpha_tau_H2O,
                final k_alpha_Qdot=k_alpha_Qdot_H2O,
                final k_alpha_Ndot=k_alpha_Ndot_H2O,
                final k_alpha_tau=k_alpha_tau_H2O,
                final k_alpha_Qdot=k_alpha_Qdot_H2O),
              redeclare Species.N2.Gas.Calibrated N2(
                final k_alpha_tau=k_alpha_tau_N2,
                final k_alpha_Qdot=k_alpha_Qdot_N2,
                final k_alpha_Ndot=k_alpha_Ndot_N2,
                final k_alpha_tau=k_alpha_tau_N2,
                final k_alpha_Qdot=k_alpha_Qdot_N2),
              redeclare Species.O2.Gas.Calibrated O2(
                final k_alpha_tau=k_alpha_tau_O2,
                final k_alpha_Qdot=k_alpha_Qdot_O2,
                final k_alpha_Ndot=k_alpha_Ndot_O2,
                final k_alpha_tau=k_alpha_tau_O2,
                final k_alpha_Qdot=k_alpha_Qdot_O2)),
            each graphite(redeclare Species.C.Graphite.Calibrated C(
                final k_alpha_tau=k_alpha_tau_C,
                final k_alpha_Qdot=k_alpha_Qdot_C,
                final k_alpha_Qdot=k_alpha_Qdot_C), redeclare
                Species.'e-'.solid.Calibrated 'e-'(
                final k_alpha_tau='k_alpha_tau_e-',
                final k_alpha_Qdot='k_alpha_Qdot_e-',
                final k_alpha_Ndot='k_alpha_Ndot_e-',
                final k_alpha_tau='k_alpha_tau_e-',
                final k_alpha_Qdot='k_alpha_Qdot_e-')),
            each ionomer(
              redeclare Species.C19HF37O5S.Solid.Calibrated C19HF37O5S(
                final k_alpha_tau=k_alpha_tau_C19HF37O5S,
                final k_alpha_Qdot=k_alpha_Qdot_C19HF37O5S,
                final k_alpha_Qdot=k_alpha_Qdot_C19HF37O5S),
              redeclare Species.H2O.Gas.Calibrated H2O(
                final k_alpha_tau=k_alpha_tau_H2O,
                final k_alpha_Qdot=k_alpha_Qdot_H2O,
                final k_alpha_Ndot=k_alpha_Ndot_H2O,
                final k_alpha_tau=k_alpha_tau_H2O,
                final k_alpha_Qdot=k_alpha_Qdot_H2O),
              redeclare Species.'H+'.Solid.Calibrated 'H+'(
                final k_alpha_tau='k_alpha_tau_H+',
                final k_alpha_Qdot='k_alpha_Qdot_H+',
                final k_alpha_Ndot='k_alpha_Ndot_H+',
                final k_alpha_tau='k_alpha_tau_H+',
                final k_alpha_Qdot='k_alpha_Qdot_H+')))),
        caGDL(subregions(each gas(
              redeclare Species.H2.Gas.Calibrated H2(
                final k_alpha_tau=k_alpha_tau_H2,
                final k_alpha_Qdot=k_alpha_Qdot_H2,
                final k_alpha_Ndot=k_alpha_Ndot_H2,
                final k_alpha_tau=k_alpha_tau_H2,
                final k_alpha_Qdot=k_alpha_Qdot_H2),
              redeclare Species.H2O.Gas.Calibrated H2O(
                final k_alpha_tau=k_alpha_tau_H2O,
                final k_alpha_Qdot=k_alpha_Qdot_H2O,
                final k_alpha_Ndot=k_alpha_Ndot_H2O,
                final k_alpha_tau=k_alpha_tau_H2O,
                final k_alpha_Qdot=k_alpha_Qdot_H2O),
              redeclare Species.N2.Gas.Calibrated N2(
                final k_alpha_tau=k_alpha_tau_N2,
                final k_alpha_Qdot=k_alpha_Qdot_N2,
                final k_alpha_Ndot=k_alpha_Ndot_N2,
                final k_alpha_tau=k_alpha_tau_N2,
                final k_alpha_Qdot=k_alpha_Qdot_N2),
              redeclare Species.O2.Gas.Calibrated O2(
                final k_alpha_tau=k_alpha_tau_O2,
                final k_alpha_Qdot=k_alpha_Qdot_O2,
                final k_alpha_Ndot=k_alpha_Ndot_O2,
                final k_alpha_tau=k_alpha_tau_O2,
                final k_alpha_Qdot=k_alpha_Qdot_O2)), each graphite(redeclare
                Species.C.Graphite.Calibrated C(
                final k_alpha_tau=k_alpha_tau_C,
                final k_alpha_Qdot=k_alpha_Qdot_C,
                final k_alpha_Qdot=k_alpha_Qdot_C), redeclare
                Species.'e-'.solid.Calibrated 'e-'(
                final k_alpha_tau='k_alpha_tau_e-',
                final k_alpha_Qdot='k_alpha_Qdot_e-',
                final k_alpha_Ndot='k_alpha_Ndot_e-',
                final k_alpha_tau='k_alpha_tau_e-',
                final k_alpha_Qdot='k_alpha_Qdot_e-')))),
        caFP(subregions(each gas(
              redeclare Species.H2.Gas.Calibrated H2(
                final k_alpha_tau=k_alpha_tau_H2,
                final k_alpha_Qdot=k_alpha_Qdot_H2,
                final k_alpha_Ndot=k_alpha_Ndot_H2,
                final k_alpha_tau=k_alpha_tau_H2,
                final k_alpha_Qdot=k_alpha_Qdot_H2),
              redeclare Species.H2O.Gas.Calibrated H2O(
                final k_alpha_tau=k_alpha_tau_H2O,
                final k_alpha_Qdot=k_alpha_Qdot_H2O,
                final k_alpha_Ndot=k_alpha_Ndot_H2O,
                final k_alpha_tau=k_alpha_tau_H2O,
                final k_alpha_Qdot=k_alpha_Qdot_H2O),
              redeclare Species.N2.Gas.Calibrated N2(
                final k_alpha_tau=k_alpha_tau_N2,
                final k_alpha_Qdot=k_alpha_Qdot_N2,
                final k_alpha_Ndot=k_alpha_Ndot_N2,
                final k_alpha_tau=k_alpha_tau_N2,
                final k_alpha_Qdot=k_alpha_Qdot_N2),
              redeclare Species.O2.Gas.Calibrated O2(
                final k_alpha_tau=k_alpha_tau_O2,
                final k_alpha_Qdot=k_alpha_Qdot_O2,
                final k_alpha_Ndot=k_alpha_Ndot_O2,
                final k_alpha_tau=k_alpha_tau_O2,
                final k_alpha_Qdot=k_alpha_Qdot_O2)), each graphite(redeclare
                Species.C.Graphite.Calibrated C(
                final k_alpha_tau=k_alpha_tau_C,
                final k_alpha_Qdot=k_alpha_Qdot_C,
                final k_alpha_Qdot=k_alpha_Qdot_C), redeclare
                Species.'e-'.solid.Calibrated 'e-'(
                final k_alpha_tau='k_alpha_tau_e-',
                final k_alpha_Qdot='k_alpha_Qdot_e-',
                final k_alpha_Ndot='k_alpha_Ndot_e-',
                final k_alpha_tau='k_alpha_tau_e-',
                final k_alpha_Qdot='k_alpha_Qdot_e-')))));

      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="cell",
        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-80,-20},{
                80,20}}), graphics),
        Documentation(info="
<html><p>For more information, see the
  <a href=\"modelica://FCSys.Assemblies.Cells.Cell\">Cell</a> model.</p></html>"));
    end CalibratedCell;

    model IntegratedCell "Baseline cell, with integrated CLs and GDLs"
      extends FCSys.BaseClasses.Icons.Cell;
      // Geometric parameters
      parameter Q.Length L_y[:]=fill(1*U.m/1, 1)
        "<html>Lengths along the channel (L<sub>y</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      parameter Q.Length L_z[:]=fill(5*U.mm/1, 1)
        "<html>Lengths across the channel (L<sub>z</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      final parameter Integer n_y=size(L_y, 1)
        "Number of regions along the channel";
      final parameter Integer n_z=size(L_z, 1)
        "Number of regions across the channel";

      // **Update these outputs based on Cell model.
      // Auxiliary variables (for analysis)
      final parameter Q.Length L[Axis]={sum(anFP.L_x) + sum(anCGDL.L_x) + sum(
          pEM.L_x) + sum(caCGDL.L_x) + sum(caFP.L_x),sum(L_y),sum(L_z)} if
        defaults.analysis "Total lengths along the x, y, and z axes";
      final parameter Q.Area A[Axis]={L[cartWrap(ax + 2)]*L[cartWrap(ax + 2)]
          for ax in 1:3} if defaults.analysis "Cross-sectional areas";
      final parameter Q.Volume V=product(L) if defaults.analysis "Volume";
      output Q.Potential Deltav_x[n_y, n_z](each stateSelect=StateSelect.never)
         = caFP.subregions[caFP.n_x, :, :].graphite.'e-'.mu_face[1, 2] - anFP.subregions[
        1, :, :].graphite.'e-'.mu_face[1, 1] if defaults.analysis
        "Electrical potential differences of the segments (x axis)";
      output Q.Power 'Wdot_e-_x'[n_y, n_z](each stateSelect=StateSelect.never)
         = anFP.subregions[1, :, :].graphite.'e-'.Ndot_face[1, 1] .* anFP.subregions[
        1, :, :].graphite.'e-'.mu_face[1, 1] + caFP.subregions[caFP.n_x, :, :].graphite.
        'e-'.Ndot_face[1, 2] .* caFP.subregions[caFP.n_x, :, :].graphite.'e-'.mu_face[
        1, 2] if defaults.analysis
        "Rates of intake of electrical energy of the segments (x axis)";
      output Q.Current I_x[n_y, n_z](each stateSelect=StateSelect.never) = (
        caFP.subregions[caFP.n_x, :, :].graphite.'e-'.Ndot_face[1, 2] - anFP.subregions[
        1, :, :].graphite.'e-'.Ndot_face[1, 1])/2 if defaults.analysis
        "Electrical currents of the segments (x axis)";
      output Q.Current I_x_tot(stateSelect=StateSelect.never) = sum(I_x) if
        defaults.analysis "Total electrical current (x axis)";
      output Q.CurrentAreic Iprimeprime_x[n_y, n_z](each stateSelect=
            StateSelect.never) = {I_x[i_y, i_z]/(L_y[i_y]*L_z[i_z]) for i_z in
        1:n_z, i_y in 1:n_y} if defaults.analysis
        "Areic electrical current of the segments (x axis)";
      output Q.CurrentAreic Iprimeprime_x_avg(stateSelect=StateSelect.never) =
        I_x_tot/A[1] if defaults.analysis
        "Average areic electrical current (x axis)";
      output Q.Current Ndot_H2(stateSelect=StateSelect.never) = sum(anFP.subregions[
        1, :, :].gas.H2.Ndot_face[1, 1]) + sum(anFP.subregions[:, 1, :].gas.H2.Ndot_face[
        2, 1] + anFP.subregions[:, n_y, :].gas.H2.Ndot_face[2, 2]) if defaults.analysis
        "Rate of hydrogen intake";
      output Q.Current Ndot_H2O(stateSelect=StateSelect.never) = sum(anFP.subregions[
        1, :, :].gas.H2O.Ndot_face[1, 1] + caFP.subregions[caFP.n_x, :, :].gas.H2O.Ndot_face[
        1, 2]) + sum(anFP.subregions[:, 1, :].gas.H2O.Ndot_face[1, 1] + anFP.subregions[
        :, n_y, :].gas.H2O.Ndot_face[1, 2]) + sum(caFP.subregions[:, 1, :].gas.H2O.Ndot_face[
        1, 1] - caFP.subregions[:, n_y, :].gas.H2O.Ndot_face[1, 2]) if defaults.analysis
        "Rate of water intake";
      output Q.Current Ndot_O2(stateSelect=StateSelect.never) = sum(caFP.subregions[
        caFP.n_x, :, :].gas.O2.Ndot_face[1, 2]) + sum(caFP.subregions[:, 1, :].gas.O2.Ndot_face[
        1, 1] - caFP.subregions[:, n_y, :].gas.O2.Ndot_face[1, 2]) if defaults.analysis
        "Rate of oxygen intake";
      output Q.Power 'Wdot_e-'(stateSelect=StateSelect.never) = sum('Wdot_e-_x')
         + sum(anFP.subregions[:, 1, :].graphite.'e-'.Ndot_face[2, 1] .* anFP.subregions[
        :, 1, :].graphite.'e-'.mu_face[2, 1] + anFP.subregions[:, n_y, :].graphite.
        'e-'.Ndot_face[2, 2] .* anFP.subregions[:, n_y, :].graphite.'e-'.mu_face[
        2, 2]) + sum(caFP.subregions[:, 1, :].graphite.'e-'.Ndot_face[2, 1] .*
        caFP.subregions[:, 1, :].graphite.'e-'.mu_face[2, 1] + caFP.subregions[
        :, n_y, :].graphite.'e-'.Ndot_face[2, 2] .* caFP.subregions[:, n_y, :].graphite.
        'e-'.mu_face[2, 2]) if defaults.analysis
        "Rate of electrical work (negative for work done)";

      FCSys.Connectors.FaceBus anFPX[n_y, n_z] "Anode plate face" annotation (
          Placement(transformation(extent={{-70,-10},{-50,10}},rotation=0),
            iconTransformation(extent={{-110,-10},{-90,10}})));
      FCSys.Connectors.FaceBus caFPX[n_y, n_z] "Cathode plate face" annotation
        (Placement(transformation(extent={{50,-10},{70,10}}, rotation=0),
            iconTransformation(extent={{90,-10},{110,10}})));
      FCSys.Connectors.FaceBus anFPPosY[anFP.n_x, n_z]
        "Positive anode flow plate face along the y axis" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-40,20}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-40,100})));
      FCSys.Connectors.FaceBus caFPPosY[caFP.n_x, n_z]
        "Positive anode flow plate face along the y axis" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={40,20}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={40,100})));
      FCSys.Connectors.FaceBus caFPNegY[caFP.n_x, n_z]
        "Negative cathode flow plate face along the y axis" annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={40,-20}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={40,-100})));
      FCSys.Connectors.FaceBus anFPNegY[anFP.n_x, n_z]
        "Negative anode flow plate face along the y axis" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-40,-20}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-40,-100})));

      replaceable FCSys.Regions.AnFPs.AnFP anFP(final L_y=L_y, final L_z=L_z)
        "Anode flow plate" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-50,-10},{-30,10}})));

      replaceable Regions.AnCLs.AnCGDL anCGDL(final L_y=L_y, final L_z=L_z)
        "Anode catalyst/gas diffusion layer" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-30,-10},{-10,10}})));

      replaceable FCSys.Regions.PEMs.PEM pEM(final L_y=L_y,final L_z=L_z) "PEM"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-10,-10},{10,10}})));
      replaceable Regions.CaCLs.CaCGDL caCGDL(final L_y=L_y, final L_z=L_z)
        "Cathode catalyst/gas diffusion layer" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{10,-10},{30,10}})));

      replaceable FCSys.Regions.CaFPs.CaFP caFP(final L_y=L_y,final L_z=L_z)
        "Cathode flow plate" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{30,-10},{50,10}})));

    protected
      outer FCSys.BCs.Defaults defaults "Environmental properties and settings";

    equation
      // Internal connections (between layers)
      connect(anFP.xPositive, anCGDL.xNegative) annotation (Line(
          points={{-30,6.10623e-16},{-30,6.10623e-16}},
          color={240,0,0},
          smooth=Smooth.None,
          thickness=0.5));
      connect(anCGDL.xPositive, pEM.xNegative) annotation (Line(
          points={{-10,6.10623e-16},{-10,6.10623e-16}},
          color={240,0,0},
          smooth=Smooth.None,
          thickness=0.5));
      connect(pEM.xPositive, caCGDL.xNegative) annotation (Line(
          points={{10,6.10623e-16},{10,6.10623e-16},{10,6.10623e-16},{10,
              6.10623e-16}},
          color={0,0,240},
          smooth=Smooth.None,
          thickness=0.5));

      connect(caCGDL.xPositive, caFP.xNegative) annotation (Line(
          points={{30,6.10623e-16},{30,6.10623e-16}},
          color={0,0,240},
          smooth=Smooth.None,
          thickness=0.5));

      // External connections
      connect(anFPX, anFP.xNegative) annotation (Line(
          points={{-60,5.55112e-16},{-60,6.10623e-16},{-50,6.10623e-16}},
          color={127,127,127},
          smooth=Smooth.None,
          thickness=0.5));
      connect(caFP.xPositive, caFPX) annotation (Line(
          points={{50,6.10623e-16},{60,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None,
          thickness=0.5));
      connect(anFP.yPositive, anFPPosY) annotation (Line(
          points={{-40,10},{-40,20}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(caFP.yPositive, caFPPosY) annotation (Line(
          points={{40,10},{40,20}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      connect(caFP.yNegative, caFPNegY) annotation (Line(
          points={{40,-10},{40,-20}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      connect(anFPNegY, anFP.yNegative) annotation (Line(
          points={{-40,-20},{-40,-10}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="cell",
        Documentation(info="

<html><p>For more information, see the
  <a href=\"modelica://FCSys.Assemblies.Cells.Cell\">Cell</a> model.</p></html>"),

        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-60,-20},{
                60,20}}), graphics),
        Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.1), graphics={Line(
                  points={{-40,100},{-40,60}},
                  color={255,128,0},
                  visible=inclY,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{-8,-1},{44,-1}},
                  color={0,128,255},
                  visible=inclX,
                  thickness=0.5,
                  origin={39,56},
                  rotation=90),Line(
                  points={{100,0},{56,0}},
                  color={127,127,127},
                  visible=inclX,
                  thickness=0.5),Line(
                  points={{-8,-1},{28,-1}},
                  color={0,128,255},
                  visible=inclX,
                  thickness=0.5,
                  origin={39,-92},
                  rotation=90),Line(
                  points={{-40,-58},{-40,-100}},
                  color={255,128,0},
                  visible=inclY,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{-66,0},{-100,0}},
                  color={127,127,127},
                  visible=inclX,
                  thickness=0.5)}),
        experimentSetupOutput);
    end IntegratedCell;

    model CellSSIC
      "Cell model with steady state initial conditions (no reaction)"

      extends Cell(
        anFP(subregions(each graphite('e-'(initMethPartNum=InitMethScalar.None)))),

        anGDL(subregions(each graphite('e-'(initMethPartNum=InitMethScalar.None)))),

        anCL(subregions(each graphite('e-'(final Ndot_IC=0, final
                  initMethPartNum=InitMethScalar.ReactionRate)), each ionomer(
                'H+'(initMethPartNum=InitMethScalar.None)))),
        pEM(subregions(each ionomer('H+'(final initMethPartNum=InitMethScalar.None)))),

        caCL(subregions(each ionomer('H+'(initMethPartNum=InitMethScalar.Amount)))),

        caGDL(subregions(each graphite('e-'(final initMethPartNum=
                    InitMethScalar.None)))),
        caFP(subregions(each graphite('e-'(final initMethPartNum=InitMethScalar.None)))));

    initial equation
      // Equipotential
      anCL.subregions.graphite.'e-'.mu = anGDL.subregions.graphite.'e-'.mu;
      anFP.subregions.graphite.'e-'.mu = anGDL.subregions.graphite.'e-'.mu;

      anCL.subregions.ionomer.'H+'.mu = pEM.subregions.ionomer.'H+'.mu;
      pEM.subregions.ionomer.'H+'.mu = caCL.subregions.ionomer.'H+'.mu;
      caCL.subregions.graphite.'e-'.mu = caGDL.subregions.graphite.'e-'.mu;
      caFP.subregions.graphite.'e-'.mu = caGDL.subregions.graphite.'e-'.mu;
    end CellSSIC;
  end Cells;
  annotation (Documentation(info="
<html>
  <p>
<b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
Copyright 2007&ndash;2012, Georgia Tech Research Corporation.
</p>
<p>
<i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.ModelicaLicense2\">
FCSys.UsersGuide.ModelicaLicense2</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
http://www.modelica.org/licenses/ModelicaLicense2</a>.</i>
</p>
</html>"));
end Assemblies;
