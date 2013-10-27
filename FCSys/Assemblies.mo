within FCSys;
package Assemblies "Combinations of regions (e.g., cells)"
  extends Modelica.Icons.Package;

  package Cells "Single-cell PEMFC models"
    extends Modelica.Icons.Package;
    extends FCSys.Icons.PackageUnderConstruction;
    package Examples "Examples"

      model Cell
        "<html>Isolated <a href=\"modelica://FCSys.Connectors.PhysicalBus\">PhysicalBus</a> or <a href=\"modelica://FCSys.Assemblies.Cells.Cell\">Cell</a> model</html>"
        extends Modelica.Icons.Example;
        // **Remove this model once Polarization is running.

        inner FCSys.Conditions.Environment environment
          annotation (Placement(transformation(extent={{20,20},{40,40}})));
        replaceable Cells.Cell cell "Fuel cell" annotation (
            __Dymola_choicesFromPackage=true, Placement(transformation(extent={
                  {-10,-10},{10,10}})));
        annotation (
          experiment(StopTime=10, Tolerance=1e-06),
          Commands(file=
                "Resources/Scripts/Dymola/Assemblies.Cells.Examples.Cell.mos"
              "Assemblies.Cells.Examples.Cell.mos"),
          __Dymola_experimentSetupOutput);

      end Cell;
      extends Modelica.Icons.ExamplesPackage;

      model Polarization "Run a cell polarization"
        extends Modelica.Icons.Example;

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
          annotation (Placement(transformation(extent={{30,32},{50,52}})));

        Conditions.TestStands.TestStand testStand(
          zJ=currentDensity.y,
          anInletRH(displayUnit="1") = 0.8,
          caInletRH(displayUnit="1") = 0.5,
          T_an=333.15*U.K,
          T_ca=333.15*U.K,
          anStoich=1.5,
          caStoich=2.0) annotation (__Dymola_choicesFromPackage=true, Placement(
              transformation(extent={{-16,-16},{16,16}})));

        replaceable Cells.Cell cell "Fuel cell" annotation (
            __Dymola_choicesFromPackage=true, Placement(transformation(extent={
                  {-10,-10},{10,10}})));
        Modelica.Blocks.Sources.Ramp currentDensity(
          height=100*U.A/U.cm^2,
          duration=100,
          startTime=0.1)
          annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
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
            points={{-4,16},{-4,14.5},{-4,14.5},{-4,13},{-4,13},{-4,10}},
            color={240,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(testStand.caNegative, cell.caNegative) annotation (Line(
            points={{4,-16},{4,-14.5},{4,-14.5},{4,-13},{4,-13},{4,-10}},
            color={0,0,240},
            thickness=0.5,
            smooth=Smooth.None));
        connect(cell.anNegative, testStand.anNegative) annotation (Line(
            points={{-4,-10},{-4,-11.5},{-4,-11.5},{-4,-13},{-4,-13},{-4,-16}},

            color={240,0,0},
            thickness=0.5,
            smooth=Smooth.None));

        connect(cell.caPositive, testStand.caPositive) annotation (Line(
            points={{4,10},{4,11.5},{4,11.5},{4,13},{4,13},{4,16}},
            color={0,0,240},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (Commands(file=
                "Resources/Scripts/Dymola/Assemblies.Cells.Examples.Polarization.mos"
              "Assemblies.Cells.Examples.Polarization.mos"), Diagram(graphics));
      end Polarization;

      model EIS "Model for electrochemical-impedance spectroscopy"
        import FCSys;
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

    end Examples;

    model Cell "Single-cell PEMFC"
      import FCSys.Utilities.average;
      extends FCSys.Icons.Cell;

      // **Add overall parameter to include or exclude liquid H2O.

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

      Connectors.FaceBus an[n_y, n_z] "Interface with the anode end plate"
        annotation (Placement(transformation(extent={{-100,-20},{-80,0}},
              rotation=0), iconTransformation(extent={{-110,-10},{-90,10}})));
      Connectors.FaceBus ca[n_y, n_z] "Interface with the cathode end plate"
        annotation (Placement(transformation(extent={{60,-20},{80,0}}, rotation
              =0), iconTransformation(extent={{90,-10},{110,10}})));
      Connectors.FaceBus anNegative[anFP.n_x, n_z] "Negative anode fluid port"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-70,-30}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-40,-100})));
      Connectors.FaceBus caNegative[caFP.n_x, n_z]
        "Negative cathode fluid port" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={50,-30}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={40,-100})));
      Connectors.FaceBus anPositive[anFP.n_x, n_z] "Positive anode fluid port"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-70,10}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-40,100})));
      Connectors.FaceBus caPositive[caFP.n_x, n_z]
        "Positive cathode fluid port" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={50,10}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={40,100})));

      replaceable Regions.AnFPs.AnFP anFP(final L_y=L_y, final L_z=L_z)
        "Anode flow plate" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-80,-20},{-60,0}})));
      replaceable Regions.AnGDLs.AnGDL anGDL(final L_y=L_y, final L_z=L_z)
        "Anode gas diffusion layer" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-60,-20},{-40,0}})));
      replaceable Regions.AnCLs.AnCL anCL(final L_y=L_y, final L_z=L_z)
        "Anode catalyst layer" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-40,-20},{-20,0}})));
      replaceable Regions.PEMs.PEM PEM(final L_y=L_y, final L_z=L_z)
        "Proton exchange membrane" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-20,-20},{0,0}})));
      replaceable Regions.CaCLs.CaCL caCL(final L_y=L_y, final L_z=L_z)
        "Cathode catalyst layer" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{0,-20},{20,0}})));
      replaceable Regions.CaGDLs.CaGDL caGDL(final L_y=L_y, final L_z=L_z)
        "Cathode gas diffusion layer" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{20,-20},{40,0}})));
      replaceable Regions.CaFPs.CaFP caFP(final L_y=L_y, final L_z=L_z)
        "Cathode flow plate" annotation (
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
        defaultComponentPrefixes="replaceable",
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

      // **Add overall parameter to include or exclude liquid H2O.

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

      Connectors.FaceBus an[n_y, n_z] "Interface with the anode end plate"
        annotation (Placement(transformation(extent={{-80,-20},{-60,0}},
              rotation=0), iconTransformation(extent={{-110,-10},{-90,10}})));
      Connectors.FaceBus ca[n_y, n_z] "Interface with the cathode end plate"
        annotation (Placement(transformation(extent={{40,-20},{60,0}}, rotation
              =0), iconTransformation(extent={{90,-10},{110,10}})));
      Connectors.FaceBus anNegative[anFP.n_x, n_z] "Negative anode fluid port"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-50,-30}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-40,-100})));
      Connectors.FaceBus caNegative[caFP.n_x, n_z]
        "Negative cathode fluid port" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={30,-30}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={40,-100})));
      Connectors.FaceBus anPositive[anFP.n_x, n_z] "Positive anode fluid port"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-50,10}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-40,100})));
      Connectors.FaceBus caPositive[caFP.n_x, n_z]
        "Positive cathode fluid port" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={30,10}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={40,100})));

      replaceable Regions.AnFPs.AnFP anFP(final L_y=L_y, final L_z=L_z)
        "Anode flow plate" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-60,-20},{-40,0}})));
      FCSys.Regions.AnCLs.AnCGDL anCGDL(final L_y=L_y, final L_z=L_z)
        "Anode catalyst and gas diffusion layer" annotation (Dialog(group=
              "Layers"), Placement(transformation(extent={{-40,-20},{-20,0}})));

      replaceable Regions.PEMs.PEM PEM(final L_y=L_y, final L_z=L_z)
        "Proton exchange membrane" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Layers"),
        Placement(transformation(extent={{-20,-20},{0,0}})));
      FCSys.Regions.CaCLs.CaCGDL caCGDL(final L_y=L_y, final L_z=L_z)
        "Cathode catalyst and gas diffusion layer" annotation (Dialog(group=
              "Layers"), Placement(transformation(extent={{0,-20},{20,0}})));

      replaceable Regions.CaFPs.CaFP caFP(final L_y=L_y, final L_z=L_z)
        "Cathode flow plate" annotation (
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
        defaultComponentPrefixes="replaceable",
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
