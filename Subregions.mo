within FCSys;
package Subregions
  "Control volumes with multi-species transport, exchange, and storage"
  extends Modelica.Icons.Package;
  package Examples "Examples and tests"
    extends Modelica.Icons.ExamplesPackage;
    model SubregionH2 "Test a subregion"
      extends Modelica.Icons.Example;
      parameter Boolean inclC=false "Carbon (C)" annotation (choices(
            __Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclC19HF37O5S=false
        "<html>Nafion sulfonate (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean 'incle-'=false "<html>Electrons (e<sup>-</sup>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclH2=true "<html>Hydrogen (H<sub>2</sub>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean 'inclH+'=false "<html>Protons (H<sup>+</sup>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclN2=false "<html>Nitrogen (N<sub>2</sub>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclO2=false "<html>Oxygen (O<sub>2</sub>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));

      FCSys.Subregions.Subregion subregion(
        L={1,1,1}*U.cm,
        graphite(
          final inclC=inclC,
          final 'incle-'='incle-',
          C(V_IC=subregion.V/4)),
        ionomer(
          final inclC19HF37O5S=inclC19HF37O5S,
          final 'inclH+'='inclH+',
          C19HF37O5S(V_IC=subregion.V/4)),
        inclLinY=false,
        inclLinZ=false,
        inclYFaces=false,
        inclZFaces=false,
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2))
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      inner FCSys.BCs.Defaults defaults(analysis=true)
        annotation (Placement(transformation(extent={{30,20},{50,40}})));
      replaceable FCSys.BCs.Face.SubregionFlow bC1(
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        graphite(final inclC=inclC, final 'incle-'='incle-'),
        ionomer(final inclC19HF37O5S=inclC19HF37O5S, final 'inclH+'='inclH+'))
        constrainedby FCSys.BCs.Face.Subregion annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-20,0})));

      replaceable FCSys.BCs.Face.SubregionFlow bC2(
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        graphite(
          final inclC=inclC,
          final 'incle-'='incle-',
          'e-'(redeclare FCSys.BCs.Face.Species.Material.Pressure materialBC)),

        ionomer(
          final inclC19HF37O5S=inclC19HF37O5S,
          final 'inclH+'='inclH+',
          'H+'(redeclare FCSys.BCs.Face.Species.Material.Pressure materialBC)))
        constrainedby FCSys.BCs.Face.Subregion annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={20,0})));

    equation
      connect(bC2.face, subregion.xPositive) annotation (Line(
          points={{16,1.23436e-15},{15,1.22125e-15},{15,6.10623e-16},{10,
              6.10623e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(bC1.face, subregion.xNegative) annotation (Line(
          points={{-16,3.65701e-16},{-15,1.11022e-15},{-15,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        Placement(transformation(extent={{70,70},{90,90}})),
        Diagram(graphics),
        experiment(StopTime=10),
        experimentSetupOutput,
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionH2.mos"));
    end SubregionH2;

    model SubregionHOR
      "<html>Test a subregion with the hydrogen oxidation reaction and the essential species for it (C, C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S, e<sup>-</sup>, H<sub>2</sub>, and H<sup>+</sup>)</html>"

      extends SubregionH2(
        inclC=true,
        inclC19HF37O5S=true,
        'incle-'=true,
        inclH2=true,
        'inclH+'=true,
        bC2(
          gas(H2(redeclare FCSys.BCs.Face.Species.Material.Pressure materialBC,
                materialSpec(k=1*U.atm))),
          graphite(C(redeclare FCSys.BCs.Face.Species.Entropy.Temperature
                entropyBC, entropySpec(k=defaults.T)), 'e-'(redeclare
                FCSys.BCs.Face.Species.Material.Current materialBC, redeclare
                Modelica.Blocks.Sources.Ramp materialSpec(duration=1000, height
                  =2*U.A))),
          ionomer('H+'(redeclare FCSys.BCs.Face.Species.Material.Pressure
                materialBC, materialSpec(k(start=1*U.atm))))));

      annotation (
        experiment(StopTime=1000, Tolerance=1e-06),
        experimentSetupOutput,
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionHOR.mos"));
    end SubregionHOR;

    model SubregionORR
      "<html>Test a subregion with the oxygen reduction reaction and the essential species for it (C, C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S, e<sup>-</sup>, H<sub>2</sub>O, H<sup>+</sup>, and O<sub>2</sub>)</html>"
      extends SubregionH2(
        inclC=true,
        inclC19HF37O5S=true,
        'incle-'=true,
        inclH2=false,
        inclH2O=true,
        'inclH+'=true,
        inclO2=true,
        bC2(
          gas(H2O(redeclare FCSys.BCs.Face.Species.Material.Pressure materialBC,
                materialSpec(k=1*U.atm)), O2(redeclare
                FCSys.BCs.Face.Species.Material.Pressure materialBC,
                materialSpec(k=1*U.atm))),
          graphite(C(redeclare FCSys.BCs.Face.Species.Entropy.Temperature
                entropyBC, entropySpec(k=defaults.T)), 'e-'(redeclare
                FCSys.BCs.Face.Species.Material.Current materialBC, redeclare
                Modelica.Blocks.Sources.Ramp materialSpec(duration=1000, height
                  =-2*U.A))),
          ionomer('H+'(redeclare FCSys.BCs.Face.Species.Material.Pressure
                materialBC, materialSpec(k(start=1*U.atm))))),
        defaults(analysis=false));

      annotation (
        experiment(StopTime=1000, Tolerance=1e-06),
        experimentSetupOutput,
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionORR.mos"));
    end SubregionORR;

    model SubregionsH2 "Test a one-dimensional array of subregions"
      extends Modelica.Icons.Example;

      parameter Integer n_x=0
        "Number of discrete subregions along the x axis, besides the 2 side subregions";
      parameter Boolean inclC=false "Carbon (C)" annotation (choices(
            __Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclC19HF37O5S=false
        "<html>Nafion sulfonate (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean 'incle-'=false "<html>Electrons (e<sup>-</sup>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclH2=true "<html>Hydrogen (H<sub>2</sub>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean 'inclH+'=false "<html>Protons (H<sup>+</sup>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclN2=false "<html>Nitrogen (N<sub>2</sub>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclO2=false "<html>Oxygen (O<sub>2</sub>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));

      FCSys.Subregions.Subregion subregion1(
        L={1,1,1}*U.cm,
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2O(
            p_IC=1.1*defaults.p,
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false)),
          N2(
            p_IC=1.1*defaults.p,
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false)),
          O2(
            p_IC=1.1*defaults.p,
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false)),
          H2(
            p_IC=1.1*defaults.p,
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false))),
        graphite(
          final inclC=inclC,
          final 'incle-'='incle-',
          C(V_IC=subregion1.V/4),
          'e-'(
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false))),
        ionomer(
          final inclC19HF37O5S=inclC19HF37O5S,
          final 'inclH+'='inclH+',
          C19HF37O5S(V_IC=subregion1.V/4),
          'H+'(
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false))),
        inclLinY=false,
        inclLinZ=false,
        inclYFaces=false,
        inclZFaces=false)
        annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));

      FCSys.Subregions.Subregion subregions[n_x](
        each L={1,1,1}*U.cm,
        gas(
          each final inclH2=inclH2,
          each final inclH2O=inclH2O,
          each final inclN2=inclN2,
          each final inclO2=inclO2,
          H2(
            each p_IC=defaults.p,
            each yNegative(viscousZ=false,viscousX=false),
            each yPositive(viscousZ=false,viscousX=false),
            each zNegative(viscousX=false,viscousY=false),
            each zPositive(viscousX=false,viscousY=false)),
          H2O(
            each p_IC=defaults.p,
            each yNegative(viscousZ=false,viscousX=false),
            each yPositive(viscousZ=false,viscousX=false),
            each zNegative(viscousX=false,viscousY=false),
            each zPositive(viscousX=false,viscousY=false)),
          N2(
            each p_IC=defaults.p,
            each yNegative(viscousZ=false, viscousX=false),
            each yPositive(viscousZ=false,viscousX=false),
            each zNegative(viscousX=false,viscousY=false),
            each zPositive(viscousX=false,viscousY=false)),
          O2(
            each p_IC=defaults.p,
            each yNegative(viscousZ=false,viscousX=false),
            each yPositive(viscousZ=false,viscousX=false),
            each zNegative(viscousX=false,viscousY=false),
            each zPositive(viscousX=false,viscousY=false))),
        graphite(
          each final inclC=inclC,
          each final 'incle-'='incle-',
          C(each V_IC=subregions[1].V/4),
          'e-'(
            each yNegative(viscousZ=false,viscousX=false),
            each yPositive(viscousZ=false,viscousX=false),
            each zNegative(viscousX=false,viscousY=false),
            each zPositive(viscousX=false,viscousY=false))),
        ionomer(
          each final inclC19HF37O5S=inclC19HF37O5S,
          each final 'inclH+'='inclH+',
          C19HF37O5S(each V_IC=subregions[1].V/4),
          'H+'(
            each yNegative(viscousZ=false,viscousX=false),
            each yPositive(viscousZ=false,viscousX=false),
            each zNegative(viscousX=false,viscousY=false),
            each zPositive(viscousX=false,viscousY=false))),
        each inclLinY=false,
        each inclLinZ=false,
        each inclYFaces=false,
        each inclZFaces=false) if n_x > 0
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      FCSys.Subregions.Subregion subregion2(
        L={1,1,1}*U.cm,
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2(
            p_IC=0.9*defaults.p,
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false)),
          H2O(
            p_IC=0.9*defaults.p,
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false)),
          N2(
            p_IC=0.9*defaults.p,
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false)),
          O2(
            p_IC=0.9*defaults.p,
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false))),
        graphite(
          final inclC=inclC,
          final 'incle-'='incle-',
          C(V_IC=subregion2.V/4),
          'e-'(
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false))),
        ionomer(
          final inclC19HF37O5S=inclC19HF37O5S,
          final 'inclH+'='inclH+',
          C19HF37O5S(V_IC=subregion2.V/4),
          'H+'(
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false))),
        inclLinY=false,
        inclLinZ=false,
        inclYFaces=false,
        inclZFaces=false)
        annotation (Placement(transformation(extent={{20,-10},{40,10}})));

      inner FCSys.BCs.Defaults defaults
        annotation (Placement(transformation(extent={{60,20},{80,40}})));
      replaceable FCSys.BCs.Face.Subregion0Current bC1 constrainedby
        FCSys.BCs.Face.Subregion(
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        graphite(final inclC=inclC, final 'incle-'='incle-'),
        ionomer(final inclC19HF37O5S=inclC19HF37O5S, final 'inclH+'='inclH+'))
        annotation (__Dymola_choicesFromPackage=true,Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-50,0})));
      replaceable FCSys.BCs.Face.Subregion0Current bC2 constrainedby
        FCSys.BCs.Face.Subregion(
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        graphite(
          final inclC=inclC,
          final 'incle-'='incle-',
          'e-'(redeclare FCSys.BCs.Face.Species.Material.Pressure materialBC)),

        ionomer(
          final inclC19HF37O5S=inclC19HF37O5S,
          final 'inclH+'='inclH+',
          'H+'(redeclare FCSys.BCs.Face.Species.Material.Pressure materialBC)))
        annotation (__Dymola_choicesFromPackage=true, Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={50,0})));

    equation
      connect(bC1.face, subregion1.xNegative) annotation (Line(
          points={{-46,3.65701e-16},{-40,6.10623e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));
      connect(subregion1.xPositive, subregions[1].xNegative) annotation (Line(
          points={{-20,6.10623e-16},{-16,-3.36456e-22},{-16,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      for i in 1:n_x - 1 loop
        connect(subregions[i].xPositive, subregions[i + 1].xNegative)
          "Not shown on the diagram";
      end for;
      if n_x == 0 then
        connect(subregion1.xPositive, subregion2.xNegative)
          "Not shown on the diagram";
      end if;
      connect(subregions[n_x].xPositive, subregion2.xNegative) annotation (Line(
          points={{10,6.10623e-16},{20,-3.36456e-22},{20,6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(bC2.face, subregion2.xPositive) annotation (Line(
          points={{46,1.23436e-15},{46,6.10623e-16},{40,6.10623e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        Placement(transformation(extent={{70,70},{90,90}})),
        Diagram(graphics),
        experiment(
          StopTime=4,
          NumberOfIntervals=5000,
          Tolerance=1e-06,
          Algorithm="Dassl"),
        experimentSetupOutput,
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionsH2.mos"));
    end SubregionsH2;

    model SubregionsCAndH2
      "<html>Test a one-dimensional array of subregions with C and H<sub>2</sub></html>"
      extends SubregionsH2(inclC=true);
      annotation (
        experiment(StopTime=3, NumberOfIntervals=5000),
        experimentSetupOutput,
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionsCAndH2.mos"));
    end SubregionsCAndH2;

    model SubregionsCAndH2AndH2O
      "<html>Test a one-dimensional array of subregions with C, H<sub>2</sub>, and H<sub>2</sub>O</html>"
      extends SubregionsH2(
        inclC=true,
        inclH2=true,
        inclH2O=true);

      annotation (
        experiment(StopTime=8, Tolerance=1e-06),
        experimentSetupOutput,
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionsCAndH2AndH2O.mos"));
    end SubregionsCAndH2AndH2O;

    model SubregionsCAndN2
      "<html>Test a one-dimensional array of subregions with C and N<sub>2</sub></html>"
      extends SubregionsH2(
        inclH2=false,
        inclC=true,
        inclN2=true);

      annotation (
        experiment(StopTime=40, Tolerance=1e-06),
        experimentSetupOutput,
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionsCAndN2.mos"));
    end SubregionsCAndN2;

    model SubregionsCAndO2
      "<html>Test a one-dimensional array of subregions with C and O<sub>2</sub></html>"
      extends SubregionsH2(
        inclH2=false,
        inclC=true,
        inclO2=true);

      annotation (
        experiment(
          StopTime=25,
          NumberOfIntervals=5000,
          Tolerance=1e-06,
          Algorithm="Dassl"),
        experimentSetupOutput,
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionsCAndO2.mos"));
    end SubregionsCAndO2;

    model SubregionsC19HF37O5SAndH2OAndHplus
      "<html>Test a one-dimensional array of subregions with C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S, H<sub>2</sub>O, and H<sup>+</sup></html>"
      extends SubregionsH2(
        inclH2=false,
        inclH2O=true,
        'inclH+'=true,
        inclC19HF37O5S=true);
      annotation (
        experiment(StopTime=150, Tolerance=1e-06),
        experimentSetupOutput,
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionsC19HF37O5SAndH2OAndHplus.mos"));
    end SubregionsC19HF37O5SAndH2OAndHplus;

    model SubregionsCAndeminus
      "<html>Test a one-dimensional array of subregions with e<sup>-</sup></html>"
      extends SubregionsH2(
        inclC=true,
        inclH2=false,
        'incle-'=true,
        defaults(analysis=false),
        subregion1(setVolume=false),
        subregions(each setVolume=false),
        subregion2(setVolume=false));

      annotation (
        experiment(Tolerance=1e-06),
        experimentSetupOutput,
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionsCAndeminus.mos"));
    end SubregionsCAndeminus;

    model ThermalConduction "Test thermal conduction (through solid)"
      extends SubregionsH2(
        n_x=8,
        inclC=true,
        inclH2=false,
        subregion1(graphite(C(
              initMethPartNum=InitMethScalar.Pressure,
              V_IC=0.99*subregion1.V,
              T_IC=1.1*defaults.T,
              T(displayUnit="degC"),
              beta_S=Modelica.Constants.inf))),
        subregions(graphite(C(
              each initMethPartNum=InitMethScalar.Pressure,
              each V_IC=0.99*subregions[1].V,
              each T(displayUnit="degC")))),
        subregion2(graphite(C(
              initMethPartNum=InitMethScalar.Pressure,
              V_IC=0.99*subregion2.V,
              T(displayUnit="degC")))),
        redeclare FCSys.BCs.Face.Subregion0Current0Power bC1,
        redeclare FCSys.BCs.Face.Subregion0Current0Power bC2);

      annotation (
        Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.ThermalConduction.mos"),

        experiment(StopTime=298.15, Algorithm="Dassl"),
        experimentSetupOutput);

    end ThermalConduction;

    model ThermalConductionConvection
      "Test combined thermal conduction and convection"
      extends SubregionsH2(
        n_x=8,
        inclC=true,
        inclH2=false,
        inclN2=true,
        subregion1(gas(N2(
              T_IC=1.1*defaults.T,
              p_IC=defaults.p,
              phi(displayUnit="mm/s"))), graphite(C(
              V_IC=0.5*subregion1.V,
              T_IC=1.1*defaults.T,
              T(displayUnit="degC"),
              beta_S=Modelica.Constants.inf))),
        subregions(gas(N2(each p_IC=defaults.p, phi(each displayUnit="mm/s"))),
            graphite(C(each V_IC=0.5*subregions[1].V, each T(displayUnit="degC")))),

        subregion2(gas(N2(p_IC=defaults.p, phi(displayUnit="mm/s"))), graphite(
              C(V_IC=0.5*subregion2.V, T(displayUnit="degC")))),
        redeclare FCSys.BCs.Face.Subregion0Current0Power bC1,
        redeclare FCSys.BCs.Face.Subregion0Current0Power bC2);

      annotation (
        Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.ThermalConductionConvection.mos"),

        experiment(StopTime=200, Algorithm="Dassl"),
        experimentSetupOutput);

    end ThermalConductionConvection;

    model ReactionRamp
      "Test chemical reaction with reaction rate ramped over time"

      extends Modelica.Icons.Example;

      FCSys.Subregions.Reaction reaction(n_spec=3)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      FCSys.BCs.Chemical.Species.Species species1(redeclare
          FCSys.Characteristics.'e-'.Graphite Data, materialBC=FCSys.BCs.Chemical.Species.BaseClasses.BCTypeMaterial.PotentialElectrochemicalPerTemperature)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-30,-24})));

      FCSys.BCs.Chemical.Species.Species species2(redeclare
          FCSys.Characteristics.'H+'.Gas Data, materialBC=FCSys.BCs.Chemical.Species.BaseClasses.BCTypeMaterial.PotentialElectrochemicalPerTemperature)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={0,-24})));

      FCSys.BCs.Chemical.Species.Species species3(
        redeclare FCSys.Characteristics.H2.Gas Data,
        materialBC=FCSys.BCs.Chemical.Species.BaseClasses.BCTypeMaterial.Current,

        redeclare Modelica.Blocks.Sources.Ramp materialSpec(height=100*U.A,
            duration=3600e2)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={30,-24})));

      inner FCSys.BCs.Defaults defaults(analysis=true)
        annotation (Placement(transformation(extent={{-90,70},{-70,90}})));
    equation
      connect(species1.chemical, reaction.chemical[1]) annotation (Line(
          points={{-30,-20},{-30,-10},{5.55112e-16,-10},{5.55112e-16,-0.666667}},

          color={208,104,0},
          smooth=Smooth.None));

      connect(species2.chemical, reaction.chemical[2]) annotation (Line(
          points={{-1.11528e-16,-20},{0,-20},{0,5.55112e-16},{5.55112e-16,
              5.55112e-16}},
          color={208,104,0},
          smooth=Smooth.None));

      connect(species3.chemical, reaction.chemical[3]) annotation (Line(
          points={{30,-20},{30,-10},{5.55112e-16,-10},{5.55112e-16,0.666667}},
          color={208,104,0},
          smooth=Smooth.None));

      annotation (
        Diagram(graphics),
        experiment(StopTime=36000),
        experimentSetupOutput,
        Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.ReactionRamp.mos"));
    end ReactionRamp;

    model Reaction "Test an electrochemical reaction"

      extends Modelica.Icons.Example;

      parameter Integer n_lin(
        final min=1,
        final max=3) = 1
        "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>"
        annotation (HideResult=true);

      FCSys.Subregions.Reaction reaction(final n_lin=n_lin, n_spec=3)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      FCSys.BCs.Chemical.Species.Species 'e-'(
        materialBC=FCSys.BCs.Chemical.Species.BaseClasses.BCTypeMaterial.PotentialElectrochemicalPerTemperature,

        redeclare FCSys.Characteristics.'e-'.Graphite Data,
        final n_lin=n_lin) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-30,-24})));

      FCSys.BCs.Chemical.Species.Species 'H+'(
        materialBC=FCSys.BCs.Chemical.Species.BaseClasses.BCTypeMaterial.PotentialElectrochemicalPerTemperature,

        redeclare FCSys.Characteristics.'H+'.Solid Data,
        final n_lin=n_lin) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={0,-24})));

      FCSys.BCs.Chemical.Species.Species H2(
        materialBC=FCSys.BCs.Chemical.Species.BaseClasses.BCTypeMaterial.Current,

        redeclare FCSys.Characteristics.H2.Gas Data,
        final n_lin=n_lin,
        redeclare Modelica.Blocks.Sources.Ramp materialSpec(height=100*U.A,
            duration=100)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={30,-24})));

      inner FCSys.BCs.Defaults defaults(analysis=true)
        annotation (Placement(transformation(extent={{-90,70},{-70,90}})));
    equation
      connect('e-'.chemical, reaction.chemical[1]) annotation (Line(
          points={{-30,-20},{-30,-10},{5.55112e-16,-10},{5.55112e-16,-0.666667}},

          color={208,104,0},
          smooth=Smooth.None));

      connect('H+'.chemical, reaction.chemical[2]) annotation (Line(
          points={{-1.11528e-16,-20},{0,-20},{0,5.55112e-16},{5.55112e-16,
              5.55112e-16}},
          color={208,104,0},
          smooth=Smooth.None));

      connect(H2.chemical, reaction.chemical[3]) annotation (Line(
          points={{30,-20},{30,-10},{5.55112e-16,-10},{5.55112e-16,0.666667}},
          color={208,104,0},
          smooth=Smooth.None));

      annotation (
        Diagram(graphics),
        experiment(StopTime=100),
        experimentSetupOutput,
        Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.Reaction.mos"));
    end Reaction;

    model SpeciesH2 "Test a species"
      extends Modelica.Icons.Example;

      // Geometric parameters
      inner parameter Q.Length L[Axis](each min=Modelica.Constants.small) =
        ones(3)*U.cm "<html>Length (<b>L</b>)</html>"
        annotation (Dialog(group="Geometry"));
      final inner parameter Q.Area A[Axis]={L[cartWrap(ax + 1)]*L[cartWrap(ax
           + 2)] for ax in 1:3} "Cross-sectional area";
      final parameter Q.Volume V=product(L) "Volume";

      inner FCSys.BCs.Defaults defaults(analysis=false)
        annotation (Placement(transformation(extent={{20,20},{40,40}})));
      replaceable Species.H2.Gas.Fixed species constrainedby
        FCSys.Subregions.Species.Species
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      FCSys.BCs.InertDalton.Species inertBC annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=225,
            origin={20,-20})));
    equation
      connect(inertBC.inert, species.inert) annotation (Line(
          points={{17.1716,-17.1716},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      annotation (
        Placement(transformation(extent={{70,70},{90,90}})),
        Diagram(graphics),
        experiment(StopTime=10),
        experimentSetupOutput);
    end SpeciesH2;

    model Specieseminus "Test a species"
      extends Modelica.Icons.Example;
      extends SpeciesH2(redeclare FCSys.Subregions.Species.'e-'.Graphite.Fixed
          species);
      FCSys.BCs.Face.Species.SpeciesX faceBC(redeclare
          FCSys.BCs.Face.Species.Material.Pressure materialBC) annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-24,0})));
    equation
      connect(faceBC.face, species.xNegative) annotation (Line(
          points={{-20,3.65701e-16},{-14,3.65701e-16},{-14,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      annotation (
        Placement(transformation(extent={{70,70},{90,90}})),
        Diagram(graphics),
        experiment(StopTime=10),
        experimentSetupOutput,
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionH2.mos"));
    end Specieseminus;

    model SubregionsCell "Test a one-dimensional array of subregions"
      extends Modelica.Icons.Example;

      parameter Integer n_x=1
        "Number of discrete regions through the PEM along the x axis";

      FCSys.Subregions.Subregion subregion1(
        L={1,1,1}*U.cm,
        gas(
          inclH2=true,
          inclH2O=false,
          inclN2=false,
          inclO2=false,
          H2O(
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false)),
          H2(
            p_IC=1.05*defaults.p,
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false))),
        graphite(
          inclC=true,
          'incle-'=true,
          C(V_IC=subregion1.V/4),
          'e-'(
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false))),
        ionomer(
          inclC19HF37O5S=true,
          'inclH+'=true,
          C19HF37O5S(V_IC=subregion1.V/4),
          'H+'(
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false))),
        inclLinY=false,
        inclLinZ=false,
        inclYFaces=false,
        inclZFaces=false)
        annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));

      FCSys.Subregions.Subregion subregions[n_x](
        each L={1,1,1}*U.cm,
        gas(each inclH2O=false),
        ionomer(
          each inclC19HF37O5S=true,
          each 'inclH+'=true,
          C19HF37O5S(each V_IC=subregions[1].V/4),
          'H+'(
            each yNegative(viscousZ=false,viscousX=false),
            each yPositive(viscousZ=false,viscousX=false),
            each zNegative(viscousX=false,viscousY=false),
            each zPositive(viscousX=false,viscousY=false))),
        each inclLinY=false,
        each inclLinZ=false,
        each inclYFaces=false,
        each inclZFaces=false,
        each setVolume=false) if n_x > 0
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      FCSys.Subregions.Subregion subregion2(
        L={1,1,1}*U.cm,
        gas(
          inclH2=false,
          inclH2O=true,
          inclN2=false,
          inclO2=true,
          H2O(
            p_IC=0.95*defaults.p,
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false)),
          N2(
            p_IC=0.95*defaults.p,
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false)),
          O2(
            p_IC=0.95*defaults.p,
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false))),
        graphite(
          inclC=true,
          'incle-'=true,
          C(V_IC=subregion2.V/4),
          'e-'(
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false))),
        ionomer(
          inclC19HF37O5S=true,
          'inclH+'=true,
          C19HF37O5S(V_IC=subregion2.V/4),
          'H+'(
            yNegative(viscousZ=false, viscousX=false),
            yPositive(viscousZ=false, viscousX=false),
            zNegative(viscousX=false, viscousY=false),
            zPositive(viscousX=false, viscousY=false))),
        inclLinY=false,
        inclLinZ=false,
        inclYFaces=false,
        inclZFaces=false)
        annotation (Placement(transformation(extent={{20,-10},{40,10}})));

      inner FCSys.BCs.Defaults defaults(analysis=true)
        annotation (Placement(transformation(extent={{60,20},{80,40}})));
      replaceable FCSys.BCs.Face.SubregionFlow bC1 constrainedby
        FCSys.BCs.Face.Subregion(gas(
          inclH2=true,
          inclH2O=false,
          inclN2=false,
          inclO2=false), graphite(
          inclC=true,
          'incle-'=true,
          C(redeclare FCSys.BCs.Face.Species.Entropy.Temperature entropyBC,
              entropySpec(k=defaults.T)),
          'e-'(redeclare FCSys.BCs.Face.Species.Material.Current materialBC,
              redeclare Modelica.Blocks.Sources.Ramp materialSpec(duration=1000,
                height=-2*U.A)))) annotation (__Dymola_choicesFromPackage=true,
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-50,0})));
      replaceable FCSys.BCs.Face.SubregionFlow bC2 constrainedby
        FCSys.BCs.Face.Subregion(gas(
          inclH2=false,
          inclH2O=true,
          inclN2=false,
          inclO2=true), graphite(
          inclC=true,
          'incle-'=true,
          C(redeclare FCSys.BCs.Face.Species.Entropy.Temperature entropyBC,
              entropySpec(k=defaults.T)),
          'e-'(redeclare FCSys.BCs.Face.Species.Material.Current materialBC,
              redeclare Modelica.Blocks.Sources.Ramp materialSpec(duration=1000,
                height=2*U.A)))) annotation (__Dymola_choicesFromPackage=true,
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={50,0})));

      replaceable FCSys.BCs.Face.SubregionFlow ground(graphite('incle-'=true,
            'e-'(redeclare FCSys.BCs.Face.Species.Material.Pressure materialBC,
              materialSpec(k=0)))) constrainedby FCSys.BCs.Face.Subregion
        annotation (__Dymola_choicesFromPackage=true, Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-50,30})));

    equation
      connect(bC1.face, subregion1.xNegative) annotation (Line(
          points={{-46,3.65701e-16},{-40,6.10623e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));
      connect(subregion1.xPositive, subregions[1].xNegative) annotation (Line(
          points={{-20,6.10623e-16},{-16,-3.36456e-22},{-16,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      for i in 1:n_x - 1 loop
        connect(subregions[i].xPositive, subregions[i + 1].xNegative)
          "Not shown on the diagram";
      end for;
      if n_x == 0 then
        connect(subregion1.xPositive, subregion2.xNegative)
          "Not shown on the diagram";
      end if;
      connect(subregions[n_x].xPositive, subregion2.xNegative) annotation (Line(
          points={{10,6.10623e-16},{20,-3.36456e-22},{20,6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(bC2.face, subregion2.xPositive) annotation (Line(
          points={{46,1.23436e-15},{46,6.10623e-16},{40,6.10623e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(ground.face, subregion1.xNegative) annotation (Line(
          points={{-46,30},{-40,30},{-40,6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        Placement(transformation(extent={{70,70},{90,90}})),
        Diagram(graphics),
        experiment(
          StopTime=4,
          NumberOfIntervals=5000,
          Tolerance=1e-06,
          Algorithm="Dassl"),
        experimentSetupOutput,
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionsH2.mos"));
    end SubregionsCell;
  end Examples;

  model Subregion "Subregion with all phases included"

    parameter Boolean inclReact=true "Include reaction(s), as appropriate"
      annotation (
      Evaluate=true,
      HideResult=true,
      Dialog(tab="Assumptions"),
      choices(__Dymola_checkBox=true));
    // Note:  This is listed above the extends clause so that it is listed
    // first in the parameter dialog.
    extends BaseClasses.PartialSubregion;

    replaceable Phases.Gas gas(inclH2O=true, final inclLin={inclLinX,inclLinY,
          inclLinZ}) "Gas" annotation (Dialog(group="Phases"), Placement(
          transformation(extent={{-10,-10},{10,10}})));

    replaceable Phases.Graphite graphite(final inclLin={inclLinX,inclLinY,
          inclLinZ}) "Graphite" annotation (Dialog(group="Phases"), Placement(
          transformation(extent={{-10,-10},{10,10}})));

    replaceable Phases.Ionomer ionomer(final inclLin={inclLinX,inclLinY,
          inclLinZ}) "Ionomer" annotation (Dialog(group="Phases"), Placement(
          transformation(extent={{-10,-10},{10,10}})));
    /*
  replaceable Phases.Liquid liquid(final inclLin={inclLinX,
        inclLinY,inclLinZ})  "Liquid" annotation (

    Dialog(group="Phases"),
    Placement(transformation(extent={{-10,-10},{10,10}})));
  */

    FCSys.Subregions.Reaction HOR(final n_lin=n_lin, n_spec=3) if inclReact
       and (graphite.'incle-' and ionomer.'inclH+' and gas.inclH2 and not (gas.inclO2
       and gas.inclH2O)) "Hydrogen oxidation reaction"
      annotation (Placement(transformation(extent={{-50,30},{-30,50}})));
    FCSys.Subregions.Reaction ORR(final n_lin=n_lin, n_spec=4) if inclReact
       and (graphite.'incle-' and ionomer.'inclH+' and gas.inclO2 and gas.inclH2O
       and not gas.inclH2) "Oxygen reduction reaction"
      annotation (Placement(transformation(extent={{-50,30},{-30,50}})));

  protected
    Connectors.ChemicalBusInternal chemical
      "Internal connector to route electrochemical interactions"
      annotation (Placement(transformation(extent={{-30,10},{-10,30}})));

  equation
    // Chemical interactions
    connect(HOR.chemical[1], chemical.'H+') annotation (Line(
        points={{-40,39.3333},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(HOR.chemical[2], chemical.'e-') annotation (Line(
        points={{-40,40},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(HOR.chemical[3], chemical.H2) annotation (Line(
        points={{-40,40.6667},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));

    connect(ORR.chemical[1], chemical.'e-') annotation (Line(
        points={{-40,39.25},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(ORR.chemical[2], chemical.'H+') annotation (Line(
        points={{-40,39.75},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(ORR.chemical[3], chemical.O2) annotation (Line(
        points={{-40,40.25},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(ORR.chemical[4], chemical.H2O) annotation (Line(
        points={{-40,40.75},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));

    // Gas
    connect(gas.chemical, chemical) annotation (Line(
        points={{-5,5},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(gas.inert, volume.inert) annotation (Line(
        points={{8,-8},{11,-11}},
        color={72,90,180},
        smooth=Smooth.None));
    connect(gas.xNegative, xNegative.gas) annotation (Line(
        points={{-8,6.10623e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(gas.xPositive, xPositive.gas) annotation (Line(
        points={{8,6.10623e-16},{24,-4.87687e-22},{24,5.55112e-16},{40,
            5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(gas.yNegative, yNegative.gas) annotation (Line(
        points={{6.10623e-16,-8.4},{-4.87687e-22,-8.4},{-4.87687e-22,-40},{
            5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.yPositive, yPositive.gas) annotation (Line(
        points={{6.10623e-16,10},{5.55112e-16,10},{5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(gas.zNegative, zNegative.gas) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.zPositive, zPositive.gas) annotation (Line(
        points={{-8,-8},{-20,-20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    // Graphite
    connect(graphite.chemical, chemical) annotation (Line(
        points={{-5,5},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));

    connect(graphite.inert, volume.inert) annotation (Line(
        points={{8,-8},{11,-11}},
        color={72,90,180},
        smooth=Smooth.None));
    connect(graphite.xNegative, xNegative.graphite) annotation (Line(
        points={{-8,6.10623e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(graphite.xPositive, xPositive.graphite) annotation (Line(
        points={{8,6.10623e-16},{24,-4.87687e-22},{24,5.55112e-16},{40,
            5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(graphite.yNegative, yNegative.graphite) annotation (Line(
        points={{6.10623e-16,-8.4},{-4.87687e-22,-8.4},{-4.87687e-22,-40},{
            5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.yPositive, yPositive.graphite) annotation (Line(
        points={{6.10623e-16,10},{5.55112e-16,10},{5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(graphite.zNegative, zNegative.graphite) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.zPositive, zPositive.graphite) annotation (Line(
        points={{-8,-8},{-20,-20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    // Ionomer
    connect(ionomer.chemical, chemical) annotation (Line(
        points={{-5,5},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));

    connect(ionomer.inert, volume.inert) annotation (Line(
        points={{8,-8},{11,-11}},
        color={72,90,180},
        smooth=Smooth.None));
    connect(ionomer.xNegative, xNegative.ionomer) annotation (Line(
        points={{-8,6.10623e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(ionomer.xPositive, xPositive.ionomer) annotation (Line(
        points={{8,6.10623e-16},{24,-4.87687e-22},{24,5.55112e-16},{40,
            5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(ionomer.yNegative, yNegative.ionomer) annotation (Line(
        points={{6.10623e-16,-8.4},{-4.87687e-22,-8.4},{-4.87687e-22,-40},{
            5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.yPositive, yPositive.ionomer) annotation (Line(
        points={{6.10623e-16,10},{5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.zNegative, zNegative.ionomer) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.zPositive, zPositive.ionomer) annotation (Line(
        points={{-8,-8},{-20,-20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    /*
  // Liquid
  connect(liquid.inert, volume.interaction) annotation (Line(
      points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
      color={0,180,0},
      smooth=Smooth.None,
      thickness=0.5));
  connect(liquid.xNegative, xNegative.liquid) annotation (Line(
      points={{-10,6.10623e-16},{-10,1.16573e-15},{-25,1.16573e-15},{-25,5.55112e-16},
          {-40,5.55112e-16}},
      color={127,127,127},
      pattern=LinePattern.None,
      thickness=0.5,
      smooth=Smooth.None));
  connect(liquid.xPositive, xPositive.liquid) annotation (Line(
      points={{10,6.10623e-16},{10,5.55112e-16},{40,5.55112e-16}},
      color={127,127,127},
      pattern=LinePattern.None,
      thickness=0.5,
      smooth=Smooth.None));
  connect(liquid.yNegative, yNegative.liquid) annotation (Line(
      points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
      color={127,127,127},
      pattern=LinePattern.None,
      thickness=0.5,
      smooth=Smooth.None));
  connect(liquid.yPositive, yPositive.liquid) annotation (Line(
      points={{6.10623e-16,10},{-4.87687e-22,10},{-4.87687e-22,40},{5.55112e-16,
          40}},
      color={127,127,127},
      pattern=LinePattern.None,
      thickness=0.5,
      smooth=Smooth.None));
  connect(liquid.zNegative, zNegative.liquid) annotation (Line(
      points={{5,5},{20,20}},
      color={127,127,127},
      pattern=LinePattern.None,
      thickness=0.5,
      smooth=Smooth.None));
  connect(liquid.zPositive, zPositive.liquid) annotation (Line(
      points={{-5,-5},{-20,-20}},
      color={127,127,127},
      pattern=LinePattern.None,
      thickness=0.5,
      smooth=Smooth.None));
  */

    annotation (
      defaultComponentPrefixes="replaceable",
      Documentation(info="<html><p>Notes:<ul>
    <li>The <code>inclReact</code> parameter may be set to
    <code>false</code>
    to eliminate unnecessary equations.</li>
  <li>H<sub>2</sub>O is included by default, since at least one species
  must be included.</li></ul>
</p><p>For more information, see the
   <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSubregion\">PartialSubregion</a> model.</p></html>"),

      Icon(graphics),
      Diagram(graphics));
  end Subregion;

  model SubregionNoGraphite "Subregion with all phases except graphite"
    parameter Boolean inclReact=false "Include reaction(s), as appropriate"
      annotation (
      Evaluate=true,
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(tab="Assumptions"));
    // Note:  This is listed above the extension clause so that it is listed
    // first in the parameter dialog.
    extends BaseClasses.PartialSubregion;

    replaceable Phases.Gas gas(
      inclH2O=true,
      inclReact=inclReact,
      final inclLin={inclLinX,inclLinY,inclLinZ}) "Gas" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(group="Phases"),
      Placement(transformation(extent={{-10,-10},{10,10}})));

    replaceable Phases.Ionomer ionomer(final inclLin={inclLinX,inclLinY,
          inclLinZ}) "Ionomer" annotation (Dialog(group="Phases"), Placement(
          transformation(extent={{-10,-10},{10,10}})));

  equation
    // Gas
    connect(gas.inert, volume.inert) annotation (Line(
        points={{8,-8},{11,-11}},
        color={72,90,180},
        smooth=Smooth.None));
    connect(gas.xNegative, xNegative.gas) annotation (Line(
        points={{-8,6.10623e-16},{-8,1.16573e-15},{-25,1.16573e-15},{-25,
            5.55112e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.xPositive, xPositive.gas) annotation (Line(
        points={{8,6.10623e-16},{8,5.55112e-16},{40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.yNegative, yNegative.gas) annotation (Line(
        points={{6.10623e-16,-8.4},{6.10623e-16,-40},{5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.yPositive, yPositive.gas) annotation (Line(
        points={{6.10623e-16,10},{-4.87687e-22,10},{-4.87687e-22,40},{
            5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.zNegative, zNegative.gas) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.zPositive, zPositive.gas) annotation (Line(
        points={{-8,-8},{-20,-20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    // Ionomer
    connect(ionomer.inert, volume.inert) annotation (Line(
        points={{8,-8},{11,-11}},
        color={72,90,180},
        smooth=Smooth.None));
    connect(ionomer.xNegative, xNegative.ionomer) annotation (Line(
        points={{-8,6.10623e-16},{-8,1.16573e-15},{-25,1.16573e-15},{-25,
            5.55112e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.xPositive, xPositive.ionomer) annotation (Line(
        points={{8,6.10623e-16},{8,5.55112e-16},{40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.yNegative, yNegative.ionomer) annotation (Line(
        points={{6.10623e-16,-8.4},{6.10623e-16,-40},{5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.yPositive, yPositive.ionomer) annotation (Line(
        points={{6.10623e-16,10},{-4.87687e-22,10},{-4.87687e-22,40},{
            5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.zNegative, zNegative.ionomer) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.zPositive, zPositive.ionomer) annotation (Line(
        points={{-8,-8},{-20,-20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    annotation (
      defaultComponentPrefixes="replaceable",
      Documentation(info="<html><p>Notes:<ul>
    <li>The <code>inclReact</code> parameter may be set to
    <code>false</code>
    to eliminate unnecessary equations.</li>
  <li>H<sub>2</sub>O is included by default, since at least one species
  must be included.</li></ul>
</p><p>For more information, see the
   <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSubregion\">PartialSubregion</a> model.</p></html>"),

      Icon(graphics),
      Diagram(graphics));
  end SubregionNoGraphite;

  model SubregionNoIonomer "Subregion with all phases except ionomer"
    parameter Boolean inclReact=false "Include reaction(s), as appropriate"
      annotation (
      Evaluate=true,
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(tab="Assumptions"));
    // Note:  This is listed above the extends clause so that it is listed
    // first in the parameter dialog.
    extends BaseClasses.PartialSubregion;

    replaceable Phases.Gas gas(
      inclH2O=true,
      inclReact=inclReact,
      final inclLin={inclLinX,inclLinY,inclLinZ}) "Gas" annotation (Dialog(
          group="Phases"), Placement(transformation(extent={{-10,-10},{10,10}})));

    replaceable Phases.Graphite graphite(final inclLin={inclLinX,inclLinY,
          inclLinZ}) "Graphite" annotation (Dialog(group="Phases"), Placement(
          transformation(extent={{-10,-10},{10,10}})));

  equation
    // Gas
    connect(gas.inert, volume.inert) annotation (Line(
        points={{8,-8},{11,-11}},
        color={72,90,180},
        smooth=Smooth.None));
    connect(gas.xNegative, xNegative.gas) annotation (Line(
        points={{-8,6.10623e-16},{-8,1.16573e-15},{-25,1.16573e-15},{-25,
            5.55112e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.xPositive, xPositive.gas) annotation (Line(
        points={{8,6.10623e-16},{8,5.55112e-16},{40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.yNegative, yNegative.gas) annotation (Line(
        points={{6.10623e-16,-8.4},{6.10623e-16,-40},{5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.yPositive, yPositive.gas) annotation (Line(
        points={{6.10623e-16,10},{-4.87687e-22,10},{-4.87687e-22,40},{
            5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.zNegative, zNegative.gas) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.zPositive, zPositive.gas) annotation (Line(
        points={{-8,-8},{-20,-20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    // Graphite
    connect(graphite.inert, volume.inert) annotation (Line(
        points={{8,-8},{11,-11}},
        color={72,90,180},
        smooth=Smooth.None));
    connect(graphite.xNegative, xNegative.graphite) annotation (Line(
        points={{-8,6.10623e-16},{-8,1.16573e-15},{-25,1.16573e-15},{-25,
            5.55112e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.xPositive, xPositive.graphite) annotation (Line(
        points={{8,6.10623e-16},{8,5.55112e-16},{40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.yNegative, yNegative.graphite) annotation (Line(
        points={{6.10623e-16,-8.4},{6.10623e-16,-40},{5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.yPositive, yPositive.graphite) annotation (Line(
        points={{6.10623e-16,10},{-4.87687e-22,10},{-4.87687e-22,40},{
            5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.zNegative, zNegative.graphite) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.zPositive, zPositive.graphite) annotation (Line(
        points={{-8,-8},{-20,-20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    annotation (
      defaultComponentPrefixes="replaceable",
      Documentation(info="<html><p>Notes:<ul>
    <li>The <code>inclReact</code> parameter may be set to
    <code>false</code>
    to eliminate unnecessary equations.</li>
  <li>H<sub>2</sub>O is included by default, since at least one species
  must be included.</li></ul>
</p><p>For more information, see the
   <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSubregion\">PartialSubregion</a> model.</p></html>"),

      Icon(graphics),
      Diagram(graphics));
  end SubregionNoIonomer;

  package Phases "Phases or mixtures of species"
    extends Modelica.Icons.Package;
    model Phase "Phase with all species conditionally included"
      extends BaseClasses.NullPhase(final n_spec=countTrue({inclC,
            inclC19HF37O5S,'incle-',inclH2,inclH2O,'inclH+',inclN2,inclO2}),
          common(phi(fixed=if inclC or inclC19HF37O5S then fill(false, n_lin)
                 else {initXVel,initYVel,initZVel}[cartAxes])));

      parameter Boolean inclReact=true "Include reaction(s), as appropriate"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(__Dymola_descriptionLabel=true));

      // Conditionally include species.
      parameter Boolean inclC=false "Carbon (C)" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.C.Graphite.Fixed C if inclC constrainedby
        Species.Species(
        final k=k,
        final inclLin=inclLin,
        initMethTemp=if initTemp and reduceFinal then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceFinal then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceFinal then StateSelect.default else StateSelect.prefer))
        "C model" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclC),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      parameter Boolean inclC19HF37O5S=false
        "<html>Nafion sulfonate (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.C19HF37O5S.Solid.Fixed C19HF37O5S if inclC19HF37O5S
        constrainedby Species.Species(
        final k=k,
        final inclLin=inclLin,
        initMethTemp=if initTemp and reduceFinal then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceFinal then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceFinal then StateSelect.default else StateSelect.prefer))
        "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S model</html>"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclC19HF37O5S),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      parameter Boolean 'incle-'=false "<html>Electrons (e<sup>-</sup>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.'e-'.Graphite.Fixed 'e-' if 'incle-' constrainedby
        Species.Species(
        final k=k,
        final inclLin=inclLin,
        initMethX=if (initXVel or inclC or inclC19HF37O5S) and reduceFinal
             then InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethY=if (initYVel or inclC or inclC19HF37O5S) and reduceFinal
             then InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethZ=if (initZVel or inclC or inclC19HF37O5S) and reduceFinal
             then InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethTemp=if initTemp and reduceFinal then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceFinal then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceFinal then StateSelect.default else StateSelect.prefer))
        "<html>'e-' model</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable='incle-'),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      parameter Boolean inclH2=false "<html>Hydrogen (H<sub>2</sub>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.H2.Gas.Fixed H2 if inclH2 constrainedby
        Species.Species(
        final k=k,
        final inclLin=inclLin,
        initMethX=if (initXVel or inclC or inclC19HF37O5S) and reduceFinal
             then InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethY=if (initYVel or inclC or inclC19HF37O5S) and reduceFinal
             then InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethZ=if (initZVel or inclC or inclC19HF37O5S) and reduceFinal
             then InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethTemp=if initTemp and reduceFinal then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceFinal then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceFinal then StateSelect.default else StateSelect.prefer))
        "<html>H<sub>2</sub> model</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclH2),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.H2O.Gas.Fixed H2O if inclH2O constrainedby
        Species.Species(
        final k=k,
        final inclLin=inclLin,
        initMethX=if (initXVel or inclC or inclC19HF37O5S) and reduceFinal
             then InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethY=if (initYVel or inclC or inclC19HF37O5S) and reduceFinal
             then InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethZ=if (initZVel or inclC or inclC19HF37O5S) and reduceFinal
             then InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethTemp=if initTemp and reduceFinal then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceFinal then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceFinal then StateSelect.default else StateSelect.prefer))
        "<html>H<sub>2</sub>O model</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclH2O),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      parameter Boolean 'inclH+'=false "<html>Protons (H<sup>+</sup>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.'H+'.Solid.Fixed 'H+' if 'inclH+' constrainedby
        Species.Species(
        final k=k,
        final inclLin=inclLin,
        initMethX=if initXVel and reduceFinal then InitMethLinMom.None else
            InitMethLinMom.Velocity,
        initMethY=if initYVel and reduceFinal then InitMethLinMom.None else
            InitMethLinMom.Velocity,
        initMethZ=if initZVel and reduceFinal then InitMethLinMom.None else
            InitMethLinMom.Velocity,
        initMethTemp=if initTemp and reduceFinal then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceFinal then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceFinal then StateSelect.default else StateSelect.prefer))
        "<html>H<sup>+</sup> model</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable='inclH+'),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      parameter Boolean inclN2=false "<html>Nitrogen (N<sub>2</sub>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.N2.Gas.Fixed N2 if inclN2 constrainedby
        Species.Species(
        final k=k,
        final inclLin=inclLin,
        initMethX=if (initXVel or inclC or inclC19HF37O5S) and reduceFinal
             then InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethY=if (initYVel or inclC or inclC19HF37O5S) and reduceFinal
             then InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethZ=if (initZVel or inclC or inclC19HF37O5S) and reduceFinal
             then InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethTemp=if initTemp and reduceFinal then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceFinal then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceFinal then StateSelect.default else StateSelect.prefer))
        "<html>N<sub>2</sub> model</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclN2),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      parameter Boolean inclO2=false "<html>Oxygen (O<sub>2</sub>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.O2.Gas.Fixed O2 if inclO2 constrainedby
        Species.Species(
        final k=k,
        final inclLin=inclLin,
        initMethX=if (initXVel or inclC or inclC19HF37O5S) and reduceFinal
             then InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethY=if (initYVel or inclC or inclC19HF37O5S) and reduceFinal
             then InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethZ=if (initZVel or inclC or inclC19HF37O5S) and reduceFinal
             then InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethTemp=if initTemp and reduceFinal then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceFinal then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceFinal then StateSelect.default else StateSelect.prefer))
        "<html>O<sub>2</sub> model</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclO2),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      Reaction '2e-+2H+=H2'(final n_lin=n_lin) if inclReact and (('incle-' and
        'inclH+') or ('incle-' and (inclH2 or (inclO2 and inclH2O))) or (
        'inclH+' and (inclH2 or (inclO2 and inclH2O))))
        "Hydrogen oxidation reaction"
        annotation (Placement(transformation(extent={{-50,30},{-30,50}})));
      Reaction '2H2+O2=2H2O'(final n_lin=n_lin) if inclReact and ((inclH2 and
        inclO2) or (inclH2 and inclH2O) or (inclO2 and inclH2O))
        annotation (Placement(transformation(extent={{-50,30},{-30,50}})));

    equation
      // Chemical interactions
      connect('2e-+2H+=H2'.chemical[1], chemical.'e-') annotation (Line(
          points={{-40,40},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect('2e-+2H+=H2'.chemical[2], chemical.'H+') annotation (Line(
          points={{-40,40},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect('2e-+2H+=H2'.chemical[3], chemical.H2) annotation (Line(
          points={{-40,40},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect('2H2+O2=2H2O'.chemical[1], chemical.H2) annotation (Line(
          points={{-40,40},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect('2H2+O2=2H2O'.chemical[2], chemical.O2) annotation (Line(
          points={{-40,40},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect('2H2+O2=2H2O'.chemical[3], chemical.H2O) annotation (Line(
          points={{-40,40},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));

      // C
      // -
      // Exchange
      connect(C.chemical, chemical.C) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(C.common, common) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(C.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(C.xNegative.material, xNegative.C.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.xPositive.material, xPositive.C.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.yNegative.material, yNegative.C.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.yPositive.material, yPositive.C.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.zNegative.material, zNegative.C.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.zPositive.material, zPositive.C.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(C.xNegative.momentumY, xNegative.C.momentumY) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.xPositive.momentumY, xPositive.C.momentumY) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.yNegative.momentumZ, yNegative.C.momentumZ) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.yPositive.momentumZ, yPositive.C.momentumZ) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.zNegative.momentumX, zNegative.C.momentumX) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.zPositive.momentumX, zPositive.C.momentumX) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(C.xNegative.momentumZ, xNegative.C.momentumZ) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.xPositive.momentumZ, xPositive.C.momentumZ) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.yNegative.momentumX, yNegative.C.momentumX) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.yPositive.momentumX, yPositive.C.momentumX) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.zNegative.momentumY, zNegative.C.momentumY) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.zPositive.momentumY, zPositive.C.momentumY) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(C.xNegative.entropy, xNegative.C.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.xPositive.entropy, xPositive.C.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.yNegative.entropy, yNegative.C.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.yPositive.entropy, yPositive.C.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.zNegative.entropy, zNegative.C.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.zPositive.entropy, zPositive.C.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Note:  It is necessary to connect the subconnectors individually
      // since they are conditional.  Currently (Dymola 7.4 and Modelica 3.2),
      // it is necessary to use conditional subconnectors.  It is not possible
      // to refer to components of expandable connectors using connect
      // statements that are conditional within the equation section (e.g.,
      // if ... then connect(C.xNegative.material, xNegative.C.material);
      // end if;).

      // C19HF37O5S
      // ----------
      // Exchange
      connect(C19HF37O5S.chemical, chemical.C19HF37O5S) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(C19HF37O5S.common, common) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(C19HF37O5S.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(C19HF37O5S.xNegative.material, xNegative.C19HF37O5S.material)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.xPositive.material, xPositive.C19HF37O5S.material)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.yNegative.material, yNegative.C19HF37O5S.material)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.yPositive.material, yPositive.C19HF37O5S.material)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.zNegative.material, zNegative.C19HF37O5S.material)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.zPositive.material, zPositive.C19HF37O5S.material)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(C19HF37O5S.xNegative.momentumY, xNegative.C19HF37O5S.momentumY)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.xPositive.momentumY, xPositive.C19HF37O5S.momentumY)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.yNegative.momentumZ, yNegative.C19HF37O5S.momentumZ)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.yPositive.momentumZ, yPositive.C19HF37O5S.momentumZ)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.zNegative.momentumX, zNegative.C19HF37O5S.momentumX)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.zPositive.momentumX, zPositive.C19HF37O5S.momentumX)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(C19HF37O5S.xNegative.momentumZ, xNegative.C19HF37O5S.momentumZ)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.xPositive.momentumZ, xPositive.C19HF37O5S.momentumZ)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.yNegative.momentumX, yNegative.C19HF37O5S.momentumX)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.yPositive.momentumX, yPositive.C19HF37O5S.momentumX)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.zNegative.momentumY, zNegative.C19HF37O5S.momentumY)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.zPositive.momentumY, zPositive.C19HF37O5S.momentumY)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(C19HF37O5S.xNegative.entropy, xNegative.C19HF37O5S.entropy)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.xPositive.entropy, xPositive.C19HF37O5S.entropy)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.yNegative.entropy, yNegative.C19HF37O5S.entropy)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.yPositive.entropy, yPositive.C19HF37O5S.entropy)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.zNegative.entropy, zNegative.C19HF37O5S.entropy)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.zPositive.entropy, zPositive.C19HF37O5S.entropy)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // e-
      // --
      // Exchange
      connect('e-'.chemical, chemical.'e-') annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect('e-'.common, common) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect('e-'.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect('e-'.xNegative.material, xNegative.'e-'.material) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.xPositive.material, xPositive.'e-'.material) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.yNegative.material, yNegative.'e-'.material) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.yPositive.material, yPositive.'e-'.material) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.zNegative.material, zNegative.'e-'.material) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.zPositive.material, zPositive.'e-'.material) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect('e-'.xNegative.momentumY, xNegative.'e-'.momentumY) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.xPositive.momentumY, xPositive.'e-'.momentumY) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.yNegative.momentumZ, yNegative.'e-'.momentumZ) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.yPositive.momentumZ, yPositive.'e-'.momentumZ) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.zNegative.momentumX, zNegative.'e-'.momentumX) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.zPositive.momentumX, zPositive.'e-'.momentumX) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect('e-'.xNegative.momentumZ, xNegative.'e-'.momentumZ) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.xPositive.momentumZ, xPositive.'e-'.momentumZ) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.yNegative.momentumX, yNegative.'e-'.momentumX) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.yPositive.momentumX, yPositive.'e-'.momentumX) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.zNegative.momentumY, zNegative.'e-'.momentumY) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.zPositive.momentumY, zPositive.'e-'.momentumY) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect('e-'.xNegative.entropy, xNegative.'e-'.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.xPositive.entropy, xPositive.'e-'.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.yNegative.entropy, yNegative.'e-'.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.yPositive.entropy, yPositive.'e-'.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.zNegative.entropy, zNegative.'e-'.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.zPositive.entropy, zPositive.'e-'.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // H2
      // --
      // Exchange
      connect(H2.chemical, chemical.H2) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(H2.common, common) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(H2.xNegative.material, xNegative.H2.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.xPositive.material, xPositive.H2.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.yNegative.material, yNegative.H2.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.yPositive.material, yPositive.H2.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.zNegative.material, zNegative.H2.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.zPositive.material, zPositive.H2.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(H2.xNegative.momentumY, xNegative.H2.momentumY) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.xPositive.momentumY, xPositive.H2.momentumY) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.yNegative.momentumZ, yNegative.H2.momentumZ) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.yPositive.momentumZ, yPositive.H2.momentumZ) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.zNegative.momentumX, zNegative.H2.momentumX) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.zPositive.momentumX, zPositive.H2.momentumX) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(H2.xNegative.momentumZ, xNegative.H2.momentumZ) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.xPositive.momentumZ, xPositive.H2.momentumZ) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.yNegative.momentumX, yNegative.H2.momentumX) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.yPositive.momentumX, yPositive.H2.momentumX) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.zNegative.momentumY, zNegative.H2.momentumY) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.zPositive.momentumY, zPositive.H2.momentumY) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(H2.xNegative.entropy, xNegative.H2.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.xPositive.entropy, xPositive.H2.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.yNegative.entropy, yNegative.H2.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.yPositive.entropy, yPositive.H2.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.zNegative.entropy, zNegative.H2.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.zPositive.entropy, zPositive.H2.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // H2O
      // ---
      // Exchange
      connect(H2O.chemical, chemical.H2O) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(H2O.common, common) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2O.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(H2O.xNegative.material, xNegative.H2O.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.material, xPositive.H2O.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.material, yNegative.H2O.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.material, yPositive.H2O.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.material, zNegative.H2O.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.material, zPositive.H2O.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(H2O.xNegative.momentumY, xNegative.H2O.momentumY) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.momentumY, xPositive.H2O.momentumY) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.momentumZ, yNegative.H2O.momentumZ) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.momentumZ, yPositive.H2O.momentumZ) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.momentumX, zNegative.H2O.momentumX) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.momentumX, zPositive.H2O.momentumX) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(H2O.xNegative.momentumZ, xNegative.H2O.momentumZ) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.momentumZ, xPositive.H2O.momentumZ) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.momentumX, yNegative.H2O.momentumX) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.momentumX, yPositive.H2O.momentumX) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.momentumY, zNegative.H2O.momentumY) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.momentumY, zPositive.H2O.momentumY) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(H2O.xNegative.entropy, xNegative.H2O.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.entropy, xPositive.H2O.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.entropy, yNegative.H2O.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.entropy, yPositive.H2O.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.entropy, zNegative.H2O.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.entropy, zPositive.H2O.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // 'H+'
      // ----
      // Exchange
      connect('H+'.chemical, chemical.'H+') annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect('H+'.common, common) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect('H+'.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect('H+'.xNegative.material, xNegative.'H+'.material) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.xPositive.material, xPositive.'H+'.material) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.yNegative.material, yNegative.'H+'.material) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.yPositive.material, yPositive.'H+'.material) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.zNegative.material, zNegative.'H+'.material) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.zPositive.material, zPositive.'H+'.material) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect('H+'.xNegative.momentumY, xNegative.'H+'.momentumY) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.xPositive.momentumY, xPositive.'H+'.momentumY) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.yNegative.momentumZ, yNegative.'H+'.momentumZ) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.yPositive.momentumZ, yPositive.'H+'.momentumZ) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.zNegative.momentumX, zNegative.'H+'.momentumX) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.zPositive.momentumX, zPositive.'H+'.momentumX) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect('H+'.xNegative.momentumZ, xNegative.'H+'.momentumZ) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.xPositive.momentumZ, xPositive.'H+'.momentumZ) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.yNegative.momentumX, yNegative.'H+'.momentumX) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.yPositive.momentumX, yPositive.'H+'.momentumX) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.zNegative.momentumY, zNegative.'H+'.momentumY) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.zPositive.momentumY, zPositive.'H+'.momentumY) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect('H+'.xNegative.entropy, xNegative.'H+'.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.xPositive.entropy, xPositive.'H+'.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.yNegative.entropy, yNegative.'H+'.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.yPositive.entropy, yPositive.'H+'.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.zNegative.entropy, zNegative.'H+'.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.zPositive.entropy, zPositive.'H+'.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // N2
      // --
      // Exchange
      connect(N2.chemical, chemical.N2) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(N2.common, common) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(N2.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(N2.xNegative.material, xNegative.N2.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.xPositive.material, xPositive.N2.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.yNegative.material, yNegative.N2.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.yPositive.material, yPositive.N2.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.zNegative.material, zNegative.N2.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.zPositive.material, zPositive.N2.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(N2.xNegative.momentumY, xNegative.N2.momentumY) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.xPositive.momentumY, xPositive.N2.momentumY) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.yNegative.momentumZ, yNegative.N2.momentumZ) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.yPositive.momentumZ, yPositive.N2.momentumZ) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.zNegative.momentumX, zNegative.N2.momentumX) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.zPositive.momentumX, zPositive.N2.momentumX) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(N2.xNegative.momentumZ, xNegative.N2.momentumZ) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.xPositive.momentumZ, xPositive.N2.momentumZ) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.yNegative.momentumX, yNegative.N2.momentumX) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.yPositive.momentumX, yPositive.N2.momentumX) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.zNegative.momentumY, zNegative.N2.momentumY) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.zPositive.momentumY, zPositive.N2.momentumY) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(N2.xNegative.entropy, xNegative.N2.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.xPositive.entropy, xPositive.N2.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.yNegative.entropy, yNegative.N2.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.yPositive.entropy, yPositive.N2.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.zNegative.entropy, zNegative.N2.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.zPositive.entropy, zPositive.N2.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // O2
      // --
      // Exchange
      connect(O2.chemical, chemical.O2) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(O2.common, common) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(O2.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{8,-8},{8,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(O2.xNegative.material, xNegative.O2.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.xPositive.material, xPositive.O2.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.yNegative.material, yNegative.O2.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.yPositive.material, yPositive.O2.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.zNegative.material, zNegative.O2.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.zPositive.material, zPositive.O2.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(O2.xNegative.momentumY, xNegative.O2.momentumY) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.xPositive.momentumY, xPositive.O2.momentumY) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.yNegative.momentumZ, yNegative.O2.momentumZ) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.yPositive.momentumZ, yPositive.O2.momentumZ) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.zNegative.momentumX, zNegative.O2.momentumX) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.zPositive.momentumX, zPositive.O2.momentumX) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(O2.xNegative.momentumZ, xNegative.O2.momentumZ) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.xPositive.momentumZ, xPositive.O2.momentumZ) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.yNegative.momentumX, yNegative.O2.momentumX) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.yPositive.momentumX, yPositive.O2.momentumX) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.zNegative.momentumY, zNegative.O2.momentumY) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.zPositive.momentumY, zPositive.O2.momentumY) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(O2.xNegative.entropy, xNegative.O2.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.xPositive.entropy, xPositive.O2.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.yNegative.entropy, yNegative.O2.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.yPositive.entropy, yPositive.O2.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.zNegative.entropy, zNegative.O2.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.zPositive.entropy, zPositive.O2.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      annotation (
        defaultComponentPrefixes="replaceable",
        Documentation(info="<html>**If C or C19HF37O5S is included, then initXVel, initYVel, and initZVel are
    ignored.  Velocity is not initialized because it is set by C or C19HF37O5S.

    **If both C and C19HF37O5S are included, then reduce must be set to false (to avoid singularity, since
    both of those species set velocity).

    <p>Notes:<ul>
    <li>Only one constant-volume species may be enabled at once.  A failure will occur if
    <code>inclC</code> and <code>inclC19HF37O5S</code> are both <code>true</code>.</li>
    <li>The <code>inclReact</code> parameter may be set to
    <code>false</code>
    to eliminate unnecessary equations.</li>
    </ul></p>

    <p>For more information, see the
 <a href=\"modelica://FCSys.Subregions.Phases.BaseClasses.NullPhase\">NullPhase</a> model.</p></html>"),

        Diagram(graphics),
        Icon(graphics));
    end Phase;

    model Gas "Phase to represent gas"

      extends BaseClasses.NullPhase(final n_spec=countTrue({inclH2,inclH2O,
            inclN2,inclO2}));

      parameter Boolean inclReact=true "Include reaction(s), as appropriate"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(tab="Assumptions"));

      // Conditionally include species.
      parameter Boolean inclH2=false "<html>Hydrogen (H<sub>2</sub>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.H2.Gas.Fixed H2 if inclH2 constrainedby
        Species.Species(
        final k=k,
        final inclLin=inclLin,
        initMethX=if initXVel and reduceFinal then InitMethLinMom.None else
            InitMethLinMom.Velocity,
        initMethY=if initYVel and reduceFinal then InitMethLinMom.None else
            InitMethLinMom.Velocity,
        initMethZ=if initZVel and reduceFinal then InitMethLinMom.None else
            InitMethLinMom.Velocity,
        initMethTemp=if initTemp and reduceFinal then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceFinal then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceFinal then StateSelect.default else StateSelect.prefer))
        "<html>H<sub>2</sub> model</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclH2),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.H2O.Gas.Fixed H2O if inclH2O constrainedby
        Species.Species(
        final k=k,
        final inclLin=inclLin,
        initMethX=if initXVel and reduceFinal then InitMethLinMom.None else
            InitMethLinMom.Velocity,
        initMethY=if initYVel and reduceFinal then InitMethLinMom.None else
            InitMethLinMom.Velocity,
        initMethZ=if initZVel and reduceFinal then InitMethLinMom.None else
            InitMethLinMom.Velocity,
        initMethTemp=if initTemp and reduceFinal then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceFinal then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceFinal then StateSelect.default else StateSelect.prefer))
        "<html>H<sub>2</sub>O model</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclH2O),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      parameter Boolean inclN2=false "<html>Nitrogen (N<sub>2</sub>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.N2.Gas.Fixed N2 if inclN2 constrainedby
        Species.Species(
        final k=k,
        final inclLin=inclLin,
        initMethX=if initXVel and reduceFinal then InitMethLinMom.None else
            InitMethLinMom.Velocity,
        initMethY=if initYVel and reduceFinal then InitMethLinMom.None else
            InitMethLinMom.Velocity,
        initMethZ=if initZVel and reduceFinal then InitMethLinMom.None else
            InitMethLinMom.Velocity,
        initMethTemp=if initTemp and reduceFinal then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceFinal then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceFinal then StateSelect.default else StateSelect.prefer))
        "<html>N<sub>2</sub> model</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclN2),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      parameter Boolean inclO2=false "<html>Oxygen (O<sub>2</sub>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.O2.Gas.Fixed O2 if inclO2 constrainedby
        Species.Species(
        final k=k,
        final inclLin=inclLin,
        initMethX=if initXVel and reduceFinal then InitMethLinMom.None else
            InitMethLinMom.Velocity,
        initMethY=if initYVel and reduceFinal then InitMethLinMom.None else
            InitMethLinMom.Velocity,
        initMethZ=if initZVel and reduceFinal then InitMethLinMom.None else
            InitMethLinMom.Velocity,
        initMethTemp=if initTemp and reduceFinal then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceFinal then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceFinal then StateSelect.default else StateSelect.prefer))
        "<html>O<sub>2</sub> model</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclO2),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      Reaction '2H2+O2=2H2O'(final n_lin=n_lin, n_spec=3) if inclReact and (
        inclH2 and inclH2O and inclO2)
        annotation (Placement(transformation(extent={{-50,30},{-30,50}})));

    equation
      // Chemical interactions
      connect('2H2+O2=2H2O'.chemical[1], H2.chemical) annotation (Line(
          points={{-40,39.3333},{-7,7}},
          color={208,104,0},
          smooth=Smooth.None));
      connect('2H2+O2=2H2O'.chemical[2], O2.chemical) annotation (Line(
          points={{-40,40},{-7,7}},
          color={208,104,0},
          smooth=Smooth.None));
      connect('2H2+O2=2H2O'.chemical[3], H2O.chemical) annotation (Line(
          points={{-40,40.6667},{-7,7}},
          color={208,104,0},
          smooth=Smooth.None));

      // H2
      // --
      // Exchange
      connect(H2.chemical, chemical.H2) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(H2.common, common) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(H2.xNegative.material, xNegative.H2.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.xPositive.material, xPositive.H2.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.yNegative.material, yNegative.H2.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.yPositive.material, yPositive.H2.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.zNegative.material, zNegative.H2.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.zPositive.material, zPositive.H2.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(H2.xNegative.momentumY, xNegative.H2.momentumY) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.xPositive.momentumY, xPositive.H2.momentumY) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.yNegative.momentumZ, yNegative.H2.momentumZ) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.yPositive.momentumZ, yPositive.H2.momentumZ) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.zNegative.momentumX, zNegative.H2.momentumX) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.zPositive.momentumX, zPositive.H2.momentumX) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(H2.xNegative.momentumZ, xNegative.H2.momentumZ) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.xPositive.momentumZ, xPositive.H2.momentumZ) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.yNegative.momentumX, yNegative.H2.momentumX) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.yPositive.momentumX, yPositive.H2.momentumX) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.zNegative.momentumY, zNegative.H2.momentumY) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.zPositive.momentumY, zPositive.H2.momentumY) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(H2.xNegative.entropy, xNegative.H2.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.xPositive.entropy, xPositive.H2.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.yNegative.entropy, yNegative.H2.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.yPositive.entropy, yPositive.H2.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.zNegative.entropy, zNegative.H2.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.zPositive.entropy, zPositive.H2.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // H2O
      // ---
      // Exchange
      connect(H2O.chemical, chemical.H2O) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(H2O.common, common) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2O.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(H2O.xNegative.material, xNegative.H2O.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.material, xPositive.H2O.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.material, yNegative.H2O.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.material, yPositive.H2O.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.material, zNegative.H2O.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.material, zPositive.H2O.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(H2O.xNegative.momentumY, xNegative.H2O.momentumY) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.momentumY, xPositive.H2O.momentumY) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.momentumZ, yNegative.H2O.momentumZ) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.momentumZ, yPositive.H2O.momentumZ) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.momentumX, zNegative.H2O.momentumX) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.momentumX, zPositive.H2O.momentumX) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(H2O.xNegative.momentumZ, xNegative.H2O.momentumZ) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.momentumZ, xPositive.H2O.momentumZ) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.momentumX, yNegative.H2O.momentumX) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.momentumX, yPositive.H2O.momentumX) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.momentumY, zNegative.H2O.momentumY) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.momentumY, zPositive.H2O.momentumY) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(H2O.xNegative.entropy, xNegative.H2O.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.entropy, xPositive.H2O.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.entropy, yNegative.H2O.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.entropy, yPositive.H2O.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.entropy, zNegative.H2O.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.entropy, zPositive.H2O.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // N2
      // --
      // Exchange
      connect(N2.chemical, chemical.N2) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(N2.common, common) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(N2.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(N2.xNegative.material, xNegative.N2.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.xPositive.material, xPositive.N2.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.yNegative.material, yNegative.N2.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.yPositive.material, yPositive.N2.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.zNegative.material, zNegative.N2.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.zPositive.material, zPositive.N2.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(N2.xNegative.momentumY, xNegative.N2.momentumY) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.xPositive.momentumY, xPositive.N2.momentumY) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.yNegative.momentumZ, yNegative.N2.momentumZ) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.yPositive.momentumZ, yPositive.N2.momentumZ) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.zNegative.momentumX, zNegative.N2.momentumX) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.zPositive.momentumX, zPositive.N2.momentumX) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(N2.xNegative.momentumZ, xNegative.N2.momentumZ) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.xPositive.momentumZ, xPositive.N2.momentumZ) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.yNegative.momentumX, yNegative.N2.momentumX) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.yPositive.momentumX, yPositive.N2.momentumX) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.zNegative.momentumY, zNegative.N2.momentumY) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.zPositive.momentumY, zPositive.N2.momentumY) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(N2.xNegative.entropy, xNegative.N2.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.xPositive.entropy, xPositive.N2.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.yNegative.entropy, yNegative.N2.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.yPositive.entropy, yPositive.N2.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.zNegative.entropy, zNegative.N2.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.zPositive.entropy, zPositive.N2.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // O2
      // --
      // Exchange
      connect(O2.chemical, chemical.O2) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(O2.common, common) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(O2.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(O2.xNegative.material, xNegative.O2.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.xPositive.material, xPositive.O2.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.yNegative.material, yNegative.O2.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.yPositive.material, yPositive.O2.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.zNegative.material, zNegative.O2.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.zPositive.material, zPositive.O2.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(O2.xNegative.momentumY, xNegative.O2.momentumY) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.xPositive.momentumY, xPositive.O2.momentumY) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.yNegative.momentumZ, yNegative.O2.momentumZ) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.yPositive.momentumZ, yPositive.O2.momentumZ) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.zNegative.momentumX, zNegative.O2.momentumX) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.zPositive.momentumX, zPositive.O2.momentumX) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(O2.xNegative.momentumZ, xNegative.O2.momentumZ) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.xPositive.momentumZ, xPositive.O2.momentumZ) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.yNegative.momentumX, yNegative.O2.momentumX) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.yPositive.momentumX, yPositive.O2.momentumX) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.zNegative.momentumY, zNegative.O2.momentumY) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.zPositive.momentumY, zPositive.O2.momentumY) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(O2.xNegative.entropy, xNegative.O2.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.xPositive.entropy, xPositive.O2.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.yNegative.entropy, yNegative.O2.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.yPositive.entropy, yPositive.O2.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.zNegative.entropy, zNegative.O2.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.zPositive.entropy, zPositive.O2.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      annotation (
        defaultComponentPrefixes="replaceable",
        Documentation(info="<html><p>Notes:<ul><li>The <code>inclReact</code> parameter may be set to
    <code>false</code>
    to eliminate unnecessary equations.</li></ul>
</p>
<p>For more information, see the
 <a href=\"modelica://FCSys.Subregions.Phases.BaseClasses.NullPhase\">NullPhase</a> model.</p></html>"),

        Diagram(graphics));
    end Gas;

    model Graphite "Phase to represent graphite"
      extends BaseClasses.NullPhase(final n_spec=countTrue({inclC,'incle-'}),
          common(phi(fixed=if inclC then fill(false, n_lin) else {initXVel,
                initYVel,initZVel}[cartAxes])));

      // Conditionally include species.
      parameter Boolean inclC=false "Carbon (C)" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.C.Graphite.Fixed C if inclC constrainedby
        Species.Species(
        final k=k,
        final inclLin=inclLin,
        initMethTemp=if initTemp and reduceFinal then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceFinal then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceFinal then StateSelect.default else StateSelect.prefer))
        "C model" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclC),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      parameter Boolean 'incle-'=false "<html>Electrons (e<sup>-</sup>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.'e-'.Graphite.Fixed 'e-' if 'incle-' constrainedby
        Species.Species(
        final k=k,
        final inclLin=inclLin,
        initMethX=if (initXVel or inclC) and reduceFinal then InitMethLinMom.None
             else InitMethLinMom.Velocity,
        initMethY=if (initYVel or inclC) and reduceFinal then InitMethLinMom.None
             else InitMethLinMom.Velocity,
        initMethZ=if (initZVel or inclC) and reduceFinal then InitMethLinMom.None
             else InitMethLinMom.Velocity,
        initMethTemp=if initTemp and reduceFinal then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceFinal then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceFinal then StateSelect.default else StateSelect.prefer))
        "<html>e<sup>-</sup> model</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable='incle-'),
        Placement(transformation(extent={{-10,-10},{10,10}})));

    equation
      // C
      // -
      // Exchange
      connect(C.chemical, chemical.C) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(C.common, common) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(C.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(C.xNegative.material, xNegative.C.material) annotation (Line(
          points={{-10,6.10623e-16},{-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.xPositive.material, xPositive.C.material) annotation (Line(
          points={{10,6.10623e-16},{40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.yNegative.material, yNegative.C.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.yPositive.material, yPositive.C.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,26},{5.55112e-16,26},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.zNegative.material, zNegative.C.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.zPositive.material, zPositive.C.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(C.xNegative.momentumY, xNegative.C.momentumY) annotation (Line(
          points={{-10,6.10623e-16},{-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.xPositive.momentumY, xPositive.C.momentumY) annotation (Line(
          points={{10,6.10623e-16},{40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.yNegative.momentumZ, yNegative.C.momentumZ) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.yPositive.momentumZ, yPositive.C.momentumZ) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,26},{5.55112e-16,26},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.zNegative.momentumX, zNegative.C.momentumX) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.zPositive.momentumX, zPositive.C.momentumX) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(C.xNegative.momentumZ, xNegative.C.momentumZ) annotation (Line(
          points={{-10,6.10623e-16},{-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.xPositive.momentumZ, xPositive.C.momentumZ) annotation (Line(
          points={{10,6.10623e-16},{40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.yNegative.momentumX, yNegative.C.momentumX) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.yPositive.momentumX, yPositive.C.momentumX) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,26},{5.55112e-16,26},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.zNegative.momentumY, zNegative.C.momentumY) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.zPositive.momentumY, zPositive.C.momentumY) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(C.xNegative.entropy, xNegative.C.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.xPositive.entropy, xPositive.C.entropy) annotation (Line(
          points={{10,6.10623e-16},{40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.yNegative.entropy, yNegative.C.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.yPositive.entropy, yPositive.C.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,26},{5.55112e-16,26},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.zNegative.entropy, zNegative.C.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.zPositive.entropy, zPositive.C.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // e-
      // --
      // Exchange
      connect('e-'.chemical, chemical.'e-') annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect('e-'.common, common) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect('e-'.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect('e-'.xNegative.material, xNegative.'e-'.material) annotation (
          Line(
          points={{-10,6.10623e-16},{-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.xPositive.material, xPositive.'e-'.material) annotation (
          Line(
          points={{10,6.10623e-16},{40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.yNegative.material, yNegative.'e-'.material) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.yPositive.material, yPositive.'e-'.material) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,26},{5.55112e-16,26},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.zNegative.material, zNegative.'e-'.material) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.zPositive.material, zPositive.'e-'.material) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect('e-'.xNegative.momentumY, xNegative.'e-'.momentumY) annotation (
          Line(
          points={{-10,6.10623e-16},{-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.xPositive.momentumY, xPositive.'e-'.momentumY) annotation (
          Line(
          points={{10,6.10623e-16},{40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.yNegative.momentumZ, yNegative.'e-'.momentumZ) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.yPositive.momentumZ, yPositive.'e-'.momentumZ) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,26},{5.55112e-16,26},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.zNegative.momentumX, zNegative.'e-'.momentumX) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.zPositive.momentumX, zPositive.'e-'.momentumX) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect('e-'.xNegative.momentumZ, xNegative.'e-'.momentumZ) annotation (
          Line(
          points={{-10,6.10623e-16},{-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.xPositive.momentumZ, xPositive.'e-'.momentumZ) annotation (
          Line(
          points={{10,6.10623e-16},{40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.yNegative.momentumX, yNegative.'e-'.momentumX) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.yPositive.momentumX, yPositive.'e-'.momentumX) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,26},{5.55112e-16,26},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.zNegative.momentumY, zNegative.'e-'.momentumY) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.zPositive.momentumY, zPositive.'e-'.momentumY) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect('e-'.xNegative.entropy, xNegative.'e-'.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.xPositive.entropy, xPositive.'e-'.entropy) annotation (Line(
          points={{10,6.10623e-16},{40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.yNegative.entropy, yNegative.'e-'.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.yPositive.entropy, yPositive.'e-'.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,26},{5.55112e-16,26},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.zNegative.entropy, zNegative.'e-'.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.zPositive.entropy, zPositive.'e-'.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      annotation (
        defaultComponentPrefixes="replaceable",
        Documentation(info="<html>**If C is included, then initXVel, initYVel, and initZVel are
    ignored.  Velocity is not initialized because it is set by C.

    <p>For more information, see the
 <a href=\"modelica://FCSys.Subregions.Phases.BaseClasses.NullPhase\">NullPhase</a> model.</p></html>"),

        Diagram(graphics));
    end Graphite;

    model Ionomer "Phase to represent ionomer"
      extends BaseClasses.NullPhase(final n_spec=countTrue({inclC19HF37O5S,
            inclH2O,'inclH+'}), common(phi(fixed=if inclC19HF37O5S then fill(
                false, n_lin) else {initXVel,initYVel,initZVel}[cartAxes])));

      // Conditionally include species.
      parameter Boolean inclC19HF37O5S=false
        "<html>Nafion sulfonate (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.C19HF37O5S.Solid.Fixed C19HF37O5S if inclC19HF37O5S
        constrainedby Species.Species(
        final k=k,
        final inclLin=inclLin,
        initMethTemp=if initTemp and reduceFinal then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceFinal then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceFinal then StateSelect.default else StateSelect.prefer))
        "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S model</html>"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclC19HF37O5S),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.H2O.Gas.Fixed H2O if inclH2O constrainedby
        Species.Species(
        final k=k,
        final inclLin=inclLin,
        initMethX=if (initXVel or inclC19HF37O5S) and reduceFinal then
            InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethY=if (initYVel or inclC19HF37O5S) and reduceFinal then
            InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethZ=if (initZVel or inclC19HF37O5S) and reduceFinal then
            InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethTemp=if initTemp then InitMethScalar.None else InitMethScalar.Temperature,

        phi(each stateSelect=if reduceFinal then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceFinal then StateSelect.default else StateSelect.prefer))
        "<html>H<sub>2</sub>O model</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclH2O),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      parameter Boolean 'inclH+'=false "<html>Protons (H<sup>+</sup>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.'H+'.Solid.Fixed 'H+' if 'inclH+' constrainedby
        Species.Species(
        final k=k,
        final inclLin=inclLin,
        initMethX=if (initXVel or inclC19HF37O5S) and reduceFinal then
            InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethY=if (initYVel or inclC19HF37O5S) and reduceFinal then
            InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethZ=if (initZVel or inclC19HF37O5S) and reduceFinal then
            InitMethLinMom.None else InitMethLinMom.Velocity,
        initMethTemp=if initTemp and reduceFinal then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceFinal then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceFinal then StateSelect.default else StateSelect.prefer))
        "<html>H<sup>+</sup> model</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable='inclH+'),
        Placement(transformation(extent={{-10,-10},{10,10}})));

    equation
      // C19HF37O5S
      // ----------
      // Exchange
      connect(C19HF37O5S.chemical, chemical.C19HF37O5S) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(C19HF37O5S.common, common) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(C19HF37O5S.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(C19HF37O5S.xNegative.material, xNegative.C19HF37O5S.material)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.xPositive.material, xPositive.C19HF37O5S.material)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.yNegative.material, yNegative.C19HF37O5S.material)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.yPositive.material, yPositive.C19HF37O5S.material)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.zNegative.material, zNegative.C19HF37O5S.material)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.zPositive.material, zPositive.C19HF37O5S.material)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(C19HF37O5S.xNegative.momentumY, xNegative.C19HF37O5S.momentumY)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.xPositive.momentumY, xPositive.C19HF37O5S.momentumY)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.yNegative.momentumZ, yNegative.C19HF37O5S.momentumZ)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.yPositive.momentumZ, yPositive.C19HF37O5S.momentumZ)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.zNegative.momentumX, zNegative.C19HF37O5S.momentumX)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.zPositive.momentumX, zPositive.C19HF37O5S.momentumX)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(C19HF37O5S.xNegative.momentumZ, xNegative.C19HF37O5S.momentumZ)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.xPositive.momentumZ, xPositive.C19HF37O5S.momentumZ)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.yNegative.momentumX, yNegative.C19HF37O5S.momentumX)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.yPositive.momentumX, yPositive.C19HF37O5S.momentumX)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.zNegative.momentumY, zNegative.C19HF37O5S.momentumY)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.zPositive.momentumY, zPositive.C19HF37O5S.momentumY)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(C19HF37O5S.xNegative.entropy, xNegative.C19HF37O5S.entropy)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.xPositive.entropy, xPositive.C19HF37O5S.entropy)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.yNegative.entropy, yNegative.C19HF37O5S.entropy)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.yPositive.entropy, yPositive.C19HF37O5S.entropy)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.zNegative.entropy, zNegative.C19HF37O5S.entropy)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.zPositive.entropy, zPositive.C19HF37O5S.entropy)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // H+
      // --
      // Exchange
      connect('H+'.chemical, chemical.'H+') annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect('H+'.common, common) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect('H+'.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect('H+'.xNegative.material, xNegative.'H+'.material) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.xPositive.material, xPositive.'H+'.material) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.yNegative.material, yNegative.'H+'.material) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.yPositive.material, yPositive.'H+'.material) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.zNegative.material, zNegative.'H+'.material) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.zPositive.material, zPositive.'H+'.material) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect('H+'.xNegative.momentumY, xNegative.'H+'.momentumY) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.xPositive.momentumY, xPositive.'H+'.momentumY) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.yNegative.momentumZ, yNegative.'H+'.momentumZ) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.yPositive.momentumZ, yPositive.'H+'.momentumZ) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.zNegative.momentumX, zNegative.'H+'.momentumX) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.zPositive.momentumX, zPositive.'H+'.momentumX) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect('H+'.xNegative.momentumZ, xNegative.'H+'.momentumZ) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.xPositive.momentumZ, xPositive.'H+'.momentumZ) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.yNegative.momentumX, yNegative.'H+'.momentumX) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.yPositive.momentumX, yPositive.'H+'.momentumX) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.zNegative.momentumY, zNegative.'H+'.momentumY) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.zPositive.momentumY, zPositive.'H+'.momentumY) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect('H+'.xNegative.entropy, xNegative.'H+'.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.xPositive.entropy, xPositive.'H+'.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.yNegative.entropy, yNegative.'H+'.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.yPositive.entropy, yPositive.'H+'.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.zNegative.entropy, zNegative.'H+'.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.zPositive.entropy, zPositive.'H+'.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // H2O
      // ---
      // Exchange
      connect(H2O.chemical, chemical.H2O) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(H2O.common, common) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));

      connect(H2O.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(H2O.xNegative.material, xNegative.H2O.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.material, xPositive.H2O.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.material, yNegative.H2O.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.material, yPositive.H2O.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.material, zNegative.H2O.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.material, zPositive.H2O.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(H2O.xNegative.momentumY, xNegative.H2O.momentumY) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.momentumY, xPositive.H2O.momentumY) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.momentumZ, yNegative.H2O.momentumZ) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.momentumZ, yPositive.H2O.momentumZ) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.momentumX, zNegative.H2O.momentumX) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.momentumX, zPositive.H2O.momentumX) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(H2O.xNegative.momentumZ, xNegative.H2O.momentumZ) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.momentumZ, xPositive.H2O.momentumZ) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.momentumX, yNegative.H2O.momentumX) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.momentumX, yPositive.H2O.momentumX) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.momentumY, zNegative.H2O.momentumY) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.momentumY, zPositive.H2O.momentumY) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(H2O.xNegative.entropy, xNegative.H2O.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.entropy, xPositive.H2O.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.entropy, yNegative.H2O.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.entropy, yPositive.H2O.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.entropy, zNegative.H2O.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.entropy, zPositive.H2O.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      annotation (
        defaultComponentPrefixes="replaceable",
        Documentation(info="<html>**If C19HF37O5S is included, then initXVel, initYVel, and initZVel are
    ignored.  Velocity is not initialized because it is set by C19HF37O5S.

    <p>For more information, see the
 <a href=\"modelica://FCSys.Subregions.Phases.BaseClasses.NullPhase\">NullPhase</a> model.</p></html>"),

        Diagram(graphics));
    end Ionomer;

    model Liquid "Phase to represent liquid"

      extends BaseClasses.NullPhase(final n_spec=if inclH2O then 1 else 0);

      // Conditionally include species.
      parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.H2O.Gas.Fixed H2O if inclH2O constrainedby
        Species.Species(
        final k=k,
        final inclLin=inclLin,
        initMethX=if initXVel and reduceFinal then InitMethLinMom.None else
            InitMethLinMom.Velocity,
        initMethY=if initYVel and reduceFinal then InitMethLinMom.None else
            InitMethLinMom.Velocity,
        initMethZ=if initZVel and reduceFinal then InitMethLinMom.None else
            InitMethLinMom.Velocity,
        initMethTemp=if initTemp and reduceFinal then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceFinal then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceFinal then StateSelect.default else StateSelect.prefer))
        "<html>H<sub>2</sub>O model</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclH2O),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      // TODO:  Add and use a model for H2O liquid.
    equation
      // H2O
      // ---
      // Exchange
      connect(H2O.chemical, chemical.H2O) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(H2O.common, common) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2O.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(H2O.xNegative.material, xNegative.H2O.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.material, xPositive.H2O.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.material, yNegative.H2O.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.material, yPositive.H2O.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.material, zNegative.H2O.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.material, zPositive.H2O.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(H2O.xNegative.momentumY, xNegative.H2O.momentumY) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.momentumY, xPositive.H2O.momentumY) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.momentumZ, yNegative.H2O.momentumZ) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.momentumZ, yPositive.H2O.momentumZ) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.momentumX, zNegative.H2O.momentumX) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.momentumX, zPositive.H2O.momentumX) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(H2O.xNegative.momentumZ, xNegative.H2O.momentumZ) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.momentumZ, xPositive.H2O.momentumZ) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.momentumX, yNegative.H2O.momentumX) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.momentumX, yPositive.H2O.momentumX) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.momentumY, zNegative.H2O.momentumY) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.momentumY, zPositive.H2O.momentumY) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(H2O.xNegative.entropy, xNegative.H2O.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.entropy, xPositive.H2O.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.entropy, yNegative.H2O.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.entropy, yPositive.H2O.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.entropy, zNegative.H2O.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.entropy, zPositive.H2O.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      annotation (
        defaultComponentPrefixes="replaceable",
        Documentation(info="<html><p>For more information, see the
 <a href=\"modelica://FCSys.Subregions.Phases.BaseClasses.NullPhase\">NullPhase</a> model.</p></html>"),

        Diagram(graphics),
        Icon(graphics));
    end Liquid;

    package BaseClasses "Base classes (not for direct use)"
      extends Modelica.Icons.BasesPackage;
      model NullPhase "Model for a phase with no species or reactions"
        //extends FCSys.BaseClasses.Icons.Names.Middle;

        // Geometric parameters
        parameter Q.NumberAbsolute k[3](
          each min=Modelica.Constants.small,
          each final nominal=1) = {1,1,1}
          "<html>Adjustment factor for transport (<b>k</b>)</html>"
          annotation (Dialog(group="Geometry"));
        outer parameter Q.Length L[Axis](each final min=Modelica.Constants.small)
          "Length" annotation (HideResult=true,missingInnerMessage=
              "This model should be used within the Subregion model.");
        outer parameter Q.Area A[Axis] "Cross-sectional area" annotation (
            HideResult=true, missingInnerMessage=
              "This model should be used within the Subregion model.");

        // Assumptions
        parameter Boolean reduce=true
          "Same velocity and temperature for all species" annotation (
          HideResult=true,
          Dialog(tab="Assumptions", enable=n_spec > 1),
          choices(__Dymola_checkBox=true));

        // Initialization
        parameter Boolean initXVel=true "Initialize the x component"
          annotation (
          compact=true,
          Dialog(
            tab="Initialization",
            group="Velocity",
            enable=reduceFinal),
          __Dymola_joinNext=true,
          choices(__Dymola_checkBox=true));
        parameter Boolean initYVel=true "Initialize the y component"
          annotation (
          compact=true,
          Dialog(
            tab="Initialization",
            group="Velocity",
            enable=reduceFinal),
          __Dymola_joinNext=true,
          choices(__Dymola_checkBox=true));
        parameter Boolean initZVel=true "Initialize the z component"
          annotation (
          compact=true,
          Dialog(
            tab="Initialization",
            group="Velocity",
            enable=reduceFinal),
          __Dymola_joinNext=true,
          choices(__Dymola_checkBox=true));
        parameter Q.Velocity phi_IC[Axis]={0,0,0}
          "<html>Initial velocity (<b>&phi;</b><sub>IC</sub>)</html>"
          annotation (Dialog(tab="Initialization", group="Velocity"));
        // This is always enabled in the dialog since it is used as a guess value.
        parameter Boolean initTemp=true "Initialize" annotation (
          compact=true,
          Dialog(
            tab="Initialization",
            group="Temperature",
            enable=reduceFinal),
          choices(__Dymola_checkBox=true));

        parameter Q.TemperatureAbsolute T_IC(nominal=298.15*U.K, start=defaults.T)
          "<html>Initial temperature (<i>T</i><sub>IC</sub>)</html>"
          annotation (Dialog(tab="Initialization", group="Temperature"));
        // This is always enabled in the dialog since it is used as a guess value.

        parameter Boolean inclLin[3]={true,false,false}
          "true, if each component of linear momentum is included"
          annotation (Evaluate=true,Dialog(tab="Assumptions"));

        FCSys.Connectors.InertAmagat inert(final n_lin=n_lin) if n_spec > 0
          annotation (Placement(transformation(extent={{10,-30},{30,-10}}),
              iconTransformation(extent={{70,-90},{90,-70}})));
        FCSys.Connectors.FaceBus xNegative if n_spec > 0
          "Negative face along the x axis" annotation (Placement(transformation(
                extent={{-50,-10},{-30,10}}), iconTransformation(extent={{-90,-10},
                  {-70,10}})));
        FCSys.Connectors.FaceBus xPositive if n_spec > 0
          "Positive face along the x axis" annotation (Placement(transformation(
                extent={{30,-10},{50,10}}), iconTransformation(extent={{70,-10},
                  {90,10}})));
        FCSys.Connectors.FaceBus yNegative if n_spec > 0
          "Negative face along the y axis" annotation (Placement(transformation(
                extent={{-10,-50},{10,-30}}), iconTransformation(extent={{-10,-94},
                  {10,-74}})));
        FCSys.Connectors.FaceBus yPositive if n_spec > 0
          "Positive face along the y axis" annotation (Placement(transformation(
                extent={{-10,30},{10,50}}), iconTransformation(extent={{-10,90},
                  {10,110}})));
        FCSys.Connectors.FaceBus zNegative if n_spec > 0
          "Negative face along the z axis" annotation (Placement(transformation(
                extent={{10,10},{30,30}}), iconTransformation(extent={{40,40},{
                  60,60}})));
        FCSys.Connectors.FaceBus zPositive if n_spec > 0
          "Positive face along the z axis" annotation (Placement(transformation(
                extent={{-30,-30},{-10,-10}}), iconTransformation(extent={{-90,
                  -90},{-70,-70}})));

        PhaseBoundary phaseBoundary(final n_lin=n_lin) if n_spec > 0
          "Phase boundary" annotation (Placement(transformation(
              extent={{-18,-18},{18,18}},
              rotation=0,
              origin={0,0})));
        // This component is conditional because if two or more empty phases
        // (without any species included) were connected within a subregion, there
        // would be a mathematical singularity.

        FCSys.Connectors.ChemicalBus chemical if n_spec > 0 annotation (
            Placement(transformation(extent={{-30,10},{-10,30}}),
              iconTransformation(extent={{-60,40},{-40,60}})));

      protected
        parameter Integer n_spec "Number of species";
        final parameter Boolean reduceFinal=reduce and n_spec > 1
          "true, if same velocity and temperature for all (more than 1) species ";
        final parameter Integer n_lin=countTrue(inclLin)
          "Number of components of linear momentum"
          annotation (Evaluate=true, HideResult=true);
        final parameter Integer cartAxes[n_lin]=index(inclLin)
          "Cartesian-axis indices of the axes of linear momentum";

        Connectors.InertInternal common(
          n_lin=n_lin,
          phi(
            each stateSelect=StateSelect.prefer,
            final start=phi_IC[cartAxes],
            fixed={initXVel,initYVel,initZVel}[cartAxes]),
          T(
            stateSelect=StateSelect.prefer,
            final start=T_IC,
            fixed=initTemp)) if reduceFinal
          "Internal connector to directly couple velocity and temperature"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={0,0}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={0,0})));

        outer BCs.Defaults defaults "Default settings" annotation (Placement(
              transformation(extent={{40,40},{60,60}}), iconTransformation(
                extent={{-10,90},{10,110}})));

      equation
        // Inert interactions
        connect(phaseBoundary.inertA, inert) annotation (Line(
            points={{12,-12},{20,-20}},
            color={72,90,180},
            smooth=Smooth.None));

        annotation (
          Documentation(info="<html>**if one of the species has setTemp=true, then initTemp should be set to false.
    The temperature should not be not initialized here since it there would be too many initialization equations.  same for initXVel, initYVel, initZVel.

    <p>Assumptions:
           <ol>
           <li>The factors that may cause anisotropic transport (<b><i>k</i></b>)
          are common to all of the species within the phase.</li>
</ol></p>

<p>Notes:<ul>
  <li>The x-axis component of linear momentum is included by default.  At least one component must be included.</li></ul></html>"),

          Icon(graphics={Ellipse(
                      extent={{-40,100},{40,20}},
                      lineColor={127,127,127},
                      startAngle=30,
                      endAngle=149,
                      pattern=LinePattern.Dash,
                      fillPattern=FillPattern.Solid,
                      fillColor={225,225,225}),Ellipse(
                      extent={{20,-4},{100,-84}},
                      lineColor={127,127,127},
                      startAngle=270,
                      endAngle=390,
                      pattern=LinePattern.Dash,
                      fillPattern=FillPattern.Solid,
                      fillColor={225,225,225}),Ellipse(
                      extent={{-100,-4},{-20,-84}},
                      lineColor={127,127,127},
                      startAngle=149,
                      endAngle=270,
                      pattern=LinePattern.Dash,
                      fillPattern=FillPattern.Solid,
                      fillColor={225,225,225}),Polygon(
                      points={{60,-84},{-60,-84},{-94.5,-24},{-34.5,80},{34.5,
                  80},{94.5,-24},{60,-84}},
                      pattern=LinePattern.None,
                      fillPattern=FillPattern.Sphere,
                      smooth=Smooth.None,
                      fillColor={225,225,225},
                      lineColor={0,0,0}),Line(
                      points={{-60,-84},{60,-84}},
                      color={127,127,127},
                      pattern=LinePattern.Dash,
                      smooth=Smooth.None),Line(
                      points={{34.5,80},{94.5,-24}},
                      color={127,127,127},
                      pattern=LinePattern.Dash,
                      smooth=Smooth.None),Line(
                      points={{-34.5,80},{-94.5,-24}},
                      color={127,127,127},
                      pattern=LinePattern.Dash,
                      smooth=Smooth.None),Text(
                      extent={{-100,-20},{100,20}},
                      textString="%name",
                      lineColor={0,0,0})}),
          Diagram(graphics));
      end NullPhase;
    end BaseClasses;
  end Phases;

  package Species
    "Models for single-species storage, transport, and exchange of material, volume, and linear momentum"
    extends Modelica.Icons.Package;
    package C "C"
      extends Modelica.Icons.Package;
      package Graphite "C graphite"
        extends Modelica.Icons.Package;
        model Calibrated "Ideal correlations, with adjustment factors"
          extends SpeciesInertStagnant(
            redeclare FCSys.Characteristics.C.Graphite Data,
            redeclare Q.Resistivity beta_Phi=k_beta_Phi*Data.gamma(T),
            redeclare Q.ResistivityThermal beta_S=k_beta_S*Data.gamma(T),
            redeclare Q.ResistivityThermal gamma_S=k_gamma_S*Data.gamma(T));

          parameter Q.NumberAbsolute k_beta_Phi(final nominal=1) = 1
            "<html>Adjustment factor for exchange of linear momentum (<i>k</i><sub>&beta;&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal exchange (<i>k</i><sub>&beta;<i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_gamma_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal transport (<i>k</i>&gamma;<sub><i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C",
            Documentation(info="<html><p>For more information, see the
    <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model.</p></html>"),

            Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}), graphics),
            Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}), graphics={Text(
                          extent={{-150,90},{-118,52}},
                          lineColor={0,0,255},
                          textString="%t.test")}));

        end Calibrated;

        model Correlated "Ideal correlations from kinetic theory"
          extends SpeciesInertStagnant(
            redeclare FCSys.Characteristics.C.Graphite Data,
            redeclare Q.Resistivity beta_Phi=Data.gamma(T),
            redeclare Q.ResistivityThermal beta_S=Data.gamma(T),
            redeclare Q.ResistivityThermal gamma_S=Data.gamma(T));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C",
            Documentation(info=
                  "<html><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model.</p></html>"),

            Diagram(graphics),
            Icon(graphics));

        end Correlated;

        model Fixed "Fixed properties"
          extends SpeciesInertStagnant(
            redeclare FCSys.Characteristics.C.Graphite Data(
              Deltah0_f=0,
              Deltah0=0,
              specHeatCapPow=0,
              T_lim_c0={0,Modelica.Constants.inf},
              b_c0=[935*U.J*Data.m/(U.kg*U.K)],
              B_c0=[-300*U.K*935*U.J*Data.m/(U.kg*U.K) + Data.Deltah0_f, 0]),
            redeclare parameter Q.Resistivity beta_Phi=Data.gamma(),
            redeclare parameter Q.ResistivityThermal beta_S=Data.gamma(),
            redeclare parameter Q.ResistivityThermal gamma_S=U.m*U.K/(11.1*U.W));
          // See the documentation for a table of values.
          // Note:  Parameter expressions (e.g., involving defaults.T) are not used
          // here since they would render the parameters unadjustable in Dymola 7.4.
          // A similar note applies to the other species.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Fixed specific heat capacity (independent of temperature)
    <li>Transport and exchange properties (&beta;<sub><i>S</i></sub>, &gamma;<sub><i>S</i></sub>, etc.) are fixed (e.g., independent of temperature)</li>
    </ol></p>

    <p>The default specific heat capacity at constant pressure (<code>b_c0=[0, 935*U.J*Data.m/(U.kg*U.K)]</code>) and thermal
   resistivity (<code>gamma_S=U.m*U.K/(11.1*U.W)</code>) is based on data of graphite fiber epoxy (25% vol)<br>composite at 300 K from
   Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 909].
   Related data is listed in Table 1.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"bottom\"><b>Table 1:</b> Properties of forms of C [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 909].</caption>
    <tr>
      <th rowspan=3 valign=\"middle\"><code>T/U.K</code></th>
      <th rowspan=1 colspan=2 width=1 valign=\"middle\">Diamond (type IIa)</th>
      <th rowspan=1 colspan=1 width=1 valign=\"middle\">Amorphous<br>carbon</th>
      <th rowspan=1 colspan=3 width=1 valign=\"middle\">Graphite (pyrolytic)</th>
      <th rowspan=1 colspan=3 width=1>Graphite fiber epoxy (25% vol)<br>composite</th>
    </tr>
    <tr>
      <th rowspan=2 valign=\"middle\"><code>c0*U.kg<br>*U.K<br>/(U.J<br>*Data.m)</code></th>
      <th rowspan=2 valign=\"middle\"><code>gamma_S<br>*U.K<br>/(U.m<br>*U.K)</code></th>
      <th rowspan=2 valign=\"middle\"><code>gamma_S<br>*U.W<br>/(U.m<br>*U.K)</code></th>
      <th rowspan=2 valign=\"middle\"><code>c0*U.kg<br>*U.W<br>/(U.J<br>*Data.m)</code></th>
      <th rowspan=1 colspan=2><code>gamma_S*U.W/(U.m*U.K)</code></th>
      <th rowspan=2 valign=\"middle\"><code>c0*U.kg<br>*U.K<br>/(U.J<br>*Data.m)</code></th>
      <th rowspan=1 colspan=2><code>gamma_S*U.W/(U.m*U.K)</code></th>
    </tr>
    <tr>
      <th valign=\"middle\">Parallel<br>to layers</th>
      <th valign=\"middle\">Perpendicular<br>to layers</th>
      <th valign=\"middle\">Parallel<br>to layers</th>
      <th valign=\"middle\">Perpendicular<br>to layers</th>
    </tr>
<tr><td>100</td><td>21</td><td>1/10000</td><td>1/0.67</td><td>136</td><td>1/4970</td><td>1/16.8</td><td>337</td><td>1/5.7</td><td>1/0.46</td></tr>
<tr><td>200</td><td>194</td><td>1/4000</td><td>1/1.18</td><td>411</td><td>1/3230</td><td>1/9.23</td><td>642</td><td>1/8.7</td><td>1/0.68</td></tr>
<tr><td>300</td><td>509</td><td>1/2300</td><td>1/1.89</td><td>709</td><td>1/1950</td><td>1/5.70</td><td>935</td><td>1/11.1</td><td>1/0.87</td></tr>
<tr><td>400</td><td>853</td><td>1/1540</td><td>1/2.19</td><td>992</td><td>1/1390</td><td>1/4.09</td><td>1216</td><td>1/13.0</td><td>1/1.1</td></tr>
<tr><td>600</td><td>-</td><td>-</td><td>1/2.37</td><td>1406</td><td>1/892</td><td>1/2.68</td><td>-</td><td>-</td><td>-</td></tr>
<tr><td>800</td><td>-</td><td>-</td><td>1/2.53</td><td>1650</td><td>1/667</td><td>1/2.01</td><td>-</td><td>-</td><td>-</td></tr>
<tr><td>1000</td><td>-</td><td>-</td><td>1/2.84</td><td>1793</td><td>1/534</td><td>1/1.60</td><td>-</td><td>-</td><td>-</td></tr>
<tr><td>1200</td><td>-</td><td>-</td><td>1/3.48</td><td>1890</td><td>1/448</td><td>1/1.34</td><td>-</td><td>-</td><td>-</td></tr>
<tr><td>1500</td><td>-</td><td>-</td><td>-</td><td>1974</td><td>1/357</td><td>1/1.08</td><td>-</td><td>-</td><td>-</td></tr>
<tr><td>2000</td><td>-</td><td>-</td><td>-</td><td>2043</td><td>1/262</td><td>1/0.81</td><td>-</td><td>-</td><td>-</td></tr>
  </table>

  <p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model.</p></html>"));
        end Fixed;
      end Graphite;
    end C;

    package C19HF37O5S
      "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S</html>"
      extends Modelica.Icons.Package;
      package Solid
        "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S solid</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Ideal correlations, with adjustment factors"
          extends SpeciesInertStagnant(
            redeclare FCSys.Characteristics.C19HF37O5S.Solid Data,
            redeclare Q.Resistivity beta_Phi=k_beta_Phi*Data.gamma(T),
            redeclare Q.ResistivityThermal gamma_S=k_gamma_S*Data.gamma(T));

          parameter Q.NumberAbsolute k_beta_Phi(final nominal=1) = 1
            "<html>Adjustment factor for exchange of linear momentum (<i>k</i><sub>&beta;&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal exchange (<i>k</i><sub>&beta;<i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_gamma_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal transport (<i>k</i>&gamma;<sub><i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C19HF37O5S",
            Documentation(info=
                  "<html><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model.</p></html>"),

            Icon(graphics),
            Diagram(graphics));

        end Calibrated;

        model Correlated "Ideal correlations from kinetic theory"
          extends SpeciesInertStagnant(
            redeclare FCSys.Characteristics.C19HF37O5S.Solid Data,
            redeclare Q.Resistivity beta_Phi=Data.gamma(T),
            redeclare Q.ResistivityThermal beta_S=Data.gamma(T),
            redeclare Q.ResistivityThermal gamma_S=Data.gamma(T));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C19HF37O5S",
            Documentation(info=
                  "<html><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends SpeciesInertStagnant(
            redeclare FCSys.Characteristics.C19HF37O5S.Solid Data,
            redeclare parameter Q.Resistivity beta_Phi=Data.gamma(),
            redeclare parameter Q.ResistivityThermal beta_S=Data.gamma(),
            redeclare parameter Q.ResistivityThermal gamma_S=U.m*U.K/(0.16*U.W));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C19HF37O5S",
            Documentation(info="<html><p>Assumptions:
    <ol>
    <li>Transport and exchange properties (&beta;<sub><i>S</i></sub>, &gamma;<sub><i>S</i></sub>, etc.) are fixed (e.g., independent of temperature)</li>
    </ol></p>

    <p>Notes:
    <ul><li>The default thermal transport resistivity (<code>gamma_S=U.m*U.K/(0.16*U.W)</code>) is of dry
  Nafion 115 [<a href=\"modelica://FCSys.UsersGuide.References\">Kandlikar2009</a>, p. 1277.</li>
  </ul>
  </p><p>For more information, see the
    <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model.</p></html>"));
        end Fixed;
      end Solid;
    end C19HF37O5S;

    package 'e-' "<html>e<sup>-</sup></html>"
      extends Modelica.Icons.Package;
      package Graphite "<html>e<sup>-</sup> in graphite</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Ideal correlations, with adjustment factors"
          extends Species0Amount(
            redeclare FCSys.Characteristics.'e-'.Graphite Data,
            initMethPartNum=InitMethScalar.Amount,
            redeclare Q.Resistivity beta_Phi=k_beta_Phi*Data.gamma(T),
            redeclare Q.Resistivity gamma_N=k_gamma_N*Data.gamma(T),
            redeclare Q.Resistivity gamma_Phi=k_gamma_Phi*Data.gamma(T),
            redeclare Q.ResistivityThermal gamma_S=k_gamma_S*Data.gamma(T));

          // Note:  initMethPartNum may not be Pressure (which is default) since the
          // EOS doesn't involve pressure.

          parameter Q.NumberAbsolute k_beta_Phi(final nominal=1) = 1
            "<html>Adjustment factor for exchange of linear momentum (<i>k</i><sub>&beta;&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_gamma_N(final nominal=1) = 1
            "<html>Adjustment factor for material transport (<i>k</i>&gamma;<sub><i>N</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_gamma_Phi(final nominal=1) = 1
            "<html>Adjustment factor for transport of linear momentum (<i>k</i>&gamma;<sub>&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_gamma_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal transport (<i>k</i>&gamma;<sub><i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'e-'",
            Documentation(info=
                  "<html><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model.</p></html>"));
        end Calibrated;

        model Correlated "Ideal correlations from kinetic theory"
          extends Species0Amount(
            redeclare FCSys.Characteristics.'e-'.Graphite Data,
            initMethPartNum=InitMethScalar.Amount,
            redeclare Q.Resistivity beta_Phi=Data.gamma(T),
            redeclare Q.Resistivity gamma_N=Data.gamma(T),
            redeclare Q.Resistivity gamma_Phi=Data.gamma(T),
            redeclare Q.ResistivityThermal gamma_S=Data.gamma(T));
          // Note:  initMethPartNum may not be Pressure (which is default) since the
          // EOS doesn't involve pressure.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'e-'",
            Documentation(info=
                  "<html><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends Species0Amount(
            redeclare FCSys.Characteristics.'e-'.Graphite Data,
            initMethPartNum=InitMethScalar.Amount,
            redeclare parameter Q.Resistivity beta_Phi=Data.gamma(),
            redeclare parameter Q.ResistivityThermal beta_S=Data.gamma(),
            redeclare parameter Q.Resistivity gamma_N=Data.gamma(),
            redeclare parameter Q.Resistivity gamma_Phi=Data.gamma(),
            redeclare parameter Q.ResistivityThermal gamma_S=Data.gamma());

          // Note:  initMethPartNum may not be Pressure (which is default) since the
          // EOS doesn't involve pressure.

          annotation (
            group="Material properties",
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'e-'",
            Documentation(info="<html>
    <p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model.</p></html>"),

            Diagram(graphics));

        end Fixed;
      end Graphite;
    end 'e-';

    package H2 "<html>H<sub>2</sub></html>"
      extends Modelica.Icons.Package;
      package Gas "<html>H<sub>2</sub> gas</html>"
        extends Modelica.Icons.Package;
        model Calibrated
          "Ideal or empirical correlations, w/ adjustment factors"
          extends Species(
            redeclare FCSys.Characteristics.H2.Gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare Q.Resistivity beta_Phi=k_beta_Phi*Data.gamma(T),
            redeclare Q.ResistivityThermal beta_S=k_beta_S*Data.gamma(T),
            redeclare Q.Resistivity gamma_N=k_gamma_N*Data.gamma(T),
            redeclare Q.Resistivity gamma_Phi=k_gamma_Phi*(if
                empiricalTransverse then Data.gamma_Phi(T) else Data.gamma(T)),

            redeclare Q.ResistivityThermal gamma_S=k_gamma_S*(if
                empiricalThermal then Data.gamma_S(T) else Data.gamma(T)));

          parameter Q.NumberAbsolute k_beta_Phi(final nominal=1) = 1
            "<html>Adjustment factor for exchange of linear momentum (<i>k</i><sub>&beta;&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal exchange (<i>k</i><sub>&beta;<i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_gamma_N(final nominal=1) = 1
            "<html>Adjustment factor for material transport (<i>k</i>&gamma;<sub><i>N</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_gamma_Phi(final nominal=1) = 1
            "<html>Adjustment factor for transport of linear momentum (<i>k</i>&gamma;<sub>&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_gamma_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal transport (<i>k</i>&gamma;<sub><i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Boolean empiricalTransverse=true
            "Use empirical correlation for transverse resistivity (otherwise, theoretical)"
            annotation (
            Evaluate=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Correlations",
              compact=true));
          parameter Boolean empiricalThermal=true
            "Use empirical correlation for thermal resistivity (otherwise, theoretical)"
            annotation (
            Evaluate=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Correlations",
              compact=true));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"));
        end Calibrated;

        model Correlated "Ideal (kinetic) or empirical correlations"
          extends Species(
            redeclare FCSys.Characteristics.H2.Gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare Q.Resistivity beta_Phi=Data.gamma(T),
            redeclare Q.ResistivityThermal beta_S=Data.gamma(T),
            redeclare Q.Resistivity gamma_N=Data.gamma(T),
            redeclare Q.Resistivity gamma_Phi=if empiricalTransverse then
                Data.gamma_Phi(T) else Data.gamma(T),
            redeclare Q.ResistivityThermal gamma_S=if empiricalThermal then
                Data.gamma_S(T) else Data.gamma(T));

          parameter Boolean empiricalTransverse=true
            "Use empirical correlation for transverse resistivity (otherwise, theoretical)"
            annotation (
            Evaluate=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Correlations",
              compact=true));
          parameter Boolean empiricalThermal=true
            "Use empirical correlation for thermal resistivity (otherwise, theoretical)"
            annotation (
            Evaluate=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Correlations",
              compact=true));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare FCSys.Characteristics.H2.Gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare parameter Q.Resistivity beta_Phi=Data.gamma(),
            redeclare parameter Q.ResistivityThermal beta_S=U.m*U.K/(183e-3*U.W),

            redeclare parameter Q.Resistivity gamma_N=Data.gamma(),
            redeclare parameter Q.Resistivity gamma_Phi=Data.m/(89.6e-7*U.Pa*U.s),

            redeclare parameter Q.ResistivityThermal gamma_S=U.m*U.K/(183e-3*U.W)
                /100);

          // See the documentation for a table of values.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
    <li>Transport and exchange properties (&beta;<sub>&Phi;</sub>, &gamma;<sub>&Phi;</sub>, etc.) are fixed (e.g., independent of temperature)</li>
</ol></p>

<p>Additional notes:<ul>
<li>
  The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

<p>The default transport resistivities (<code>gamma_Phi=Data.m/(89.6e-7*U.Pa*U.s)</code>
and <code>gamma_S=U.m*U.K/(183e-3*U.W)</code>) are based on data of H<sub>2</sub> gas at 1 atm and
  300 K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 919&ndash;920].
  Table 1 lists the properties at  other temperatures. </p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"bottom\"><b>Table 1:</b> Properties of H<sub>2</sub> gas at 1 atm (<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 919&ndash;920).</caption>
  <tr>
      <th rowspan=2 valign=\"middle\"><code>T/U.K</code></th>
      <th rowspan=2 width=1><code>c*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th rowspan=2 width=1 ><code>gamma_Phi<br>*U.Pa*U.s/Data.m</code></th>
      <th rowspan=2 width=1 ><code>gamma_S*U.W<br>/(U.m*U.K)</code></th>
    </tr>
    <tr>
<tr><td>100</td><td>11.23e3</td><td>1/42.1e-7</td><td>1/67.0e-3</td></tr>
<tr><td>150</td><td>12.60e3</td><td>1/56.0e-7</td><td>1/101e-3</td></tr>
<tr><td>200</td><td>13.54e3</td><td>1/68.1e-7</td><td>1/131e-3</td></tr>
<tr><td>250</td><td>14.06e3</td><td>1/78.9e-7</td><td>1/157e-3</td></tr>
<tr><td>300</td><td>14.31e3</td><td>1/89.6e-7</td><td>1/183e-3</td></tr>
<tr><td>350</td><td>14.43e3</td><td>1/98.8e-7</td><td>1/204e-3</td></tr>
<tr><td>400</td><td>14.48e3</td><td>1/108.2e-7</td><td>1/226e-3</td></tr>
<tr><td>450</td><td>14.50e3</td><td>1/117.2e-7</td><td>1/247e-3</td></tr>
<tr><td>500</td><td>14.52e3</td><td>1/126.4e-7</td><td>1/266e-3</td></tr>
<tr><td>550</td><td>14.53e3</td><td>1/134.3e-7</td><td>1/285e-3</td></tr>
<tr><td>600</td><td>14.55e3</td><td>1/142.4e-7</td><td>1/305e-3</td></tr>
<tr><td>700</td><td>14.61e3</td><td>1/157.8e-7</td><td>1/342e-3</td></tr>
<tr><td>800</td><td>14.70e3</td><td>1/172.4e-7</td><td>1/378e-3</td></tr>
<tr><td>900</td><td>14.83e3</td><td>1/186.5e-7</td><td>1/412e-3</td></tr>
<tr><td>1000</td><td>14.99e3</td><td>1/201.3e-7</td><td>1/448e-3</td></tr>
<tr><td>1100</td><td>15.17e3</td><td>1/213.0e-7</td><td>1/488e-3</td></tr>
<tr><td>1200</td><td>15.37e3</td><td>1/226.2e-7</td><td>1/528e-3</td></tr>
<tr><td>1300</td><td>15.59e3</td><td>1/238.5e-7</td><td>1/568e-3</td></tr>
<tr><td>1400</td><td>15.81e3</td><td>1/250.7e-7</td><td>1/610e-3</td></tr>
<tr><td>1500</td><td>16.02e3</td><td>1/262.7e-7</td><td>1/655e-3</td></tr>
<tr><td>1600</td><td>16.28e3</td><td>1/273.7e-7</td><td>1/697e-3</td></tr>
<tr><td>1700</td><td>16.58e3</td><td>1/284.9e-7</td><td>1/742e-3</td></tr>
<tr><td>1800</td><td>16.96e3</td><td>1/296.1e-7</td><td>1/786e-3</td></tr>
<tr><td>1900</td><td>17.49e3</td><td>1/307.2e-7</td><td>1/835e-3</td></tr>
<tr><td>2000</td><td>18.25e3</td><td>1/318.2e-7</td><td>1/878e-3</td></tr>
    </tr>
  </table>
<p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"),

            Diagram(graphics));

        end Fixed;
      end Gas;
    end H2;

    package H2O "<html>H<sub>2</sub>O</html>"
      extends Modelica.Icons.Package;
      package Gas "<html>H<sub>2</sub>O gas</html>"
        extends Modelica.Icons.Package;
        model Calibrated
          "Ideal or empirical correlations, w/ adjustment factors"
          extends Species(
            redeclare FCSys.Characteristics.H2O.Gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare Q.Resistivity beta_Phi=k_beta_Phi*Data.gamma(T),
            redeclare Q.ResistivityThermal beta_S=k_beta_S*Data.gamma(T),
            redeclare Q.Resistivity gamma_N=k_gamma_N*Data.gamma(T),
            redeclare Q.Resistivity gamma_Phi=k_gamma_Phi*(if
                empiricalTransverse then Data.gamma_Phi(T) else Data.gamma(T)),

            redeclare Q.ResistivityThermal gamma_S=k_gamma_S*(if
                empiricalThermal then Data.gamma_S(T) else Data.gamma(T)));

          parameter Q.NumberAbsolute k_beta_Phi(final nominal=1) = 1
            "<html>Adjustment factor for exchange of linear momentum (<i>k</i><sub>&beta;&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal exchange (<i>k</i><sub>&beta;<i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_gamma_N(final nominal=1) = 1
            "<html>Adjustment factor for material transport (<i>k</i>&gamma;<sub><i>N</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_gamma_Phi(final nominal=1) = 1
            "<html>Adjustment factor for transport of linear momentum (<i>k</i>&gamma;<sub>&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_gamma_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal transport (<i>k</i>&gamma;<sub><i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Boolean empiricalTransverse=true
            "Use empirical correlation for transverse resistivity (otherwise, theoretical)"
            annotation (
            Evaluate=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Correlations",
              compact=true));
          parameter Boolean empiricalThermal=true
            "Use empirical correlation for thermal resistivity (otherwise, theoretical)"
            annotation (
            Evaluate=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Correlations",
              compact=true));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"));
        end Calibrated;

        model Correlated "Ideal (kinetic) or empirical correlations"
          extends Species(
            redeclare FCSys.Characteristics.H2O.Gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare Q.Resistivity beta_Phi=Data.gamma(T),
            redeclare Q.ResistivityThermal beta_S=Data.gamma(T),
            redeclare Q.Resistivity gamma_N=Data.gamma(T),
            redeclare Q.Resistivity gamma_Phi=if empiricalTransverse then
                Data.gamma_Phi(T) else Data.gamma(T),
            redeclare Q.ResistivityThermal gamma_S=if empiricalThermal then
                Data.gamma_S(T) else Data.gamma(T));

          parameter Boolean empiricalTransverse=true
            "Use empirical correlation for transverse resistivity (otherwise, theoretical)"
            annotation (
            Evaluate=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Correlations",
              compact=true));
          parameter Boolean empiricalThermal=true
            "Use empirical correlation for thermal resistivity (otherwise, theoretical)"
            annotation (
            Evaluate=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Correlations",
              compact=true));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare FCSys.Characteristics.H2O.Gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare parameter Q.Resistivity beta_Phi=Data.gamma(),
            redeclare parameter Q.ResistivityThermal beta_S=U.m*U.K/(19.6e-3*U.W),

            redeclare parameter Q.Resistivity gamma_N=Data.gamma(),
            redeclare parameter Q.Resistivity gamma_Phi=Data.m/(9.09e-6*U.Pa*U.s),

            redeclare parameter Q.ResistivityThermal gamma_S=U.m*U.K/(19.6e-3*U.W));

          // See the documentation for tables of values.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
    <li>Transport and exchange properties (&beta;<sub>&Phi;</sub>, &gamma;<sub>&Phi;</sub>, etc.) are fixed (e.g., independent of temperature)</li>
  </ol></p>

          <p>Notes:<ul>
<ul><li>
  The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

<p>The default transport resistivities (<code>gamma_Phi=Data.m/(9.09e-6*U.Pa*U.s)</code>
and <code>gamma_S=U.m*U.K/(19.6e-3*U.W)</code>) are of H<sub>2</sub>O gas at saturation pressure and
  300 K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 921].  Table 1 lists the properties at
  saturation pressure and other temperatures.  Table 2 lists the properties of H<sub>2</sub>O gas at 1 atm.
  Table 3 lists resistivity to transport of linear momentum (inverse of bulk viscosity) based on its ratio to
  transverse resistivity (inverse of shear viscosity). See also
  <a href=\"http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html\">http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html</a>.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"bottom\"><b>Table 1:</b> Properties of H<sub>2</sub>O at saturation pressure (<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 924&ndash;925).</caption>
  <tr>
      <th rowspan=2 valign=\"middle\"><code>T/U.K</code></th>
      <th rowspan=1 colspan=2 width=1><code>c*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th rowspan=1 colspan=2 width=1><code>gamma_Phi<br>*U.Pa*U.s/Data.m</code></th>
      <th rowspan=1 colspan=2 width=1><code>gamma_S*U.W<br>/(U.m*U.K)</code></th>
    </tr>
    <tr>
      <th>Gas</th>
      <th>Liquid</th>
      <th>Gas</th>
      <th>Liquid</th>
      <th>Gas</th>
      <th>Liquid</th>
    </tr>
<tr><td>273.15</td><td>1.854e3</td><td>4.217e3</td><td>1/8.02e-6</td><td>1/1750e-6</td><td>1/18.2e-3</td><td>1/569e-3</td></tr>
<tr><td>275</td><td>1.855e3</td><td>4.211e3</td><td>1/8.09e-6</td><td>1/1652e-6</td><td>1/18.3e-3</td><td>1/574e-3</td></tr>
<tr><td>280</td><td>1.858e3</td><td>4.198e3</td><td>1/8.29e-6</td><td>1/1422e-6</td><td>1/18.6e-3</td><td>1/582e-3</td></tr>
<tr><td>285</td><td>1.861e3</td><td>4.189e3</td><td>1/8.49e-6</td><td>1/1225e-6</td><td>1/18.9e-3</td><td>1/590e-3</td></tr>
<tr><td>290</td><td>1.864e3</td><td>4.184e3</td><td>1/8.69e-6</td><td>1/1080e-6</td><td>1/19.3e-3</td><td>1/598e-3</td></tr>
<tr><td>295</td><td>1.868e3</td><td>4.181e3</td><td>1/8.89e-6</td><td>1/959e-6</td><td>1/19.5e-3</td><td>1/606e-3</td></tr>
<tr><td>300</td><td>1.872e3</td><td>4.179e3</td><td>1/9.09e-6</td><td>1/855e-6</td><td>1/19.6e-3</td><td>1/613e-3</td></tr>
<tr><td>305</td><td>1.877e3</td><td>4.178e3</td><td>1/9.29e-6</td><td>1/769e-6</td><td>1/20.1e-3</td><td>1/620e-3</td></tr>
<tr><td>310</td><td>1.882e3</td><td>4.178e3</td><td>1/9.49e-6</td><td>1/695e-6</td><td>1/20.4e-3</td><td>1/628e-3</td></tr>
<tr><td>315</td><td>1.888e3</td><td>4.179e3</td><td>1/9.69e-6</td><td>1/631e-6</td><td>1/20.7e-3</td><td>1/634e-3</td></tr>
<tr><td>320</td><td>1.895e3</td><td>4.180e3</td><td>1/9.89e-6</td><td>1/577e-6</td><td>1/21.0e-3</td><td>1/640e-3</td></tr>
<tr><td>325</td><td>1.903e3</td><td>4.182e3</td><td>1/10.09e-6</td><td>1/528e-6</td><td>1/21.3e-3</td><td>1/645e-3</td></tr>
<tr><td>330</td><td>1.911e3</td><td>4.184e3</td><td>1/10.29e-6</td><td>1/489e-6</td><td>1/21.7e-3</td><td>1/650e-3</td></tr>
<tr><td>335</td><td>1.920e3</td><td>4.186e3</td><td>1/10.49e-6</td><td>1/453e-6</td><td>1/22.0e-3</td><td>1/656e-3</td></tr>
<tr><td>340</td><td>1.930e3</td><td>4.188e3</td><td>1/10.69e-6</td><td>1/420e-6</td><td>1/22.3e-3</td><td>1/660e-3</td></tr>
<tr><td>345</td><td>1.941e3</td><td>4.191e3</td><td>1/10.89e-6</td><td>1/389e-6</td><td>1/22.6e-3</td><td>1/668e-3</td></tr>
<tr><td>350</td><td>1.954e3</td><td>4.195e3</td><td>1/11.09e-6</td><td>1/365e-6</td><td>1/23.0e-3</td><td>1/668e-3</td></tr>
<tr><td>355</td><td>1.968e3</td><td>4.199e3</td><td>1/11.29e-6</td><td>1/343e-6</td><td>1/23.3e-3</td><td>1/671e-3</td></tr>
<tr><td>360</td><td>1.983e3</td><td>4.203e3</td><td>1/11.49e-6</td><td>1/324e-6</td><td>1/23.7e-3</td><td>1/674e-3</td></tr>
<tr><td>365</td><td>1.999e3</td><td>4.209e3</td><td>1/11.69e-6</td><td>1/306e-6</td><td>1/24.1e-3</td><td>1/677e-3</td></tr>
<tr><td>370</td><td>2.017e3</td><td>4.214e3</td><td>1/11.89e-6</td><td>1/289e-6</td><td>1/24.5e-3</td><td>1/679e-3</td></tr>
<tr><td>373.15</td><td>2.029e3</td><td>4.217e3</td><td>1/12.02e-6</td><td>1/279e-6</td><td>1/24.8e-3</td><td>1/680e-3</td></tr>
<tr><td>375</td><td>2.036e3</td><td>4.220e3</td><td>1/12.09e-6</td><td>1/274e-6</td><td>1/24.9e-3</td><td>1/681e-3</td></tr>
<tr><td>380</td><td>2.057e3</td><td>4.226e3</td><td>1/12.29e-6</td><td>1/260e-6</td><td>1/25.4e-3</td><td>1/683e-3</td></tr>
<tr><td>385</td><td>2.080e3</td><td>4.232e3</td><td>1/12.49e-6</td><td>1/248e-6</td><td>1/25.8e-3</td><td>1/685e-3</td></tr>
<tr><td>390</td><td>2.104e3</td><td>4.239e3</td><td>1/12.69e-6</td><td>1/237e-6</td><td>1/26.3e-3</td><td>1/686e-3</td></tr>
<tr><td>400</td><td>2.158e3</td><td>4.256e3</td><td>1/13.05e-6</td><td>1/217e-6</td><td>1/27.2e-3</td><td>1/688e-3</td></tr>
<tr><td>410</td><td>2.221e3</td><td>4.278e3</td><td>1/13.42e-6</td><td>1/200e-6</td><td>1/28.2e-3</td><td>1/688e-3</td></tr>
<tr><td>420</td><td>2.291e3</td><td>4.302e3</td><td>1/13.79e-6</td><td>1/185e-6</td><td>1/29.8e-3</td><td>1/688e-3</td></tr>
<tr><td>430</td><td>2.369e3</td><td>4.331e3</td><td>1/14.14e-6</td><td>1/173e-6</td><td>1/30.4e-3</td><td>1/685e-3</td></tr>
<tr><td>440</td><td>2.46e3</td><td>4.36e3</td><td>1/14.50e-6</td><td>1/162e-6</td><td>1/3.17e-3</td><td>1/682e-3</td></tr>
<tr><td>450</td><td>2.56e3</td><td>4.40e3</td><td>1/14.85e-6</td><td>1/152e-6</td><td>1/33.1e-3</td><td>1/678e-3</td></tr>
<tr><td>460</td><td>2.68e3</td><td>4.44e3</td><td>1/15.19e-6</td><td>1/143e-6</td><td>1/34.6e-3</td><td>1/673e-3</td></tr>
<tr><td>470</td><td>2.79e3</td><td>4.48e3</td><td>1/15.54e-6</td><td>1/136e-6</td><td>1/36.3e-3</td><td>1/667e-3</td></tr>
<tr><td>480</td><td>2.94e3</td><td>4.53e3</td><td>1/15.88e-6</td><td>1/129e-6</td><td>1/38.1e-3</td><td>1/660e-3</td></tr>
<tr><td>490</td><td>3.10e3</td><td>4.59e3</td><td>1/16.23e-6</td><td>1/124e-6</td><td>1/40.1e-3</td><td>1/651e-3</td></tr>
<tr><td>500</td><td>3.27e3</td><td>4.66e3</td><td>1/16.59e-6</td><td>1/118e-6</td><td>1/42.3e-3</td><td>1/642e-3</td></tr>
<tr><td>510</td><td>3.47e3</td><td>4.74e3</td><td>1/16.95e-6</td><td>1/113e-6</td><td>1/44.7e-3</td><td>1/631e-3</td></tr>
<tr><td>520</td><td>3.70e3</td><td>4.84e3</td><td>1/17.33e-6</td><td>1/108e-6</td><td>1/47.5e-3</td><td>1/621e-3</td></tr>
<tr><td>530</td><td>3.96e3</td><td>4.95e3</td><td>1/17.72e-6</td><td>1/104e-6</td><td>1/50.6e-3</td><td>1/608e-3</td></tr>
<tr><td>540</td><td>4.27e3</td><td>5.08e3</td><td>1/18.1e-6</td><td>1/101e-6</td><td>1/54.0e-3</td><td>1/594e-3</td></tr>
<tr><td>550</td><td>4.64e3</td><td>5.24e3</td><td>1/18.6e-6</td><td>1/97e-6</td><td>1/58.3e-3</td><td>1/580e-3</td></tr>
<tr><td>560</td><td>5.09e3</td><td>5.43e3</td><td>1/19.1e-6</td><td>1/94e-6</td><td>1/63.7e-3</td><td>1/563e-3</td></tr>
<tr><td>570</td><td>5.67e3</td><td>5.68e3</td><td>1/19.7e-6</td><td>1/91e-6</td><td>1/76.7e-3</td><td>1/548e-3</td></tr>
<tr><td>580</td><td>6.40e3</td><td>6.00e3</td><td>1/20.4e-6</td><td>1/88e-6</td><td>1/76.7e-3</td><td>1/528e-3</td></tr>
<tr><td>590</td><td>7.35e3</td><td>6.41e3</td><td>1/21.5e-6</td><td>1/84e-6</td><td>1/84.1e-3</td><td>1/513e-3</td></tr>
<tr><td>600</td><td>8.75e3</td><td>7.00e3</td><td>1/22.7e-6</td><td>1/81e-6</td><td>1/92.9e-3</td><td>1/497e-3</td></tr>
<tr><td>610</td><td>11.1e3</td><td>7.85e3</td><td>1/24.1e-6</td><td>1/77e-6</td><td>1/103e-3</td><td>1/467e-3</td></tr>
<tr><td>620</td><td>15.4e3</td><td>9.35e3</td><td>1/25.9e-6</td><td>1/72e-6</td><td>1/114e-3</td><td>1/444e-3</td></tr>
<tr><td>635</td><td>18.3e3</td><td>10.6e3</td><td>1/27.0e-6</td><td>1/70e-6</td><td>1/121e-3</td><td>1/430e-3</td></tr>
<tr><td>630</td><td>22.1e3</td><td>12.6e3</td><td>1/28.0e-6</td><td>1/67e-6</td><td>1/130e-3</td><td>1/412e-3</td></tr>
<tr><td>635</td><td>27.6e3</td><td>16.4e3</td><td>1/30.0e-6</td><td>1/64e-6</td><td>1/141e-3</td><td>1/392e-3</td></tr>
<tr><td>640</td><td>42e3</td><td>26e3</td><td>1/32.0e-6</td><td>1/59e-6</td><td>1/155e-3</td><td>1/367e-3</td></tr>
    </tr>
  </table>

<br>

    <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"bottom\"><b>Table 2:</b> Properties of H<sub>2</sub>O gas at 1 atm (<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 921).</caption>
  <tr>
      <th valign=\"middle\"><code>T/U.K</code></th>
      <th width=1><code>c*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1 ><code>gamma_Phi<br>*U.Pa*U.s/Data.m</code></th>
      <th width=1 ><code>gamma_S*U.W<br>/(U.m*U.K)</code></th>
    </tr>
<tr><td>380</td><td>2.060e3</td><td>1/127.1e-7</td><td>1/24.6e-3</td></tr>
<tr><td>400</td><td>2.014e3</td><td>1/134.4e-7</td><td>1/26.1e-3</td></tr>
<tr><td>450</td><td>1.980e3</td><td>1/152.5e-7</td><td>1/29.9e-3</td></tr>
<tr><td>500</td><td>1.985e3</td><td>1/170.4e-7</td><td>1/33.9e-3</td></tr>
<tr><td>550</td><td>1.997e3</td><td>1/188.4e-7</td><td>1/37.9e-3</td></tr>
<tr><td>600</td><td>2.206e3</td><td>1/206.7e-7</td><td>1/42.2e-3</td></tr>
<tr><td>650</td><td>2.056e3</td><td>1/224.7e-7</td><td>1/46.4e-3</td></tr>
<tr><td>700</td><td>2.085e3</td><td>1/242.6e-7</td><td>1/50.5e-3</td></tr>
<tr><td>750</td><td>2.119e3</td><td>1/260.4e-7</td><td>1/54.9e-3</td></tr>
<tr><td>800</td><td>2.152e3</td><td>1/278.6e-7</td><td>1/59.2e-3</td></tr>
<tr><td>850</td><td>2.186e3</td><td>1/296.9e-7</td><td>1/63.7e-3</td></tr>
  </table></ul>
<br>

  </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"));
        end Fixed;
      end Gas;
    end H2O;

    package 'H+' "<html>H<sup>+</sup></html>"
      extends Modelica.Icons.Package;
      package Solid "<html>H<sup>+</sup> in solid</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Ideal correlations, with adjustment factors"
          extends Species0Amount(
            redeclare FCSys.Characteristics.'H+'.Solid Data,
            initMethPartNum=InitMethScalar.Amount,
            redeclare Q.Resistivity beta_Phi=k_beta_Phi*Data.gamma(T),
            redeclare Q.Resistivity gamma_N=k_gamma_N*Data.gamma(T),
            redeclare Q.Resistivity gamma_Phi=k_gamma_Phi*Data.gamma(T),
            redeclare Q.ResistivityThermal gamma_S=k_gamma_S*Data.gamma(T));

          // Note:  initMethPartNum may not be Pressure (which is default) since the
          // EOS doesn't involve pressure.

          parameter Q.NumberAbsolute k_beta_Phi(final nominal=1) = 1
            "<html>Adjustment factor for exchange of linear momentum (<i>k</i><sub>&beta;&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_gamma_N(final nominal=1) = 1
            "<html>Adjustment factor for material transport (<i>k</i>&gamma;<sub><i>N</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_gamma_Phi(final nominal=1) = 1
            "<html>Adjustment factor for transport of linear momentum (<i>k</i>&gamma;<sub>&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_gamma_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal transport (<i>k</i>&gamma;<sub><i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'H+'",
            Documentation(info=
                  "<html><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model.</p></html>"));
        end Calibrated;

        model Correlated "Ideal correlations from kinetic theory"
          extends Species0Amount(
            redeclare FCSys.Characteristics.'H+'.Solid Data,
            initMethPartNum=InitMethScalar.Amount,
            redeclare Q.Resistivity beta_Phi=Data.gamma(T),
            redeclare Q.Resistivity gamma_N=Data.gamma(T),
            redeclare Q.Resistivity gamma_Phi=Data.gamma(T),
            redeclare Q.ResistivityThermal gamma_S=Data.gamma(T));
          // Note:  initMethPartNum may not be Pressure (which is default) since the
          // EOS doesn't involve pressure.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'H+'",
            Documentation(info=
                  "<html><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends Species0Amount(
            redeclare FCSys.Characteristics.'H+'.Solid Data,
            initMethPartNum=InitMethScalar.Amount,
            redeclare parameter Q.Resistivity beta_Phi=Data.gamma(),
            redeclare parameter Q.ResistivityThermal beta_S=U.m*U.K/(0.1661*U.W),

            redeclare parameter Q.Resistivity gamma_N=Data.gamma(),
            redeclare parameter Q.Resistivity gamma_Phi=Data.m/(5.3e-6*U.Pa*U.s),

            redeclare parameter Q.ResistivityThermal gamma_S=U.m*U.K/(0.1661*U.W));

          // Note:  initMethPartNum may not be Pressure (which is default) since
          // overrideEOS is true.

          // See the documentation for a table of values.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'H+'",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Transport and exchange properties (&beta;<sub>&Phi;</sub>, &gamma;<sub>&Phi;</sub>, etc.) are fixed 
    (e.g., independent of temperature)</li>
    </ol></p>

  <p>The default transport resistivities (<code>gamma_Phi=Data.m/(5.3e-6*U.Pa*U.s)</code> and <code>gamma_S=U.m*U.K/(0.1661*U.W)</code>) are of H gas
  (rather than H<sup>+</sup>) at 300 K from [<a href=\"modelica://FCSys.UsersGuide.References\">Schetz1996</a>, p. 139].
  Table 1 lists the properties at other temperatures.</p>

    <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"bottom\"><b>Table 1:</b> Properties of H gas (not H<sup>+</sup>) (<a href=\"modelica://FCSys.UsersGuide.References\">Schetz and Fuhs, 1996</a>, p. 139)</caption>
<tr>
      <th valign=\"middle\"><code>T/U.K</code></th>
      <th width=1><code>beta_Phi<br>*U.Pa*U.s/Data.m</code></th>
      <th width=1><code>beta_S*U.W<br>/(U.m*U.K)</code></th>
    </tr>
<tr><td>200</td><td>1/3.8e-6</td><td>1/0.1197</td></tr>
<tr><td>300</td><td>1/5.3e-6</td><td>1/0.1661</td></tr>
<tr><td>400</td><td>1/6.7e-6</td><td>1/0.2094</td></tr>
<tr><td>500</td><td>1/8.1e-6</td><td>1/0.2507</td></tr>
<tr><td>600</td><td>1/9.3e-6</td><td>1/0.2906</td></tr>
<tr><td>700</td><td>1/10.6e-6</td><td>1/0.3292</td></tr>
<tr><td>800</td><td>1/11.8e-6</td><td>1/0.3670</td></tr>
<tr><td>900</td><td>1/13.0e-6</td><td>1/0.4040</td></tr>
<tr><td>1000</td><td>1/14.2e-6</td><td>1/0.4403</td></tr>
<tr><td>1100</td><td>1/15.3e-6</td><td>1/0.4761</td></tr>
<tr><td>1200</td><td>1/16.5e-6</td><td>1/0.5114</td></tr>
<tr><td>1300</td><td>1/17.6e-6</td><td>1/0.5462</td></tr>
<tr><td>1400</td><td>1/18.7e-6</td><td>1/0.5807</td></tr>
<tr><td>1500</td><td>1/19.8e-6</td><td>1/0.6149</td></tr>
<tr><td>1600</td><td>1/20.9e-6</td><td>1/0.6487</td></tr>
<tr><td>1700</td><td>1/22.0e-6</td><td>1/0.6823</td></tr>
<tr><td>1800</td><td>1/23.1e-6</td><td>1/0.7156</td></tr>
<tr><td>1900</td><td>1/24.2e-6</td><td>1/0.7488</td></tr>
<tr><td>2000</td><td>1/25.2e-6</td><td>1/0.7817</td></tr>
<tr><td>2100</td><td>1/26.3e-6</td><td>1/0.8144</td></tr>
<tr><td>2200</td><td>1/27.3e-6</td><td>1/0.8470</td></tr>
<tr><td>2300</td><td>1/28.4e-6</td><td>1/0.8794</td></tr>
<tr><td>2400</td><td>1/29.4e-6</td><td>1/0.9117</td></tr>
<tr><td>2500</td><td>1/30.5e-6</td><td>1/0.9438</td></tr>
<tr><td>2600</td><td>1/31.5e-6</td><td>1/0.9758</td></tr>
<tr><td>2700</td><td>1/32.5e-6</td><td>1/1.0077</td></tr>
<tr><td>2800</td><td>1/33.6e-6</td><td>1/1.0395</td></tr>
<tr><td>2900</td><td>1/34.6e-6</td><td>1/1.0711</td></tr>
<tr><td>3000</td><td>1/35.6e-6</td><td>1/1.1027</td></tr>
<tr><td>3100</td><td>1/36.6e-6</td><td>1/1.1347</td></tr>
<tr><td>3200</td><td>1/37.7e-6</td><td>1/1.1664</td></tr>
<tr><td>3300</td><td>1/38.7e-6</td><td>1/1.1978</td></tr>
<tr><td>3400</td><td>1/39.7e-6</td><td>1/1.2288</td></tr>
<tr><td>3500</td><td>1/40.7e-6</td><td>1/1.2592</td></tr>
<tr><td>3600</td><td>1/41.6e-6</td><td>1/1.2884</td></tr>
<tr><td>3700</td><td>1/42.5e-6</td><td>1/1.3171</td></tr>
<tr><td>3800</td><td>1/43.4e-6</td><td>1/1.3455</td></tr>
<tr><td>3900</td><td>1/44.4e-6</td><td>1/1.3735</td></tr>
<tr><td>4000</td><td>1/45.2e-6</td><td>1/1.4012</td></tr>
<tr><td>4100</td><td>1/46.1e-6</td><td>1/1.4290</td></tr>
<tr><td>4200</td><td>1/47.0e-6</td><td>1/1.4566</td></tr>
<tr><td>4300</td><td>1/47.9e-6</td><td>1/1.4842</td></tr>
<tr><td>4400</td><td>1/48.8e-6</td><td>1/1.5116</td></tr>
<tr><td>4500</td><td>1/49.7e-6</td><td>1/1.5389</td></tr>
<tr><td>4600</td><td>1/50.6e-6</td><td>1/1.5661</td></tr>
<tr><td>4700</td><td>1/51.5e-6</td><td>1/1.5933</td></tr>
<tr><td>4800</td><td>1/52.3e-6</td><td>1/1.6204</td></tr>
<tr><td>4900</td><td>1/53.2e-6</td><td>1/1.6477</td></tr>
<tr><td>5000</td><td>1/54.1e-6</td><td>1/1.6750</td></tr>
  </table>

</p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model.</p></html>"));
        end Fixed;
      end Solid;
    end 'H+';

    package N2 "<html>N<sub>2</sub></html>"
      extends Modelica.Icons.Package;
      package Gas "<html>N<sub>2</sub> gas</html>"
        extends Modelica.Icons.Package;
        model Calibrated
          "Ideal or empirical correlations, w/ adjustment factors"
          extends Species(
            redeclare FCSys.Characteristics.N2.Gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare Q.Resistivity beta_Phi=k_beta_Phi*Data.gamma(T),
            redeclare Q.ResistivityThermal beta_S=k_beta_S*Data.gamma(T),
            redeclare Q.Resistivity gamma_N=k_gamma_N*Data.gamma(T),
            redeclare Q.Resistivity gamma_Phi=k_gamma_Phi*(if
                empiricalTransverse then Data.gamma_Phi(T) else Data.gamma(T)),

            redeclare Q.ResistivityThermal gamma_S=k_gamma_S*(if
                empiricalThermal then Data.gamma_S(T) else Data.gamma(T)));

          parameter Q.NumberAbsolute k_beta_Phi(final nominal=1) = 1
            "<html>Adjustment factor for exchange of linear momentum (<i>k</i><sub>&beta;&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal exchange (<i>k</i><sub>&beta;<i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_gamma_N(final nominal=1) = 1
            "<html>Adjustment factor for material transport (<i>k</i>&gamma;<sub><i>N</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_gamma_Phi(final nominal=1) = 1
            "<html>Adjustment factor for transport of linear momentum (<i>k</i>&gamma;<sub>&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_gamma_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal transport (<i>k</i>&gamma;<sub><i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Boolean empiricalTransverse=true
            "Use empirical correlation for transverse resistivity (otherwise, theoretical)"
            annotation (
            Evaluate=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Correlations",
              compact=true));
          parameter Boolean empiricalThermal=true
            "Use empirical correlation for thermal resistivity (otherwise, theoretical)"
            annotation (
            Evaluate=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Correlations",
              compact=true));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="N2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"));
        end Calibrated;

        model Correlated "Ideal (kinetic) or empirical correlations"
          extends Species(
            redeclare FCSys.Characteristics.N2.Gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare Q.Resistivity beta_Phi=Data.gamma(T),
            redeclare Q.ResistivityThermal beta_S=Data.gamma(T),
            redeclare Q.Resistivity gamma_N=Data.gamma(T),
            redeclare Q.Resistivity gamma_Phi=if empiricalTransverse then
                Data.gamma_Phi(T) else Data.gamma(T),
            redeclare Q.ResistivityThermal gamma_S=if empiricalThermal then
                Data.gamma_S(T) else Data.gamma(T));

          parameter Boolean empiricalTransverse=true
            "Use empirical correlation for transverse resistivity (otherwise, theoretical)"
            annotation (
            Evaluate=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Correlations",
              compact=true));
          parameter Boolean empiricalThermal=true
            "Use empirical correlation for thermal resistivity (otherwise, theoretical)"
            annotation (
            Evaluate=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Correlations",
              compact=true));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="N2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare FCSys.Characteristics.N2.Gas Data(
              b_v=[1],
              specVolPow={-1,0},
              specHeatCapPow=0,
              T_lim_c0={0,Modelica.Constants.inf},
              b_c0=[1.041e3*U.J*Data.m/(U.kg*U.K)],
              B_c0=[-300*U.K*1.041e3*U.J*Data.m/(U.kg*U.K) + Data.Deltah0_f, 0]),

            redeclare parameter Q.Resistivity beta_Phi=Data.gamma(),
            redeclare parameter Q.ResistivityThermal beta_S=U.m*U.K/(25.9e-3*U.W),

            redeclare parameter Q.Resistivity gamma_N=Data.gamma(),
            redeclare parameter Q.Resistivity gamma_Phi=Data.m/(178.2e-7*U.Pa*U.s),

            redeclare parameter Q.ResistivityThermal gamma_S=U.m*U.K/(25.9e-3*U.W));

          // See the documentation for a table of values.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="N2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
    <li>Fixed specific heat capacity (independent of temperature)</ol>
    <li>Transport and exchange properties (&beta;<sub>&Phi;</sub>, &gamma;<sub>&Phi;</sub>, etc.) are fixed (e.g., independent of temperature)</li></p>

<p>The default specific heat capacity (<code>b_c0=[1.041e3*U.J*Data.m/(U.kg*U.K)]</code>) and transport resistivities (<code>gamma_Phi=Data.m/(178.2e-7*U.Pa*U.s)</code> and <code>gamma_S=U.m*U.K/(25.9e-3*U.W))</code>) are based on data of gas at 1 atm and
  300 K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 920].
  Table 1 lists the properties at  other temperatures. Note that the value for specific heat capacity at constant pressure at
  800 K (<code>c=1.22e3*U.J*Data.m/(U.kg*U.K)</code>) seems unusual, but it matches the
  reference.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"bottom\"><b>Table 1:</b> Properties of N<sub>2</sub> gas at 1 atm (<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 920)</caption>
  <tr>
      <th valign=\"middle\"><code>T/U.K</code></th>
      <th width=1><code>c*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1><code>gamma_Phi<br>*U.Pa*U.s/Data.m</code></th>
      <th width=1><code>gamma_S*U.W<br>/(U.m*U.K)</code></th>
    </tr>
<tr><td>100</td><td>1.070e3</td><td>1/68.8e-7</td><td>1/9.58e-3</td></tr>
<tr><td>150</td><td>1.050e3</td><td>1/100.6e-7</td><td>1/13.9e-3</td></tr>
<tr><td>200</td><td>1.043e3</td><td>1/129.2e-7</td><td>1/18.3e-3</td></tr>
<tr><td>250</td><td>1.042e3</td><td>1/154.9e-7</td><td>1/22.2e-3</td></tr>
<tr><td>300</td><td>1.041e3</td><td>1/178.2e-7</td><td>1/25.9e-3</td></tr>
<tr><td>350</td><td>1.042e3</td><td>1/200.0e-7</td><td>1/29.3e-3</td></tr>
<tr><td>400</td><td>1.045e3</td><td>1/220.4e-7</td><td>1/32.7e-3</td></tr>
<tr><td>450</td><td>1.050e3</td><td>1/239.6e-7</td><td>1/35.8e-3</td></tr>
<tr><td>500</td><td>1.056e3</td><td>1/257.7e-7</td><td>1/38.9e-3</td></tr>
<tr><td>550</td><td>1.065e3</td><td>1/274.7e-7</td><td>1/41.7e-3</td></tr>
<tr><td>600</td><td>1.075e3</td><td>1/290.8e-7</td><td>1/44.6e-3</td></tr>
<tr><td>700</td><td>1.098e3</td><td>1/320.1e-7</td><td>1/49.9e-3</td></tr>
<tr><td>800</td><td>1.220e3</td><td>1/349.1e-7</td><td>1/54.8e-3</td></tr>
<tr><td>900</td><td>1.146e3</td><td>1/375.3e-7</td><td>1/59.7e-3</td></tr>
<tr><td>1000</td><td>1.167e3</td><td>1/399.9e-7</td><td>1/64.7e-3</td></tr>
<tr><td>1100</td><td>1.187e3</td><td>1/423.2e-7</td><td>1/70.0e-3</td></tr>
<tr><td>1200</td><td>1.204e3</td><td>1/445.3e-7</td><td>1/75.8e-3</td></tr>
<tr><td>1300</td><td>1.219e3</td><td>1/466.2e-7</td><td>1/81.0e-3</td></tr>
  </table>

  <p>The transverse resistivity of air at 15.0 &deg;C and 1 atm is given by
       <code>gamma_Phi=Data.m*U.s/(178e-7*U.Pa)</code>
   (<a href=\"http://en.wikipedia.org/wiki/Viscosity\">http://en.wikipedia.org/wiki/Viscosity</a>).</p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"));
        end Fixed;
      end Gas;
    end N2;

    package O2 "<html>O<sub>2</sub></html>"
      extends Modelica.Icons.Package;
      package Gas "<html>O<sub>2</sub> gas</html>"
        extends Modelica.Icons.Package;
        model Calibrated
          "Ideal or empirical correlations, w/ adjustment factors"
          extends Species(
            redeclare FCSys.Characteristics.O2.Gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare Q.Resistivity beta_Phi=k_beta_Phi*Data.gamma(T),
            redeclare Q.ResistivityThermal beta_S=k_beta_S*Data.gamma(T),
            redeclare Q.Resistivity gamma_N=k_gamma_N*Data.gamma(T),
            redeclare Q.Resistivity gamma_Phi=k_gamma_Phi*(if
                empiricalTransverse then Data.gamma_Phi(T) else Data.gamma(T)),

            redeclare Q.ResistivityThermal gamma_S=k_gamma_S*(if
                empiricalThermal then Data.gamma_S(T) else Data.gamma(T)));

          parameter Q.NumberAbsolute k_beta_Phi(final nominal=1) = 1
            "<html>Adjustment factor for exchange of linear momentum (<i>k</i><sub>&beta;&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal exchange (<i>k</i><sub>&beta;<i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_gamma_N(final nominal=1) = 1
            "<html>Adjustment factor for material transport (<i>k</i>&gamma;<sub><i>N</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_gamma_Phi(final nominal=1) = 1
            "<html>Adjustment factor for transport of linear momentum (<i>k</i>&gamma;<sub>&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_gamma_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal transport (<i>k</i>&gamma;<sub><i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Boolean empiricalTransverse=true
            "Use empirical correlation for transverse resistivity (otherwise, theoretical)"
            annotation (
            Evaluate=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Correlations",
              compact=true));
          parameter Boolean empiricalThermal=true
            "Use empirical correlation for thermal resistivity (otherwise, theoretical)"
            annotation (
            Evaluate=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Correlations",
              compact=true));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="O2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"));
        end Calibrated;

        model Correlated "Ideal (kinetic) or empirical correlations"
          extends Species(
            redeclare FCSys.Characteristics.O2.Gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare Q.Resistivity beta_Phi=Data.gamma(T),
            redeclare Q.ResistivityThermal beta_S=Data.gamma(T),
            redeclare Q.Resistivity gamma_N=Data.gamma(T),
            redeclare Q.Resistivity gamma_Phi=if empiricalTransverse then
                Data.gamma_Phi(T) else Data.gamma(T),
            redeclare Q.ResistivityThermal gamma_S=if empiricalThermal then
                Data.gamma_S(T) else Data.gamma(T));

          parameter Boolean empiricalTransverse=true
            "Use empirical correlation for transverse resistivity (otherwise, theoretical)"
            annotation (
            Evaluate=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Correlations",
              compact=true));
          parameter Boolean empiricalThermal=true
            "Use empirical correlation for thermal resistivity (otherwise, theoretical)"
            annotation (
            Evaluate=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              tab="Assumptions",
              group="Correlations",
              compact=true));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="O2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare FCSys.Characteristics.O2.Gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare parameter Q.Resistivity beta_Phi=Data.gamma(),
            redeclare parameter Q.ResistivityThermal beta_S=U.m*U.K/(26.8e-3*U.W),

            redeclare parameter Q.Resistivity gamma_N=Data.gamma(),
            redeclare parameter Q.Resistivity gamma_Phi=Data.m/(207.2e-7*U.Pa*U.s),

            redeclare parameter Q.ResistivityThermal gamma_S=U.m*U.K/(26.8e-3*U.W));

          // See the documentation for a table of values.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="O2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
    <li>Transport and exchange properties (&beta;<sub>&Phi;</sub>, &gamma;<sub>&Phi;</sub>, etc.) are fixed (e.g., independent of temperature)</li>
    </ol></p>

<p>Additional notes:
<ul>
          <li>
  The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

  <p>The default transport resistivities (<code>gamma_Phi=Data.m/(207.2e-7*U.Pa*U.s)</code> and <code>gamma_S=U.m*U.K/(26.8e-3*U.W)</code>) are based on data of gas at 1 atm and
  300 K from Incropera and DeWitt  [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 920&ndash;921].
  Table 1 lists the properties at other temperatures.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"bottom\"><b>Table 1:</b> Properties of O<sub>2</sub> gas at 1 atm [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002, pp. 920&ndash;921</a>]</caption>
  <tr>
      <th valign=\"middle\"><code>T/U.K</code></th>
      <th width=1><code>c*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1><code>gamma_Phi<br>*U.Pa*U.s/Data.m</code></th>
      <th width=1><code>gamma_S*U.W<br>/(U.m*U.K)</code></th>
    </tr>
<tr><td>100</td><td>0.962e3</td><td>1/76.4e-7</td><td>1/9.25e-3</td></tr>
<tr><td>150</td><td>0.921e3</td><td>1/114.8e-7</td><td>1/13.8e-3</td></tr>
<tr><td>200</td><td>0.915e3</td><td>1/147.5e-7</td><td>1/18.3e-3</td></tr>
<tr><td>250</td><td>0.915e3</td><td>1/178.6e-7</td><td>1/22.6e-3</td></tr>
<tr><td>300</td><td>0.920e3</td><td>1/207.2e-7</td><td>1/26.8e-3</td></tr>
<tr><td>350</td><td>0.929e3</td><td>1/233.5e-7</td><td>1/29.6e-3</td></tr>
<tr><td>400</td><td>0.942e3</td><td>1/258.2e-7</td><td>1/33.0e-3</td></tr>
<tr><td>450</td><td>0.956e3</td><td>1/281.4e-7</td><td>1/36.3e-3</td></tr>
<tr><td>500</td><td>0.972e3</td><td>1/303.3e-7</td><td>1/41.2e-3</td></tr>
<tr><td>550</td><td>0.988e3</td><td>1/324.0e-7</td><td>1/44.1e-3</td></tr>
<tr><td>600</td><td>1.003e3</td><td>1/343.7e-7</td><td>1/47.3e-3</td></tr>
<tr><td>700</td><td>1.031e3</td><td>1/380.8e-7</td><td>1/52.8e-3</td></tr>
<tr><td>800</td><td>1.054e3</td><td>1/415.2e-7</td><td>1/58.9e-3</td></tr>
<tr><td>900</td><td>1.074e3</td><td>1/447.2e-7</td><td>1/64.9e-3</td></tr>
<tr><td>1000</td><td>1.090e3</td><td>1/477.0e-7</td><td>1/71.0e-3</td></tr>
<tr><td>1100</td><td>1.103e3</td><td>1/505.5e-7</td><td>1/75.8e-3</td></tr>
<tr><td>1200</td><td>1.115e3</td><td>1/532.5e-7</td><td>1/81.9e-3</td></tr>
<tr><td>1300</td><td>1.125e3</td><td>1/588.4e-7</td><td>1/87.1e-3</td></tr>
  </table>
</p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"));
        end Fixed;
      end Gas;
    end O2;

    model Species0Amount "Species with zero particle number"
      extends Species(
        final overrideEOS=true,
        final rho_IC=0,
        xNegative(viscousY=false, viscousZ=false),
        xPositive(viscousY=false, viscousZ=false),
        yNegative(viscousZ=false, viscousX=false),
        yPositive(viscousZ=false, viscousX=false),
        zNegative(viscousX=false, viscousY=false),
        zPositive(viscousX=false, viscousY=false),
        N(stateSelect=StateSelect.never),
        phi(each stateSelect=StateSelect.never),
        T(stateSelect=StateSelect.never));
      // Note:  StateSelect.never is necessary to avoid dynamic state selection
      // in Dymola 7.4.
      annotation (Documentation(info="<html>**e.g., purely electrical species without inductance
    <p>Assumptions:<ol>
  <li>**</li>
  </p>

  <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));
    end Species0Amount;

    model SpeciesInertStagnant "Inert and stagnant species"
      extends Species(
        redeclare final parameter Q.Resistivity gamma_N=Modelica.Constants.inf,

        redeclare final parameter Q.Resistivity gamma_Phi=0,
        final upstreamX,
        final upstreamY,
        final upstreamZ,
        final setPartNum=true,
        final setXVel=true,
        final setYVel=true,
        final setZVel=true,
        initMethPartNum=InitMethScalar.Volume,
        final initMethX=InitMethLinMom.Velocity,
        final initMethY=InitMethLinMom.Velocity,
        final initMethZ=InitMethLinMom.Velocity,
        final phi_IC=zeros(3),
        final derphi_IC,
        final I_IC,
        final derI_IC,
        xNegative(
          thermoOpt=ThermoOpt.ClosedDiabatic,
          final viscousY=false,
          final viscousZ=false),
        xPositive(
          thermoOpt=ThermoOpt.ClosedDiabatic,
          final viscousY=false,
          final viscousZ=false),
        yNegative(
          thermoOpt=ThermoOpt.ClosedDiabatic,
          final viscousZ=false,
          final viscousX=false),
        yPositive(
          thermoOpt=ThermoOpt.ClosedDiabatic,
          final viscousZ=false,
          final viscousX=false),
        zNegative(
          thermoOpt=ThermoOpt.ClosedDiabatic,
          final viscousX=false,
          final viscousY=false),
        zPositive(
          thermoOpt=ThermoOpt.ClosedDiabatic,
          final viscousX=false,
          final viscousY=false));

      // Note:  gamma_N and gamma_Phi don't matter since material and linear
      // momentum aren't transported.  upstreamX, upstreamY, and upstreamZ don't
      // matter since bulk current is zero.

      annotation (Documentation(info="<html><p>Assumptions:<ol>
  <li>Zero velocity</li>
  <li>No material exchange or transport</li</ol>
  </p>

  <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));
    end SpeciesInertStagnant;

    model Species
      "Model for single-species exchange, transport, and storage of material, linear momentum, and entropy"
      //extends FCSys.BaseClasses.Icons.Names.Top1;

      // Geometric parameters
      outer parameter Q.Length L[Axis](each min=Modelica.Constants.small)
        "Length" annotation (HideResult=true, missingInnerMessage=
            "This model should be used within the Subregion model.");
      outer parameter Q.Area A[Axis] "Cross-sectional area" annotation (
          HideResult=true, missingInnerMessage=
            "This model should be used within the Subregion model.");
      parameter Q.NumberAbsolute k[Axis](
        each min=Modelica.Constants.small,
        each final nominal=1) = {1,1,1}
        "<html>Anisotropic factor for transport (<b>k</b>)</html>"
        annotation (HideResult=true, Dialog(group="Geometry"));
      parameter Q.Length Lstar(
        min=Modelica.Constants.small,
        nominal=10*U.m,
        start=1e3*product(L)^(1/3))
        "<html>Characteristic length for exchange (<i>L</i><sup>&#9733;</sup>)</html>"
        annotation (Dialog(group="Geometry"));

      // Material properties
      replaceable FCSys.Characteristics.BaseClasses.Characteristic Data
        "Characteristic data of the species" annotation (
        Dialog(group="Material properties"),
        __Dymola_choicesAllMatching=true,
        Placement(transformation(extent={{-60,40},{-40,60}}),
            iconTransformation(extent={{-10,90},{10,110}})));

      // General assumptions
      parameter Boolean inclLin[Axis]={true,false,false}
        "true, if each component of linear momentum is included"
        annotation (Evaluate=true,Dialog(tab="Assumptions"));
      parameter Boolean overrideEOS=false
        "<html>Override the equation of state with the value of &rho;<sub>IC</sub></html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        Dialog(tab="Assumptions", compact=true),
        choices(__Dymola_checkBox=true));

      // Assumptions about upstream discretization
      parameter Boolean upstreamX=true "X" annotation (
        Evaluate=true,
        HideResult=true,
        Dialog(
          tab="Assumptions",
          group="Axes with upstream discretization",
          compact=true),
        choices(__Dymola_checkBox=true));
      parameter Boolean upstreamY=true "Y" annotation (
        Evaluate=true,
        HideResult=true,
        Dialog(
          tab="Assumptions",
          group="Axes with upstream discretization",
          compact=true),
        choices(__Dymola_checkBox=true));
      parameter Boolean upstreamZ=true "Z" annotation (
        Evaluate=true,
        HideResult=true,
        Dialog(
          tab="Assumptions",
          group="Axes with upstream discretization",
          compact=true),
        choices(__Dymola_checkBox=true));

      // Assumptions about dynamics
      parameter Boolean setPartNum=false "Particle number" annotation (
        Evaluate=true,
        Dialog(
          tab="Assumptions",
          group="Prescribed states (via initialization parameters)",
          compact=true),
        choices(__Dymola_checkBox=true));
      parameter Boolean setXVel=false "X-axis component of velocity"
        annotation (
        Evaluate=true,
        Dialog(
          tab="Assumptions",
          group="Prescribed states (via initialization parameters)",
          enable=inclLin[1],
          compact=true),
        choices(__Dymola_checkBox=true));
      parameter Boolean setYVel=false "Y-axis component of velocity"
        annotation (
        Evaluate=true,
        Dialog(
          tab="Assumptions",
          group="Prescribed states (via initialization parameters)",
          enable=inclLin[2],
          compact=true),
        choices(__Dymola_checkBox=true));
      parameter Boolean setZVel=false "Z-axis component of velocity"
        annotation (
        Evaluate=true,
        Dialog(
          tab="Assumptions",
          group="Prescribed states (via initialization parameters)",
          enable=inclLin[3],
          compact=true),
        choices(__Dymola_checkBox=true));
      parameter Boolean setTemp=false "Temperature" annotation (
        Evaluate=true,
        Dialog(
          tab="Assumptions",
          group="Prescribed states (via initialization parameters)",
          compact=true),
        choices(__Dymola_checkBox=true));

      // Initialization parameters for scalar properties
      parameter BaseClasses.InitMethScalar initMethPartNum=InitMethScalar.Pressure
        "Method of initializing the particle number" annotation (Evaluate=true,
          Dialog(tab="Initialization", group="Scalar properties"));
      parameter BaseClasses.InitMethScalar initMethTemp=InitMethScalar.Temperature
        "Method of initializing the temperature" annotation (Evaluate=true,
          Dialog(tab="Initialization", group="Scalar properties"));
      parameter Q.Amount N_IC(start=V_IC*rho_IC)
        "<html>Initial particle number (<i>N</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Scalar properties"));
      // Note:  This parameter is left enabled even if it isn't used to
      // explicitly initialize any states, since it is used as a guess value.
      // Similar notes apply to some other initial conditions below.
      parameter Q.Current derN_IC=0
        "<html>Initial rate of particle number ((&part;<i>N</i>/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Scalar properties",
          enable=initMethPartNum == 3 or initMethTemp == 3));
      // Note:  Dymola 7.4 doesn't recognize enumerations in the dialog enable
      // option, e.g.,
      //     enable=initMethPartNum == InitMethScalar.AmountRate.
      // Therefore, the values of the enumerations are specified numerically for
      // this initial condition and some others below.
      parameter Q.AmountVolumic rho_IC(min=if overrideEOS then 0 else Modelica.Constants.small,
          start=1/Data.v_pT(p_IC, T_IC))
        "<html>Initial volumic amount (&rho;<sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Scalar properties"));
      parameter Q.AmountVolumicRate derrho_IC=0
        "<html>Initial rate of volumic amount ((&part;&rho;/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Scalar properties",
          enable=initMethPartNum == 5 or initMethTemp == 5));
      parameter Q.Volume V_IC(start=product(L))
        "<html>Initial volume (<i>V</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Scalar properties"));
      parameter Q.VolumeRate derV_IC=0
        "<html>Initial rate of volume ((&part;<i>V</i>/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Scalar properties",
          enable=initMethPartNum == 7 or initMethTemp == 7));
      parameter Q.PressureAbsolute p_IC(start=defaults.p)
        "<html>Initial pressure (<i>p</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Scalar properties"));
      parameter Q.PressureRate derp_IC=0
        "<html>Initial rate of pressure ((&part;<i>p</i>/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Scalar properties",
          enable=initMethPartNum == 9 or initMethTemp == 9));
      parameter Q.TemperatureAbsolute T_IC(nominal=298.15*U.K, start=defaults.T)
        "<html>Initial temperature (<i>T</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Scalar properties"));
      parameter Q.TemperatureRate derT_IC=0
        "<html>Initial rate of temperature (&part;<i>T</i>/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Scalar properties",
          enable=initMethPartNum == 11 or initMethTemp == 11));
      parameter Q.NumberAbsolute s_IC(min=Modelica.Constants.small, start=
            Data.s_pT(p_IC, T_IC))
        "<html>Initial specific entropy (<i>s</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization", group="Scalar properties"));
      parameter Q.NumberRate ders_IC=0
        "<html>Initial rate of specific entropy ((&part;<i>s</i>/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Scalar properties",
          enable=initMethPartNum == 13 or initMethTemp == 13));
      parameter Q.Potential h_IC(start=Data.h0_T(T_IC))
        "<html>Initial specific enthalpy (<i>h</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization", group="Scalar properties"));
      parameter Q.PotentialRate derh_IC=0
        "<html>Initial rate of specific enthalpy ((&part;<i>h</i>/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Scalar properties",
          enable=initMethPartNum == 15 or initMethTemp == 15));
      parameter Q.Potential mu_IC(start=Data.g_pT(p_IC, T_IC))
        "<html>Initial electrochemical potential (&mu;<sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization", group="Scalar properties"));
      parameter Q.PotentialRate dermu_IC=0
        "<html>Initial rate of electrochemical potential ((&part;&mu;/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Scalar properties",
          enable=initMethPartNum == 17 or initMethTemp == 17));
      parameter Q.Current Ndot_IC=0
        "<html>Initial reaction rate (<i>N&#775;</i><sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Scalar properties",
          enable=initMethPartNum == 18 or initMethTemp == 18));

      // Initialization parameters for linear momentum
      parameter BaseClasses.InitMethLinMom initMethX=InitMethLinMom.Velocity
        "Method of initializing the x-axis component" annotation (Dialog(
          tab="Initialization",
          group="Linear momentum",
          enable=inclLin[1]));
      parameter BaseClasses.InitMethLinMom initMethY=InitMethLinMom.Velocity
        "Method of initializing the y-axis component" annotation (Dialog(
          tab="Initialization",
          group="Linear momentum",
          enable=inclLin[2]));
      parameter BaseClasses.InitMethLinMom initMethZ=InitMethLinMom.Velocity
        "Method of initializing the z-axis component" annotation (Dialog(
          tab="Initialization",
          group="Linear momentum",
          enable=inclLin[3]));
      // Note:  Dymola 7.4 doesn't provide pull-down lists for arrays of
      // enumerations; therefore, a parameter is used for each axis.
      parameter Q.Velocity phi_IC[Axis]={0,0,0}
        "<html>Initial velocity (<b>&phi;</b><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Linear momentum"));
      parameter Q.Acceleration derphi_IC[Axis]={0,0,0}
        "<html>Initial acceleration ((&part;<b>&phi;</b>/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Linear momentum",
          enable=initMethX == 3 or initMethY == 3 or initMethZ == 3));
      parameter Q.Current I_IC[Axis]={0,0,0}
        "<html>Initial current (<i><b>I</b></i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Linear momentum"));
      parameter Q.CurrentRate derI_IC[Axis]={0,0,0}
        "<html>Initial rate of current ((&part;<i><b>I</b></i>/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Linear momentum",
          enable=initMethX == 5 or initMethY == 5 or initMethZ == 5));

      // Preferred states
      Q.Amount N(
        nominal=1*U.mol,
        final start=N_IC,
        final fixed=false,
        stateSelect=StateSelect.prefer) "Particle number";
      // Note:  The start value for this variable (and others below) isn't fixed
      // because the related initial condition is applied in the initial
      // equation section.
      Q.Velocity phi[n_lin](
        each nominal=1*U.cm/U.s,
        final start=phi_IC[cartAxes],
        each final fixed=false,
        each stateSelect=StateSelect.prefer) "Velocity";
      Q.TemperatureAbsolute T(
        nominal=298.15*U.K,
        final start=T_IC,
        final fixed=false,
        stateSelect=StateSelect.prefer) "Temperature";

      // Aliases (for common terms)
      Q.Mass M(nominal=1*U.g, start=Data.m*N_IC) "Mass";
      Q.Volume V(
        min=0,
        final start=V_IC,
        final fixed=false) "Volume";
      Q.PressureAbsolute p(
        final start=p_IC,
        final fixed=false,
        stateSelect=StateSelect.never) "Pressure";
      // Note:  In Dymola 7.4, StateSelect.never is necessary to avoid dynamic
      // state selection.
      Q.Potential h(
        nominal=1*U.V,
        final start=h_IC,
        final fixed=false) "Specific enthalpy";
      Q.NumberAbsolute s(
        nominal=10,
        final start=s_IC,
        final fixed=false) "Specific entropy";
      Q.Potential mu(
        nominal=1*U.V,
        final start=mu_IC,
        final fixed=false) "Electrochemical potential";
      Q.Current I[n_lin](
        each nominal=1*U.A,
        final start=I_IC[cartAxes],
        each final fixed=false) "Current";

      // Material properties
      input Q.Resistivity beta_Phi(nominal=10*U.cm/U.A, start=Data.gamma(T_IC))
        "<html>Exchange resistivity for linear momentum (&beta;<sub>&Phi;</sub>)</html>"
        annotation (Dialog(group="Material properties"));
      input Q.ResistivityThermal beta_S(nominal=10*U.cm/U.A, start=Data.gamma(
            T_IC))
        "<html>Thermal exchange resistivity (&beta;<sub><i>S</i></sub>)</html>"
        annotation (Dialog(group="Material properties"));
      input Q.Resistivity gamma_N(nominal=10*U.cm/U.A, start=Data.gamma(T_IC))
        "<html>Material transport resistivity (&gamma;<sub><i>N</i></sub>)</html>"
        annotation (Dialog(group="Material properties"));
      input Q.Resistivity gamma_Phi(nominal=10*U.cm/U.A, start=Data.gamma(T_IC))
        "<html>Transport resistivity for linear momentum (&gamma;<sub>&Phi;</sub>)</html>"
        annotation (Dialog(group="Material properties"));
      input Q.ResistivityThermal gamma_S(nominal=10*U.cm/U.A, start=Data.gamma(
            T_IC))
        "<html>Thermal transport resistivity (&gamma;<sub><i>S</i></sub>)</html>"
        annotation (Dialog(group="Material properties"));
      // Note:  The properties are defined as inputs so that they can be
      // redeclared as parameters or acausal time-varying variables.

      // Auxiliary variables (for analysis)
      // ----------------------------------
      // General properties
      output Q.AmountVolumic rho(stateSelect=StateSelect.never) = N/V if
        defaults.analysis "Molar density";
      // Note:  The reciprocal, specific volume (v), isn't provided because
      // particle number (N) can be zero.
      output Q.CapacityThermalSpecific c_p(stateSelect=StateSelect.never) =
        Data.c0_T(T) if defaults.analysis "Specific heat capacity";
      output Q.CapacityThermal C_p(stateSelect=StateSelect.never) = N*c_p if
        defaults.analysis "Heat capacity";
      output Q.PressureAbsolute q[n_lin](each stateSelect=StateSelect.never) =
        Data.m*phi .* I ./ A[cartAxes] if defaults.analysis "Dynamic pressure";
      //
      // Time constants
      output Q.Time tau_exch_Phi(stateSelect=StateSelect.never) = beta_Phi*N/
        Lstar if defaults.analysis
        "Time constant for exchange of linear momentum";
      output Q.Time tau_exch_S(stateSelect=StateSelect.never) = beta_S*C_p/
        Lstar if defaults.analysis "Time constant for thermal exchange";
      output Q.Time tau_trans_N[Axis](each stateSelect=StateSelect.never) = -
        fill(gamma_N*N*p/Data.dp(
            v=1/rho,
            T=T,
            dv=1,
            dT=0), 3) ./ Lstar_trans if defaults.analysis and Data.isCompressible
        "Time constants for material transport";
      output Q.Time tau_trans_Phi[Axis](each stateSelect=StateSelect.never) =
        fill(gamma_Phi*N, 3) ./ Lstar_trans if defaults.analysis
        "Time constants for transport of linear momentum";
      output Q.Time tau_trans_S[Axis](each stateSelect=StateSelect.never) =
        fill(gamma_S*N*c_p, 3) ./ Lstar_trans if defaults.analysis
        "Time constants for thermal transport";
      //
      // Peclet numbers (only for the axes with linear momentum included; others
      // are zero)
      output Q.Number Pe_N[n_lin](each stateSelect=StateSelect.never) = I*
        gamma_N ./ Lstar_trans[cartAxes] if defaults.analysis
        "Material Peclet numbers";
      output Q.Number Pe_Phi[n_lin](each stateSelect=StateSelect.never) = I*
        gamma_Phi ./ Lstar_trans[cartAxes] if defaults.analysis
        "Peclet numbers for linear momentum";
      output Q.Number Pe_S[n_lin](each stateSelect=StateSelect.never) = I*
        gamma_S ./ Lstar_trans[cartAxes] if defaults.analysis
        "Thermal Peclet numbers";
      //
      // Bulk flow rates
      output Q.Force mphiI[n_lin, Orientation](each stateSelect=StateSelect.never)
         = {(if inclLin[cartWrap(cartAxes[axis] + orientation)] then Data.m*phi[
        linAxes[cartWrap(cartAxes[axis] + orientation)]]*I[axis] else 0) for
        orientation in Orientation, axis in 1:n_lin} if n_lin > 0 and defaults.analysis
        "Bulk rate of advection of 1st and 2nd transverse linear momentum";
      output Q.Force TsI[n_lin](each stateSelect=StateSelect.never) = T*s*I if
        defaults.analysis "Bulk rate of thermal advection";
      //
      // Linear momentum balance
      output Q.Force Ma[n_lin](each stateSelect=StateSelect.never) = M*der(phi)
        /U.s if defaults.analysis "Acceleration force (constant mass)";
      output Q.Force f_exch_adv[n_lin](each stateSelect=StateSelect.never) =
        chemical.mPhidot - Data.m*phi*chemical.Ndot if defaults.analysis
        "Acceleration force due to chemical (advective) exchange";
      output Q.Force f_exch_diff[n_lin](each stateSelect=StateSelect.never) =
        common.mPhidot + inert.mPhidot if defaults.analysis
        "Friction from other species (diffusive exchange)";
      output Q.Force f_trans_adv[n_lin](each stateSelect=StateSelect.never) =
        Data.m*({sum(if {{xNegative.thermoOpt == ThermoOpt.OpenDiabatic,
        xPositive.thermoOpt == ThermoOpt.OpenDiabatic},{yNegative.thermoOpt ==
        ThermoOpt.OpenDiabatic,yPositive.thermoOpt == ThermoOpt.OpenDiabatic},{
        zNegative.thermoOpt == ThermoOpt.OpenDiabatic,zPositive.thermoOpt ==
        ThermoOpt.OpenDiabatic}}[cartAxes[axis], side] then inSign(side)*
        Data.v_pT(p_face[cartAxes[axis], side], T_face[cartAxes[axis], side])*
        Ndot_face[cartAxes[axis], side]^2/A[cartAxes[axis]] else 0 for side in
        Side) + sum(phi_face[cartWrap(cartAxes[axis] - orientation), :,
        orientation]*Ndot_face[cartWrap(cartAxes[axis] - orientation), :] for
        orientation in Orientation) for axis in 1:n_lin} - phi*sum(Ndot_face))
        if defaults.analysis and not overrideEOS
        "Acceleration force due to material transport (dynamic pressure)";
      output Q.Force f_trans_diff[n_lin](each stateSelect=StateSelect.never) =
        {sum(if {{xNegative.thermoOpt == ThermoOpt.OpenDiabatic,xPositive.thermoOpt
         == ThermoOpt.OpenDiabatic},{yNegative.thermoOpt == ThermoOpt.OpenDiabatic,
        yPositive.thermoOpt == ThermoOpt.OpenDiabatic},{zNegative.thermoOpt ==
        ThermoOpt.OpenDiabatic,zPositive.thermoOpt == ThermoOpt.OpenDiabatic}}[
        cartAxes[axis], side] then inSign(side)*(p_face[cartAxes[axis], side]
         - p)*A[cartAxes[axis]] else 0 for side in Side) for axis in 1:n_lin}
         + {sum(Sigma(mPhidot_face[cartWrap(cartAxes[axis] - orientation), :,
        orientation]) for orientation in Orientation) for axis in 1:n_lin} if
        defaults.analysis
        "Friction from other subregions (diffusive transport; includes volume viscosity)";
      //
      // Energy balance
      output Q.Power derE(stateSelect=StateSelect.never) = (N*Data.c0_T(T)*der(
        T) - (if overrideEOS then 0 else V*der(Data.p_vT(V/N, T))) + Data.m*der(
        phi*phi)/2)/U.s if defaults.analysis
        "Rate of energy storage (internal and kinetic) at constant mass";
      output Q.Power Wdot_exch(stateSelect=StateSelect.never) = -chemical.phi*
        chemical.mPhidot - (Data.m*(chemical.hbar - phi*phi) - h)*chemical.Ndot
        if defaults.analysis
        "Rate of work (internal, flow, and kinetic) done by chemical exchange (advection)";
      output Q.Power Qdot_gen_exch(stateSelect=StateSelect.never) = common.phi*
        common.mPhidot + inert.phi*inert.mPhidot if defaults.analysis
        "Rate of heat generation due to friction with other species";
      output Q.Power Qdot_exch(stateSelect=StateSelect.never) = common.T*common.Sdot
         + inert.T*inert.Sdot if defaults.analysis
        "Rate of thermal conduction from other species";
      output Q.Power Wdot_trans(stateSelect=StateSelect.never) = -sum(sum((
        Data.h0_T(T_face[axis, side]) + Data.m*(Data.v_pT(p_face[axis, side],
        T_face[axis, side])*Ndot_face[axis, side]/A[axis])^2 + phi_face[axis,
        side, :]*phi_face[axis, side, :] - phi*phi - h)*Ndot_face[axis, side]
        for side in Side) for axis in Axis) if defaults.analysis
        "Rate of work (internal, flow, and kinetic) done by material transport (advection)";
      output Q.Power Qdot_gen_trans(stateSelect=StateSelect.never) = sum(
        phi_face .* mPhidot_face) if defaults.analysis
        "Rate of heat generation due to friction with other subregions";
      output Q.Power Qdot_trans(stateSelect=StateSelect.never) = sum(T_face .*
        Sdot_face) if defaults.analysis
        "Rate of thermal conduction from other subregions";
      // Note:  These auxiliary variables should not be used as states; the
      // structure of the problem should not change if they are included.
      FCSys.Connectors.ChemicalOutput chemical(
        final n_lin=n_lin,
        final m=Data.m,
        final formula=Data.formula,
        muPerT(final start=mu_IC/T_IC),
        phi(final start=phi_IC[cartAxes]),
        Ndot(final start=Ndot_IC,final fixed=false),
        hbar(final start=Data.h0_T(T_IC)/Data.m,final fixed=false))
        "Connector to exchange material with advection of linear momentum"
        annotation (Placement(transformation(extent={{-30,10},{-10,30}}),
            iconTransformation(extent={{-80,60},{-60,80}})));
      Connectors.Inert common(
        final n_lin=n_lin,
        phi(start=phi_IC[cartAxes]),
        T(start=T_IC),
        Sdot(start=0)) "Connector to couple velocity and temperature"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}}),
            iconTransformation(extent={{-10,-10},{10,10}})));
      FCSys.Connectors.InertDalton inert(
        final n_lin=n_lin,
        V(
          min=0,
          final start=V_IC,
          final fixed=false),
        p(final start=p_IC, final fixed=false),
        phi(start=phi_IC[cartAxes]),
        T(start=T_IC),
        Sdot(start=0))
        "Connector to add pressure and exchange linear momentum and entropy by diffusion"
        annotation (Placement(transformation(extent={{10,-30},{30,-10}}),
            iconTransformation(extent={{60,-80},{80,-60}})));

      replaceable FCSys.Connectors.FaceX xNegative(
        thermoOpt=ThermoOpt.OpenDiabatic,
        viscousY=inclLin[Axis.y],
        viscousZ=inclLin[Axis.z],
        material(final p(start=p_IC) = p_face[Axis.x, Side.n], final Ndot(start
              =I_IC[Axis.x]) = Ndot_face[Axis.x, Side.n]),
        momentumY(final phi(start=phi_IC[Axis.y]) = phi_face[Axis.x, Side.n,
            Orientation.following], final mPhidot=mPhidot_face[Axis.x, Side.n,
              Orientation.following]),
        momentumZ(final phi(start=phi_IC[Axis.z]) = phi_face[Axis.x, Side.n,
            Orientation.preceding], final mPhidot=mPhidot_face[Axis.x, Side.n,
              Orientation.preceding]),
        entropy(final T(start=T_IC) = T_face[Axis.x, Side.n],final Sdot(start=0)
             = Sdot_face[Axis.x, Side.n])) "Negative face along the x axis"
        annotation (Dialog(
          tab="Assumptions",
          group="Characteristics of the faces (click to edit)",
          __Dymola_label="xNegative",
          __Dymola_descriptionLabel=true),Placement(transformation(extent={{-50,
                -10},{-30,10}}), iconTransformation(extent={{-110,-10},{-90,10}})));
      replaceable FCSys.Connectors.FaceX xPositive(
        thermoOpt=xNegative.thermoOpt,
        viscousY=inclLin[Axis.y],
        viscousZ=inclLin[Axis.z],
        material(final p(start=p_IC) = p_face[Axis.x, Side.p], final Ndot(start
              =-I_IC[Axis.x]) = Ndot_face[Axis.x, Side.p]),
        momentumY(final phi(start=phi_IC[Axis.y]) = phi_face[Axis.x, Side.p,
            Orientation.following], final mPhidot=mPhidot_face[Axis.x, Side.p,
              Orientation.following]),
        momentumZ(final phi(start=phi_IC[Axis.z]) = phi_face[Axis.x, Side.p,
            Orientation.preceding], final mPhidot=mPhidot_face[Axis.x, Side.p,
              Orientation.preceding]),
        entropy(final T(start=T_IC) = T_face[Axis.x, Side.p],final Sdot(start=0)
             = Sdot_face[Axis.x, Side.p])) "Positive face along the x axis"
        annotation (Dialog(
          tab="Assumptions",
          group="Characteristics of the faces (click to edit)",
          __Dymola_label="xPositive",
          __Dymola_descriptionLabel=true),Placement(transformation(extent={{30,
                -10},{50,10}}), iconTransformation(extent={{90,-10},{110,10}})));

      replaceable FCSys.Connectors.FaceY yNegative(
        viscousZ=inclLin[Axis.z],
        viscousX=inclLin[Axis.x],
        material(final p(start=p_IC) = p_face[Axis.y, Side.n], final Ndot(start
              =I_IC[Axis.y]) = Ndot_face[Axis.y, Side.n]),
        momentumZ(final phi(start=phi_IC[Axis.z]) = phi_face[Axis.y, Side.n,
            Side.n], final mPhidot=mPhidot_face[Axis.y, Side.n, Orientation.following]),

        momentumX(final phi(start=phi_IC[Axis.x]) = phi_face[Axis.y, Side.n,
            Side.p], final mPhidot=mPhidot_face[Axis.y, Side.n, Orientation.preceding]),

        entropy(final T(start=T_IC) = T_face[Axis.y, Side.n],final Sdot(start=0)
             = Sdot_face[Axis.y, Side.n])) "Negative face along the y axis"
        annotation (Dialog(
          tab="Assumptions",
          group="Characteristics of the faces (click to edit)",
          __Dymola_label="yNegative",
          __Dymola_descriptionLabel=true),Placement(transformation(extent={{-10,
                -50},{10,-30}}), iconTransformation(extent={{-10,-110},{10,-90}})));

      replaceable FCSys.Connectors.FaceY yPositive(
        thermoOpt=yNegative.thermoOpt,
        viscousZ=inclLin[Axis.z],
        viscousX=inclLin[Axis.x],
        material(final p(start=p_IC) = p_face[Axis.y, Side.p], final Ndot(start
              =-I_IC[Axis.y]) = Ndot_face[Axis.y, Side.p]),
        momentumZ(final phi(start=phi_IC[Axis.z]) = phi_face[Axis.y, Side.p,
            Orientation.following], final mPhidot=mPhidot_face[Axis.y, Side.p,
              Orientation.following]),
        momentumX(final phi(start=phi_IC[Axis.x]) = phi_face[Axis.y, Side.p,
            Orientation.preceding], final mPhidot=mPhidot_face[Axis.y, Side.p,
              Orientation.preceding]),
        entropy(final T(start=T_IC) = T_face[Axis.y, Side.p],final Sdot(start=0)
             = Sdot_face[Axis.y, Side.p])) "Positive face along the y axis"
        annotation (Dialog(
          tab="Assumptions",
          group="Characteristics of the faces (click to edit)",
          __Dymola_label="yPositive",
          __Dymola_descriptionLabel=true),Placement(transformation(extent={{-10,
                30},{10,50}}), iconTransformation(extent={{-10,90},{10,110}})));

      replaceable FCSys.Connectors.FaceZ zNegative(
        viscousX=inclLin[Axis.x],
        viscousY=inclLin[Axis.y],
        material(final p(start=p_IC) = p_face[Axis.z, Side.n], final Ndot(start
              =I_IC[Axis.z]) = Ndot_face[Axis.z, Side.n]),
        momentumX(final phi(start=phi_IC[Axis.x]) = phi_face[Axis.z, Side.n,
            Orientation.following], final mPhidot=mPhidot_face[Axis.z, Side.n,
              Orientation.following]),
        momentumY(final phi(start=phi_IC[Axis.y]) = phi_face[Axis.z, Side.n,
            Orientation.preceding], final mPhidot=mPhidot_face[Axis.z, Side.n,
              Orientation.preceding]),
        entropy(final T(start=T_IC) = T_face[Axis.z, Side.n],final Sdot(start=0)
             = Sdot_face[Axis.z, Side.n])) "Negative face along the z axis"
        annotation (Dialog(
          tab="Assumptions",
          group="Characteristics of the faces (click to edit)",
          __Dymola_label="zNegative",
          __Dymola_descriptionLabel=true),Placement(transformation(extent={{10,
                10},{30,30}}), iconTransformation(extent={{60,60},{80,80}})));
      replaceable FCSys.Connectors.FaceZ zPositive(
        thermoOpt=zNegative.thermoOpt,
        viscousX=inclLin[Axis.x],
        viscousY=inclLin[Axis.y],
        material(final p(start=p_IC) = p_face[Axis.z, Side.p], final Ndot(start
              =-I_IC[Axis.z]) = Ndot_face[Axis.z, Side.p]),
        momentumX(final phi(start=phi_IC[Axis.x]) = phi_face[Axis.z, Side.p,
            Orientation.following], final mPhidot=mPhidot_face[Axis.z, Side.p,
              Orientation.following]),
        momentumY(final phi(start=phi_IC[Axis.y]) = phi_face[Axis.z, Side.p,
            Orientation.preceding], final mPhidot=mPhidot_face[Axis.z, Side.p,
              Orientation.preceding]),
        entropy(final T(start=T_IC) = T_face[Axis.z, Side.p],final Sdot(start=0)
             = Sdot_face[Axis.z, Side.p])) "Positive face along the z axis"
        annotation (Dialog(
          tab="Assumptions",
          group="Characteristics of the faces (click to edit)",
          __Dymola_label="zPositive",
          __Dymola_descriptionLabel=true),Placement(transformation(extent={{-30,
                -30},{-10,-10}}), iconTransformation(extent={{-80,-80},{-60,-60}})));
      // Note:  These connectors are replaceable so that their parameters can be
      // edited directly in the parameter dialog.

      // Geometric parameters
    protected
      final parameter Q.Length Lstar_trans[Axis]=k .* A ./ L
        "Effective cross-sectional area per length";
      final parameter Integer n_lin=countTrue(inclLin)
        "Number of components of linear momentum";
      final parameter Integer cartAxes[n_lin]=index(inclLin)
        "Cartesian-axis indices of the axes of linear momentum";
      final parameter Integer linAxes[Axis]=enumerate(inclLin)
        "Linear momentum indices of the Cartesian axes";
      final parameter Boolean upstream[Axis]={upstreamX,upstreamY,upstreamZ}
        "true, if each Cartesian axis uses upstream discretization";
      final parameter Boolean setVel[Axis]={setXVel,setYVel,setZVel}
        "true, if each component of velocity is prescribed";
      final parameter FCSys.Subregions.Species.BaseClasses.InitMethLinMom
        initMethLin[Axis]={initMethX,initMethY,initMethZ}
        "Initialization methods for linear momentum";

      // Efforts and flows of the conditional faces
      Q.Potential p_face[Axis, Side](each start=p_IC)
        "Static pressures at the faces";
      Q.Current Ndot_face[Axis, Side](start=outerProduct(I_IC, {1,-1}))
        "Currents into the faces";
      Q.Velocity phi_face[Axis, Side, Orientation](start={fill({phi_IC[cartWrap(
            axis + orientation)] for orientation in Orientation}, 2) for axis
             in Axis}) "Transverse velocities at the faces";
      Q.Force mPhidot_face[Axis, Side, Orientation]
        "Transverse forces on the faces";
      Q.TemperatureAbsolute T_face[Axis, Side](each start=T_IC)
        "Temperatures at the faces";
      Q.Current Sdot_face[Axis, Side] "Entropy flow rates into the faces";

      outer FCSys.BCs.Defaults defaults "Default settings" annotation (
          missingInnerMessage="Your model is using an outer \"defaults\" record, but an inner \"defaults\" record is not defined.
For simulation, specify global default settings by dragging FCSys.BCs.Defaults into your model.
The default global default settings will be used for the current simulation.",
          Placement(transformation(extent={{40,40},{60,60}}),
            iconTransformation(extent={{-10,90},{10,110}})));
      // Note:  In Dymola 7.4, it is necessary to add the missing inner message
      // here to give a warning message, even though it is included in the
      // Defaults model too.

    initial equation
      // Check that the initialization methods are valid.
      assert(initMethPartNum <> initMethTemp or initMethPartNum ==
        InitMethScalar.None,
        "The initialization methods for particle number and temperature cannot be the same (unless None).");
      assert(not (overrideEOS and (initMethPartNum == InitMethScalar.AmountVolumic
         or initMethTemp == InitMethScalar.AmountVolumic)),
        "Volumic amount cannot be used as an initial or fixed condition since it is used to override the equation of state (overrideEOS = true).");
      if not Data.isCompressible then
        assert(initMethPartNum <> InitMethScalar.Pressure and initMethPartNum
           <> InitMethScalar.PressureRate or setPartNum, "The material is incompressible,
      yet the initialization method for particle number involves pressure.");
        assert(initMethTemp <> InitMethScalar.Pressure and initMethTemp <>
          InitMethScalar.PressureRate or setTemp, "The material is incompressible,
      yet the initialization method for temperature involves pressure.");
        if not Data.hasThermalExpansion then
          assert(initMethPartNum <> InitMethScalar.AmountVolumic and
            initMethPartNum <> InitMethScalar.AmountVolumicRate or setPartNum, "The material has constant density,
      yet the initialization method for particle number involves density.");
          assert(initMethTemp <> InitMethScalar.AmountVolumic and initMethTemp
             <> InitMethScalar.AmountVolumicRate or setPartNum, "The material has constant density,
      yet the initialization method for temperature involves density.");
        end if;
      end if;

      /* This is commented out because it may be annoying.
  // Warn when index reduction may be necessary.
  if abs(beta_Phi) > Modelica.Constants.small then
    Modelica.Utilities.Streams.print("Warning: The resistivity to exchange of linear momentum is zero.
    This may directly couple the velocities of species within a subregion.
    Consider setting the value of beta_Phi as final (if not already) so that index reduction may be performed.");
  end if;
  if abs(beta_S) > Modelica.Constants.small then
    Modelica.Utilities.Streams.print("Warning: The thermal resistance to exchange is zero.
    This may directly couple the temperatures of species within a subregion.
    Consider setting the value of beta_S as final (if not already) so that index reduction may be performed.");
  end if;
  if abs(gamma_N) < Modelica.Constants.small then
    Modelica.Utilities.Streams.print("Warning: The material resistance to transport is zero.
    This may directly couple the density within neighboring subregions.\nConsider setting the value of gamma_N as final (if not already) so that index reduction may be performed.");
  end if;
  if abs(gamma_Phi) > Modelica.Constants.small then
    Modelica.Utilities.Streams.print("Warning: The resistance to transport of linear momentum is zero.
    This may directly couple the velocity within neighboring subregions.\nConsider setting the value of gamma_Phi as final (if not already) so that index reduction may be performed.");
  end if;
  if abs(gamma_S) > Modelica.Constants.small then
    Modelica.Utilities.Streams.print("Warning: The thermal resistance to transport is zero.
    This may directly couple the temperature within neighboring subregions.\nConsider setting the value of gamma_S as final (if not already) so that index reduction may be performed.");
  end if;
  // Note:  According to the Modelica 3.0 specification (and later), these
  // checks should be possible using the assert() command with
  // level=AssertionLevel.warning.  However, this isn't supported in
  // Dymola 7.4.
  */

      // Material
      if setPartNum then
        // Ensure that a condition is selected, since the state is prescribed.
        assert(initMethPartNum <> InitMethScalar.None, "The state for particle number is prescribed,
    yet its condition is not defined.\nChoose a condition besides None.");
      elseif not overrideEOS or rho_IC > 0 then
        // Initialize since there's a time-varying state.
        if initMethPartNum == InitMethScalar.Amount then
          N = N_IC;
        elseif initMethPartNum == InitMethScalar.AmountRate then
          der(N)/U.s = derN_IC;
        elseif initMethPartNum == InitMethScalar.AmountVolumic then
          N/V = rho_IC;
        elseif initMethPartNum == InitMethScalar.AmountVolumicRate then
          der(N/V)/U.s = derrho_IC;
        elseif initMethPartNum == InitMethScalar.Volume then
          V = V_IC;
        elseif initMethPartNum == InitMethScalar.VolumeRate then
          der(V)/U.s = derV_IC;
        elseif initMethPartNum == InitMethScalar.Pressure then
          p = p_IC;
        elseif initMethPartNum == InitMethScalar.PressureRate then
          der(p)/U.s = derp_IC;
        elseif initMethPartNum == InitMethScalar.Temperature then
          T = T_IC;
        elseif initMethPartNum == InitMethScalar.TemperatureRate then
          der(T)/U.s = derT_IC;
        elseif initMethPartNum == InitMethScalar.SpecificEntropy then
          s = s_IC;
        elseif initMethPartNum == InitMethScalar.SpecificEntropyRate then
          der(s)/U.s = ders_IC;
        elseif initMethPartNum == InitMethScalar.SpecificEnthalpy then
          h = h_IC;
        elseif initMethPartNum == InitMethScalar.SpecificEnthalpyRate then
          der(h)/U.s = derh_IC;
        elseif initMethPartNum == InitMethScalar.PotentialElectrochemical then
          mu = mu_IC;
        elseif initMethPartNum == InitMethScalar.PotentialElectrochemicalRate
             then
          der(mu) = dermu_IC;
        elseif initMethPartNum == InitMethScalar.ReactionRate then
          chemical.Ndot = Ndot_IC;
          // Else, initMethPartNum == InitMethScalar.None; then, there are no
          // initial equations.
        end if;
      end if;

      // Linear momentum
      for axis in Axis loop
        if inclLin[axis] then
          if setVel[axis] then
            // Ensure that a condition is selected, since the state is
            // prescribed.
            assert(initMethLin[axis] <> InitMethLinMom.None,
              "The state for the " + {"x","y","z"}[axis] + "-axis component of linear momentum is prescribed,
        yet its condition is not defined.\nChoose any condition besides None.");
          elseif not overrideEOS or rho_IC > 0 then
            // Initialize since there's a time-varying state.
            if initMethLin[axis] == InitMethLinMom.Velocity then
              phi[linAxes[axis]] = phi_IC[axis];
            elseif initMethLin[axis] == InitMethLinMom.Acceleration then
              der(phi[linAxes[axis]])/U.s = derphi_IC[axis];
            elseif initMethX == InitMethLinMom.Current then
              I[linAxes[axis]] = I_IC[axis];
            elseif initMethLin[axis] == InitMethLinMom.CurrentRate then
              der(I[linAxes[axis]])/U.s = derI_IC[axis];
              // Else, initMethLin[axis] == InitMethLinMom.None; then, there are
              // no initial equations.
            end if;
          end if;
        end if;
      end for;

      // Temperature
      if setTemp then
        // Ensure that a condition is selected, since the state is prescribed.
        assert(initMethTemp <> InitMethScalar.None, "The state for temperature is prescribed,
    yet its condition is not defined.\nChoose a condition besides None.");
      elseif not overrideEOS or rho_IC > 0 then
        // Initialize since there's a time-varying state.
        if initMethTemp == InitMethScalar.Amount then
          N = N_IC;
        elseif initMethTemp == InitMethScalar.AmountRate then
          der(N)/U.s = derN_IC;
        elseif initMethPartNum == InitMethScalar.AmountVolumic then
          N/V = rho_IC;
        elseif initMethPartNum == InitMethScalar.AmountVolumicRate then
          der(N/V)/U.s = derrho_IC;
        elseif initMethTemp == InitMethScalar.Volume then
          V = V_IC;
        elseif initMethTemp == InitMethScalar.VolumeRate then
          der(V)/U.s = derV_IC;
        elseif initMethTemp == InitMethScalar.Pressure then
          p = p_IC;
        elseif initMethTemp == InitMethScalar.PressureRate then
          der(p)/U.s = derp_IC;
        elseif initMethTemp == InitMethScalar.Temperature then
          T = T_IC;
        elseif initMethTemp == InitMethScalar.TemperatureRate then
          der(T)/U.s = derT_IC;
        elseif initMethTemp == InitMethScalar.SpecificEntropy then
          s = s_IC;
        elseif initMethTemp == InitMethScalar.SpecificEntropyRate then
          der(s)/U.s = ders_IC;
        elseif initMethTemp == InitMethScalar.SpecificEnthalpy then
          h = h_IC;
        elseif initMethTemp == InitMethScalar.SpecificEnthalpyRate then
          der(h)/U.s = derh_IC;
        elseif initMethTemp == InitMethScalar.PotentialElectrochemical then
          mu = mu_IC;
        elseif initMethTemp == InitMethScalar.PotentialElectrochemicalRate then
          der(mu) = dermu_IC;
        elseif initMethTemp == InitMethScalar.ReactionRate then
          chemical.Ndot = Ndot_IC;
          // Else, initMethTemp == InitMethScalar.None; then, there are no
          // initial equations.
        end if;
      end if;

    equation
      // Aliases (only for clarity)
      phi = common.phi;
      T = common.T;
      p = inert.p;
      V = inert.V;
      mu = T*chemical.muPerT;
      h = mu + T*s;
      N*phi = L[cartAxes] .* I;
      M = Data.m*N;

      // Thermodynamic correlations
      if overrideEOS then
        N = rho_IC*V;
      elseif Data.isCompressible then
        p = Data.p_vT(V/N, T);
      else
        V = N*Data.v_pT(p, T);
      end if;
      h = Data.h0_T(T);
      s = Data.s_pT(p, T);

      // Exchange
      // --------
      // Linear momentum
      chemical.mPhidot = semiLinear(
            Data.m*chemical.Ndot,
            chemical.phi,
            phi) "Advection";
      beta_Phi*inert.mPhidot = 2*Lstar*Data.m*(inert.phi - phi) "Diffusion";
      //
      // Energy
      chemical.Hdot = semiLinear(
            chemical.Ndot,
            chemical.hbar*Data.m,
            h) "Advection";
      beta_S*inert.Sdot = 2*Lstar*(inert.T/T - 1) "Diffusion";

      // Transport
      for axis in Axis loop
        for side in Side loop
          // Material
          if [xNegative.thermoOpt == ThermoOpt.OpenDiabatic, xPositive.thermoOpt
               == ThermoOpt.OpenDiabatic; yNegative.thermoOpt == ThermoOpt.OpenDiabatic,
              yPositive.thermoOpt == ThermoOpt.OpenDiabatic; zNegative.thermoOpt
               == ThermoOpt.OpenDiabatic, zPositive.thermoOpt == ThermoOpt.OpenDiabatic]
              [axis, side] then
            p*gamma_N*(Ndot_face[axis, side] - inSign(side)*(if inclLin[axis]
               then I[linAxes[axis]] else 0)) = Lstar_trans[axis]*(p_face[axis,
              side] - p)*(if upstream[axis] and inclLin[axis] then (exp(inSign(
              side)*I[linAxes[axis]]*gamma_N/(2*Lstar_trans[axis])) + 1) else 2);
          else
            p_face[axis, side] = 0;
            Ndot_face[axis, side] = 0;
          end if;

          // Linear momentum
          for orientation in Orientation loop
            if {{{xNegative.viscousY,xNegative.viscousZ},{xPositive.viscousY,
                xPositive.viscousZ}},{{yNegative.viscousZ,yNegative.viscousX},{
                yPositive.viscousZ,yPositive.viscousX}},{{zNegative.viscousX,
                zNegative.viscousY},{zPositive.viscousX,zPositive.viscousY}}}[
                axis, side, orientation] then
              gamma_Phi*mPhidot_face[axis, side, orientation]/Data.m =
                Lstar_trans[axis]*(phi_face[axis, side, orientation] - (if
                inclLin[cartWrap(axis + orientation)] then phi[linAxes[cartWrap(
                axis + orientation)]] else 0))*(if upstream[axis] and inclLin[
                axis] then (exp(inSign(side)*I[linAxes[axis]]*gamma_Phi/(2*
                Lstar_trans[axis])) + 1) else 2);
            else
              phi_face[axis, side, orientation] = 0;
              mPhidot_face[axis, side, orientation] = 0;
            end if;
          end for;

          // Entropy
          if [xNegative.thermoOpt == ThermoOpt.OpenDiabatic or xNegative.thermoOpt
               == ThermoOpt.ClosedDiabatic, xPositive.thermoOpt == ThermoOpt.OpenDiabatic
               or xPositive.thermoOpt == ThermoOpt.ClosedDiabatic; yNegative.thermoOpt
               == ThermoOpt.OpenDiabatic or yNegative.thermoOpt == ThermoOpt.ClosedDiabatic,
              yPositive.thermoOpt == ThermoOpt.OpenDiabatic or yPositive.thermoOpt
               == ThermoOpt.ClosedDiabatic; zNegative.thermoOpt == ThermoOpt.OpenDiabatic
               or zNegative.thermoOpt == ThermoOpt.ClosedDiabatic, zPositive.thermoOpt
               == ThermoOpt.OpenDiabatic or zPositive.thermoOpt == ThermoOpt.ClosedDiabatic]
              [axis, side] then
            T*gamma_S*Sdot_face[axis, side] = 2*Lstar_trans[axis]*(T_face[axis,
              side] - T)*(if upstream[axis] and inclLin[axis] then (exp(inSign(
              side)*I[linAxes[axis]]*gamma_S/(2*Lstar_trans[axis])) + 1) else 2);
          else
            T_face[axis, side] = 0;
            Sdot_face[axis, side] = 0;
          end if;
        end for;
      end for;

      // Material dynamics
      if setPartNum then
        // Apply the IC for all time (material not conserved).
        if initMethPartNum == InitMethScalar.Amount then
          N = N_IC;
        elseif initMethPartNum == InitMethScalar.AmountRate then
          der(N)/U.s = derN_IC;
        elseif initMethPartNum == InitMethScalar.AmountVolumic then
          N/V = rho_IC;
        elseif initMethPartNum == InitMethScalar.AmountVolumicRate then
          der(N/V)/U.s = derrho_IC;
        elseif initMethPartNum == InitMethScalar.Volume then
          V = V_IC;
        elseif initMethPartNum == InitMethScalar.VolumeRate then
          der(V)/U.s = derV_IC;
        elseif initMethPartNum == InitMethScalar.Pressure then
          p = p_IC;
        elseif initMethPartNum == InitMethScalar.PressureRate then
          der(p)/U.s = derp_IC;
        elseif initMethPartNum == InitMethScalar.Temperature then
          T = T_IC;
        elseif initMethPartNum == InitMethScalar.TemperatureRate then
          der(T)/U.s = derT_IC;
        elseif initMethPartNum == InitMethScalar.SpecificEntropy then
          s = s_IC;
        elseif initMethPartNum == InitMethScalar.SpecificEntropyRate then
          der(s)/U.s = ders_IC;
        elseif initMethPartNum == InitMethScalar.SpecificEnthalpy then
          h = h_IC;
        elseif initMethPartNum == InitMethScalar.SpecificEnthalpyRate then
          der(h)/U.s = derh_IC;
        elseif initMethPartNum == InitMethScalar.PotentialElectrochemical then
          mu = mu_IC;
        elseif initMethPartNum == InitMethScalar.PotentialElectrochemicalRate
             then
          der(mu) = dermu_IC;
        else
          //if initMethPartNum == InitMethScalar.ReactionRate then
          chemical.Ndot = Ndot_IC;
          // Note:  initMethPartNum == InitMethScalar.None can't occur due to an
          // assertion.
        end if;
      else
        der(N)/U.s = chemical.Ndot + sum(Ndot_face) "Material conservation";
      end if;

      // Dynamics of linear momentum
      for axis in 1:n_lin loop
        if setVel[cartAxes[axis]] then
          // Apply the IC for all time (linear momentum isn't conserved along
          // this axis).
          if initMethLin[cartAxes[axis]] == InitMethLinMom.Velocity then
            phi[axis] = phi_IC[cartAxes[axis]];
          elseif initMethLin[cartAxes[axis]] == InitMethLinMom.Acceleration
               then
            der(phi[axis])/U.s = derphi_IC[cartAxes[axis]];
          elseif initMethX == InitMethLinMom.Current then
            I[axis] = I_IC[cartAxes[axis]];
          elseif initMethLin[cartAxes[axis]] == InitMethLinMom.CurrentRate then
            der(I[axis])/U.s = derI_IC[cartAxes[axis]];
            // Note:  initMethLin[cartAxes[axis]] == InitMethLinMom.None can't
            // occur due to an assertion.
          end if;
        else
          der(M*phi[axis])/U.s = chemical.mPhidot[axis] + common.mPhidot[axis]
             + inert.mPhidot[axis] + sum(if {{xNegative.thermoOpt == ThermoOpt.OpenDiabatic,
            xPositive.thermoOpt == ThermoOpt.OpenDiabatic},{yNegative.thermoOpt
             == ThermoOpt.OpenDiabatic,yPositive.thermoOpt == ThermoOpt.OpenDiabatic},
            {zNegative.thermoOpt == ThermoOpt.OpenDiabatic,zPositive.thermoOpt
             == ThermoOpt.OpenDiabatic}}[cartAxes[axis], side] then inSign(side)
            *((p_face[cartAxes[axis], side] - p)*A[cartAxes[axis]] + Data.m*
            Data.v_pT(p_face[cartAxes[axis], side], T_face[cartAxes[axis], side])
            *Ndot_face[cartAxes[axis], side]^2/A[cartAxes[axis]]) else 0 for
            side in Side) + sum(Data.m*phi_face[cartWrap(cartAxes[axis] -
            orientation), :, orientation]*Ndot_face[cartWrap(cartAxes[axis] -
            orientation), :] + Sigma(mPhidot_face[cartWrap(cartAxes[axis] -
            orientation), :, orientation]) for orientation in Orientation)
            "Conservation of linear momentum";
        end if;
      end for;

      // Thermal dynamics
      if setTemp then
        // Apply the IC for all time (energy not conserved).
        if initMethTemp == InitMethScalar.Amount then
          N = N_IC;
        elseif initMethTemp == InitMethScalar.AmountRate then
          der(N)/U.s = derN_IC;
        elseif initMethTemp == InitMethScalar.Volume then
          V = V_IC;
        elseif initMethPartNum == InitMethScalar.AmountVolumic then
          N/V = rho_IC;
        elseif initMethPartNum == InitMethScalar.AmountVolumicRate then
          der(N/V)/U.s = derrho_IC;
        elseif initMethTemp == InitMethScalar.VolumeRate then
          der(V)/U.s = derV_IC;
        elseif initMethTemp == InitMethScalar.Pressure then
          p = p_IC;
        elseif initMethTemp == InitMethScalar.PressureRate then
          der(p)/U.s = derp_IC;
        elseif initMethTemp == InitMethScalar.Temperature then
          T = T_IC;
        elseif initMethTemp == InitMethScalar.TemperatureRate then
          der(T)/U.s = derT_IC;
        elseif initMethTemp == InitMethScalar.SpecificEntropy then
          s = s_IC;
        elseif initMethTemp == InitMethScalar.SpecificEntropyRate then
          der(s)/U.s = ders_IC;
        elseif initMethTemp == InitMethScalar.SpecificEnthalpy then
          h = h_IC;
        elseif initMethTemp == InitMethScalar.SpecificEnthalpyRate then
          der(h)/U.s = derh_IC;
        elseif initMethTemp == InitMethScalar.PotentialElectrochemical then
          mu = mu_IC;
        elseif initMethTemp == InitMethScalar.PotentialElectrochemicalRate then
          der(mu) = dermu_IC;
        else
          //if initMethTemp == InitMethScalar.ReactionRate then
          chemical.Ndot = Ndot_IC;
          // Note:  initMethTemp == InitMethScalar.None can't occur due to an
          // assertion.
        end if;
      else
        (N*Data.c0_T(T)*der(T) - (if overrideEOS then 0 else V*der(Data.p_vT(V/
          N, T))) + der(M*phi*phi)/2)/U.s = chemical.phi*chemical.mPhidot + (
          Data.m*chemical.hbar - h)*chemical.Ndot + common.phi*common.mPhidot
           + common.T*common.Sdot + inert.phi*inert.mPhidot + inert.T*inert.Sdot
           + sum(sum(if [xNegative.thermoOpt == ThermoOpt.OpenDiabatic,
          xPositive.thermoOpt == ThermoOpt.OpenDiabatic; yNegative.thermoOpt
           == ThermoOpt.OpenDiabatic, yPositive.thermoOpt == ThermoOpt.OpenDiabatic;
          zNegative.thermoOpt == ThermoOpt.OpenDiabatic, zPositive.thermoOpt
           == ThermoOpt.OpenDiabatic][axis, side] then (Data.h0_T(T_face[axis,
          side]) + Data.m*((Data.v_pT(p_face[axis, side], T_face[axis, side])*
          Ndot_face[axis, side]/A[axis])^2 + phi_face[axis, side, :]*phi_face[
          axis, side, :]) - h)*Ndot_face[axis, side] else 0 for side in Side)
          for axis in Axis) + sum(phi_face .* mPhidot_face) + sum(T_face .*
          Sdot_face) "Conservation of energy";
        // Note:  Although it is mathematically equivalent,
        // der(Data.p_vT(V/N, T)) is used instead of der(inert.p) or der(p)
        // so that the term can be expanded to avoid dynamic state selection
        // in Dymola 7.4.
      end if;
      annotation (
        defaultComponentPrefixes="replaceable",
        Documentation(info="<html>
    <p>This model is based on the following fixed assumptions.  Other assumptions are optional via the parameters.
    <ol>
       <li>All faces are rectangular.
       <li>The material is orthorhombic.  This implies that a
          gradient which induces diffusion along an axis does not induce
          diffusion along axes orthogonal to it [<a href=\"modelica://FCSys.UsersGuide.References\">Bejan2006</a>,
          pp. 691&ndash;692].</li>
       <li>The coordinate system (x, y, z) is aligned with the principle
          axes of transport.  For example, if the species is stratified, the
          layers must be parallel to one of the planes in the rectilinear
          grid.</li>
       <li>The factors that may cause anisotropic behavior (<b><i>k</i></b>)
          are common to the transport of material, linear momentum, and
          entropy.</li>
       <li>There are no body or inertial forces (e.g., gravity).</li>
    </ol>
    </p>

    <p>Figure 1 shows the manner in which instances of
    <a href=\"modelica://FCSys.Subregions.Species\">Species</a> models (derived from this model) are
    connected within a <a href=\"modelica://FCSys.Subregions\">Subregion</a>.  The
    exchange resistances
    are internal to each instance.  In the diagram,
    <i>Q</i> is a generic quantity (particle number <i>N</i>, linear momentum <i>m</i> &Phi;, or entropy <i>S</i> ), <i>Q&#775;</i> is the flow rate of that quantity,
    and <i>q</i> is the associated effort (&mu;, &phi;, or <i>T</i> ).
    The connection of the <a href=\"modelica://FCSys.Connectors.Material\">material connectors</a> is through a
    <a href=\"modelica://FCSys.Subregions.Reaction\">Reaction</a> model.  The
    <a href=\"modelica://FCSys.Subregions.Reaction\">Reaction</a> model advects linear momentum
    and entropy between the reactant and product species; diffusion is handled separately.
    </p>

    <p align=center><img src=\"modelica://FCSys/resources/documentation/Subregions/Species/Species/exchange.png\">
<br><b>Figure 1:</b>  Exchange of a quantity (particle number, linear momentum, or entropy) among species (A, B, and C) within a subregion.</p>

    <p>Figure 2 shows the manner in which <a href=\"modelica://FCSys.Subregions.Species\">Species</a>
    instances
    of the same type are connected between neighboring <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> instances.
    Transport is similar to exchange
    except that advection and diffusion are directly coupled; both effects
    are included in the same connections.
    Upstream discretization is applied if it is enabled (via the <code>upstreamX</code>, etc. parameters).

    <p align=center><img src=\"modelica://FCSys/resources/documentation/Subregions/Species/Species/transport.png\">
<br><b>Figure 2:</b>  Transport of a quantity associated with the same chemical species between subregions (1 and 2).</p>

    <p>Within a phase, <a href=\"modelica://FCSys.Subregions.Species\">Species</a> instances are combined
    by Dalton's law (see the <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector), as shown in
    Figure 3a.  The pressures are additive, and each species is
    assumed to exist at the volume of the phase.
    Within a subregion, phases are combined by Amagat's law (see the <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector), as shown in
    Figure 3b.  The volumes are additive, and each species is
    assumed to exist at the pressure of the subregion.</p>

    <table border=\"0\" cellspacing=\"0\" cellpadding=\"2\" align=center>
      <tr align=center>
        <td align=center>
          <img src=\"modelica://FCSys/resources/documentation/Subregions/Species/Species/share_pressure.png\">
<br><b>a:</b>  Pressures of species (A, B, and C) are additive within a phase.
        </td>
        <td align=center>
          <img src=\"modelica://FCSys/resources/documentation/Subregions/Species/Species/share_volume.png\">
<br><b>b:</b>  Volumes of phases (I, II, and III) are additive within a subregion.
        </td>
      </tr>
      <tr align=center>
        <td colspan=2 align=center><b>Figure 3:</b> Methods of sharing pressure and volume.</td>
      </tr>
    </table>

    <p>**Update: The following variables reflect the actual properties of the
    species: <code>chemical.mphi</code> (specific mass times velocity), <code>chemical.Ts</code> (specific entropy times temperature),
     <code>V</code> or <code>inert.V</code> (volume), and <code>p</code> or <code>inert.p</code> (pressure).  However, due to exchange losses
    <code>chemical.mu</code> (electrochemical potential), <code>inert.phi</code> (velocity), <code>inert.T</code> (temperature)
    are generally not the bulk properties of the species.</p>

    <p> The following notes apply to the parameters:
    <ul>
    <li>The \"specific\" adjective is taken to mean a quantity divided by particle
    number.  (\"Massic\" would indicate a quantity divided by mass.)</li>
    <li>The term \"resistivity\" indicates a generalized resistivity.  Its dimension is L.T/N (where \"L\" is length, \"T\" is time,
    and \"N\" is particle number), regardless of the quantity
    being transported or exchanged.</li>
    <li>In general, if a resistivity is zero, then it should be set as <code>final</code>
    so that index reduction may be performed.  If two <a href=\"modelica://FCSys.Subregions.Species\">Species</a> instances
    are connected through their exchange connectors
    (<code>chemical</code> or <code>inert</code>) or faces (<code>xNegative</code>, <code>xPositive</code>, etc.) and both have zero resistivities for a
    quantity, then index reduction is necessary.</li>
    <li>Even if an initialization parameter is not selected for explicit use,
    it may be used a guess value.</li>
    <li>The <b><i>k</i></b> factor can be used to account for reduced cross-sectional area
    (e.g., porosity).  It affects all the transport coefficients (for material, linear momentum, and entropy) equally.
    Its components are unity by default.</li>
    <li>By default, only the x-axis component of linear momentum is included.  Also by default,
    only material and thermal transport are included through the x-axis faces and only
    x-direction momentum is included through the y- and z-axis faces.</li>
    <li>If a state is prescribed, then the
    associated initial condition (IC) will be applied for all time.  The
    corresponding conservation equation will not be imposed.
    If <code>setPartNum</code>, <code>setXVel</code>, <code>setYVel</code>, or <code>setZVel</code> is
    <code>true</code>, then there will generally be a secondary effect on the energy conservation equation
    and thus temperature.
    In that case, it may be helpful to set <code>setTemp</code> to <code>true</code> so that
    the energy conservation equation is not imposed.</li>
    <li>If a subregion does not contain any compressible species, then pressure must be prescribed.
    Set <code>setPartNum</code> to <code>true</code> and <code>initMethPartNum</code>
    to <code>InitMethScalar.Pressure</code> for the species.  In general, only one incompressible
    species can be included if there are no incompressible species.</li>
    <li>The <code>start</code> values of the initial conditions for pressure and temperature
    (<i>p</i><sub>IC</sub> and <i>T</i><sub>IC</sub>) are the global default pressure and
    temperature (via the <code>outer</code> instance of the <a href=\"modelica://FCSys.BCs.Defaults\">Defaults</a> model).
    The <code>start</code> values of the initial conditions for
    other intensive properties (&rho;<sub>IC</sub>, <i>s</i><sub>IC</sub>, <i>h</i><sub>IC</sub>, and
    &mu;<sub>IC</sub>) are related to the initial pressure and temperature
    by the characteristics of the species.  The <code>start</code> value of the
    initial condition for the extensive volume (<i>V</i><sub>IC</sub>) is the volume of the
    subregion, and the <code>start</code> value for particle number (<i>N</i><sub>IC</sub>)
    is related to it via the characteristics (in <code>Data</code>) and the initial pressure and temperature.
    In order to apply other values for any of these initial conditions,
    it may be necessary to do so before translating the model.</li>
    <li>**Update If the species has charged (i.e., is ionic) and permittivity (<code>epsilon</code>) is
    zero, then it should be set as <code>final</code> to eliminate
    the associated state.  Otherwise, errors may occur.</li>
    <li>With the <code>overrideEOS</code> parameter, it is possible to specify that
    the volumic amount (i.e., molar concentration) and thus the amount is zero.
    Set <code>overrideEOS = true</code> and <code>rho_IC = 0</code>; then,
    the states for material, linear momentum, and energy will be eliminated (leaving
    only the electrochemical double-layer state, &Delta;&mu;, if applicable).  If a species
    is included with this setting, then there must be an external reference
    for electrochemical potential (i.e., ground).  There must be at least one other
    species in the subregion or the velocity must be set (e.g., <code>setXVel = true</code>).</li>
    </p>

    <p>In order to reduce numerical error during simulation, enthalpy of formation
    (<code>Data.Deltah0_f</code>) is excluded
    from the face connectors (e.g., <code>xNegative.material.mu</code>).  There is
    no mathematical effect since the linear momentum and energy balances are adjusted accordingly.
    **Update: However, enthalpy of formation is included in the chemical connector
    (<code>chemical.mu</code>); it is necessary for the proper chemical equilibrium.</p>

    <p>In evaluating the dynamics of a phase, it is usually assumed
    that all of the species exist at the same temperature.
    The time constants that govern the temperatures/heat capacities of the species and entropy flow rates among them
    are usually
    much shorter than the time span of interest.
    This assumption can be applied in the model by setting the thermal exchange resistivities
    (&beta;<sub><i>S</i></sub>) of the species as <code>final</code> parameters equal to zero.
    Then, the translator can perform index reduction and retain only one
    state associated with temperature.  However, this will likely lead to nonlinear systems of equations
    and may reduce the performance of the simulation. Likewise, if the reaction rates are very fast with respect to the
    observed time span,
    then the reaction resistivities (&beta;<sub><i>N</i></sub>) of the species may be redeclared as
    <code>final</code> parameters and set to zero.  A similar situation applies to
    momentum exchange (&beta;<sub>&Phi;</sub>), material transport (&gamma;<sub><i>N</i></sub>),
    compressive momentum transport
    (&gamma;<sub>&#8214;</sub>), transverse momentum transport (&gamma;<sub>&Phi;</sub>),
    and thermal transport (&gamma;<sub><i>S</i></sub>).</p>

    <p>In the variables that relate to transport,
    the first index is the axis and the second index is the side.  The sides
    are ordered from negative to positive, according to the
    <a href=\"modelica://FCSys.BaseClasses.Side\">Side</a> enumeration.
    Linear momentum is additionally indexed by
    the orientation of the momentum with respect to the face.
    The orientations are ordered in Cartesian space starting with the axis after the
    normal face, according to the
    <a href=\"modelica://FCSys.BaseClasses.Orientation\">Orientation</a> enumeration.</p>
    </html>"),
        Diagram(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.1), graphics),
        Icon(graphics={Ellipse(
              extent={{-100,100},{100,-100}},
              lineColor={127,127,127},
              pattern=LinePattern.Dash,
              fillColor={225,225,225},
              fillPattern=FillPattern.Solid), Text(
              extent={{-100,20},{100,60}},
              textString="%name",
              lineColor={0,0,0})}));
    end Species;

    package BaseClasses "Base classes (not for direct use)"

      extends Modelica.Icons.BasesPackage;

      type InitMethScalar = enumeration(
          None "Do not explicitly initialize.",
          Amount "Initialize the amount.",
          AmountRate "Initialize the rate of ditto.",
          AmountVolumic "Initialize the volumic amount.",
          AmountVolumicRate "Initialize the rate of ditto.",
          Volume "Initialize the volume.",
          VolumeRate "Initialize the rate of ditto.",
          Pressure "Initialize the pressure.",
          PressureRate "Initialize the rate of ditto.",
          Temperature "Initialize the temperature.",
          TemperatureRate "Initialize the rate of ditto.",
          SpecificEntropy "Initialize the specific entropy.",
          SpecificEntropyRate "Initialize the rate of ditto.",
          SpecificEnthalpy "Initialize the specific enthalpy.",
          SpecificEnthalpyRate "Initialize the rate of ditto.",
          PotentialElectrochemical "Initialize the electrochemical potential.",

          PotentialElectrochemicalRate "Initialize the rate of ditto.",
          ReactionRate "Initialize the reaction rate.")
        "Methods of initializing scalar properties (particle number and temperature)";

      type InitMethLinMom = enumeration(
          None "Do not explicitly initialize.",
          Velocity "Initialize the velocity.",
          Acceleration "Initialize the acceleration.",
          Current "Initialize the current.",
          CurrentRate "Initialize the rate of ditto.")
        "Methods of initializing linear momentum";
    end BaseClasses;
  end Species;

  model PhaseBoundary
    "Phase boundary (adapter between Amagat and Dalton mixtures)"
    //extends FCSys.BaseClasses.Icons.Names.Top6;
    parameter Integer n_lin(
      final min=0,
      final max=3) = 1
      "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>"
      annotation (
      Evaluate=true,
      HideResult=true,
      Dialog(group="Geometry"));
    FCSys.Connectors.InertAmagat inertA(final n_lin=n_lin)
      "Connector for volume, linear momentum, and entropy&mdash;with Amagat's law"
      annotation (Placement(transformation(extent={{60,-80},{80,-60}}),
          iconTransformation(extent={{110,-130},{130,-110}})));
    Connectors.InertDalton inertD(final n_lin=n_lin)
      "Connector for volume, linear momentum, and entropy&mdash;with Dalton's law"
      annotation (Placement(transformation(extent={{30,-50},{50,-30}}),
          iconTransformation(extent={{60,-80},{80,-60}})));

  equation
    // Equal properties
    inertA.phi = inertD.phi;
    inertA.T = inertD.T;

    // Static balances
    0 = inertA.p + inertD.p "Pressure";
    0 = inertA.V + inertD.V "Volume";

    // Rate balances (without storage or generation)
    zeros(n_lin) = inertA.mPhidot + inertD.mPhidot "Linear momentum";
    0 = inertA.Sdot + inertD.Sdot "Entropy";
    annotation (
      Documentation(info="<html><p>This model is essentially an
    adapter between the <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> and
    <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connectors.  Inside a phase,
    Dalton's law is applied.  Outside, Amagat's law is used.</p>

    <p>See also the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
              100,100}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-180,-180},{180,
              180}}), graphics={
          Rectangle(
            extent={{-170,120},{170,160}},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Ellipse(
            extent={{-60,188},{60,68}},
            lineColor={127,127,127},
            startAngle=30,
            endAngle=149,
            pattern=LinePattern.Dash,
            fillPattern=FillPattern.Solid,
            fillColor={255,255,255}),
          Ellipse(
            extent={{-170,-2},{-50,-122}},
            lineColor={127,127,127},
            startAngle=149,
            endAngle=270,
            pattern=LinePattern.Dash,
            fillPattern=FillPattern.Solid,
            fillColor={255,255,255}),
          Ellipse(
            extent={{50,-2},{170,-122}},
            lineColor={127,127,127},
            startAngle=270,
            endAngle=390,
            pattern=LinePattern.Dash,
            fillPattern=FillPattern.Solid,
            fillColor={255,255,255}),
          Polygon(
            points={{51.5,159},{162,-32},{110,-122},{-110,-122},{-162,-32},{-51.5,
                159},{51.5,159}},
            smooth=Smooth.None,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Line(
            points={{51.5,159},{162,-32}},
            color={127,127,127},
            smooth=Smooth.None,
            pattern=LinePattern.Dash),
          Line(
            points={{110,-122},{-110,-122}},
            color={127,127,127},
            smooth=Smooth.None,
            pattern=LinePattern.Dash),
          Line(
            points={{-162,-32},{-51.5,159}},
            color={127,127,127},
            smooth=Smooth.None,
            pattern=LinePattern.Dash),
          Text(
            extent={{-170,120},{170,160}},
            textString="%name",
            lineColor={0,0,0})}));
  end PhaseBoundary;

  model Reaction "Model for a chemical/electrochemical reaction"
    //extends FCSys.BaseClasses.Icons.Names.Top2;

    parameter Integer n_spec(min=2) = 0 "Number of chemical species"
      annotation (Dialog(connectorSizing=true));
    // Note  The minimum is 2 for a meaningful reaction, but the default
    // must be 0 to use connectorSizing.
    parameter Integer n_lin=1
      "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>"
      annotation (Evaluate=true, HideResult=true);
    parameter Q.Number alpha(
      min=0,
      max=1) = 0.5 "<html>Symmetry factor (&alpha;)</html>";
    parameter Q.Number epsilonPerT=0
      "<html>Additional activation barrier (&epsilon;/<i>T</i> )</html>";
    input Q.Resistance betaPerLstar=10/U.A
      "<html>Quotient of resistivity and characteristic length (&beta;/<i>L</i><sup>&#9733;</sup>)</html>"
      annotation (Dialog);

    Real nu[n_spec]=Chemistry.stoich(chemical.formula)
      "Stoichiometric coefficients";
    // Note 1:  As of Modelica 3.2 and Dymola 7.4, nu can't be a parameter or
    // constant even though it isn't time-varying.  The strings that represent
    // the chemical formulas can't be passed through the connectors as
    // parameters or constants.  However, the translator should recognize that
    // these equations are static.
    // Note 2:  This is a Real variable (rather than Integer) to avoid the
    // following warning in Dymola 7.4:
    //     "Cannot differentiate discrete or record variable:
    //         [...].nu[...]
    //     with respect to time."
    Q.Current Xidot(nominal=1*U.A) "Reaction rate";

    Connectors.ChemicalInput chemical[n_spec](each final n_lin=n_lin)
      "Connector for chemical species"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

  protected
    Q.Velocity phi[n_lin](each nominal=1*U.cm/U.s,each start=0) "Velocity";
    Q.Velocity2 hbar(nominal=1*U.V*U.mol/U.g, start=0) "Massic enthalpy";

  equation
    // Equilibrium
    0 = nu*chemical.muPerT;

    // Conservation
    nu[1:n_spec]*Xidot = chemical.Ndot "Chemical material (no storage)";
    zeros(n_lin) = sum(chemical[i].mPhidot for i in 1:n_spec)
      "Linear momentum (no storage)";
    0 = sum(chemical.Hdot) "Thermal energy (no storage)";

    // Ideal mixing/upstream discretization
    // Chemical species
    for i in 1:n_spec loop
      chemical[i].mPhidot = semiLinear(
          chemical[i].m*chemical[i].Ndot,
          chemical[i].phi,
          phi) "Linear momentum";
      chemical[i].Hdot = semiLinear(
          chemical[i].m*chemical[i].Ndot,
          chemical[i].hbar,
          hbar) "Energy";
    end for;
    // TODO:  Consider using stream connectors once they are better supported
    // (some errors occurred in Dymola 7.4).

    // Note:  This model is marked as structurally incomplete.  It must have
    // zero species by default (for automatic connector sizing), but at least
    // one species is mathematically required (two for a meaningful reaction).
    annotation (
      defaultComponentName="reaction",
      structurallyIncomplete=true,
      Documentation(info="<html>
    <p>The size of the chemical connector is automatically increased each time a connection is made.
    At least two species must be connected.
    The stoichiometry is determined automatically from the chemical formulas
    of the connected species.  No intermediate species are considered. Each reaction must be
    completely and uniquely defined by the connected species.  Otherwise an error message is given.
    If you suspect a bug in the library, please report it using the
    <a href=\"modelica://FCSys.UsersGuide.Contact\">contact information</a>.</p>

    <p>**Update this: For material, this model is essentially the opposite of a standard single-species connection.
    The stoichiometric sum of the efforts (&Sigma; &nu;<sub><i>i</i></sub> &mu;<sub><i>i</i></sub>)
    is zero, which is analogous to Kirchhoff's Current Law.  The flow rates divided by the
    stoichiometric coefficients (<i>N&#775;</i><sub><i>i</i></sub> /&nu;<sub><i>i</i></sub>)
    are equal&mdash;analogous to Kirchhoff's Voltage Law.</p>

    <p>Momentum and energy are advected using the <code>semiLinear()</code> operator.  There is no diffusion;
    it is included in the inert connections among species
    (see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model).<p>

    <p>The parameter &beta;/<i>L</i><sup>&#9733;</sup> is the reciprocal of exchange current density.</p>

    <p>Assumptions:<ul>
    <li>No storage of material, linear momentum, or energy</li></ul>
    </p>
    </html>"),
      Icon(graphics={
          Rectangle(
            extent={{-140,40},{140,80}},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Text(
            extent={{-140,40},{140,80}},
            textString="%name",
            lineColor={0,0,0}),
          Ellipse(
            extent={{-80,40},{80,-40}},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            lineColor={127,127,127},
            pattern=LinePattern.Dash),
          Text(
            extent={{-100,-16},{100,-40}},
            lineColor={127,127,127},
            textString="%n_spec")}),
      Diagram(graphics));
  end Reaction;

  model Volume "Model to establish a fixed volume for phases"
    //extends FCSys.BaseClasses.Icons.Names.Top7;

    // Geometric parameters
    parameter Boolean setVolume=true
      "true, if volume is specified (otherwise, pressure)"
      annotation (Dialog(group="Geometry"), choices(__Dymola_checkBox=true));
    parameter Q.Volume V(start=1*U.cm^3) "Volume"
      annotation (Dialog(enable=setVolume));
    parameter Q.Pressure p(start=defaults.p) "Pressure"
      annotation (Dialog(enable=not setVolume));
    parameter Integer n_lin(
      final min=0,
      final max=3) = 1
      "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>"
      annotation (Evaluate=true, HideResult=true);

    FCSys.Connectors.InertAmagat inert(final n_lin=n_lin)
      "Connector for linear momentum and entropy, with shared volume"
      annotation (Placement(transformation(extent={{60,-80},{80,-60}}),
          iconTransformation(extent={{100,-120},{120,-100}})));
    outer FCSys.BCs.Defaults defaults "Default settings" annotation (Placement(
          transformation(extent={{40,40},{60,60}}), iconTransformation(extent={
              {-10,90},{10,110}})));

  equation
    // Specified volume
    if setVolume then
      V = inert.V;
    else
      p = inert.p;
    end if;

    // Rate balances (without storage or generation)
    zeros(n_lin) = inert.mPhidot "Linear momentum";
    0 = inert.Sdot "Entropy";
    annotation (
      Documentation(info="<html><p>This model uses a <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector that imposes
    additivity of volume.  In order to convert to additivity of pressure, use
    the <a href=\"modelica://FCSys.Subregions.PhaseBoundary\">PhaseBoundary</a> model.</p>

<p>By default, <code>setVolume</code> is <code>true</code>, which sets the volume to a prescribed
value (<code>V</code>).  If there are no compressible species within a subregion, it is
necessary to set <code>setVolume</code> to <code>false</code>.  Then, pressure will be set to a prescribed
value (<code>p</code>).</p>

    <p>See also the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-160,-160},{160,
              160}}), graphics={
          Rectangle(
            extent={{-160,112},{160,152}},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Polygon(
            points={{-160,60},{-60,160},{160,160},{160,-60},{60,-160},{-160,-160},
                {-160,60}},
            lineColor={127,127,127},
            smooth=Smooth.None,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.Dash),
          Text(
            extent={{-160,112},{160,152}},
            textString="%name",
            lineColor={0,0,0})}),
      Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
              100,100}}), graphics),
      Icon(graphics));
  end Volume;

  package BaseClasses "Base classes (not for direct use)"

    extends Modelica.Icons.BasesPackage;
    partial model PartialSubregion
      "Partial subregion model for multi-dimensional and multi-species storage, transport, and exchange"
      extends FCSys.BaseClasses.Icons.Names.Top3;

      // Geometric parameters
      inner parameter Q.Length L[Axis](each min=Modelica.Constants.small,start=
            ones(3)*U.cm) "<html>Length (<b>L</b>)</html>"
        annotation (Dialog(group="Geometry"));
      final inner parameter Q.Area A[Axis]={L[cartWrap(axis + 1)]*L[cartWrap(
          axis + 2)] for axis in Axis} "Cross-sectional area";
      final parameter Q.Volume V=product(L) "Volume";
      parameter Boolean setVolume=true
        "true, if volume is specified (otherwise, pressure)"
        annotation (Dialog(group="Geometry"), choices(__Dymola_checkBox=true));

      // Assumptions about components of linear momentum
      parameter Boolean inclLinX=true "X" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          tab="Assumptions",
          group="Axes with linear momentum included",
          compact=true));
      parameter Boolean inclLinY=false "Y" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          tab="Assumptions",
          group="Axes with linear momentum included",
          compact=true));
      parameter Boolean inclLinZ=false "Z" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          tab="Assumptions",
          group="Axes with linear momentum included",
          compact=true));

      // Assumptions about faces
      parameter Boolean inclXFaces=true "X" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          tab="Assumptions",
          group="Axes with faces included",
          compact=true));
      parameter Boolean inclYFaces=true "Y" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          tab="Assumptions",
          group="Axes with faces included",
          compact=true));
      parameter Boolean inclZFaces=true "Z" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          tab="Assumptions",
          group="Axes with faces included",
          compact=true));

      FCSys.Connectors.FaceBus xNegative if inclXFaces
        "Negative face along the x axis" annotation (Placement(transformation(
              extent={{-50,-10},{-30,10}}), iconTransformation(extent={{-110,-10},
                {-90,10}})));
      FCSys.Connectors.FaceBus xPositive if inclXFaces
        "Positive face along the x axis" annotation (Placement(transformation(
              extent={{30,-10},{50,10}}), iconTransformation(extent={{90,-10},{
                110,10}})));
      FCSys.Connectors.FaceBus yNegative if inclYFaces
        "Negative face along the y axis" annotation (Placement(transformation(
              extent={{-10,-50},{10,-30}}), iconTransformation(extent={{-10,-110},
                {10,-90}})));
      FCSys.Connectors.FaceBus yPositive if inclYFaces
        "Positive face along the y axis" annotation (Placement(transformation(
              extent={{-10,30},{10,50}}), iconTransformation(extent={{-10,90},{
                10,110}})));
      FCSys.Connectors.FaceBus zNegative if inclZFaces
        "Negative face along the z axis" annotation (Placement(transformation(
              extent={{10,10},{30,30}}), iconTransformation(extent={{40,40},{60,
                60}})));
      FCSys.Connectors.FaceBus zPositive if inclZFaces
        "Positive face along the z axis" annotation (Placement(transformation(
              extent={{-30,-30},{-10,-10}}), iconTransformation(extent={{-60,-60},
                {-40,-40}})));

      FCSys.Subregions.Volume volume(
        final n_lin=n_lin,
        final V=V,
        final setVolume=setVolume) "Model to establish space for species"
        annotation (Placement(transformation(extent={{-16,-16},{16,16}})));
    protected
      final parameter Integer n_lin=countTrue({inclLinX,inclLinY,inclLinZ})
        "Number of components of linear momentum" annotation (Evaluate=true);

      annotation (
        Documentation(info="<html><p>Notes:
  <ul><li>This model must be be extended so that models can be added for
  relevant species, phases, and reactions.</li>
  <li>Material will be transported between two subregions only if both of the connected faces are marked
  as open (<code>thermoOpt==ThermoOpt.OpenDiabatic</code>)
  within the instances of the matched <a href=\"modelica://FCSys.Subregions.Species\">Species</a> models.
  If either or both of the faces are closed (<code>thermoOpt==ThermoOpt.ClosedAdiabatic</code> or
  <code>thermoOpt==ThermoOpt.ClosedDiabatic</code>), then the interface will be closed.
  note applies to the viscous/inviscous and diabatic/adiabatic properties.</li>
  <li>The x-axis component of linear momentum is included by default.  At least one component must be included.</li>
  <li>By default, <code>setVolume</code> is <code>true</code>, which sets the volume to a prescribed
  value (<code>V</code>).  If there are no compressible species within a subregion, it is
  necessary to set <code>setVolume</code> to <code>false</code>.  Then, pressure will be set to a prescribed
  value (<code>p</code>).</li>
  </ul></p></html>"),
        Diagram(graphics),
        Icon(graphics={
            Line(
              points={{-100,0},{-40,0}},
              color={127,127,127},
              thickness=0.5,
              visible=inclXFaces,
              smooth=Smooth.None),
            Line(
              points={{0,-40},{0,-100}},
              color={127,127,127},
              thickness=0.5,
              visible=inclYFaces,
              smooth=Smooth.None),
            Line(
              points={{40,40},{50,50}},
              color={127,127,127},
              thickness=0.5,
              visible=inclZFaces,
              smooth=Smooth.None),
            Polygon(
              points={{-40,16},{-16,40},{40,40},{40,-16},{16,-40},{-40,-40},{-40,
                  16}},
              lineColor={127,127,127},
              smooth=Smooth.None,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Line(
              points={{-40,-40},{-16,-16}},
              color={127,127,127},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),
            Line(
              points={{-16,40},{-16,-16},{40,-16}},
              color={127,127,127},
              smooth=Smooth.None,
              pattern=LinePattern.Dash),
            Line(
              points={{-40,0},{28,0}},
              color={210,210,210},
              visible=inclXFaces,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{0,28},{0,-40}},
              color={210,210,210},
              visible=inclYFaces,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{28,0},{100,0}},
              color={127,127,127},
              thickness=0.5,
              visible=inclXFaces,
              smooth=Smooth.None),
            Line(
              points={{0,100},{0,28}},
              color={127,127,127},
              thickness=0.5,
              visible=inclYFaces,
              smooth=Smooth.None),
            Line(
              points={{-12,-12},{40,40}},
              color={210,210,210},
              visible=inclZFaces,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{-40,16},{16,16},{16,-40}},
              color={127,127,127},
              smooth=Smooth.None),
            Line(
              points={{-50,-50},{-12,-12}},
              color={127,127,127},
              thickness=0.5,
              visible=inclZFaces,
              smooth=Smooth.None),
            Polygon(
              points={{-40,16},{-16,40},{40,40},{40,-16},{16,-40},{-40,-40},{-40,
                  16}},
              lineColor={127,127,127},
              smooth=Smooth.None),
            Line(
              points={{40,40},{16,16}},
              color={127,127,127},
              smooth=Smooth.None)}));
    end PartialSubregion;
  end BaseClasses;

  annotation (Documentation(info="<html>
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
end Subregions;
