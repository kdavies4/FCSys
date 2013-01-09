within FCSys;
package Subregions
  "Control volumes with multi-species transport, exchange, and storage"
  extends Modelica.Icons.Package;
  package Examples "Examples and tests"
    extends Modelica.Icons.ExamplesPackage;
    extends FCSys.BaseClasses.Icons.PackageUnderConstruction;
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
      replaceable FCSys.BCs.FaceBus.SubregionClosed bC1(
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        graphite(final inclC=inclC, final 'incle-'='incle-'),
        ionomer(final inclC19HF37O5S=inclC19HF37O5S, final 'inclH+'='inclH+'))
        constrainedby FCSys.BCs.FaceBus.Subregion annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-20,0})));

      replaceable FCSys.BCs.FaceBus.SubregionClosed bC2(
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        graphite(final inclC=inclC, final 'incle-'='incle-'),
        ionomer(final inclC19HF37O5S=inclC19HF37O5S, final 'inclH+'='inclH+'))
        constrainedby FCSys.BCs.FaceBus.Subregion annotation (Placement(
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
          gas(H2(redeclare FCSys.BCs.Face.Material.Pressure materialBC,
                materialSpec(k=1*U.atm))),
          graphite(C(redeclare FCSys.BCs.Face.Entropy.Temperature entropyBC,
                entropySpec(k=defaults.T)), 'e-'(redeclare
                FCSys.BCs.Face.Material.Current materialBC, redeclare
                Modelica.Blocks.Sources.Ramp materialSpec(duration=1000, height
                  =2*U.A))),
          ionomer('H+'(redeclare FCSys.BCs.Face.Material.Pressure materialBC,
                materialSpec(k(start=1*U.atm))))));

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
          gas(H2O(redeclare FCSys.BCs.Face.Material.Pressure materialBC,
                materialSpec(k=1*U.atm)), O2(redeclare
                FCSys.BCs.Face.Material.Pressure materialBC, materialSpec(k=1*U.atm))),

          graphite(C(redeclare FCSys.BCs.Face.Entropy.Temperature entropyBC,
                entropySpec(k=defaults.T)), 'e-'(redeclare
                FCSys.BCs.Face.Material.Current materialBC, redeclare
                Modelica.Blocks.Sources.Ramp materialSpec(duration=1000, height
                  =-2*U.A))),
          ionomer('H+'(redeclare FCSys.BCs.Face.Material.Pressure materialBC,
                materialSpec(k(start=1*U.atm))))),
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
            p_IC=defaults.p + 1*U.kPa,
            yNegative(inviscidZ=true, inviscidX=true),
            yPositive(inviscidZ=true, inviscidX=true),
            zNegative(inviscidX=true, inviscidY=true),
            zPositive(inviscidX=true, inviscidY=true)),
          N2(
            p_IC=defaults.p + 1*U.kPa,
            yNegative(inviscidZ=true, inviscidX=true),
            yPositive(inviscidZ=true, inviscidX=true),
            zNegative(inviscidX=true, inviscidY=true),
            zPositive(inviscidX=true, inviscidY=true)),
          O2(
            p_IC=defaults.p + 1*U.kPa,
            yNegative(inviscidZ=true, inviscidX=true),
            yPositive(inviscidZ=true, inviscidX=true),
            zNegative(inviscidX=true, inviscidY=true),
            zPositive(inviscidX=true, inviscidY=true)),
          H2(
            p_IC=defaults.p + 1*U.kPa,
            yNegative(inviscidZ=true, inviscidX=true),
            yPositive(inviscidZ=true, inviscidX=true),
            zNegative(inviscidX=true, inviscidY=true),
            zPositive(inviscidX=true, inviscidY=true))),
        graphite(
          final inclC=inclC,
          final 'incle-'='incle-',
          C(V_IC=subregion1.V/4),
          'e-'(
            yNegative(inviscidZ=true, inviscidX=true),
            yPositive(inviscidZ=true, inviscidX=true),
            zNegative(inviscidX=true, inviscidY=true),
            zPositive(inviscidX=true, inviscidY=true))),
        ionomer(
          final inclC19HF37O5S=inclC19HF37O5S,
          final 'inclH+'='inclH+',
          C19HF37O5S(V_IC=subregion1.V/4),
          'H+'(
            yNegative(inviscidZ=true, inviscidX=true),
            yPositive(inviscidZ=true, inviscidX=true),
            zNegative(inviscidX=true, inviscidY=true),
            zPositive(inviscidX=true, inviscidY=true))),
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
            each yNegative(inviscidZ=true,inviscidX=true),
            each yPositive(inviscidZ=true,inviscidX=true),
            each zNegative(inviscidX=true,inviscidY=true),
            each zPositive(inviscidX=true,inviscidY=true)),
          H2O(
            each p_IC=defaults.p,
            each yNegative(inviscidZ=true,inviscidX=true),
            each yPositive(inviscidZ=true,inviscidX=true),
            each zNegative(inviscidX=true,inviscidY=true),
            each zPositive(inviscidX=true,inviscidY=true)),
          N2(
            each p_IC=defaults.p,
            each yNegative(inviscidZ=true, inviscidX=true),
            each yPositive(inviscidZ=true,inviscidX=true),
            each zNegative(inviscidX=true,inviscidY=true),
            each zPositive(inviscidX=true,inviscidY=true)),
          O2(
            each p_IC=defaults.p,
            each yNegative(inviscidZ=true,inviscidX=true),
            each yPositive(inviscidZ=true,inviscidX=true),
            each zNegative(inviscidX=true,inviscidY=true),
            each zPositive(inviscidX=true,inviscidY=true))),
        graphite(
          each final inclC=inclC,
          each final 'incle-'='incle-',
          C(each V_IC=subregions[1].V/4),
          'e-'(
            each yNegative(inviscidZ=true,inviscidX=true),
            each yPositive(inviscidZ=true,inviscidX=true),
            each zNegative(inviscidX=true,inviscidY=true),
            each zPositive(inviscidX=true,inviscidY=true))),
        ionomer(
          each final inclC19HF37O5S=inclC19HF37O5S,
          each final 'inclH+'='inclH+',
          C19HF37O5S(each V_IC=subregions[1].V/4),
          'H+'(
            each yNegative(inviscidZ=true,inviscidX=true),
            each yPositive(inviscidZ=true,inviscidX=true),
            each zNegative(inviscidX=true,inviscidY=true),
            each zPositive(inviscidX=true,inviscidY=true))),
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
            p_IC=defaults.p - 1*U.kPa,
            yNegative(inviscidZ=true, inviscidX=true),
            yPositive(inviscidZ=true, inviscidX=true),
            zNegative(inviscidX=true, inviscidY=true),
            zPositive(inviscidX=true, inviscidY=true)),
          H2O(
            p_IC=defaults.p - 1*U.kPa,
            yNegative(inviscidZ=true, inviscidX=true),
            yPositive(inviscidZ=true, inviscidX=true),
            zNegative(inviscidX=true, inviscidY=true),
            zPositive(inviscidX=true, inviscidY=true)),
          N2(
            p_IC=defaults.p - 1*U.kPa,
            yNegative(inviscidZ=true, inviscidX=true),
            yPositive(inviscidZ=true, inviscidX=true),
            zNegative(inviscidX=true, inviscidY=true),
            zPositive(inviscidX=true, inviscidY=true)),
          O2(
            p_IC=defaults.p - 1*U.kPa,
            yNegative(inviscidZ=true, inviscidX=true),
            yPositive(inviscidZ=true, inviscidX=true),
            zNegative(inviscidX=true, inviscidY=true),
            zPositive(inviscidX=true, inviscidY=true))),
        graphite(
          final inclC=inclC,
          final 'incle-'='incle-',
          C(V_IC=subregion2.V/4),
          'e-'(
            yNegative(inviscidZ=true, inviscidX=true),
            yPositive(inviscidZ=true, inviscidX=true),
            zNegative(inviscidX=true, inviscidY=true),
            zPositive(inviscidX=true, inviscidY=true))),
        ionomer(
          final inclC19HF37O5S=inclC19HF37O5S,
          final 'inclH+'='inclH+',
          C19HF37O5S(V_IC=subregion2.V/4),
          'H+'(
            yNegative(inviscidZ=true, inviscidX=true),
            yPositive(inviscidZ=true, inviscidX=true),
            zNegative(inviscidX=true, inviscidY=true),
            zPositive(inviscidX=true, inviscidY=true))),
        inclLinY=false,
        inclLinZ=false,
        inclYFaces=false,
        inclZFaces=false)
        annotation (Placement(transformation(extent={{20,-10},{40,10}})));

      inner FCSys.BCs.Defaults defaults(analysis=true)
        annotation (Placement(transformation(extent={{60,20},{80,40}})));
      replaceable FCSys.BCs.FaceBus.SubregionClosed bC1 constrainedby
        FCSys.BCs.FaceBus.Subregion(
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        graphite(
          final inclC=inclC,
          final 'incle-'='incle-',
          C(isochoric=true)),
        ionomer(
          final inclC19HF37O5S=inclC19HF37O5S,
          final 'inclH+'='inclH+',
          C19HF37O5S(isochoric=true))) annotation (__Dymola_choicesFromPackage=
            true, Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-50,0})));
      replaceable FCSys.BCs.FaceBus.SubregionClosed bC2 constrainedby
        FCSys.BCs.FaceBus.Subregion(
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        graphite(
          final inclC=inclC,
          final 'incle-'='incle-',
          C(isochoric=true)),
        ionomer(
          final inclC19HF37O5S=inclC19HF37O5S,
          final 'inclH+'='inclH+',
          C19HF37O5S(isochoric=true))) annotation (__Dymola_choicesFromPackage=
            true, Placement(transformation(
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
        experiment(
          StopTime=4,
          Tolerance=1e-06,
          Algorithm="Dassl"),
        experimentSetupOutput,
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionsH2.mos"),
        Diagram(graphics));

    end SubregionsH2;

    model SubregionsCAndH2
      "<html>Test a one-dimensional array of subregions with C and H<sub>2</sub></html>"
      extends SubregionsH2(inclC=true, defaults(analysis=true));
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
        defaults(analysis=false));

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
              alpha_Qdot=Modelica.Constants.inf))),
        subregions(graphite(C(
              each initMethPartNum=InitMethScalar.Pressure,
              each V_IC=0.99*subregions[1].V,
              each T(displayUnit="degC")))),
        subregion2(graphite(C(
              initMethPartNum=InitMethScalar.Pressure,
              V_IC=0.99*subregion2.V,
              T(displayUnit="degC")))),
        redeclare FCSys.BCs.FaceBus.SubregionClosedAdiabatic bC1,
        redeclare FCSys.BCs.FaceBus.SubregionClosedAdiabatic bC2);

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
              alpha_Qdot=Modelica.Constants.inf))),
        subregions(gas(N2(each p_IC=defaults.p, phi(each displayUnit="mm/s"))),
            graphite(C(each V_IC=0.5*subregions[1].V, each T(displayUnit="degC")))),

        subregion2(gas(N2(p_IC=defaults.p, phi(displayUnit="mm/s"))), graphite(
              C(V_IC=0.5*subregion2.V, T(displayUnit="degC")))),
        redeclare FCSys.BCs.FaceBus.SubregionClosedAdiabatic bC1,
        redeclare FCSys.BCs.FaceBus.SubregionClosedAdiabatic bC2);

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
      FCSys.BCs.Chemical.Species species1(redeclare
          FCSys.Characteristics.'e-'.Graphite Data, materialBC=FCSys.BCs.Chemical.BaseClasses.BCTypeMaterial.PotentialElectrochemicalPerTemperature)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-30,-24})));

      FCSys.BCs.Chemical.Species species2(redeclare
          FCSys.Characteristics.'H+'.Gas Data, materialBC=FCSys.BCs.Chemical.BaseClasses.BCTypeMaterial.PotentialElectrochemicalPerTemperature)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={0,-24})));

      FCSys.BCs.Chemical.Species species3(
        redeclare FCSys.Characteristics.H2.Gas Data,
        materialBC=FCSys.BCs.Chemical.BaseClasses.BCTypeMaterial.Current,
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
      FCSys.BCs.Chemical.Species 'e-'(
        materialBC=FCSys.BCs.Chemical.BaseClasses.BCTypeMaterial.PotentialElectrochemicalPerTemperature,

        redeclare FCSys.Characteristics.'e-'.Graphite Data,
        final n_lin=n_lin) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-30,-24})));

      FCSys.BCs.Chemical.Species 'H+'(
        materialBC=FCSys.BCs.Chemical.BaseClasses.BCTypeMaterial.PotentialElectrochemicalPerTemperature,

        redeclare FCSys.Characteristics.'H+'.Solid Data,
        final n_lin=n_lin) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={0,-24})));

      FCSys.BCs.Chemical.Species H2(
        materialBC=FCSys.BCs.Chemical.BaseClasses.BCTypeMaterial.Current,
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
        experiment(StopTime=10),
        experimentSetupOutput);
    end SpeciesH2;

    model Specieseminus "Test a species"
      extends Modelica.Icons.Example;
      extends SpeciesH2(redeclare FCSys.Subregions.Species.'e-'.Graphite.Fixed
          species);
      FCSys.BCs.Face.BaseClasses.PartialSpecies faceBC(redeclare
          FCSys.BCs.Face.Material.Pressure materialBC) annotation (Placement(
            transformation(
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
            yNegative(slipZ=false, slipX=false),
            yPositive(slipZ=false, slipX=false),
            zNegative(slipX=false, slipY=false),
            zPositive(slipX=false, slipY=false)),
          H2(
            p_IC=1.05*defaults.p,
            yNegative(slipZ=false, slipX=false),
            yPositive(slipZ=false, slipX=false),
            zNegative(slipX=false, slipY=false),
            zPositive(slipX=false, slipY=false))),
        graphite(
          inclC=true,
          'incle-'=true,
          C(V_IC=subregion1.V/4),
          'e-'(
            yNegative(slipZ=false, slipX=false),
            yPositive(slipZ=false, slipX=false),
            zNegative(slipX=false, slipY=false),
            zPositive(slipX=false, slipY=false))),
        ionomer(
          inclC19HF37O5S=true,
          'inclH+'=true,
          C19HF37O5S(V_IC=subregion1.V/4),
          'H+'(
            yNegative(slipZ=false, slipX=false),
            yPositive(slipZ=false, slipX=false),
            zNegative(slipX=false, slipY=false),
            zPositive(slipX=false, slipY=false))),
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
            each yNegative(slipZ=false,slipX=false),
            each yPositive(slipZ=false,slipX=false),
            each zNegative(slipX=false,slipY=false),
            each zPositive(slipX=false,slipY=false))),
        each inclLinY=false,
        each inclLinZ=false,
        each inclYFaces=false,
        each inclZFaces=false) if n_x > 0
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
            yNegative(slipZ=false, slipX=false),
            yPositive(slipZ=false, slipX=false),
            zNegative(slipX=false, slipY=false),
            zPositive(slipX=false, slipY=false)),
          N2(
            p_IC=0.95*defaults.p,
            yNegative(slipZ=false, slipX=false),
            yPositive(slipZ=false, slipX=false),
            zNegative(slipX=false, slipY=false),
            zPositive(slipX=false, slipY=false)),
          O2(
            p_IC=0.95*defaults.p,
            yNegative(slipZ=false, slipX=false),
            yPositive(slipZ=false, slipX=false),
            zNegative(slipX=false, slipY=false),
            zPositive(slipX=false, slipY=false))),
        graphite(
          inclC=true,
          'incle-'=true,
          C(V_IC=subregion2.V/4),
          'e-'(
            yNegative(slipZ=false, slipX=false),
            yPositive(slipZ=false, slipX=false),
            zNegative(slipX=false, slipY=false),
            zPositive(slipX=false, slipY=false))),
        ionomer(
          inclC19HF37O5S=true,
          'inclH+'=true,
          C19HF37O5S(V_IC=subregion2.V/4),
          'H+'(
            yNegative(slipZ=false, slipX=false),
            yPositive(slipZ=false, slipX=false),
            zNegative(slipX=false, slipY=false),
            zPositive(slipX=false, slipY=false))),
        inclLinY=false,
        inclLinZ=false,
        inclYFaces=false,
        inclZFaces=false)
        annotation (Placement(transformation(extent={{20,-10},{40,10}})));

      inner FCSys.BCs.Defaults defaults(analysis=true)
        annotation (Placement(transformation(extent={{60,20},{80,40}})));
      replaceable FCSys.BCs.FaceBus.SubregionFlow bC1 constrainedby
        FCSys.BCs.FaceBus.Subregion(gas(
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
      replaceable FCSys.BCs.FaceBus.SubregionFlow bC2 constrainedby
        FCSys.BCs.FaceBus.Subregion(gas(
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

      replaceable FCSys.BCs.FaceBus.SubregionFlow ground(graphite('incle-'=true,
            'e-'(redeclare FCSys.BCs.Face.Material.Pressure materialBC,
              materialSpec(k=0)))) constrainedby FCSys.BCs.FaceBus.Subregion
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
    // Note:  This is listed above the extends clause so that it's listed
    // first in the parameter dialog.
    extends BaseClasses.PartialSubregion;

    Phases.Gas gas(inclH2O=true, final inclLin={inclLinX,inclLinY,inclLinZ})
      "Gas" annotation (Dialog(group="Phases"), Placement(transformation(extent
            ={{-10,-10},{10,10}})));
    Phases.Graphite graphite(final inclLin={inclLinX,inclLinY,inclLinZ})
      "Graphite" annotation (Dialog(group="Phases"), Placement(transformation(
            extent={{-10,-10},{10,10}})));
    Phases.Ionomer ionomer(final inclLin={inclLinX,inclLinY,inclLinZ})
      "Ionomer" annotation (Dialog(group="Phases"), Placement(transformation(
            extent={{-10,-10},{10,10}})));
    Phases.Liquid liquid(final inclLin={inclLinX,inclLinY,inclLinZ}) "Liquid"
      annotation (Dialog(group="Phases"), Placement(transformation(extent={{-10,
              -10},{10,10}})));

    FCSys.Subregions.Reaction HOR(final n_lin=n_lin, n_spec=3) if inclReact
       and (graphite.'incle-' and ionomer.'inclH+' and gas.inclH2 and not (gas.inclO2
       and gas.inclH2O)) "Hydrogen oxidation reaction"
      annotation (Placement(transformation(extent={{-50,30},{-30,50}})));
    FCSys.Subregions.Reaction ORR(final n_lin=n_lin, n_spec=4) if inclReact
       and (graphite.'incle-' and ionomer.'inclH+' and gas.inclO2 and gas.inclH2O
       and not gas.inclH2) "Oxygen reduction reaction"
      annotation (Placement(transformation(extent={{-50,30},{-30,50}})));

  protected
    FCSys.Connectors.ChemicalBusInternal chemical
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
        points={{8,6.10623e-16},{24,5.55112e-16},{40,5.55112e-16}},
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
        points={{8,6.10623e-16},{24,5.55112e-16},{40,5.55112e-16}},
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
        points={{8,6.10623e-16},{24,5.55112e-16},{40,5.55112e-16}},
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

    // Liquid
    connect(liquid.chemical, chemical) annotation (Line(
        points={{-5,5},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(liquid.inert, volume.inert) annotation (Line(
        points={{8,-8},{11,-11}},
        color={0,180,0},
        smooth=Smooth.None,
        thickness=0.5));
    connect(liquid.xNegative, xNegative.liquid) annotation (Line(
        points={{-8,6.10623e-16},{-8,1.16573e-15},{-25,1.16573e-15},{-25,
            5.55112e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.xPositive, xPositive.liquid) annotation (Line(
        points={{8,6.10623e-16},{24,5.55112e-16},{40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.yNegative, yNegative.liquid) annotation (Line(
        points={{6.10623e-16,-8.4},{6.10623e-16,-40},{5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.yPositive, yPositive.liquid) annotation (Line(
        points={{6.10623e-16,10},{-4.87687e-22,10},{-4.87687e-22,40},{
            5.55112e-16,40}},
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

      Diagram(graphics));
  end Subregion;

  model SubregionNoGraphite "Subregion with all phases except graphite"
    parameter Boolean inclReact=false "Include reaction(s), as appropriate"
      annotation (
      Evaluate=true,
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(tab="Assumptions"));
    // Note:  This is listed above the extension clause so that it's listed
    // first in the parameter dialog.
    extends BaseClasses.PartialSubregion;

    Phases.Gas gas(
      inclH2O=true,
      inclReact=inclReact,
      final inclLin={inclLinX,inclLinY,inclLinZ}) "Gas" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(group="Phases"),
      Placement(transformation(extent={{-10,-10},{10,10}})));

    Phases.Ionomer ionomer(final inclLin={inclLinX,inclLinY,inclLinZ})
      "Ionomer" annotation (Dialog(group="Phases"), Placement(transformation(
            extent={{-10,-10},{10,10}})));
    Phases.Liquid liquid(final inclLin={inclLinX,inclLinY,inclLinZ}) "Liquid"
      annotation (Dialog(group="Phases"), Placement(transformation(extent={{-10,
              -10},{10,10}})));

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

    // Liquid
    connect(liquid.inert, volume.inert) annotation (Line(
        points={{8,-8},{11,-11}},
        color={0,180,0},
        smooth=Smooth.None,
        thickness=0.5));
    connect(liquid.xNegative, xNegative.liquid) annotation (Line(
        points={{-8,6.10623e-16},{-8,1.16573e-15},{-25,1.16573e-15},{-25,
            5.55112e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.xPositive, xPositive.liquid) annotation (Line(
        points={{8,6.10623e-16},{8,5.55112e-16},{40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.yNegative, yNegative.liquid) annotation (Line(
        points={{6.10623e-16,-8.4},{6.10623e-16,-40},{5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.yPositive, yPositive.liquid) annotation (Line(
        points={{6.10623e-16,10},{-4.87687e-22,10},{-4.87687e-22,40},{
            5.55112e-16,40}},
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

      Diagram(graphics));
  end SubregionNoGraphite;

  model SubregionNoIonomer "Subregion with all phases except ionomer"
    parameter Boolean inclReact=false "Include reaction(s), as appropriate"
      annotation (
      Evaluate=true,
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(tab="Assumptions"));
    // Note:  This is listed above the extends clause so that it's listed
    // first in the parameter dialog.
    extends BaseClasses.PartialSubregion;

    Phases.Gas gas(
      inclH2O=true,
      inclReact=inclReact,
      final inclLin={inclLinX,inclLinY,inclLinZ}) "Gas" annotation (Dialog(
          group="Phases"), Placement(transformation(extent={{-10,-10},{10,10}})));

    Phases.Graphite graphite(final inclLin={inclLinX,inclLinY,inclLinZ})
      "Graphite" annotation (Dialog(group="Phases"), Placement(transformation(
            extent={{-10,-10},{10,10}})));
    Phases.Liquid liquid(final inclLin={inclLinX,inclLinY,inclLinZ}) "Liquid"
      annotation (Dialog(group="Phases"), Placement(transformation(extent={{-10,
              -10},{10,10}})));

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

    // Liquid
    connect(liquid.inert, volume.inert) annotation (Line(
        points={{8,-8},{11,-11}},
        color={0,180,0},
        smooth=Smooth.None,
        thickness=0.5));
    connect(liquid.xNegative, xNegative.liquid) annotation (Line(
        points={{-8,6.10623e-16},{-8,1.16573e-15},{-25,1.16573e-15},{-25,
            5.55112e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.xPositive, xPositive.liquid) annotation (Line(
        points={{8,6.10623e-16},{8,5.55112e-16},{40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.yNegative, yNegative.liquid) annotation (Line(
        points={{6.10623e-16,-8.4},{6.10623e-16,-40},{5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.yPositive, yPositive.liquid) annotation (Line(
        points={{6.10623e-16,10},{-4.87687e-22,10},{-4.87687e-22,40},{
            5.55112e-16,40}},
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

      Diagram(graphics));
  end SubregionNoIonomer;

  package Phases "Phases or mixtures of species"
    extends Modelica.Icons.Package;

    model Gas "Gas phase"

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
        initMethX=if initVelX and reduceVel then InitMethVelocity.None else
            InitMethVelocity.Velocity,
        initMethY=if initVelY and reduceVel then InitMethVelocity.None else
            InitMethVelocity.Velocity,
        initMethZ=if initVelZ and reduceVel then InitMethVelocity.None else
            InitMethVelocity.Velocity,
        initMethTemp=if initTemp and reduceTemp then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceVel then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceTemp then StateSelect.default else StateSelect.prefer))
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
        initMethX=if initVelX and reduceVel then InitMethVelocity.None else
            InitMethVelocity.Velocity,
        initMethY=if initVelY and reduceVel then InitMethVelocity.None else
            InitMethVelocity.Velocity,
        initMethZ=if initVelZ and reduceVel then InitMethVelocity.None else
            InitMethVelocity.Velocity,
        initMethTemp=if initTemp and reduceTemp then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceVel then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceTemp then StateSelect.default else StateSelect.prefer))
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
        initMethX=if initVelX and reduceVel then InitMethVelocity.None else
            InitMethVelocity.Velocity,
        initMethY=if initVelY and reduceVel then InitMethVelocity.None else
            InitMethVelocity.Velocity,
        initMethZ=if initVelZ and reduceVel then InitMethVelocity.None else
            InitMethVelocity.Velocity,
        initMethTemp=if initTemp and reduceTemp then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceVel then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceTemp then StateSelect.default else StateSelect.prefer))
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
        initMethX=if initVelX and reduceVel then InitMethVelocity.None else
            InitMethVelocity.Velocity,
        initMethY=if initVelY and reduceVel then InitMethVelocity.None else
            InitMethVelocity.Velocity,
        initMethZ=if initVelZ and reduceVel then InitMethVelocity.None else
            InitMethVelocity.Velocity,
        initMethTemp=if initTemp and reduceTemp then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceVel then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceTemp then StateSelect.default else StateSelect.prefer))
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
      connect(H2.common.mechanical, common.mechanical) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2.common.thermal, common.thermal) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Transport of x-axis linear momentum
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

      connect(H2.yNegative.mechanicalX, yNegative.H2.mechanicalX) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.yPositive.mechanicalX, yPositive.H2.mechanicalX) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.zNegative.mechanicalX, zNegative.H2.mechanicalX) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.zPositive.mechanicalX, zPositive.H2.mechanicalX) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Transport of y-axis linear momentum
      connect(H2.yNegative.material, yNegative.H2.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.yPositive.material, yPositive.H2.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.zNegative.mechanicalY, zNegative.H2.mechanicalY) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.zPositive.mechanicalY, zPositive.H2.mechanicalY) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.xNegative.mechanicalY, xNegative.H2.mechanicalY) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.xPositive.mechanicalY, xPositive.H2.mechanicalY) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      // Transport of z-axis linear momentum
      connect(H2.zNegative.material, zNegative.H2.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.zPositive.material, zPositive.H2.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.xNegative.mechanicalZ, xNegative.H2.mechanicalZ) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.xPositive.mechanicalZ, xPositive.H2.mechanicalZ) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.yNegative.mechanicalZ, yNegative.H2.mechanicalZ) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.yPositive.mechanicalZ, yPositive.H2.mechanicalZ) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      // Thermal transport
      connect(H2.xNegative.thermal, xNegative.H2.thermal) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.xPositive.thermal, xPositive.H2.thermal) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.yNegative.thermal, yNegative.H2.thermal) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.yPositive.thermal, yPositive.H2.thermal) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.zNegative.thermal, zNegative.H2.thermal) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.zPositive.thermal, zPositive.H2.thermal) annotation (Line(
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
      connect(H2O.common.mechanical, common.mechanical) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2O.common.thermal, common.thermal) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2O.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Transport of x-axis linear momentum
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

      connect(H2O.yNegative.mechanicalX, yNegative.H2O.mechanicalX) annotation
        (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.yPositive.mechanicalX, yPositive.H2O.mechanicalX) annotation
        (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.mechanicalX, zNegative.H2O.mechanicalX) annotation
        (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.mechanicalX, zPositive.H2O.mechanicalX) annotation
        (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Transport of y-axis linear momentum
      connect(H2O.yNegative.material, yNegative.H2O.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.yPositive.material, yPositive.H2O.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.mechanicalY, zNegative.H2O.mechanicalY) annotation
        (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.mechanicalY, zPositive.H2O.mechanicalY) annotation
        (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.xNegative.mechanicalY, xNegative.H2O.mechanicalY) annotation
        (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.mechanicalY, xPositive.H2O.mechanicalY) annotation
        (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      // Transport of z-axis linear momentum
      connect(H2O.zNegative.material, zNegative.H2O.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.material, zPositive.H2O.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xNegative.mechanicalZ, xNegative.H2O.mechanicalZ) annotation
        (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.mechanicalZ, xPositive.H2O.mechanicalZ) annotation
        (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.mechanicalZ, yNegative.H2O.mechanicalZ) annotation
        (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.mechanicalZ, yPositive.H2O.mechanicalZ) annotation
        (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      // Thermal transport
      connect(H2O.xNegative.thermal, xNegative.H2O.thermal) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.thermal, xPositive.H2O.thermal) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.thermal, yNegative.H2O.thermal) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.thermal, yPositive.H2O.thermal) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.thermal, zNegative.H2O.thermal) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.thermal, zPositive.H2O.thermal) annotation (Line(
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
      connect(N2.common.mechanical, common.mechanical) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(N2.common.thermal, common.thermal) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(N2.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Transport of x-axis linear momentum
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

      connect(N2.yNegative.mechanicalX, yNegative.N2.mechanicalX) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.yPositive.mechanicalX, yPositive.N2.mechanicalX) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.zNegative.mechanicalX, zNegative.N2.mechanicalX) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.zPositive.mechanicalX, zPositive.N2.mechanicalX) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Transport of y-axis linear momentum
      connect(N2.yNegative.material, yNegative.N2.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.yPositive.material, yPositive.N2.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.zNegative.mechanicalY, zNegative.N2.mechanicalY) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.zPositive.mechanicalY, zPositive.N2.mechanicalY) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.xNegative.mechanicalY, xNegative.N2.mechanicalY) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.xPositive.mechanicalY, xPositive.N2.mechanicalY) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      // Transport of z-axis linear momentum
      connect(N2.zNegative.material, zNegative.N2.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.zPositive.material, zPositive.N2.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.xNegative.mechanicalZ, xNegative.N2.mechanicalZ) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.xPositive.mechanicalZ, xPositive.N2.mechanicalZ) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.yNegative.mechanicalZ, yNegative.N2.mechanicalZ) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.yPositive.mechanicalZ, yPositive.N2.mechanicalZ) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      // Thermal transport
      connect(N2.xNegative.thermal, xNegative.N2.thermal) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.xPositive.thermal, xPositive.N2.thermal) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.yNegative.thermal, yNegative.N2.thermal) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.yPositive.thermal, yPositive.N2.thermal) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.zNegative.thermal, zNegative.N2.thermal) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.zPositive.thermal, zPositive.N2.thermal) annotation (Line(
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
      connect(O2.common.mechanical, common.mechanical) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(O2.common.thermal, common.thermal) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(O2.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{8,-8},{8,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Transport of x-axis linear momentum
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

      connect(O2.yNegative.mechanicalX, yNegative.O2.mechanicalX) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.yPositive.mechanicalX, yPositive.O2.mechanicalX) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.zNegative.mechanicalX, zNegative.O2.mechanicalX) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.zPositive.mechanicalX, zPositive.O2.mechanicalX) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Transport of y-axis linear momentum
      connect(O2.yNegative.material, yNegative.O2.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.yPositive.material, yPositive.O2.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.zNegative.mechanicalY, zNegative.O2.mechanicalY) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.zPositive.mechanicalY, zPositive.O2.mechanicalY) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.xNegative.mechanicalY, xNegative.O2.mechanicalY) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.xPositive.mechanicalY, xPositive.O2.mechanicalY) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      // Transport of z-axis linear momentum
      connect(O2.zNegative.material, zNegative.O2.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.zPositive.material, zPositive.O2.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.xNegative.mechanicalZ, xNegative.O2.mechanicalZ) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.xPositive.mechanicalZ, xPositive.O2.mechanicalZ) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.yNegative.mechanicalZ, yNegative.O2.mechanicalZ) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.yPositive.mechanicalZ, yPositive.O2.mechanicalZ) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      // Thermal transport
      connect(O2.xNegative.thermal, xNegative.O2.thermal) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.xPositive.thermal, xPositive.O2.thermal) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.yNegative.thermal, yNegative.O2.thermal) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.yPositive.thermal, yPositive.O2.thermal) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.zNegative.thermal, zNegative.O2.thermal) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.zPositive.thermal, zPositive.O2.thermal) annotation (Line(
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

    model Graphite "Graphite phase"
      extends BaseClasses.NullPhase(final n_spec=countTrue({inclC,'incle-'}),
          common(mechanical(phi(fixed=if inclC then fill(false, n_lin) else {
                  initVelX,initVelY,initVelZ}[cartAxes]))));

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
        initMethTemp=if initTemp and reduceTemp then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceVel then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceTemp then StateSelect.default else StateSelect.prefer))
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
        initMethX=if (initVelX or inclC) and reduceVel then InitMethVelocity.None
             else InitMethVelocity.Velocity,
        initMethY=if (initVelY or inclC) and reduceVel then InitMethVelocity.None
             else InitMethVelocity.Velocity,
        initMethZ=if (initVelZ or inclC) and reduceVel then InitMethVelocity.None
             else InitMethVelocity.Velocity,
        initMethTemp=if initTemp and reduceTemp then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceVel then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceTemp then StateSelect.default else StateSelect.prefer))
        "<html>'e-' model</html>" annotation (
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
      connect(C.common.mechanical, common.mechanical) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(C.common.thermal, common.thermal) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(C.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Transport of x-axis linear momentum
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

      connect(C.yNegative.mechanicalX, yNegative.C.mechanicalX) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.yPositive.mechanicalX, yPositive.C.mechanicalX) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.zNegative.mechanicalX, zNegative.C.mechanicalX) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.zPositive.mechanicalX, zPositive.C.mechanicalX) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Transport of y-axis linear momentum
      connect(C.yNegative.material, yNegative.C.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.yPositive.material, yPositive.C.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.zNegative.mechanicalY, zNegative.C.mechanicalY) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.zPositive.mechanicalY, zPositive.C.mechanicalY) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.xNegative.mechanicalY, xNegative.C.mechanicalY) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.xPositive.mechanicalY, xPositive.C.mechanicalY) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      // Transport of z-axis linear momentum
      connect(C.zNegative.material, zNegative.C.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.zPositive.material, zPositive.C.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.xNegative.mechanicalZ, xNegative.C.mechanicalZ) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.xPositive.mechanicalZ, xPositive.C.mechanicalZ) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.yNegative.mechanicalZ, yNegative.C.mechanicalZ) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.yPositive.mechanicalZ, yPositive.C.mechanicalZ) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      // Thermal transport
      connect(C.xNegative.thermal, xNegative.C.thermal) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.xPositive.thermal, xPositive.C.thermal) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.yNegative.thermal, yNegative.C.thermal) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.yPositive.thermal, yPositive.C.thermal) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.zNegative.thermal, zNegative.C.thermal) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.zPositive.thermal, zPositive.C.thermal) annotation (Line(
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
      connect('e-'.common.mechanical, common.mechanical) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect('e-'.common.thermal, common.thermal) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect('e-'.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Transport of x-axis linear momentum
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

      connect('e-'.yNegative.mechanicalX, yNegative.'e-'.mechanicalX)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.yPositive.mechanicalX, yPositive.'e-'.mechanicalX)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.zNegative.mechanicalX, zNegative.'e-'.mechanicalX)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.zPositive.mechanicalX, zPositive.'e-'.mechanicalX)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Transport of y-axis linear momentum
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

      connect('e-'.zNegative.mechanicalY, zNegative.'e-'.mechanicalY)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.zPositive.mechanicalY, zPositive.'e-'.mechanicalY)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.xNegative.mechanicalY, xNegative.'e-'.mechanicalY)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.xPositive.mechanicalY, xPositive.'e-'.mechanicalY)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      // Transport of z-axis linear momentum
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

      connect('e-'.xNegative.mechanicalZ, xNegative.'e-'.mechanicalZ)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.xPositive.mechanicalZ, xPositive.'e-'.mechanicalZ)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.yNegative.mechanicalZ, yNegative.'e-'.mechanicalZ)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.yPositive.mechanicalZ, yPositive.'e-'.mechanicalZ)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      // Thermal transport
      connect('e-'.xNegative.thermal, xNegative.'e-'.thermal) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.xPositive.thermal, xPositive.'e-'.thermal) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.yNegative.thermal, yNegative.'e-'.thermal) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.yPositive.thermal, yPositive.'e-'.thermal) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.zNegative.thermal, zNegative.'e-'.thermal) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.zPositive.thermal, zPositive.'e-'.thermal) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      annotation (
        defaultComponentPrefixes="replaceable",
        Documentation(info="<html><p>If C is included (<code>inclC</code>=<code>true</code>),
    then <code>initVelX</code>, <code>initVelY</code>, and <code>initVelZ</code> are
    ignored.  Velocity is not initialized because it is set by C.</p>

    <p>For more information, see the
 <a href=\"modelica://FCSys.Subregions.Phases.BaseClasses.NullPhase\">NullPhase</a> model.</p></html>"),

        Diagram(graphics));
    end Graphite;

    model Ionomer "Ionomer phase"
      extends BaseClasses.NullPhase(final n_spec=countTrue({inclC19HF37O5S,
            inclH2O,'inclH+'}), common(mechanical(phi(fixed=if inclC19HF37O5S
                   then fill(false, n_lin) else {initVelX,initVelY,initVelZ}[
                  cartAxes]))));

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
        initMethTemp=if initTemp and reduceTemp then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceVel then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceTemp then StateSelect.default else StateSelect.prefer))
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
        initMethX=if (initVelX19HF37O5S) and reduceVel then InitMethVelocity.None
             else InitMethVelocity.Velocity,
        initMethY=if (initVelY or inclC19HF37O5S) and reduceVel then
            InitMethVelocity.None else InitMethVelocity.Velocity,
        initMethZ=if (initVelZ or inclC19HF37O5S) and reduceVel then
            InitMethVelocity.None else InitMethVelocity.Velocity,
        initMethTemp=if initTemp and reduceTemp then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceVel then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceTemp then StateSelect.default else StateSelect.prefer))
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
        initMethX=if initVelX and reduceVel then InitMethVelocity.None else
            InitMethVelocity.Velocity,
        initMethY=if initVelY and reduceVel then InitMethVelocity.None else
            InitMethVelocity.Velocity,
        initMethZ=if initVelZ and reduceVel then InitMethVelocity.None else
            InitMethVelocity.Velocity,
        initMethTemp=if initTemp and reduceTemp then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceVel then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceTemp then StateSelect.default else StateSelect.prefer))
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
      connect(C19HF37O5S.common.mechanical, common.mechanical) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(C19HF37O5S.common.thermal, common.thermal) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(C19HF37O5S.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Transport of x-axis linear momentum
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

      connect(C19HF37O5S.yNegative.mechanicalX, yNegative.C19HF37O5S.mechanicalX)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.yPositive.mechanicalX, yPositive.C19HF37O5S.mechanicalX)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.zNegative.mechanicalX, zNegative.C19HF37O5S.mechanicalX)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.zPositive.mechanicalX, zPositive.C19HF37O5S.mechanicalX)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Transport of y-axis linear momentum
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

      connect(C19HF37O5S.zNegative.mechanicalY, zNegative.C19HF37O5S.mechanicalY)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.zPositive.mechanicalY, zPositive.C19HF37O5S.mechanicalY)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.xNegative.mechanicalY, xNegative.C19HF37O5S.mechanicalY)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.xPositive.mechanicalY, xPositive.C19HF37O5S.mechanicalY)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      // Transport of z-axis linear momentum
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

      connect(C19HF37O5S.xNegative.mechanicalZ, xNegative.C19HF37O5S.mechanicalZ)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.xPositive.mechanicalZ, xPositive.C19HF37O5S.mechanicalZ)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.yNegative.mechanicalZ, yNegative.C19HF37O5S.mechanicalZ)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.yPositive.mechanicalZ, yPositive.C19HF37O5S.mechanicalZ)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      // Thermal transport
      connect(C19HF37O5S.xNegative.thermal, xNegative.C19HF37O5S.thermal)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.xPositive.thermal, xPositive.C19HF37O5S.thermal)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.yNegative.thermal, yNegative.C19HF37O5S.thermal)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.yPositive.thermal, yPositive.C19HF37O5S.thermal)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.zNegative.thermal, zNegative.C19HF37O5S.thermal)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.zPositive.thermal, zPositive.C19HF37O5S.thermal)
        annotation (Line(
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
      connect('H+'.common.mechanical, common.mechanical) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect('H+'.common.thermal, common.thermal) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect('H+'.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Transport of x-axis linear momentum
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

      connect('H+'.yNegative.mechanicalX, yNegative.'H+'.mechanicalX)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.yPositive.mechanicalX, yPositive.'H+'.mechanicalX)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.zNegative.mechanicalX, zNegative.'H+'.mechanicalX)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.zPositive.mechanicalX, zPositive.'H+'.mechanicalX)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Transport of y-axis linear momentum
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

      connect('H+'.zNegative.mechanicalY, zNegative.'H+'.mechanicalY)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.zPositive.mechanicalY, zPositive.'H+'.mechanicalY)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.xNegative.mechanicalY, xNegative.'H+'.mechanicalY)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.xPositive.mechanicalY, xPositive.'H+'.mechanicalY)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      // Transport of z-axis linear momentum
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

      connect('H+'.xNegative.mechanicalZ, xNegative.'H+'.mechanicalZ)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.xPositive.mechanicalZ, xPositive.'H+'.mechanicalZ)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.yNegative.mechanicalZ, yNegative.'H+'.mechanicalZ)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.yPositive.mechanicalZ, yPositive.'H+'.mechanicalZ)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      // Thermal transport
      connect('H+'.xNegative.thermal, xNegative.'H+'.thermal) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.xPositive.thermal, xPositive.'H+'.thermal) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.yNegative.thermal, yNegative.'H+'.thermal) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.yPositive.thermal, yPositive.'H+'.thermal) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.zNegative.thermal, zNegative.'H+'.thermal) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.zPositive.thermal, zPositive.'H+'.thermal) annotation (Line(
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
      connect(H2O.common.mechanical, common.mechanical) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2O.common.thermal, common.thermal) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2O.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Transport of x-axis linear momentum
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

      connect(H2O.yNegative.mechanicalX, yNegative.H2O.mechanicalX) annotation
        (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.yPositive.mechanicalX, yPositive.H2O.mechanicalX) annotation
        (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.mechanicalX, zNegative.H2O.mechanicalX) annotation
        (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.mechanicalX, zPositive.H2O.mechanicalX) annotation
        (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Transport of y-axis linear momentum
      connect(H2O.yNegative.material, yNegative.H2O.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.yPositive.material, yPositive.H2O.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.mechanicalY, zNegative.H2O.mechanicalY) annotation
        (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.mechanicalY, zPositive.H2O.mechanicalY) annotation
        (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.xNegative.mechanicalY, xNegative.H2O.mechanicalY) annotation
        (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.mechanicalY, xPositive.H2O.mechanicalY) annotation
        (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      // Transport of z-axis linear momentum
      connect(H2O.zNegative.material, zNegative.H2O.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.material, zPositive.H2O.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xNegative.mechanicalZ, xNegative.H2O.mechanicalZ) annotation
        (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.mechanicalZ, xPositive.H2O.mechanicalZ) annotation
        (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.mechanicalZ, yNegative.H2O.mechanicalZ) annotation
        (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.mechanicalZ, yPositive.H2O.mechanicalZ) annotation
        (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      // Thermal transport
      connect(H2O.xNegative.thermal, xNegative.H2O.thermal) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.thermal, xPositive.H2O.thermal) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.thermal, yNegative.H2O.thermal) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.thermal, yPositive.H2O.thermal) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.thermal, zNegative.H2O.thermal) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.thermal, zPositive.H2O.thermal) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      annotation (
        defaultComponentPrefixes="replaceable",
        Documentation(info="<html><p>If C19HF37O5S is included (<code>inclC19HF37O5S</code>=<code>true</code>),
    then <code>initVelX</code>, <code>initVelY</code>, and <code>initVelZ</code> are
    ignored.  Velocity is not initialized because it is set by C19HF37O5S.</p>

    <p>For more information, see the
 <a href=\"modelica://FCSys.Subregions.Phases.BaseClasses.NullPhase\">NullPhase</a> model.</p></html>"),

        Diagram(graphics));
    end Ionomer;

    model Liquid "Liquid phase"

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
      replaceable Species.H2O.Liquid.Fixed H2O if inclH2O constrainedby
        Species.Species(
        final k=k,
        final inclLin=inclLin,
        initMethX=if initVelX and reduceVel then InitMethVelocity.None else
            InitMethVelocity.Velocity,
        initMethY=if initVelY and reduceVel then InitMethVelocity.None else
            InitMethVelocity.Velocity,
        initMethZ=if initVelZ and reduceVel then InitMethVelocity.None else
            InitMethVelocity.Velocity,
        initMethTemp=if initTemp and reduceTemp then InitMethScalar.None else
            InitMethScalar.Temperature,
        phi(each stateSelect=if reduceVel then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceTemp then StateSelect.default else StateSelect.prefer))
        "<html>H<sub>2</sub>O model</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclH2O),
        Placement(transformation(extent={{-10,-10},{10,10}})));

    equation
      // H2O
      // ---
      // Exchange
      connect(H2O.chemical, chemical.H2O) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(H2O.common.mechanical, common.mechanical) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2O.common.thermal, common.thermal) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2O.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Transport of x-axis linear momentum
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

      connect(H2O.yNegative.mechanicalX, yNegative.H2O.mechanicalX) annotation
        (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.yPositive.mechanicalX, yPositive.H2O.mechanicalX) annotation
        (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.mechanicalX, zNegative.H2O.mechanicalX) annotation
        (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.mechanicalX, zPositive.H2O.mechanicalX) annotation
        (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Transport of y-axis linear momentum
      connect(H2O.yNegative.material, yNegative.H2O.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.yPositive.material, yPositive.H2O.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.mechanicalY, zNegative.H2O.mechanicalY) annotation
        (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.mechanicalY, zPositive.H2O.mechanicalY) annotation
        (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.xNegative.mechanicalY, xNegative.H2O.mechanicalY) annotation
        (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.mechanicalY, xPositive.H2O.mechanicalY) annotation
        (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      // Transport of z-axis linear momentum
      connect(H2O.zNegative.material, zNegative.H2O.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.material, zPositive.H2O.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xNegative.mechanicalZ, xNegative.H2O.mechanicalZ) annotation
        (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.mechanicalZ, xPositive.H2O.mechanicalZ) annotation
        (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.mechanicalZ, yNegative.H2O.mechanicalZ) annotation
        (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.mechanicalZ, yPositive.H2O.mechanicalZ) annotation
        (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      // Thermal transport
      connect(H2O.xNegative.thermal, xNegative.H2O.thermal) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.xPositive.thermal, xPositive.H2O.thermal) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yNegative.thermal, yNegative.H2O.thermal) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.yPositive.thermal, yPositive.H2O.thermal) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.zNegative.thermal, zNegative.H2O.thermal) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.zPositive.thermal, zPositive.H2O.thermal) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      annotation (
        defaultComponentPrefixes="replaceable",
        Documentation(info="<html><p>For more information, see the
 <a href=\"modelica://FCSys.Subregions.Phases.BaseClasses.NullPhase\">NullPhase</a> model.</p></html>"),

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
          "<html>Area fill factor for transport (<b>k</b>)</html>"
          annotation (Dialog(group="Geometry"));
        outer parameter Q.Length L[Axis](each final min=Modelica.Constants.small)
          "Length" annotation (HideResult=true, missingInnerMessage=
              "This model should be used within the Subregion model.");
        outer parameter Q.Area A[Axis] "Cross-sectional area" annotation (
            HideResult=true, missingInnerMessage=
              "This model should be used within the Subregion model.");

        // Assumptions
        parameter Boolean reduceVel=true "Same velocity for all species"
          annotation (
          HideResult=true,
          Dialog(tab="Assumptions",enable=n_spec > 1),
          choices(__Dymola_checkBox=true));
        parameter Boolean reduceTemp=true "Same temperature for all species"
          annotation (
          HideResult=true,
          Dialog(tab="Assumptions",enable=n_spec > 1),
          choices(__Dymola_checkBox=true));

        // Initialization
        parameter Boolean initVelX=true "Initialize the x component"
          annotation (
          compact=true,
          Dialog(
            tab="Initialization",
            group="Velocity",
            enable=reduceVel),
          __Dymola_joinNext=true,
          choices(__Dymola_checkBox=true));
        parameter Boolean initVelY=true "Initialize the y component"
          annotation (
          compact=true,
          Dialog(
            tab="Initialization",
            group="Velocity",
            enable=reduceVel),
          __Dymola_joinNext=true,
          choices(__Dymola_checkBox=true));
        parameter Boolean initVelZ=true "Initialize the z component"
          annotation (
          compact=true,
          Dialog(
            tab="Initialization",
            group="Velocity",
            enable=reduceVel),
          __Dymola_joinNext=true,
          choices(__Dymola_checkBox=true));
        parameter Q.Velocity phi_IC[Axis]={0,0,0}
          "<html>Initial velocity (<b>&phi;</b><sub>IC</sub>)</html>"
          annotation (Dialog(tab="Initialization",group="Velocity"));
        // This is always enabled in the dialog since it's used as a guess value.
        parameter Boolean initTemp=true "Initialize" annotation (
          compact=true,
          Dialog(
            tab="Initialization",
            group="Temperature",
            enable=reduceTemp),
          choices(__Dymola_checkBox=true));

        parameter Q.TemperatureAbsolute T_IC(nominal=298.15*U.K,start=defaults.T)
          "<html>Initial temperature (<i>T</i><sub>IC</sub>)</html>"
          annotation (Dialog(tab="Initialization",group="Temperature"));
        // This is always enabled in the dialog since it's used as a guess value.

        parameter Boolean inclLin[3]={true,false,false}
          "true, if each component of linear momentum is included"
          annotation (Evaluate=true, Dialog(tab="Assumptions"));

        FCSys.Connectors.InertAmagat inert(final n_lin=n_lin) if n_spec > 0
          annotation (Placement(transformation(extent={{10,-30},{30,-10}}),
              iconTransformation(extent={{70,-90},{90,-70}})));
        FCSys.Connectors.FaceBus xNegative if n_spec > 0
          "Negative face along the x axis" annotation (Placement(transformation(
                extent={{-50,-10},{-30,10}}),iconTransformation(extent={{-90,-10},
                  {-70,10}})));
        FCSys.Connectors.FaceBus xPositive if n_spec > 0
          "Positive face along the x axis" annotation (Placement(transformation(
                extent={{30,-10},{50,10}}),iconTransformation(extent={{70,-10},
                  {90,10}})));
        FCSys.Connectors.FaceBus yNegative if n_spec > 0
          "Negative face along the y axis" annotation (Placement(transformation(
                extent={{-10,-50},{10,-30}}),iconTransformation(extent={{-10,-94},
                  {10,-74}})));
        FCSys.Connectors.FaceBus yPositive if n_spec > 0
          "Positive face along the y axis" annotation (Placement(transformation(
                extent={{-10,30},{10,50}}),iconTransformation(extent={{-10,90},
                  {10,110}})));
        FCSys.Connectors.FaceBus zNegative if n_spec > 0
          "Negative face along the z axis" annotation (Placement(transformation(
                extent={{10,10},{30,30}}),iconTransformation(extent={{40,40},{
                  60,60}})));
        FCSys.Connectors.FaceBus zPositive if n_spec > 0
          "Positive face along the z axis" annotation (Placement(transformation(
                extent={{-30,-30},{-10,-10}}),iconTransformation(extent={{-90,-90},
                  {-70,-70}})));

        PhaseBoundary phaseBoundary(final n_lin=n_lin) if n_spec > 0
          "Phase boundary" annotation (Placement(transformation(
              extent={{-18,-18},{18,18}},
              rotation=0,
              origin={0,0})));
        // This component is conditional because if two or more empty phases
        // (without any species included) were connected within a subregion, there
        // would be a mathematical singularity.

        FCSys.Connectors.ChemicalBus chemical if n_spec > 0
          "Connector for chemical exchange among species" annotation (Placement(
              transformation(extent={{-30,10},{-10,30}}),iconTransformation(
                extent={{-60,40},{-40,60}})));

      protected
        parameter Integer n_spec "Number of species";
        final parameter Integer n_lin=countTrue(inclLin)
          "Number of components of linear momentum"
          annotation (Evaluate=true,HideResult=true);
        final parameter Integer cartAxes[n_lin]=index(inclLin)
          "Cartesian-axis indices of the axes of linear momentum";

        FCSys.Connectors.InertInternal common(
          n_lin=n_lin,
          final uniformVelocity=reduceVel,
          final isothermal=reduceTemp,
          mechanical(phi(
              each stateSelect=StateSelect.prefer,
              final start=phi_IC[cartAxes],
              fixed={initVelX,initVelY,initVelZ}[cartAxes])),
          thermal(T(
              stateSelect=StateSelect.prefer,
              final start=T_IC,
              final fixed=initTemp))) if n_spec > 0
          "Internal connector to directly couple velocities and/or temperatures"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={0,0}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={0,0})));

        outer BCs.Defaults defaults "Default settings" annotation (Placement(
              transformation(extent={{40,40},{60,60}}),iconTransformation(
                extent={{-10,90},{10,110}})));

      equation
        // Inert interactions
        connect(phaseBoundary.inertA, inert) annotation (Line(
            points={{12,-12},{20,-20}},
            color={72,90,180},
            smooth=Smooth.None));

        annotation (
          Documentation(info="<html><p>If one of the species has <code>setTemp = true</code>, then
    <code>initTemp</code> should be set to <code>false</code>.
    Likewise, if one of the species has <code>setVelX = true</code>,
    <code>setVelY = true</code>, or <code>setVelZ = true</code>, then
    <code>initVelX</code>, <code>initVelY</code>, or <code>initVelZ</code> should
    be set to <code>false</code> (respectively).</p>

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
    "Models for single-species storage, transport, and exchange of material, linear momentum, and energy"
    extends Modelica.Icons.Package;
    package C "C"
      extends Modelica.Icons.Package;
      package Graphite "C graphite"
        extends Modelica.Icons.Package;
        model Calibrated "Correlations with adjustment factors"
          extends SpeciesInertStagnant(redeclare replaceable package Data =
                FCSys.Characteristics.C.Graphite, R=k_R*Data.R(T));

          parameter Q.NumberAbsolute k_R(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub><i>R</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C",
            Documentation(info="<html><p>For more information, see the
    <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"),

            Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}), graphics),
            Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}), graphics={Text(
                          extent={{-150,90},{-118,52}},
                          lineColor={0,0,255},
                          textString="%t.test")}));

        end Calibrated;

        model Correlated "Correlated properties"
          extends SpeciesInertStagnant(redeclare replaceable package Data =
                FCSys.Characteristics.C.Graphite, R=Data.R(T));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C",
            Documentation(info=
                  "<html><p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"),

            Icon(graphics));

        end Correlated;

        model Fixed "Fixed properties"
          extends SpeciesInertStagnant(redeclare replaceable package Data =
                FCSys.Characteristics.C.Graphite (
                Deltah0_f=0,
                Deltah0=0,
                specHeatCapPow=0,
                T_lim_c={0,Modelica.Constants.inf},
                b_c=[935*U.J*Data.m/(U.kg*U.K)],
                B_c=[-298.15*U.K*935*U.J*Data.m/(U.kg*U.K) + Data.Deltah0_f,
                    Polynomial.F(
                            298.15*U.K,
                            FCSys.Characteristics.C.Graphite.b_c[1, :],
                            -3) + FCSys.Characteristics.C.Graphite.B_c[1, 2] -
                    Data.b_c[1, 1]*ln(298.15*U.K)]), redeclare parameter
              Q.Resistivity R=U.m*U.K/(5.70*U.W));
          // See the documentation for a table of values.
          // Note:  Parameter expressions (e.g., involving defaults.T) are not used
          // here since they would render the parameters unadjustable in Dymola 7.4.
          // A similar note applies to the other species.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Fixed specific heat capacity (independent of temperature)
    <li>Thermal resistivity <i>r</i><sub>th</sub> is fixed (e.g., independent of temperature)</li>
    </ol></p>

    <p>The default isobaric specific heat capacity (<code>b_c=[0, 935*U.J*Data.m/(U.kg*U.K)]</code>) and thermal
   resistivity (<code>alpha_Qdot=U.m*U.K/(11.1*U.W)</code>) is based on data of graphite fiber epoxy (25% vol)<br>composite at 300 K from
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
      <th rowspan=2 valign=\"middle\"><code>c_p*U.kg<br>*U.K<br>/(U.J<br>*Data.m)</code></th>
      <th rowspan=2 valign=\"middle\"><code>R<br>*U.W<br>/(U.m<br>*U.K)</code></th>
      <th rowspan=2 valign=\"middle\"><code>R<br>*U.W<br>/(U.m<br>*U.K)</code></th>
      <th rowspan=2 valign=\"middle\"><code>c_p*U.kg<br>*U.W<br>/(U.J<br>*Data.m)</code></th>
      <th rowspan=1 colspan=2><code>R*U.W/(U.m*U.K)</code></th>
      <th rowspan=2 valign=\"middle\"><code>c_p*U.kg<br>*U.K<br>/(U.J<br>*Data.m)</code></th>
      <th rowspan=1 colspan=2><code>R*U.W/(U.m*U.K)</code></th>
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

  <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));
        end Fixed;
      end Graphite;
    end C;

    package C19HF37O5S
      "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S</html>"
      extends Modelica.Icons.Package;
      package Solid
        "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S solid</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Correlations with adjustment factors"
          extends SpeciesInertStagnant(redeclare replaceable package Data =
                FCSys.Characteristics.C19HF37O5S.Solid, R=k_R*Data.R(T));

          parameter Q.NumberAbsolute k_R(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub><i>R</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C19HF37O5S",
            Documentation(info=
                  "<html><p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"),

            Diagram(graphics));

        end Calibrated;

        model Correlated "Correlated properties"
          extends SpeciesInertStagnant(redeclare replaceable package Data =
                FCSys.Characteristics.C19HF37O5S.Solid, R=Data.R(T));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C19HF37O5S",
            Documentation(info=
                  "<html><p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends SpeciesInertStagnant(redeclare replaceable package Data =
                FCSys.Characteristics.C19HF37O5S.Solid, redeclare parameter
              Q.ResistivityThermal R=Data.R());

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C19HF37O5S",
            Documentation(info="<html><p>Assumptions:
    <ol>
    <li>Thermal resistivity <i>r</i><sub>th</sub> is fixed (e.g., independent of temperature)</li>
    </ol></p>

    <p>Notes:
    <ul><li>The default thermal resistivity (<code>R=U.m*U.K/(0.16*U.W)</code>) is of dry
  Nafion 115 [<a href=\"modelica://FCSys.UsersGuide.References\">Kandlikar2009</a>, p. 1277].</li>
  </ul>
  </p><p>For more information, see the
    <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));
        end Fixed;
      end Solid;
    end C19HF37O5S;

    package 'e-' "<html>e<sup>-</sup></html>"
      extends Modelica.Icons.Package;
      package Graphite "<html>e<sup>-</sup> in graphite</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Correlations with adjustment factors"
          extends Species0Amount(
            redeclare replaceable package Data =
                FCSys.Characteristics.'e-'.Graphite,
            Xi=k_Xi*Data.Xi(T),
            F=k_F*Data.F(T),
            R=k_R*Data.R(T));

          parameter Q.NumberAbsolute k_Xi(final nominal=1) = 1
            "<html>Adjustment factor for self diffusivity (<i>k</i><sub>&Xi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_F(final nominal=1) = 1
            "<html>Adjustment factor for fluidity (<i>k</i><sub><i>F</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_R(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&alpha; <i>S&#775;</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'e-'",
            Documentation(info=
                  "<html><p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));
        end Calibrated;

        model Correlated "Correlated properties"
          extends Species0Amount(
            redeclare replaceable package Data =
                FCSys.Characteristics.'e-'.Graphite,
            Xi=Data.Xi(T),
            F=Data.F(T),
            R=Data.R(T));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'e-'",
            Documentation(info=
                  "<html><p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends Species0Amount(
            redeclare replaceable package Data =
                FCSys.Characteristics.'e-'.Graphite,
            redeclare parameter Q.Fusivity Xi=Data.Xi(),
            redeclare parameter Q.Fluidity F=Data.F(),
            redeclare parameter Q.ResistivityThermal R=Data.R());

          annotation (
            group="Material properties",
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'e-'",
            Documentation(info="<html>
    <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"),

            Diagram(graphics));

        end Fixed;
      end Graphite;
    end 'e-';

    package 'H+' "<html>H<sup>+</sup></html>"
      extends Modelica.Icons.Package;
      package Solid "<html>H<sup>+</sup> in solid</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Correlations with adjustment factors"
          extends Species0Amount(
            redeclare replaceable package Data =
                FCSys.Characteristics.'H+'.Solid,
            initMethPartNum=InitMethScalar.Amount,
            Xi=k_Xi*Data.Xi(T),
            F=k_F*Data.F(T),
            R=k_R*Data.R(T));

          // Note:  initMethPartNum may not be Pressure (which is default) since the
          // EOS doesn't involve pressure.

          parameter Q.NumberAbsolute k_Xi(final nominal=1) = 1
            "<html>Adjustment factor for self diffusivity (<i>k</i><sub>&Xi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_F(final nominal=1) = 1
            "<html>Adjustment factor for fluidity (<i>k</i><sub><i>F</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_R(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&alpha; <i>S&#775;</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'H+'",
            Documentation(info=
                  "<html><p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));
        end Calibrated;

        model Correlated "Correlated properties"
          extends Species0Amount(
            redeclare replaceable package Data =
                FCSys.Characteristics.'H+'.Solid,
            initMethPartNum=InitMethScalar.Amount,
            Xi=Data.Xi(T),
            F=Data.F(T),
            R=Data.R(T));

          // Note:  initMethPartNum may not be Pressure (which is default) since the
          // EOS doesn't involve pressure.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'H+'",
            Documentation(info=
                  "<html><p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends Species0Amount(
            redeclare replaceable package Data =
                FCSys.Characteristics.'H+'.Solid,
            initMethPartNum=InitMethScalar.Amount,
            redeclare parameter Q.Fusivity Xi=Data.Xi(),
            redeclare parameter Q.Fluidity F=Data.F(),
            redeclare parameter Q.ResistivityThermal R=U.m*U.K/(0.1661*U.W));

          // Note:  initMethPartNum may not be Pressure (which is default) since
          // overrideEOS is true.

          // See the documentation for a table of values.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'H+'",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Resistivities (<i>f</i><sub>0</sub>, <i>r</i><sub>th</sub>, etc.) are fixed (e.g., independent of temperature)</li>
    </ol></p>

  <p>The default resistivities (<code>F=1/(5.3e-6*U.Pa*U.s)</code> and <code>R=U.m*U.K/(0.1661*U.W)</code>) are of H gas
  (rather than H<sup>+</sup>) at 300 K from [<a href=\"modelica://FCSys.UsersGuide.References\">Schetz1996</a>, p. 139].
  Table 1 lists the properties at other temperatures.</p>

    <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"bottom\"><b>Table 1:</b> Properties of H gas (not H<sup>+</sup>) (<a href=\"modelica://FCSys.UsersGuide.References\">Schetz and Fuhs, 1996</a>, p. 139)</caption>
<tr>
      <th valign=\"middle\"><code>T/U.K</code></th>
      <th width=1><code>F<br>*U.Pa*U.s</code></th>
      <th width=1><code>R*U.W<br>/(U.m*U.K)</code></th>
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

</p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));
        end Fixed;
      end Solid;
    end 'H+';

    package H2 "<html>H<sub>2</sub></html>"
      extends Modelica.Icons.Package;
      package Gas "<html>H<sub>2</sub> gas</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Correlations with adjustment factors"
          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.H2.Gas (
                  b_v=[1], specVolPow={-1,0}),
            Xi=k_Xi*Data.Xi(T),
            F=k_F*Data.F(T),
            R=k_R*Data.R(T));

          parameter Q.NumberAbsolute k_Xi(final nominal=1) = 1
            "<html>Adjustment factor for self diffusivity (<i>k</i><sub>&Xi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_F(final nominal=1) = 1
            "<html>Adjustment factor for fluidity (<i>k</i><sub><i>F</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_R(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub><i>R</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));
        end Calibrated;

        model Correlated "Correlated properties"
          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.H2.Gas (
                  b_v=[1], specVolPow={-1,0}),
            Xi=Data.Xi(T),
            F=Data.F(T),
            R=Data.R(T));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.H2.Gas (
                  b_v=[1], specVolPow={-1,0}),
            redeclare parameter Q.Fusivity Xi=Data.Xi(),
            redeclare parameter Q.Fluidity F=1/(89.6e-7*U.Pa*U.s),
            redeclare parameter Q.ResistivityThermal R=U.m*U.K/(183e-3*U.W));

          // See the documentation for a table of values.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
    <li>Resistivities (<i>f</i><sub>0</sub>, <i>r</i><sub>th</sub>, etc.) are fixed (e.g., independent of temperature)</li>
    </ol></p>

<p>Additional notes:<ul>
<li>
  The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

<p>The default resistivities (<code>F=1/(89.6e-7*U.Pa*U.s)</code>
and <code>R=U.m*U.K/(183e-3*U.W)</code>) are based on data of H<sub>2</sub> gas at 1 atm and
  300 K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 919&ndash;920].
  Table 1 lists the properties at  other temperatures. </p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"bottom\"><b>Table 1:</b> Properties of H<sub>2</sub> gas at 1 atm (<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 919&ndash;920).</caption>
  <tr>
      <th rowspan=2 valign=\"middle\"><code>T/U.K</code></th>
      <th rowspan=2 width=1><code>c_p*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th rowspan=2 width=1 ><code>F<br>*U.Pa*U.s</code></th>
      <th rowspan=2 width=1 ><code>R*U.W<br>/(U.m*U.K)</code></th>
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
<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"),

            Diagram(graphics));

        end Fixed;
      end Gas;
    end H2;

    package H2O "<html>H<sub>2</sub>O</html>"
      extends Modelica.Icons.Package;
      package Gas "<html>H<sub>2</sub>O gas</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Correlations with adjustment factors"
          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.H2O.Gas
                (b_v=[1], specVolPow={-1,0}),
            Xi=k_Xi*Data.Xi(T),
            F=k_F*Data.F(T),
            R=k_R*Data.R(T));

          parameter Q.NumberAbsolute k_Xi(final nominal=1) = 1
            "<html>Adjustment factor for self diffusivity (<i>k</i><sub>&Xi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_F(final nominal=1) = 1
            "<html>Adjustment factor for fluidity (<i>k</i><sub><i>F</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_R(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub><i>R</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));
        end Calibrated;

        model Correlated "Correlated properties"
          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.H2O.Gas
                (b_v=[1], specVolPow={-1,0}),
            Xi=Data.Xi(T),
            F=Data.F(T),
            R=Data.R(T));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.H2O.Gas
                (b_v=[1], specVolPow={-1,0}),
            redeclare parameter Q.Fusivity Xi=Data.Xi(),
            redeclare parameter Q.Fluidity F=1/(9.09e-6*U.Pa*U.s),
            redeclare parameter Q.ResistivityThermal R=U.m*U.K/(19.6e-3*U.W));

          // See the documentation for tables of values.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
    <li>Resistivities (<i>f</i><sub>0</sub>, <i>r</i><sub>th</sub>, etc.) are fixed (e.g., independent of temperature)</li>
    </ol></p>

          <p>Notes:<ul>
<ul><li>
  The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

<p>The default resistivities (<code>F=1/(9.09e-6*U.Pa*U.s)</code>
and <code>R=U.m*U.K/(19.6e-3*U.W)</code>) are of H<sub>2</sub>O gas at saturation pressure and
  300 K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 921].  Table 1 lists the properties at
  saturation pressure and other temperatures.  Table 2 lists the properties of H<sub>2</sub>O gas at 1 atm.
  See also
  <a href=\"http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html\">http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html</a>.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"bottom\"><b>Table 1:</b> Properties of H<sub>2</sub>O at saturation pressure (<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 924&ndash;925).</caption>
  <tr>
      <th valign=\"middle\"><code>T/U.K</code></th>
      <th width=1><code>c_p*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1><code>F<br>*U.Pa*U.s</code></th>
      <th width=1><code>R*U.W<br>/(U.m*U.K)</code></th>
    </tr>
<tr><td>273.15</td><td>1854</td><td>1/8.02e-6</td><td>1/18.2e-3</td></tr>
<tr><td>275</td><td>1855</td><td>1/8.09e-6</td><td>1/18.3e-3</td></tr>
<tr><td>280</td><td>1858</td><td>1/8.29e-6</td><td>1/18.6e-3</td></tr>
<tr><td>285</td><td>1861</td><td>1/8.49e-6</td><td>1/18.9e-3</td></tr>
<tr><td>290</td><td>1864</td><td>1/8.69e-6</td><td>1/19.3e-3</td></tr>
<tr><td>295</td><td>1868</td><td>1/8.89e-6</td><td>1/19.5e-3</td></tr>
<tr><td>300</td><td>1872</td><td>1/9.09e-6</td><td>1/19.6e-3</td></tr>
<tr><td>305</td><td>1877</td><td>1/9.29e-6</td><td>1/20.1e-3</td></tr>
<tr><td>310</td><td>1882</td><td>1/9.49e-6</td><td>1/20.4e-3</td></tr>
<tr><td>315</td><td>1888</td><td>1/9.69e-6</td><td>1/20.7e-3</td></tr>
<tr><td>320</td><td>1895</td><td>1/9.89e-6</td><td>1/21.0e-3</td></tr>
<tr><td>325</td><td>1903</td><td>1/10.09e-6</td><td>1/21.3e-3</td></tr>
<tr><td>330</td><td>1911</td><td>1/10.29e-6</td><td>1/21.7e-3</td></tr>
<tr><td>335</td><td>1920</td><td>1/10.49e-6</td><td>1/22.0e-3</td></tr>
<tr><td>340</td><td>1930</td><td>1/10.69e-6</td><td>1/22.3e-3</td></tr>
<tr><td>345</td><td>1941</td><td>1/10.89e-6</td><td>1/22.6e-3</td></tr>
<tr><td>350</td><td>1954</td><td>1/11.09e-6</td><td>1/23.0e-3</td></tr>
<tr><td>355</td><td>1968</td><td>1/11.29e-6</td><td>1/23.3e-3</td></tr>
<tr><td>360</td><td>1983</td><td>1/11.49e-6</td><td>1/23.7e-3</td></tr>
<tr><td>365</td><td>1999</td><td>1/11.69e-6</td><td>1/24.1e-3</td></tr>
<tr><td>370</td><td>2017</td><td>1/11.89e-6</td><td>1/24.5e-3</td></tr>
<tr><td>373.15</td><td>2029</td><td>1/12.02e-6</td><td>1/24.8e-3</td></tr>
<tr><td>375</td><td>2036</td><td>1/12.09e-6</td><td>1/24.9e-3</td></tr>
<tr><td>380</td><td>2057</td><td>1/12.29e-6</td><td>1/25.4e-3</td></tr>
<tr><td>385</td><td>2080</td><td>1/12.49e-6</td><td>1/25.8e-3</td></tr>
<tr><td>390</td><td>2104</td><td>1/12.69e-6</td><td>1/26.3e-3</td></tr>
<tr><td>400</td><td>2158</td><td>1/13.05e-6</td><td>1/27.2e-3</td></tr>
<tr><td>410</td><td>2221</td><td>1/13.42e-6</td><td>1/28.2e-3</td></tr>
<tr><td>420</td><td>2291</td><td>1/13.79e-6</td><td>1/29.8e-3</td></tr>
<tr><td>430</td><td>2369</td><td>1/14.14e-6</td><td>1/30.4e-3</td></tr>
<tr><td>440</td><td>2460</td><td>1/14.50e-6</td><td>1/3.17e-3</td></tr>
<tr><td>450</td><td>2560</td><td>1/14.85e-6</td><td>1/33.1e-3</td></tr>
<tr><td>460</td><td>2680</td><td>1/15.19e-6</td><td>1/34.6e-3</td></tr>
<tr><td>470</td><td>2790</td><td>1/15.54e-6</td><td>1/36.3e-3</td></tr>
<tr><td>480</td><td>2940</td><td>1/15.88e-6</td><td>1/38.1e-3</td></tr>
<tr><td>490</td><td>3100</td><td>1/16.23e-6</td><td>1/40.1e-3</td></tr>
<tr><td>500</td><td>3270</td><td>1/16.59e-6</td><td>1/42.3e-3</td></tr>
<tr><td>510</td><td>3470</td><td>1/16.95e-6</td><td>1/44.7e-3</td></tr>
<tr><td>520</td><td>3700</td><td>1/17.33e-6</td><td>1/47.5e-3</td></tr>
<tr><td>530</td><td>3960</td><td>1/17.72e-6</td><td>1/50.6e-3</td></tr>
<tr><td>540</td><td>4270</td><td>1/18.1e-6</td><td>1/54.0e-3</td></tr>
<tr><td>550</td><td>4640</td><td>1/18.6e-6</td><td>1/58.3e-3</td></tr>
<tr><td>560</td><td>5090</td><td>1/19.1e-6</td><td>1/63.7e-3</td></tr>
<tr><td>570</td><td>5670</td><td>1/19.7e-6</td><td>1/76.7e-3</td></tr>
<tr><td>580</td><td>6400</td><td>1/20.4e-6</td><td>1/76.7e-3</td></tr>
<tr><td>590</td><td>7350</td><td>1/21.5e-6</td><td>1/84.1e-3</td></tr>
<tr><td>600</td><td>8750</td><td>1/22.7e-6</td><td>1/92.9e-3</td></tr>
<tr><td>610</td><td>11100</td><td>1/24.1e-6</td><td>1/103e-3</td></tr>
<tr><td>620</td><td>15400</td><td>1/25.9e-6</td><td>1/114e-3</td></tr>
<tr><td>635</td><td>18300</td><td>1/27.0e-6</td><td>1/121e-3</td></tr>
<tr><td>630</td><td>22100</td><td>1/28.0e-6</td><td>1/130e-3</td></tr>
<tr><td>635</td><td>27600</td><td>1/30.0e-6</td><td>1/141e-3</td></tr>
<tr><td>640</td><td>42000</td><td>1/32.0e-6</td><td>1/155e-3</td></tr>
  </table>

<br>

    <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"bottom\"><b>Table 2:</b> Properties of H<sub>2</sub>O gas at 1 atm (<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 921).</caption>
  <tr>
      <th valign=\"middle\"><code>T/U.K</code></th>
      <th width=1><code>c_p*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1 ><code>F<br>*U.Pa*U.s</code></th>
      <th width=1 ><code>R*U.W<br>/(U.m*U.K)</code></th>
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

  </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));
        end Fixed;
      end Gas;

      package Liquid "<html>H<sub>2</sub>O liquid</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Correlations with adjustment factors"
          extends Species(
            redeclare replaceable package Data =
                FCSys.Characteristics.H2O.Liquid,
            Xi=k_Xi*Data.Xi(T),
            F=k_F*Data.F(T),
            R=k_R*Data.R(T));

          parameter Q.NumberAbsolute k_Xi(final nominal=1) = 1
            "<html>Adjustment factor for self diffusivity (<i>k</i><sub>&Xi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_F(final nominal=1) = 1
            "<html>Adjustment factor for fluidity (<i>k</i><sub><i>F</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_R(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub><i>R</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html>
    <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));
        end Calibrated;

        model Correlated "Correlated properties"
          extends Species(
            redeclare replaceable package Data =
                FCSys.Characteristics.H2O.Liquid,
            Xi=Data.Xi(T),
            F=Data.F(T),
            R=Data.R(T));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html>
         <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare replaceable package Data =
                FCSys.Characteristics.H2O.Liquid,
            redeclare parameter Q.Fusivity Xi=Data.Xi(),
            redeclare parameter Q.Fluidity F=1/(855e-6*U.Pa*U.s),
            redeclare parameter Q.ResistivityThermal R=U.m*U.K/(613e-3*U.W));

          // See the documentation for tables of values.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Resistivities (<i>f</i><sub>0</sub>, <i>r</i><sub>th</sub>, etc.) are fixed (e.g., independent of temperature)</li>
    </ol></p>

          <p>Notes:<ul>
  <li>
  The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

<p>The default resistivities (<code>F=1/(855e-6*U.Pa*U.s)</code>
and <code>R=U.m*U.K/(613e-3*U.W)</code>) are of H<sub>2</sub>O liquid at saturation pressure and
  300 K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 921].  Table 1 lists the properties at
  saturation pressure and other temperatures.
  See also
  <a href=\"http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html\">http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html</a>.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"bottom\"><b>Table 1:</b> Properties of H<sub>2</sub>O at saturation pressure (<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 924&ndash;925).</caption>
  <tr>
      <th valign=\"middle\"><code>T/U.K</code></th>
      <th width=1><code>c_p*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1><code>F<br>*U.Pa*U.s</code></th>
      <th width=1><code>R*U.W<br>/(U.m*U.K)</code></th>
    </tr>
<tr><td>273.15</td><td>4217</td><td>1/1750e-6</td><td>1/569e-3</td></tr>
<tr><td>275</td><td>4211</td><td>1/1652e-6</td><td>1/574e-3</td></tr>
<tr><td>280</td><td>4198</td><td>1/1422e-6</td><td>1/582e-3</td></tr>
<tr><td>285</td><td>4189</td><td>1/1225e-6</td><td>1/590e-3</td></tr>
<tr><td>290</td><td>4184</td><td>1/1080e-6</td><td>1/598e-3</td></tr>
<tr><td>295</td><td>4181</td><td>1/959e-6</td><td>1/606e-3</td></tr>
<tr><td>300</td><td>4179</td><td>1/855e-6</td><td>1/613e-3</td></tr>
<tr><td>305</td><td>4178</td><td>1/769e-6</td><td>1/620e-3</td></tr>
<tr><td>310</td><td>4178</td><td>1/695e-6</td><td>1/628e-3</td></tr>
<tr><td>315</td><td>4179</td><td>1/631e-6</td><td>1/634e-3</td></tr>
<tr><td>320</td><td>4180</td><td>1/577e-6</td><td>1/640e-3</td></tr>
<tr><td>325</td><td>4182</td><td>1/528e-6</td><td>1/645e-3</td></tr>
<tr><td>330</td><td>4184</td><td>1/489e-6</td><td>1/650e-3</td></tr>
<tr><td>335</td><td>4186</td><td>1/453e-6</td><td>1/656e-3</td></tr>
<tr><td>340</td><td>4188</td><td>1/420e-6</td><td>1/660e-3</td></tr>
<tr><td>345</td><td>4191</td><td>1/389e-6</td><td>1/668e-3</td></tr>
<tr><td>350</td><td>4195</td><td>1/365e-6</td><td>1/668e-3</td></tr>
<tr><td>355</td><td>4199</td><td>1/343e-6</td><td>1/671e-3</td></tr>
<tr><td>360</td><td>4203</td><td>1/324e-6</td><td>1/674e-3</td></tr>
<tr><td>365</td><td>4209</td><td>1/306e-6</td><td>1/677e-3</td></tr>
<tr><td>370</td><td>4214</td><td>1/289e-6</td><td>1/679e-3</td></tr>
<tr><td>373.15</td><td>4217</td><td>1/279e-6</td><td>1/680e-3</td></tr>
<tr><td>375</td><td>4220</td><td>1/274e-6</td><td>1/681e-3</td></tr>
<tr><td>380</td><td>4226</td><td>1/260e-6</td><td>1/683e-3</td></tr>
<tr><td>385</td><td>4232</td><td>1/248e-6</td><td>1/685e-3</td></tr>
<tr><td>390</td><td>4239</td><td>1/237e-6</td><td>1/686e-3</td></tr>
<tr><td>400</td><td>4256</td><td>1/217e-6</td><td>1/688e-3</td></tr>
<tr><td>410</td><td>4278</td><td>1/200e-6</td><td>1/688e-3</td></tr>
<tr><td>420</td><td>4302</td><td>1/185e-6</td><td>1/688e-3</td></tr>
<tr><td>430</td><td>4331</td><td>1/173e-6</td><td>1/685e-3</td></tr>
<tr><td>440</td><td>4360</td><td>1/162e-6</td><td>1/682e-3</td></tr>
<tr><td>450</td><td>4400</td><td>1/152e-6</td><td>1/678e-3</td></tr>
<tr><td>460</td><td>4440</td><td>1/143e-6</td><td>1/673e-3</td></tr>
<tr><td>470</td><td>4480</td><td>1/136e-6</td><td>1/667e-3</td></tr>
<tr><td>480</td><td>4530</td><td>1/129e-6</td><td>1/660e-3</td></tr>
<tr><td>490</td><td>4590</td><td>1/124e-6</td><td>1/651e-3</td></tr>
<tr><td>500</td><td>4660</td><td>1/118e-6</td><td>1/642e-3</td></tr>
<tr><td>510</td><td>4740</td><td>1/113e-6</td><td>1/631e-3</td></tr>
<tr><td>520</td><td>4840</td><td>1/108e-6</td><td>1/621e-3</td></tr>
<tr><td>530</td><td>4950</td><td>1/104e-6</td><td>1/608e-3</td></tr>
<tr><td>540</td><td>5080</td><td>1/101e-6</td><td>1/594e-3</td></tr>
<tr><td>550</td><td>5240</td><td>1/97e-6</td><td>1/580e-3</td></tr>
<tr><td>560</td><td>5430</td><td>1/94e-6</td><td>1/563e-3</td></tr>
<tr><td>570</td><td>5680</td><td>1/91e-6</td><td>1/548e-3</td></tr>
<tr><td>580</td><td>6000</td><td>1/88e-6</td><td>1/528e-3</td></tr>
<tr><td>590</td><td>6410</td><td>1/84e-6</td><td>1/513e-3</td></tr>
<tr><td>600</td><td>7000</td><td>1/81e-6</td><td>1/497e-3</td></tr>
<tr><td>610</td><td>7850</td><td>1/77e-6</td><td>1/467e-3</td></tr>
<tr><td>620</td><td>9350</td><td>1/72e-6</td><td>1/444e-3</td></tr>
<tr><td>635</td><td>10600</td><td>1/70e-6</td><td>1/430e-3</td></tr>
<tr><td>630</td><td>12600</td><td>1/67e-6</td><td>1/412e-3</td></tr>
<tr><td>635</td><td>16400</td><td>1/64e-6</td><td>1/392e-3</td></tr>
<tr><td>640</td><td>26000</td><td>1/59e-6</td><td>1/367e-3</td></tr>
  </table>
<br>

  <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));
        end Fixed;
      end Liquid;
    end H2O;

    package N2 "<html>N<sub>2</sub></html>"
      extends Modelica.Icons.Package;
      package Gas "<html>N<sub>2</sub> gas</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Correlations with adjustment factors"
          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.N2.Gas (
                  b_v=[1], specVolPow={-1,0}),
            Xi=k_Xi*Data.Xi(T),
            F=k_F*Data.F(T),
            R=k_R*Data.R(T));

          parameter Q.NumberAbsolute k_Xi(final nominal=1) = 1
            "<html>Adjustment factor for self diffusivity (<i>k</i><sub>&Xi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_F(final nominal=1) = 1
            "<html>Adjustment factor for fluidity (<i>k</i><sub><i>F</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_R(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub><i>R</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="N2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));
        end Calibrated;

        model Correlated "Correlated properties"
          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.N2.Gas (
                  b_v=[1], specVolPow={-1,0}),
            Xi=Data.Xi(T),
            F=Data.F(T),
            R=Data.R(T));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="N2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.N2.Gas (
                b_v=[1],
                specVolPow={-1,0},
                specHeatCapPow=0,
                T_lim_c={0,Modelica.Constants.inf},
                b_c=[1.041e3*U.J*Data.m/(U.kg*U.K)],
                B_c=[-298.15*U.K*1.041e3*U.J*Data.m/(U.kg*U.K) + Data.Deltah0_f,
                    Polynomial.F(
                            298.15*U.K,
                            FCSys.Characteristics.N2.Gas.b_c[1, :],
                            -3) + FCSys.Characteristics.N2.Gas.B_c[1, 2] - Data.b_c[
                    1, 1]*ln(298.15*U.K)]),
            redeclare parameter Q.Fusivity Xi=Data.Xi(),
            redeclare parameter Q.Fluidity F=1/(178.2e-7*U.Pa*U.s),
            redeclare parameter Q.ResistivityThermal R=U.m*U.K/(25.9e-3*U.W));

          // See the documentation for a table of values.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="N2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
    <li>Fixed specific heat capacity (independent of temperature)</ol>
    <li>Resistivities (<i>f</i><sub>0</sub>, <i>r</i><sub>th</sub>, etc.) are fixed (e.g., independent of temperature)</li>
    </ol></p>

<p>The default specific heat capacity (<code>b_c=[1.041e3*U.J*Data.m/(U.kg*U.K)]</code>) and resistivities (<code>F=1/(178.2e-7*U.Pa*U.s)</code> and <code>R=U.m*U.K/(25.9e-3*U.W))</code>) are based on data of gas at 1 atm and
  300 K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 920].
  Table 1 lists the properties at  other temperatures. Note that the value for isobaric specific heat capacity at
  800 K (<code>c_p=1.22e3*U.J*Data.m/(U.kg*U.K)</code>) seems unusual, but it matches the
  reference.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"bottom\"><b>Table 1:</b> Properties of N<sub>2</sub> gas at 1 atm (<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 920)</caption>
  <tr>
      <th valign=\"middle\"><code>T/U.K</code></th>
      <th width=1><code>c_p*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1><code>F<br>*U.Pa*U.s</code></th>
      <th width=1><code>R*U.W<br>/(U.m*U.K)</code></th>
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

  <p>The fluidity of air at 15.0 &deg;C and 1 atm is given by
       <code>F=1/(178e-7*U.Pa*U.s)</code>
   (<a href=\"http://en.wikipedia.org/wiki/Viscosity\">http://en.wikipedia.org/wiki/Viscosity</a>).</p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));
        end Fixed;
      end Gas;
    end N2;

    package O2 "<html>O<sub>2</sub></html>"
      extends Modelica.Icons.Package;
      package Gas "<html>O<sub>2</sub> gas</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Correlations with adjustment factors"
          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.O2.Gas (
                  b_v=[1], specVolPow={-1,0}),
            Xi=k_Xi*Data.Xi(T),
            F=k_F*Data.F(T),
            R=k_R*Data.R(T));

          parameter Q.NumberAbsolute k_Xi(final nominal=1) = 1
            "<html>Adjustment factor for self diffusivity (<i>k</i><sub>&Xi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_F(final nominal=1) = 1
            "<html>Adjustment factor for fluidity (<i>k</i><sub><i>F</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_R(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub><i>R</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="O2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));
        end Calibrated;

        model Correlated "Correlated properties"
          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.O2.Gas (
                  b_v=[1], specVolPow={-1,0}),
            Xi=Data.Xi(T),
            F=Data.F(T),
            R=Data.R(T));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="O2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.O2.Gas (
                  b_v=[1], specVolPow={-1,0}),
            redeclare parameter Q.Fusivity Xi=Data.Xi(),
            redeclare parameter Q.Fluidity F=1/(207.2e-7*U.Pa*U.s),
            redeclare parameter Q.ResistivityThermal R=U.m*U.K/(26.8e-3*U.W));

          // See the documentation for a table of values.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="O2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
    <li>Resistivities (<i>f</i><sub>0</sub>, <i>r</i><sub>th</sub>, etc.) are fixed (e.g., independent of temperature)</li>
    </ol></p>

<p>Additional notes:
<ul>
          <li>
  The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

  <p>The default resistivities (<code>F=1/(207.2e-7*U.Pa*U.s)</code> and <code>R=U.m*U.K/(26.8e-3*U.W)</code>) are based on data of gas at 1 atm and
  300 K from Incropera and DeWitt  [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 920&ndash;921].
  Table 1 lists the properties at other temperatures.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"bottom\"><b>Table 1:</b> Properties of O<sub>2</sub> gas at 1 atm [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002, pp. 920&ndash;921</a>]</caption>
  <tr>
      <th valign=\"middle\"><code>T/U.K</code></th>
      <th width=1><code>c_p*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1><code>F<br>*U.Pa*U.s</code></th>
      <th width=1><code>R*U.W<br>/(U.m*U.K)</code></th>
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
</p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));
        end Fixed;
      end Gas;
    end O2;

    model Species0Amount "Species with zero particle number"
      extends Species(
        final overrideEOS=true,
        final rho_IC=0,
        final derrho_IC=0,
        N(stateSelect=StateSelect.never),
        phi(each stateSelect=StateSelect.never),
        T(stateSelect=StateSelect.never));
      // Note:  StateSelect.never is necessary to avoid dynamic state selection
      // in Dymola 7.4.
      annotation (Documentation(info="<html>
  <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));
    end Species0Amount;

    model SpeciesInertStagnant "Inert and stagnant species"
      extends Species(
        final Xi=Modelica.Constants.inf,
        final F=0,
        final upstreamX,
        final upstreamY,
        final upstreamZ,
        final setPartNum=true,
        final setVelX=true,
        final setVelY=true,
        final setVelZ=true,
        initMethPartNum=InitMethScalar.Volume,
        final initMethX=InitMethVelocity.Velocity,
        final initMethY=InitMethVelocity.Velocity,
        final initMethZ=InitMethVelocity.Velocity,
        final Ndot_IC=0,
        final phi_IC=zeros(3),
        final derphi_IC,
        final I_IC,
        final derI_IC,
        xNegative(
          final isochoric=true,
          final inviscidY=true,
          final inviscidZ=true),
        xPositive(
          final isochoric=true,
          final inviscidY=true,
          final inviscidZ=true),
        yNegative(
          final isochoric=true,
          final inviscidZ=true,
          final inviscidX=true),
        yPositive(
          final isochoric=true,
          final inviscidZ=true,
          final inviscidX=true),
        zNegative(
          final isochoric=true,
          final inviscidX=true,
          final inviscidY=true),
        zPositive(
          final isochoric=true,
          final inviscidX=true,
          final inviscidY=true));

      // Note:  upstreamX, upstreamY, and upstreamZ don't matter since bulk
      // current is zero.

      annotation (Documentation(info="<html><p>Assumptions:<ol>
  <li>Zero velocity</li>
  <li>No material exchange or transport</li</ol>
  </p>

  <p><code>F</code> only affects exchange&mdash;not transport.  Since a no-slip
  boundary condition is applied and the bulk velocity is set to zero, there are
  never any shear forces, regardless of the value of <code>F</code>.</p>

  <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));
    end SpeciesInertStagnant;

    model Species
      "Model for single-species exchange, transport, and storage of material, linear momentum, and energy"
      //extends FCSys.BaseClasses.Icons.Names.Top1;

      // Geometric parameters
      outer parameter Q.Length L[Axis](each min=Modelica.Constants.small)
        "Length" annotation (HideResult=true, missingInnerMessage=
            "This model should be used within the Subregion model.");
      outer parameter Q.Area A[Axis] "Cross-sectional area" annotation (
          HideResult=true, missingInnerMessage=
            "This model should be used within the Subregion model.");
      parameter Q.Length Lstar(
        min=Modelica.Constants.small,
        nominal=10*U.m,
        start=1e3*product(L)^(1/3))
        "<html>Characteristic length for exchange (<i>L</i><sup>&#9733;</sup>)</html>"
        annotation (Dialog(group="Geometry"));
      parameter Q.NumberAbsolute k[Axis](
        each min=Modelica.Constants.small,
        each final nominal=1) = {1,1,1}
        "<html>Area fill factor for transport (<b>k</b>)</html>"
        annotation (HideResult=true, Dialog(group="Geometry"));

      // Material properties
      replaceable package Data =
          FCSys.Characteristics.BaseClasses.Characteristic constrainedby
        FCSys.Characteristics.BaseClasses.Characteristic
        "Characteristic data of the species" annotation (
        Dialog(group="Material properties"),
        __Dymola_choicesFromPackage=true,
        Placement(transformation(extent={{-60,40},{-40,60}}),
            iconTransformation(extent={{-10,90},{10,110}})));

      // Assumptions
      // -----------
      // General
      parameter Boolean overrideEOS=false
        "<html>Override the equation of state with the value of &rho;<sub>IC</sub></html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        Dialog(tab="Assumptions", compact=true),
        choices(__Dymola_checkBox=true));
      parameter Boolean inclLin[Axis]={true,false,false}
        "<html>true, if each component of linear momentum is included (<i>Do not adjust here.</i>)</html>"
        annotation (
        HideResult=true,
        Evaluate=true,
        Dialog(tab="Assumptions", enable=false));
      // Even though this parameter is set as final within the constrainedby
      // clauses of the Phase models, Dymola 7.4 still shows it in the
      // parameter dialog (hence the "Do not adjust").
      //
      // Upstream discretization
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
      //
      // Dynamics
      parameter Boolean setPartNum=false "Particle number" annotation (
        Evaluate=true,
        Dialog(
          tab="Assumptions",
          group="Prescribed states (via initialization parameters)",
          compact=true),
        choices(__Dymola_checkBox=true));
      parameter Boolean setVelX=false "X-axis component of velocity"
        annotation (
        Evaluate=true,
        Dialog(
          tab="Assumptions",
          group="Prescribed states (via initialization parameters)",
          enable=inclLin[1],
          compact=true),
        choices(__Dymola_checkBox=true));
      parameter Boolean setVelY=false "Y-axis component of velocity"
        annotation (
        Evaluate=true,
        Dialog(
          tab="Assumptions",
          group="Prescribed states (via initialization parameters)",
          enable=inclLin[2],
          compact=true),
        choices(__Dymola_checkBox=true));
      parameter Boolean setVelZ=false "Z-axis component of velocity"
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
      // Note:  This parameter is left enabled even it isn't used to
      // explicitly initialize any states, since it's used as a guess value.
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
          start=1/Data.v_Tp(T_IC, p_IC))
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
            Data.s(T=T_IC, p=p_IC))
        "<html>Initial specific entropy (<i>s</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization", group="Scalar properties"));
      parameter Q.NumberRate ders_IC=0
        "<html>Initial rate of specific entropy ((&part;<i>s</i>/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Scalar properties",
          enable=initMethPartNum == 13 or initMethTemp == 13));
      parameter Q.Potential h_IC(start=Data.h(T_IC, p_IC))
        "<html>Initial specific enthalpy (<i>h</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization", group="Scalar properties"));
      parameter Q.PotentialRate derh_IC=0
        "<html>Initial rate of specific enthalpy ((&part;<i>h</i>/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Scalar properties",
          enable=initMethPartNum == 15 or initMethTemp == 15));
      parameter Q.Potential mu_IC(start=Data.g(T_IC, p_IC))
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

      // Initialization parameters for velocity
      parameter BaseClasses.InitMethVelocity initMethX=InitMethVelocity.Velocity
        "Method of initializing the x-axis component" annotation (Dialog(
          tab="Initialization",
          group="Velocity",
          enable=inclLin[1]));
      parameter BaseClasses.InitMethVelocity initMethY=InitMethVelocity.Velocity
        "Method of initializing the y-axis component" annotation (Dialog(
          tab="Initialization",
          group="Velocity",
          enable=inclLin[2]));
      parameter BaseClasses.InitMethVelocity initMethZ=InitMethVelocity.Velocity
        "Method of initializing the z-axis component" annotation (Dialog(
          tab="Initialization",
          group="Velocity",
          enable=inclLin[3]));
      // Note:  Dymola 7.4 doesn't provide pull-down lists for arrays of
      // enumerations; therefore, a parameter is used for each axis.
      parameter Q.Velocity phi_IC[Axis]={0,0,0}
        "<html>Initial velocity (<b>&phi;</b><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Velocity"));
      parameter Q.Acceleration derphi_IC[Axis]={0,0,0}
        "<html>Initial acceleration ((&part;<b>&phi;</b>/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Velocity",
          enable=initMethX == 3 or initMethY == 3 or initMethZ == 3));
      parameter Q.Current I_IC[Axis]={0,0,0}
        "<html>Initial current (<i><b>I</b></i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Velocity"));
      parameter Q.CurrentRate derI_IC[Axis]={0,0,0}
        "<html>Initial rate of current ((&part;<i><b>I</b></i>/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Velocity",
          enable=initMethX == 5 or initMethY == 5 or initMethZ == 5));

      // Material properties
      Q.Fusivity Xi(nominal=1*U.s/U.cm^2) = Data.Xi(T, 1/rho)
        "<html>Fusivity (&Xi;)</html>"
        annotation (Dialog(group="Material properties"));
      Q.Fluidity F(nominal=10*U.cm*U.s/U.g) = Data.F(T, 1/rho) "Fluidity"
        annotation (Dialog(group="Material properties"));
      Q.ResistivityThermal R(nominal=10*U.cm/U.A) = Data.R(T, 1/rho)
        "Thermal resistivity" annotation (Dialog(group="Material properties"));

      // Preferred states
      Q.Amount N(
        nominal=1*U.mol,
        final start=N_IC,
        final fixed=false,
        stateSelect=StateSelect.default) "Particle number";
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
        nominal=1*U.cm^3,
        final start=V_IC,
        final fixed=false) "Volume";
      Q.PressureAbsolute p(
        nominal=1*U.atm,
        final start=p_IC,
        final fixed=false) "Pressure";
      // Note:  In Dymola 7.4 StateSelect.never is necessary to avoid dynamic
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
      Q.CapacityThermalSpecific c_V "Isochoric specific heat capacity";

      // Auxiliary variables (for analysis)
      // ----------------------------------
      // Misc. properties and conditions
      output Q.AmountVolumic rho(stateSelect=StateSelect.never) = N/V if
        defaults.analysis "Molar density";
      // Note:  The reciprocal, specific volume (v), isn't included because
      // particle number (N) can be zero.
      output Q.PressureAbsolute q[n_lin](each stateSelect=StateSelect.never) =
        Data.m*phi .* I ./ (2*A[cartAxes]) if defaults.analysis
        "Dynamic pressure";
      //
      // Capacitances
      output Q.Capacitance C(stateSelect=StateSelect.never) = if Data.isCompressible
         and not overrideEOS then -N*rho^2*(U.m^3/U.C)/Data.dp(
            v=1/rho,
            T=T,
            dv=U.m^3/U.C,
            dT=0) else 0 if defaults.analysis "Chemical capacitance";
      // Note:  This is delN/delg at constant T and V.
      output Q.CapacityThermal C_V(stateSelect=StateSelect.never) = N*c_V if
        defaults.analysis "Isochoric heat capacity";
      //
      // Time constants
      output Q.Time tau_exch_mechanical(stateSelect=StateSelect.never) =
        alpha_F*N/Lstar if defaults.analysis
        "Time constant for mechanical exchange";
      output Q.Time tau_exch_thermal(stateSelect=StateSelect.never) = alpha_R*N
        /Lstar if defaults.analysis "Time constant for thermal exchange";
      output Q.Time tau_trans_mat[Axis](each stateSelect=StateSelect.never) =
        fill(alpha_Xi*N, 3) ./ Lstar_trans if defaults.analysis
        "Time constants for material transport";
      output Q.Time tau_trans_mech[Axis](each stateSelect=StateSelect.never) =
        fill(alpha_F*N, 3) ./ Lstar_trans if defaults.analysis
        "Time constants for mechanical transport";
      output Q.Time tau_trans_therm[Axis](each stateSelect=StateSelect.never)
         = fill(alpha_R*N, 3) ./ Lstar_trans if defaults.analysis
        "Time constants for thermal transport";
      //
      // Peclet numbers (only for the axes with linear momentum included; others are
      // zero)
      output Q.Number Pe_mat[n_lin](each stateSelect=StateSelect.never) = I*
        alpha_Xi ./ Lstar_trans[cartAxes] if defaults.analysis
        "Material Peclet numbers";
      output Q.Number Pe_mech[n_lin](each stateSelect=StateSelect.never) = I*
        alpha_F ./ Lstar_trans[cartAxes] if defaults.analysis
        "Mechanical Peclet numbers";
      output Q.Number Pe_therm[n_lin](each stateSelect=StateSelect.never) = I*
        alpha_R ./ Lstar_trans[cartAxes] if defaults.analysis
        "Thermal Peclet numbers";
      //
      // Bulk flow rates
      output Q.Force mphiI[n_lin, Orientation](each stateSelect=StateSelect.never)
         = {(if inclLin[cartWrap(cartAxes[axis] + orientation)] then Data.m*phi[
        linAxes[cartWrap(cartAxes[axis] + orientation)]]*I[axis] else 0) for
        orientation in Orientation, axis in 1:n_lin} if n_lin > 0 and defaults.analysis
        "Bulk rate of mechanical advection";
      output Q.Power TsI[n_lin](each stateSelect=StateSelect.never) = T*s*I if
        defaults.analysis "Bulk rate of thermal advection";
      //
      // Linear momentum balance
      output Q.Force Ma[n_lin](each stateSelect=StateSelect.never) = M*der(phi)
        /U.s if defaults.analysis "Acceleration force (constant mass)";
      output Q.Force f_exch_adv[n_lin](each stateSelect=StateSelect.never) =
        chemical.mPhidot - Data.m*phi*chemical.Ndot if defaults.analysis
        "Acceleration force due to material (advective) exchange";
      output Q.Force f_exch_diff[n_lin](each stateSelect=StateSelect.never) =
        common.mechanical.mPhidot + inert.mPhidot if defaults.analysis
        "Friction from other species (diffusive exchange)";
      output Q.Force f_trans_adv[n_lin](each stateSelect=StateSelect.never) = {
        Data.m*Delta(Data.v_Tp(T_face[cartAxes[axis], :], p_face[cartAxes[axis],
        :]) .* Ndot_face[cartAxes[axis], :] .^ 2)/A[cartAxes[axis]] + sum(Data.m
        *Sigma(phi_face[cartWrap(cartAxes[axis] - orientation), :, orientation]
         .* Ndot_face[cartWrap(cartAxes[axis] - orientation), :]) for
        orientation in Orientation) for axis in 1:n_lin} if defaults.analysis
        "Acceleration force due to material (advective) transport";
      output Q.Force f_trans_diff[n_lin](each stateSelect=StateSelect.never) =
        {A[cartAxes[axis]]*Delta(p_face[cartAxes[axis], :]) + sum(Sigma(
        mPhidot_face[cartWrap(cartAxes[axis] - orientation), :, orientation])
        for orientation in Orientation) for axis in 1:n_lin} if defaults.analysis
        "Friction from other subregions (diffusive transport; includes volume viscosity)";
      //
      // Energy balance
      output Q.Power Ndere(stateSelect=StateSelect.never) = (N*(der(h) + Data.m
        *der(phi*phi)/2) - V*der(p))/U.s if defaults.analysis
        "Rate of energy storage (internal and kinetic) at constant mass";
      output Q.Power Wdot_exch(stateSelect=StateSelect.never) = -((Data.m*(
        chemical.hbar - phi*phi/2) - h)*chemical.Ndot + chemical.phi*chemical.mPhidot
        /2) if defaults.analysis
        "Relative rate of work (internal, flow, and kinetic) done by chemical exchange (advection)";
      output Q.Power Qdot_gen_exch(stateSelect=StateSelect.never) = phi*common.mechanical.mPhidot
         + inert.phi*inert.mPhidot if defaults.analysis
        "Rate of heat generation due to friction with other species";
      output Q.Power Qdot_exch(stateSelect=StateSelect.never) = common.thermal.Qdot
         + inert.Qdot if defaults.analysis
        "Rate of thermal conduction from other species";
      output Q.Power Wdot_trans(stateSelect=StateSelect.never) = -sum(sum((
        Data.h(T_face[axis, side], p_face[axis, side]) + Data.m*((Data.v_Tp(
        T_face[axis, side], p_face[axis, side])*Ndot_face[axis, side]/A[axis])^
        2 + phi_face[axis, side, :]*phi_face[axis, side, :])/2 - h - Data.m*phi
        *phi/2)*Ndot_face[axis, side] for side in Side) for axis in Axis) if
        defaults.analysis
        "Relative rate of work (internal, flow, and kinetic) done by material transport (advection)";
      output Q.Power Qdot_gen_trans(stateSelect=StateSelect.never) = sum(
        phi_face .* mPhidot_face) if defaults.analysis
        "Rate of heat generation due to friction with other subregions";
      output Q.Power Qdot_trans(stateSelect=StateSelect.never) = sum(Qdot_face)
        if defaults.analysis "Rate of thermal conduction from other subregions";
      // Note:  These auxiliary variables should not be used as states; the
      // structure of the problem should not change if they are included.

      FCSys.Connectors.ChemicalOutput chemical(
        final n_lin=n_lin,
        final m=Data.m,
        final formula=Data.formula,
        muPerT(final start=mu_IC/T_IC),
        phi(final start=phi_IC[cartAxes]),
        Ndot(final start=Ndot_IC,final fixed=false),
        hbar(final start=Data.h(T_IC, p_IC)/Data.m, final fixed=false))
        "Connector to exchange material while advecting linear momentum and energy"
        annotation (Placement(transformation(extent={{-30,10},{-10,30}}),
            iconTransformation(extent={{-80,60},{-60,80}})));
      FCSys.Connectors.Inert common(
        final n_lin=n_lin,
        mechanical(phi(final start=phi_IC[cartAxes], each final fixed=false)),
        thermal(T(final start=T_IC,final fixed=false)))
        "Connector for direct mechanical and thermal coupling of multiple species"
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
        T(start=T_IC))
        "Connector to add pressure and exchange linear momentum and heat by diffusion"
        annotation (Placement(transformation(extent={{10,-30},{30,-10}}),
            iconTransformation(extent={{60,-80},{80,-60}})));
      FCSys.Connectors.FaceX xNegative(
        isochoric=not inclLin[Axis.x],
        inviscidY=not inclLin[Axis.y],
        inviscidZ=not inclLin[Axis.z],
        material(final rho(start=rho_IC) = rho_face[Axis.x, Side.n], final Ndot(
              start=I_IC[Axis.x]) = Ndot_face[Axis.x, Side.n]),
        mechanicalY(final phi(start=phi_IC[Axis.y]) = phi_face[Axis.x, Side.n,
            Orientation.following], final mPhidot=mPhidot_face[Axis.x, Side.n,
              Orientation.following]),
        mechanicalZ(final phi(start=phi_IC[Axis.z]) = phi_face[Axis.x, Side.n,
            Orientation.preceding], final mPhidot=mPhidot_face[Axis.x, Side.n,
              Orientation.preceding]),
        thermal(final T(start=T_IC) = T_face[Axis.x, Side.n],final Qdot(start=0)
             = Qdot_face[Axis.x, Side.n])) "Negative face along the x axis"
        annotation (Dialog(
          tab="Assumptions",
          group="Characteristics of the faces (click to edit)",
          __Dymola_label="xNegative",
          __Dymola_descriptionLabel=true),Placement(transformation(extent={{-50,
                -10},{-30,10}}), iconTransformation(extent={{-110,-10},{-90,10}})));
      FCSys.Connectors.FaceX xPositive(
        isochoric=not inclLin[Axis.x],
        inviscidY=not inclLin[Axis.y],
        inviscidZ=not inclLin[Axis.z],
        material(final rho(start=rho_IC) = rho_face[Axis.x, Side.p], final Ndot(
              start=-I_IC[Axis.x]) = Ndot_face[Axis.x, Side.p]),
        mechanicalY(final phi(start=phi_IC[Axis.y]) = phi_face[Axis.x, Side.p,
            Orientation.following], final mPhidot=mPhidot_face[Axis.x, Side.p,
              Orientation.following]),
        mechanicalZ(final phi(start=phi_IC[Axis.z]) = phi_face[Axis.x, Side.p,
            Orientation.preceding], final mPhidot=mPhidot_face[Axis.x, Side.p,
              Orientation.preceding]),
        thermal(final T(start=T_IC) = T_face[Axis.x, Side.p],final Qdot(start=0)
             = Qdot_face[Axis.x, Side.p])) "Positive face along the x axis"
        annotation (Dialog(
          tab="Assumptions",
          group="Characteristics of the faces (click to edit)",
          __Dymola_label="xPositive",
          __Dymola_descriptionLabel=true),Placement(transformation(extent={{30,
                -10},{50,10}}), iconTransformation(extent={{90,-10},{110,10}})));
      FCSys.Connectors.FaceY yNegative(
        isochoric=not inclLin[Axis.y],
        inviscidZ=not inclLin[Axis.z],
        inviscidX=not inclLin[Axis.x],
        material(final rho(start=rho_IC) = rho_face[Axis.y, Side.n], final Ndot(
              start=I_IC[Axis.y]) = Ndot_face[Axis.y, Side.n]),
        mechanicalZ(final phi(start=phi_IC[Axis.z]) = phi_face[Axis.y, Side.n,
            Side.n], final mPhidot=mPhidot_face[Axis.y, Side.n, Orientation.following]),

        mechanicalX(final phi(start=phi_IC[Axis.x]) = phi_face[Axis.y, Side.n,
            Side.p], final mPhidot=mPhidot_face[Axis.y, Side.n, Orientation.preceding]),

        thermal(final T(start=T_IC) = T_face[Axis.y, Side.n],final Qdot(start=0)
             = Qdot_face[Axis.y, Side.n])) "Negative face along the y axis"
        annotation (Dialog(
          tab="Assumptions",
          group="Characteristics of the faces (click to edit)",
          __Dymola_label="yNegative",
          __Dymola_descriptionLabel=true),Placement(transformation(extent={{-10,
                -50},{10,-30}}), iconTransformation(extent={{-10,-110},{10,-90}})));

      FCSys.Connectors.FaceY yPositive(
        isochoric=not inclLin[Axis.y],
        inviscidZ=not inclLin[Axis.z],
        inviscidX=not inclLin[Axis.x],
        material(final rho(start=rho_IC) = rho_face[Axis.y, Side.p],final Ndot(
              start=-I_IC[Axis.y]) = Ndot_face[Axis.y, Side.p]),
        mechanicalZ(final phi(start=phi_IC[Axis.z]) = phi_face[Axis.y, Side.p,
            Orientation.following], final mPhidot=mPhidot_face[Axis.y, Side.p,
              Orientation.following]),
        mechanicalX(final phi(start=phi_IC[Axis.x]) = phi_face[Axis.y, Side.p,
            Orientation.preceding], final mPhidot=mPhidot_face[Axis.y, Side.p,
              Orientation.preceding]),
        thermal(final T(start=T_IC) = T_face[Axis.y, Side.p],final Qdot(start=0)
             = Qdot_face[Axis.y, Side.p])) "Positive face along the y axis"
        annotation (Dialog(
          tab="Assumptions",
          group="Characteristics of the faces (click to edit)",
          __Dymola_label="yPositive",
          __Dymola_descriptionLabel=true),Placement(transformation(extent={{-10,
                30},{10,50}}), iconTransformation(extent={{-10,90},{10,110}})));
      FCSys.Connectors.FaceZ zNegative(
        isochoric=not inclLin[Axis.z],
        inviscidX=not inclLin[Axis.x],
        inviscidY=not inclLin[Axis.y],
        material(final rho(start=rho_IC) = rho_face[Axis.z, Side.n], final Ndot(
              start=I_IC[Axis.z]) = Ndot_face[Axis.z, Side.n]),
        mechanicalX(final phi(start=phi_IC[Axis.x]) = phi_face[Axis.z, Side.n,
            Orientation.following], final mPhidot=mPhidot_face[Axis.z, Side.n,
              Orientation.following]),
        mechanicalY(final phi(start=phi_IC[Axis.y]) = phi_face[Axis.z, Side.n,
            Orientation.preceding], final mPhidot=mPhidot_face[Axis.z, Side.n,
              Orientation.preceding]),
        thermal(final T(start=T_IC) = T_face[Axis.z, Side.n],final Qdot(start=0)
             = Qdot_face[Axis.z, Side.n])) "Negative face along the z axis"
        annotation (Dialog(
          tab="Assumptions",
          group="Characteristics of the faces (click to edit)",
          __Dymola_label="zNegative",
          __Dymola_descriptionLabel=true),Placement(transformation(extent={{10,
                10},{30,30}}), iconTransformation(extent={{60,60},{80,80}})));
      FCSys.Connectors.FaceZ zPositive(
        isochoric=not inclLin[Axis.z],
        inviscidX=not inclLin[Axis.x],
        inviscidY=not inclLin[Axis.y],
        material(final rho(start=rho_IC) = rho_face[Axis.z, Side.p], final Ndot(
              start=-I_IC[Axis.z]) = Ndot_face[Axis.z, Side.p]),
        mechanicalX(final phi(start=phi_IC[Axis.x]) = phi_face[Axis.z, Side.p,
            Orientation.following], final mPhidot=mPhidot_face[Axis.z, Side.p,
              Orientation.following]),
        mechanicalY(final phi(start=phi_IC[Axis.y]) = phi_face[Axis.z, Side.p,
            Orientation.preceding], final mPhidot=mPhidot_face[Axis.z, Side.p,
              Orientation.preceding]),
        thermal(final T(start=T_IC) = T_face[Axis.z, Side.p],final Qdot(start=0)
             = Qdot_face[Axis.z, Side.p])) "Positive face along the z axis"
        annotation (Dialog(
          tab="Assumptions",
          group="Characteristics of the faces (click to edit)",
          __Dymola_label="zPositive",
          __Dymola_descriptionLabel=true),Placement(transformation(extent={{-30,
                -30},{-10,-10}}), iconTransformation(extent={{-80,-80},{-60,-60}})));

      // Geometric parameters
    protected
      final parameter Q.Length Lstar_trans[Axis]=k .* A ./ L
        "Effective cross-sectional area per length";
      final parameter Integer n_lin=countTrue(inclLin)
        "Number of components of linear momentum";
      final parameter Integer cartAxes[n_lin]=index(inclLin)
        "Cartesian-axis indices of the axes of linear momentum";
      final parameter Integer linAxes[Axis]=enumerate(inclLin)
        "Linear momentum component indices of the Cartesian axes";
      final parameter Boolean upstream[Axis]={upstreamX,upstreamY,upstreamZ}
        "true, if each Cartesian axis uses upstream discretization";
      final parameter Boolean setVel[Axis]={setVelX,setVelY,setVelZ}
        "true, if each component of linear momentum is prescribed";
      final parameter BaseClasses.InitMethVelocity initMethVel[Axis]={initMethX,
          initMethY,initMethZ} "Initialization methods for velocity";

      // Base resistivity factors
      Q.Resistivity alpha_Xi(nominal=10*U.cm/U.A) = Xi/rho
        "Base resistivity factor for fusivity";
      Q.Resistivity alpha_F(nominal=10*U.cm/U.A) = F*Data.m
        "Base resistivity factor for fluidity";
      Q.Resistivity alpha_R(nominal=10*U.cm/U.A) = R*c_V
        "Base resistivity factor for thermal resistivity";

      // Efforts and flows of the conditional connectors
      Q.Force rho_face[Axis, Side](start=fill(
                rho_IC,
                3,
                2)) "Densities at the faces";
      Q.Pressure p_face[Axis, Side](start=fill(
                p_IC,
                3,
                2)) "Pressures at the faces";
      Q.Current Ndot_face[Axis, Side](start=outerProduct(I_IC, {1,-1}))
        "Currents into the faces";
      Q.Velocity phi_face[Axis, Side, Orientation](start={{{if inclLin[cartWrap(
            axis + orientation)] then phi_IC[cartWrap(axis + orientation)]
             else 0 for orientation in Orientation} for side in Side} for axis
             in Axis}) "Shear velocities at the faces";
      Q.Force mPhidot_face[Axis, Side, Orientation](start={fill({phi_IC[
            cartWrap(axis + orientation)] for orientation in Orientation}, 2)
            for axis in Axis}) "Shear forces on the faces";
      Q.TemperatureAbsolute T_face[Axis, Side](each start=T_IC)
        "Temperatures at the faces";
      Q.Power Qdot_face[Axis, Side] "Heat flow rates into the faces";

      outer FCSys.BCs.Defaults defaults "Default settings" annotation (
          missingInnerMessage="Your model is using an outer \"defaults\" record, but an inner \"defaults\" record is not defined.
For simulation, specify global default settings by dragging FCSys.BCs.Defaults into your model.
The default global default settings will be used for the current simulation.",
          Placement(transformation(extent={{40,40},{60,60}}),
            iconTransformation(extent={{-10,90},{10,110}})));
      // Note:  In Dymola 7.4 it's necessary to add the missing inner message
      // here to give a warning message, even though it's included in the
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
  if abs(D) < Modelica.Constants.small then
    Modelica.Utilities.Streams.print("Warning: The self diffusivity is zero.
    This may directly couple the currents within neighboring subregions.\nConsider setting the value of D as final (if not already) so that index reduction may be performed.");
  end if;
  if abs(F) > Modelica.Constants.small then
    Modelica.Utilities.Streams.print("Warning: The fluidity is zero.
    This may directly couple the velocity of this species with others within the subregion or with the same species within neighboring subregions.\nConsider setting the value of F as final (if not already) so that index reduction may be performed.");
  end if;
  if abs(R) > Modelica.Constants.small then
    Modelica.Utilities.Streams.print("Warning: The thermal resistance is zero.
    This may directly couple the temperature of this species with others within the subregion or with the same species within neighboring subregions.\nConsider setting the value of R as final (if not already) so that index reduction may be performed.");
  end if;
  // Note:  According to the Modelica >=3.0 specification, these
  // checks should be possible using the assert() command with
  // level=AssertionLevel.warning.  However, this isn't supported in
  // Dymola 7.4 or FD2012.
  */

      // Particle number
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
          der(mu)/U.s = dermu_IC;
        elseif initMethPartNum == InitMethScalar.ReactionRate then
          chemical.Ndot = Ndot_IC;
          // Else, initMethPartNum == InitMethScalar.None; then, there are no
          // initial equations.
        end if;
      end if;

      // Velocity
      for axis in Axis loop
        if inclLin[axis] then
          if setVel[axis] then
            // Ensure that a condition is selected, since the state is
            // prescribed.
            assert(initMethVel[axis] <> InitMethVelocity.None,
              "The state for the " + {"x","y","z"}[axis] + "-axis component of linear momentum is prescribed,
        yet its condition is not defined.\nChoose any condition besides None.");
          elseif not overrideEOS or rho_IC > 0 then
            // Initialize since there's a time-varying state.
            if initMethVel[axis] == InitMethVelocity.Velocity then
              phi[linAxes[axis]] = phi_IC[axis];
            elseif initMethVel[axis] == InitMethVelocity.Acceleration then
              der(phi[linAxes[axis]])/U.s = derphi_IC[axis];
            elseif initMethX == InitMethVelocity.Current then
              I[linAxes[axis]] = I_IC[axis];
            elseif initMethVel[axis] == InitMethVelocity.CurrentRate then
              der(I[linAxes[axis]])/U.s = derI_IC[axis];
              // Else, initMethVel[axis] == InitMethVelocity.None; then, there are
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
          der(mu)/U.s = dermu_IC;
        elseif initMethTemp == InitMethScalar.ReactionRate then
          chemical.Ndot = Ndot_IC;
          // Else, initMethTemp == InitMethScalar.None; then, there are no
          // initial equations.
        end if;
      end if;

    equation
      // Aliases (only for clarity)
      p = inert.p;
      V = inert.V;
      T = common.thermal.T;
      phi = common.mechanical.phi;
      mu = T*chemical.muPerT;
      h = mu + T*s;
      N*phi = L[cartAxes] .* I;
      M = Data.m*N;
      p_face = {Data.p_Tv(T_face[axis, :], {1,1} ./ rho_face[axis, :]) for axis
         in Axis};
      c_V = Data.c_V(T, p);

      // Thermodynamic correlations
      if overrideEOS then
        N = rho_IC*V;
      elseif Data.isCompressible then
        p = Data.p_Tv(T, V/N);
      else
        V = N*Data.v_Tp(T, p);
      end if;
      h = Data.h(T, p);
      s = Data.s(T, p);

      // Exchange
      // --------
      // Material
      chemical.mPhidot = semiLinear(
            Data.m*chemical.Ndot,
            chemical.phi,
            phi) "Advection";
      F*inert.mPhidot = 2*Lstar*(inert.phi - phi) "Diffusion";
      //
      // Fluid/thermal
      chemical.Hdot = semiLinear(
            chemical.Ndot,
            chemical.hbar*Data.m,
            h) "Advection";
      R*inert.Qdot = 2*Lstar*(inert.T - T) "Diffusion";

      // Transport
      for axis in Axis loop
        for side in Side loop
          // Material
          Xi*(Ndot_face[axis, side] - (if inclLin[axis] then inSign(side)*I[
            linAxes[axis]] else 0)) = Lstar_trans[axis]*(rho_face[axis, side]
             - rho)*(if upstream[axis] and inclLin[axis] then (exp(inSign(side)
            *I[linAxes[axis]]*alpha_Xi/(2*Lstar_trans[axis])) + 1) else 2);
          if [xNegative.isochoric, xPositive.isochoric; yNegative.isochoric,
              yPositive.isochoric; zNegative.isochoric, zPositive.isochoric][
              axis, side] then
            rho_face[axis, side] = rho "Isochoric BC";
          end if;

          // Mechanical
          for orientation in Orientation loop
            F*mPhidot_face[axis, side, orientation] = 4*Lstar_trans[axis]*(
              phi_face[axis, side, orientation] - (if inclLin[cartWrap(axis +
              orientation)] then phi[linAxes[cartWrap(axis + orientation)]]
               else 0))*(if upstream[axis] and inclLin[axis] then (exp(inSign(
              side)*I[linAxes[axis]]*alpha_F/(2*Lstar_trans[axis])) + 1) else 2);
            if {{{xNegative.inviscidY,xNegative.inviscidZ},{xPositive.inviscidY,
                xPositive.inviscidZ}},{{yNegative.inviscidZ,yNegative.inviscidX},
                {yPositive.inviscidZ,yPositive.inviscidX}},{{zNegative.inviscidX,
                zNegative.inviscidY},{zPositive.inviscidX,zPositive.inviscidY}}}
                [axis, side, orientation] then
              mPhidot_face[axis, side, orientation] = 0 "Inviscid BC";
            end if;
          end for;

          // Thermal
          R*Qdot_face[axis, side] = Lstar_trans[axis]*(T_face[axis, side] - T)*
            (if upstream[axis] and inclLin[axis] then (exp(inSign(side)*I[
            linAxes[axis]]*alpha_R/(2*Lstar_trans[axis])) + 1) else 2);
          if [xNegative.adiabatic, xPositive.adiabatic; yNegative.adiabatic,
              yPositive.adiabatic; zNegative.adiabatic, zPositive.adiabatic][
              axis, side] then
            Qdot_face[axis, side] = 0 "Adiabatic BC";
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
          der(mu)/U.s = dermu_IC;
        else
          //if initMethPartNum == InitMethScalar.ReactionRate then
          chemical.Ndot = Ndot_IC;
          // Note:  initMethPartNum == InitMethScalar.None can't occur due to an
          // assertion.
        end if;
      else
        der(N)/U.s = chemical.Ndot + sum(Ndot_face) "Material conservation";
      end if;

      // Mechanical dynamics
      for axis in 1:n_lin loop
        if setVel[cartAxes[axis]] then
          // Apply the IC for all time (linear momentum isn't conserved along
          // this axis).
          if initMethVel[cartAxes[axis]] == InitMethVelocity.Velocity then
            phi[axis] = phi_IC[cartAxes[axis]];
          elseif initMethVel[cartAxes[axis]] == InitMethVelocity.Acceleration
               then
            der(phi[axis])/U.s = derphi_IC[cartAxes[axis]];
          elseif initMethX == InitMethVelocity.Current then
            I[axis] = I_IC[cartAxes[axis]];
          elseif initMethVel[cartAxes[axis]] == InitMethVelocity.CurrentRate
               then
            der(I[axis])/U.s = derI_IC[cartAxes[axis]];
            // Note:  initMethVel[cartAxes[axis]] == InitMethVelocity.None can't
            // occur due to an assertion.
          end if;
        else
          der(M*phi[axis])/U.s = chemical.mPhidot[axis] + common.mechanical.mPhidot[
            axis] + inert.mPhidot[axis] + A[cartAxes[axis]]*Delta(p_face[
            cartAxes[axis], :]) + Data.m*Delta(Data.v_Tp(T_face[cartAxes[axis],
            :], p_face[cartAxes[axis], :]) .* Ndot_face[cartAxes[axis], :] .^ 2)
            /A[cartAxes[axis]] + sum(Data.m*Sigma(phi_face[cartWrap(cartAxes[
            axis] - orientation), :, orientation] .* Ndot_face[cartWrap(
            cartAxes[axis] - orientation), :]) + Sigma(mPhidot_face[cartWrap(
            cartAxes[axis] - orientation), :, orientation]) for orientation in
            Orientation) "Conservation of linear momentum";
        end if;
      end for;

      // Thermal dynamics
      if setTemp then
        // Apply the IC for all time (energy not conserved).
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
          der(mu)/U.s = dermu_IC;
        else
          //if initMethTemp == InitMethScalar.ReactionRate then
          chemical.Ndot = Ndot_IC;
          // Note:  initMethTemp == InitMethScalar.None can't occur due to an
          // assertion.
        end if;
      else
        (der(N*h) + der(M*phi*phi)/2 - V*der(p))/U.s = chemical.phi*chemical.mPhidot
          /2 + Data.m*chemical.hbar*chemical.Ndot + phi*common.mechanical.mPhidot
           + common.thermal.Qdot + inert.phi*inert.mPhidot + inert.Qdot + sum(
          sum((Data.h(T_face[axis, side], p_face[axis, side]) + Data.m*((
          Data.v_Tp(T_face[axis, side], p_face[axis, side])*Ndot_face[axis,
          side]/A[axis])^2 + phi_face[axis, side, :]*phi_face[axis, side, :])/2)
          *Ndot_face[axis, side] for side in Side) for axis in Axis) + sum(
          phi_face .* mPhidot_face) + sum(Qdot_face) "Energy conservation";
      end if;
      annotation (
        defaultComponentPrefixes="replaceable",
        Documentation(info="<html>
    <p>This model is based on the following fixed assumptions.  Other assumptions are
    optional via the parameters.
    <ol>
       <li>All faces are rectangular.
       <li>The material is orthorhombic.  This implies that a
          gradient which induces diffusion along an axis does not induce
          diffusion along axes orthogonal to it
          [<a href=\"modelica://FCSys.UsersGuide.References\">Bejan2006</a>,
          pp. 691&ndash;692].</li>
       <li>The coordinate system (x, y, z) is aligned with the principle
          axes of transport.  For example, if the species is stratified, the
          layers must be parallel to one of the planes in the rectilinear
          grid.</li>
       <li>The factors that may cause anisotropic behavior (<b><i>k</i></b>)
          are common to material, mechanical, and thermal transport.</li>
       <li>There are no body or inertial forces (e.g., gravity).</li>
    </ol>
    </p>

    <p>Figure 1 shows the manner in which instances of
    <a href=\"modelica://FCSys.Subregions.Species\">Species</a> models (derived from this
    model) are
    connected within a <a href=\"modelica://FCSys.Subregions\">Subregion</a>.  The
    generalized resistances (<i>R</i>) affect the flow rates of linear momentum and
    heat associated with differences in velocity and temperature (respectively) between
    each species and a common node.  This exchange is diffusive.

    <p>Linear momentum and enthalpy are advected as material is exchanged in a chemical
    reaction.  This occurs at the velocity and massic enthalpy of the reactants (source
    species).  Resistance is not involved directly.  When the species are connected through
    a <a href=\"modelica://FCSys.Subregions.Reaction\">Reaction</a> model, the material
    states are directly coupled according to the chemical equilibrium of the reaction.
    This reduces the DAE index by one.  The reaction rate is determined solely by
    the dynamics of the transport equations.</p>

    <p align=center><img src=\"modelica://FCSys/resources/documentation/Subregions/Species/Species/exchange.png\">
<br><b>Figure 1:</b>  Exchange of a quantity (linear momentum or heat) among species
    (A, B, and C) within a subregion.</p>

    <p>Figure 2 shows how <a href=\"modelica://FCSys.Subregions.Species\">Species</a>
    instances of the same type are connected between neighboring
    <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> instances.  Material,
    linear momentum, and heat are transported by both advection and diffusion.
    Upstream discretization is applied if it is enabled (via the <code>upstreamX</code>,
    etc. parameters).</p>

    <p align=center><img src=\"modelica://FCSys/resources/documentation/Subregions/Species/Species/transport.png\">
<br><b>Figure 2:</b>  Transport of a quantity associated with the same species
    between subregions (1 and 2).</p>

    <p>Within a phase, <a href=\"modelica://FCSys.Subregions.Species\">Species</a> instances
    are combined by Dalton's law (see the
    <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector), as shown
    in Figure 3a.  The pressures are additive, and each species is assumed to exist at the
    volume of the phase.  Within a subregion, phases are combined by Amagat's law (see the
    <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector), as shown
    in Figure 3b.  The volumes are additive, and each species is assumed to exist at the
    pressure of the subregion.</p>

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

    <p> The following notes apply to the parameters:
    <ul>
    <li>The \"specific\" adjective is taken to mean a quantity divided by particle
    number.  (\"Massic\" would indicate a quantity divided by mass.)</li>
    <li>In general, if fluidity, self diffusivity, or thermal resistivity is zero, then
    it should be set as <code>final</code> so that index reduction may be performed.
    If two <a href=\"modelica://FCSys.Subregions.Species\">Species</a> instances
    are connected through their exchange connectors
    (<code>chemical</code> or <code>inert</code>) or faces (<code>xNegative</code>, <code>xPositive</code>, etc.) and both have zero resistivities for a
    quantity, then index reduction is necessary.</li>
    <li>Even if an initialization parameter is not selected for explicit use,
    it may be used a guess value.</li>
    <li>The <b><i>k</i></b> factor can be used to account for the effects of porosity and tortousity
    on transport.
    It should be changed directly with effective area and inversely with effective length.
    The factor may reflect anisotropic properties; it is a vector with independent components
    for each axis. It affects all of the diffusive transport rates (material, mechanical, and
    thermal) by the same factor.  By default, its components are unity.</li>
    <li>By default, only the x-axis component of linear momentum is included.  Also by default,
    only material and thermal transport are included through the x-axis faces and only
    x-axis displacement/shear stress is included through the y- and z-axis faces.</li>
    <li>By default, faces are assumed to be isobaric and invisicid if the corresponding 
    component of linear momenum is disabled.</li>
    <li>If a state is prescribed, then the
    associated initial condition (IC) will be applied for all time.  The
    corresponding conservation equation will not be imposed.
    If <code>setPartNum</code>, <code>setVelX</code>, <code>setVelY</code>, or <code>setVelZ</code> is
    <code>true</code>, then there may be a secondary effect on the energy conservation equation
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
    <li>With the <code>overrideEOS</code> parameter, it is possible to specify that
    the volumic amount (i.e., molar concentration) and thus the amount or particle number is
    zero.
    Set <code>overrideEOS = true</code> and <code>rho_IC = 0</code>; then,
    the states for material, linear momentum, and energy will be eliminated.  If a species
    is included with this setting, then there must be an external reference
    for electrochemical potential (i.e., ground).  There must be at least one other
    species in the subregion or the velocity must be set (e.g.,
    <code>setVelX = true</code>).</li>
    </p>

    <p>In evaluating the dynamics of a phase, it is usually assumed that all of the species
    exist at the same temperature and bulk velocity.  The time constants that govern the
    temperatures/heat capacities of the species and heat flow rates among them are usually
    much shorter than the time span of interest.
    This assumption can be applied in the model by connecting the <code>common</code>
    connectors of the species.  It will cause index reduction during translation.</p>

    <p>In the variables that relate to transport,
    the first index is the axis and the second index is the side.  The sides
    are ordered from negative to positive, according to the
    <a href=\"modelica://FCSys.BaseClasses.Side\">Side</a> enumeration.
    Shear velocity and force are additionally indexed by
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

      type InitMethVelocity = enumeration(
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
      "Connector for volume, linear momentum, and heat&mdash;with Amagat's law"
      annotation (Placement(transformation(extent={{60,-80},{80,-60}}),
          iconTransformation(extent={{110,-130},{130,-110}})));
    FCSys.Connectors.InertDalton inertD(final n_lin=n_lin)
      "Connector for volume, linear momentum, and heat&mdash;with Dalton's law"
      annotation (Placement(transformation(extent={{30,-50},{50,-30}}),
          iconTransformation(extent={{60,-80},{80,-60}})));

  equation
    // Equal intensive properties
    inertA.phi = inertD.phi;
    inertA.T = inertD.T;

    // Static balances
    0 = inertA.p + inertD.p "Pressure";
    0 = inertA.V + inertD.V "Volume";

    // Rate balances (without storage or generation)
    zeros(n_lin) = inertA.mPhidot + inertD.mPhidot "Linear momentum";
    0 = inertA.Qdot + inertD.Qdot "Energy";
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
    Q.Current Ndot(nominal=1*U.A) "Reaction rate";

    FCSys.Connectors.ChemicalInput chemical[n_spec](each final n_lin=n_lin)
      "Connector for chemical species"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

  protected
    Q.Velocity phi[n_lin](each nominal=1*U.cm/U.s,each start=0) "Velocity";
    Q.Velocity2 hbar(nominal=1*U.V*U.mol/U.g, start=0) "Massic enthalpy";

  initial equation
    assert(abs(nu*chemical.m) < 1e-4*U.g/U.mol, "The mass balance is incorrect.
Check the chemical formulas and the specific masses of the species.");

  equation
    // Chemical equilibrium
    0 = nu*chemical.muPerT;

    // Conservation (without storage)
    nu[1:n_spec]*Ndot = chemical.Ndot "Material";
    zeros(n_lin) = sum(chemical[i].mPhidot for i in 1:n_spec) "Linear momentum";
    0 = sum(chemical.Hdot) "Energy";

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
    // TODO:  Use stream connectors once they are better supported (some
    // errors occurred in Dymola 7.4).

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

  <p>Linear momentum and energy are advected using the <code>semiLinear</code> operator.
  The rate of advection of linear momentum is the
  product of the velocity of the source and the rate of mass
  (<i>m</i> &phi; <i>N&#775;</i>).  The rate of thermal advection is the
  product of the massic enthalpy of the source and the rate of mass
  (<i>m</i> <i>h&#772;</i> <i>N&#775;</i>).  If there multiple sources, then
  their contributions are additive.  If there are multiple sinks, then
  the flow is split on a mass basis.</p>

  <p>At uniform temperature, the stoichiometrically weighted sum of the chemical
  potentials zero.  The <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model
  specifies that the chemical potential is the Gibbs potential (<i>g</i> = <i>h</i> - <i>sT</i>).
  Therefore, assuming uniform temperature, the rate of advection of enthalpy (<i>H&#775;</i>) from reactants
  to products (or vice versa) is also temperature times the rate of advection of entropy (<i>TS&#775;</i>).</p>

    <p>For material, this model is essentially the opposite of a standard single-species connection.
    The stoichiometric sum of the efforts (&Sigma; &nu;<sub><i>i</i></sub> &mu;<sub><i>i</i></sub>)
    is zero, which is analogous to Kirchhoff's Current Law.  The flow rates divided by the
    stoichiometric coefficients (<i>N&#775;</i><sub><i>i</i></sub> /&nu;<sub><i>i</i></sub>)
    are equal&mdash;analogous to Kirchhoff's Voltage Law.</p>

    <p>Linear momentum and enthalpy are advected using the <code>semiLinear()</code> operator.  There is no diffusion;
    it is included in the inert connections among species
    (see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model).<p>

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
    parameter Q.Volume V(start=1*U.cm^3) "Volume";
    parameter Integer n_lin(
      final min=0,
      final max=3) = 1
      "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>"
      annotation (Evaluate=true, HideResult=true);

    FCSys.Connectors.InertAmagat inert(final n_lin=n_lin)
      "Connector for linear momentum and heat,  with shared volume" annotation
      (Placement(transformation(extent={{60,-80},{80,-60}}), iconTransformation(
            extent={{100,-120},{120,-100}})));
    outer FCSys.BCs.Defaults defaults "Default settings" annotation (Placement(
          transformation(extent={{40,40},{60,60}}), iconTransformation(extent={
              {-10,90},{10,110}})));

  equation
    // Specified volume
    V = inert.V;

    // Rate balances (without storage or generation)
    zeros(n_lin) = inert.mPhidot "Linear momentum";
    0 = inert.Qdot "Energy";
    annotation (
      Documentation(info="<html><p>This model uses a <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector that imposes
    additivity of volume.  In order to convert to additivity of pressure, use
    the <a href=\"modelica://FCSys.Subregions.PhaseBoundary\">PhaseBoundary</a> model.</p>

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

      FCSys.Subregions.Volume volume(final n_lin=n_lin, final V=V)
        "Model to establish space for species"
        annotation (Placement(transformation(extent={{-16,-16},{16,16}})));
    protected
      final parameter Integer n_lin=countTrue({inclLinX,inclLinY,inclLinZ})
        "Number of components of linear momentum" annotation (Evaluate=true);

      annotation (Documentation(info="<html><p>Notes:
  <ul><li>This model must be be extended so that models can be added for
  relevant species, phases, and reactions.</li>
  <li>Material will be transported between two subregions only if both of the connected faces are marked
  as open (<code>thermoOpt==ThermoOpt.OpenDiabatic</code>)
  within the instances of the matched <a href=\"modelica://FCSys.Subregions.Species\">Species</a> models.
  If either or both of the faces are closed (<code>thermoOpt==ThermoOpt.ClosedAdiabatic</code> or
  <code>thermoOpt==ThermoOpt.ClosedDiabatic</code>), then the interface will be closed.
  note applies to the viscous/inviscous and diabatic/adiabatic properties.</li>
  <li>The x-axis component of linear momentum is included by default.  At least one component must be included.</li>
  </ul></p></html>"), Icon(graphics={
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
