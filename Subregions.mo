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
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        graphite(
          final inclC=inclC,
          final 'incle-'='incle-',
          C(V_IC=subregion.V/4)),
        ionomer(
          final inclC19HF37O5S=inclC19HF37O5S,
          final 'inclH+'='inclH+',
          C19HF37O5S(V_IC=subregion.V/4)),
        inclYMom=false,
        inclZMom=false,
        inclYFaces=false,
        inclZFaces=false)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      inner FCSys.BCs.Defaults defaults(analysis=true)
        annotation (Placement(transformation(extent={{30,20},{50,40}})));
      FCSys.BCs.Face.SubregionFlow bC1(
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        graphite(final inclC=inclC, final 'incle-'='incle-'),
        ionomer(final inclC19HF37O5S=inclC19HF37O5S, final 'inclH+'='inclH+'))
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-20,0})));

      FCSys.BCs.Face.SubregionFlow bC2(
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        graphite(final inclC=inclC, final 'incle-'='incle-'),
        ionomer(final inclC19HF37O5S=inclC19HF37O5S, final 'inclH+'='inclH+'))
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={20,0})));

    equation
      connect(bC2.face, subregion.positiveX) annotation (Line(
          points={{16,1.23436e-15},{17.5,1.23436e-15},{17.5,1.22125e-15},{15,
              1.22125e-15},{15,6.10623e-16},{10,6.10623e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));
      connect(bC1.face, subregion.negativeX) annotation (Line(
          points={{-16,3.65701e-16},{-17.5,3.65701e-16},{-17.5,1.11022e-15},{-15,
              1.11022e-15},{-15,6.10623e-16},{-10,6.10623e-16}},
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

    model SubregionCAndH2
      "<html>Test a subregion with C and H<sub>2</sub></html>"
      extends SubregionH2(inclC=true, subregion(gas(H2(T_IC=310*U.K))));

      annotation (
        experiment(StopTime=0.2, Algorithm="Dassl"),
        experimentSetupOutput,
        Diagram(graphics),
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionCAndH2.mos"));
    end SubregionCAndH2;

    model SubregionHOR
      "<html>Test a subregion with the hydrogen oxidation reaction and the essential species for it (C, C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S, e<sup>-</sup>, H<sub>2</sub>, and H<sup>+</sup>)</html>"
      extends SubregionH2(
        inclC=true,
        inclC19HF37O5S=true,
        'incle-'=true,
        inclH2=true,
        'inclH+'=true);

      annotation (
        experiment(StopTime=5, Tolerance=1e-06),
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
        subregion(graphite('e-'(partNumInitMeth=InitMethScalar.ParticleNumber,
                N_IC=1*U.q)), ionomer('H+'(partNumInitMeth=InitMethScalar.ParticleNumber,
                N_IC=1*U.q))));
      // Electrons and protons are initially nearly depleted in order
      // to mostly balance the enthalpy of formation of water.  Otherwise,
      // the reaction will be tremendously fast (consuming electrons, protons,
      // and oxygen) and the simulation may fail.
      annotation (
        experiment(StopTime=0.01, Tolerance=1e-06),
        experimentSetupOutput,
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionORR .mos"));
      // Note:  The tolerance is currently very important.  A tolerance of
      // 1e-4 will likely fail except at very small stop times.
    end SubregionORR;

    model SubregionsH2 "Test a one-dimensional array of subregions"
      extends Modelica.Icons.Example;

      parameter Integer n_x=0
        "Number of discrete regions along the x axis, besides the 2 face regions";
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
        graphite(
          final inclC=inclC,
          final 'incle-'='incle-',
          C(V_IC=subregion1.V/4),
          'e-'(
            p_IC=1.05*defaults.p,
            negativeY(viscousZ=false, viscousX=false),
            positiveY(viscousZ=false, viscousX=false),
            negativeZ(viscousX=false, viscousY=false),
            positiveZ(viscousX=false, viscousY=false))),
        ionomer(
          final inclC19HF37O5S=inclC19HF37O5S,
          final 'inclH+'='inclH+',
          C19HF37O5S(V_IC=subregion1.V/4),
          'H+'(
            p_IC=1.05*defaults.p,
            negativeY(viscousZ=false, viscousX=false),
            positiveY(viscousZ=false, viscousX=false),
            negativeZ(viscousX=false, viscousY=false),
            positiveZ(viscousX=false, viscousY=false))),
        inclYMom=false,
        inclZMom=false,
        inclYFaces=false,
        inclZFaces=false,
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2O(
            p_IC=1.05*defaults.p,
            negativeY(viscousZ=false, viscousX=false),
            positiveY(viscousZ=false, viscousX=false),
            negativeZ(viscousX=false, viscousY=false),
            positiveZ(viscousX=false, viscousY=false)),
          N2(
            p_IC=1.05*defaults.p,
            negativeY(viscousZ=false, viscousX=false),
            positiveY(viscousZ=false, viscousX=false),
            negativeZ(viscousX=false, viscousY=false),
            positiveZ(viscousX=false, viscousY=false)),
          O2(
            p_IC=1.05*defaults.p,
            negativeY(viscousZ=false, viscousX=false),
            positiveY(viscousZ=false, viscousX=false),
            negativeZ(viscousX=false, viscousY=false),
            positiveZ(viscousX=false, viscousY=false)),
          H2(
            p_IC=1.05*defaults.p,
            negativeY(viscousZ=false, viscousX=false),
            positiveY(viscousZ=false, viscousX=false),
            negativeZ(viscousX=false, viscousY=false),
            positiveZ(viscousX=false, viscousY=false))))
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
            negativeY(viscousZ=false, viscousX=false),
            positiveY(viscousZ=false, viscousX=false),
            negativeZ(viscousX=false, viscousY=false),
            positiveZ(viscousX=false, viscousY=false)),
          H2O(
            each p_IC=defaults.p,
            negativeY(viscousZ=false, viscousX=false),
            positiveY(viscousZ=false, viscousX=false),
            negativeZ(viscousX=false, viscousY=false),
            positiveZ(viscousX=false, viscousY=false)),
          each N2(
            each p_IC=defaults.p,
            negativeY(viscousZ=false, viscousX=false),
            positiveY(viscousZ=false, viscousX=false),
            negativeZ(viscousX=false, viscousY=false),
            positiveZ(viscousX=false, viscousY=false)),
          O2(
            each p_IC=defaults.p,
            negativeY(viscousZ=false, viscousX=false),
            positiveY(viscousZ=false, viscousX=false),
            negativeZ(viscousX=false, viscousY=false),
            positiveZ(viscousX=false, viscousY=false))),
        graphite(
          each final inclC=inclC,
          each final 'incle-'='incle-',
          C(each V_IC=subregions[1].V/4),
          'e-'(
            each p_IC=defaults.p,
            negativeY(viscousZ=false, viscousX=false),
            positiveY(viscousZ=false, viscousX=false),
            negativeZ(viscousX=false, viscousY=false),
            positiveZ(viscousX=false, viscousY=false))),
        ionomer(
          each final inclC19HF37O5S=inclC19HF37O5S,
          each final 'inclH+'='inclH+',
          C19HF37O5S(each V_IC=subregions[1].V/4),
          'H+'(
            each p_IC=defaults.p,
            negativeY(viscousZ=false, viscousX=false),
            positiveY(viscousZ=false, viscousX=false),
            negativeZ(viscousX=false, viscousY=false),
            positiveZ(viscousX=false, viscousY=false))),
        each inclYMom=false,
        each inclZMom=false,
        inclYFaces=false,
        inclZFaces=false) if n_x > 0
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      FCSys.Subregions.Subregion subregion2(
        L={1,1,1}*U.cm,
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2(
            p_IC=0.95*defaults.p,
            negativeY(viscousZ=false, viscousX=false),
            positiveY(viscousZ=false, viscousX=false),
            negativeZ(viscousX=false, viscousY=false),
            positiveZ(viscousX=false, viscousY=false)),
          H2O(
            p_IC=0.95*defaults.p,
            negativeY(viscousZ=false, viscousX=false),
            positiveY(viscousZ=false, viscousX=false),
            negativeZ(viscousX=false, viscousY=false),
            positiveZ(viscousX=false, viscousY=false)),
          N2(
            p_IC=0.95*defaults.p,
            negativeY(viscousZ=false, viscousX=false),
            positiveY(viscousZ=false, viscousX=false),
            negativeZ(viscousX=false, viscousY=false),
            positiveZ(viscousX=false, viscousY=false)),
          O2(
            p_IC=0.95*defaults.p,
            negativeY(viscousZ=false, viscousX=false),
            positiveY(viscousZ=false, viscousX=false),
            negativeZ(viscousX=false, viscousY=false),
            positiveZ(viscousX=false, viscousY=false))),
        graphite(
          final inclC=inclC,
          final 'incle-'='incle-',
          C(V_IC=subregion2.V/4),
          'e-'(
            p_IC=0.95*defaults.p,
            negativeY(viscousZ=false, viscousX=false),
            positiveY(viscousZ=false, viscousX=false),
            negativeZ(viscousX=false, viscousY=false),
            positiveZ(viscousX=false, viscousY=false))),
        ionomer(
          final inclC19HF37O5S=inclC19HF37O5S,
          final 'inclH+'='inclH+',
          C19HF37O5S(V_IC=subregion2.V/4),
          'H+'(
            p_IC=0.95*defaults.p,
            negativeY(viscousZ=false, viscousX=false),
            positiveY(viscousZ=false, viscousX=false),
            negativeZ(viscousX=false, viscousY=false),
            positiveZ(viscousX=false, viscousY=false))),
        inclYMom=false,
        inclZMom=false,
        inclYFaces=false,
        inclZFaces=false)
        annotation (Placement(transformation(extent={{20,-10},{40,10}})));

      inner FCSys.BCs.Defaults defaults(analysis=true)
        annotation (Placement(transformation(extent={{70,20},{90,40}})));
      replaceable FCSys.BCs.Face.Subregion0Current bC1 constrainedby
        'e-'.FCSys.BCs.Face.Subregion(
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
            origin={-60,0})));
      replaceable FCSys.BCs.Face.Subregion0Current bC2 constrainedby
        'e-'.FCSys.BCs.Face.Subregion(
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        graphite(final inclC=inclC, final 'incle-'='incle-'),
        ionomer(final inclC19HF37O5S=inclC19HF37O5S, final 'inclH+'='inclH+'))
        annotation (__Dymola_choicesFromPackage=true,Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={60,0})));

    equation
      connect(bC1.face, subregion1.negativeX) annotation (Line(
          points={{-56,3.65701e-16},{-40,6.10623e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));
      connect(subregion1.positiveX, subregions[1].negativeX) annotation (Line(
          points={{-20,6.10623e-16},{-16,-3.36456e-22},{-16,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      for i in 1:n_x - 1 loop
        connect(subregions[i].positiveX, subregions[i + 1].negativeX)
          "Not shown on the diagram";
      end for;
      if n_x == 0 then
        connect(subregion1.positiveX, subregion2.negativeX)
          "Not shown on the diagram";
      end if;
      connect(subregions[n_x].positiveX, subregion2.negativeX) annotation (Line(
          points={{10,6.10623e-16},{20,-3.36456e-22},{20,6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(bC2.face, subregion2.positiveX) annotation (Line(
          points={{56,1.23436e-15},{50,-3.36456e-22},{50,6.10623e-16},{40,
              6.10623e-16}},
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
        experiment(StopTime=5, Tolerance=1e-06),
        experimentSetupOutput,
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionsCAndH2AndH2O.mos"));
    end SubregionsCAndH2AndH2O;

    model SubregionsCAndH2O
      "<html>Test a one-dimensional array of subregions with C and H<sub>2</sub>O</html>"
      extends SubregionsH2(
        inclC=true,
        inclH2=false,
        inclH2O=true);

      annotation (
        experiment(StopTime=20),
        experimentSetupOutput,
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionsCAndH2O.mos"));
    end SubregionsCAndH2O;

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
          StopTime=20,
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

    model Subregionseminus
      "<html>Test a one-dimensional array of subregions with e<sup>-</sup></html>"
      extends SubregionsH2(inclH2=false, 'incle-'=true);
      annotation (
        experiment(StopTime=6e-10, Tolerance=1e-06),
        experimentSetupOutput,
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.Subregionseminus.mos"));
    end Subregionseminus;

    model ThermalConduction "Test thermal conduction (through solid)"
      extends SubregionsH2(
        n_x=8,
        inclC=true,
        inclH2=false,
        subregion1(graphite(C(
              partNumInitMeth=InitMethScalar.Pressure,
              V_IC=0.99*subregion1.V,
              T_IC=1.1*defaults.T,
              T(displayUnit="degC"),
              alpha_S=Modelica.Constants.inf))),
        subregions(graphite(C(
              each partNumInitMeth=InitMethScalar.Pressure,
              each V_IC=0.99*subregions[1].V,
              each T(displayUnit="degC")))),
        subregion2(graphite(C(
              partNumInitMeth=InitMethScalar.Pressure,
              V_IC=0.99*subregion2.V,
              T(displayUnit="degC")))),
        redeclare FCSys.BCs.Face.Subregion0Current0Power bC1,
        redeclare FCSys.BCs.Face.Subregion0Current0Power bC2);
      extends Modelica.Icons.UnderConstruction;
      annotation (
        Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.ThermalConduction.mos"),

        experiment(StopTime=300, Algorithm="Dassl"),
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
              alpha_S=Modelica.Constants.inf))),
        subregions(gas(N2(each p_IC=defaults.p, phi(each displayUnit="mm/s"))),
            graphite(C(each V_IC=0.5*subregions[1].V, each T(displayUnit="degC")))),

        subregion2(gas(N2(p_IC=defaults.p, phi(displayUnit="mm/s"))), graphite(
              C(V_IC=0.5*subregion2.V, T(displayUnit="degC")))),
        redeclare FCSys.BCs.Face.Subregion0Current0Power bC1,
        redeclare FCSys.BCs.Face.Subregion0Current0Power bC2);

      extends Modelica.Icons.UnderConstruction;

      //partNumInitMeth=InitMethScalar.Pressure,
      // each partNumInitMeth=InitMethScalar.Pressure,
      // partNumInitMeth=InitMethScalar.Pressure,

      annotation (
        Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.ThermalConductionConvection.mos"),

        experiment(StopTime=200, Algorithm="Dassl"),
        experimentSetupOutput);

    end ThermalConductionConvection;

    model ReactionRamp
      "Test chemical reaction with reaction rate ramped over time"

      extends Modelica.Icons.Example;
      extends Modelica.Icons.UnderConstruction;
      FCSys.Subregions.Reaction reaction
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      FCSys.BCs.Chemical.Species.Species species1(redeclare
          FCSys.BCs.Chemical.Species.Material.PotentialElectrochemical
          materialBC) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-30,-24})));
      FCSys.BCs.Chemical.Species.Species species2(redeclare
          FCSys.BCs.Chemical.Species.Material.PotentialElectrochemical
          materialBC) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={0,-24})));
      FCSys.BCs.Chemical.Species.Species species3(redeclare
          FCSys.BCs.Chemical.Species.Material.Current materialBC, redeclare
          Modelica.Blocks.Sources.Ramp materialSpec(height=100*U.A, duration=
              3600e2)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={30,-24})));
      inner FCSys.BCs.Defaults defaults(analysis=true)
        annotation (Placement(transformation(extent={{-90,70},{-70,90}})));
    equation
      connect(species1.chemical, reaction.chemical[1]) annotation (Line(
          points={{-30,-20},{-30,-10},{5.55112e-16,-10},{5.55112e-16,
              5.55112e-16}},
          color={170,0,0},
          smooth=Smooth.None));

      connect(species2.chemical, reaction.chemical[2]) annotation (Line(
          points={{-1.11528e-16,-20},{0,-20},{0,5.55112e-16},{5.55112e-16,
              5.55112e-16}},
          color={170,0,0},
          smooth=Smooth.None));

      connect(species3.chemical, reaction.chemical[3]) annotation (Line(
          points={{30,-20},{30,-10},{5.55112e-16,-10},{5.55112e-16,5.55112e-16}},

          color={170,0,0},
          smooth=Smooth.None));

      annotation (
        Diagram(graphics),
        experiment(StopTime=36000),
        experimentSetupOutput,
        Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.ReactionRamp.mos"));
    end ReactionRamp;

    model H2_O2_H2ODynamic
      "<html>Test the 2H<sub>2</sub> + O<sub>2</sub> &#8652; 2H<sub>2</sub>O reaction dynamically</html>"

      extends Modelica.Icons.Example;
      extends Modelica.Icons.UnderConstruction;
      parameter Q.PotentialAbsolute T(displayUnit="K") = 298.15*U.K
        "Temperature";
      parameter Q.Volume V=1*U.cm^3 "Volume";

      FCSys.Subregions.Reactions.'Reaction2H2+O2<=>2H2O' reaction
        "Chemical reaction"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      SimpleSpecies species1(
        final V=V,
        final T=T,
        p_IC=0.4*U.atm,
        termDepleted=false,
        redeclare FCSys.Characteristics.H2.gas Data(b_v=[1], specVolPow={-1,0}))
        "1st species" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-30,-24})));
      SimpleSpecies species2(
        final V=V,
        final T=T,
        p_IC=1*U.atm,
        redeclare FCSys.Characteristics.O2.gas Data(b_v=[1], specVolPow={-1,0}))
        "2nd species" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={0,-24})));
      SimpleSpecies species3(
        final V=V,
        final T=T,
        p_IC=0.4*U.atm,
        redeclare FCSys.Characteristics.H2O.gas Data(b_v=[1], specVolPow={-1,0}))
        "3rd species" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={30,-24})));

      inner FCSys.BCs.Defaults defaults(analysis=true)
        annotation (Placement(transformation(extent={{-90,70},{-70,90}})));
    equation

      connect(species1.chemical, reaction.chemical[1]) annotation (Line(
          points={{-30,-20},{-30,-10},{5.55112e-16,-10},{5.55112e-16,
              5.55112e-16}},
          color={170,0,0},
          smooth=Smooth.None));

      connect(species2.chemical, reaction.chemical[2]) annotation (Line(
          points={{9.89443e-16,-20},{0,-20},{0,5.55112e-16},{5.55112e-16,
              5.55112e-16}},
          color={170,0,0},
          smooth=Smooth.None));

      connect(species3.chemical, reaction.chemical[3]) annotation (Line(
          points={{30,-20},{30,-10},{5.55112e-16,-10},{5.55112e-16,5.55112e-16}},

          color={170,0,0},
          smooth=Smooth.None));

      annotation (
        Diagram(graphics),
        experiment(
          StopTime=0.01,
          NumberOfIntervals=5000,
          Tolerance=1e-06),
        experimentSetupOutput,
        Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.H2_O2_H2ODynamic.mos"));
    end H2_O2_H2ODynamic;

    model H2O_H2ODynamic
      "<html>Test the H<sub>2</sub>O evaporation/condensation reaction dynamically</html>"

      extends Modelica.Icons.Example;
      extends Modelica.Icons.UnderConstruction;
      parameter Q.PotentialAbsolute T(displayUnit="K") = 298.15*U.K
        "Temperature";
      parameter Q.Volume V=1*U.cm^3 "Volume";

      FCSys.Subregions.Reactions.'ReactionH2O<=>H2O' reaction
        "Chemical reaction"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      SimpleSpecies species1(
        final V=V,
        final T=T,
        redeclare FCSys.Characteristics.H2O.liquid Data) "1st species"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-30,-24})));
      SimpleSpecies species2(
        final V=V,
        final T=T,
        p_IC=0.4*U.atm,
        redeclare FCSys.Characteristics.H2O.gas Data(b_v=[1], specVolPow={-1,0}))
        "2nd species" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={0,-24})));

      inner FCSys.BCs.Defaults defaults(analysis=true)
        annotation (Placement(transformation(extent={{-90,70},{-70,90}})));
    equation

      connect(species1.chemical, reaction.chemical[1]) annotation (Line(
          points={{-30,-20},{-30,-10},{5.55112e-16,-10},{5.55112e-16,
              5.55112e-16}},
          color={170,0,0},
          smooth=Smooth.None));

      connect(species2.chemical, reaction.chemical[2]) annotation (Line(
          points={{9.89443e-16,-20},{0,-20},{0,5.55112e-16},{5.55112e-16,
              5.55112e-16}},
          color={170,0,0},
          smooth=Smooth.None));

      annotation (
        Documentation(info="<html>Note that the pressure of liquid H<sub>2</sub>O shows as zero.
  Since it is assumed to be incompressible, its pressure cannot be determined by the specific volume and temperature.
  In that case, the <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.p_vT\">p_vT</a> function
  returns zero.</html>"),
        Diagram(graphics),
        experiment(
          StopTime=2,
          NumberOfIntervals=5000,
          Tolerance=1e-06),
        experimentSetupOutput,
        Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.H2O_H2ODynamic.mos"));
    end H2O_H2ODynamic;

    model SimpleSpecies "Simple species model to test reactions"
      //extends FCSys.BaseClasses.Icons.Names.Middle;
      extends Modelica.Icons.UnderConstruction;
      replaceable FCSys.Characteristics.BaseClasses.Characteristic Data
        "Characteristic data of the species"
        annotation (Dialog(group="Material properties"));
      parameter Q.Volume V=1*U.cm^3 "Volume";
      parameter Q.PotentialAbsolute T(displayUnit="K") = 298.15*U.K
        "Temperature";
      parameter Q.PressureAbsolute p_IC=1*U.atm
        "<html>Initial pressure (<i>p</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization", group="Scalar properties"));
      parameter Q.Resistance alphaPerLstar=1/U.A
        "<html>Material exchange resistance (&alpha;/<i>L</i><sup>&#9733;</sup>)</html>";
      parameter Boolean termDepleted=true "Terminate when depleted" annotation
        (
        Evaluate=true,
        Dialog(tab="Assumptions"),
        choices(__Dymola_checkBox=true));

      Q.ParticleNumber N(nominal=1*U.C, start=1*U.C) "Particle number";
      Q.VolumeSpecific v(nominal=4*U.C/U.cm^3, start=4*U.C/U.cm^3)
        "Specific volume";
      Q.Pressure p(
        nominal=1*U.atm,
        start=p_IC,
        fixed=true) "Pressure";
      Q.Potential h(nominal=1*U.V) "Specific enthalpy";
      Q.Potential mu(nominal=1*U.V) "Electrochemical potential";

      FCSys.Connectors.ChemicalOutput chemical annotation (Placement(
            transformation(extent={{-10,-50},{10,-30}}), iconTransformation(
              extent={{-10,-50},{10,-30}})));

    protected
      outer FCSys.BCs.Defaults defaults "Default settings" annotation (
          Placement(transformation(extent={{-10,-10},{10,10}}),
            iconTransformation(extent={{-10,90},{10,110}})));
      Boolean depleted "true, if pressure is at or below the minimum";

    equation
      // Properties
      v*N = V;
      p = Data.p_vT(v, T);
      h = Data.h_pT(
            p,
            T,
            referenceEnthalpy=ReferenceEnthalpy.ZeroAt25degC);
      //    referenceEnthalpy=ReferenceEnthalpy.EnthalpyOfFormationAt25degC);
      mu = h - T*Data.s0_T(T);
      depleted = p <= Data.p_min;
      if termDepleted then
        when depleted then
          terminate("The " + Data.name + " species is depleted (N = " + String(
            N/U.C) + " C and p = " + String(p/U.Pa) + " Pa).");
        end when;
      end if;

      // Material exchange
      T*alphaPerLstar*chemical.Ndot = if depleted and not termDepleted then 0
         else 2*(chemical.mu - mu - Data.Deltah0_f);

      // Conservation of material
      der(N)/U.s = chemical.Ndot;
      annotation (defaultComponentName="species", Icon(graphics={Rectangle(
                  extent={{-100,40},{100,-40}},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Line(
                  points={{-100,-40},{100,-40}},
                  color={0,0,0},
                  smooth=Smooth.None,
                  pattern=LinePattern.Dash),Line(
                  points={{-100,-40},{-100,40},{100,40},{100,-40}},
                  pattern=LinePattern.None,
                  smooth=Smooth.None),Text(
                  extent={{-100,-20},{100,20}},
                  textString="%name",
                  lineColor={0,0,0})}));
    end SimpleSpecies;

    model Cell "Test both half reactions of a cell"
      extends Modelica.Icons.UnderConstruction;
      extends Modelica.Icons.Example;
      output Q.Number DeltamuPerT='e-_an'.chemical.mu - 'e-_ca'.chemical.mu
        "Voltage of cathode w.r.t. anode";
      // The negative factor is due to the charge negative charge of electrons.

      FCSys.Subregions.Reaction HOR(nu={2,2,1}) "Hydrogen oxidation reaction"
        annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));

      FCSys.BCs.Chemical.Species.Species 'e-_an'(redeclare
          FCSys.BCs.Chemical.Species.Material.PotentialElectrochemical
          materialBC) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-52,-24})));
      FCSys.BCs.Chemical.Species.Species H2(redeclare
          FCSys.BCs.Chemical.Species.Material.PotentialElectrochemical
          materialBC) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-28,-24})));
      inner FCSys.BCs.Defaults defaults(analysis=true)
        annotation (Placement(transformation(extent={{-90,70},{-70,90}})));
      FCSys.Subregions.Reaction ORR(nu={4,4,1,-2}) "Oxygen reduction reaction"
        annotation (Placement(transformation(extent={{30,-10},{50,10}})));
      FCSys.BCs.Chemical.Species.Species 'e-_ca'(redeclare
          FCSys.BCs.Chemical.Species.Material.Current materialBC, redeclare
          Modelica.Blocks.Sources.Ramp materialSpec(duration=3600e2, height=100
              *U.A)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={14,-24})));
      FCSys.BCs.Chemical.Species.Species O2(redeclare
          FCSys.BCs.Chemical.Species.Material.PotentialElectrochemical
          materialBC) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={40,-24})));
      FCSys.BCs.Chemical.Species.Species H2O(redeclare
          FCSys.BCs.Chemical.Species.Material.PotentialElectrochemical
          materialBC, redeclare Modelica.Blocks.Sources.Constant materialSpec(k
            =-2*1.20646*U.V/(298.15*U.K))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={64,-24})));
      // Note:  The OCV of H2/O2 cell is 1.20646 V at 300 K, assuming the product is
      // gaseous.
    equation
      connect('e-_an'.chemical, HOR.chemical[1]) annotation (Line(
          points={{-52,-20},{-52,-10},{-40,-10},{-40,5.55112e-16}},
          color={170,0,0},
          smooth=Smooth.None));
      connect(H2.chemical, HOR.chemical[3]) annotation (Line(
          points={{-28,-20},{-28,-10},{-40,-10},{-40,5.55112e-16}},
          color={170,0,0},
          smooth=Smooth.None));
      connect('e-_ca'.chemical, ORR.chemical[1]) annotation (Line(
          points={{14,-20},{14,-10},{40,-10},{40,5.55112e-16}},
          color={170,0,0},
          smooth=Smooth.None));
      connect(O2.chemical, ORR.chemical[3]) annotation (Line(
          points={{40,-20},{40,5.55112e-16}},
          color={170,0,0},
          smooth=Smooth.None));
      connect(H2O.chemical, ORR.chemical[4]) annotation (Line(
          points={{64,-20},{64,-10},{40,-10},{40,5.55112e-16}},
          color={170,0,0},
          smooth=Smooth.None));
      connect(HOR.chemical[2], ORR.chemical[2]) annotation (Line(
          points={{-40,5.55112e-16},{8,-4.87687e-22},{8,5.55112e-16},{40,
              5.55112e-16}},
          color={170,0,0},
          smooth=Smooth.None));

      annotation (
        Diagram(graphics),
        experiment(StopTime=360000),
        experimentSetupOutput,
        Commands(file="resources/scripts/Dymola/Subregions.Examples.Cell.mos"));
    end Cell;

    model Reaction "Test an electrochemical reaction"
      extends Modelica.Icons.Example;

      parameter Integer n_lin(
        final min=1,
        final max=3) = 1
        "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>"
        annotation (HideResult=true);

      FCSys.Subregions.Reaction reaction(final n_lin=n_lin,n_spec=3)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      FCSys.BCs.Chemical.Species.Species 'e-'(
        materialBC=FCSys.BCs.Chemical.Species.BaseClasses.BCTypeMaterial.PotentialElectrochemical,

        formula="e-",
        final n_lin=n_lin) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-30,-24})));

      FCSys.BCs.Chemical.Species.Species 'H+'(
        materialBC=FCSys.BCs.Chemical.Species.BaseClasses.BCTypeMaterial.PotentialElectrochemical,

        formula="H+",
        final n_lin=n_lin) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={0,-24})));

      FCSys.BCs.Chemical.Species.Species H2(
        materialBC=FCSys.BCs.Chemical.Species.BaseClasses.BCTypeMaterial.Current,

        formula="H2",
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

    model SubregionH22 "Test a subregion"
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
      inner FCSys.BCs.Defaults defaults(analysis=false)
        annotation (Placement(transformation(extent={{30,20},{50,40}})));

      FCSys.Subregions.Subregion subregion(
        L={1,1,1}*U.cm,
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        graphite(
          final inclC=inclC,
          final 'incle-'='incle-',
          C(V_IC=subregion.V/4)),
        ionomer(
          final inclC19HF37O5S=inclC19HF37O5S,
          final 'inclH+'='inclH+',
          C19HF37O5S(V_IC=subregion.V/4)),
        inclYMom=false,
        inclZMom=false,
        inclYFaces=false,
        inclZFaces=false)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      annotation (
        Placement(transformation(extent={{70,70},{90,90}})),
        Diagram(graphics),
        experiment(StopTime=10),
        experimentSetupOutput,
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionH2.mos"));
    end SubregionH22;

    model Species "Test a species"
      import FCSys;
      extends Modelica.Icons.Example;
      // Geometric parameters
      inner parameter Q.Length L[3](each final min=Modelica.Constants.small) =
        ones(3)*U.m "<html>Length (<b>L</b>)</html>"
        annotation (Dialog(group="Geometry"));
      final inner parameter Q.Area A[3]={L[cartWrap(ax + 1)]*L[cartWrap(ax + 2)]
          for ax in 1:3} "Cross-sectional area";
      final parameter Q.Volume V=product(L) "Volume";

      // Species
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
      inner FCSys.BCs.Defaults defaults(analysis=true)
        annotation (Placement(transformation(extent={{30,20},{50,40}})));

      replaceable FCSys.Subregions.Species.H2.gas.Fixed H2
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      FCSys.BCs.InertDalton.Species speciesInertBC annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=225,
            origin={20,-20})));
    equation
      connect(speciesInertBC.inert, H2.inert) annotation (Line(
          points={{17.1716,-17.1716},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      annotation (
        Placement(transformation(extent={{70,70},{90,90}})),
        Diagram(graphics),
        experiment(StopTime=10),
        experimentSetupOutput,
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionH2.mos"));
    end Species;
  end Examples;

  model Subregion "Subregion with all phases included"

    parameter Boolean inclReact=true "Include reaction(s), as appropriate"
      annotation (
      Evaluate=true,
      HideResult=true,
      Dialog(tab="Assumptions"),
      choices(__Dymola_checkBox=true));
    // Note:  This is listed above the extends clause so that it is
    // listed first in the parameter dialog.
    extends BaseClasses.PartialSubregion;

    replaceable FCSys.Subregions.Phases.Gas gas(inclH2O=true, final inclLin={
          inclXMom,inclYMom,inclZMom}) "Gas" annotation (Dialog(group="Phases"),
        Placement(transformation(extent={{-10,-10},{10,10}})));

    // **Currently, both reactions must be included at this level.

    replaceable FCSys.Subregions.Phases.Graphite graphite(final inclLin={
          inclXMom,inclYMom,inclZMom}) "Graphite" annotation (Dialog(group=
            "Phases"), Placement(transformation(extent={{-10,-10},{10,10}})));

    replaceable FCSys.Subregions.Phases.Ionomer ionomer(final inclLin={inclXMom,
          inclYMom,inclZMom}) "Ionomer" annotation (Dialog(group="Phases"),
        Placement(transformation(extent={{-10,-10},{10,10}})));

    /*
  replaceable FCSys.Subregions.Phases.Liquid liquid(final inclLin={inclXMom,
        inclYMom,inclZMom})  "Liquid" annotation (

    Dialog(group="Phases"),
    Placement(transformation(extent={{-10,-10},{10,10}})));
  */

    FCSys.Subregions.Reaction HOR(final n_lin=n_lin,n_spec=3) if inclReact and
      (graphite.'incle-' and ionomer.'inclH+' and gas.inclH2 and not (gas.inclO2
       and gas.inclH2O)) "Hydrogen oxidation reaction"
      annotation (Placement(transformation(extent={{-50,30},{-30,50}})));
    FCSys.Subregions.Reaction ORR(final n_lin=n_lin,n_spec=4) if inclReact and
      (graphite.'incle-' and ionomer.'inclH+' and gas.inclO2 and gas.inclH2O
       and not gas.inclH2) "Oxygen reduction reaction"
      annotation (Placement(transformation(extent={{-50,30},{-30,50}})));

    /* **Note: Multiple reactions cannot be included at once due to the following error:
Error: Failed to expand the variable HOR.chemical[1].mphi
Error: Failed to expand the variable HOR.chemical[2].mphi
Error: Failed to expand the variable ORR.chemical[1].mphi
Error: Failed to expand the variable ORR.chemical[2].mphi
  */

  protected
    Connectors.ChemicalBusInternal chemical
      "Internal connector to route chemical reaction"
      annotation (Placement(transformation(extent={{-30,10},{-10,30}})));
  equation
    // Chemical interactions
    connect(HOR.chemical[1], chemical.'e-') annotation (Line(
        points={{-40,39.3333},{-20,20}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(HOR.chemical[2], chemical.'H+') annotation (Line(
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
    connect(gas.negativeX, negativeX.gas) annotation (Line(
        points={{-8,6.10623e-16},{-8,5.55112e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.positiveX, positiveX.gas) annotation (Line(
        points={{8,6.10623e-16},{8,5.55112e-16},{40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.negativeY, negativeY.gas) annotation (Line(
        points={{6.10623e-16,-8.4},{6.10623e-16,-40},{5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.positiveY, positiveY.gas) annotation (Line(
        points={{6.10623e-16,10},{0,10},{0,10},{-4.87687e-22,40},{-4.87687e-22,
            40},{5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.negativeZ, negativeZ.gas) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.positiveZ, positiveZ.gas) annotation (Line(
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
    connect(graphite.negativeX, negativeX.graphite) annotation (Line(
        points={{-8,6.10623e-16},{-8,5.55112e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.positiveX, positiveX.graphite) annotation (Line(
        points={{8,6.10623e-16},{8,5.55112e-16},{40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.negativeY, negativeY.graphite) annotation (Line(
        points={{6.10623e-16,-8.4},{6.10623e-16,-40},{5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.positiveY, positiveY.graphite) annotation (Line(
        points={{6.10623e-16,10},{0,10},{0,10},{-4.87687e-22,40},{-4.87687e-22,
            40},{5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.negativeZ, negativeZ.graphite) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.positiveZ, positiveZ.graphite) annotation (Line(
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
    connect(ionomer.negativeX, negativeX.ionomer) annotation (Line(
        points={{-8,6.10623e-16},{-8,5.55112e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.positiveX, positiveX.ionomer) annotation (Line(
        points={{8,6.10623e-16},{8,5.55112e-16},{40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.negativeY, negativeY.ionomer) annotation (Line(
        points={{6.10623e-16,-8.4},{6.10623e-16,-40},{5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.positiveY, positiveY.ionomer) annotation (Line(
        points={{6.10623e-16,10},{0,10},{0,10},{-4.87687e-22,40},{-4.87687e-22,
            40},{5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.negativeZ, negativeZ.ionomer) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.positiveZ, positiveZ.ionomer) annotation (Line(
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
  connect(liquid.negativeX, negativeX.liquid) annotation (Line(
      points={{-10,6.10623e-16},{-10,1.16573e-15},{-25,1.16573e-15},{-25,5.55112e-16},
          {-40,5.55112e-16}},
      color={127,127,127},
      pattern=LinePattern.None,
      thickness=0.5,
      smooth=Smooth.None));
  connect(liquid.positiveX, positiveX.liquid) annotation (Line(
      points={{10,6.10623e-16},{10,5.55112e-16},{40,5.55112e-16}},
      color={127,127,127},
      pattern=LinePattern.None,
      thickness=0.5,
      smooth=Smooth.None));
  connect(liquid.negativeY, negativeY.liquid) annotation (Line(
      points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
      color={127,127,127},
      pattern=LinePattern.None,
      thickness=0.5,
      smooth=Smooth.None));
  connect(liquid.positiveY, positiveY.liquid) annotation (Line(
      points={{6.10623e-16,10},{-4.87687e-22,10},{-4.87687e-22,40},{5.55112e-16,
          40}},
      color={127,127,127},
      pattern=LinePattern.None,
      thickness=0.5,
      smooth=Smooth.None));
  connect(liquid.negativeZ, negativeZ.liquid) annotation (Line(
      points={{5,5},{20,20}},
      color={127,127,127},
      pattern=LinePattern.None,
      thickness=0.5,
      smooth=Smooth.None));
  connect(liquid.positiveZ, positiveZ.liquid) annotation (Line(
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
    // Note:  This is listed above the extension clause so that it is
    // listed first in the parameter dialog.
    extends BaseClasses.PartialSubregion;

    replaceable FCSys.Subregions.Phases.Gas gas(
      inclH2O=true,
      inclReact=inclReact,
      final inclLin={inclXMom,inclYMom,inclZMom}) "Gas" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(group="Phases"),
      Placement(transformation(extent={{-10,-10},{10,10}})));

    replaceable FCSys.Subregions.Phases.Ionomer ionomer(final inclLin={inclXMom,
          inclYMom,inclZMom}) "Ionomer" annotation (Dialog(group="Phases"),
        Placement(transformation(extent={{-10,-10},{10,10}})));

  equation
    // Gas
    connect(gas.inert, volume.inert) annotation (Line(
        points={{8,-8},{11,-11}},
        color={72,90,180},
        smooth=Smooth.None));
    connect(gas.negativeX, negativeX.gas) annotation (Line(
        points={{-8,6.10623e-16},{-8,1.16573e-15},{-25,1.16573e-15},{-25,
            5.55112e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.positiveX, positiveX.gas) annotation (Line(
        points={{8,6.10623e-16},{8,5.55112e-16},{40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.negativeY, negativeY.gas) annotation (Line(
        points={{6.10623e-16,-8.4},{6.10623e-16,-40},{5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.positiveY, positiveY.gas) annotation (Line(
        points={{6.10623e-16,10},{-4.87687e-22,10},{-4.87687e-22,40},{
            5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.negativeZ, negativeZ.gas) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.positiveZ, positiveZ.gas) annotation (Line(
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
    connect(ionomer.negativeX, negativeX.ionomer) annotation (Line(
        points={{-8,6.10623e-16},{-8,1.16573e-15},{-25,1.16573e-15},{-25,
            5.55112e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.positiveX, positiveX.ionomer) annotation (Line(
        points={{8,6.10623e-16},{8,5.55112e-16},{40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.negativeY, negativeY.ionomer) annotation (Line(
        points={{6.10623e-16,-8.4},{6.10623e-16,-40},{5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.positiveY, positiveY.ionomer) annotation (Line(
        points={{6.10623e-16,10},{-4.87687e-22,10},{-4.87687e-22,40},{
            5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.negativeZ, negativeZ.ionomer) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.positiveZ, positiveZ.ionomer) annotation (Line(
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
    // Note:  This is listed above the extends clause so that it is
    // listed first in the parameter dialog.
    extends BaseClasses.PartialSubregion;

    replaceable FCSys.Subregions.Phases.Gas gas(
      inclH2O=true,
      inclReact=inclReact,
      final inclLin={inclXMom,inclYMom,inclZMom}) "Gas" annotation (Dialog(
          group="Phases"), Placement(transformation(extent={{-10,-10},{10,10}})));

    replaceable FCSys.Subregions.Phases.Graphite graphite(final inclLin={
          inclXMom,inclYMom,inclZMom}) "Graphite" annotation (Dialog(group=
            "Phases"), Placement(transformation(extent={{-10,-10},{10,10}})));

  equation
    // Gas
    connect(gas.inert, volume.inert) annotation (Line(
        points={{8,-8},{11,-11}},
        color={72,90,180},
        smooth=Smooth.None));
    connect(gas.negativeX, negativeX.gas) annotation (Line(
        points={{-8,6.10623e-16},{-8,1.16573e-15},{-25,1.16573e-15},{-25,
            5.55112e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.positiveX, positiveX.gas) annotation (Line(
        points={{8,6.10623e-16},{8,5.55112e-16},{40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.negativeY, negativeY.gas) annotation (Line(
        points={{6.10623e-16,-8.4},{6.10623e-16,-40},{5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.positiveY, positiveY.gas) annotation (Line(
        points={{6.10623e-16,10},{-4.87687e-22,10},{-4.87687e-22,40},{
            5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.negativeZ, negativeZ.gas) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.positiveZ, positiveZ.gas) annotation (Line(
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
    connect(graphite.negativeX, negativeX.graphite) annotation (Line(
        points={{-8,6.10623e-16},{-8,1.16573e-15},{-25,1.16573e-15},{-25,
            5.55112e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.positiveX, positiveX.graphite) annotation (Line(
        points={{8,6.10623e-16},{8,5.55112e-16},{40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.negativeY, negativeY.graphite) annotation (Line(
        points={{6.10623e-16,-8.4},{6.10623e-16,-40},{5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.positiveY, positiveY.graphite) annotation (Line(
        points={{6.10623e-16,10},{-4.87687e-22,10},{-4.87687e-22,40},{
            5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.negativeZ, negativeZ.graphite) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.positiveZ, positiveZ.graphite) annotation (Line(
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
      extends BaseClasses.NullPhase(final hasSpecies=inclC or inclC19HF37O5S
             or 'incle-' or inclH2 or inclH2O or 'inclH+' or inclN2 or inclO2);

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
      replaceable FCSys.Subregions.Species.C.graphite.Fixed C(final k=k, final
          inclLin=inclLin) if inclC constrainedby
        FCSys.Subregions.Species.Species "C model" annotation (
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
      replaceable FCSys.Subregions.Species.C19HF37O5S.solid.Fixed C19HF37O5S(
          final k=k, final inclLin=inclLin) if inclC19HF37O5S constrainedby
        FCSys.Subregions.Species.Species
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
      replaceable FCSys.Subregions.Species.'e-'.graphite.Fixed 'e-'(final k=k,
          final inclLin=inclLin) if 'incle-' constrainedby
        FCSys.Subregions.Species.Species "<html>'e-' model</html>" annotation (
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
      replaceable FCSys.Subregions.Species.H2.gas.Fixed H2(final k=k, final
          inclLin=inclLin) if inclH2 constrainedby
        FCSys.Subregions.Species.Species "<html>H<sub>2</sub> model</html>"
        annotation (
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
      replaceable FCSys.Subregions.Species.H2O.gas.Fixed H2O(final k=k, final
          inclLin=inclLin) if inclH2O constrainedby
        FCSys.Subregions.Species.Species "<html>H<sub>2</sub>O model</html>"
        annotation (
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
      replaceable FCSys.Subregions.Species.'H+'.solid.Fixed 'H+'(final k=k,
          final inclLin=inclLin) if 'inclH+' constrainedby
        FCSys.Subregions.Species.Species "<html>H<sup>+</sup> model</html>"
        annotation (
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
      replaceable FCSys.Subregions.Species.N2.gas.Fixed N2(final k=k, final
          inclLin=inclLin) if inclN2 constrainedby
        FCSys.Subregions.Species.Species "<html>N<sub>2</sub> model</html>"
        annotation (
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
      replaceable FCSys.Subregions.Species.O2.gas.Fixed O2(final k=k, final
          inclLin=inclLin) if inclO2 constrainedby
        FCSys.Subregions.Species.Species "<html>O<sub>2</sub> model</html>"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclO2),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      FCSys.Subregions.Reaction '2e-+2H+=H2'(final n_lin=n_lin) if inclReact
         and (('incle-' and 'inclH+') or ('incle-' and (inclH2 or (inclO2 and
        inclH2O))) or ('inclH+' and (inclH2 or (inclO2 and inclH2O))))
        "Hydrogen oxidation reaction"
        annotation (Placement(transformation(extent={{-50,30},{-30,50}})));
      FCSys.Subregions.Reaction '2H2+O2=2H2O'(final n_lin=n_lin) if inclReact
         and ((inclH2 and inclO2) or (inclH2 and inclH2O) or (inclO2 and
        inclH2O))
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
      // Exchange
      connect(C.chemical, chemical.C) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(C.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(C.negativeX.material, negativeX.C.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.positiveX.material, positiveX.C.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.negativeY.material, negativeY.C.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.positiveY.material, positiveY.C.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.negativeZ.material, negativeZ.C.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.positiveZ.material, positiveZ.C.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(C.negativeX.momentumY, negativeX.C.momentumY) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.positiveX.momentumY, positiveX.C.momentumY) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.negativeY.momentumZ, negativeY.C.momentumZ) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.positiveY.momentumZ, positiveY.C.momentumZ) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.negativeZ.momentumX, negativeZ.C.momentumX) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.positiveZ.momentumX, positiveZ.C.momentumX) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(C.negativeX.momentumZ, negativeX.C.momentumZ) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.positiveX.momentumZ, positiveX.C.momentumZ) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.negativeY.momentumX, negativeY.C.momentumX) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.positiveY.momentumX, positiveY.C.momentumX) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.negativeZ.momentumY, negativeZ.C.momentumY) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.positiveZ.momentumY, positiveZ.C.momentumY) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(C.negativeX.entropy, negativeX.C.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.positiveX.entropy, positiveX.C.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.negativeY.entropy, negativeY.C.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.positiveY.entropy, positiveY.C.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.negativeZ.entropy, negativeZ.C.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.positiveZ.entropy, positiveZ.C.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // C19HF37O5S
      // Exchange
      connect(C19HF37O5S.chemical, chemical.C19HF37O5S) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(C19HF37O5S.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(C19HF37O5S.negativeX.material, negativeX.C19HF37O5S.material)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.positiveX.material, positiveX.C19HF37O5S.material)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.negativeY.material, negativeY.C19HF37O5S.material)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.positiveY.material, positiveY.C19HF37O5S.material)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.negativeZ.material, negativeZ.C19HF37O5S.material)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.positiveZ.material, positiveZ.C19HF37O5S.material)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(C19HF37O5S.negativeX.momentumY, negativeX.C19HF37O5S.momentumY)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.positiveX.momentumY, positiveX.C19HF37O5S.momentumY)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.negativeY.momentumZ, negativeY.C19HF37O5S.momentumZ)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.positiveY.momentumZ, positiveY.C19HF37O5S.momentumZ)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.negativeZ.momentumX, negativeZ.C19HF37O5S.momentumX)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.positiveZ.momentumX, positiveZ.C19HF37O5S.momentumX)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(C19HF37O5S.negativeX.momentumZ, negativeX.C19HF37O5S.momentumZ)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.positiveX.momentumZ, positiveX.C19HF37O5S.momentumZ)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.negativeY.momentumX, negativeY.C19HF37O5S.momentumX)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.positiveY.momentumX, positiveY.C19HF37O5S.momentumX)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.negativeZ.momentumY, negativeZ.C19HF37O5S.momentumY)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.positiveZ.momentumY, positiveZ.C19HF37O5S.momentumY)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(C19HF37O5S.negativeX.entropy, negativeX.C19HF37O5S.entropy)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.positiveX.entropy, positiveX.C19HF37O5S.entropy)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.negativeY.entropy, negativeY.C19HF37O5S.entropy)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.positiveY.entropy, positiveY.C19HF37O5S.entropy)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.negativeZ.entropy, negativeZ.C19HF37O5S.entropy)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.positiveZ.entropy, positiveZ.C19HF37O5S.entropy)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // e-
      // Exchange
      connect('e-'.chemical, chemical.'e-') annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect('e-'.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect('e-'.negativeX.material, negativeX.'e-'.material) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.positiveX.material, positiveX.'e-'.material) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.negativeY.material, negativeY.'e-'.material) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.positiveY.material, positiveY.'e-'.material) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.negativeZ.material, negativeZ.'e-'.material) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.positiveZ.material, positiveZ.'e-'.material) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect('e-'.negativeX.momentumY, negativeX.'e-'.momentumY) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.positiveX.momentumY, positiveX.'e-'.momentumY) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.negativeY.momentumZ, negativeY.'e-'.momentumZ) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.positiveY.momentumZ, positiveY.'e-'.momentumZ) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.negativeZ.momentumX, negativeZ.'e-'.momentumX) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.positiveZ.momentumX, positiveZ.'e-'.momentumX) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect('e-'.negativeX.momentumZ, negativeX.'e-'.momentumZ) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.positiveX.momentumZ, positiveX.'e-'.momentumZ) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.negativeY.momentumX, negativeY.'e-'.momentumX) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.positiveY.momentumX, positiveY.'e-'.momentumX) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.negativeZ.momentumY, negativeZ.'e-'.momentumY) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.positiveZ.momentumY, positiveZ.'e-'.momentumY) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect('e-'.negativeX.entropy, negativeX.'e-'.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.positiveX.entropy, positiveX.'e-'.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.negativeY.entropy, negativeY.'e-'.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.positiveY.entropy, positiveY.'e-'.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.negativeZ.entropy, negativeZ.'e-'.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.positiveZ.entropy, positiveZ.'e-'.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // H2
      // Exchange
      connect(H2.chemical, chemical.H2) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(H2.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(H2.negativeX.material, negativeX.H2.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.positiveX.material, positiveX.H2.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.negativeY.material, negativeY.H2.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.positiveY.material, positiveY.H2.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.negativeZ.material, negativeZ.H2.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.positiveZ.material, positiveZ.H2.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(H2.negativeX.momentumY, negativeX.H2.momentumY) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.positiveX.momentumY, positiveX.H2.momentumY) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.negativeY.momentumZ, negativeY.H2.momentumZ) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.positiveY.momentumZ, positiveY.H2.momentumZ) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.negativeZ.momentumX, negativeZ.H2.momentumX) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.positiveZ.momentumX, positiveZ.H2.momentumX) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(H2.negativeX.momentumZ, negativeX.H2.momentumZ) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.positiveX.momentumZ, positiveX.H2.momentumZ) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.negativeY.momentumX, negativeY.H2.momentumX) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.positiveY.momentumX, positiveY.H2.momentumX) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.negativeZ.momentumY, negativeZ.H2.momentumY) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.positiveZ.momentumY, positiveZ.H2.momentumY) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(H2.negativeX.entropy, negativeX.H2.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.positiveX.entropy, positiveX.H2.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.negativeY.entropy, negativeY.H2.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.positiveY.entropy, positiveY.H2.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.negativeZ.entropy, negativeZ.H2.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.positiveZ.entropy, positiveZ.H2.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // H2O
      // Exchange
      connect(H2O.chemical, chemical.H2O) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(H2O.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(H2O.negativeX.material, negativeX.H2O.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveX.material, positiveX.H2O.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeY.material, negativeY.H2O.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveY.material, positiveY.H2O.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeZ.material, negativeZ.H2O.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.positiveZ.material, positiveZ.H2O.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(H2O.negativeX.momentumY, negativeX.H2O.momentumY) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveX.momentumY, positiveX.H2O.momentumY) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeY.momentumZ, negativeY.H2O.momentumZ) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveY.momentumZ, positiveY.H2O.momentumZ) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeZ.momentumX, negativeZ.H2O.momentumX) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.positiveZ.momentumX, positiveZ.H2O.momentumX) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(H2O.negativeX.momentumZ, negativeX.H2O.momentumZ) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveX.momentumZ, positiveX.H2O.momentumZ) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeY.momentumX, negativeY.H2O.momentumX) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveY.momentumX, positiveY.H2O.momentumX) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeZ.momentumY, negativeZ.H2O.momentumY) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.positiveZ.momentumY, positiveZ.H2O.momentumY) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(H2O.negativeX.entropy, negativeX.H2O.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveX.entropy, positiveX.H2O.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeY.entropy, negativeY.H2O.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveY.entropy, positiveY.H2O.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeZ.entropy, negativeZ.H2O.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.positiveZ.entropy, positiveZ.H2O.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // 'H+'
      // Exchange
      connect('H+'.chemical, chemical.'H+') annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect('H+'.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect('H+'.negativeX.material, negativeX.'H+'.material) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.positiveX.material, positiveX.'H+'.material) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.negativeY.material, negativeY.'H+'.material) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.positiveY.material, positiveY.'H+'.material) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.negativeZ.material, negativeZ.'H+'.material) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.positiveZ.material, positiveZ.'H+'.material) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect('H+'.negativeX.momentumY, negativeX.'H+'.momentumY) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.positiveX.momentumY, positiveX.'H+'.momentumY) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.negativeY.momentumZ, negativeY.'H+'.momentumZ) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.positiveY.momentumZ, positiveY.'H+'.momentumZ) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.negativeZ.momentumX, negativeZ.'H+'.momentumX) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.positiveZ.momentumX, positiveZ.'H+'.momentumX) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect('H+'.negativeX.momentumZ, negativeX.'H+'.momentumZ) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.positiveX.momentumZ, positiveX.'H+'.momentumZ) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.negativeY.momentumX, negativeY.'H+'.momentumX) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.positiveY.momentumX, positiveY.'H+'.momentumX) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.negativeZ.momentumY, negativeZ.'H+'.momentumY) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.positiveZ.momentumY, positiveZ.'H+'.momentumY) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect('H+'.negativeX.entropy, negativeX.'H+'.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.positiveX.entropy, positiveX.'H+'.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.negativeY.entropy, negativeY.'H+'.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.positiveY.entropy, positiveY.'H+'.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.negativeZ.entropy, negativeZ.'H+'.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.positiveZ.entropy, positiveZ.'H+'.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // N2
      // Exchange
      connect(N2.chemical, chemical.N2) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(N2.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(N2.negativeX.material, negativeX.N2.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.positiveX.material, positiveX.N2.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.negativeY.material, negativeY.N2.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.positiveY.material, positiveY.N2.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.negativeZ.material, negativeZ.N2.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.positiveZ.material, positiveZ.N2.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(N2.negativeX.momentumY, negativeX.N2.momentumY) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.positiveX.momentumY, positiveX.N2.momentumY) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.negativeY.momentumZ, negativeY.N2.momentumZ) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.positiveY.momentumZ, positiveY.N2.momentumZ) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.negativeZ.momentumX, negativeZ.N2.momentumX) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.positiveZ.momentumX, positiveZ.N2.momentumX) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(N2.negativeX.momentumZ, negativeX.N2.momentumZ) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.positiveX.momentumZ, positiveX.N2.momentumZ) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.negativeY.momentumX, negativeY.N2.momentumX) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.positiveY.momentumX, positiveY.N2.momentumX) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.negativeZ.momentumY, negativeZ.N2.momentumY) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.positiveZ.momentumY, positiveZ.N2.momentumY) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(N2.negativeX.entropy, negativeX.N2.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.positiveX.entropy, positiveX.N2.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.negativeY.entropy, negativeY.N2.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.positiveY.entropy, positiveY.N2.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.negativeZ.entropy, negativeZ.N2.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.positiveZ.entropy, positiveZ.N2.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // O2
      // Exchange
      connect(O2.chemical, chemical.O2) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(O2.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(O2.negativeX.material, negativeX.O2.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.positiveX.material, positiveX.O2.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.negativeY.material, negativeY.O2.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.positiveY.material, positiveY.O2.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.negativeZ.material, negativeZ.O2.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.positiveZ.material, positiveZ.O2.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(O2.negativeX.momentumY, negativeX.O2.momentumY) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.positiveX.momentumY, positiveX.O2.momentumY) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.negativeY.momentumZ, negativeY.O2.momentumZ) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.positiveY.momentumZ, positiveY.O2.momentumZ) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.negativeZ.momentumX, negativeZ.O2.momentumX) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.positiveZ.momentumX, positiveZ.O2.momentumX) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(O2.negativeX.momentumZ, negativeX.O2.momentumZ) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.positiveX.momentumZ, positiveX.O2.momentumZ) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.negativeY.momentumX, negativeY.O2.momentumX) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.positiveY.momentumX, positiveY.O2.momentumX) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.negativeZ.momentumY, negativeZ.O2.momentumY) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.positiveZ.momentumY, positiveZ.O2.momentumY) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(O2.negativeX.entropy, negativeX.O2.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.positiveX.entropy, positiveX.O2.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.negativeY.entropy, negativeY.O2.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.positiveY.entropy, positiveY.O2.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.negativeZ.entropy, negativeZ.O2.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.positiveZ.entropy, positiveZ.O2.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      annotation (
        defaultComponentPrefixes="replaceable",
        Documentation(info="<html>
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

      extends BaseClasses.NullPhase(final hasSpecies=inclH2 or inclH2O or
            inclN2 or inclO2);

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
      replaceable FCSys.Subregions.Species.H2.gas.Fixed H2 if inclH2
        constrainedby FCSys.Subregions.Species.Species(final k=k, final inclLin
          =inclLin) "<html>H<sub>2</sub> model</html>" annotation (
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
      replaceable FCSys.Subregions.Species.H2O.gas.Fixed H2O if inclH2O
        constrainedby FCSys.Subregions.Species.Species(final k=k, final inclLin
          =inclLin) "<html>H<sub>2</sub>O model</html>" annotation (
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
      replaceable FCSys.Subregions.Species.N2.gas.Fixed N2 if inclN2
        constrainedby FCSys.Subregions.Species.Species(final k=k, final inclLin
          =inclLin) "<html>N<sub>2</sub> model</html>" annotation (
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
      replaceable FCSys.Subregions.Species.O2.gas.Fixed O2 if inclO2
        constrainedby FCSys.Subregions.Species.Species(final k=k, final inclLin
          =inclLin) "<html>O<sub>2</sub> model</html>" annotation (
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
      // Exchange
      connect(H2.chemical, chemical.H2) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(H2.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(H2.negativeX.material, negativeX.H2.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.positiveX.material, positiveX.H2.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.negativeY.material, negativeY.H2.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.positiveY.material, positiveY.H2.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.negativeZ.material, negativeZ.H2.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.positiveZ.material, positiveZ.H2.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(H2.negativeX.momentumY, negativeX.H2.momentumY) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.positiveX.momentumY, positiveX.H2.momentumY) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.negativeY.momentumZ, negativeY.H2.momentumZ) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.positiveY.momentumZ, positiveY.H2.momentumZ) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.negativeZ.momentumX, negativeZ.H2.momentumX) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.positiveZ.momentumX, positiveZ.H2.momentumX) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(H2.negativeX.momentumZ, negativeX.H2.momentumZ) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.positiveX.momentumZ, positiveX.H2.momentumZ) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.negativeY.momentumX, negativeY.H2.momentumX) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.positiveY.momentumX, positiveY.H2.momentumX) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.negativeZ.momentumY, negativeZ.H2.momentumY) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.positiveZ.momentumY, positiveZ.H2.momentumY) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(H2.negativeX.entropy, negativeX.H2.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.positiveX.entropy, positiveX.H2.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.negativeY.entropy, negativeY.H2.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.positiveY.entropy, positiveY.H2.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.negativeZ.entropy, negativeZ.H2.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.positiveZ.entropy, positiveZ.H2.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // H2O
      // Exchange
      connect(H2O.chemical, chemical.H2O) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(H2O.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(H2O.negativeX.material, negativeX.H2O.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveX.material, positiveX.H2O.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeY.material, negativeY.H2O.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveY.material, positiveY.H2O.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeZ.material, negativeZ.H2O.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.positiveZ.material, positiveZ.H2O.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(H2O.negativeX.momentumY, negativeX.H2O.momentumY) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveX.momentumY, positiveX.H2O.momentumY) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeY.momentumZ, negativeY.H2O.momentumZ) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveY.momentumZ, positiveY.H2O.momentumZ) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeZ.momentumX, negativeZ.H2O.momentumX) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.positiveZ.momentumX, positiveZ.H2O.momentumX) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(H2O.negativeX.momentumZ, negativeX.H2O.momentumZ) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveX.momentumZ, positiveX.H2O.momentumZ) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeY.momentumX, negativeY.H2O.momentumX) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveY.momentumX, positiveY.H2O.momentumX) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeZ.momentumY, negativeZ.H2O.momentumY) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.positiveZ.momentumY, positiveZ.H2O.momentumY) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(H2O.negativeX.entropy, negativeX.H2O.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveX.entropy, positiveX.H2O.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeY.entropy, negativeY.H2O.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveY.entropy, positiveY.H2O.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeZ.entropy, negativeZ.H2O.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.positiveZ.entropy, positiveZ.H2O.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // N2
      // Exchange
      connect(N2.chemical, chemical.N2) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(N2.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(N2.negativeX.material, negativeX.N2.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.positiveX.material, positiveX.N2.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.negativeY.material, negativeY.N2.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.positiveY.material, positiveY.N2.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.negativeZ.material, negativeZ.N2.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.positiveZ.material, positiveZ.N2.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(N2.negativeX.momentumY, negativeX.N2.momentumY) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.positiveX.momentumY, positiveX.N2.momentumY) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.negativeY.momentumZ, negativeY.N2.momentumZ) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.positiveY.momentumZ, positiveY.N2.momentumZ) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.negativeZ.momentumX, negativeZ.N2.momentumX) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.positiveZ.momentumX, positiveZ.N2.momentumX) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(N2.negativeX.momentumZ, negativeX.N2.momentumZ) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.positiveX.momentumZ, positiveX.N2.momentumZ) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.negativeY.momentumX, negativeY.N2.momentumX) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.positiveY.momentumX, positiveY.N2.momentumX) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.negativeZ.momentumY, negativeZ.N2.momentumY) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.positiveZ.momentumY, positiveZ.N2.momentumY) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(N2.negativeX.entropy, negativeX.N2.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.positiveX.entropy, positiveX.N2.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.negativeY.entropy, negativeY.N2.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.positiveY.entropy, positiveY.N2.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.negativeZ.entropy, negativeZ.N2.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.positiveZ.entropy, positiveZ.N2.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // O2
      // Exchange
      connect(O2.chemical, chemical.O2) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(O2.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(O2.negativeX.material, negativeX.O2.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.positiveX.material, positiveX.O2.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.negativeY.material, negativeY.O2.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.positiveY.material, positiveY.O2.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.negativeZ.material, negativeZ.O2.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.positiveZ.material, positiveZ.O2.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(O2.negativeX.momentumY, negativeX.O2.momentumY) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.positiveX.momentumY, positiveX.O2.momentumY) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.negativeY.momentumZ, negativeY.O2.momentumZ) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.positiveY.momentumZ, positiveY.O2.momentumZ) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.negativeZ.momentumX, negativeZ.O2.momentumX) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.positiveZ.momentumX, positiveZ.O2.momentumX) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(O2.negativeX.momentumZ, negativeX.O2.momentumZ) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.positiveX.momentumZ, positiveX.O2.momentumZ) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.negativeY.momentumX, negativeY.O2.momentumX) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.positiveY.momentumX, positiveY.O2.momentumX) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.negativeZ.momentumY, negativeZ.O2.momentumY) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.positiveZ.momentumY, positiveZ.O2.momentumY) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(O2.negativeX.entropy, negativeX.O2.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.positiveX.entropy, positiveX.O2.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.negativeY.entropy, negativeY.O2.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.positiveY.entropy, positiveY.O2.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.negativeZ.entropy, negativeZ.O2.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.positiveZ.entropy, positiveZ.O2.entropy) annotation (Line(
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
      extends BaseClasses.NullPhase(final hasSpecies=inclC or 'incle-');

      // Conditionally include species.
      parameter Boolean inclC=false "Carbon (C)" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable FCSys.Subregions.Species.C.graphite.Fixed C if inclC
        constrainedby FCSys.Subregions.Species.Species(final k=k, final inclLin
          =inclLin) "C model" annotation (
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
      replaceable FCSys.Subregions.Species.'e-'.graphite.Fixed 'e-' if 'incle-'
        constrainedby FCSys.Subregions.Species.Species(final k=k, final inclLin
          =inclLin) "<html>e<sup>-</sup> model</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable='incle-'),
        Placement(transformation(extent={{-10,-10},{10,10}})));

    equation
      // C
      // Exchange
      connect(C.chemical, chemical.C) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(C.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(C.negativeX.material, negativeX.C.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.positiveX.material, positiveX.C.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.negativeY.material, negativeY.C.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.positiveY.material, positiveY.C.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.negativeZ.material, negativeZ.C.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.positiveZ.material, positiveZ.C.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(C.negativeX.momentumY, negativeX.C.momentumY) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.positiveX.momentumY, positiveX.C.momentumY) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.negativeY.momentumZ, negativeY.C.momentumZ) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.positiveY.momentumZ, positiveY.C.momentumZ) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.negativeZ.momentumX, negativeZ.C.momentumX) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.positiveZ.momentumX, positiveZ.C.momentumX) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(C.negativeX.momentumZ, negativeX.C.momentumZ) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.positiveX.momentumZ, positiveX.C.momentumZ) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.negativeY.momentumX, negativeY.C.momentumX) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.positiveY.momentumX, positiveY.C.momentumX) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.negativeZ.momentumY, negativeZ.C.momentumY) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.positiveZ.momentumY, positiveZ.C.momentumY) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(C.negativeX.entropy, negativeX.C.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.positiveX.entropy, positiveX.C.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.negativeY.entropy, negativeY.C.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.positiveY.entropy, positiveY.C.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C.negativeZ.entropy, negativeZ.C.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C.positiveZ.entropy, positiveZ.C.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // e-
      // Exchange
      connect('e-'.chemical, chemical.'e-') annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect('e-'.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect('e-'.negativeX.material, negativeX.'e-'.material) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.positiveX.material, positiveX.'e-'.material) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.negativeY.material, negativeY.'e-'.material) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.positiveY.material, positiveY.'e-'.material) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.negativeZ.material, negativeZ.'e-'.material) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.positiveZ.material, positiveZ.'e-'.material) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect('e-'.negativeX.momentumY, negativeX.'e-'.momentumY) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.positiveX.momentumY, positiveX.'e-'.momentumY) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.negativeY.momentumZ, negativeY.'e-'.momentumZ) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.positiveY.momentumZ, positiveY.'e-'.momentumZ) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.negativeZ.momentumX, negativeZ.'e-'.momentumX) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.positiveZ.momentumX, positiveZ.'e-'.momentumX) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect('e-'.negativeX.momentumZ, negativeX.'e-'.momentumZ) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.positiveX.momentumZ, positiveX.'e-'.momentumZ) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.negativeY.momentumX, negativeY.'e-'.momentumX) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.positiveY.momentumX, positiveY.'e-'.momentumX) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.negativeZ.momentumY, negativeZ.'e-'.momentumY) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.positiveZ.momentumY, positiveZ.'e-'.momentumY) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect('e-'.negativeX.entropy, negativeX.'e-'.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.positiveX.entropy, positiveX.'e-'.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.negativeY.entropy, negativeY.'e-'.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.positiveY.entropy, positiveY.'e-'.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.negativeZ.entropy, negativeZ.'e-'.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.positiveZ.entropy, positiveZ.'e-'.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      annotation (defaultComponentPrefixes="replaceable", Documentation(info="<html><p>For more information, see the
 <a href=\"modelica://FCSys.Subregions.Phases.BaseClasses.NullPhase\">NullPhase</a> model.</p></html>"));
    end Graphite;

    model Ionomer "Phase to represent ionomer"
      extends BaseClasses.NullPhase(final hasSpecies=inclC19HF37O5S or inclH2O
             or 'inclH+');

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
      replaceable FCSys.Subregions.Species.C19HF37O5S.solid.Fixed C19HF37O5S
        if inclC19HF37O5S constrainedby FCSys.Subregions.Species.Species(final
          k=k, final inclLin=inclLin)
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
      replaceable FCSys.Subregions.Species.H2O.gas.Fixed H2O if inclH2O
        constrainedby FCSys.Subregions.Species.Species(final k=k, final inclLin
          =inclLin) "<html>H<sub>2</sub>O model</html>" annotation (
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
      replaceable FCSys.Subregions.Species.'H+'.solid.Fixed 'H+' if 'inclH+'
        constrainedby FCSys.Subregions.Species.Species(final k=k, final inclLin
          =inclLin) "<html>H<sup>+</sup> model</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable='inclH+'),
        Placement(transformation(extent={{-10,-10},{10,10}})));

    equation
      // C19HF37O5S
      // Exchange
      connect(C19HF37O5S.chemical, chemical.C19HF37O5S) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(C19HF37O5S.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(C19HF37O5S.negativeX.material, negativeX.C19HF37O5S.material)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.positiveX.material, positiveX.C19HF37O5S.material)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.negativeY.material, negativeY.C19HF37O5S.material)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.positiveY.material, positiveY.C19HF37O5S.material)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.negativeZ.material, negativeZ.C19HF37O5S.material)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.positiveZ.material, positiveZ.C19HF37O5S.material)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(C19HF37O5S.negativeX.momentumY, negativeX.C19HF37O5S.momentumY)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.positiveX.momentumY, positiveX.C19HF37O5S.momentumY)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.negativeY.momentumZ, negativeY.C19HF37O5S.momentumZ)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.positiveY.momentumZ, positiveY.C19HF37O5S.momentumZ)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.negativeZ.momentumX, negativeZ.C19HF37O5S.momentumX)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.positiveZ.momentumX, positiveZ.C19HF37O5S.momentumX)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(C19HF37O5S.negativeX.momentumZ, negativeX.C19HF37O5S.momentumZ)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.positiveX.momentumZ, positiveX.C19HF37O5S.momentumZ)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.negativeY.momentumX, negativeY.C19HF37O5S.momentumX)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.positiveY.momentumX, positiveY.C19HF37O5S.momentumX)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.negativeZ.momentumY, negativeZ.C19HF37O5S.momentumY)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.positiveZ.momentumY, positiveZ.C19HF37O5S.momentumY)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(C19HF37O5S.negativeX.entropy, negativeX.C19HF37O5S.entropy)
        annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.positiveX.entropy, positiveX.C19HF37O5S.entropy)
        annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.negativeY.entropy, negativeY.C19HF37O5S.entropy)
        annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.positiveY.entropy, positiveY.C19HF37O5S.entropy)
        annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.negativeZ.entropy, negativeZ.C19HF37O5S.entropy)
        annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.positiveZ.entropy, positiveZ.C19HF37O5S.entropy)
        annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // H+
      // Exchange
      connect('H+'.chemical, chemical.'H+') annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect('H+'.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect('H+'.negativeX.material, negativeX.'H+'.material) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.positiveX.material, positiveX.'H+'.material) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.negativeY.material, negativeY.'H+'.material) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.positiveY.material, positiveY.'H+'.material) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.negativeZ.material, negativeZ.'H+'.material) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.positiveZ.material, positiveZ.'H+'.material) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect('H+'.negativeX.momentumY, negativeX.'H+'.momentumY) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.positiveX.momentumY, positiveX.'H+'.momentumY) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.negativeY.momentumZ, negativeY.'H+'.momentumZ) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.positiveY.momentumZ, positiveY.'H+'.momentumZ) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.negativeZ.momentumX, negativeZ.'H+'.momentumX) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.positiveZ.momentumX, positiveZ.'H+'.momentumX) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect('H+'.negativeX.momentumZ, negativeX.'H+'.momentumZ) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.positiveX.momentumZ, positiveX.'H+'.momentumZ) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.negativeY.momentumX, negativeY.'H+'.momentumX) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.positiveY.momentumX, positiveY.'H+'.momentumX) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.negativeZ.momentumY, negativeZ.'H+'.momentumY) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.positiveZ.momentumY, positiveZ.'H+'.momentumY) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect('H+'.negativeX.entropy, negativeX.'H+'.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.positiveX.entropy, positiveX.'H+'.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.negativeY.entropy, negativeY.'H+'.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.positiveY.entropy, positiveY.'H+'.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.negativeZ.entropy, negativeZ.'H+'.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.positiveZ.entropy, positiveZ.'H+'.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // H2O
      // Exchange

      connect(H2O.chemical, chemical.H2O) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));

      connect(H2O.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(H2O.negativeX.material, negativeX.H2O.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveX.material, positiveX.H2O.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeY.material, negativeY.H2O.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveY.material, positiveY.H2O.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeZ.material, negativeZ.H2O.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.positiveZ.material, positiveZ.H2O.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(H2O.negativeX.momentumY, negativeX.H2O.momentumY) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveX.momentumY, positiveX.H2O.momentumY) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeY.momentumZ, negativeY.H2O.momentumZ) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveY.momentumZ, positiveY.H2O.momentumZ) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeZ.momentumX, negativeZ.H2O.momentumX) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.positiveZ.momentumX, positiveZ.H2O.momentumX) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(H2O.negativeX.momentumZ, negativeX.H2O.momentumZ) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveX.momentumZ, positiveX.H2O.momentumZ) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeY.momentumX, negativeY.H2O.momentumX) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveY.momentumX, positiveY.H2O.momentumX) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeZ.momentumY, negativeZ.H2O.momentumY) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.positiveZ.momentumY, positiveZ.H2O.momentumY) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(H2O.negativeX.entropy, negativeX.H2O.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveX.entropy, positiveX.H2O.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeY.entropy, negativeY.H2O.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveY.entropy, positiveY.H2O.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeZ.entropy, negativeZ.H2O.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.positiveZ.entropy, positiveZ.H2O.entropy) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      annotation (
        defaultComponentPrefixes="replaceable",
        Documentation(info="<html><p>For more information, see the
 <a href=\"modelica://FCSys.Subregions.Phases.BaseClasses.NullPhase\">NullPhase</a> model.</p></html>"),

        Diagram(graphics));
    end Ionomer;

    model Liquid "Phase to represent liquid"

      extends BaseClasses.NullPhase(final hasSpecies=inclH2O);

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
      replaceable FCSys.Subregions.Species.H2O.gas.Fixed H2O if inclH2O
        constrainedby FCSys.Subregions.Species.Species(final k=k, final inclLin
          =inclLin) "<html>H<sub>2</sub>O model</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclH2O),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      // TODO: Add and use a model for H2O liquid.
    equation
      // H2O
      // Exchange
      connect(H2O.chemical, chemical.H2O) annotation (Line(
          points={{-7,7},{-20,20}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(H2O.inert, phaseBoundary.inertD) annotation (Line(
          points={{7,-7},{7,-7},{7,-7}},
          color={72,90,180},
          smooth=Smooth.None));
      // Material transport
      connect(H2O.negativeX.material, negativeX.H2O.material) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveX.material, positiveX.H2O.material) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeY.material, negativeY.H2O.material) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveY.material, positiveY.H2O.material) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeZ.material, negativeZ.H2O.material) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.positiveZ.material, positiveZ.H2O.material) annotation (Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 1st transverse transport of linear momentum
      connect(H2O.negativeX.momentumY, negativeX.H2O.momentumY) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveX.momentumY, positiveX.H2O.momentumY) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeY.momentumZ, negativeY.H2O.momentumZ) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveY.momentumZ, positiveY.H2O.momentumZ) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeZ.momentumX, negativeZ.H2O.momentumX) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.positiveZ.momentumX, positiveZ.H2O.momentumX) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // 2nd transverse transport of linear momentum
      connect(H2O.negativeX.momentumZ, negativeX.H2O.momentumZ) annotation (
          Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveX.momentumZ, positiveX.H2O.momentumZ) annotation (
          Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeY.momentumX, negativeY.H2O.momentumX) annotation (
          Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveY.momentumX, positiveY.H2O.momentumX) annotation (
          Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeZ.momentumY, negativeZ.H2O.momentumY) annotation (
          Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.positiveZ.momentumY, positiveZ.H2O.momentumY) annotation (
          Line(
          points={{-7,-7},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // Thermal transport
      connect(H2O.negativeX.entropy, negativeX.H2O.entropy) annotation (Line(
          points={{-10,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},{-40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveX.entropy, positiveX.H2O.entropy) annotation (Line(
          points={{10,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeY.entropy, negativeY.H2O.entropy) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.positiveY.entropy, positiveY.H2O.entropy) annotation (Line(
          points={{6.10623e-16,10},{-4.87687e-22,20},{5.55112e-16,20},{
              5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.negativeZ.entropy, negativeZ.H2O.entropy) annotation (Line(
          points={{7,7},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.positiveZ.entropy, positiveZ.H2O.entropy) annotation (Line(
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
        outer parameter Q.Length L[3](each final min=Modelica.Constants.small)
          "Length" annotation (HideResult=true,missingInnerMessage=
              "This model should be used within the Subregion model.");
        outer parameter Q.Area A[3] "Cross-sectional area" annotation (
            HideResult=true, missingInnerMessage=
              "This model should be used within the Subregion model.");

        // Assumptions
        parameter Boolean inclLin[3]={true,false,false}
          "true, if each component of linear momentum is included"
          annotation (Evaluate=true,Dialog(tab="Assumptions"));
        final parameter Integer n_lin=countTrue(inclLin)
          "Number of components of linear momentum"
          annotation (Evaluate=true, HideResult=true);

        FCSys.Connectors.InertAmagat inert(final n_lin=n_lin) if hasSpecies
          annotation (Placement(transformation(extent={{10,-30},{30,-10}}),
              iconTransformation(extent={{70,-90},{90,-70}})));
        FCSys.Connectors.FaceBus negativeX "Negative face along the x axis"
          annotation (Placement(transformation(extent={{-50,-10},{-30,10}}),
              iconTransformation(extent={{-90,-10},{-70,10}})));
        FCSys.Connectors.FaceBus positiveX "Positive face along the x axis"
          annotation (Placement(transformation(extent={{30,-10},{50,10}}),
              iconTransformation(extent={{70,-10},{90,10}})));
        FCSys.Connectors.FaceBus negativeY "Negative face along the y axis"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}}),
              iconTransformation(extent={{-10,-94},{10,-74}})));
        FCSys.Connectors.FaceBus positiveY "Positive face along the y axis"
          annotation (Placement(transformation(extent={{-10,30},{10,50}}),
              iconTransformation(extent={{-10,90},{10,110}})));
        FCSys.Connectors.FaceBus negativeZ "Negative face along the z axis"
          annotation (Placement(transformation(extent={{10,10},{30,30}}),
              iconTransformation(extent={{40,40},{60,60}})));
        FCSys.Connectors.FaceBus positiveZ "Positive face along the z axis"
          annotation (Placement(transformation(extent={{-30,-30},{-10,-10}}),
              iconTransformation(extent={{-90,-90},{-70,-70}})));

        PhaseBoundary phaseBoundary(final n_lin=n_lin) if hasSpecies
          "Phase boundary" annotation (Placement(transformation(
              extent={{-18,-18},{18,18}},
              rotation=0,
              origin={0,0})));
        // This component is conditional because if two or more empty phases
        // (without any species included) are connected within a subregion, there
        // would be a mathematical singularity.

        FCSys.Connectors.ChemicalBus chemical annotation (Placement(
              transformation(extent={{-30,10},{-10,30}}), iconTransformation(
                extent={{-60,40},{-40,60}})));

      protected
        parameter Boolean hasSpecies "true, if any species are included";

      equation
        // Inert interactions
        connect(phaseBoundary.inertA, inert) annotation (Line(
            points={{12,-12},{20,-20}},
            color={72,90,180},
            smooth=Smooth.None));

        annotation (
          Documentation(info="<html><p>Assumptions:
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
      package graphite "C graphite"
        extends Modelica.Icons.Package;
        model Calibrated "Ideal correlations, with adjustment factors"
          extends SpeciesInertStagnant(
            redeclare FCSys.Characteristics.C.graphite Data,
            redeclare Q.Resistivity alpha_Phi,
            redeclare Q.Resistivity alpha_S,
            redeclare Q.Resistivity beta_S,
            redeclare final parameter Q.Permittivity epsilon=0,
            final Deltamu_IC=0);
          // Note:  Permittivity and double layer electrochemical potential don't
          // apply since the species is electrically neutral.

          parameter Q.NumberAbsolute k_alpha_Phi(final nominal=1) = 1
            "<html>Adjustment factor for exchange of linear momentum (<i>k</i><sub>&alpha;&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_alpha_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal exchange (<i>k</i><sub>&alpha;<i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal transport (<i>k</i>&beta;<sub><i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
        equation
          alpha_Phi = k_alpha_Phi*Data.beta(T);
          alpha_S = k_alpha_S*Data.beta(T);
          beta_S = k_beta_S*Data.beta(T);
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
            redeclare FCSys.Characteristics.C.graphite Data,
            redeclare Q.Resistivity alpha_Phi,
            redeclare Q.Resistivity alpha_S,
            redeclare Q.Resistivity beta_S,
            redeclare final parameter Q.Permittivity epsilon=0,
            final Deltamu_IC=0);
          // Note:  Permittivity and double layer electrochemical potential don't
          // apply since the species is electrically neutral.

        equation
          alpha_Phi = Data.beta(T);
          alpha_S = Data.beta(T);
          beta_S = Data.beta(T);
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
            redeclare FCSys.Characteristics.C.graphite Data(
              Deltah0_f=0,
              Deltah0=0,
              specHeatCapPow=0,
              T_lim_c0={0,Modelica.Constants.inf},
              b_c0=[935*U.J*Data.m/(U.kg*U.K)],
              B_c0=[-298.15*U.K*935*U.J*Data.m/(U.kg*U.K) + Data.Deltah0_f, 0]),

            redeclare parameter Q.Resistivity alpha_Phi=Data.beta(),
            redeclare parameter Q.Resistivity alpha_S=Data.beta(),
            redeclare parameter Q.Resistivity beta_S=U.m*U.K/(11.1*U.W),
            redeclare final parameter Q.Permittivity epsilon=0,
            final Deltamu_IC=0);

          // See the documentation for a table of values.
          // Note:  Permittivity and double layer electrochemical potential don't
          // apply since the species is electrically neutral.
          // Note:  Parameter expressions (e.g., involving defaults.T) are not used
          // here since they would render the parameters unadjustable in Dymola 7.4.
          // A similar note applies to the other species.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Fixed specific heat capacity (independent of temperature)
    <li>Transport and exchange properties (&alpha;<sub><i>S</i></sub>, &beta;<sub><i>S</i></sub>, etc.) are fixed (e.g., independent of temperature)</li>
    </ol></p>

    <p>The default specific heat capacity at constant pressure (<code>b_c0=[0, 935*U.J*Data.m/(U.kg*U.K)]</code>) and thermal
   resistivity (<code>beta_S=U.m*U.K/(11.1*U.W)</code>) is based on data of graphite fiber epoxy (25% vol)<br>composite at 300 K from
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
      <th rowspan=2 valign=\"middle\"><code>beta_S<br>*U.K<br>/(U.m<br>*U.K)</code></th>
      <th rowspan=2 valign=\"middle\"><code>beta_S<br>*U.W<br>/(U.m<br>*U.K)</code></th>
      <th rowspan=2 valign=\"middle\"><code>c0*U.kg<br>*U.W<br>/(U.J<br>*Data.m)</code></th>
      <th rowspan=1 colspan=2><code>beta_S*U.W/(U.m*U.K)</code></th>
      <th rowspan=2 valign=\"middle\"><code>c0*U.kg<br>*U.K<br>/(U.J<br>*Data.m)</code></th>
      <th rowspan=1 colspan=2><code>beta_S*U.W/(U.m*U.K)</code></th>
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
      end graphite;
    end C;

    package C19HF37O5S
      "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S</html>"
      extends Modelica.Icons.Package;
      package solid
        "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S solid</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Ideal correlations, with adjustment factors"
          extends SpeciesInertStagnant(
            redeclare FCSys.Characteristics.C19HF37O5S.solid Data,
            redeclare Q.Resistivity alpha_Phi,
            redeclare Q.Resistivity alpha_S,
            redeclare Q.Resistivity beta_S,
            redeclare final parameter Q.Permittivity epsilon=0,
            final Deltamu_IC=0);
          // Note:  Permittivity and double layer electrochemical potential don't
          // apply since the species is electrically neutral.

          parameter Q.NumberAbsolute k_alpha_Phi(final nominal=1) = 1
            "<html>Adjustment factor for exchange of linear momentum (<i>k</i><sub>&alpha;&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_alpha_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal exchange (<i>k</i><sub>&alpha;<i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal transport (<i>k</i>&beta;<sub><i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
        equation
          alpha_Phi = k_alpha_Phi*Data.beta(T);
          alpha_S = k_alpha_S*Data.beta(T);
          beta_S = k_beta_S*Data.beta(T);
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
            redeclare FCSys.Characteristics.C19HF37O5S.solid Data,
            redeclare Q.Resistivity alpha_Phi,
            redeclare Q.Resistivity alpha_S,
            redeclare Q.Resistivity beta_S,
            redeclare final parameter Q.Permittivity epsilon=0,
            final Deltamu_IC=0);
          // Note:  Permittivity and double layer electrochemical potential don't
          // apply since the species is electrically neutral.

        equation
          alpha_Phi = Data.beta(T);
          alpha_S = Data.beta(T);
          beta_S = Data.beta(T);
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C19HF37O5S",
            Documentation(info=
                  "<html><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends SpeciesInertStagnant(
            redeclare FCSys.Characteristics.C19HF37O5S.solid Data,
            redeclare parameter Q.Resistivity alpha_Phi=Data.beta(),
            redeclare parameter Q.Resistivity alpha_S=Data.beta(),
            redeclare parameter Q.Resistivity beta_S=U.m*U.K/(0.16*U.W),
            redeclare final parameter Q.Permittivity epsilon=0,
            final Deltamu_IC=0);
          // Note:  Permittivity and double layer electrochemical potential don't
          // apply since the species is electrically neutral.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C19HF37O5S",
            Documentation(info="<html><p>Assumptions:
    <ol>
    <li>Transport and exchange properties (&alpha;<sub><i>S</i></sub>, &beta;<sub><i>S</i></sub>, etc.) are fixed (e.g., independent of temperature)</li>
    </ol></p>

    <p>Notes:
    <ul><li>The default thermal transport resistivity (<code>beta_S=U.m*U.K/(0.16*U.W)</code>) is of dry
  Nafion 115 [<a href=\"modelica://FCSys.UsersGuide.References\">Kandlikar2009</a>, p. 1277.</li>
  </ul>
  </p><p>For more information, see the
    <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model.</p></html>"));
        end Fixed;
      end solid;
    end C19HF37O5S;

    package 'e-' "<html>e<sup>-</sup></html>"
      extends Modelica.Icons.Package;
      package graphite "<html>e<sup>-</sup> in graphite</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Ideal correlations, with adjustment factors"
          extends Species(
            redeclare FCSys.Characteristics.'e-'.graphite Data,
            redeclare Q.Resistivity alpha_N,
            redeclare Q.Resistivity alpha_Phi,
            redeclare Q.Resistivity alpha_S,
            redeclare Q.Resistivity beta_N,
            redeclare Q.Resistivity beta_Phi,
            redeclare Q.Resistivity beta_S,
            partNumInitMeth=InitMethScalar.PotentialElectrochemical,
            Deltamu(fixed=true));

          parameter Q.Capacitance C=1*U.micro*U.F "Electrical capacitance"
            annotation (group="Material properties");

          parameter Q.NumberAbsolute k_alpha_N(final nominal=1) = 1
            "<html>Adjustment factor for reaction (<i>k</i><sub>&alpha;<i>N</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_alpha_Phi(final nominal=1) = 1
            "<html>Adjustment factor for exchange of linear momentum (<i>k</i><sub>&alpha;&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_alpha_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal exchange (<i>k</i><sub>&alpha;<i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_N(final nominal=1) = 1
            "<html>Adjustment factor for material transport (<i>k</i>&beta;<sub><i>N</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_Phi(final nominal=1) = 1
            "<html>Adjustment factor for transport of linear momentum (<i>k</i>&beta;<sub>&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal transport (<i>k</i>&beta;<sub><i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
        equation
          alpha_N = k_alpha_N*Data.beta(T);
          alpha_Phi = k_alpha_Phi*Data.beta(T);
          alpha_S = k_alpha_S*Data.beta(T);
          beta_N = k_beta_N*Data.beta(T);
          beta_Phi = k_beta_Phi*Data.beta(T);
          beta_S = k_beta_S*Data.beta(T);
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'e-'",
            Documentation(info=
                  "<html><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model.</p></html>"));
        end Calibrated;

        model Correlated "Ideal correlations from kinetic theory"
          extends Species(
            redeclare FCSys.Characteristics.'e-'.graphite Data,
            redeclare Q.Resistivity alpha_N,
            redeclare Q.Resistivity alpha_Phi,
            redeclare Q.Resistivity alpha_S,
            redeclare Q.Resistivity beta_N,
            redeclare Q.Resistivity beta_Phi,
            redeclare Q.Resistivity beta_S,
            partNumInitMeth=InitMethScalar.PotentialElectrochemical,
            Deltamu(fixed=true));

        equation
          alpha_N = Data.beta(T);
          alpha_Phi = Data.beta(T);
          alpha_S = Data.beta(T);
          beta_N = Data.beta(T);

          beta_Phi = Data.beta(T);
          beta_S = Data.beta(T);
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'e-'",
            Documentation(info=
                  "<html><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare FCSys.Characteristics.'e-'.graphite Data,
            redeclare parameter Q.Resistivity alpha_N=1e23*Data.beta(),
            redeclare parameter Q.Resistivity alpha_Phi=Data.beta(),
            redeclare parameter Q.Resistivity alpha_S=Data.beta(),
            redeclare parameter Q.Resistivity beta_N=Data.beta(),
            redeclare parameter Q.Resistivity beta_Phi=Data.beta(),
            redeclare parameter Q.Resistivity beta_S=Data.beta(),
            redeclare parameter Q.Permittivity epsilon=1e3*U.epsilon_0,
            partNumInitMeth=InitMethScalar.PotentialElectrochemical,
            Deltamu(fixed=true),
            Lstar=1e6*product(L)^(1/3));
          // **temp Lstar
          annotation (
            group="Material properties",
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'e-'",
            Documentation(info="<html>
          <p>Notes:<ul>
          <li>The default value for <code>alpha_N</code> (<code>1e4*Data.beta()</code>)
          gives reasonable reaction rates for the HOR and ORR when used in conjunction
          with the defaults for H<sub>2</sub>, H<sub>2</sub>O, H<sup>+</sup>, and O<sub>2</sub>.
          The default <code>alpha_N</code> is greater (higher resistance) for H<sub>2</sub>O
          and O<sub>2</sub> to slow the ORR.</li>
          </ul></p>
    <p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model.</p></html>"),

            Diagram(graphics));

        end Fixed;
      end graphite;
    end 'e-';

    package H2 "<html>H<sub>2</sub></html>"
      extends Modelica.Icons.Package;
      package gas "<html>H<sub>2</sub> gas</html>"
        extends Modelica.Icons.Package;
        model Calibrated
          "Ideal or empirical correlations, w/ adjustment factors"
          extends SpeciesNeutral(
            redeclare FCSys.Characteristics.H2.gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare Q.Resistivity alpha_N,
            redeclare Q.Resistivity alpha_Phi,
            redeclare Q.Resistivity alpha_S,
            redeclare Q.Resistivity beta_N,
            redeclare Q.Resistivity beta_Phi,
            redeclare Q.Resistivity beta_S);

          parameter Q.NumberAbsolute k_alpha_N(final nominal=1) = 1
            "<html>Adjustment factor for reaction (<i>k</i><sub>&alpha;<i>N</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_alpha_Phi(final nominal=1) = 1
            "<html>Adjustment factor for exchange of linear momentum (<i>k</i><sub>&alpha;&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_alpha_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal exchange (<i>k</i><sub>&alpha;<i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_N(final nominal=1) = 1
            "<html>Adjustment factor for material transport (<i>k</i>&beta;<sub><i>N</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_Phi(final nominal=1) = 1
            "<html>Adjustment factor for transport of linear momentum (<i>k</i>&beta;<sub>&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal transport (<i>k</i>&beta;<sub><i>S</i></sub>)</html>"
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
        equation
          alpha_N = k_alpha_N*Data.beta(T);
          alpha_Phi = k_alpha_Phi*Data.beta(T);
          alpha_S = k_alpha_S*Data.beta(T);
          beta_N = k_beta_N*Data.beta(T);

          beta_Phi = k_beta_Phi*(if empiricalTransverse then Data.beta_Phi(T)
             else Data.beta(T));
          beta_S = k_beta_S*(if empiricalThermal then Data.beta_S(T) else
            Data.beta(T));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"));
        end Calibrated;

        model Correlated "Ideal (kinetic) or empirical correlations"
          extends SpeciesNeutral(
            redeclare FCSys.Characteristics.H2.gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare Q.Resistivity alpha_N,
            redeclare Q.Resistivity alpha_Phi,
            redeclare Q.Resistivity alpha_S,
            redeclare Q.Resistivity beta_N,
            redeclare Q.Resistivity beta_Phi,
            redeclare Q.Resistivity beta_S);

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
        equation
          alpha_N = Data.beta(T);
          alpha_Phi = Data.beta(T);
          alpha_S = Data.beta(T);
          beta_N = Data.beta(T);

          beta_Phi = (if empiricalTransverse then Data.beta_Phi(T) else
            Data.beta(T));
          beta_S = (if empiricalThermal then Data.beta_S(T) else Data.beta(T));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends SpeciesNeutral(
            redeclare FCSys.Characteristics.H2.gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare parameter Q.Resistivity alpha_N=1e4*Data.beta(),
            redeclare parameter Q.Resistivity alpha_Phi=Data.beta(),
            redeclare parameter Q.Resistivity alpha_S=U.m*U.K/(183e-3*U.W),
            redeclare parameter Q.Resistivity beta_N=Data.beta(),
            redeclare parameter Q.Resistivity beta_Phi=Data.m/(89.6e-7*U.Pa*U.s),

            redeclare parameter Q.Resistivity beta_S=U.m*U.K/(183e-3*U.W)/100);

          // See the documentation for a table of values.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
    <li>Transport and exchange properties (&alpha;<sub>&Phi;</sub>, &beta;<sub>&Phi;</sub>, etc.) are fixed (e.g., independent of temperature)</li>
</ol></p>

<p>Additional notes:<ul>
          <li>The default value for <code>alpha_N</code> (<code>1e4*Data.beta()</code>)
          gives reasonable reaction rates for the HOR and ORR when used in conjunction
          with the defaults for e<sup>-</sup>, H<sub>2</sub>O, H<sup>+</sup>, and O<sub>2</sub>.
          The default <code>alpha_N</code> is greater (higher resistance) for H<sub>2</sub>O
          and O<sub>2</sub> to slow the ORR.</li>
<li>
  The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

<p>The default transport resistivities (<code>beta_Phi=Data.m/(89.6e-7*U.Pa*U.s)</code>
and <code>beta_S=U.m*U.K/(183e-3*U.W)</code>) are based on data of H<sub>2</sub> gas at 1 atm and
  300 K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 919&ndash;920].
  Table 1 lists the properties at  other temperatures. </p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"bottom\"><b>Table 1:</b> Properties of H<sub>2</sub> gas at 1 atm (<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 919&ndash;920).</caption>
  <tr>
      <th rowspan=2 valign=\"middle\"><code>T/U.K</code></th>
      <th rowspan=2 width=1><code>c*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th rowspan=2 width=1 ><code>beta_Phi<br>*U.Pa*U.s/Data.m</code></th>
      <th rowspan=2 width=1 ><code>beta_S*U.W<br>/(U.m*U.K)</code></th>
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
      end gas;
    end H2;

    package H2O "<html>H<sub>2</sub>O</html>"
      extends Modelica.Icons.Package;
      package gas "<html>H<sub>2</sub>O gas</html>"
        extends Modelica.Icons.Package;
        model Calibrated
          "Ideal or empirical correlations, w/ adjustment factors"
          extends SpeciesNeutral(
            redeclare FCSys.Characteristics.H2O.gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare Q.Resistivity alpha_N,
            redeclare Q.Resistivity alpha_Phi,
            redeclare Q.Resistivity alpha_S,
            redeclare Q.Resistivity beta_N,
            redeclare Q.Resistivity beta_Phi,
            redeclare Q.Resistivity beta_S);

          parameter Q.NumberAbsolute k_alpha_N(final nominal=1) = 1
            "<html>Adjustment factor for reaction (<i>k</i><sub>&alpha;<i>N</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_alpha_Phi(final nominal=1) = 1
            "<html>Adjustment factor for exchange of linear momentum (<i>k</i><sub>&alpha;&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_alpha_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal exchange (<i>k</i><sub>&alpha;<i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_N(final nominal=1) = 1
            "<html>Adjustment factor for material transport (<i>k</i>&beta;<sub><i>N</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_Phi(final nominal=1) = 1
            "<html>Adjustment factor for transport of linear momentum (<i>k</i>&beta;<sub>&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal transport (<i>k</i>&beta;<sub><i>S</i></sub>)</html>"
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
        equation
          alpha_N = k_alpha_N*Data.beta(T);
          alpha_Phi = k_alpha_Phi*Data.beta(T);
          alpha_S = k_alpha_S*Data.beta(T);
          beta_N = k_beta_N*Data.beta(T);

          beta_Phi = k_beta_Phi*(if empiricalTransverse then Data.beta_Phi(T)
             else Data.beta(T));
          beta_S = k_beta_S*(if empiricalThermal then Data.beta_S(T) else
            Data.beta(T));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"));
        end Calibrated;

        model Correlated "Ideal (kinetic) or empirical correlations"
          extends SpeciesNeutral(
            redeclare FCSys.Characteristics.H2O.gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare Q.Resistivity alpha_N,
            redeclare Q.Resistivity alpha_Phi,
            redeclare Q.Resistivity alpha_S,
            redeclare Q.Resistivity beta_N,
            redeclare Q.Resistivity beta_Phi,
            redeclare Q.Resistivity beta_S);

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
        equation
          alpha_N = Data.beta(T);
          alpha_Phi = Data.beta(T);
          alpha_S = Data.beta(T);
          beta_N = Data.beta(T);

          beta_Phi = (if empiricalTransverse then Data.beta_Phi(T) else
            Data.beta(T));
          beta_S = (if empiricalThermal then Data.beta_S(T) else Data.beta(T));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends SpeciesNeutral(
            redeclare FCSys.Characteristics.H2O.gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare parameter Q.Resistivity alpha_N=1e10*Data.beta(),
            redeclare parameter Q.Resistivity alpha_Phi=Data.beta(),
            redeclare parameter Q.Resistivity alpha_S=U.m*U.K/(19.6e-3*U.W),
            redeclare parameter Q.Resistivity beta_N=Data.beta(),
            redeclare parameter Q.Resistivity beta_Phi=Data.m/(9.09e-6*U.Pa*U.s),

            redeclare parameter Q.Resistivity beta_S=U.m*U.K/(19.6e-3*U.W));

          // See the documentation for tables of values.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
    <li>Transport and exchange properties (&alpha;<sub>&Phi;</sub>, &beta;<sub>&Phi;</sub>, etc.) are fixed (e.g., independent of temperature)</li>
  </ol></p>

          <p>Notes:<ul>
          <li>The default value for <code>alpha_N</code> (<code>1e10*Data.beta()</code>)
          gives reasonable reaction rates for the HOR and ORR when used in conjunction
          with the defaults for e<sup>-</sup>, H<sub>2</sub>, H<sub>2</sub>, H<sup>+</sup>, and O<sub>2</sub>.
          The default <code>alpha_N</code> is greater (higher resistance) for H<sub>2</sub>O
          and O<sub>2</sub> to slow the ORR.</li>
<ul><li>
  The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

<p>The default transport resistivities (<code>beta_Phi=Data.m/(9.09e-6*U.Pa*U.s)</code>
and <code>beta_S=U.m*U.K/(19.6e-3*U.W)</code>) are of H<sub>2</sub>O gas at saturation pressure and
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
      <th rowspan=1 colspan=2 width=1><code>beta_Phi<br>*U.Pa*U.s/Data.m</code></th>
      <th rowspan=1 colspan=2 width=1><code>beta_S*U.W<br>/(U.m*U.K)</code></th>
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
      <th width=1 ><code>beta_Phi<br>*U.Pa*U.s/Data.m</code></th>
      <th width=1 ><code>beta_S*U.W<br>/(U.m*U.K)</code></th>
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
      end gas;
    end H2O;

    package 'H+' "<html>H<sup>+</sup></html>"
      extends Modelica.Icons.Package;
      package solid "<html>H<sup>+</sup> in solid</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Ideal correlations, with adjustment factors"
          extends Species(
            redeclare FCSys.Characteristics.'H+'.solid Data,
            redeclare Q.Resistivity alpha_N,
            redeclare Q.Resistivity alpha_Phi,
            redeclare Q.Resistivity alpha_S,
            redeclare Q.Resistivity beta_N,
            redeclare Q.Resistivity beta_Phi,
            redeclare Q.Resistivity beta_S,
            partNumInitMeth=InitMethScalar.PotentialElectrochemical,
            Deltamu(fixed=true));

          parameter Q.NumberAbsolute k_alpha_N(final nominal=1) = 1
            "<html>Adjustment factor for reaction (<i>k</i><sub>&alpha;<i>N</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_alpha_Phi(final nominal=1) = 1
            "<html>Adjustment factor for exchange of linear momentum (<i>k</i><sub>&alpha;&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_alpha_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal exchange (<i>k</i><sub>&alpha;<i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_N(final nominal=1) = 1
            "<html>Adjustment factor for material transport (<i>k</i>&beta;<sub><i>N</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_Phi(final nominal=1) = 1
            "<html>Adjustment factor for transport of linear momentum (<i>k</i>&beta;<sub>&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal transport (<i>k</i>&beta;<sub><i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
        equation
          alpha_N = k_alpha_N*Data.beta(T);
          alpha_Phi = k_alpha_Phi*Data.beta(T);
          alpha_S = k_alpha_S*Data.beta(T);
          beta_N = k_beta_N*Data.beta(T);

          beta_Phi = k_beta_Phi*Data.beta(T);
          beta_S = k_beta_S*Data.beta(T);
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'H+'",
            Documentation(info=
                  "<html><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model.</p></html>"));
        end Calibrated;

        model Correlated "Ideal correlations from kinetic theory"
          extends Species(
            redeclare FCSys.Characteristics.'H+'.solid Data,
            redeclare Q.Resistivity alpha_N,
            redeclare Q.Resistivity alpha_Phi,
            redeclare Q.Resistivity alpha_S,
            redeclare Q.Resistivity beta_N,
            redeclare Q.Resistivity beta_Phi,
            redeclare Q.Resistivity beta_S,
            partNumInitMeth=InitMethScalar.PotentialElectrochemical,
            Deltamu(fixed=true));

        equation
          alpha_N = Data.beta(T);
          alpha_Phi = Data.beta(T);
          alpha_S = Data.beta(T);
          beta_N = Data.beta(T);

          beta_Phi = Data.beta(T);
          beta_S = Data.beta(T);
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'H+'",
            Documentation(info=
                  "<html><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare FCSys.Characteristics.'H+'.solid Data,
            redeclare parameter Q.Resistivity alpha_N=1e5*Data.beta(),
            redeclare parameter Q.Resistivity alpha_Phi=Data.beta(),
            redeclare parameter Q.Resistivity alpha_S=U.m*U.K/(0.1661*U.W),
            redeclare parameter Q.Resistivity beta_N=Data.beta(),
            redeclare parameter Q.Resistivity beta_Phi=Data.m/(5.3e-6*U.Pa*U.s),

            redeclare parameter Q.Resistivity beta_S=U.m*U.K/(0.1661*U.W),
            redeclare parameter Q.Permittivity epsilon=1e10*U.epsilon_0,
            partNumInitMeth=InitMethScalar.PotentialElectrochemical,
            Deltamu(fixed=true),
            Lstar=1e6*product(L)^(1/3));

          // **temp Lstar

          // See the documentation for a table of values.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'H+'",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Transport and exchange properties (&alpha;<sub>&Phi;</sub>, &beta;<sub>&Phi;</sub>, etc.) are fixed (e.g., independent of temperature)</li>
    </ol></p>

              <p>Additional notes:<ul>
          <li>The default value for <code>alpha_N</code> (<code>1e4*Data.beta()</code>)
          gives reasonable reaction rates for the HOR and ORR when used in conjunction
          with the defaults for e<sup>-</sup>, H<sub>2</sub>, H<sub>2</sub>O, and O<sub>2</sub>.
          The default <code>alpha_N</code> is greater (higher resistance) for H<sub>2</sub>O
          and O<sub>2</sub> to slow the ORR.</li>
          </ul></p>

  <p>The default transport resistivities (<code>beta_Phi=Data.m/(5.3e-6*U.Pa*U.s)</code> and <code>beta_S=U.m*U.K/(0.1661*U.W)</code>) are of H gas
  (rather than H<sup>+</sup>) at 300 K from [<a href=\"modelica://FCSys.UsersGuide.References\">Schetz1996</a>, p. 139].
  Table 1 lists the properties at other temperatures.</p>

    <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"bottom\"><b>Table 1:</b> Properties of H gas (not H<sup>+</sup>) (<a href=\"modelica://FCSys.UsersGuide.References\">Schetz and Fuhs, 1996</a>, p. 139)</caption>
<tr>
      <th valign=\"middle\"><code>T/U.K</code></th>
      <th width=1><code>alpha_Phi<br>*U.Pa*U.s/Data.m</code></th>
      <th width=1><code>alpha_S*U.W<br>/(U.m*U.K)</code></th>
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
      end solid;
    end 'H+';

    package N2 "<html>N<sub>2</sub></html>"
      extends Modelica.Icons.Package;
      package gas "<html>N<sub>2</sub> gas</html>"
        extends Modelica.Icons.Package;
        model Calibrated
          "Ideal or empirical correlations, w/ adjustment factors"
          extends SpeciesNeutral(
            redeclare FCSys.Characteristics.N2.gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare Q.Resistivity alpha_N,
            redeclare Q.Resistivity alpha_Phi,
            redeclare Q.Resistivity alpha_S,
            redeclare Q.Resistivity beta_N,
            redeclare Q.Resistivity beta_Phi,
            redeclare Q.Resistivity beta_S);

          parameter Q.NumberAbsolute k_alpha_N(final nominal=1) = 1
            "<html>Adjustment factor for reaction (<i>k</i><sub>&alpha;<i>N</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_alpha_Phi(final nominal=1) = 1
            "<html>Adjustment factor for exchange of linear momentum (<i>k</i><sub>&alpha;&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_alpha_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal exchange (<i>k</i><sub>&alpha;<i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_N(final nominal=1) = 1
            "<html>Adjustment factor for material transport (<i>k</i>&beta;<sub><i>N</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_Phi(final nominal=1) = 1
            "<html>Adjustment factor for transport of linear momentum (<i>k</i>&beta;<sub>&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal transport (<i>k</i>&beta;<sub><i>S</i></sub>)</html>"
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
        equation
          alpha_N = k_alpha_N*Data.beta(T);
          alpha_Phi = k_alpha_Phi*Data.beta(T);
          alpha_S = k_alpha_S*Data.beta(T);
          beta_N = k_beta_N*Data.beta(T);

          beta_Phi = k_beta_Phi*(if empiricalTransverse then Data.beta_Phi(T)
             else Data.beta(T));
          beta_S = k_beta_S*(if empiricalThermal then Data.beta_S(T) else
            Data.beta(T));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="N2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"));
        end Calibrated;

        model Correlated "Ideal (kinetic) or empirical correlations"
          extends SpeciesNeutral(
            redeclare FCSys.Characteristics.N2.gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare Q.Resistivity alpha_N,
            redeclare Q.Resistivity alpha_Phi,
            redeclare Q.Resistivity alpha_S,
            redeclare Q.Resistivity beta_N,
            redeclare Q.Resistivity beta_Phi,
            redeclare Q.Resistivity beta_S);

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
        equation
          alpha_N = Data.beta(T);
          alpha_Phi = Data.beta(T);
          alpha_S = Data.beta(T);
          beta_N = Data.beta(T);

          beta_Phi = (if empiricalTransverse then Data.beta_Phi(T) else
            Data.beta(T));
          beta_S = (if empiricalThermal then Data.beta_S(T) else Data.beta(T));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="N2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends SpeciesNeutral(
            redeclare FCSys.Characteristics.N2.gas Data(
              b_v=[1],
              specVolPow={-1,0},
              specHeatCapPow=0,
              T_lim_c0={0,Modelica.Constants.inf},
              b_c0=[1.041e3*U.J*Data.m/(U.kg*U.K)],
              B_c0=[-298.15*U.K*1.041e3*U.J*Data.m/(U.kg*U.K) + Data.Deltah0_f,
                  0]),
            redeclare parameter Q.Resistivity alpha_N=Data.beta(),
            redeclare parameter Q.Resistivity alpha_Phi=Data.beta(),
            redeclare parameter Q.Resistivity alpha_S=U.m*U.K/(25.9e-3*U.W),
            redeclare parameter Q.Resistivity beta_N=Data.beta(),
            redeclare parameter Q.Resistivity beta_Phi=Data.m/(178.2e-7*U.Pa*U.s),

            redeclare parameter Q.Resistivity beta_S=U.m*U.K/(25.9e-3*U.W));

          // See the documentation for a table of values.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="N2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
    <li>Fixed specific heat capacity (independent of temperature)</ol>
    <li>Transport and exchange properties (&alpha;<sub>&Phi;</sub>, &beta;<sub>&Phi;</sub>, etc.) are fixed (e.g., independent of temperature)</li></p>

<p>The default specific heat capacity (<code>b_c0=[1.041e3*U.J*Data.m/(U.kg*U.K)]</code>) and transport resistivities (<code>beta_Phi=Data.m/(178.2e-7*U.Pa*U.s)</code> and <code>beta_S=U.m*U.K/(25.9e-3*U.W))</code>) are based on data of gas at 1 atm and
  300 K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 920].
  Table 1 lists the properties at  other temperatures. Note that the value for specific heat capacity at constant pressure at
  800 K (<code>c=1.22e3*U.J*Data.m/(U.kg*U.K)</code>) seems unusual, but it matches the
  reference.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"bottom\"><b>Table 1:</b> Properties of N<sub>2</sub> gas at 1 atm (<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 920)</caption>
  <tr>
      <th valign=\"middle\"><code>T/U.K</code></th>
      <th width=1><code>c*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1><code>beta_Phi<br>*U.Pa*U.s/Data.m</code></th>
      <th width=1><code>beta_S*U.W<br>/(U.m*U.K)</code></th>
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
       <code>beta_Phi=Data.m*U.s/(178e-7*U.Pa)</code>
   (<a href=\"http://en.wikipedia.org/wiki/Viscosity\">http://en.wikipedia.org/wiki/Viscosity</a>).</p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"));
        end Fixed;
      end gas;
    end N2;

    package O2 "<html>O<sub>2</sub></html>"
      extends Modelica.Icons.Package;
      package gas "<html>O<sub>2</sub> gas</html>"
        extends Modelica.Icons.Package;
        model Calibrated
          "Ideal or empirical correlations, w/ adjustment factors"
          extends SpeciesNeutral(
            redeclare FCSys.Characteristics.O2.gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare Q.Resistivity alpha_N,
            redeclare Q.Resistivity alpha_Phi,
            redeclare Q.Resistivity alpha_S,
            redeclare Q.Resistivity beta_N,
            redeclare Q.Resistivity beta_Phi,
            redeclare Q.Resistivity beta_S);

          parameter Q.NumberAbsolute k_alpha_N(final nominal=1) = 1
            "<html>Adjustment factor for reaction (<i>k</i><sub>&alpha;<i>N</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_alpha_Phi(final nominal=1) = 1
            "<html>Adjustment factor for exchange of linear momentum (<i>k</i><sub>&alpha;&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_alpha_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal exchange (<i>k</i><sub>&alpha;<i>S</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_N(final nominal=1) = 1
            "<html>Adjustment factor for material transport (<i>k</i>&beta;<sub><i>N</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_Phi(final nominal=1) = 1
            "<html>Adjustment factor for transport of linear momentum (<i>k</i>&beta;<sub>&Phi;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_beta_S(final nominal=1) = 1
            "<html>Adjustment factor for thermal transport (<i>k</i>&beta;<sub><i>S</i></sub>)</html>"
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
        equation
          alpha_N = k_alpha_N*Data.beta(T);
          alpha_Phi = k_alpha_Phi*Data.beta(T);
          alpha_S = k_alpha_S*Data.beta(T);
          beta_N = k_beta_N*Data.beta(T);

          beta_Phi = k_beta_Phi*(if empiricalTransverse then Data.beta_Phi(T)
             else Data.beta(T));
          beta_S = k_beta_S*(if empiricalThermal then Data.beta_S(T) else
            Data.beta(T));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="O2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"));
        end Calibrated;

        model Correlated "Ideal (kinetic) or empirical correlations"
          extends SpeciesNeutral(
            redeclare FCSys.Characteristics.O2.gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare Q.Resistivity alpha_N,
            redeclare Q.Resistivity alpha_Phi,
            redeclare Q.Resistivity alpha_S,
            redeclare Q.Resistivity beta_N,
            redeclare Q.Resistivity beta_Phi,
            redeclare Q.Resistivity beta_S);

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
        equation
          alpha_N = Data.beta(T);
          alpha_Phi = Data.beta(T);
          alpha_S = Data.beta(T);
          beta_N = Data.beta(T);

          beta_Phi = (if empiricalTransverse then Data.beta_Phi(T) else
            Data.beta(T));
          beta_S = (if empiricalThermal then Data.beta_S(T) else Data.beta(T));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="O2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol>
          </p><p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpeciesAmagat\">PartialSpeciesAmagat</a> model.</p></html>"));
        end Correlated;

        model Fixed "Fixed properties"
          extends SpeciesNeutral(
            redeclare FCSys.Characteristics.O2.gas Data(b_v=[1], specVolPow={-1,
                  0}),
            redeclare parameter Q.Resistivity alpha_N=1e10*Data.beta(),
            redeclare parameter Q.Resistivity alpha_Phi=Data.beta(),
            redeclare parameter Q.Resistivity alpha_S=U.m*U.K/(26.8e-3*U.W),
            redeclare parameter Q.Resistivity beta_N=Data.beta(),
            redeclare parameter Q.Resistivity beta_Phi=Data.m/(207.2e-7*U.Pa*U.s),

            redeclare parameter Q.Resistivity beta_S=U.m*U.K/(26.8e-3*U.W));

          // See the documentation for a table of values.

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="O2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
    <li>Transport and exchange properties (&alpha;<sub>&Phi;</sub>, &beta;<sub>&Phi;</sub>, etc.) are fixed (e.g., independent of temperature)</li>
    </ol></p>

<p>Additional notes:
<ul>
          <li>The default value for <code>alpha_N</code> (<code>1e10*Data.beta()</code>)
          gives reasonable reaction rates for the HOR and ORR when used in conjunction
          with the defaults for e<sup>-</sup>, H<sub>2</sub>, H<sub>2</sub>O, and H<sup>+</sup>.
          The default <code>alpha_N</code> is greater (higher resistance) for H<sub>2</sub>O
          and O<sub>2</sub> to slow the ORR.</li>
          <li>
  The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

  <p>The default transport resistivities (<code>beta_Phi=Data.m/(207.2e-7*U.Pa*U.s)</code> and <code>beta_S=U.m*U.K/(26.8e-3*U.W)</code>) are based on data of gas at 1 atm and
  300 K from Incropera and DeWitt  [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 920&ndash;921].
  Table 1 lists the properties at other temperatures.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"bottom\"><b>Table 1:</b> Properties of O<sub>2</sub> gas at 1 atm [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002, pp. 920&ndash;921</a>]</caption>
  <tr>
      <th valign=\"middle\"><code>T/U.K</code></th>
      <th width=1><code>c*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1><code>beta_Phi<br>*U.Pa*U.s/Data.m</code></th>
      <th width=1><code>beta_S*U.W<br>/(U.m*U.K)</code></th>
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
      end gas;
    end O2;

    model SpeciesInertStagnant "Inert and stagnant species model"
      extends Species(
        Data(final p_min=0),
        redeclare final constant Q.Resistivity alpha_N=Modelica.Constants.inf,
        redeclare final constant Q.Resistivity beta_N=Data.beta(),
        redeclare final constant Q.Resistivity beta_Phi=Data.beta(),
        negativeX(
          matEntOpt=MaterialEntropyOpt.ClosedDiabatic,
          final viscousY=false,
          final viscousZ=false),
        positiveX(
          matEntOpt=MaterialEntropyOpt.ClosedDiabatic,
          final viscousY=false,
          final viscousZ=false),
        negativeY(final viscousZ=false,final viscousX=false),
        positiveY(final viscousZ=false,final viscousX=false),
        negativeZ(final viscousX=false,final viscousY=false),
        positiveZ(final viscousX=false,final viscousY=false),
        final termDepleted=false,
        final upstreamX=false,
        final upstreamY=false,
        final upstreamZ=false,
        final setPartNum=true,
        final setXVel=true,
        final setYVel=true,
        final setZVel=true,
        partNumInitMeth=InitMethScalar.Volume,
        final xInitMeth=InitMethLinear.Velocity,
        final yInitMeth=InitMethLinear.Velocity,
        final zInitMeth=InitMethLinear.Velocity,
        final phi_IC=zeros(3),
        final derphi_IC=zeros(3),
        final I_IC=zeros(3),
        final derI_IC=zeros(3));

      // Note:  beta_N and beta_Phi don't matter since material and linear
      // momentum isn't transported.

      annotation (Documentation(info="<html><p>Assumptions:<ol>
  <li>Zero velocity</li>
  <li>No material exchange or transport</li</ol>
  </p>

  <p>For more information, see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model.</p></html>"));
    end SpeciesInertStagnant;

    model SpeciesNeutral "Partial electrically neutral species model"
      import FCSys;
      extends Species(final Deltamu_IC=0, redeclare final parameter
          Q.Permittivity epsilon=0);
      // Note:  These variables don't apply since the species is electrically
      // neutral.

    initial equation
      assert(Data.z == 0,
        "The PartialSpeciesNeutral model is being used with a species that is not electrically neutral.");
    end SpeciesNeutral;

    model Species
      "Model for single-species exchange, transport, and storage of material, linear momentum, and entropy"
      //extends FCSys.BaseClasses.Icons.Names.Middle;

      // Geometric parameters
      outer parameter Q.Length L[3](each final min=Modelica.Constants.small)
        "Length" annotation (HideResult=true,missingInnerMessage=
            "This model should be used within the Subregion model.");
      outer parameter Q.Area A[3] "Cross-sectional area" annotation (HideResult
          =true, missingInnerMessage=
            "This model should be used within the Subregion model.");
      parameter Q.NumberAbsolute k[3](
        each min=Modelica.Constants.small,
        each final nominal=1) = {1,1,1}
        "<html>Anisotropic factor for transport (<b>k</b>)</html>"
        annotation (HideResult=true, Dialog(group="Geometry"));
      parameter Q.Length Lstar(
        min=0,
        nominal=1e5*U.m,
        start=1e7*product(L)^(1/3))
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
      parameter Boolean termDepleted=false
        "Terminate when material is depleted" annotation (
        Evaluate=true,
        Dialog(tab="Assumptions", compact=true),
        choices(__Dymola_checkBox=true));
      parameter Boolean inclLin[3]={true,false,false}
        "true, if each component of linear momentum is included"
        annotation (Evaluate=true,Dialog(tab="Assumptions"));
      final parameter Integer n_lin=countTrue(inclLin)
        "Number of components of linear momentum"
        annotation (Evaluate=true, HideResult=true);

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
      parameter FCSys.Subregions.BaseClasses.InitMethScalar partNumInitMeth=
          InitMethScalar.Pressure "Method of initializing the particle number"
        annotation (Evaluate=true, Dialog(tab="Initialization", group=
              "Scalar properties"));
      parameter FCSys.Subregions.BaseClasses.InitMethScalar tempInitMeth=
          InitMethScalar.Temperature "Method of initializing the temperature"
        annotation (Evaluate=true, Dialog(tab="Initialization", group=
              "Scalar properties"));
      parameter Q.ParticleNumber N_IC(start=V_IC/v_IC)
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
          enable=partNumInitMeth == 3 or tempInitMeth == 3));
      // Note:  Dymola 7.4 doesn't recognize enumerations in the dialog enable
      // option, e.g.,
      //     enable=partNumInitMeth == InitMethScalar.ParticleNumberRate.
      // Therefore, the values of the enumerations are specified numerically for
      // this initial condition and some others below.
      parameter Q.Volume V_IC(start=product(L))
        "<html>Initial volume (<i>V</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Scalar properties"));
      parameter Q.VolumeRate derV_IC=0
        "<html>Initial rate of volume ((&part;<i>V</i>/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Scalar properties",
          enable=partNumInitMeth == 5 or tempInitMeth == 5));
      parameter Q.VolumeSpecific v_IC(min=Modelica.Constants.small, start=
            Data.v_pT(p_IC, T_IC))
        "<html>Initial specific volume (<i>v</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Scalar properties"));
      parameter Q.VolumeSpecificRate derv_IC=0
        "<html>Initial rate of specific volume ((&part;<i>v</i>/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Scalar properties",
          enable=partNumInitMeth == 7 or tempInitMeth == 7));
      parameter Q.PressureAbsolute p_IC(start=defaults.p)
        "<html>Initial pressure (<i>p</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Scalar properties"));
      parameter Q.PressureRate derp_IC=0
        "<html>Initial rate of pressure ((&part;<i>p</i>/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Scalar properties",
          enable=partNumInitMeth == 9 or tempInitMeth == 9));
      parameter Q.PotentialAbsolute T_IC(
        nominal=298.15*U.K,
        displayUnit="K",
        start=defaults.T)
        "<html>Initial temperature (<i>T</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Scalar properties"));
      parameter Q.PotentialRate derT_IC(displayUnit="K/s") = 0
        "<html>Initial rate of temperature (&part;<i>T</i>/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Scalar properties",
          enable=partNumInitMeth == 11 or tempInitMeth == 11));
      parameter Q.NumberAbsolute s_IC(min=Modelica.Constants.small, start=
            Data.s_pT(p_IC, T_IC))
        "<html>Initial specific entropy (<i>s</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization", group="Scalar properties"));
      parameter Q.NumberRate ders_IC=0
        "<html>Initial rate of specific entropy ((&part;<i>s</i>/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Scalar properties",
          enable=partNumInitMeth == 13 or tempInitMeth == 13));
      parameter Q.Potential h_IC(start=Data.h0_T(T_IC))
        "<html>Initial specific enthalpy (<i>h</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization", group="Scalar properties"));
      parameter Q.PotentialRate derh_IC=0
        "<html>Initial rate of specific enthalpy ((&part;<i>h</i>/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Scalar properties",
          enable=partNumInitMeth == 15 or tempInitMeth == 15));
      parameter Q.Potential mu_IC(start=Data.g_pT(p_IC, T_IC))
        "<html>Initial electrochemical potential (&mu;<sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization", group="Scalar properties"));
      parameter Q.PotentialRate dermu_IC=0
        "<html>Initial rate of electrochemical potential ((&part;&mu;/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Scalar properties",
          enable=partNumInitMeth == 17 or tempInitMeth == 17));
      parameter Q.Current Ndot_IC=0
        "<html>Initial reaction rate (<i>N&#775;</i><sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Scalar properties",
          enable=partNumInitMeth == 18 or tempInitMeth == 18));

      // Electrochemical initialization parameter
      parameter Q.Potential Deltamu_IC=0
        "<html>Initial electrochemical potential of reaction relative to bulk ((&Delta;&mu;)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Electrical",
          enable=Data.z <> 0));

      // Initialization parameters for linear momentum
      parameter FCSys.Subregions.BaseClasses.InitMethLinear xInitMeth=
          InitMethLinear.Velocity "Method of initializing the x-axis component"
        annotation (Dialog(
          tab="Initialization",
          group="Linear momentum",
          enable=inclLin[1]));
      parameter FCSys.Subregions.BaseClasses.InitMethLinear yInitMeth=
          InitMethLinear.Velocity "Method of initializing the y-axis component"
        annotation (Dialog(
          tab="Initialization",
          group="Linear momentum",
          enable=inclLin[2]));
      parameter FCSys.Subregions.BaseClasses.InitMethLinear zInitMeth=
          InitMethLinear.Velocity "Method of initializing the z-axis component"
        annotation (Dialog(
          tab="Initialization",
          group="Linear momentum",
          enable=inclLin[3]));
      // Note:  Dymola 7.4 doesn't provide pull-down lists for arrays of
      // enumerations; therefore, a parameter is used for each axis.
      parameter Q.Velocity phi_IC[3]={0,0,0}
        "<html>Initial velocity (<b>&phi;</b><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Linear momentum"));
      parameter Q.Acceleration derphi_IC[3]={0,0,0}
        "<html>Initial acceleration ((&part;<b>&phi;</b>/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Linear momentum",
          enable=xInitMeth == 3 or yInitMeth == 3 or zInitMeth == 3));
      parameter Q.Current I_IC[3]={0,0,0}
        "<html>Initial current (<i><b>I</b></i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Linear momentum"));
      parameter Q.CurrentRate derI_IC[3]={0,0,0}
        "<html>Initial rate of current ((&part;<i><b>I</b></i>/&part;<i>t</i>)<sub>IC</sub>)</html>"
        annotation (Dialog(
          tab="Initialization",
          group="Linear momentum",
          enable=xInitMeth == 5 or yInitMeth == 5 or zInitMeth == 5));

      // Preferred states
      Q.ParticleNumber N(
        nominal=1*U.mol,
        min=Modelica.Constants.small,
        final start=N_IC,
        final fixed=false,
        stateSelect=StateSelect.prefer) "Particle number";
      // Note:  The start value for this variable (and others below) isn't
      // fixed because the related initial condition is applied in the initial
      // equation section.
      Q.Velocity phi[n_lin](
        each nominal=1*U.cm/U.s,
        final start=phi_IC[cartAxes],
        each final fixed=false,
        each stateSelect=StateSelect.prefer) "Velocity";
      Q.PotentialAbsolute T(
        nominal=298.15*U.K,
        final start=T_IC,
        final fixed=false,
        stateSelect=StateSelect.prefer,
        displayUnit="K") "Temperature";

      // Aliases (for common terms)
      Q.VolumeSpecific v(
        nominal=Data.v_pT(1*U.atm, 298.15*U.K),
        final start=v_IC,
        final fixed=false) "Specific volume";
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
        each nominal=100*U.A,
        final start=I_IC[cartAxes],
        each final fixed=false) "Current";
      Q.Potential ke(final start=if n_lin == 0 then 0 else Data.m*phi_IC[
            cartAxes]*phi_IC[cartAxes]/2) "Specific macroscopic kinetic energy";
      Q.Potential Deltamu(final start=Deltamu_IC)
        "Electrochemical potential of exchange relative to bulk";
      // Note:  The initial condition isn't fixed because the state
      // only appears if the species is charged and it is involved in a
      // reaction.

      // Material properties
      input Q.Resistivity alpha_N(nominal=10*U.cm/U.A, start=Data.beta(T_IC))
        "<html>Reaction resistivity (&alpha;<sub><i>N</i></sub>)</html>"
        annotation (Dialog(group="Material properties"));
      input Q.Resistivity alpha_Phi(nominal=10*U.cm/U.A, start=Data.beta(T_IC))
        "<html>Exchange resistivity for linear momentum (&alpha;<sub>&Phi;</sub>)</html>"
        annotation (Dialog(group="Material properties"));
      input Q.Resistivity alpha_S(nominal=10*U.cm/U.A, start=Data.beta(T_IC))
        "<html>Thermal exchange resistivity (&alpha;<sub><i>S</i></sub>)</html>"
        annotation (Dialog(group="Material properties"));
      input Q.Resistivity beta_N(nominal=10*U.cm/U.A, start=Data.beta(T_IC))
        "<html>Material transport resistivity (&beta;<sub><i>N</i></sub>)</html>"
        annotation (Dialog(group="Material properties"));
      input Q.Resistivity beta_Phi(nominal=10*U.cm/U.A, start=Data.beta(T_IC))
        "<html>Transport resistivity for linear momentum (&beta;<sub>&Phi;</sub>)</html>"
        annotation (Dialog(group="Material properties"));
      input Q.Resistivity beta_S(nominal=10*U.cm/U.A, start=Data.beta(T_IC))
        "<html>Thermal transport resistivity (&beta;<sub><i>S</i></sub>)</html>"
        annotation (Dialog(group="Material properties"));
      input Q.Permittivity epsilon(nominal=1e6*U.epsilon_0, start=1e6*U.epsilon_0)
        "<html>Permittivity (if charged; &epsilon;)</html>"
        annotation (Dialog(group="Material properties", enable=Data.z <> 0));
      // Note:  The properties are defined as inputs so that they can be
      // redeclared as parameters or acausal time-varying variables.

      output Q.Time t=beta_Phi*Lstar/N if defaults.analysis "**";

      // Nonessential variables (for analysis)
      // General
      output Q.Mass M=Data.m*N if defaults.analysis "Mass";
      output Q.Volume V=inert.V if defaults.analysis "Volume";
      output Q.PressureAbsolute p=inert.p if defaults.analysis "Pressure";
      output Q.ParticleNumberVolumic rho=N/V if defaults.analysis
        "Molar density";
      output Q.Capacitance C=Lstar*epsilon if defaults.analysis and (
        Modelica.Utilities.Strings.find(Data.formula, "+") <> 0 or
        Modelica.Utilities.Strings.find(Data.formula, "-") <> 0)
        "Double layer capacitance";

      // Note:  The string find functions are used instead of
      // Data.z <> 0 due to the following error in Dymola 7.4:
      //     "Current version of Dymola can only handle conditional
      //      components with fixed condition."
      // even though Data.z is calculated directly from constants.
      // Peclet numbers
      output Q.Time RC=C*alpha_N/Lstar if defaults.analysis and (
        Modelica.Utilities.Strings.find(Data.formula, "+") <> 0 or
        Modelica.Utilities.Strings.find(Data.formula, "-") <> 0)
        "**Electrochemical time constant";
      output Q.Number Pe_N[n_lin]=I*beta_N ./ Lstar_trans[cartAxes] if defaults.analysis
        "Material Peclet numbers";
      output Q.Number Pe_Phi[n_lin]=I*beta_Phi ./ Lstar_trans[cartAxes] if
        defaults.analysis "Peclet numbers for linear momentum";
      output Q.Number Pe_S[n_lin]=I*beta_S ./ Lstar_trans[cartAxes] if defaults.analysis
        "Thermal Peclet numbers";
      // Advection
      output Q.PressureAbsolute q[n_lin]=Data.m*phi .* I ./ A[cartAxes] if
        defaults.analysis "Bulk dynamic pressure";
      output Q.Force mphiI[n_lin, 2]={(if inclLin[cartWrap(cartAxes[axis] +
          orientation)] then Data.m*phi[linAxes[cartWrap(cartAxes[axis] +
          orientation)]]*I[axis] else 0) for orientation in 1:2, axis in 1:
          n_lin} if n_lin > 0 and defaults.analysis
        "Bulk rate of advection of 1st and 2nd transverse linear momentum";
      output Q.Force sI[n_lin]=s*I if defaults.analysis
        "Bulk rate of advection of entropy";
      // Linear momentum balance
      output Q.Force Mderphi[n_lin]=M*der(phi)/U.s if defaults.analysis
        "Rate of storage of linear momentum at constant mass";
      output Q.Force mPhidot_exch_adv[n_lin]=chemical.mPhidot - Data.m*phi*
          chemical.Ndot if defaults.analysis
        "Acceleration force due to chemical exchange";
      output Q.Force mPhidot_exch_diff[n_lin]=inert.mPhidot
        "Force due to friction with other species";
      output Q.Force mPhidot_parallel[n_lin]={sum(if [negativeX.matEntOpt ==
          MaterialEntropyOpt.OpenDiabatic, positiveX.matEntOpt ==
          MaterialEntropyOpt.OpenDiabatic; negativeY.matEntOpt ==
          MaterialEntropyOpt.OpenDiabatic, positiveY.matEntOpt ==
          MaterialEntropyOpt.OpenDiabatic; negativeZ.matEntOpt ==
          MaterialEntropyOpt.OpenDiabatic, positiveZ.matEntOpt ==
          MaterialEntropyOpt.OpenDiabatic][cartAxes[axis], side] then (3 - 2*
          side)*(mu_face[cartAxes[axis], side] - mu + Data.Deltah0_f)*N/L[
          cartAxes[axis]] else 0 for side in 1:2) for axis in 1:n_lin} - Data.m
          *sum(Ndot_face)*phi if defaults.analysis
        "Acceleration force through normal faces";
      output Q.Force mPhidot_perpendicular[n_lin]={sum(Sigma(mPhidot_face[
          cartWrap(cartAxes[axis] - orientation), :, orientation]) - Data.m*phi[
          axis]*Sigma(Ndot_face[cartWrap(cartAxes[axis] - orientation), :])
          for orientation in 1:2) for axis in 1:n_lin} if defaults.analysis
        "Acceleration force through transverse faces";
      // Energy balance
      output Q.Power derE=(N*(Data.c0_T(T)*der(T) + der(ke)) + inert.V*der(
          Data.p_vT(inert.V/N, T)))/U.s if defaults.analysis
        "Rate of energy storage at constant particle number";
      output Q.Power Edot_exch_N=(chemical.mu - h - ke)*chemical.Ndot +
          semiLinear(
              chemical.Ndot,
              chemical.mphi*chemical.mphi/Data.m,
              Data.m*phi*phi) + chemical.TSdot if defaults.analysis
        "Rate of energy intake due to chemical exchange";
      output Q.Power Edot_exch_Phi=inert.phi*inert.mPhidot if defaults.analysis
        "Rate of energy intake due to diffusion of linear momentum from other species";
      output Q.Power Qdot_exch=inert.T*inert.Sdot if defaults.analysis
        "Rate of thermal conduction from other species";
      output Q.Power Edot_trans_N=sum((mu_face - fill(
              h + ke - Data.Deltah0_f,
              3,
              2)) .* Ndot_face) if defaults.analysis
        "Rate of energy intake due to material transport";
      output Q.Power Edot_trans_Phi=sum(phi_face .* mPhidot_face) if defaults.analysis
        "Rate of energy intake due to shear force";
      output Q.Power Qdot_trans=sum(T_face .* Sdot_face) if defaults.analysis
        "Rate of thermal convection from other subregions";
      Connectors.ChemicalOutput chemical(
        final n_lin=n_lin,
        final formula=Data.formula,
        mu(start=mu_IC),
        mphi(final start=Data.m*phi_IC[cartAxes]),
        Ndot(final start=Ndot_IC,final fixed=false))
        "Connector to exchange material with advection of linear momentum"
        annotation (Placement(transformation(extent={{-30,10},{-10,30}}),
            iconTransformation(extent={{-80,60},{-60,80}})));
      FCSys.Connectors.InertDalton inert(
        final n_lin=n_lin,
        V(
          min=0,
          final start=V_IC,
          final fixed=false),
        p(final start=p_IC,final fixed=false),
        phi(start=phi_IC[cartAxes]),
        T(start=T_IC))
        "Connector to add pressure and exchange linear momentum and entropy by diffusion"
        annotation (Placement(transformation(extent={{10,-30},{30,-10}}),
            iconTransformation(extent={{60,-80},{80,-60}})));
      replaceable FCSys.Connectors.FaceX negativeX(
        matEntOpt=MaterialEntropyOpt.OpenDiabatic,
        viscousY=inclLin[2],
        viscousZ=inclLin[3],
        material(final mu(start=mu_IC - Data.Deltah0_f) = mu_face[1, 1], final
            Ndot(start=I_IC[1]) = Ndot_face[1, 1]),
        momentumY(final phi(start=phi_IC[2]) = phi_face[1, 1, 1],final mPhidot=
              mPhidot_face[1, 1, 1]),
        momentumZ(final phi(start=phi_IC[3]) = phi_face[1, 1, 2],final mPhidot=
              mPhidot_face[1, 1, 2]),
        entropy(final T(start=T_IC) = T_face[1, 1],final Sdot=Sdot_face[1, 1]))
        "Negative face along the x axis" annotation (Dialog(
          tab="Assumptions",
          group="Characteristics of the faces (click to edit)",
          __Dymola_label="negativeX",
          __Dymola_descriptionLabel=true),Placement(transformation(extent={{-50,
                -10},{-30,10}}), iconTransformation(extent={{-110,-10},{-90,10}})));
      replaceable FCSys.Connectors.FaceX positiveX(
        matEntOpt=negativeX.matEntOpt,
        viscousY=inclLin[2],
        viscousZ=inclLin[3],
        material(final mu(start=mu_IC - Data.Deltah0_f) = mu_face[1, 2], final
            Ndot(start=-I_IC[1]) = Ndot_face[1, 2]),
        momentumY(final phi(start=phi_IC[2]) = phi_face[1, 2, 1],final mPhidot=
              mPhidot_face[1, 2, 1]),
        momentumZ(final phi(start=phi_IC[3]) = phi_face[1, 2, 2],final mPhidot=
              mPhidot_face[1, 2, 2]),
        entropy(final T(start=T_IC) = T_face[1, 2],final Sdot=Sdot_face[1, 2]))
        "Positive face along the x axis" annotation (Dialog(
          tab="Assumptions",
          group="Characteristics of the faces (click to edit)",
          __Dymola_label="positiveX",
          __Dymola_descriptionLabel=true),Placement(transformation(extent={{30,
                -10},{50,10}}), iconTransformation(extent={{90,-10},{110,10}})));

      replaceable FCSys.Connectors.FaceY negativeY(
        viscousZ=inclLin[3],
        viscousX=inclLin[1],
        material(final mu(start=mu_IC - Data.Deltah0_f) = mu_face[2, 1], final
            Ndot(start=I_IC[2]) = Ndot_face[2, 1]),
        momentumZ(final phi(start=phi_IC[3]) = phi_face[2, 1, 1],final mPhidot=
              mPhidot_face[2, 1, 1]),
        momentumX(final phi(start=phi_IC[1]) = phi_face[2, 1, 2],final mPhidot=
              mPhidot_face[2, 1, 2]),
        entropy(final T(start=T_IC) = T_face[2, 1],final Sdot=Sdot_face[2, 1]))
        "Negative face along the y axis" annotation (Dialog(
          tab="Assumptions",
          group="Characteristics of the faces (click to edit)",
          __Dymola_label="negativeY",
          __Dymola_descriptionLabel=true),Placement(transformation(extent={{-10,
                -50},{10,-30}}), iconTransformation(extent={{-10,-110},{10,-90}})));
      replaceable FCSys.Connectors.FaceY positiveY(
        matEntOpt=negativeY.matEntOpt,
        viscousZ=inclLin[3],
        viscousX=inclLin[1],
        material(final mu(start=mu_IC - Data.Deltah0_f) = mu_face[2, 2], final
            Ndot(start=-I_IC[2]) = Ndot_face[2, 2]),
        momentumZ(final phi(start=phi_IC[3]) = phi_face[2, 2, 1],final mPhidot=
              mPhidot_face[2, 2, 1]),
        momentumX(final phi(start=phi_IC[1]) = phi_face[2, 2, 2],final mPhidot=
              mPhidot_face[2, 2, 2]),
        entropy(final T(start=T_IC) = T_face[2, 2],final Sdot=Sdot_face[2, 2]))
        "Positive face along the y axis" annotation (Dialog(
          tab="Assumptions",
          group="Characteristics of the faces (click to edit)",
          __Dymola_label="positiveY",
          __Dymola_descriptionLabel=true),Placement(transformation(extent={{-10,
                30},{10,50}}), iconTransformation(extent={{-10,90},{10,110}})));

      replaceable FCSys.Connectors.FaceZ negativeZ(
        viscousX=inclLin[1],
        viscousY=inclLin[2],
        material(final mu(start=mu_IC - Data.Deltah0_f) = mu_face[3, 1], final
            Ndot(start=I_IC[3]) = Ndot_face[3, 1]),
        momentumX(final phi(start=phi_IC[1]) = phi_face[3, 1, 1],final mPhidot=
              mPhidot_face[3, 1, 1]),
        momentumY(final phi(start=phi_IC[2]) = phi_face[3, 1, 2],final mPhidot=
              mPhidot_face[3, 1, 2]),
        entropy(final T(start=T_IC) = T_face[3, 1],final Sdot=Sdot_face[3, 1]))
        "Negative face along the z axis" annotation (Dialog(
          tab="Assumptions",
          group="Characteristics of the faces (click to edit)",
          __Dymola_label="negativeZ",
          __Dymola_descriptionLabel=true),Placement(transformation(extent={{10,
                10},{30,30}}), iconTransformation(extent={{60,60},{80,80}})));
      replaceable FCSys.Connectors.FaceZ positiveZ(
        matEntOpt=negativeZ.matEntOpt,
        viscousX=inclLin[1],
        viscousY=inclLin[2],
        material(final mu(start=mu_IC - Data.Deltah0_f) = mu_face[3, 2], final
            Ndot(start=-I_IC[3]) = Ndot_face[3, 2]),
        momentumX(final phi(start=phi_IC[1]) = phi_face[3, 2, 1],final mPhidot=
              mPhidot_face[3, 2, 1]),
        momentumY(final phi(start=phi_IC[2]) = phi_face[3, 2, 2],final mPhidot=
              mPhidot_face[3, 2, 2]),
        entropy(final T(start=T_IC) = T_face[3, 2],final Sdot=Sdot_face[3, 2]))
        "Positive face along the z axis" annotation (Dialog(
          tab="Assumptions",
          group="Characteristics of the faces (click to edit)",
          __Dymola_label="positiveZ",
          __Dymola_descriptionLabel=true),Placement(transformation(extent={{-30,
                -30},{-10,-10}}), iconTransformation(extent={{-80,-80},{-60,-60}})));
      // Note:  These connectors are replaceable so that their parameters can be
      // edited directly in the parameter dialog.

      // Geometric parameters
    protected
      final parameter Integer cartAxes[n_lin]=index(inclLin)
        "Cartesian-axis indices of the axes of linear momentum";
      final parameter Integer linAxes[3]=enumerate(inclLin)
        "Linear momentum indices of the Cartesian axes";
      final parameter Boolean upstream[3]={upstreamX,upstreamY,upstreamZ}
        "true, if each Cartesian axis uses upstream discretization";
      final parameter Q.Length Lstar_trans[3]=k .* A ./ L
        "Effective cross-sectional area per length";

      Boolean depleted "true, if nearly no material remains";

      // Efforts and flows of the conditional faces
      Q.Potential mu_face[3, 2](each start=mu_IC - Data.Deltah0_f)
        "Electrochemical potentials at the faces";
      Q.Current Ndot_face[3, 2](start=outerProduct(I_IC, {1,-1}))
        "Currents into the faces";
      Q.Velocity phi_face[3, 2, 2](start={fill({phi_IC[cartWrap(axis +
            orientation)] for orientation in 1:2}, 2) for axis in 1:3})
        "Transverse velocities at the faces";
      Q.Force mPhidot_face[3, 2, 2] "Transverse forces on the faces";
      Q.Potential T_face[3, 2](each start=T_IC,each displayUnit="K")
        "Temperatures at the faces";
      Q.Current Sdot_face[3, 2] "Entropy flow rates into the faces";

      outer FCSys.BCs.Defaults defaults "Default settings" annotation (
          missingInnerMessage="Your model is using an outer \"defaults\" record, but an inner \"defaults\" record is not defined.
For simulation, specify global default settings by dragging FCSys.BCs.Defaults into your model.
The default global default settings will be used for the current simulation.",
          Placement(transformation(extent={{40,40},{60,60}}),
            iconTransformation(extent={{-10,90},{10,110}})));
      // Note:  In Dymola 7.4, it is necessary to add the missing inner message
      // here to give a warning message, even though it is included in the Defaults
      // model too.
    initial equation
      // Check that the initialization methods are valid.
      assert(partNumInitMeth <> tempInitMeth or partNumInitMeth ==
        InitMethScalar.None,
        "The initialization methods for particle number and temperature cannot be the same (unless None).");
      if not Data.isCompressible then
        assert(partNumInitMeth <> InitMethScalar.Pressure and partNumInitMeth
           <> InitMethScalar.PressureRate or setPartNum, "The material is incompressible,
      yet the initialization method for particle number involves pressure.");
        assert(tempInitMeth <> InitMethScalar.Pressure and tempInitMeth <>
          InitMethScalar.PressureRate or setTemp, "The material is incompressible,
      yet the initialization method for temperature involves pressure.");
        assert(Data.hasThermalExpansion or (partNumInitMeth <> InitMethScalar.VolumeSpecific
           and partNumInitMeth <> InitMethScalar.VolumeSpecificRate or
          setPartNum), "The material has constant specific volume,
      yet the initialization method for particle number involves specific volume.");
        assert(Data.hasThermalExpansion or (tempInitMeth <> InitMethScalar.VolumeSpecific
           and tempInitMeth <> InitMethScalar.VolumeSpecificRate or setPartNum),
          "The material has constant specific volume,
      yet the initialization method for temperature involves specific volume.");
      end if;

      /* This is commented out because it may be annoying.
  // Warn when index reduction may be necessary.
  if abs(alpha_N) < Modelica.Constants.small then
    Modelica.Utilities.Streams.print("Warning: The material resistivity to exchange is zero.
    This may directly couple the chemical potentials of species reacting within a subregion.
    Consider setting the value of alpha_N as final (if not already) so that index reduction may be performed.");
  end if;
  if abs(alpha_Phi) > Modelica.Constants.small then
    Modelica.Utilities.Streams.print("Warning: The resistivity to exchange of linear momentum is zero.
    This may directly couple the velocities of species within a subregion.
    Consider setting the value of alpha_Phi as final (if not already) so that index reduction may be performed.");
  end if;
  if abs(alpha_S) > Modelica.Constants.small then
    Modelica.Utilities.Streams.print("Warning: The thermal resistivity to exchange is zero.
    This may directly couple the temperatures of species within a subregion.
    Consider setting the value of alpha_S as final (if not already) so that index reduction may be performed.");
  end if;
  if abs(beta_N) < Modelica.Constants.small then
    Modelica.Utilities.Streams.print("Warning: The material resistivity to transport is zero.
    This may directly couple the density within neighboring subregions.\nConsider setting the value of beta_N as final (if not already) so that index reduction may be performed.");
  end if;
  if abs(beta_Phi) > Modelica.Constants.small then
    Modelica.Utilities.Streams.print("Warning: The resistivity to transport of linear momentum is zero.
    This may directly couple the velocity within neighboring subregions.\nConsider setting the value of beta_Phi as final (if not already) so that index reduction may be performed.");
  end if;
  if abs(beta_S) > Modelica.Constants.small then
    Modelica.Utilities.Streams.print("Warning: The thermal resistivity to transport is zero.
    This may directly couple the temperature within neighboring subregions.\nConsider setting the value of beta_S as final (if not already) so that index reduction may be performed.");
  end if;
  // Note:  According to the Modelica 3.0 specification (and later), these
  // checks should be possible using the assert() command with
  // level=AssertionLevel.warning.  However, this isn't supported in
  // Dymola 7.4.
  */

      // Material
      if setPartNum then
        // Ensure that a condition is selected, since the state is prescribed.
        assert(partNumInitMeth <> InitMethScalar.None, "The state for particle number is prescribed,
    yet its condition is not defined.\nChoose a condition besides None.");
      else
        // Initialize only if there is a time-varying state.
        if partNumInitMeth == InitMethScalar.ParticleNumber then
          N = N_IC;
        elseif partNumInitMeth == InitMethScalar.ParticleNumberRate then
          der(N)/U.s = derN_IC;
        elseif partNumInitMeth == InitMethScalar.Volume then
          inert.V = V_IC;
        elseif partNumInitMeth == InitMethScalar.VolumeRate then
          der(inert.V)/U.s = derV_IC;
        elseif partNumInitMeth == InitMethScalar.VolumeSpecific then
          v = v_IC;
        elseif partNumInitMeth == InitMethScalar.VolumeSpecificRate then
          der(v)/U.s = derv_IC;
        elseif partNumInitMeth == InitMethScalar.Pressure then
          inert.p = p_IC;
        elseif partNumInitMeth == InitMethScalar.PressureRate then
          der(inert.p)/U.s = derp_IC;
        elseif partNumInitMeth == InitMethScalar.Temperature then
          T = T_IC;
        elseif partNumInitMeth == InitMethScalar.TemperatureRate then
          der(T)/U.s = derT_IC;
        elseif partNumInitMeth == InitMethScalar.SpecificEntropy then
          s = s_IC;
        elseif partNumInitMeth == InitMethScalar.SpecificEntropyRate then
          der(s)/U.s = ders_IC;
        elseif partNumInitMeth == InitMethScalar.SpecificEnthalpy then
          h = h_IC;
        elseif partNumInitMeth == InitMethScalar.SpecificEnthalpyRate then
          der(h)/U.s = derh_IC;
        elseif partNumInitMeth == InitMethScalar.PotentialElectrochemical then
          mu = mu_IC;
        elseif partNumInitMeth == InitMethScalar.PotentialElectrochemicalRate
             then
          der(mu) = dermu_IC;
        elseif partNumInitMeth == InitMethScalar.ReactionRate then
          chemical.Ndot = Ndot_IC;
          // Else, partNumInitMeth == InitMethScalar.None; then, there are no
          // initial equations.
        end if;
      end if;

      // Temperature
      if setTemp then
        // Ensure that a condition is selected, since the state is prescribed.
        assert(tempInitMeth <> InitMethScalar.None, "The state for temperature is prescribed,
    yet its condition is not defined.\nChoose a condition besides None.");
      else
        // Initialize only if there is a time-varying state.
        if tempInitMeth == InitMethScalar.ParticleNumber then
          N = N_IC;
        elseif tempInitMeth == InitMethScalar.ParticleNumberRate then
          der(N)/U.s = derN_IC;
        elseif tempInitMeth == InitMethScalar.Volume then
          inert.V = V_IC;
        elseif tempInitMeth == InitMethScalar.VolumeRate then
          der(inert.V)/U.s = derV_IC;
        elseif tempInitMeth == InitMethScalar.VolumeSpecific then
          v = v_IC;
        elseif tempInitMeth == InitMethScalar.VolumeSpecificRate then
          der(v)/U.s = derv_IC;
        elseif tempInitMeth == InitMethScalar.Pressure then
          inert.p = p_IC;
        elseif tempInitMeth == InitMethScalar.PressureRate then
          der(inert.p)/U.s = derp_IC;
        elseif tempInitMeth == InitMethScalar.Temperature then
          T = T_IC;
        elseif tempInitMeth == InitMethScalar.TemperatureRate then
          der(T)/U.s = derT_IC;
        elseif tempInitMeth == InitMethScalar.SpecificEntropy then
          s = s_IC;
        elseif tempInitMeth == InitMethScalar.SpecificEntropyRate then
          der(s)/U.s = ders_IC;
        elseif tempInitMeth == InitMethScalar.SpecificEnthalpy then
          h = h_IC;
        elseif tempInitMeth == InitMethScalar.SpecificEnthalpyRate then
          der(h)/U.s = derh_IC;
        elseif tempInitMeth == InitMethScalar.PotentialElectrochemical then
          mu = mu_IC;
        elseif tempInitMeth == InitMethScalar.PotentialElectrochemicalRate then
          der(mu) = dermu_IC;
        elseif tempInitMeth == InitMethScalar.ReactionRate then
          chemical.Ndot = Ndot_IC;
          // Else, tempInitMeth == InitMethScalar.None; then, there are no
          // initial equations.
        end if;
      end if;

      // Linear momentum
      for axis in 1:n_lin loop
        if cartAxes[axis] == 1 then
          // X axis
          if setXVel then
            // Ensure that a condition is selected, since the state is
            // prescribed.
            assert(xInitMeth <> InitMethLinear.None, "The state for the x-axis component of linear momentum is prescribed,
        yet its condition is not defined.\nChoose any condition besides None.");
          else
            // Initialize only if there is a time-varying state.
            if xInitMeth == InitMethLinear.Velocity then
              phi[axis] = phi_IC[1];
            elseif xInitMeth == InitMethLinear.Acceleration then
              der(phi[axis])/U.s = derphi_IC[1];
            elseif xInitMeth == InitMethLinear.Current then
              I[axis] = I_IC[1];
            elseif xInitMeth == InitMethLinear.CurrentRate then
              der(I[axis])/U.s = derI_IC[1];
              // Else, xInitMeth == InitMethLinear.None; then, there are no initial
              // equations.
            end if;
          end if;
        elseif cartAxes[axis] == 2 then
          // Y axis
          if setYVel then
            // Ensure that a condition is selected, since the state is
            // prescribed.
            assert(yInitMeth <> InitMethLinear.None, "The state for the y-axis component of linear momentum is prescribed,
        yet its condition is not defined.\nChoose any condition besides None.");
          else
            // Initialize only if there is a time-varying state.
            if yInitMeth == InitMethLinear.Velocity then
              phi[axis] = phi_IC[2];
            elseif yInitMeth == InitMethLinear.Acceleration then
              der(phi[axis])/U.s = derphi_IC[2];
            elseif yInitMeth == InitMethLinear.Current then
              I[axis] = I_IC[2];
            elseif yInitMeth == InitMethLinear.CurrentRate then
              der(I[axis])/U.s = derI_IC[2];
              // Else, yInitMeth == InitMethLinear.None; then, there are no initial
              // equations.
            end if;
          end if;
        elseif cartAxes[axis] == 3 then
          // Z axis
          if setZVel then
            // Ensure that a condition is selected, since the state is
            // prescribed.
            assert(zInitMeth <> InitMethLinear.None, "The state for the z-axis component of linear momentum is prescribed,
        yet its condition is not defined.\nChoose any condition besides None.");
          else
            // Initialize only if there is a time-varying state.
            if zInitMeth == InitMethLinear.Velocity then
              phi[axis] = phi_IC[3];
            elseif zInitMeth == InitMethLinear.Acceleration then
              der(phi[axis])/U.s = derphi_IC[3];
            elseif zInitMeth == InitMethLinear.Current then
              I[axis] = I_IC[3];
            elseif zInitMeth == InitMethLinear.CurrentRate then
              der(I[axis])/U.s = derI_IC[3];
              // Else, zInitMeth == InitMethLinear.None; then, there are no initial
              // equations.
            end if;
          end if;
        end if;
      end for;

    equation
      // Empirical thermodynamic correlations
      if Data.isCompressible then
        inert.p = Data.p_vT(v, T);
      else
        v = Data.v_pT(inert.p, T);
      end if;
      h = Data.h0_T(T, ReferenceEnthalpy.EnthalpyOfFormationAt25degC);
      s = Data.s_pT(inert.p, T);

      // Aliases (for clarity and simplification)
      v*N = inert.V;
      N*phi = I .* L[cartAxes];
      //  mu = h - chemical.Ts;
      mu = h - T*s;
      Deltamu = chemical.mu - mu;
      2*ke = if n_lin == 0 then 0 else Data.m*phi*phi;

      // Protection against depletion
      depleted = (Data.p_min > 0 and inert.p <= Data.p_min) or N < 1*U.q;
      if termDepleted then
        when depleted then
          terminate("The " + Data.formula + "  species (" + Data.phase +
            " phase) is depleted (N = " + String(N/U.C) + " C and p = " +
            String(inert.p/U.kPa) + " kPa).");
        end when;
      end if;

      // Exchange
      /*
    alpha_N*T*(chemical.Ndot/(2*Lstar) - (if Data.z == 0 then 0 else epsilon*der(
    Deltamu)/U.s)) = if termDepleted or not depleted then Deltamu else 0
    "Material";

    **Is depletion protection necessary?  If not, remove it in Characteristic record.  If
    so, re-add this code, eliminate below.
  */
      alpha_N*T*(chemical.Ndot/(2*Lstar) - (if Data.z == 0 then 0 else epsilon*
        der(Deltamu)/U.s)) = Deltamu "Material";
      chemical.mPhidot = semiLinear(
            chemical.Ndot,
            chemical.mphi,
            Data.m*phi) "Advection of linear momentum";
      alpha_Phi*inert.mPhidot = 2*Lstar*Data.m*(inert.phi - phi)
        "Diffusion of linear momentum";
      chemical.TSdot = semiLinear(
            chemical.Ndot,
            chemical.Ts,
            T*s) "Thermal advection";
      alpha_S*inert.Sdot = 2*Lstar*(inert.T/T - 1) "Thermal diffusion";

      // Transport
      for axis in 1:3 loop
        for side in 1:2 loop
          // Material
          if [negativeX.matEntOpt == MaterialEntropyOpt.OpenDiabatic, positiveX.matEntOpt
               == MaterialEntropyOpt.OpenDiabatic; negativeY.matEntOpt ==
              MaterialEntropyOpt.OpenDiabatic, positiveY.matEntOpt ==
              MaterialEntropyOpt.OpenDiabatic; negativeZ.matEntOpt ==
              MaterialEntropyOpt.OpenDiabatic, positiveZ.matEntOpt ==
              MaterialEntropyOpt.OpenDiabatic][axis, side] then
            T*beta_N*(Ndot_face[axis, side] + (2*side - 3)*(if inclLin[axis]
               then I[linAxes[axis]] else 0)) = Lstar_trans[axis]*(mu_face[axis,
              side] - mu + Data.Deltah0_f)*(if upstream[axis] and inclLin[axis]
               then (exp((3 - 2*side)*I[linAxes[axis]]*beta_N/(2*Lstar_trans[
              axis])) + 1) else 2) "Advection and diffusion";
          else
            mu_face[axis, side] = 0;
            Ndot_face[axis, side] = 0;
          end if;

          // Linear momentum
          for orientation in 1:2 loop
            if {{{negativeX.viscousY,negativeX.viscousZ},{positiveX.viscousY,
                positiveX.viscousZ}},{{negativeY.viscousZ,negativeY.viscousX},{
                positiveY.viscousZ,positiveY.viscousX}},{{negativeZ.viscousX,
                negativeZ.viscousY},{positiveZ.viscousX,positiveZ.viscousY}}}[
                axis, side, orientation] then
              beta_Phi*(mPhidot_face[axis, side, orientation]/Data.m + (if
                inclLin[axis] then (2*side - 3)*I[linAxes[axis]] else 0)*(if
                inclLin[cartWrap(axis + orientation)] then phi[linAxes[cartWrap(
                axis + orientation)]] else 0)) = Lstar_trans[axis]*(phi_face[
                axis, side, orientation] - (if inclLin[cartWrap(axis +
                orientation)] then phi[linAxes[cartWrap(axis + orientation)]]
                 else 0))*(if upstream[axis] and inclLin[axis] then (exp((3 - 2
                *side)*I[linAxes[axis]]*beta_Phi/(2*Lstar_trans[axis])) + 1)
                 else 2) "Diffusion";
            else
              phi_face[axis, side, orientation] = 0;
              mPhidot_face[axis, side, orientation] = 0;
            end if;
          end for;

          // Entropy
          if [negativeX.matEntOpt == MaterialEntropyOpt.OpenDiabatic or
              negativeX.matEntOpt == MaterialEntropyOpt.ClosedDiabatic,
              positiveX.matEntOpt == MaterialEntropyOpt.OpenDiabatic or
              positiveX.matEntOpt == MaterialEntropyOpt.ClosedDiabatic;
              negativeY.matEntOpt == MaterialEntropyOpt.OpenDiabatic or
              negativeY.matEntOpt == MaterialEntropyOpt.ClosedDiabatic,
              positiveY.matEntOpt == MaterialEntropyOpt.OpenDiabatic or
              positiveY.matEntOpt == MaterialEntropyOpt.ClosedDiabatic;
              negativeZ.matEntOpt == MaterialEntropyOpt.OpenDiabatic or
              negativeZ.matEntOpt == MaterialEntropyOpt.ClosedDiabatic,
              positiveZ.matEntOpt == MaterialEntropyOpt.OpenDiabatic or
              positiveZ.matEntOpt == MaterialEntropyOpt.ClosedDiabatic][axis,
              side] then
            beta_S*(Sdot_face[axis, side] + (2*side - 3)*s*(if inclLin[axis]
               then I[linAxes[axis]] else 0)) = 2*Lstar_trans[axis]*(T_face[
              axis, side]/T - 1)*(if upstream[axis] and inclLin[axis] then (exp(
              (3 - 2*side)*I[linAxes[axis]]*beta_S/(2*Lstar_trans[axis])) + 1)
               else 2) "Diffusion";
            // Note:  Advection is included in material transport.
          else
            T_face[axis, side] = 0;
            Sdot_face[axis, side] = 0;
          end if;
        end for;
      end for;

      // Material dynamics
      if setPartNum then
        // Apply the IC for all time (material not conserved).
        if partNumInitMeth == InitMethScalar.ParticleNumber then
          N = N_IC;
        elseif partNumInitMeth == InitMethScalar.ParticleNumberRate then
          der(N)/U.s = derN_IC;
        elseif partNumInitMeth == InitMethScalar.Volume then
          inert.V = V_IC;
        elseif partNumInitMeth == InitMethScalar.VolumeRate then
          der(inert.V)/U.s = derV_IC;
        elseif partNumInitMeth == InitMethScalar.VolumeSpecific then
          v = v_IC;
        elseif partNumInitMeth == InitMethScalar.VolumeSpecificRate then
          der(v)/U.s = derv_IC;
        elseif partNumInitMeth == InitMethScalar.Pressure then
          inert.p = p_IC;
        elseif partNumInitMeth == InitMethScalar.PressureRate then
          der(inert.p)/U.s = derp_IC;
        elseif partNumInitMeth == InitMethScalar.Temperature then
          T = T_IC;
        elseif partNumInitMeth == InitMethScalar.TemperatureRate then
          der(T)/U.s = derT_IC;
        elseif partNumInitMeth == InitMethScalar.SpecificEntropy then
          s = s_IC;
        elseif partNumInitMeth == InitMethScalar.SpecificEntropyRate then
          der(s)/U.s = ders_IC;
        elseif partNumInitMeth == InitMethScalar.SpecificEnthalpy then
          h = h_IC;
        elseif partNumInitMeth == InitMethScalar.SpecificEnthalpyRate then
          der(h)/U.s = derh_IC;
        elseif partNumInitMeth == InitMethScalar.PotentialElectrochemical then
          mu = mu_IC;
        elseif partNumInitMeth == InitMethScalar.PotentialElectrochemicalRate
             then
          der(mu) = dermu_IC;
        else
          //if partNumInitMeth == InitMethScalar.ReactionRate then
          chemical.Ndot = Ndot_IC;
          // Note:  partNumInitMeth == InitMethScalar.None can't occur due to
          // an assertion.
        end if;
      else
        der(N)/U.s = chemical.Ndot + sum(Ndot_face) "Material conservation";
      end if;

      // Dynamics of linear momentum
      for axis in 1:n_lin loop
        if cartAxes[axis] == 1 and setXVel then
          // Apply the IC for all time (x-axis component not conserved).
          if xInitMeth == InitMethLinear.Velocity then
            phi[axis] = phi_IC[1];
          elseif xInitMeth == InitMethLinear.Acceleration then
            der(phi[axis])/U.s = derphi_IC[1];
          elseif xInitMeth == InitMethLinear.Current then
            I[axis] = I_IC[1];
          else
            //if xInitMeth == InitMethLinear.CurrentRate then
            der(I[axis])/U.s = derI_IC[1];
            // Note:  xInitMeth == InitMethLinear.None can't occur due to an
            // assertion.
          end if;
        elseif cartAxes[axis] == 2 and setYVel then
          // Apply the IC for all time (y-axis component not conserved).
          if yInitMeth == InitMethLinear.Velocity then
            phi[axis] = phi_IC[2];
          elseif yInitMeth == InitMethLinear.Acceleration then
            der(phi[axis])/U.s = derphi_IC[2];
          elseif yInitMeth == InitMethLinear.Current then
            I[axis] = I_IC[2];
          else
            //if yInitMeth == InitMethLinear.CurrentRate then
            der(I[axis])/U.s = derI_IC[2];
            // Note:  yInitMeth == InitMethLinear.None can't occur due to an
            // assertion.
          end if;
        elseif cartAxes[axis] == 3 and setZVel then
          // Apply the IC for all time (z-axis component not conserved).
          if zInitMeth == InitMethLinear.Velocity then
            phi[axis] = phi_IC[3];
          elseif zInitMeth == InitMethLinear.Acceleration then
            der(phi[axis])/U.s = derphi_IC[3];
          elseif zInitMeth == InitMethLinear.Current then
            I[axis] = I_IC[3];
          else
            //if zInitMeth == InitMethLinear.CurrentRate then
            der(I[axis])/U.s = derI_IC[3];
            // Note:  zInitMeth == InitMethLinear.None can't occur due to an
            // assertion.
          end if;
        else
          der(Data.m*N*phi[axis])/U.s = chemical.mPhidot[axis] + inert.mPhidot[
            axis] + sum(if {{negativeX.matEntOpt == MaterialEntropyOpt.OpenDiabatic,
            positiveX.matEntOpt == MaterialEntropyOpt.OpenDiabatic},{negativeY.matEntOpt
             == MaterialEntropyOpt.OpenDiabatic,positiveY.matEntOpt ==
            MaterialEntropyOpt.OpenDiabatic},{negativeZ.matEntOpt ==
            MaterialEntropyOpt.OpenDiabatic,positiveZ.matEntOpt ==
            MaterialEntropyOpt.OpenDiabatic}}[cartAxes[axis], side] then (3 - 2
            *side)*(mu_face[cartAxes[axis], side] - mu + Data.Deltah0_f)*N/L[
            cartAxes[axis]] else 0 for side in 1:2) + sum(Sigma(mPhidot_face[
            cartWrap(cartAxes[axis] - orientation), :, orientation]) for
            orientation in 1:2) "Conservation of linear momentum";
        end if;
      end for;

      // Thermal dynamics
      if setTemp then
        // Apply the IC for all time (energy not conserved).
        if tempInitMeth == InitMethScalar.ParticleNumber then
          N = N_IC;
        elseif tempInitMeth == InitMethScalar.ParticleNumberRate then
          der(N)/U.s = derN_IC;
        elseif tempInitMeth == InitMethScalar.Volume then
          inert.V = V_IC;
        elseif tempInitMeth == InitMethScalar.VolumeRate then
          der(inert.V)/U.s = derV_IC;
        elseif tempInitMeth == InitMethScalar.VolumeSpecific then
          v = v_IC;
        elseif tempInitMeth == InitMethScalar.VolumeSpecificRate then
          der(v)/U.s = derv_IC;
        elseif tempInitMeth == InitMethScalar.Pressure then
          inert.p = p_IC;
        elseif tempInitMeth == InitMethScalar.PressureRate then
          der(inert.p)/U.s = derp_IC;
        elseif tempInitMeth == InitMethScalar.Temperature then
          T = T_IC;
        elseif tempInitMeth == InitMethScalar.TemperatureRate then
          der(T)/U.s = derT_IC;
        elseif tempInitMeth == InitMethScalar.SpecificEntropy then
          s = s_IC;
        elseif tempInitMeth == InitMethScalar.SpecificEntropyRate then
          der(s)/U.s = ders_IC;
        elseif tempInitMeth == InitMethScalar.SpecificEnthalpy then
          h = h_IC;
        elseif tempInitMeth == InitMethScalar.SpecificEnthalpyRate then
          der(h)/U.s = derh_IC;
        elseif tempInitMeth == InitMethScalar.PotentialElectrochemical then
          mu = mu_IC;
        elseif tempInitMeth == InitMethScalar.PotentialElectrochemicalRate then
          der(mu) = dermu_IC;
        else
          //if tempInitMeth == InitMethScalar.ReactionRate then
          chemical.Ndot = Ndot_IC;
          // Note:  tempInitMeth == InitMethScalar.None can't occur due to an
          // assertion.
        end if;
      else
        (N*(Data.c0_T(T)*der(T) + der(ke)) + inert.V*der(Data.p_vT(inert.V/N, T)))
          /U.s = (chemical.mu - h - ke)*chemical.Ndot + semiLinear(
              chemical.Ndot,
              chemical.mphi*chemical.mphi/Data.m,
              Data.m*phi*phi) + chemical.TSdot + inert.phi*inert.mPhidot +
          inert.T*inert.Sdot + sum((mu_face - fill(
              h + ke - Data.Deltah0_f,
              3,
              2)) .* Ndot_face) + sum(phi_face .* mPhidot_face) + sum(T_face
           .* Sdot_face) "Conservation of energy";
        // Note:  Although it is mathematically equivalent,
        // der(Data.p_vT(inert.V/N, T)) is used instead of der(Data.p_vT(v, T)),
        // der(inert.p), or der(inert.p) so that the term can be expanded to
        // avoid dynamic state selection.
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

    <p align=center><img src=\"modelica://FCSys/resources/images/Subregions/BaseClasses/PartialSpecies/exchange.png\">
    <br><b>Figure 1:</b>  Exchange of a quantity (particle number, linear momentum, or entropy) among species (A, B, and C) within a subregion.</p>

    <p>Figure 2 shows the manner in which <a href=\"modelica://FCSys.Subregions.Species\">Species</a>
    instances
    of the same type are connected between neighboring <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> instances.
    Transport is similar to exchange
    except that advection and diffusion are directly coupled; both effects
    are included in the same connections.
    Upstream discretization is applied if it is enabled (via the <code>upstreamX</code>, etc. parameters).

    <p align=center><img src=\"modelica://FCSys/resources/images/Subregions/BaseClasses/PartialSpecies/transport.png\">
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
          <img src=\"modelica://FCSys/resources/images/Subregions/BaseClasses/PartialSpecies/share_pressure.png\">
          <br><b>a:</b>  Pressures of species (A, B, and C) are additive within a phase.
        </td>
        <td align=center>
          <img src=\"modelica://FCSys/resources/images/Subregions/BaseClasses/PartialSpecies/share_volume.png\">
          <br><b>b:</b>  Volumes of phases (I, II, and III) are additive within a subregion.
        </td>
      </tr>
      <tr align=center>
        <td colspan=2 align=center><b>Figure 3:</b> Methods of sharing pressure and volume.</td>
      </tr>
    </table>

    <p>The following variables reflect the actual properties of the
    species: <code>chemical.mphi</code> (specific mass times velocity), <code>chemical.Ts</code> (specific entropy times temperature),
     <code>inert.V</code> (volume), and <code>inert.p</code> (pressure).  However, due to exchange losses
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
    (<code>chemical</code> or <code>inert</code>) or faces (<code>negativeX</code>, <code>positiveX</code>, etc.) and both have zero resistivities for a
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
    Set <code>setPartNum</code> to <code>true</code> and <code>partNumInitMeth</code>
    to <code>InitMethScalar.Pressure</code> for the species.  In general, only one incompressible
    species can be included if there are no incompressible species.</li>
    <li>If <code>termDepleted</code> is <code>true</code>, then the simulation will be terminated when
    pressure is less than or equal to <code>Data.p_min</code>
    or the particle number is less than or equal to zero. <code>Data</code> is an instance
    of the <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</li>
    <li>The <code>start</code> values of the initial conditions for pressure and temperature
    (<i>p</i><sub>IC</sub> and <i>T</i><sub>IC</sub>) are the global default pressure and
    temperature (via the <code>outer</code> instance of the <a href=\"modelica://FCSys.BCs.Defaults\">Defaults</a> model).
    The <code>start</code> values of the initial conditions for
    other intensive properties (<i>v</i><sub>IC</sub>, <i>s</i><sub>IC</sub>, <i>h</i><sub>IC</sub>, and
    &mu;<sub>IC</sub>) are related to the initial pressure and temperature
    by the characteristics of the species.  The <code>start</code> value of the
    initial condition for the extensive volume (<i>V</i><sub>IC</sub>) is the volume of the
    subregion, and the <code>start</code> value for particle number (<i>N</i><sub>IC</sub>)
    is related to it via the characteristics (in <code>Data</code>) and the initial pressure and temperature.
    In order to apply other values for any of these initial conditions,
    it may be necessary to do so before translating the model.</li>
    <li>If the species has charged (i.e., is ionic) and permittivity (<code>epsilon</code>) is
    zero, then it should be set as <code>final</code> to eliminate
    the associated state.  Otherwise, errors may occur.</li>
    </p>

    <p>In order to reduce numerical error during simulation, enthalpy of formation
    (<code>Data.Deltah0_f</code>) is excluded
    from the face connectors (e.g., <code>negativeX.material.mu</code>).  There is
    no mathematical effect since the linear momentum and energy balances are adjusted accordingly.
    However, enthalpy of formation is included in the chemical connector
    (<code>chemical.mu</code>); it is necessary for the proper chemical equilibrium.</p>

    <p>In evaluating the dynamics of a phase, it is usually assumed
    that all of the species exist at the same temperature.
    The time constants that govern the temperatures/heat capacities of the species and entropy flow rates among them
    are usually
    much shorter than the time span of interest.
    This assumption can be applied in the model by setting the thermal exchange resistivities
    (&alpha;<sub><i>S</i></sub>) of the species as <code>final</code> parameters equal to zero.
    Then, the translator can perform index reduction and retain only one
    state associated with temperature.  However, this will likely lead to nonlinear systems of equations
    and may reduce the performance of the simulation. Likewise, if the reaction rates are very fast with respect to the
    observed time span,
    then the reaction resistivities (&alpha;<sub><i>N</i></sub>) of the species may be redeclared as
    <code>final</code> parameters and set to zero.  A similar situation applies to
    momentum exchange (&alpha;<sub>&Phi;</sub>), material transport (&beta;<sub><i>N</i></sub>),
    compressive momentum transport
    (&beta;<sub>&#8214;</sub>), transverse momentum transport (&beta;<sub>&Phi;</sub>),
    and thermal transport (&beta;<sub><i>S</i></sub>).</p>

    <p>In the variables that relate to transport,
    the first index is the axis and the second index is the side.  The sides
    are indexed from negative (1) to positive (2).  Linear momentum is additionally indexed by
    the orientation of the momentum with respect to the face.
    The first index corresponds to momentum oriented along the axis following the face-normal axis
    in Cartesian space (x, y, z).
    The second index is of linear momentum oriented along the axis twice after the face-normal axis in
    Cartesian space
    (or, equivalently, the axis preceding it).</p>
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
                  fillPattern=FillPattern.Solid),Text(
                  extent={{-100,-20},{100,20}},
                  textString="%name",
                  lineColor={0,0,0})}));
    end Species;
  end Species;

  model PhaseBoundary
    "Phase boundary (adapter between Amagat and Dalton mixtures)"
    extends FCSys.BaseClasses.Icons.Names.Top6;
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

  model Reaction "Model for a chemical reaction"
    //extends FCSys.BaseClasses.Icons.Names.Top2;

    parameter Integer n_lin=1
      "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>"
      annotation (Evaluate=true,HideResult=true);
    parameter Integer n_spec(min=2) = 0 "Number of species"
      annotation (Dialog(connectorSizing=true));
    // Note:  The minimum is 2 for a meaningful reaction, but the default must
    // be 0 to use connectorSizing.
    final Integer nu[n_spec]=Chemistry.stoich(chemical.formula)
      "Stoichiometric coefficients";
    //  Q.MassSpecific m[n_spec]=chemical.m "Specific masses";
    // Note:  As of Modelica 3.2 and Dymola 7.4, this can't be a parameter or
    // constant even though it isn't time-varying.  The strings that
    // represent the chemical formulas cannnot be passed through the
    // connectors with parameter or constant prefixes.  However, the
    // translator should recognize that these equations are static.
    Q.VelocityMassSpecific mphi[n_lin] "Specific mass times velocity";
    Q.Potential Ts "Temperature times specific entropy";
    Q.Current Xidot(nominal=1*U.A) "Reaction rate";
    Connectors.ChemicalInput chemical[n_spec](each final n_lin=n_lin)
      "Chemical connector"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

  equation
    // Chemical equilibrium
    nu*chemical.mu = 0;

    // Ideal mixing/upstream discretization
    for i in 1:n_spec loop
      chemical[i].mPhidot = semiLinear(
          chemical[i].Ndot,
          chemical[i].mphi,
          mphi);
    end for;
    chemical.TSdot = semiLinear(
        chemical.Ndot,
        chemical.Ts,
        Ts);

    // Conservation (without storage)
    nu*Xidot = chemical.Ndot "Material";
    sum(chemical[i].mPhidot for i in 1:n_spec) = zeros(n_lin) "Linear momentum";
    sum(chemical.TSdot) = 0 "Energy";

    // Note:  This model is marked as structurally incomplete.  It must have
    // no species by default (for automatic connector sizing), but at least
    // one species is mathematically required (and two for a meaningful
    // reaction).
    annotation (
      structurallyIncomplete=true,
      Documentation(info="<html>
    <p>The size of the chemical connector is automatically increased each time a connection is made.
    At least two species must be connected.
    The stoichiometry is determined automatically from the chemical formulas
    of the connected species.  No intermediate species are considered. Each reaction must be
    completely and uniquely defined by the connected species.  Otherwise an error message is given.
    If you suspect a bug in the library, please report it using the
    <a href=\"modelica://FCSys.UsersGuide.Contact\">contact information</a>.</p>

    <p>For material, this model is essentially the opposite of a standard single-species connection.
    The stoichiometric sum of the efforts (&Sigma; &nu;<sub><i>i</i></sub> &mu;<sub><i>i</i></sub>)
    is zero, which is analogous to Kirchhoff's Current Law.  The flow rates divided by the
    stoichiometric coefficients (<i>N&#775;</i><sub><i>i</i></sub> /&nu;<sub><i>i</i></sub>)
    are equal&mdash;analogous to Kirchhoff's Voltage Law.</p>

    <p>Momentum and energy are advected using <code>stream</code> variables.  There is no diffusion;
    it is included in the inert connections among species
    (see the <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSpecies\">PartialSpecies</a> model).<p>

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
            extent={{-100,40},{100,-40}},
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
    parameter Q.Volume V=1*U.cm^3 "Volume" annotation (Dialog(group="Geometry"));
    parameter Integer n_lin(
      final min=0,
      final max=3) = 1
      "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>"
      annotation (
      Evaluate=true,
      HideResult=true,
      Dialog(group="Geometry"));

    FCSys.Connectors.InertAmagat inert(final n_lin=n_lin)
      "Connector for linear momentum and entropy, with shared volume"
      annotation (Placement(transformation(extent={{60,-80},{80,-60}}),
          iconTransformation(extent={{100,-120},{120,-100}})));

  equation
    // Specified volume
    V = inert.V;

    // Rate balances (without storage or generation)
    zeros(n_lin) = inert.mPhidot "Linear momentum";
    0 = inert.Sdot "Entropy";
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
      inner parameter Q.Length L[3](each final min=Modelica.Constants.small) =
        ones(3)*U.m "<html>Length (<b>L</b>)</html>"
        annotation (Dialog(group="Geometry"));
      final inner parameter Q.Area A[3]={L[cartWrap(ax + 1)]*L[cartWrap(ax + 2)]
          for ax in 1:3} "Cross-sectional area";
      final parameter Q.Volume V=product(L) "Volume";

      // Assumptions about components of linear momentum
      parameter Boolean inclXMom=true "X" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          tab="Assumptions",
          group="Axes with linear momentum included",
          compact=true));
      parameter Boolean inclYMom=false "Y" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          tab="Assumptions",
          group="Axes with linear momentum included",
          compact=true));
      parameter Boolean inclZMom=false "Z" annotation (
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

      FCSys.Connectors.FaceBus negativeX if inclXFaces
        "Negative face along the x axis" annotation (Placement(transformation(
              extent={{-50,-10},{-30,10}}), iconTransformation(extent={{-110,-10},
                {-90,10}})));
      FCSys.Connectors.FaceBus positiveX if inclXFaces
        "Positive face along the x axis" annotation (Placement(transformation(
              extent={{30,-10},{50,10}}), iconTransformation(extent={{90,-10},{
                110,10}})));
      FCSys.Connectors.FaceBus negativeY if inclYFaces
        "Negative face along the y axis" annotation (Placement(transformation(
              extent={{-10,-50},{10,-30}}), iconTransformation(extent={{-10,-110},
                {10,-90}})));
      FCSys.Connectors.FaceBus positiveY if inclYFaces
        "Positive face along the y axis" annotation (Placement(transformation(
              extent={{-10,30},{10,50}}), iconTransformation(extent={{-10,90},{
                10,110}})));
      FCSys.Connectors.FaceBus negativeZ if inclZFaces
        "Negative face along the z axis" annotation (Placement(transformation(
              extent={{10,10},{30,30}}), iconTransformation(extent={{40,40},{60,
                60}})));
      FCSys.Connectors.FaceBus positiveZ if inclZFaces
        "Positive face along the z axis" annotation (Placement(transformation(
              extent={{-30,-30},{-10,-10}}), iconTransformation(extent={{-60,-60},
                {-40,-40}})));

      FCSys.Subregions.Volume volume(final n_lin=n_lin,final V=V)
        "Model to establish space for species"
        annotation (Placement(transformation(extent={{-16,-16},{16,16}})));
    protected
      final parameter Integer n_lin=countTrue({inclXMom,inclYMom,inclZMom})
        "Number of components of linear momentum" annotation (Evaluate=true);

      annotation (
        Documentation(info="<html><p>Notes:
<ul><li>This model must be be extended so that models can be added for
  relevant species, phases, and reactions.</li>
  <li>Material will be transported between two subregions only if both of the connected faces are marked
  as open (<code>matEntOpt==MaterialEntropyOpt.OpenDiabatic</code>)
  within the instances of the matched <a href=\"modelica://FCSys.Subregions.Species\">Species</a> models.
  If either or both of the faces are closed (<code>matEntOpt==MaterialEntropyOpt.ClosedAdiabatic</code> or
  <code>matEntOpt==MaterialEntropyOpt.ClosedDiabatic</code>), then the interface will be closed.
  note applies to the viscous/inviscous and diabatic/adiabatic properties.</li>
  <li>The x-axis component of linear momentum is included by default.  At least one component must be included.</li></ul></p></html>"),

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

    type InitMethScalar = enumeration(
        None "Do not explicitly initialize.",
        ParticleNumber "Initialize the particle number.",
        ParticleNumberRate "Initialize the rate of ditto.",
        Volume "Initialize the volume.",
        VolumeRate "Initialize the rate of ditto.",
        VolumeSpecific "Initialize the specific volume.",
        VolumeSpecificRate "Initialize the rate of ditto.",
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
    type InitMethLinear = enumeration(
        None "Do not explicitly initialize.",
        Velocity "Initialize the velocity.",
        Acceleration "Initialize the acceleration.",
        Current "Initialize the current.",
        CurrentRate "Initialize the rate of ditto.")
      "Methods of initializing linear momentum";
  end BaseClasses;

  annotation (Documentation(info="<html>
<p>
<b>Licensed by Kevin Davies under the Modelica License 2</b><br>
Copyright 2007&ndash;2012, Kevin Davies.
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
