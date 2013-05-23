within FCSys;
package Subregions
  "Control volumes with multi-species transport, exchange, and storage"
  extends Modelica.Icons.Package;
  package Examples "Examples"
    extends Modelica.Icons.ExamplesPackage;

    model Subregion
      "<html>Evaluate a single subregion, with H<sub>2</sub> by default</html>"

      extends Modelica.Icons.Example;

      parameter Boolean 'inclC+'=false
        "<html>Carbon plus (C<sup>+</sup>)</html>" annotation (choices(
            __Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean 'inclC19HF37O5S-'=false
        "<html>Nafion sulfonate (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean 'incle-'=false "<html>Electrons (e<sup>-</sup>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean 'inclH+'=false "<html>Protons (H<sup>+</sup>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclH2=true "<html>Hydrogen (H<sub>2</sub>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclH2O=false
        "<html>Water vapor (H<sub>2</sub>O)</html>" annotation (choices(
            __Dymola_checkBox=true), Dialog(group="Species",
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
          final 'inclC+'='inclC+',
          final 'incle-'='incle-',
          'C+'(V_IC=subregion.V/4)),
        ionomer(
          final 'inclC19HF37O5S-'='inclC19HF37O5S-',
          final 'inclH+'='inclH+',
          'C19HF37O5S-'(V_IC=subregion.V/4)),
        liquid(H2O(V_IC=subregion.V/4)),
        inclTransY=false,
        inclTransZ=false,
        inclFacesX=false,
        inclFacesY=false,
        inclFacesZ=false)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      inner FCSys.Conditions.Environment environment(analysis=true)
        annotation (Placement(transformation(extent={{30,20},{50,40}})));
      annotation (
        Placement(transformation(extent={{70,70},{90,90}})),
        experiment(StopTime=10),
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.Subregion.mos"));

    end Subregion;

    model SubregionHOR
      "<html>Test a subregion with the hydrogen oxidation reaction and the essential species for it (C+, C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S, e<sup>-</sup>, H<sub>2</sub>, and H<sup>+</sup>)</html>"

      extends Examples.Subregion(
        'inclC+'=true,
        'inclC19HF37O5S-'=true,
        'incle-'=true,
        'inclH+'=true,
        inclH2=true,
        subregion(
          inclFacesX=true,
          graphite('e-'(redeclare package Data = FCSys.Characteristics.'e-'.Gas)),

          ionomer('H+'(redeclare package Data = FCSys.Characteristics.'H+'.Gas))));

      Conditions.FaceBus.SubregionFlows reactionBC(
        graphite(
          final 'inclC+'='inclC+',
          final 'incle-'='incle-',
          'e-'(redeclare Conditions.Face.Material.Current material(redeclare
                Modelica.Blocks.Sources.Ramp source(duration=1000, height=0.001
                    *U.A)))),
        ionomer(
          final 'inclC19HF37O5S-'='inclC19HF37O5S-',
          final 'inclH+'='inclH+',
          'H+'(redeclare Conditions.Face.Material.Pressure material(source(k=U.atm)))),

        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={20,0})));

      extends Modelica.Icons.UnderConstruction;
      // **fails sim

    equation
      connect(subregion.xPositive, reactionBC.face) annotation (Line(
          points={{10,6.10623e-16},{14,6.10623e-16},{14,1.23436e-15},{16,
              1.23436e-15}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (experiment(StopTime=1000, Tolerance=1e-06), Commands(file(
              ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionHOR.mos"));
    end SubregionHOR;

    model SubregionORR
      "<html>Test a subregion with the oxygen reduction reaction and the essential species for it (C+, C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S, e<sup>-</sup>, H<sub>2</sub>O, H<sup>+</sup>, and O<sub>2</sub>)</html>"

      extends Examples.SubregionHOR(
        inclH2=false,
        inclH2O=true,
        inclO2=true);
      annotation (experiment(StopTime=1000, Tolerance=1e-06), Commands(file(
              ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionORR.mos"));

    end SubregionORR;

    model Subregions
      "<html>Test a one-dimensional array of subregions with an initial pressure gradient (H<sub>2</sub> included by default)</html>"
      extends Modelica.Icons.Example;

      parameter Integer n_x=0
        "Number of discrete subregions along the x axis, besides the 2 side subregions";
      parameter Boolean 'inclC+'=false
        "<html>Carbon plus (C<sup>+</sup>)</html>" annotation (choices(
            __Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean 'inclC19HF37O5S-'=false
        "<html>Nafion sulfonate (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean 'incle-'=false "<html>Electrons (e<sup>-</sup>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclH2=true "<html>Hydrogen (H<sub>2</sub>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclH2O=false
        "<html>Water vapor (H<sub>2</sub>O)</html>" annotation (choices(
            __Dymola_checkBox=true), Dialog(group="Species",
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
      parameter Q.Pressure Deltap_IC=10*U.kPa
        "<html>Initial pressure difference (&Delta;<i>p</i><sub>IC</sub>)</html>";

      FCSys.Subregions.Subregion subregion1(
        L={10,1,1}*U.mm,
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2O(p_IC=environment.p + Deltap_IC/2),
          N2(p_IC=environment.p + Deltap_IC/2),
          O2(p_IC=environment.p + Deltap_IC/2),
          H2(p_IC=environment.p + Deltap_IC/2)),
        graphite(
          final 'inclC+'='inclC+',
          final 'incle-'='incle-',
          'C+'(V_IC=subregion1.V/4)),
        ionomer(
          final 'inclC19HF37O5S-'='inclC19HF37O5S-',
          final 'inclH+'='inclH+',
          'C19HF37O5S-'(V_IC=subregion1.V/4)),
        liquid(H2O(V_IC=subregion1.V/4)),
        inclFacesY=false,
        inclFacesZ=false,
        inclTransY=false,
        inclTransZ=false)
        annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));

      FCSys.Subregions.Subregion subregions[n_x](
        each L={10,1,1}*U.mm,
        gas(
          each final inclH2=inclH2,
          each final inclH2O=inclH2O,
          each final inclN2=inclN2,
          each final inclO2=inclO2,
          H2(p_IC={environment.p + Deltap_IC/2 - i*Deltap_IC/(n_x + 1) for i
                 in 1:n_x}),
          H2O(p_IC=subregions.gas.H2.p_IC),
          N2(p_IC=subregions.gas.H2.p_IC),
          O2(p_IC=subregions.gas.H2.p_IC)),
        graphite(
          each final 'inclC+'='inclC+',
          each final 'incle-'='incle-',
          'C+'(each V_IC=subregions[1].V/4)),
        ionomer(
          each final 'inclC19HF37O5S-'='inclC19HF37O5S-',
          each final 'inclH+'='inclH+',
          'C19HF37O5S-'(each V_IC=subregions[1].V/4)),
        each inclTransY=false,
        each inclTransZ=false,
        each inclFacesY=false,
        each inclFacesZ=false) if n_x > 0
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      FCSys.Subregions.Subregion subregion2(
        L={10,1,1}*U.mm,
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2(p_IC=environment.p - Deltap_IC/2),
          H2O(p_IC=environment.p - Deltap_IC/2),
          N2(p_IC=environment.p - Deltap_IC/2),
          O2(p_IC=environment.p - Deltap_IC/2)),
        graphite(
          final 'inclC+'='inclC+',
          final 'incle-'='incle-',
          'C+'(V_IC=subregion2.V/4)),
        ionomer(
          final 'inclC19HF37O5S-'='inclC19HF37O5S-',
          final 'inclH+'='inclH+',
          'C19HF37O5S-'(V_IC=subregion2.V/4)),
        inclFacesY=false,
        inclFacesZ=false,
        inclTransY=false,
        inclTransZ=false)
        annotation (Placement(transformation(extent={{20,-10},{40,10}})));

      inner FCSys.Conditions.Environment environment(analysis=true)
        annotation (Placement(transformation(extent={{60,20},{80,40}})));

    equation
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
      annotation (
        Placement(transformation(extent={{70,70},{90,90}})),
        experiment(
          StopTime=250,
          Tolerance=1e-06,
          Algorithm="Dassl"),
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.Subregions.mos"));
    end Subregions;

    model SubregionsCAndH2
      "<html>Test a one-dimensional array of subregions with C and H<sub>2</sub></html>"
      extends Examples.Subregions(
        'inclC+'=true,
        subregion1(graphite('C+'(
              setVelX=false,
              setVelY=false,
              setVelZ=false))),
        environment(analysis=false));
      annotation (experiment(StopTime=3, NumberOfIntervals=5000), Commands(file(
              ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionsCAndH2.mos"));

    end SubregionsCAndH2;

    model SubregionsCAndH2AndH2O
      "<html>Test a one-dimensional array of subregions with C, H<sub>2</sub>, and H<sub>2</sub>O</html>"
      extends Examples.Subregions(
        'inclC+'=true,
        inclH2=true,
        inclH2O=true);
      annotation (experiment(StopTime=8, Tolerance=1e-06), Commands(file(
              ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionsCAndH2AndH2O.mos"));

    end SubregionsCAndH2AndH2O;

    model SubregionsCAndN2
      "<html>Test a one-dimensional array of subregions with C and N<sub>2</sub></html>"
      extends Examples.Subregions(
        inclH2=false,
        'inclC+'=true,
        inclN2=true);
      annotation (experiment(StopTime=40, Tolerance=1e-06), Commands(file(
              ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionsCAndN2.mos"));

    end SubregionsCAndN2;

    model SubregionsCAndO2
      "<html>Test a one-dimensional array of subregions with C and O<sub>2</sub></html>"
      extends Examples.Subregions(
        inclH2=false,
        'inclC+'=true,
        inclO2=true);
      annotation (experiment(
          StopTime=25,
          NumberOfIntervals=5000,
          Tolerance=1e-06,
          Algorithm="Dassl"), Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionsCAndO2.mos"));

    end SubregionsCAndO2;

    model SubregionsC19HF37O5SminusAndH2OAndHplus
      "<html>Test a one-dimensional array of subregions with C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S, H<sub>2</sub>O, and H<sup>+</sup></html>"
      extends Examples.Subregions(
        inclH2=false,
        inclH2O=true,
        'inclH+'=true,
        'inclC19HF37O5S-'=true);
      extends Modelica.Icons.UnderConstruction;
      annotation (experiment(StopTime=150, Tolerance=1e-06), Commands(file(
              ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionsC19HF37O5SminusAndH2OAndHplus.mos"));

    end SubregionsC19HF37O5SminusAndH2OAndHplus;

    model SubregionsCAndeminus
      "<html>Test a one-dimensional array of subregions with e<sup>-</sup></html>"
      extends Examples.Subregions(
        'inclC+'=true,
        inclH2=false,
        'incle-'=true,
        environment(analysis=false));
      extends Modelica.Icons.UnderConstruction;
      annotation (experiment(Tolerance=1e-06), Commands(file(ensureSimulated=
                true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionsCAndeminus.mos"));

    end SubregionsCAndeminus;

    model ThermalConduction "Test thermal conduction (through solid)"
      extends Examples.Subregions(
        n_x=8,
        'inclC+'=true,
        inclH2=false,
        subregion1(graphite('C+'(
              initMaterial=InitScalar.Pressure,
              V_IC=0.99*subregion1.V,
              T_IC=1.1*environment.T,
              T(displayUnit="degC")))),
        subregions(graphite('C+'(
              each initMaterial=InitScalar.Pressure,
              each V_IC=0.99*subregions[1].V,
              each T(displayUnit="degC")))),
        subregion2(graphite('C+'(
              initMaterial=InitScalar.Pressure,
              V_IC=0.99*subregion2.V,
              T(displayUnit="degC")))),
        redeclare FCSys.Conditions.FaceBus.SubregionClosedAdiabatic condition1,

        redeclare FCSys.Conditions.FaceBus.SubregionClosedAdiabatic condition2);

      annotation (Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.ThermalConduction.mos"),
          experiment(StopTime=298.15, Algorithm="Dassl"));

    end ThermalConduction;

    model ThermalConductionConvection
      "Test combined thermal conduction and convection"
      extends Examples.Subregions(
        n_x=8,
        'inclC+'=true,
        inclH2=false,
        inclN2=true,
        subregion1(gas(N2(
              T_IC=1.1*environment.T,
              p_IC=environment.p,
              phi(displayUnit="mm/s"))), graphite('C+'(
              V_IC=0.5*subregion1.V,
              T_IC=1.1*environment.T,
              T(displayUnit="degC")))),
        subregions(gas(N2(each p_IC=environment.p, phi(each displayUnit="mm/s"))),
            graphite('C+'(each V_IC=0.5*subregions[1].V, each T(displayUnit=
                    "degC")))),
        subregion2(gas(N2(p_IC=environment.p, phi(displayUnit="mm/s"))),
            graphite('C+'(V_IC=0.5*subregion2.V, T(displayUnit="degC")))),
        redeclare FCSys.Conditions.FaceBus.SubregionClosedAdiabatic condition1,

        redeclare FCSys.Conditions.FaceBus.SubregionClosedAdiabatic condition2);

      annotation (Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.ThermalConductionConvection.mos"),
          experiment(StopTime=200, Algorithm="Dassl"));

    end ThermalConductionConvection;

    model ReactionRamp
      "Test chemical reaction with reaction rate ramped over time"

      extends Modelica.Icons.Example;
      extends Modelica.Icons.UnderConstruction;

      Reaction reaction(n_spec=3)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Conditions.Chemical.Species species1(redeclare
          Characteristics.'e-'.Graphite Data, material=FCSys.Conditions.Chemical.BaseClasses.ConditionTypeMaterial.PotentialPerTemperature)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-30,-24})));

      Conditions.Chemical.Species species2(redeclare Characteristics.'H+'.Gas
          Data, material=FCSys.Conditions.Chemical.BaseClasses.ConditionTypeMaterial.PotentialPerTemperature)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={0,-24})));

      Conditions.Chemical.Species species3(
        redeclare FCSys.Characteristics.H2.Gas Data,
        material=FCSys.Conditions.Chemical.BaseClasses.ConditionTypeMaterial.Current,

        redeclare Modelica.Blocks.Sources.Ramp materialSpec(height=100*U.A,
            duration=3600e2)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={30,-24})));

      inner FCSys.Conditions.Environment environment(analysis=true)
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
      annotation (experiment(StopTime=36000), Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.ReactionRamp.mos"));
    end ReactionRamp;

    model Reaction "Test an electrochemical reaction"

      extends Modelica.Icons.Example;
      extends Modelica.Icons.UnderConstruction;

      parameter Integer n_trans(
        final min=1,
        final max=3) = 1
        "<html>Number of components of translational momentum (<i>n</i><sub>trans</sub>)</html>";

      Reaction reaction(n_spec=3)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Conditions.Chemical.Species 'e-'(
        material=FCSys.Conditions.Chemical.BaseClasses.ConditionTypeMaterial.PotentialPerTemperature,

        redeclare FCSys.Characteristics.'e-'.Graphite Data,
        final n_trans=n_trans) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={-30,-24})));

      Conditions.Chemical.Species 'H+'(
        material=FCSys.Conditions.Chemical.BaseClasses.ConditionTypeMaterial.PotentialPerTemperature,

        redeclare FCSys.Characteristics.'H+'.Ionomer Data,
        final n_trans=n_trans) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={0,-24})));

      Conditions.Chemical.Species H2(
        material=FCSys.Conditions.Chemical.BaseClasses.ConditionTypeMaterial.Current,

        redeclare FCSys.Characteristics.H2.Gas Data,
        final n_trans=n_trans,
        redeclare Modelica.Blocks.Sources.Ramp materialSpec(height=100*U.A,
            duration=100)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=180,
            origin={30,-24})));

      inner FCSys.Conditions.Environment environment(analysis=true)
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
      annotation (experiment(StopTime=100), Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.Reaction.mos"));
    end Reaction;

    model TestSpecies "Test the Species model"
      extends Modelica.Icons.Example;

      import FCSys.BaseClasses.Utilities.cartWrap;
      import FCSys.BaseClasses.Utilities.countTrue;
      import FCSys.BaseClasses.Utilities.enumerate;
      import FCSys.BaseClasses.Utilities.index;
      // extends FCSys.BaseClasses.Icons.Names.Top3;

      // Geometric parameters
      inner parameter Q.Length L[Axis](each min=Modelica.Constants.small) = {U.cm,
        U.cm,U.cm} "<html>Length (<b>L</b>)</html>"
        annotation (Dialog(group="Geometry"));
      final inner parameter Q.Volume V=product(L) "Volume";
      parameter Q.NumberAbsolute k[Axis](
        each min=Modelica.Constants.small,
        each final nominal=1) = {1,1,1}
        "<html>Area fill factor (<b>k</b>)</html>"
        annotation (Dialog(group="Geometry"));

      // Assumptions
      // -----------
      // Included components of translational momentum
      parameter Boolean inclTransX=true "X" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));
      parameter Boolean inclTransY=true "Y" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));
      parameter Boolean inclTransZ=true "Z" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));

      replaceable FCSys.Subregions.Species.Species species(redeclare package
          Data = FCSys.Characteristics.H2.Gas, n_react=1) constrainedby
        FCSys.Subregions.Species.Species
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

    protected
      final inner parameter Q.Length Lprime[:]=k .* A ./ L
        "Effective cross-sectional area per length";
      final inner parameter Q.Area A[Axis]={L[cartWrap(axis + 1)]*L[cartWrap(
          axis + 2)] for axis in Axis} "Cross-sectional areas";
      final inner parameter Boolean inclTrans[Axis]={inclTransX,inclTransY,
          inclTransZ}
        "true, if each component of translational momentum is included";
      final inner parameter Integer n_trans=countTrue(inclTrans)
        "Number of components of translational momentum";
      final inner parameter Integer cartAxes[:]=index(inclTrans)
        "Cartesian-axis indices of the components of translational momentum";
      final inner parameter Integer transAxes[Axis]=enumerate(inclTrans)
        "Translational-momentum-component indices of the Cartesian axes";

      Volume volume "Model to establish a fixed total volume"
        annotation (Placement(transformation(extent={{-16,-92},{16,-60}})));
      PhaseBoundary phaseBoundary "Phase boundary" annotation (Placement(
            transformation(
            extent={{-18,-18},{18,18}},
            rotation=0,
            origin={0,-36})));
    public
      inner Conditions.Environment environment(analysis=true)
        annotation (Placement(transformation(extent={{70,72},{90,92}})));
    equation
      //species.chemical[1].axis = 1;
      connect(phaseBoundary.inertAmagat, volume.inert) annotation (Line(
          points={{6.2,-48.2},{6.2,-56.1},{11,-56.1},{11,-87}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(species.inertDalton, phaseBoundary.inertDalton) annotation (Line(
          points={{3.578,-7.155},{3.578,-16.155},{3.578,-16.155},{3.578,-25.155},
              {3.578,-43.155},{3.578,-43.155}},
          color={72,90,180},
          smooth=Smooth.None));
      annotation (Diagram(graphics));
    end TestSpecies;

    model Specieseminus "Test a species"
      extends Modelica.Icons.Example;
      extends Modelica.Icons.UnderConstruction;

      extends FCSys.Subregions.Examples.TestSpecies(redeclare
          Species.'e-'.Graphite.Fixed species);
      Conditions.Face.BaseClasses.PartialSpecies faceCondition(redeclare
          Conditions.Face.Material.Pressure material) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-24,0})));

    equation
      connect(faceCondition.face, species.xNegative) annotation (Line(
          points={{-20,3.65701e-16},{-14,3.65701e-16},{-14,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      annotation (
        Placement(transformation(extent={{70,70},{90,90}})),
        experiment(StopTime=10),
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionH2.mos"));
    end Specieseminus;

    model SubregionsCell "Test a one-dimensional array of subregions"
      extends Modelica.Icons.Example;
      extends Modelica.Icons.UnderConstruction;

      parameter Integer n_x=1
        "Number of discrete regions through the PEM along the x axis";

      Subregion subregion1(
        L={1,1,1}*U.cm,
        gas(
          inclH2=true,
          inclH2O=false,
          inclN2=false,
          inclO2=false,
          H2O(
            yNegative(inviscidX=true),
            yPositive(inviscidX=true),
            zNegative(inviscidX=true),
            zPositive(inviscidX=true)),
          H2(
            p_IC=1.05*environment.p,
            yNegative(inviscidX=true),
            yPositive(inviscidX=true),
            zNegative(inviscidX=true),
            zPositive(inviscidX=true))),
        graphite(
          'inclC+'=true,
          'incle-'=true,
          'C+'(V_IC=subregion1.V/4),
          'e-'(
            yNegative(inviscidX=true),
            yPositive(inviscidX=true),
            zNegative(inviscidX=true),
            zPositive(inviscidX=true))),
        ionomer(
          'inclC19HF37O5S-'=true,
          'inclH+'=true,
          'C19HF37O5S-'(V_IC=subregion1.V/4),
          'H+'(
            yNegative(inviscidX=true),
            yPositive(inviscidX=true),
            zNegative(inviscidX=true),
            zPositive(inviscidX=true))),
        inclTransY=false,
        inclTransZ=false,
        inclFacesY=false,
        inclFacesZ=false)
        annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));

      Subregion subregions[n_x](
        each L={1,1,1}*U.cm,
        gas(each inclH2O=false),
        ionomer(
          each 'inclC19HF37O5S-'=true,
          each 'inclH+'=true,
          'C19HF37O5S-'(each V_IC=subregions[1].V/4),
          'H+'(
            each yNegative(inviscidX=true),
            each yPositive(inviscidX=true),
            each zNegative(inviscidX=true),
            each zPositive(inviscidX=true))),
        each inclTransY=false,
        each inclTransZ=false,
        each inclFacesY=false,
        each inclFacesZ=false) if n_x > 0
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      Subregion subregion2(
        L={1,1,1}*U.cm,
        gas(
          inclH2=false,
          inclH2O=true,
          inclN2=false,
          inclO2=true,
          H2O(
            p_IC=0.95*environment.p,
            yNegative(inviscidX=true),
            yPositive(inviscidX=true),
            zNegative(inviscidX=true),
            zPositive(inviscidX=true)),
          N2(
            p_IC=0.95*environment.p,
            yNegative(inviscidX=true),
            yPositive(inviscidX=true),
            zNegative(inviscidX=true),
            zPositive(inviscidX=true)),
          O2(
            p_IC=0.95*environment.p,
            yNegative(inviscidX=true),
            yPositive(inviscidX=true),
            zNegative(inviscidX=true),
            zPositive(inviscidX=true))),
        graphite(
          'inclC+'=true,
          'incle-'=true,
          'C+'(V_IC=subregion2.V/4),
          'e-'(
            yNegative(inviscidX=true),
            yPositive(inviscidX=true),
            zNegative(inviscidX=true),
            zPositive(inviscidX=true))),
        ionomer(
          'inclC19HF37O5S-'=true,
          'inclH+'=true,
          'C19HF37O5S-'(V_IC=subregion2.V/4),
          'H+'(
            yNegative(inviscidX=true),
            yPositive(inviscidX=true),
            zNegative(inviscidX=true),
            zPositive(inviscidX=true))),
        inclTransY=false,
        inclTransZ=false,
        inclFacesY=false,
        inclFacesZ=false)
        annotation (Placement(transformation(extent={{20,-10},{40,10}})));

      inner FCSys.Conditions.Environment environment(analysis=true)
        annotation (Placement(transformation(extent={{60,20},{80,40}})));
      replaceable FCSys.Conditions.FaceBus.SubregionFlows condition1
        constrainedby FCSys.Conditions.FaceBus.Subregion(gas(
          inclH2=true,
          inclH2O=false,
          inclN2=false,
          inclO2=false), graphite(
          'inclC+'=true,
          'incle-'=true,
          'C+'(redeclare FCSys.Conditions.Face.Thermal.Temperature thermal(
                source(k=environment.T))),
          'e-'(redeclare FCSys.Conditions.Face.Material.Current normal(
                redeclare Modelica.Blocks.Sources.Ramp source(duration=1000,
                  height=-2*U.A))))) annotation (__Dymola_choicesFromPackage=
            true, Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-50,0})));

      replaceable FCSys.Conditions.FaceBus.SubregionFlows condition2
        constrainedby FCSys.Conditions.FaceBus.Subregion(gas(
          inclH2=false,
          inclH2O=true,
          inclN2=false,
          inclO2=true), graphite(
          'inclC+'=true,
          'incle-'=true,
          'C+'(redeclare FCSys.Conditions.Face.Thermal.Temperature thermal(
                source(k=environment.T))),
          'e-'(redeclare FCSys.Conditions.Face.Material.Current normal(
                redeclare Modelica.Blocks.Sources.Ramp source(duration=1000,
                  height=2*U.A))))) annotation (__Dymola_choicesFromPackage=
            true, Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={50,0})));

      replaceable FCSys.Conditions.FaceBus.SubregionFlows ground(graphite(
            'incle-'=true, 'e-'(redeclare Conditions.Face.Material.Pressure
              material, materialSpec(k=0)))) constrainedby
        Conditions.FaceBus.Subregion annotation (__Dymola_choicesFromPackage=
            true, Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-50,30})));

    equation
      connect(condition1.face, subregion1.xNegative) annotation (Line(
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
      connect(condition2.face, subregion2.xPositive) annotation (Line(
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
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionsH2.mos"));
    end SubregionsCell;

    model SubregionH2PipeTest
      extends Examples.Subregion(subregion(L={100,1,1}*U.mm), condition1(gas(H2
              (redeclare FCSys.Conditions.Face.Material.Density normal(
                  redeclare Modelica.Blocks.Sources.Ramp source(height=U.C/U.cc,
                    offset=298.15*U.K/U.atm))))));

    end SubregionH2PipeTest;

    model SubregionCplusAlone
      "<html>Test a subregion with H<sub>2</sub>, without external boundary conditions</html>"

      extends Modelica.Icons.Example;

      parameter Boolean 'inclC+'=true
        "<html>Carbon plus (C<sup>+</sup>)</html>" annotation (choices(
            __Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean 'inclC19HF37O5S-'=false
        "<html>Nafion sulfonate (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean 'incle-'=false "<html>Electrons (e<sup>-</sup>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean 'inclH+'=false "<html>Protons (H<sup>+</sup>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclH2=true "<html>Hydrogen (H<sub>2</sub>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclH2O=false
        "<html>Water vapor (H<sub>2</sub>O)</html>" annotation (choices(
            __Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclN2=false "<html>Nitrogen (N<sub>2</sub>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclO2=false "<html>Oxygen (O<sub>2</sub>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));

      Subregion subregion(
        L={1,1,1}*U.cm,
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        ionomer(
          final 'inclC19HF37O5S-'='inclC19HF37O5S-',
          final 'inclH+'='inclH+',
          'C19HF37O5S-'(V_IC=subregion.V/4)),
        inclFacesX=false,
        inclFacesY=false,
        inclFacesZ=false,
        graphite(
          final 'inclC+'='inclC+',
          final 'incle-'='incle-',
          'C+'(
            V_IC=subregion.V/4,
            setVelX=true,
            setVelY=true,
            setVelZ=true)))
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      inner FCSys.Conditions.Environment environment(analysis=false)
        annotation (Placement(transformation(extent={{30,20},{50,40}})));
      annotation (
        Placement(transformation(extent={{70,70},{90,90}})),
        experiment(StopTime=10),
        Commands(file(ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionH2Alone.mos"),

        experimentSetupOutput);

    end SubregionCplusAlone;

    model SubregionCAndH2
      extends Examples.Subregion('inclC+'=true, subregion(liquid(inclH2O=true)));

    end SubregionCAndH2;

    model SubregionEvaporation
      "<html>**Test a subregion with the hydrogen oxidation reaction and the essential species for it (C+, C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S, e<sup>-</sup>, H<sub>2</sub>, and H<sup>+</sup>)</html>"

      extends Examples.Subregion(
        inclH2O=true,
        inclH2=false,
        subregion(
          inclFacesX=true,
          gas(H2O(initMaterial=InitScalar.None)),
          liquid(inclH2O=inclH2O)));

      Conditions.FaceBus.SubregionFlows reactionBC(
        graphite(
          final 'inclC+'='inclC+',
          final 'incle-'='incle-',
          'e-'(redeclare Conditions.Face.Material.Current material(redeclare
                Modelica.Blocks.Sources.Ramp source(duration=1000, height=0.001
                    *U.A)))),
        ionomer(
          final 'inclC19HF37O5S-'='inclC19HF37O5S-',
          final 'inclH+'='inclH+',
          'H+'(redeclare Conditions.Face.Material.Pressure material(source(k=U.atm)))),

        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        liquid(inclH2O=inclH2O)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={20,0})));

      extends Modelica.Icons.UnderConstruction;
      // **fails sim

    equation
      connect(subregion.xPositive, reactionBC.face) annotation (Line(
          points={{10,6.10623e-16},{14,6.10623e-16},{14,1.23436e-15},{16,
              1.23436e-15}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (experiment(StopTime=1000, Tolerance=1e-06), Commands(file(
              ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionHOR.mos"));
    end SubregionEvaporation;

  end Examples;

  model Subregion "Subregion with all phases"

    parameter Boolean inclReact=true "Include reactions (as appropriate)"
      annotation (
      HideResult=true,
      Dialog(tab="Assumptions"),
      choices(__Dymola_checkBox=true));
    // Note:  This is listed above the extends clause so that it's listed
    // first in the parameter dialog.
    extends BaseClasses.PartialSubregion;

    Phases.Gas gas(inclH2O=true,final inclReact=inclReact) "Gas" annotation (
        Dialog(group="Phases (click to edit)"), Placement(transformation(extent
            ={{-10,-10},{10,10}})));
    Phases.Graphite graphite('e-'(initMaterial=if inclHOR or inclORR then
            InitScalar.None else InitScalar.Pressure)) "Graphite" annotation (
        Dialog(group="Phases (click to edit)"), Placement(transformation(extent
            ={{-10,-10},{10,10}})));
    Phases.Ionomer ionomer('C19HF37O5S-'(
        setVelX=not graphite.'inclC+',
        setVelY=not graphite.'inclC+',
        setVelZ=not graphite.'inclC+',
        initTransX=if graphite.'inclC+' then InitTranslational.None else
            InitTranslational.Velocity,
        initTransY=if graphite.'inclC+' then InitTranslational.None else
            InitTranslational.Velocity,
        initTransZ=if graphite.'inclC+' then InitTranslational.None else
            InitTranslational.Velocity)) "Ionomer" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-10,-10},
              {10,10}})));

    Phases.Liquid liquid "Liquid" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-10,-10},
              {10,10}})));

  protected
    final parameter Boolean inclHOR=inclReact and (graphite.'incle-' and
        ionomer.'inclH+' and gas.inclH2 and not (gas.inclO2 and gas.inclH2O))
      "true, if HOR is included" annotation (HideResult=true);
    final parameter Boolean inclORR=inclReact and (graphite.'incle-' and
        ionomer.'inclH+' and gas.inclO2 and gas.inclH2O and not gas.inclH2)
      "true, if ORR is included" annotation (HideResult=true);
    Reaction hOR(n_neut=1) if inclHOR "Hydrogen oxidation reaction" annotation
      (Placement(transformation(extent={{-66.66,33.33},{-46.66,53.33}})));
    Reaction oRR(n_neut=2) if inclORR "Oxygen reduction reaction"
      annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
    // Note:  The additional condition (not gas.inclH2O) prevents a singularity
    // if water is included as gas, liquid, and in ionomer.

    Connectors.ChemicalBusInternal CondEvap
      "Connector for H2O condensation and evaporation" annotation (Placement(
          transformation(extent={{-30,50},{-10,70}}), iconTransformation(extent
            ={{-76,16},{-56,36}})));
    Connectors.ChemicalBusInternal Hydration
      "Connector for ionomer hydration and drying" annotation (Placement(
          transformation(extent={{-43.33,36.66},{-23.33,56.66}}),
          iconTransformation(extent={{-52,56},{-32,76}})));
    Connectors.ChemicalBusInternal HOR
      "Connector for hydrogen oxidation reaction" annotation (Placement(
          transformation(extent={{-56.66,23.33},{-36.66,43.33}}),
          iconTransformation(extent={{-64,36},{-44,56}})));
    Connectors.ChemicalBusInternal ORR
      "Connector for oxygen reduction reaction" annotation (Placement(
          transformation(extent={{-70,10},{-50,30}}), iconTransformation(extent
            ={{-40,76},{-20,96}})));

  equation
    // Chemical reactions
    // ------------------
    // Condensation/evaporation
    connect(CondEvap, gas.CondEvap) annotation (Line(
        points={{-20,60},{-3,8.6}},
        color={208,104,0},
        pattern=LinePattern.Solid,
        thickness=0.5,
        smooth=Smooth.None));
    connect(CondEvap, liquid.CondEvap) annotation (Line(
        points={{-20,60},{-3,8.6}},
        color={208,104,0},
        pattern=LinePattern.Solid,
        thickness=0.5,
        smooth=Smooth.None));
    // HOR
    connect(HOR, gas.HOR) annotation (Line(
        points={{-46.66,33.33},{-4.9,5.7}},
        color={208,104,0},
        pattern=LinePattern.Solid,
        thickness=0.5,
        smooth=Smooth.None));
    connect(HOR, graphite.HOR) annotation (Line(
        points={{-46.66,33.33},{-4.9,5.7}},
        color={208,104,0},
        pattern=LinePattern.Solid,
        thickness=0.5,
        smooth=Smooth.None));
    connect(HOR, ionomer.HOR) annotation (Line(
        points={{-46.66,33.33},{-4.9,5.7}},
        color={208,104,0},
        pattern=LinePattern.Solid,
        thickness=0.5,
        smooth=Smooth.None));
    connect(hOR.positive, HOR.'H+') annotation (Line(
        points={{-46.66,43.33},{-46.66,33.33}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(hOR.chemical[2], HOR.H2) annotation (Line(
        points={{-56.66,43.33},{-46.66,33.33}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(hOR.negative, HOR.'e-') annotation (Line(
        points={{-66.66,43.33},{-46.66,33.33}},
        color={208,104,0},
        smooth=Smooth.None));
    // Hydration/drying
    connect(Hydration, gas.Hydration) annotation (Line(
        points={{-33.33,46.66},{-4.2,6.6}},
        color={208,104,0},
        pattern=LinePattern.Solid,
        thickness=0.5,
        smooth=Smooth.None));
    connect(Hydration, ionomer.Hydration) annotation (Line(
        points={{-33.33,46.66},{-4.2,6.6}},
        color={208,104,0},
        pattern=LinePattern.Solid,
        thickness=0.5,
        smooth=Smooth.None));
    // ORR
    connect(ORR, gas.ORR) annotation (Line(
        points={{-60,20},{-6.4,2.8}},
        color={208,104,0},
        pattern=LinePattern.Solid,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ORR, graphite.ORR) annotation (Line(
        points={{-60,20},{-6.4,2.8}},
        color={208,104,0},
        pattern=LinePattern.Solid,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ORR, ionomer.ORR) annotation (Line(
        points={{-60,20},{-6.4,2.8}},
        color={208,104,0},
        pattern=LinePattern.Solid,
        thickness=0.5,
        smooth=Smooth.None));
    connect(oRR.positive, ORR.'e-') annotation (Line(
        points={{-60,30},{-60,20}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(oRR.chemical[1], ORR.O2) annotation (Line(
        points={{-70,29.5},{-60,20}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(oRR.chemical[2], ORR.H2O) annotation (Line(
        points={{-70,30.5},{-60,20}},
        color={208,104,0},
        smooth=Smooth.None));
    connect(oRR.negative, ORR.'H+') annotation (Line(
        points={{-80,30},{-60,20}},
        color={208,104,0},
        smooth=Smooth.None));

    // Gas
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
        points={{8,6.10623e-16},{8,-4.87687e-22},{40,-4.87687e-22},{40,
            5.55112e-16}},
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
        points={{6.10623e-16,10},{0,10},{0,40},{5.55112e-16,40}},
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
        points={{-8,6.10623e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.xPositive, xPositive.graphite) annotation (Line(
        points={{8,6.10623e-16},{8,-4.87687e-22},{40,-4.87687e-22},{40,
            5.55112e-16}},
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
        points={{6.10623e-16,10},{0,10},{0,40},{5.55112e-16,40}},
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
        points={{8,6.10623e-16},{8,-4.87687e-22},{40,-4.87687e-22},{40,
            5.55112e-16}},
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
    connect(liquid.inert, volume.inert) annotation (Line(
        points={{8,-8},{11,-11}},
        color={0,180,0},
        smooth=Smooth.None,
        thickness=0.5));
    connect(liquid.xNegative, xNegative.liquid) annotation (Line(
        points={{-8,6.10623e-16},{-8,6.10623e-16},{-8,5.55112e-16},{-40,
            5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(liquid.xPositive, xPositive.liquid) annotation (Line(
        points={{8,6.10623e-16},{8,-4.87687e-22},{40,-4.87687e-22},{40,
            5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(liquid.yNegative, yNegative.liquid) annotation (Line(
        points={{6.10623e-16,-8.4},{-4.87687e-22,-8.4},{-4.87687e-22,-40},{
            5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.yPositive, yPositive.liquid) annotation (Line(
        points={{6.10623e-16,10},{0,10},{0,10},{-4.87687e-22,40},{5.55112e-16,
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
  <li>H<sub>2</sub>O vapor is included by default, since at least one species
  must be included.</li></ul></p>

<p>For more information, see the
   <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSubregion\">PartialSubregion</a> model.</p></html>"),

      Diagram(graphics));
  end Subregion;

  model SubregionIonomerOnly "Subregion with only the ionomer phase"
    extends BaseClasses.PartialSubregion;

    Phases.Ionomer ionomer(inclH2O=true, final inclTrans={inclTransX,inclTransY,
          inclTransZ}) "Ionomer" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-10,-10},
              {10,10}})));

  equation
    // Ionomer
    connect(ionomer.chemical, chemical.ionomer);
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
      defaultComponentName="subregion",
      Documentation(info="<html><p>Notes:<ul>
  <li>H<sub>2</sub>O is included by default, since at least one species
  must be included.</li></ul></p>

<p>For more information, see the
   <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSubregion\">PartialSubregion</a> model.</p></html>"),

      Diagram(graphics));
  end SubregionIonomerOnly;

  model SubregionNoIonomer "Subregion with all phases except ionomer"
    parameter Boolean inclReact=false "Include reactions (as appropriate)"
      annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(tab="Assumptions"));

    // Note:  This is listed above the extends clause so that it's listed
    // first in the parameter dialog.
    extends BaseClasses.PartialSubregion;

    Phases.Gas gas(inclH2O=true, inclReact=inclReact) "Gas" annotation (Dialog(
          group="Phases (click to edit)"), Placement(transformation(extent={{-10,
              -10},{10,10}})));

    Phases.Graphite graphite "Graphite" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-10,-10},
              {10,10}})));
    Phases.Liquid liquid "Liquid" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-10,-10},
              {10,10}})));

  protected
    Connectors.ChemicalBusInternal CondEvap
      "Connector for H2O condensation and evaporation" annotation (Placement(
          transformation(extent={{-30,10},{-10,30}}), iconTransformation(extent
            ={{-76,16},{-56,36}})));

  equation
    // Chemical interactions (not shown graphically)

    // Gas
    connect(gas.chemical, chemical.gas);

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
    connect(graphite.chemical, chemical.graphite);

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
    connect(liquid.chemical, chemical.liquid);

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

    connect(CondEvap, liquid.CondEvap) annotation (Line(
        points={{-20,20},{-3,8.6}},
        color={208,104,0},
        pattern=LinePattern.Solid,
        thickness=0.5,
        smooth=Smooth.None));
    annotation (
      defaultComponentPrefixes="replaceable",
      defaultComponentName="subregion",
      Documentation(info="<html><p>Notes:<ul>
    <li>The <code>inclReact</code> parameter may be set to
    <code>false</code>
    to eliminate unnecessary equations.</li>
  <li>H<sub>2</sub>O vapor is included by default, since at least one species
  must be included.</li></ul></p>

<p>For more information, see the
   <a href=\"modelica://FCSys.Subregions.BaseClasses.PartialSubregion\">PartialSubregion</a> model.</p></html>"),

      Diagram(graphics));
  end SubregionNoIonomer;

  package Phases "Phases or mixtures of species"
    extends Modelica.Icons.Package;

    model Gas "Gas phase"
      import FCSys.BaseClasses.Utilities.countTrue;

      extends BaseClasses.NullPhase(final n_spec=countTrue({inclH2,inclH2O,
            inclN2,inclO2}));

      parameter Boolean inclReact=true "Include reactions (as appropriate)"
        annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(tab="Assumptions"));

      // Conditionally include species.
      parameter Boolean inclH2=true "<html>Hydrogen (H<sub>2</sub>)</html>"
        annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));

      replaceable Species.H2.Gas.Fixed H2 if inclH2 constrainedby
        Species.Species(
        n_react=1,
        initTransX=if initVelX and reduceVel then InitTranslational.None else
            InitTranslational.Velocity,
        initTransY=if initVelY and reduceVel then InitTranslational.None else
            InitTranslational.Velocity,
        initTransZ=if initVelZ and reduceVel then InitTranslational.None else
            InitTranslational.Velocity,
        initEnergy=if initTemp and reduceTemp then InitScalar.None else
            InitScalar.Temperature,
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
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));

      replaceable Species.H2O.Gas.Fixed H2O if inclH2O constrainedby
        Species.Species(
        n_react=2,
        initTransX=if initVelX and reduceVel then InitTranslational.None else
            InitTranslational.Velocity,
        initTransY=if initVelY and reduceVel then InitTranslational.None else
            InitTranslational.Velocity,
        initTransZ=if initVelZ and reduceVel then InitTranslational.None else
            InitTranslational.Velocity,
        initEnergy=if initTemp and reduceTemp then InitScalar.None else
            InitScalar.Temperature,
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
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));

      replaceable Species.N2.Gas.Fixed N2 if inclN2 constrainedby
        Species.Species(
        n_react=0,
        initTransX=if initVelX and reduceVel then InitTranslational.None else
            InitTranslational.Velocity,
        initTransY=if initVelY and reduceVel then InitTranslational.None else
            InitTranslational.Velocity,
        initTransZ=if initVelZ and reduceVel then InitTranslational.None else
            InitTranslational.Velocity,
        initEnergy=if initTemp and reduceTemp then InitScalar.None else
            InitScalar.Temperature,
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
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.O2.Gas.Fixed O2 if inclO2 constrainedby
        Species.Species(
        n_react=1,
        initTransX=if initVelX and reduceVel then InitTranslational.None else
            InitTranslational.Velocity,
        initTransY=if initVelY and reduceVel then InitTranslational.None else
            InitTranslational.Velocity,
        initTransZ=if initVelZ and reduceVel then InitTranslational.None else
            InitTranslational.Velocity,
        initEnergy=if initTemp and reduceTemp then InitScalar.None else
            InitScalar.Temperature,
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

      Connectors.ChemicalBus PC "Phase change" annotation (Placement(
            transformation(extent={{-30,50},{-10,70}}), iconTransformation(
              extent={{-42,74},{-22,94}})));
      Connectors.ChemicalBus HOR "Hydrogen oxidation reaction" annotation (
          Placement(transformation(extent={{-49.66,30.33},{-29.66,50.33}}),
            iconTransformation(extent={{-59,47},{-39,67}})));
      Connectors.ChemicalBus ORR "Oxygen reduction reaction" annotation (
          Placement(transformation(extent={{-70,10},{-50,30}}),
            iconTransformation(extent={{-74,18},{-54,38}})));
    protected
      Connectors.ChemicalSpecies chemical4(n_trans=n_trans)
        annotation (Placement(transformation(extent={{-36,18},{-16,38}})));
      Connectors.ChemicalSpecies chemical1(n_trans=n_trans)
        annotation (Placement(transformation(extent={{-16,46},{4,66}})));
      Connectors.ChemicalSpecies chemical2(n_trans=n_trans)
        annotation (Placement(transformation(extent={{-30,34},{-10,54}})));
      Connectors.ChemicalSpecies chemical3(n_trans=n_trans)
        annotation (Placement(transformation(extent={{-50,0},{-30,20}})));
    equation
      // Phase change
      connect(chemical1, H2O.chemical[1]) annotation (Line(
          points={{-6,56},{-5.578,5.555}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(PC.H2O, chemical1) annotation (Line(
          points={{-20,60},{-6,56}},
          color={208,104,0},
          smooth=Smooth.None));
      // HOR
      connect(chemical2, H2.chemical[1]) annotation (Line(
          points={{-20,44},{-5.578,5.555}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(HOR.H2, chemical2) annotation (Line(
          points={{-39.66,40.33},{-20,44}},
          color={208,104,0},
          smooth=Smooth.None));
      // ORR
      connect(chemical3, H2O.chemical[2]) annotation (Line(
          points={{-40,10},{-5.578,5.555}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(ORR.H2O, chemical3) annotation (Line(
          points={{-60,20},{-40,10}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(chemical4, O2.chemical[1]) annotation (Line(
          points={{-26,28},{-5.578,5.555}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(ORR.O2, chemical4) annotation (Line(
          points={{-60,20},{-26,28}},
          color={208,104,0},
          smooth=Smooth.None));

      // H2
      // --
      // Exchange
      connect(H2.inert.translational, inert.translational) annotation (Line(
          points={{7.155,-3.578},{26.67,-13.33}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2.inert.thermal, inert.thermal) annotation (Line(
          points={{7.155,-3.578},{26.67,-13.33}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2.inertDalton, phaseBoundary.inertDalton) annotation (Line(
          points={{3.578,-7.155},{4,-8},{4,-7.155},{3.578,-7.155}},
          color={72,90,180},
          smooth=Smooth.None));
      // Transport
      connect(H2.faces[Axis.x, Side.n], xNegative.H2) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.faces[Axis.x, Side.p], xPositive.H2) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},
              {40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.faces[Axis.y, Side.n], yNegative.H2) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{6.10623e-16,-40},{5.55112e-16,-40}},

          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.faces[Axis.y, Side.p], yPositive.H2) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.faces[Axis.z, Side.n], zNegative.H2) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.faces[Axis.z, Side.p], zPositive.H2) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // H2O
      // ---
      // Exchange
      connect(H2O.inert.translational, inert.translational) annotation (Line(
          points={{7.155,-3.578},{26.67,-13.33}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2O.inert.thermal, inert.thermal) annotation (Line(
          points={{7.155,-3.578},{26.67,-13.33}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2O.inertDalton, phaseBoundary.inertDalton) annotation (Line(
          points={{3.578,-7.155},{4,-8},{4,-7.155},{3.578,-7.155}},
          color={72,90,180},
          smooth=Smooth.None));
      // Transport
      connect(H2O.faces[Axis.x, Side.n], xNegative.H2O) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[Axis.x, Side.p], xPositive.H2O) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},
              {40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[Axis.y, Side.n], yNegative.H2O) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{6.10623e-16,-40},{5.55112e-16,-40}},

          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[Axis.y, Side.p], yPositive.H2O) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[Axis.z, Side.n], zNegative.H2O) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.faces[Axis.z, Side.p], zPositive.H2O) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // N2
      // --
      // Exchange
      connect(N2.inert.translational, inert.translational) annotation (Line(
          points={{7.155,-3.578},{26.67,-13.33}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(N2.inert.thermal, inert.thermal) annotation (Line(
          points={{7.155,-3.578},{26.67,-13.33}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(N2.inertDalton, phaseBoundary.inertDalton) annotation (Line(
          points={{3.578,-7.155},{4,-8},{4,-7.155},{3.578,-7.155}},
          color={72,90,180},
          smooth=Smooth.None));
      // Transport
      connect(N2.faces[Axis.x, Side.n], xNegative.N2) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.faces[Axis.x, Side.p], xPositive.N2) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},
              {40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.faces[Axis.y, Side.n], yNegative.N2) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{6.10623e-16,-40},{5.55112e-16,-40}},

          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.faces[Axis.y, Side.p], yPositive.N2) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.faces[Axis.z, Side.n], zNegative.N2) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.faces[Axis.z, Side.p], zPositive.N2) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // O2
      // --
      // Exchange
      connect(O2.inert.translational, inert.translational) annotation (Line(
          points={{7.155,-3.578},{26.67,-13.33}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(O2.inert.thermal, inert.thermal) annotation (Line(
          points={{7.155,-3.578},{26.67,-13.33}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(O2.inertDalton, phaseBoundary.inertDalton) annotation (Line(
          points={{3.578,-7.155},{3.578,-7.155},{3.578,-7.155},{3.578,-7.155},{
              3.578,-7.155},{3.578,-7.155}},
          color={72,90,180},
          smooth=Smooth.None));

      // Transport
      connect(O2.faces[Axis.x, Side.n], xNegative.O2) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.faces[Axis.x, Side.p], xPositive.O2) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.faces[Axis.y, Side.n], yNegative.O2) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.82867e-16,6.10623e-16},{
              5.82867e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.faces[Axis.y, Side.p], yPositive.O2) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.faces[Axis.z, Side.n], zNegative.O2) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.faces[Axis.z, Side.p], zPositive.O2) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      annotation (
        defaultComponentPrefixes="replaceable",
        Documentation(info="<html><p>Notes:<ul><li>The <code>inclReact</code> parameter may be set to
    <code>false</code>
    to eliminate unnecessary equations.</li></ul></p>

<p>For more information, see the
 <a href=\"modelica://FCSys.Subregions.Phases.BaseClasses.NullPhase\">NullPhase</a> model.</p></html>"),

        Icon(graphics),
        Diagram(graphics));
    end Gas;

    model Graphite "Graphite phase"
      import FCSys.BaseClasses.Utilities.countTrue;
      extends BaseClasses.NullPhase(
        final n_spec=countTrue({'inclC+','incle-'}),
        initVelX=not 'inclC+',
        initVelY=not 'inclC+',
        initVelZ=not 'inclC+');

      // Conditionally include species.
      parameter Boolean 'inclC+'=false
        "<html>Carbon plus (C<sup>+</sup>)</html>" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));

      replaceable Species.'C+'.Graphite.Fixed 'C+' if 'inclC+' constrainedby
        Species.Species(
        n_react=0,
        initEnergy=if initTemp and reduceTemp then InitScalar.None else
            InitScalar.Temperature,
        phi(each stateSelect=if reduceVel then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceTemp then StateSelect.default else StateSelect.prefer))
        "<html>C<sup>+</sup> model</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable='inclC+'),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      parameter Boolean 'incle-'=false "<html>Electrons (e<sup>-</sup>)</html>"
        annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));

      replaceable Species.'e-'.Graphite.Fixed 'e-' if 'incle-' constrainedby
        Species.Species(
        n_react=2,
        initTransX=if (initTransX or 'inclC+') and reduceVel then
            InitTranslational.None else InitTranslational.Velocity,
        initTransY=if (initTransY or 'inclC+') and reduceVel then
            InitTranslational.None else InitTranslational.Velocity,
        initTransZ=if (initTransZ or 'inclC+') and reduceVel then
            InitTranslational.None else InitTranslational.Velocity,
        initEnergy=if initTemp and reduceTemp then InitScalar.None else
            InitScalar.Temperature,
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

      Connectors.ChemicalBus HOR "Hydrogen oxidation reaction" annotation (
          Placement(transformation(extent={{-49.66,30.33},{-29.66,50.33}}),
            iconTransformation(extent={{-59,47},{-39,67}})));
      Connectors.ChemicalBus ORR "Oxygen reduction reaction" annotation (
          Placement(transformation(extent={{-70,10},{-50,30}}),
            iconTransformation(extent={{-74,18},{-54,38}})));

    equation
      // HOR
      connect(HOR.'e-', 'e-'.chemical[1]) annotation (Line(
          points={{-39.66,40.33},{-5.578,5.555}},
          color={208,104,0},
          smooth=Smooth.None));
      // ORR
      connect(ORR.'e-', 'e-'.chemical[2]) annotation (Line(
          points={{-60,20},{-5.578,5.555}},
          color={208,104,0},
          smooth=Smooth.None));

      // C+
      // --
      // Exchange
      connect('C+'.inert.translational, inert.translational) annotation (Line(
          points={{7.155,-3.578},{26.67,-13.33}},
          color={72,90,180},
          smooth=Smooth.None));
      connect('C+'.inert.thermal, inert.thermal) annotation (Line(
          points={{7.155,-3.578},{26.67,-13.33}},
          color={72,90,180},
          smooth=Smooth.None));
      connect('C+'.inertDalton, phaseBoundary.inertD) annotation (Line(
          points={{7.155,-3.578},{3.578,-7.155},{3.578,-7.155}},
          color={72,90,180},
          smooth=Smooth.None));
      // Transport
      connect('C+'.faces[Axis.x, Side.n], xNegative.'C+') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('C+'.faces[Axis.x, Side.p], xPositive.'C+') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{6.10623e-16,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('C+'.faces[Axis.y, Side.n], yNegative.'C+') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-4.87687e-22,6.10623e-16},{-4.87687e-22,
              -40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('C+'.faces[Axis.y, Side.p], yPositive.'C+') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('C+'.faces[Axis.z, Side.n], zNegative.'C+') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('C+'.faces[Axis.z, Side.p], zPositive.'C+') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // e-
      // --
      // Exchange
      connect('e-'.inert.translational, inert.translational) annotation (Line(
          points={{7.155,-3.578},{26.67,-13.33}},
          color={72,90,180},
          smooth=Smooth.None));
      connect('e-'.inert.thermal, inert.thermal) annotation (Line(
          points={{7.155,-3.578},{26.67,-13.33}},
          color={72,90,180},
          smooth=Smooth.None));
      connect('e-'.inertDalton, phaseBoundary.inertD) annotation (Line(
          points={{7.155,-3.578},{3.578,-7.155},{3.578,-7.155}},
          color={72,90,180},
          smooth=Smooth.None));
      // Transport
      connect('e-'.faces[Axis.x, Side.n], xNegative.'e-') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.faces[Axis.x, Side.p], xPositive.'e-') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{6.10623e-16,5.55112e-16},{40,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.faces[Axis.y, Side.n], yNegative.'e-') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-4.87687e-22,6.10623e-16},{-4.87687e-22,
              -40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.faces[Axis.y, Side.p], yPositive.'e-') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.faces[Axis.z, Side.n], zNegative.'e-') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.faces[Axis.z, Side.p], zPositive.'e-') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      annotation (
        defaultComponentPrefixes="replaceable",
        Documentation(info="<html><p>If C<sup>+</sup> is included (<code>'inclC+'</code>=<code>true</code>),
    then <code>initTransX</code>, <code>initTransY</code>, and <code>initTransZ</code> are
    ignored.  Velocity is not initialized because it is set by C<sup>+</sup>. **update this</p>

    <p>Assumptions:
    <ol>
    <li>The density of e<sup>-</sup> is equal to that of C<sup>+</sup>.</li>
    <li>All of the C has been ionized to C<sup>+</sup>.</li>
    </ol></p>

    <p>For more information, see the
 <a href=\"modelica://FCSys.Subregions.Phases.BaseClasses.NullPhase\">NullPhase</a> model.</p></html>"),

        Diagram(graphics),
        Icon(graphics));
    end Graphite;

    model Ionomer "Ionomer phase"
      import FCSys.BaseClasses.Utilities.countTrue;
      extends BaseClasses.NullPhase(
        final n_spec=countTrue({'inclC19HF37O5S-','inclH+',inclH2O}),
        initVelX=not 'inclC19HF37O5S-',
        initVelY=not 'inclC19HF37O5S-',
        initVelZ=not 'inclC19HF37O5S-');

      // Conditionally include species.
      parameter Boolean 'inclC19HF37O5S-'=false
        "<html>Nafion sulfonate minus (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>)</html>"
        annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));

      replaceable Species.'C19HF37O5S-'.Ionomer.Fixed 'C19HF37O5S-' if
        'inclC19HF37O5S-' constrainedby Species.Species(
        n_react=0,
        initEnergy=if initTemp and reduceTemp then InitScalar.None else
            InitScalar.Temperature,
        phi(each stateSelect=if reduceVel then StateSelect.default else
              StateSelect.prefer),
        T(stateSelect=if reduceTemp then StateSelect.default else StateSelect.prefer))
        "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup> model</html>"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable='inclC19HF37O5S-'),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      parameter Boolean 'inclH+'=false "<html>Protons (H<sup>+</sup>)</html>"
        annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));

      replaceable Species.'H+'.Ionomer.Fixed 'H+' if 'inclH+' constrainedby
        Species.Species(
        n_react=2,
        initTransX=if initVelX and reduceVel then InitTranslational.None else
            InitTranslational.Velocity,
        initTransY=if initVelY and reduceVel then InitTranslational.None else
            InitTranslational.Velocity,
        initTransZ=if initVelZ and reduceVel then InitTranslational.None else
            InitTranslational.Velocity,
        initEnergy=if initTemp and reduceTemp then InitScalar.None else
            InitScalar.Temperature,
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

      parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
        annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));

      replaceable Species.H2O.Ionomer.Fixed H2O if inclH2O constrainedby
        Species.Species(
        n_react=1,
        initTransX=if (initTransX or 'inclC19HF37O5S-') and reduceVel then
            InitTranslational.None else InitTranslational.Velocity,
        initTransY=if (initTransY or 'inclC19HF37O5S-') and reduceVel then
            InitTranslational.None else InitTranslational.Velocity,
        initTransZ=if (initTransZ or 'inclC19HF37O5S-') and reduceVel then
            InitTranslational.None else InitTranslational.Velocity,
        initEnergy=if initTemp and reduceTemp then InitScalar.None else
            InitScalar.Temperature,
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

      Connectors.ChemicalBus PC "Phase change" annotation (Placement(
            transformation(extent={{-30,50},{-10,70}}), iconTransformation(
              extent={{-42,74},{-22,94}})));
      Connectors.ChemicalBus HOR "Hydrogen oxidation reaction" annotation (
          Placement(transformation(extent={{-49.66,30.33},{-29.66,50.33}}),
            iconTransformation(extent={{-59,47},{-39,67}})));
      Connectors.ChemicalBus ORR "Oxygen reduction reaction" annotation (
          Placement(transformation(extent={{-70,10},{-50,30}}),
            iconTransformation(extent={{-74,18},{-54,38}})));

    equation
      // Phase change
      connect(PC.H2O, H2O.chemical[1]) annotation (Line(
          points={{-20,60},{-5.578,5.555}},
          color={208,104,0},
          smooth=Smooth.None));
      // HOR
      connect(HOR.'H+', 'H+'.chemical[1]) annotation (Line(
          points={{-39.66,40.33},{-5.578,5.555}},
          color={208,104,0},
          smooth=Smooth.None));
      // ORR
      connect(ORR.'H+', 'H+'.chemical[2]) annotation (Line(
          points={{-60,20},{-5.578,5.555}},
          color={208,104,0},
          smooth=Smooth.None));

      // C19HF37O5S-
      // -----------
      // Exchange
      connect('C19HF37O5S-'.inert.translational, inert.translational)
        annotation (Line(
          points={{7.155,-3.578},{26.67,-13.33}},
          color={72,90,180},
          smooth=Smooth.None));
      connect('C19HF37O5S-'.inert.thermal, inert.thermal) annotation (Line(
          points={{7.155,-3.578},{26.67,-13.33}},
          color={72,90,180},
          smooth=Smooth.None));
      connect('C19HF37O5S-'.inertDalton, phaseBoundary.inertDalton) annotation
        (Line(
          points={{3.578,-7.155},{7,-7},{7,-7.155},{3.578,-7.155}},
          color={72,90,180},
          smooth=Smooth.None));
      // Transport
      connect('C19HF37O5S-'.faces[Axis.x, Side.n], xNegative.'C19HF37O5S-')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},
              {-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('C19HF37O5S-'.faces[Axis.x, Side.p], xPositive.'C19HF37O5S-')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},
              {40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('C19HF37O5S-'.faces[Axis.y, Side.n], yNegative.'C19HF37O5S-')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{6.10623e-16,-40},{5.55112e-16,-40}},

          color={127,127,127},
          smooth=Smooth.None));

      connect('C19HF37O5S-'.faces[Axis.y, Side.p], yPositive.'C19HF37O5S-')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-4.87687e-22,20},{5.55112e-16,20},
              {5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('C19HF37O5S-'.faces[Axis.z, Side.n], zNegative.'C19HF37O5S-')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('C19HF37O5S-'.faces[Axis.z, Side.p], zPositive.'C19HF37O5S-')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // 'H+'
      // ----
      // Exchange
      connect('H+'.inert.translational, inert.translational) annotation (Line(
          points={{7.155,-3.578},{26.67,-13.33}},
          color={72,90,180},
          smooth=Smooth.None));
      connect('H+'.inert.thermal, inert.thermal) annotation (Line(
          points={{7.155,-3.578},{26.67,-13.33}},
          color={72,90,180},
          smooth=Smooth.None));
      connect('H+'.inertDalton, phaseBoundary.inertDalton) annotation (Line(
          points={{3.578,-7.155},{7,-7},{7,-7.155},{3.578,-7.155}},
          color={72,90,180},
          smooth=Smooth.None));
      // Transport
      connect('H+'.faces[Axis.x, Side.n], xNegative.'H+') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},
              {-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.faces[Axis.x, Side.p], xPositive.'H+') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},
              {40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.faces[Axis.y, Side.n], yNegative.'H+') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{6.10623e-16,-40},{5.55112e-16,-40}},

          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.faces[Axis.y, Side.p], yPositive.'H+') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-4.87687e-22,20},{5.55112e-16,20},
              {5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.faces[Axis.z, Side.n], zNegative.'H+') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.faces[Axis.z, Side.p], zPositive.'H+') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // H2O
      // ---
      // Exchange
      connect(H2O.inert.translational, inert.translational) annotation (Line(
          points={{7.155,-3.578},{26.67,-13.33}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2O.inert.thermal, inert.thermal) annotation (Line(
          points={{7.155,-3.578},{26.67,-13.33}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2O.inertDalton, phaseBoundary.inertDalton) annotation (Line(
          points={{3.578,-7.155},{7,-7},{7,-7.155},{3.578,-7.155}},
          color={72,90,180},
          smooth=Smooth.None));
      // Transport
      connect(H2O.faces[Axis.x, Side.n], xNegative.H2O) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},
              {-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[Axis.x, Side.p], xPositive.H2O) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},
              {40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[Axis.y, Side.n], yNegative.H2O) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{6.10623e-16,-40},{5.55112e-16,-40}},

          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[Axis.y, Side.p], yPositive.H2O) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-4.87687e-22,20},{5.55112e-16,20},
              {5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[Axis.z, Side.n], zNegative.H2O) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.faces[Axis.z, Side.p], zPositive.H2O) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      annotation (
        defaultComponentPrefixes="replaceable",
        Documentation(info="<html><p>If C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup> is included
    (<code>'inclC19HF37O5S-'</code>=<code>true</code>),
    then <code>initTransX</code>, <code>initTransY</code>, and <code>initTransZ</code> are
    ignored.  Velocity is not initialized because it is set by C19HF37O5S.**update this</p>

    <p>Assumptions:
    <ol>
    <li>The density of H<sup>+</sup> is equal to that of C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>.</li>
    <li>All of the C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S has been ionized to
    C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>.  Note that this is not true in practice
    (see <a href=\"modelica://FCSys.Characteristics.'H+'.Ionomer\">Characteristics.'H+'.Ionomer</a>).</li>
    </ol></p>

    <p>For more information, see the
 <a href=\"modelica://FCSys.Subregions.Phases.BaseClasses.NullPhase\">NullPhase</a> model.</p></html>"),

        Diagram(graphics),
        Icon(graphics));
    end Ionomer;

    model Liquid "Liquid phase"

      extends BaseClasses.NullPhase(final n_spec=if inclH2O then 1 else 0);

      // Conditionally include species.
      parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
        annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.H2O.Liquid.Fixed H2O if inclH2O constrainedby
        Species.Species(
        n_react=1,
        initTransX=if initVelX and reduceVel then InitTranslational.None else
            InitTranslational.Velocity,
        initTransY=if initVelY and reduceVel then InitTranslational.None else
            InitTranslational.Velocity,
        initTransZ=if initVelZ and reduceVel then InitTranslational.None else
            InitTranslational.Velocity,
        initEnergy=if initTemp and reduceTemp then InitScalar.None else
            InitScalar.Temperature,
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

      Connectors.ChemicalBus PC "Phase change" annotation (Placement(
            transformation(extent={{-30,50},{-10,70}}), iconTransformation(
              extent={{-42,74},{-22,94}})));

    equation
      // Phase change
      connect(PC.H2O, H2O.chemical[1]) annotation (Line(
          points={{-20,60},{-5.578,5.555}},
          color={208,104,0},
          smooth=Smooth.None));

      // H2O
      // ---
      // Exchange
      connect(H2O.inert.translational, inert.translational) annotation (Line(
          points={{7.155,-3.578},{26.67,-13.33}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2O.inert.thermal, inert.thermal) annotation (Line(
          points={{7.155,-3.578},{26.67,-13.33}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(H2O.inertDalton, phaseBoundary.inertD) annotation (Line(
          points={{3.578,-7.155},{3.8,-7},{3.8,-7.155},{3.578,-7.155}},
          color={72,90,180},
          smooth=Smooth.None));
      // Transport
      connect(H2O.faces[Axis.x, Side.n], xNegative.H2O) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},
              {-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[Axis.x, Side.p], xPositive.H2O) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},
              {40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[Axis.y, Side.n], yNegative.H2O) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{6.10623e-16,-40},{5.55112e-16,-40}},

          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[Axis.y, Side.p], yPositive.H2O) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-4.87687e-22,20},{5.55112e-16,20},
              {5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[Axis.z, Side.n], zNegative.H2O) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.faces[Axis.z, Side.p], zPositive.H2O) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      annotation (
        defaultComponentPrefixes="replaceable",
        Documentation(info="<html><p>See the information in the
 <a href=\"modelica://FCSys.Subregions.Phases.BaseClasses.NullPhase\">NullPhase</a> model.</p></html>"),

        Icon(graphics),
        Diagram(graphics));
    end Liquid;

    package BaseClasses "Base classes (not generally for direct use)"
      extends Modelica.Icons.BasesPackage;
      model NullPhase "Model for a phase with no species or reactions"
        import FCSys.BaseClasses.Utilities.index;
        // extends FCSys.BaseClasses.Icons.Names.Middle;

        // Geometric parameters
        parameter Q.NumberAbsolute k[Axis](
          each min=Modelica.Constants.small,
          each final nominal=1) = {1,1,1} if n_spec > 0
          "<html>Area fill factor (<b><i>k</i></b>)</html>"
          annotation (Dialog(group="Geometry"));

        // Assumptions
        parameter Boolean reduceVel=true if n_spec > 0
          "Same velocity for all species" annotation (Dialog(tab="Assumptions",
              enable=n_spec > 1), choices(__Dymola_checkBox=true));
        parameter Boolean reduceTemp=true if n_spec > 0
          "Same temperature for all species" annotation (Dialog(tab=
                "Assumptions", enable=n_spec > 1), choices(__Dymola_checkBox=
                true));

        // Initialization
        parameter Boolean initVelX=true if n_spec > 0
          "Initialize the x component" annotation (
          compact=true,
          Dialog(
            tab="Initialization",
            group="Velocity",
            enable=reduceVel),
          __Dymola_joinNext=true,
          choices(__Dymola_checkBox=true));
        parameter Boolean initVelY=true if n_spec > 0
          "Initialize the y component" annotation (
          compact=true,
          Dialog(
            tab="Initialization",
            group="Velocity",
            enable=reduceVel),
          __Dymola_joinNext=true,
          choices(__Dymola_checkBox=true));
        parameter Boolean initVelZ=true if n_spec > 0
          "Initialize the z component" annotation (
          compact=true,
          Dialog(
            tab="Initialization",
            group="Velocity",
            enable=reduceVel),
          __Dymola_joinNext=true,
          choices(__Dymola_checkBox=true));
        parameter Q.Velocity phi_IC[Axis]={0,0,0} if n_spec > 0
          "<html>Initial velocity (<b>&phi;</b><sub>IC</sub>)</html>"
          annotation (Dialog(
            tab="Initialization",
            group="Velocity",
            enable=reduceVel));
        // This is always enabled in the dialog since it's used as a guess value.
        parameter Boolean initTemp=true if n_spec > 0 "Initialize" annotation (
          compact=true,
          Dialog(
            tab="Initialization",
            group="Temperature",
            enable=reduceTemp),
          choices(__Dymola_checkBox=true));

        parameter Q.TemperatureAbsolute T_IC(nominal=300*U.K, start=environment.T)
          if n_spec > 0
          "<html>Initial temperature (<i>T</i><sub>IC</sub>)</html>"
          annotation (Dialog(
            tab="Initialization",
            group="Temperature",
            enable=reduceTemp));
        // This is always enabled in the dialog since it's used as a guess value.

        Connectors.InertAmagat inertAmagat(final n_trans=n_trans) if n_spec > 0
          annotation (Placement(transformation(extent={{3.33,-36.67},{23.33,-16.67}}),
              iconTransformation(extent={{70,-90},{90,-70}})));
        Connectors.FaceBus xPositive if n_spec > 0
          "Positive face along the x axis" annotation (Placement(transformation(
                extent={{30,-10},{50,10}}), iconTransformation(extent={{70,-10},
                  {90,10}})));
        Connectors.FaceBus xNegative if n_spec > 0
          "Negative face along the x axis" annotation (Placement(transformation(
                extent={{-50,-10},{-30,10}}), iconTransformation(extent={{-90,-10},
                  {-70,10}})));
        Connectors.FaceBus yPositive if n_spec > 0
          "Positive face along the y axis" annotation (Placement(transformation(
                extent={{-10,30},{10,50}}), iconTransformation(extent={{-10,90},
                  {10,110}})));
        Connectors.FaceBus yNegative if n_spec > 0
          "Negative face along the y axis" annotation (Placement(transformation(
                extent={{-10,-50},{10,-30}}), iconTransformation(extent={{-10,-94},
                  {10,-74}})));
        Connectors.FaceBus zPositive if n_spec > 0
          "Positive face along the z axis" annotation (Placement(transformation(
                extent={{-30,-30},{-10,-10}}), iconTransformation(extent={{-90,
                  -90},{-70,-70}})));
        Connectors.FaceBus zNegative if n_spec > 0
          "Negative face along the z axis" annotation (Placement(transformation(
                extent={{10,10},{30,30}}), iconTransformation(extent={{40,40},{
                  60,60}})));

        outer parameter Q.Length L[Axis] if n_spec > 0 "Length" annotation (
            HideResult=true,missingInnerMessage=
              "This model should be used within a subregion model.");
        outer parameter Q.Area A[Axis] if n_spec > 0 "Cross-sectional area"
          annotation (HideResult=true,missingInnerMessage=
              "This model should be used within a subregion model.");
        // Note:  These must be public in Dymola 7.4, so HideResult is set true
        // instead.

      protected
        parameter Integer n_spec(start=0) "Number of species";
        final inner parameter Q.Length Lprime[Axis]=k .* A ./ L if n_spec > 0
          "Effective cross-sectional area per length";
        outer parameter Integer n_trans
          "Number of components of translational momentum" annotation (
            missingInnerMessage=
              "This model should be used within a subregion model.");
        outer parameter Integer cartAxes[:]
          "Cartesian-axis indices of the components of translational momentum"
          annotation (missingInnerMessage=
              "This model should be used within a subregion model.");
        outer parameter Boolean inclTrans[Axis]
          "true, if each component of translational momentum is included"
          annotation (missingInnerMessage=
              "This model should be used within a subregion model.");

        outer Conditions.Environment environment "Environmental conditions";
        PhaseBoundary phaseBoundary if n_spec > 0 "Phase boundary" annotation (
            Placement(transformation(
              extent={{-18,-18},{18,18}},
              rotation=0,
              origin={0,0})));
        // This component is conditional to prevent a mathematical singularity
        // when two or more empty phases (without any species included) are
        // connected.

        Connectors.InertInternal inert(
          n_trans=n_trans,
          final inclTranslational=reduceVel,
          final inclThermal=reduceTemp,
          translational(phi(
              each stateSelect=StateSelect.prefer,
              final start=phi_IC[cartAxes],
              final fixed={initVelX,initVelY,initVelZ}[index(inclTrans)])),
          thermal(T(
              stateSelect=StateSelect.prefer,
              final start=T_IC,
              final fixed=initTemp))) if n_spec > 0 and (reduceVel or
          reduceTemp)
          "Internal connector to directly couple velocities and/or temperatures"
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
                origin={26.67,-13.33}), iconTransformation(extent={{-10,-10},{
                  10,10}}, origin={26,-14})));
        // Note:  It would be simpler to use {initVelX, initVelY, initVelZ}[cartAxes]
        // for the fixed attribute of inert.translational, but Dymola 7.4 refuses to
        // accept it.
      equation
        // Inert interactions
        connect(phaseBoundary.inertAmagat, inertAmagat) annotation (Line(
            points={{6.2,-12.2},{13.33,-26.67}},
            color={72,90,180},
            smooth=Smooth.None));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="phase",
          Documentation(info="<html><p>If one of the species has <code>consEnergy = Conservation.IC</code>, then
    <code>initTemp</code> should be set to <code>false</code>.
    Likewise, if one of the species has <code>consTransX = Conservation.IC</code>,
    <code>consTransY = Conservation.IC</code>, or <code>consTransZ = Conservation.IC</code>, then
    <code>initVelX</code>, <code>initVelY</code>, or <code>initVelZ</code> should
    be set to <code>false</code> (respectively).</p>

    <p>The area fill factor (<b><i>k</i></b>) is a vector which inversely scales all
    the transport coefficients (&beta;, &zeta;, &eta;, and &theta;) of all of the species
    within the phase.  It can be used to introduce minor head loss or the effects of
    porosity or tortousity.  These effects may be anisotropic.</p>

    <p>Porosity is often quoted in material data sheets (e.g.,
    [<a href=\"modelica://FCSys.UsersGuide.References\">Toray2010</a>]) as volumetric porosity.  Using the
    Bruggeman correction factor [<a href=\"modelica://FCSys.UsersGuide.References\">Weber2004</a>, p. 4696],
    the area fill factor for the solid should be set to (1 - &epsilon;)<sup>3/2</sup>
    along each axis, where &epsilon; is the volumetric porosity (or volumetric fill factor
    of the gas).<sup><a href=\"#fn1\" id=\"ref1\">1</a></sup></p>
    
    <hr>

    <small>
    <p id=\"fn1\">1. Note that the Bruggeman correction contradicts what one would
    expect based on geometry&mdash;that the area fill factor would be the volumetric fill factor (1 - &epsilon;)
    raised to the two-thirds power (not three halfs).<a href=\"#ref1\" title=\"Jump back to footnote 1 in the text.\">&#8629;</a></p>

</html>"),
          Icon(graphics={
              Ellipse(
                extent={{-40,100},{40,20}},
                lineColor={127,127,127},
                startAngle=30,
                endAngle=149,
                pattern=LinePattern.Dash,
                fillPattern=FillPattern.Solid,
                fillColor={225,225,225}),
              Ellipse(
                extent={{20,-4},{100,-84}},
                lineColor={127,127,127},
                startAngle=270,
                endAngle=390,
                pattern=LinePattern.Dash,
                fillPattern=FillPattern.Solid,
                fillColor={225,225,225}),
              Ellipse(
                extent={{-100,-4},{-20,-84}},
                lineColor={127,127,127},
                startAngle=149,
                endAngle=270,
                pattern=LinePattern.Dash,
                fillPattern=FillPattern.Solid,
                fillColor={225,225,225}),
              Polygon(
                points={{60,-84},{-60,-84},{-94.5,-24},{-34.5,80},{34.5,80},{
                    94.5,-24},{60,-84}},
                pattern=LinePattern.None,
                fillPattern=FillPattern.Sphere,
                smooth=Smooth.None,
                fillColor={225,225,225},
                lineColor={0,0,0}),
              Line(
                points={{-60,-84},{60,-84}},
                color={127,127,127},
                pattern=LinePattern.Dash,
                smooth=Smooth.None),
              Line(
                points={{34.5,80},{94.5,-24}},
                color={127,127,127},
                pattern=LinePattern.Dash,
                smooth=Smooth.None),
              Line(
                points={{-34.5,80},{-94.5,-24}},
                color={127,127,127},
                pattern=LinePattern.Dash,
                smooth=Smooth.None),
              Text(
                extent={{-100,-20},{100,20}},
                textString="%name",
                lineColor={0,0,0})}),
          Diagram(graphics));
      end NullPhase;

    end BaseClasses;

  end Phases;

  package Species
    "Models for single-species storage, transport, and exchange of material, translational momentum, and energy"
    extends Modelica.Icons.Package;
    package 'C+' "C"
      extends Modelica.Icons.Package;
      package Graphite "<html>C<sup>+</sup> graphite</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Correlations with adjustment factors"
          extends SpeciesSolid(redeclare replaceable package Data =
                Characteristics.'C+'.Graphite, theta=k_theta*Data.theta(T, v));
          // Note:  In Dymola 7.4,
          // "redeclare replaceable package Data = FCSys.Characteristics.C.Graphite"
          // must be used instead of
          // "redeclare replaceable FCSys.Characteristics.C.Graphite Data" in
          // order for this model to pass its check.  This applies to the other
          // species models too.

          parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C",
            Documentation(info="<html><p>See the information in the
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
          extends SpeciesSolid(redeclare replaceable package Data =
                Characteristics.'C+'.Graphite);
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C",
            Documentation(info=
                  "<html><p>See the information in the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"),

            Icon(graphics));

        end Correlated;

        model Fixed "Fixed properties"
          import Modelica.Math.log;
          import FCSys.BaseClasses.Utilities.Polynomial;

          extends SpeciesSolid(redeclare replaceable package Data =
                Characteristics.'C+'.Graphite,redeclare parameter
              Q.ResistivityThermal theta=Data.theta());
          /* **
(
        Deltah0_f=0,
        Deltah0=0,
        n_c=0,
        T_lim_c={0,Modelica.Constants.inf,Modelica.Constants.inf},
        b_c=[935*U.J*Data.m/(U.kg*U.K);0],
        B_c=[0, 0;0,0])
        B_c=[-298.15*U.K*935*U.J*Data.m/(U.kg*U.K) + Data.Deltah0_f,
            Polynomial.F(
            298.15*U.K,
            Characteristics.'C+'.Graphite.b_c[1, :],
            -3) + FCSys.Characteristics.'C+'.Graphite.B_c[1, 2] - Data.b_c[1, 1]
            *lnog298.15*U.K)]), redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(5.70
        *U.W)*/
          // See the documentation for a table of values.
          // Note:  Parameter expressions (e.g., involving environment.T) are not
          // used here since they would render the parameters unadjustable in Dymola
          // 7.4.  A similar note applies to the other species.
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Fixed specific heat capacity (independent of temperature)
    <li>Thermal resistivity &theta; is fixed (e.g., independent of temperature)</li>
    </ol></p>

    <p>The default isobaric specific heat capacity (<code>b_c=[0, 935*U.J*Data.m/(U.kg*U.K)]</code>) and thermal
   resistivity (<code>theta=U.m*U.K/(11.1*U.W)</code>) is based on data of graphite fiber epoxy (25% vol)<br>composite at 300 K from
   Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 909].
   Related data is listed in <a href=\"#Tab1\">Table 1</a>.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"top\" id=\"Tab1\">Table 1: Properties of forms of C [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 909].</caption>
    <tr>
      <th rowspan=3 valign=\"middle\"><code>T<br>/U.K</code></th>
      <th rowspan=1 colspan=2 width=1 valign=\"middle\">Diamond (type IIa)</th>
      <th rowspan=1 colspan=1 width=1 valign=\"middle\">Amorphous<br>carbon</th>
      <th rowspan=1 colspan=3 width=1 valign=\"middle\">Graphite (pyrolytic)</th>
      <th rowspan=1 colspan=3 width=1>Graphite fiber epoxy (25% vol)<br>composite</th>
    </tr>
    <tr>
      <th rowspan=2 valign=\"middle\"><code>c_p*U.kg<br>*U.K<br>/(U.J<br>*Data.m)</code></th>
      <th rowspan=2 valign=\"middle\"><code>theta<br>*U.W<br>/(U.m<br>*U.K)</code></th>
      <th rowspan=2 valign=\"middle\"><code>theta<br>*U.W<br>/(U.m<br>*U.K)</code></th>
      <th rowspan=2 valign=\"middle\"><code>c_p*U.kg<br>*U.W<br>/(U.J<br>*Data.m)</code></th>
      <th rowspan=1 colspan=2><code>theta*U.W/(U.m*U.K)</code></th>
      <th rowspan=2 valign=\"middle\"><code>c_p*U.kg<br>*U.K<br>/(U.J<br>*Data.m)</code></th>
      <th rowspan=1 colspan=2><code>theta*U.W/(U.m*U.K)</code></th>
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

    end 'C+';

    package 'C19HF37O5S-'
      "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup></html>"
      extends Modelica.Icons.Package;
      package Ionomer
        "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S ionomer</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Correlations with adjustment factors"
          extends SpeciesSolid(redeclare replaceable package Data =
                Characteristics.'C19HF37O5S-'.Ionomer, theta=k_theta*Data.theta(
                T, v));

          parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C19HF37O5S",
            Documentation(info=
                  "<html><p>See the information in the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"),

            Diagram(graphics));

        end Calibrated;

        model Correlated "Correlated properties"
          extends SpeciesSolid(redeclare replaceable package Data =
                Characteristics.'C19HF37O5S-'.Ionomer);
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C19HF37O5S",
            Documentation(info=
                  "<html><p>See the information in the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Correlated;

        model Fixed "Fixed properties"
          extends SpeciesSolid(redeclare replaceable package Data =
                Characteristics.'C19HF37O5S-'.Ionomer, redeclare parameter
              Q.ResistivityThermal theta=Data.theta());
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="C19HF37O5S",
            Documentation(info="<html><p>Assumptions:
    <ol>
    <li>Thermal resistivity &theta; is fixed (e.g., independent of temperature)</li>
    </ol></p>

    <p>Notes:
    <ul><li>The default thermal resistivity (<code>theta=U.m*U.K/(0.16*U.W)</code>) is of dry
  Nafion 115 [<a href=\"modelica://FCSys.UsersGuide.References\">Kandlikar2009</a>, p. 1277].</li>
  </ul></p>

<p>For more information, see the
    <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Fixed;

      end Ionomer;

    end 'C19HF37O5S-';

    package 'e-' "<html>e<sup>-</sup></html>"
      extends Modelica.Icons.Package;
      package Graphite "<html>e<sup>-</sup> in graphite</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Correlations with adjustment factors"
          extends Species(
            redeclare replaceable package Data = Characteristics.'e-'.Graphite,

            beta=k_beta*Data.beta(T, v),
            zeta=k_zeta*Data.zeta(T, v),
            theta=k_theta*Data.theta(T, v));

          parameter Q.NumberAbsolute k_beta(final nominal=1) = 1
            "<html>Adjustment factor for dynamic compressibility (<i>k</i><sub>&beta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_zeta(final nominal=1) = 1
            "<html>Adjustment factor for fluidity (<i>k</i><sub>&zeta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'e-'",
            Documentation(info=
                  "<html><p>See the information in the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));
          // **Add parameters k_mu, k_nu, k_eta for this and all calibrated species.
        end Calibrated;

        model Correlated "Correlated properties"
          extends Species(redeclare replaceable package Data =
                Characteristics.'e-'.Graphite);

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'e-'",
            Documentation(info=
                  "<html><p>See the information in the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare replaceable package Data = Characteristics.'e-'.Graphite,

            redeclare parameter Q.Fluidity beta=Data.beta(),
            redeclare parameter Q.Fluidity zeta=Data.zeta(),
            redeclare parameter Q.ResistivityThermal theta=Data.theta());

          annotation (
            group="Material properties",
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'e-'",
            Documentation(info="<html>
    <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"),

            Diagram(graphics));

          // **set theta=0 (final), zeta=0 (final), use small Lprime factor to reduce thermal and translational coupling with solid

        end Fixed;

      end Graphite;

    end 'e-';

    package 'H+' "<html>H<sup>+</sup></html>"
      extends Modelica.Icons.Package;
      package Ionomer "<html>H<sup>+</sup> in ionomer</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Correlations with adjustment factors"
          extends Species(
            redeclare replaceable package Data = Characteristics.'H+'.Ionomer (
                  n_v=FCSys.Characteristics.'C19HF37O5S-'.Ionomer.n_v, b_v=
                    FCSys.Characteristics.'C19HF37O5S-'.Ionomer.b_v),
            beta=k_beta*Data.beta(T, v),
            zeta=k_zeta*Data.zeta(T, v),
            theta=k_theta*Data.theta(T, v));

          parameter Q.NumberAbsolute k_beta(final nominal=1) = 1
            "<html>Adjustment factor for dynamic compressibility (<i>k</i><sub>&beta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_zeta(final nominal=1) = 1
            "<html>Adjustment factor for fluidity (<i>k</i><sub>&zeta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</i></sub>)</html>"
            annotation (Dialog(group="Material properties"));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'H+'",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>The density of H<sup>+</sup> is equal to that of
  C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup> or approximately 1.912 M.  Note that
  this is greater than that measured by Spry and Fayer (see
  <a href=\"modelica://FCSys.Characteristics.'H+'.Ionomer\">Characteristics.'H+'.Ionomer</a>), but it
  simplifies the model by requiring only C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>
  (not C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S) for charge neutrality.
    </li>
    </ol></p>

  <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Calibrated;

        model Correlated "Correlated properties"
          extends Species(redeclare replaceable package Data =
                Characteristics.'H+'.Ionomer (n_v=FCSys.Characteristics.
                    'C19HF37O5S-'.Ionomer.n_v, b_v=FCSys.Characteristics.
                    'C19HF37O5S-'.Ionomer.b_v));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'H+'",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>The density of H<sup>+</sup> is equal to that of
  C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup> or approximately 1.912 M.  Note that
  this is greater than that measured by Spry and Fayer (see
  <a href=\"modelica://FCSys.Characteristics.'H+'.Ionomer\">Characteristics.'H+'.Ionomer</a>), but it
  simplifies the model by requiring only C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>
  (not C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S) for charge neutrality.
    </li>
    </ol></p>

    <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare replaceable package Data = Characteristics.'H+'.Ionomer,
            initMaterial=InitScalar.None,
            redeclare parameter Q.Fluidity beta=Data.beta(),
            redeclare parameter Q.Fluidity zeta=Data.zeta(),
            redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(0.1661*U.W));

          /*
    *Are the trivial Data modifications necessary?
    (
          n_v=FCSys.Characteristics.'C19HF37O5S-'.Ionomer.n_v,
          b_v=FCSys.Characteristics.'C19HF37O5S-'.Ionomer.b_v)*/
          // **temp initmeth (if keep it, copy to other models)
          // *
          // See the documentation for a table of values.
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'H+'",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>The generalized resistivities (&beta;, &zeta;, &theta;) are fixed (e.g., independent of temperature).</li>
    <li>The density of H<sup>+</sup> is equal to that of
  C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup> or approximately 1.912 M.  Note that
  this is greater than that measured by Spry and Fayer (see
  <a href=\"modelica://FCSys.Characteristics.'H+'.Ionomer\">Characteristics.'H+'.Ionomer</a>), but it
  simplifies the model by requiring only C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>
  (not C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S) for charge neutrality.
    </li>
    </ol></p>

  <p>The default resistivities (<code>zeta=1/(5.3e-6*U.Pa*U.s)</code> and <code>theta=U.m*U.K/(0.1661*U.W)</code>) are of H gas
  (rather than H<sup>+</sup>) at 300 K from [<a href=\"modelica://FCSys.UsersGuide.References\">Schetz1996</a>, p. 139].
  <a href=\"#Tab1\">Table 1</a> lists the properties at other temperatures.</p>

    <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"top\" id=\"Tab1\">Table 1: Properties of H gas (not H<sup>+</sup>) [<a href=\"modelica://FCSys.UsersGuide.References\">Schetz1996</a>, p. 139]</caption>
<tr>
      <th valign=\"middle\"><code>T<br>/U.K</code></th>
      <th width=1><code>zeta<br>*U.Pa*U.s</code></th>
      <th width=1><code>theta*U.W<br>/(U.m*U.K)</code></th>
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
  </table></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Fixed;

      end Ionomer;

    end 'H+';

    package H2 "<html>H<sub>2</sub></html>"
      extends Modelica.Icons.Package;
      package Gas "<html>H<sub>2</sub> gas</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Correlations with adjustment factors"
          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.H2.Gas (
                  b_v=[1],n_v={-1,0}),
            beta=k_beta*Data.beta(T, v),
            zeta=k_zeta*Data.zeta(T, v),
            theta=k_theta*Data.theta(T, v));

          parameter Q.NumberAbsolute k_beta(final nominal=1) = 1
            "<html>Adjustment factor for dynamic compressibility (<i>k</i><sub>&beta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_zeta(final nominal=1) = 1
            "<html>Adjustment factor for fluidity (<i>k</i><sub>&zeta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));

        end Calibrated;

        model Correlated "Correlated properties"
          extends Species(redeclare replaceable package Data =
                FCSys.Characteristics.H2.Gas (b_v=[1], n_v={-1,0}));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));

        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.H2.Gas (
                  b_v=[1],n_v={-1,0}),
            redeclare parameter Q.Fluidity beta=Data.beta(),
            redeclare parameter Q.Fluidity zeta=1/(89.6e-7*U.Pa*U.s),
            redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(183e-3*U.W));
          // See the documentation for a table of values.
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
        <li>The generalized resistivities (&beta;, &zeta;, &theta;) are fixed (e.g., independent of temperature).</li>

    </ol></p>

<p>Additional notes:<ul>
<li>
  The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

<p>The default resistivities (<code>zeta=1/(89.6e-7*U.Pa*U.s)</code>
and <code>theta=U.m*U.K/(183e-3*U.W)</code>) are based on data of H<sub>2</sub> gas at 1 atm and
  300 K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 919&ndash;920].
  <a href=\"#Tab1\">Table 1</a> lists the properties at other temperatures.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"top\" id=\"Tab1\">Table 1: Properties of H<sub>2</sub> gas at 1 atm [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 919&ndash;920].</caption>
  <tr>
      <th valign=\"middle\"><code>T<br>/U.K</code></th>
      <th width=1><code>c_p*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1 ><code>zeta<br>*U.Pa*U.s</code></th>
      <th width=1 ><code>theta*U.W<br>/(U.m*U.K)</code></th>
    </tr>
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
                (b_v=[1], n_v={-1,0}),
            beta=k_beta*Data.beta(T, v),
            zeta=k_zeta*Data.zeta(T, v),
            theta=k_theta*Data.theta(T, v));

          parameter Q.NumberAbsolute k_beta(final nominal=1) = 1
            "<html>Adjustment factor for dynamic compressibility (<i>k</i><sub>&beta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_zeta(final nominal=1) = 1
            "<html>Adjustment factor for fluidity (<i>k</i><sub>&zeta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));

        end Calibrated;

        model Correlated "Correlated properties"
          extends Species(redeclare replaceable package Data =
                FCSys.Characteristics.H2O.Gas (b_v=[1], n_v={-1,0}));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));

        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.H2O.Gas
                (b_v=[1], n_v={-1,0}),
            redeclare parameter Q.Fluidity beta=Data.beta(),
            redeclare parameter Q.Fluidity zeta=1/(9.09e-6*U.Pa*U.s),
            redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(19.6e-3*U.W));

          // See the documentation for tables of values.
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
        <li>The generalized resistivities (&beta;, &zeta;, &theta;) are fixed (e.g., independent of temperature).</li>
    </ol></p>

  <p>Notes:<ul>
  <li>The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

<p>The default resistivities (<code>zeta=1/(9.09e-6*U.Pa*U.s)</code>
and <code>theta=U.m*U.K/(19.6e-3*U.W)</code>) are of H<sub>2</sub>O gas at saturation pressure and
  300 K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 921].  <a href=\"#Tab1\">Table 1</a> lists the properties at
  saturation pressure and other temperatures.  <a href=\"#Tab2\">Table 2</a> lists the properties of H<sub>2</sub>O gas at 1 atm.
  See also
  <a href=\"http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html\">http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html</a>.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"top\" id=\"Tab1\">Table 1: Properties of H<sub>2</sub>O gas at saturation pressure [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 924&ndash;925].</caption>
  <tr>
      <th valign=\"middle\"><code>T<br>/U.K</code></th>
      <th width=1><code>c_p*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1><code>zeta<br>*U.Pa*U.s</code></th>
      <th width=1><code>theta*U.W<br>/(U.m*U.K)</code></th>
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
<tr><td>440</td><td>2460</td><td>1/14.50e-6</td><td>1/31.7e-3</td></tr>
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
    <caption align=\"top\" id=\"Tab2\">Table 2: Properties of H<sub>2</sub>O gas at 1 atm [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 921].</caption>
  <tr>
      <th valign=\"middle\"><code>T<br>/U.K</code></th>
      <th width=1><code>c_p*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1 ><code>zeta<br>*U.Pa*U.s</code></th>
      <th width=1 ><code>theta*U.W<br>/(U.m*U.K)</code></th>
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
<br></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));

        end Fixed;

      end Gas;

      package Ionomer "<html>H<sub>2</sub>O in ionomer</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Correlations with adjustment factors"
          extends Species(
            redeclare replaceable package Data = Characteristics.H2O.Ionomer (
                  b_v=[1], n_v={-1,0}),
            beta=k_beta*Data.beta(T, v),
            zeta=k_zeta*Data.zeta(T, v),
            theta=k_theta*Data.theta(T, v));

          parameter Q.NumberAbsolute k_beta(final nominal=1) = 1
            "<html>Adjustment factor for dynamic compressibility (<i>k</i><sub>&beta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_zeta(final nominal=1) = 1
            "<html>Adjustment factor for fluidity (<i>k</i><sub>&zeta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));

        end Calibrated;

        model Correlated "Correlated properties"
          extends Species(redeclare replaceable package Data =
                Characteristics.H2O.Ionomer (b_v=[1], n_v={-1,0}));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));

        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare replaceable package Data = Characteristics.H2O.Ionomer (
                  b_v=[1], n_v={-1,0}),
            redeclare parameter Q.Fluidity beta=Data.beta(),
            redeclare parameter Q.Fluidity zeta=Data.zeta(),
            redeclare parameter Q.ResistivityThermal theta=Data.theta());
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
        <li>The generalized resistivities (&beta;, &zeta;, &theta;) are fixed (e.g., independent of temperature).</li>
    </ol></p></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));

        end Fixed;

      end Ionomer;

      package Liquid "<html>H<sub>2</sub>O liquid</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Correlations with adjustment factors"
          extends SpeciesIncompressible(
            redeclare replaceable package Data = Characteristics.H2O.Liquid,
            beta=k_beta*Data.beta(T, v),
            zeta=k_zeta*Data.zeta(T, v),
            theta=k_theta*Data.theta(T, v));

          parameter Q.NumberAbsolute k_beta(final nominal=1) = 1
            "<html>Adjustment factor for dynamic compressibility (<i>k</i><sub>&beta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_zeta(final nominal=1) = 1
            "<html>Adjustment factor for fluidity (<i>k</i><sub>&zeta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html>
    <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));

        end Calibrated;

        model Correlated "Correlated properties"
          extends SpeciesIncompressible(redeclare replaceable package Data =
                Characteristics.H2O.Liquid);
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html>
         <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));

        end Correlated;

        model Fixed "Fixed properties"
          extends SpeciesIncompressible(
            redeclare replaceable package Data = Characteristics.H2O.Liquid,
            redeclare parameter Q.Fluidity beta=Data.beta(),
            redeclare parameter Q.Fluidity zeta=1/(855e-6*U.Pa*U.s),
            redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(613e-3*U.W));

          // See the documentation for tables of values.
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
        <li>The generalized resistivities (&beta;, &zeta;, &theta;) are fixed (e.g., independent of temperature).</li>

    </ol></p>

          <p>Notes:<ul>
  <li>
  The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

<p>The default resistivities (<code>zeta=1/(855e-6*U.Pa*U.s)</code>
and <code>theta=U.m*U.K/(613e-3*U.W)</code>) are of H<sub>2</sub>O liquid at saturation pressure and
  300 K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 921].  <a href=\"#Tab1\">Table 1</a> lists the properties at
  saturation pressure and other temperatures.
  See also
  <a href=\"http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html\">http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html</a>.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"top\" id=\"Tab1\">Table 1: Properties of H<sub>2</sub>O liquid at saturation pressure [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 924&ndash;925].</caption>
  <tr>
      <th valign=\"middle\"><code>T<br>/U.K</code></th>
      <th width=1><code>c_p*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1><code>zeta<br>*U.Pa*U.s</code></th>
      <th width=1><code>theta*U.W<br>/(U.m*U.K)</code></th>
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
                  b_v=[1],n_v={-1,0}),
            beta=k_beta*Data.beta(T, v),
            zeta=k_zeta*Data.zeta(T, v),
            theta=k_theta*Data.theta(T, v));

          parameter Q.NumberAbsolute k_beta(final nominal=1) = 1
            "<html>Adjustment factor for dynamic compressibility (<i>k</i><sub>&beta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_zeta(final nominal=1) = 1
            "<html>Adjustment factor for fluidity (<i>k</i><sub>&zeta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="N2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));

        end Calibrated;

        model Correlated "Correlated properties"
          extends Species(redeclare replaceable package Data =
                FCSys.Characteristics.N2.Gas (b_v=[1], n_v={-1,0}));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="N2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));

        end Correlated;

        model Fixed "Fixed properties"
          import FCSys.BaseClasses.Utilities.Polynomial;
          import Modelica.Math.log;

          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.N2.Gas (
                b_v=[1],
                n_v={-1,0},
                n_c=0,
                T_lim_c={0,Modelica.Constants.inf},
                b_c=[1.041e3*U.J*Data.m/(U.kg*U.K)],
                B_c=[-298.15*U.K*1.041e3*U.J*Data.m/(U.kg*U.K) + Data.Deltah0_f,
                    Polynomial.F(
                            298.15*U.K,
                            Characteristics.N2.Gas.b_c[1, :],
                            -3) + FCSys.Characteristics.N2.Gas.B_c[1, 2] - Data.b_c[
                    1, 1]*log(298.15*U.K)]),
            redeclare parameter Q.Fluidity beta=Data.beta(),
            redeclare parameter Q.Fluidity zeta=1/(178.2e-7*U.Pa*U.s),
            redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(25.9e-3*U.W));
          // **Add table from [Present1958, p. 263] to the documentation (see Tests.Characteristics.N2.eta).

          // See the documentation for a table of values.
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="N2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
    <li>Fixed specific heat capacity (independent of temperature)</ol>
        <li>The generalized resistivities (&beta;, &zeta;, &theta;) are fixed (e.g., independent of temperature).</li>

    </ol></p>

<p>The default specific heat capacity (<code>b_c=[1.041e3*U.J*Data.m/(U.kg*U.K)]</code>) and resistivities (<code>zeta=1/(178.2e-7*U.Pa*U.s)</code> and <code>theta=U.m*U.K/(25.9e-3*U.W))</code>) are based on data of gas at 1 atm and
  300 K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 920].
  <a href=\"#Tab1\">Table 1</a> lists the properties at other temperatures. Note that the value for isobaric specific heat capacity at
  800 K (<code>c_p=1.22e3*U.J*Data.m/(U.kg*U.K)</code>) seems unusual, but it matches the
  reference.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"top\" id=\"Tab1\">Table 1: Properties of N<sub>2</sub> gas at 1 atm [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 920]</caption>
  <tr>
      <th valign=\"middle\"><code>T<br>/U.K</code></th>
      <th width=1><code>c_p*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1><code>zeta<br>*U.Pa*U.s</code></th>
      <th width=1><code>theta*U.W<br>/(U.m*U.K)</code></th>
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
       <code>zeta=1/(178e-7*U.Pa*U.s)</code>
   (<a href=\"http://en.wikipedia.org/wiki/Viscosity\">http://en.wikipedia.org/wiki/Viscosity</a>).</p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));

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
                  b_v=[1],n_v={-1,0}),
            beta=k_beta*Data.beta(T, v),
            zeta=k_zeta*Data.zeta(T, v),
            theta=k_theta*Data.theta(T, v));

          parameter Q.NumberAbsolute k_beta(final nominal=1) = 1
            "<html>Adjustment factor for dynamic compressibility (<i>k</i><sub>&beta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_zeta(final nominal=1) = 1
            "<html>Adjustment factor for fluidity (<i>k</i><sub>&zeta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="O2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));

        end Calibrated;

        model Correlated "Correlated properties"
          extends Species(redeclare replaceable package Data =
                FCSys.Characteristics.O2.Gas (b_v=[1], n_v={-1,0}));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="O2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"));

        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.O2.Gas (
                  b_v=[1],n_v={-1,0}),
            redeclare parameter Q.Fluidity beta=Data.beta(),
            redeclare parameter Q.Fluidity zeta=1/(207.2e-7*U.Pa*U.s),
            redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(26.8e-3*U.W));
          // **Add table from Present1958 p. 263 to the documentation (see Tests.Characteristics.O2.eta).

          // See the documentation for a table of values.
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="O2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
        <li>The generalized resistivities (&beta;, &zeta;, &theta;) are fixed (e.g., independent of temperature).</li>

    </ol></p>

<p>Additional notes:
<ul>
          <li>
  The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

  <p>The default resistivities (<code>zeta=1/(207.2e-7*U.Pa*U.s)</code> and <code>theta=U.m*U.K/(26.8e-3*U.W)</code>) are based on data of gas at 1 atm and
  300 K from Incropera and DeWitt  [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 920&ndash;921].
  <a href=\"#Tab1\">Table 1</a> lists the properties at other temperatures.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"top\" id=\"Tab1\">Table 1: Properties of O<sub>2</sub> gas at 1 atm [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 920&ndash;921]</caption>
  <tr>
      <th valign=\"middle\"><code>T<br>/U.K</code></th>
      <th width=1><code>c_p*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1><code>zeta<br>*U.Pa*U.s</code></th>
      <th width=1><code>theta*U.W<br>/(U.m*U.K)</code></th>
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
  </table></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Specues,Species\">Species</a> model.</p></html>"),

            Icon(graphics));

        end Fixed;

      end Gas;

    end O2;

    model SpeciesSolid
      "<html><a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model for a solid (inert and stagnant)</html>"
      extends SpeciesIncompressible(
        final upstreamX=false,
        final upstreamY=false,
        final upstreamZ=false,
        final eta=0,
        final beta=0,
        final zeta=0,
        final phi_IC=zeros(3),
        final I_IC,
        invertEOS=false);
      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="species",
        Documentation(info="<html><p>Assumptions:<ol>
  <li>Zero dynamic compressibility (&rArr; uniform current)</li>
  <li>Zero fluidity (&rArr; no shearing)</li></ol></p>

  <p>Usually, conditions should be applied to specify the velocity (typically <b>0</b>).  In a group of connected solid species
  of a single type (instances of a model derived from this one), there should be exactly one equation to specify the velocity along each Cartesian axis.
  For example, the x-axis velocity may be given by setting <code>setVelX</code> to <code>true</code> in one of the instances
  (x-axis velocity will be <code>phi_IC[Axis.x]</code> for all time).  Alternatively, one of the faces on the outside of the
  group could be removed by setting the appropriate <code>inclFaceNegX</code>, <code>inclFacePosX</code>, etc. parameter
  to <code>false</code>, which will set the current and transverse components of velocity to zero.</p>

  <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

    end SpeciesSolid;

    model SpeciesIncompressible
      "<html><a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model for an incompressible fluid</html>"
      extends Species(initMaterial=InitScalar.Volume);
      // Note:  The default, pressure, can't be used to initialize an
      // incompressible species.
      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="species",
        Documentation(info="<html>
  <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

    end SpeciesIncompressible;

    model Species
      "Model to exchange, transport, and store the material, momentum, and energy of one species"
      import Modelica.Math.log10;
      import FCSys.BaseClasses.Utilities.cartWrap;
      import FCSys.BaseClasses.Utilities.inSign;
      import FCSys.BaseClasses.Utilities.Delta;
      import FCSys.BaseClasses.Utilities.Sigma;
      import FCSys;
      extends FCSys.BaseClasses.Icons.Names.Top4;

      // Material properties
      replaceable package Data = Characteristics.BaseClasses.Characteristic
        constrainedby Characteristics.BaseClasses.Characteristic
        "Characteristic data" annotation (
        Dialog(group="Material properties"),
        choicesAllMatching=true,
        __Dymola_choicesFromPackage=true,
        Placement(transformation(extent={{-60,40},{-40,60}}),
            iconTransformation(extent={{-10,90},{10,110}})));
      parameter Boolean invertEOS=true "Invert the equation of state"
        annotation (Dialog(Evaluate=true, group="Material properties"), choices(
            __Dymola_checkBox=true));
      Q.Mobility mu(nominal=0.1*U.C*U.s/U.kg) = Data.mu(T, v)
        "<html>Mobility (&mu;)</html>"
        annotation (Dialog(group="Material properties"));
      Q.TimeAbsolute nu(nominal=1e-9*U.s) = Data.nu(T, v)
        "<html>Thermal independity (&nu;)</html>"
        annotation (Dialog(group="Material properties"));
      Q.ResistivityMaterial eta(nominal=10e-6*U.s/U.m^2) = Data.eta(T, v)
        "<html>Material resistivity (&eta;)</html>"
        annotation (Dialog(group="Material properties"));
      Q.Fluidity beta(nominal=10*U.cm*U.s/U.g) = Data.beta(T, v)
        "<html>Dynamic compressibility (&beta;)</html>"
        annotation (Dialog(group="Material properties"));
      Q.Fluidity zeta(nominal=10*U.cm*U.s/U.g) = Data.zeta(T, v)
        "<html>Fluidity (&zeta;)</html>"
        annotation (Dialog(group="Material properties"));
      Q.ResistivityThermal theta(nominal=10*U.cm/U.A) = Data.theta(T, v)
        "<html>Thermal resistivity (&theta;)</html>"
        annotation (Dialog(group="Material properties"));

      // Reactions
      parameter Integer n_react(min=0) = 1
        "<html>Number of interactions (<i>n</i><sub>react</sub>)</html>"
        annotation (Dialog(group="Reactions and phase change"), HideResult=true);
      Q.CurrentAbsolute Ndot_0[n_react](each nominal=1e-3*U.A) = fill(1e-6*U.A,
        n_react) "<html>Exchange currents (<i>N&#775;</i><sub>0</sub>)</html>"
        annotation (Dialog(group="Reactions and phase change"));

      // Assumptions
      // -----------
      // Dynamics
      parameter Conservation consMaterial=Conservation.Dynamic "Material"
        annotation (Dialog(
          Evaluate=true,
          tab="Assumptions",
          group="Formulation of conservation equations"));
      parameter Conservation consTransX=Conservation.Dynamic
        "X-axis translational momentum" annotation (Evaluate=true, Dialog(
          tab="Assumptions",
          group="Formulation of conservation equations",
          enable=inclTrans[1]));
      parameter Conservation consTransY=Conservation.Dynamic
        "Y-axis translational momentum" annotation (Evaluate=true, Dialog(
          tab="Assumptions",
          group="Formulation of conservation equations",
          enable=inclTrans[2]));
      parameter Conservation consTransZ=Conservation.Dynamic
        "Z-axis translational momentum" annotation (Evaluate=true, Dialog(
          tab="Assumptions",
          group="Formulation of conservation equations",
          enable=inclTrans[3]));
      parameter Conservation consEnergy=Conservation.Dynamic "Energy"
        annotation (Evaluate=true, Dialog(tab="Assumptions", group=
              "Formulation of conservation equations"));
      // **If the static option isn't useful, remove it and go back
      // to Boolean setMaterial, setTransX, etc.
      //
      // Flow conditions
      Q.NumberAbsolute Nu_Phi[Axis]={4,4,4}
        "<html>Translational Nusselt numbers (<b><i>Nu</i><sub>&Phi;</sub></b>)</html>"
        annotation (Dialog(tab="Assumptions",group="Flow conditions"));
      Q.NumberAbsolute Nu_Q=3.66
        "<html>Thermal Nusselt number (<i>Nu</i><sub><i>Q</i></sub>)</html>"
        annotation (Dialog(tab="Assumptions",group="Flow conditions"));
      //
      // Upstream discretization
      parameter Boolean upstreamX=true "X" annotation (
        Evaluate=true,
        Dialog(
          tab="Assumptions",
          group="Axes with upstream discretization",
          enable=inclTrans[1],
          compact=true),
        choices(__Dymola_checkBox=true));
      parameter Boolean upstreamY=true "Y" annotation (
        Evaluate=true,
        Dialog(
          tab="Assumptions",
          group="Axes with upstream discretization",
          enable=inclTrans[2],
          compact=true),
        choices(__Dymola_checkBox=true));
      parameter Boolean upstreamZ=true "Z" annotation (
        Evaluate=true,
        Dialog(
          tab="Assumptions",
          group="Axes with upstream discretization",
          enable=inclTrans[3],
          compact=true),
        choices(__Dymola_checkBox=true));

      // Initialization parameters
      // -------------------------
      // Scalar properties
      parameter InitScalar initMaterial=InitScalar.Pressure
        "Method of initializing the material balance" annotation (Evaluate=true,
          Dialog(tab="Initialization", group="Material and energy"));
      parameter InitScalar initEnergy=InitScalar.Temperature
        "Method of initializing the energy balance" annotation (Evaluate=true,
          Dialog(tab="Initialization", group="Material and energy"));
      parameter Q.Amount N_IC(start=V_IC*rho_IC)
        "<html>Initial particle number (<i>N</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Material and energy"));
      // Note:  This parameter is left enabled even it isn't used to
      // explicitly initialize any states, since it's used as a guess value.
      // Similar notes apply to some other initial conditions below.
      parameter Q.Density rho_IC(start=1/Data.v_Tp(T_IC, p_IC))
        "<html>Initial density (&rho;<sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Material and energy"));
      parameter Q.Volume V_IC(start=product(L))
        "<html>Initial volume (<i>V</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Material and energy"));
      parameter Q.PressureAbsolute p_IC(start=environment.p)
        "<html>Initial pressure (<i>p</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Material and energy"));
      parameter Q.TemperatureAbsolute T_IC(start=environment.T)
        "<html>Initial temperature (<i>T</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Material and energy"));
      parameter Q.Potential h_IC(start=Data.h(T_IC, p_IC))
        "<html>Initial specific enthalpy (<i>h</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization", group="Material and energy"));
      parameter Q.Potential g_IC(start=Data.g(T_IC, p_IC))
        "<html>Initial Gibbs potential (<i>g</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization", group="Material and energy"));
      //
      // Velocity
      parameter InitTranslational initTransX=InitTranslational.Velocity
        "Method of initializing the x-axis balance" annotation (Evaluate=true,
          Dialog(
          tab="Initialization",
          group="Translational momentum",
          enable=inclTrans[1]));
      parameter InitTranslational initTransY=InitTranslational.Velocity
        "Method of initializing the y-axis balance" annotation (Evaluate=true,
          Dialog(
          tab="Initialization",
          group="Translational momentum",
          enable=inclTrans[2]));
      parameter InitTranslational initTransZ=InitTranslational.Velocity
        "Method of initializing the z-axis balance" annotation (Evaluate=true,
          Dialog(
          tab="Initialization",
          group="Translational momentum",
          enable=inclTrans[3]));
      // Note:  Dymola 7.4 doesn't provide pull-down lists for arrays of
      // enumerations; therefore, a parameter is used for each axis.
      parameter Q.Velocity phi_IC[Axis]={0,0,0}
        "<html>Initial velocity (<b>&phi;</b><sub>IC</sub>)</html>" annotation
        (Dialog(tab="Initialization", group="Translational momentum"));
      parameter Q.Current I_IC[Axis]={0,0,0}
        "<html>Initial current (<i><b>I</b></i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization",group="Translational momentum"));

      // Preferred states
      Q.Amount N(
        nominal=4*U.C,
        final start=N_IC,
        final fixed=false,
        stateSelect=StateSelect.prefer) "Particle number";
      // Note:  The start value for this variable (and others below) isn't fixed
      // because the related initial condition is applied in the initial
      // equation section.
      Q.Velocity phi[n_trans](
        each nominal=10*U.cm/U.s,
        final start=phi_IC[cartAxes],
        each final fixed=false,
        each stateSelect=StateSelect.prefer) "Velocity";
      Q.TemperatureAbsolute T(
        nominal=300*U.K,
        final start=T_IC,
        final fixed=false,
        stateSelect=StateSelect.prefer) "Temperature";

      // Aliases (for common terms)
      Q.PressureAbsolute p(
        nominal=U.atm,
        final start=p_IC,
        final fixed=false,
        stateSelect=StateSelect.never) "Pressure";
      // StateSelect.never is necessary to avoid dynamic state selection
      // in Dymola 7.4.
      Q.PressureAbsolute p_faces[Axis, Side](each nominal=U.atm, each final
          start=p_IC) "Pressures at the faces";
      Q.PotentialAbsolute sT(
        nominal=3000*U.K,
        final start=h_IC - g_IC,
        final fixed=false) "Product of specific entropy and temperature";
      Q.Mass M(nominal=1e-3*U.g, start=Data.m*N_IC) "Mass";
      Q.Volume V(
        nominal=U.cc,
        final start=V_IC,
        final fixed=false) "Volume";
      Q.VolumeSpecific v(
        nominal=U.cc/(4*U.C),
        final start=1/rho_IC,
        final fixed=false) "Specific volume";
      Q.Density rho(
        nominal=4*U.C/U.cc,
        final start=rho_IC,
        final fixed=false) "Density";
      Q.Potential h(
        nominal=U.V,
        final start=h_IC,
        final fixed=false,
        stateSelect=StateSelect.never) "Specific enthalpy";
      // StateSelect.never is necessary to avoid dynamic state selection
      // in Dymola 7.4.
      Q.Current I[n_trans](
        each nominal=U.A,
        final start=I_IC[cartAxes],
        each final fixed=false) "Current";

      // Auxiliary variables (for analysis)
      // ----------------------------------
      // Misc. properties and conditions
      /*
  output Q.Potential g(stateSelect=StateSelect.never) = chemical.mu if
    environment.analysis "Electrochemical potential";
  */
      /* **
  output Q.NumberAbsolute s(stateSelect=StateSelect.never) = Data.s(T, p) if
    environment.analysis "Specific entropy";
  output Q.Amount S(stateSelect=StateSelect.never) = N*s if environment.analysis
    "Entropy";
  output Q.PressureAbsolute q[n_trans](each stateSelect=StateSelect.never) = Data.m
    *phi .* I ./ (2*A[cartAxes]) if environment.analysis "Dynamic pressure";
  output Q.CapacityThermalSpecific c_v(stateSelect=StateSelect.never) =
    Data.c_v(T, p) if environment.analysis "Isochoric specific heat capacity";
  output Q.PressureReciprocal beta_T(stateSelect=StateSelect.never) =
    Data.beta_T(T, p) if environment.analysis
    "Isothermal static compressibility";
  //
  // Time constants
  output Q.Time tau_exch_translational(stateSelect=StateSelect.never) =
    halfalpha_zeta if environment.analysis
    "Time constant for translational exchange";
  output Q.Time tau_exch_thermal(stateSelect=StateSelect.never) = halfalpha_theta*N/ if environment.analysis "Time constant for thermal exchange";
  output Q.Time tau_trans_material(stateSelect=StateSelect.never) = noEvent(if
    halfalpha_beta > Modelica.Constants.small then Data.m*beta_T/halfalpha_beta
     else 0) if environment.analysis "Time constant for material transport";
  // **Note that this isn't dependent on length or area---only intensive
  // properties.
  output Q.Time tau_trans_normal[n_trans](each stateSelect=StateSelect.never) =
    fill(halfalpha_beta*N, n_trans) ./ Lprime[cartAxes] if environment.analysis
    "Time constants for normal translational transport";
  // **Note that only for the axes with translational momentum included; others are
  // infinite
  output Q.Time tau_trans_transverse[n_trans](each stateSelect=StateSelect.never)
     = fill(halfalpha_zeta*N, n_trans) ./ Lprime[cartAxes] if environment.analysis
    "Time constants for transverse translational transport";
  // **Note that only for the axes with translational momentum included; others are
  // infinite
  output Q.Time tau_trans_thermal[Axis](each stateSelect=StateSelect.never) =
    fill(halfalpha_theta*N, 3) ./ Lprime if environment.analysis
    "Time constants for thermal transport";
 /*
 /*
  output Q.NumberAbsolute eta(stateSelect=StateSelect.never) = log10(max({
    tau_exch_translational,tau_exch_thermal,max([tau_trans_normal;
    tau_trans_transverse; tau_trans_thermal])})) - log10(min({
    tau_exch_translational,tau_exch_thermal,min([tau_trans_normal;
    tau_trans_transverse; tau_trans_thermal])})) if environment.analysis
    "Range of time constants in order of magnitude";
*/
      //
      // Peclet numbers (only for the axes with translational momentum included; others
      // are zero)
      /* **
  output Q.Number Pe_0[n_trans](each stateSelect=StateSelect.never) = I*
    halfalpha_beta ./ Lprime[cartAxes] if environment.analysis
    "Normal Peclet numbers";
  output Q.Number Pe_12[n_trans](each stateSelect=StateSelect.never) = I*
    halfalpha_zeta ./ Lprime[cartAxes] if environment.analysis
    "Transverse Peclet numbers";
  output Q.Number Pe_therm[n_trans](each stateSelect=StateSelect.never) = I*
    halfalpha_theta ./ Lprime[cartAxes] if environment.analysis
    "Thermal Peclet numbers";
    */
      //
      // Bulk flow rates
      /* **indexing error
  output Q.Force mphiI[n_trans, Orientation](each stateSelect=StateSelect.never)
     = {(if inclTrans[cartWrap(cartAxes[axis] + orientation)] then Data.m*phi[
    transAxes[cartWrap(cartAxes[axis] + orientation)]]*I[axis] else 0) for
    orientation in Orientation, axis in 1:n_trans} if n_trans > 0 and environment.analysis
    "Bulk rate of translational advection";
  */
      /* **
  output Q.Power sTI[n_trans](each stateSelect=StateSelect.never) = T*s*I if
    environment.analysis "Bulk rate of thermal advection";
  //
  // Translational momentum balance
  output Q.Force Ma[n_trans](each stateSelect=StateSelect.never) = M*(der(phi)/U.s
     - environment.a[cartAxes]) if environment.analysis
    "Acceleration force relative to the frame of reference (constant mass)";
  output Q.Force f_exch_adv[n_trans](each stateSelect=StateSelect.never) =
    chemical.mPhidot - Data.m*phi*chemical.Ndot if environment.analysis
    "Acceleration force due to advective exchange";
  output Q.Force f_exch_diff[n_trans](each stateSelect=StateSelect.never) =
    inert.translational.mPhidot + inertDalton.mPhidot if environment.analysis
    "Friction from other species (diffusive exchange)";
  output Q.Force f_trans_adv[n_trans](each stateSelect=StateSelect.never) = Data.m
    *{faces.phi[1][cartAxes[axis], :]*faces[cartAxes[axis], :].Ndot + sum(faces[
    cartWrap(cartAxes[axis] - orientation), :].phi[orientation]*faces[cartWrap(
    cartAxes[axis] - orientation), :].Ndot for orientation in Orientation) for
    axis in 1:n_trans} if environment.analysis
    "Acceleration force due to advective transport";
  output Q.Force f_trans_diff[n_trans](each stateSelect=StateSelect.never) = {sum
    (Sigma(faces[cartWrap(cartAxes[axis] - orientation), :].mPhidot[orientation])
    for orientation in Orientation) - Delta(faces[cartAxes[axis], :].rho)*A[
    cartAxes[axis]] for axis in 1:n_trans} if environment.analysis
    "Friction from other subregions (diffusive transport; includes volume viscosity)";
  //
  // Energy balance
  output Q.Power Ndere(stateSelect=StateSelect.never) = (M*der(phi*phi)/2 + N*
    der(h) - V*der(p))/U.s if environment.analysis
    "Rate of energy storage (internal and kinetic) at constant mass";
  output Q.Power Wdot_exch(stateSelect=StateSelect.never) = -(chemical.mphi*
    chemical.mPhidot/(2*Data.m) + (Data.m*(chemical.hbar - phi*phi/2) - h)*chemical.Ndot)
    if environment.analysis
    "Relative rate of work (internal, flow, and kinetic) done by chemical exchange (advection)";
  output Q.Power Qdot_gen_exch(stateSelect=StateSelect.never) = phi*inert.translational.mPhidot
     + inertDalton.phi*inertDalton.mPhidot if environment.analysis
    "Rate of heat generation due to friction with other species";
  output Q.Power Qdot_exch(stateSelect=StateSelect.never) = inert.thermal.Qdot +
    inertDalton.Qdot if environment.analysis
    "Rate of thermal conduction from other species";
  output Q.Power Wdot_trans(stateSelect=StateSelect.never) = -sum(sum((Data.m*(
    faces.phi[axis, side,0]^2 + faces[axis, side].phi*faces[axis, side].phi)/2 +
    Data.h(faces[axis, side].T, faces[axis, side].rho) - Data.m*phi*phi/2 - h)*
    faces[axis, side].Ndot for side in Side) for axis in Axis) if environment.analysis
    "Relative rate of work (internal, flow, and kinetic) done by advective transport";
  output Q.Power Qdot_gen_trans(stateSelect=StateSelect.never) = sum(faces.phi .*
    faces.mPhidot) if environment.analysis
    "Rate of heat generation due to friction with other subregions";
  output Q.Power Qdot_trans(stateSelect=StateSelect.never) = sum(faces.Qdot)
    if environment.analysis "Rate of thermal conduction from other subregions";
  // Note:  These auxiliary variables should not be used as states (hence
  // StateSelect.never); the structure of the problem should not change if
  // they are included.
*/

      FCSys.Connectors.ChemicalSpecies chemical[n_react](
        each final n_trans=n_trans,
        mu(each start=g_IC),
        sT(each start=h_IC - g_IC)) "Connector for phase change and reactions"
        annotation (Placement(transformation(extent={{-24,4},{-4,24}}),
            iconTransformation(extent={{-45.78,45.55},{-65.78,65.55}})));
      // each final formula=Data.formula,
      // each final m=Data.m,

      //each final phi=phi,

      //   ** final phi(start=fill(phi_IC[cartAxes], n_react)) = chemical_phi,

      Connectors.Inert inert(
        final n_trans=n_trans,
        translational(phi(final start=phi_IC[cartAxes], each final fixed=false)),

        thermal(T(final start=T_IC, final fixed=false)))
        "Connector to directly couple velocities and temperatures of configurations"
        annotation (Placement(transformation(extent={{10,-10},{30,10}}),
            iconTransformation(extent={{61.55,-25.78},{81.55,-45.78}})));

      Connectors.InertDalton inertDalton(
        final n_trans=n_trans,
        V(
          min=0,
          final start=V_IC,
          final fixed=false),
        p(final start=p_IC, final fixed=false),
        phi(start=phi_IC[cartAxes]),
        T(start=T_IC))
        "Connector to exchange translational momentum and thermal energy by diffusion, with additivity of pressure"
        annotation (Placement(transformation(extent={{-10,-30},{10,-10}}),
            iconTransformation(extent={{25.78,-61.55},{45.78,-81.55}})));
      Connectors.Face faces[Axis, Side](
        rho(each start=rho_IC),
        Ndot(start=outerProduct(I_IC, {1,-1})),
        mPhidot(each start=0),
        T(each start=T_IC),
        Qdot(each start=0))
        "Connectors to transport material, translational momentum, and thermal energy through the boundaries"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}}),
            iconTransformation(extent={{-10,-10},{10,10}})));
      //    phi(start={{fill(phi_IC[cartWrap(axis + orientation - 1)], 2) for         orientation in Orientation} for axis in Axis}),

      // Geometric parameters
    protected
      outer parameter Q.Length L[Axis] "Length" annotation (missingInnerMessage
          ="This model should be used within a subregion model.");
      outer parameter Q.Area A[Axis] "Cross-sectional area" annotation (
          missingInnerMessage=
            "This model should be used within a subregion model.");
      outer parameter Q.Length Lprime[Axis]
        "Effective cross-sectional area per length" annotation (
          missingInnerMessage="This model should be used within a phase model.");
      outer parameter Boolean inclTrans[Axis]
        "true, if each component of translational momentum is included"
        annotation (missingInnerMessage=
            "This model should be used within a subregion model.");
      outer parameter Integer n_trans
        "Number of components of translational momentum" annotation (
          missingInnerMessage=
            "This model should be used within a subregion model.");
      outer parameter Integer transAxes[:]
        "Translational-momentum-component indices of the Cartesian axes"
        annotation (missingInnerMessage=
            "This model should be used within a subregion model.");
      outer parameter Integer cartAxes[:]
        "Cartesian-axis indices of the components of translational momentum"
        annotation (missingInnerMessage=
            "This model should be used within a subregion model.");
      // Note:  The size is n_trans, but it can't be specified here due to an error
      // in Dymola 7.4.
      final parameter Boolean upstream[Axis]={upstreamX,upstreamY,upstreamZ}
        "true, if each Cartesian axis uses upstream discretization"
        annotation (HideResult=true);
      final parameter Conservation consTrans[Axis]={consTransX,consTransY,
          consTransZ} "Formulation of the translational conservation equations"
        annotation (HideResult=true);
      final parameter InitTranslational initTrans[Axis]={initTransX,initTransY,
          initTransZ} "Initialization methods for translational momentum"
        annotation (HideResult=true);

      // Aliases to the chemical connector variables
      // These are necessary because Dymola does not fully support zero-sized
      // connectors.
      //Q.Potential chemical_mu[n_react](each start=g_IC) "Electrochemical potential";
      //Q.Current chemical_Ndot[n_react] "Diffusion current";
      //Q.Velocity chemical_phi[n_react, n_trans](start=fill(phi_IC[cartAxes],
      //      n_react)) "Velocity";
      Q.Force chemical_mPhidot[n_react, n_trans] "Advective force";
      //Q.PotentialAbsolute chemical_sT[n_react](each start=h_IC - g_IC) "Specific entropy-temperature product";
      //Q.Power chemical_Qdot[n_react] "Rate of thermal advection";

      outer Conditions.Environment environment "Environmental conditions";

    initial equation
      // Check that the initialization methods are valid.
      assert(initMaterial <> initEnergy or initMaterial == InitScalar.None or
        consMaterial == Conservation.Static or consEnergy == Conservation.Static,
        "The initialization methods for material and energy must be different (unless None).");
      if not Data.isCompressible then
        assert(initMaterial <> InitScalar.Pressure and initMaterial <>
          InitScalar.PressureSS or consMaterial == Conservation.IC, "The material is incompressible,
yet the initialization method for material involves pressure.");
        assert(initEnergy <> InitScalar.Pressure and initEnergy <> InitScalar.PressureSS
           or consEnergy == Conservation.IC, "The material is incompressible,
yet the initialization method for energy involves pressure.");
        if not Data.hasThermalExpansion then
          assert(initMaterial <> InitScalar.Density and initMaterial <>
            InitScalar.DensitySS or consMaterial == Conservation.IC, "The material is isochoric,
yet the initialization method for material involves density.");
          assert(initEnergy <> InitScalar.Density and initEnergy <> InitScalar.DensitySS
             or consMaterial == Conservation.IC, "The material is isochoric,
yet the initialization method for energy involves density.");
        end if;
      end if;

      // Material
      if consMaterial == Conservation.IC then
        // Ensure that a condition is selected, since the state is prescribed.
        assert(initMaterial <> InitScalar.None, "The material state is prescribed,
yet its condition is not defined.
Choose a condition besides None.");
      elseif consMaterial == Conservation.Dynamic then
        // Initialize since there's a time-varying state.
        if initMaterial == InitScalar.Amount then
          N = N_IC;
        elseif initMaterial == InitScalar.AmountSS then
          der(N) = 0;
        elseif initMaterial == InitScalar.Density then
          rho = rho_IC;
        elseif initMaterial == InitScalar.DensitySS then
          der(rho)/U.s = rho_IC;
        elseif initMaterial == InitScalar.Volume then
          V = V_IC;
        elseif initMaterial == InitScalar.VolumeSS then
          der(V) = 0;
        elseif initMaterial == InitScalar.Pressure then
          p = p_IC;
        elseif initMaterial == InitScalar.PressureSS then
          der(p) = 0;
        elseif initMaterial == InitScalar.Temperature then
          T = T_IC;
        elseif initMaterial == InitScalar.TemperatureSS then
          der(T) = 0;
        elseif initMaterial == InitScalar.SpecificEnthalpy then
          h = h_IC;
        elseif initMaterial == InitScalar.SpecificEnthalpySS then
          der(h) = 0;
        elseif initMaterial == InitScalar.PotentialGibbs then
          h - sT = g_IC;
        else
          //if initMaterial == InitScalar.PotentialGibbsSS then
          der(h - sT) = 0;
          // Else there is no initial equation because
          // initMaterial == InitScalar.None or
          // consMaterial == Conservation.Static.
        end if;
      end if;

      // Velocity
      for axis in Axis loop
        if inclTrans[axis] then
          if consTrans[axis] == Conservation.IC then
            // Ensure that a condition is selected, since the state is
            // prescribed.
            assert(initTrans[axis] <> InitTranslational.None,
              "The state for the " + {"x","y","z"}[axis] + "-axis component of translational momentum is prescribed,
yet its condition is not defined.
Choose any condition besides None.");
          elseif consTrans[axis] == Conservation.Dynamic then
            // Initialize since there's a time-varying state.
            if initTrans[axis] == InitTranslational.Velocity then
              phi[transAxes[axis]] = phi_IC[axis];
            elseif initTrans[axis] == InitTranslational.VelocitySS then
              der(phi[transAxes[axis]]) = 0;
            elseif initTrans[axis] == InitTranslational.Current then
              I[transAxes[axis]] = I_IC[axis];
            elseif initTrans[axis] == InitTranslational.CurrentSS then
              der(I[transAxes[axis]]) = 0;
              // Else there is no initial equation because
              // initTrans[axis] == InitTranslational.None or
              // initTrans[axis] == Conservation.Static.
            end if;
          end if;
        end if;
      end for;

      // Energy
      if consEnergy == Conservation.IC then
        // Ensure that a condition is selected, since the state is prescribed.
        assert(initEnergy <> InitScalar.None, "The energy state is prescribed,
yet its condition is not defined.
Choose a condition besides None.");
      elseif consEnergy == Conservation.Dynamic then
        // Initialize since there's a time-varying state.
        if initEnergy == InitScalar.Amount then
          N = N_IC;
        elseif initEnergy == InitScalar.AmountSS then
          der(N) = 0;
        elseif initMaterial == InitScalar.Density then
          rho = rho_IC;
        elseif initMaterial == InitScalar.DensitySS then
          der(rho) = 0;
        elseif initEnergy == InitScalar.Volume then
          V = V_IC;
        elseif initEnergy == InitScalar.VolumeSS then
          der(V) = 0;
        elseif initEnergy == InitScalar.Pressure then
          p = p_IC;
        elseif initEnergy == InitScalar.PressureSS then
          der(p) = 0;
        elseif initEnergy == InitScalar.Temperature then
          T = T_IC;
        elseif initEnergy == InitScalar.TemperatureSS then
          der(T) = 0;
        elseif initEnergy == InitScalar.SpecificEnthalpy then
          h = h_IC;
        elseif initEnergy == InitScalar.SpecificEnthalpySS then
          der(h) = 0;
        elseif initEnergy == InitScalar.PotentialGibbs then
          h - sT = g_IC;
        else
          //if initEnergy == InitScalar.PotentialGibbsSS then
          der(h - sT) = 0;
          // Else there is no initial equation because
          // initEnergy == InitScalar.None or
          // consEnergy == Conservation.Static.
        end if;
      end if;

    equation
      // Aliases (only for clarity)
      p = inertDalton.p;
      V = inertDalton.V;
      V = v*N;
      1 = rho*v;
      M = Data.m*N;
      T = inert.thermal.T;
      phi = inert.translational.phi;
      N*phi = L[cartAxes] .* I;
      p_faces = {{Data.p_Tv(faces[axis, side].T, 1/faces[axis, side].rho) for
        side in Side} for axis in Axis};

      // Thermodynamic correlations
      if invertEOS then
        p = Data.p_Tv(T, v);
      else
        v = Data.v_Tp(T, p);
      end if;
      h = Data.h(T, p);
      sT = Data.s(T, p)*T;

      // Exchange
      // --------
      // Material
      if Data.z == 0 then
        for i in 1:n_react loop
          h = chemical[i].mu + sT "Direct application of potential";
          //    chemical_Ndot[i] = 2 - 2*exp(chemical_mu[i]/T) "Diffusion";
          // **Add rate coefficient
        end for;
      else
        for i in 1:n_react loop
          Data.z*chemical[i].Ndot = 2 - 2*exp(chemical[i].mu/T) "Diffusion";
          // **Add rate coefficient
        end for;
      end if;
      //
      // Translational momentum
      for i in 1:n_react loop
        chemical[i].phi = phi "Advected property upon outflow";
        chemical_mPhidot[i, :] = Data.m*actualStream(chemical[i].phi)*chemical[
          i].Ndot "Advection";
      end for;
      mu*inertDalton.mPhidot = N*(inertDalton.phi - phi) "Diffusion";
      //
      // Thermal energy
      chemical.sT = fill(sT, n_react) "Advected property upon outflow";
      nu*inertDalton.Qdot = N*(inertDalton.T - T) "Diffusion";

      // Transport
      for axis in Axis loop
        for side in Side loop
          // Material
          eta*faces[axis, side].Ndot = Lprime[axis]*(faces[axis, side].rho -
            rho)*(if upstream[axis] and inclTrans[axis] then 1 + exp(-inSign(
            side)*I[transAxes[axis]]*eta*v/(2*Lprime[axis])) else 2);

          // Translational momentum
          beta*faces[axis, side].mPhidot[1] = Lprime[axis]*(faces[axis, side].phi[
            1] - (if inclTrans[axis] then phi[transAxes[axis]] else 0))*(if
            inclTrans[axis] and upstream[axis] then 1 + exp(-inSign(side)*I[
            transAxes[axis]]*beta*Data.m/(2*Lprime[axis])) else 2) "Normal";
          zeta*faces[axis, side].mPhidot[2] = Nu_Phi[axis]*Lprime[axis]*(faces[
            axis, side].phi[2] - (if inclTrans[cartWrap(axis + 1)] then phi[
            transAxes[cartWrap(axis + 1)]] else 0))*(if inclTrans[axis] and
            upstream[axis] then 1 + exp(-inSign(side)*I[transAxes[axis]]*zeta*
            Data.m/(2*Lprime[axis])) else 2) "1st transverse";
          zeta*faces[axis, side].mPhidot[3] = Nu_Phi[axis]*Lprime[axis]*(faces[
            axis, side].phi[3] - (if inclTrans[cartWrap(axis + 2)] then phi[
            transAxes[cartWrap(axis + 2)]] else 0))*(if inclTrans[axis] and
            upstream[axis] then 1 + exp(-inSign(side)*I[transAxes[axis]]*zeta*
            Data.m/(2*Lprime[axis])) else 2) "2nd transverse";

          // Thermal energy
          theta*faces[axis, side].Qdot = Nu_Q*Lprime[axis]*(faces[axis, side].T
             - T)*(if inclTrans[axis] and upstream[axis] then 1 + exp(-inSign(
            side)*I[transAxes[axis]]*theta*Data.c_v(T, p)/(2*Lprime[axis]))
             else 2);
        end for;
      end for;

      // Material dynamics
      if consMaterial == Conservation.IC then
        // Apply the IC forever (material not conserved).
        if initMaterial == InitScalar.Amount then
          N = N_IC;
        elseif initMaterial == InitScalar.AmountSS then
          der(N) = 0;
        elseif initMaterial == InitScalar.Density then
          rho = rho_IC;
        elseif initMaterial == InitScalar.DensitySS then
          der(rho) = 0;
        elseif initMaterial == InitScalar.Volume then
          V = V_IC;
        elseif initMaterial == InitScalar.VolumeSS then
          der(V) = 0;
        elseif initMaterial == InitScalar.Pressure then
          p = p_IC;
        elseif initMaterial == InitScalar.PressureSS then
          der(p) = 0;
        elseif initMaterial == InitScalar.Temperature then
          T = T_IC;
        elseif initMaterial == InitScalar.TemperatureSS then
          der(T) = 0;
        elseif initMaterial == InitScalar.SpecificEnthalpy then
          h = h_IC;
        elseif initMaterial == InitScalar.SpecificEnthalpySS then
          der(h) = 0;
        elseif initMaterial == InitScalar.PotentialGibbs then
          h - sT = g_IC;
        else
          //if initMaterial == InitScalar.PotentialGibbsSS then
          der(h - sT) = 0;
          // Note:  initMaterial == InitScalar.None can't occur due to an
          // assertion.
        end if;
      else
        (if consMaterial == Conservation.Dynamic then der(N)/U.s else 0) = sum(
          chemical.Ndot) + sum(faces.Ndot) "Material conservation";
      end if;

      // Translational dynamics
      for axis in 1:n_trans loop
        if consTrans[cartAxes[axis]] == Conservation.IC then
          // Apply the IC forever (translational momentum isn't conserved along
          // this axis).
          if initTrans[cartAxes[axis]] == InitTranslational.Velocity then
            phi[axis] = phi_IC[cartAxes[axis]];
          elseif initTrans[cartAxes[axis]] == InitTranslational.VelocitySS then
            der(phi[axis]) = 0;
          elseif initTransX == InitTranslational.Current then
            I[axis] = I_IC[cartAxes[axis]];
          elseif initTrans[cartAxes[axis]] == InitTranslational.CurrentSS then
            der(I[axis]) = 0;
            // Note:  initTrans[cartAxes[axis]] == InitTranslational.None can't
            // occur due to an assertion.
          end if;
        else
          (if consTrans[cartAxes[axis]] == Conservation.Dynamic then der(M*phi[
            axis])/U.s else 0) + Delta(p_faces[cartAxes[axis], :])*A[cartAxes[
            axis]] + M*environment.a[cartAxes[axis]] + N*Data.z*environment.E[
            cartAxes[axis]] = sum(chemical_mPhidot[:, axis]) + inert.translational.mPhidot[
            axis] + inertDalton.mPhidot[axis] + sum(faces[cartWrap(cartAxes[
            axis] - orientation + 1), :].phi[orientation]*faces[cartWrap(
            cartAxes[axis] - orientation + 1), :].Ndot*Data.m + Sigma(faces[
            cartWrap(cartAxes[axis] - orientation + 1), :].mPhidot[orientation])
            for orientation in Orientation)
            "Conservation of translational momentum";
        end if;
      end for;

      // Thermal dynamics
      if consEnergy == Conservation.IC then
        // Apply the IC forever (energy not conserved).
        if initEnergy == InitScalar.Amount then
          N = N_IC;
        elseif initEnergy == InitScalar.AmountSS then
          der(N) = 0;
        elseif initMaterial == InitScalar.Density then
          rho = rho_IC;
        elseif initMaterial == InitScalar.DensitySS then
          der(rho) = 0;
        elseif initEnergy == InitScalar.Volume then
          V = V_IC;
        elseif initEnergy == InitScalar.VolumeSS then
          der(V) = 0;
        elseif initEnergy == InitScalar.Pressure then
          p = p_IC;
        elseif initEnergy == InitScalar.PressureSS then
          der(p) = 0;
        elseif initEnergy == InitScalar.Temperature then
          T = T_IC;
        elseif initEnergy == InitScalar.TemperatureSS then
          der(T) = 0;
        elseif initEnergy == InitScalar.SpecificEnthalpy then
          h = h_IC;
        elseif initEnergy == InitScalar.SpecificEnthalpySS then
          der(h) = 0;
        elseif initEnergy == InitScalar.PotentialGibbs then
          h - sT = g_IC;
        else
          //if initEnergy == InitScalar.PotentialGibbsSS then
          der(h - sT) = 0;
          // Note:  initEnergy == InitScalar.None can't occur due to an
          // assertion.
        end if;
      else
        (if consEnergy == Conservation.Dynamic then (der(N*h) - V*der(p) + der(
          M*phi*phi)/2)/U.s else 0) = (chemical.mu)*chemical.Ndot + sum(
          actualStream(chemical[i].sT)*chemical[i].Ndot for i in 1:n_react) +
          inert.translational.phi*inert.translational.mPhidot + inert.thermal.Qdot
           + inertDalton.phi*inertDalton.mPhidot + inertDalton.Qdot + sum(sum((
          Data.h(faces[axis, side].T, p_faces[axis, side]) + faces[axis, side].phi
          *faces[axis, side].phi*Data.m/2)*faces[axis, side].Ndot + faces[axis,
          side].phi*faces[axis, side].mPhidot for side in Side) for axis in
          Axis) + sum(faces.Qdot) "Energy conservation";
        //**+ sum(chemical_phi[:, ax] .^ 2 for ax in 1:n_trans)*(Data.m/2)
        // sum(actualStream(chemical[i].sT)*chemical[i].Ndot for i in 1:n_react)
        // is used because Dymola 7.4 isn't able to interpret the vectorized
        // form, actualStream(chemical.sT)*chemical.Ndot.

      end if;
      annotation (
        defaultComponentPrefixes="replaceable",
        Documentation(info="<html>
    <p>This model is based on the following fixed assumptions.  Other assumptions are
    optional via the parameters.
    <ol>
       <li>All faces are rectangular.
       <li>The material is orthorhombic.  This implies that a gradient which induces diffusion
       along an axis does not induce diffusion along axes orthogonal to it
       [<a href=\"modelica://FCSys.UsersGuide.References\">Bejan2006</a>,
       pp. 691&ndash;692].</li>
       <li>The coordinate system (x, y, z) is aligned with the principle
       axes of transport.  For example, if the species is stratified, the
       layers must be parallel to one of the planes in the rectilinear
       grid.</li>
       <li>The factors that may cause anisotropic behavior (<b><i>k</i></b>)
          are common to material, translational, and thermal transport.</li>
       <li>There is no radiative heat transfer.</li>
       <li>Rotational momentum is not exchanged, transported, or stored.</li>
       <li>For the purpose of the material, translational momentum, and energy balances, the
       cross sectional areas of the faces are assumed to be the full cross-sectional
       areas of the subregion.  If multiple phases are present, then areas are
       actually smaller.</li>
       <li>Relativistic effects are negligible.</li>
    </ol></p>

    <p>Figure 1 shows how instances of
    <a href=\"modelica://FCSys.Subregions.Species\">Species</a> models (derived from this
    model) are
    connected within a <a href=\"modelica://FCSys.Subregions\">Subregion</a>.  The
    generalized resistances (<i>R</i>) affect the flow rates of translational momentum and
    thermal energy associated with differences in velocity and temperature (respectively) between
    each species and a common node.  This exchange is diffusive.

    <p>Translational momentum and enthalpy are advected as material is exchanged in a chemical
    reaction.  This occurs at the velocity and massic enthalpy of the reactants (source
    species), where the reactant/product designation depends on the current conditions.
    If species are connected through
    a <a href=\"modelica://FCSys.Subregions.Reaction\">Reaction</a> model, then the material
    states (e.g., amounts of each material) are directly coupled to impose chemical
    equilibrium.
    This reduces the DAE index by one in accordance with Gibbs' phase rule
    [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>].
    Resistance is not included directly in the reaction equations;
    the reaction rate is determined solely by
    the transport equations.</p>

    <p align=center><img src=\"modelica://FCSys/resources/documentation/Subregions/Species/Species/Exchange.png\">
<br>Figure 1:  Exchange of a quantity (translational momentum or thermal energy) among species
    (A, B, and C) within a subregion.</p>

    <p>Figure 2 shows how <a href=\"modelica://FCSys.Subregions.Species\">Species</a>
    instances of the same type are connected between neighboring
    <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> instances.
    Normal and transverse translational momentum and thermal energy are transported by both advection and diffusion.
    Upstream discretization is applied if it is enabled via the <code>upstreamX</code>,
    etc. parameters.</p>

    <p align=center><img src=\"modelica://FCSys/resources/documentation/Subregions/Species/Species/Transport.png\">
<br>Figure 2:  Transport of a quantity associated with the same species
    between subregions (1 and 2).</p>

    <p>All <a href=\"modelica://FCSys.Subregions.Species\">Species</a> instances
    within a <a href=\"modelica://FCSys.Subregions.Phases\">Phase</a> are joined by Dalton's law (see the
    <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector), as shown
    in Figure 3a.  The pressures are additive, and each species is assumed to exist at the
    total extensive volume of the phase.  Within a <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a>,
    the <a href=\"modelica://FCSys.Subregions.Phases\">Phases</a> are combined by Amagat's law (see the
    <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector), as shown
    in Figure 3b.  The volumes are additive, and each species is assumed to exist at the
    total pressure in the subregion.</p>

    <table border=0 cellspacing=0 cellpadding=2 align=center class=noBorder style=\"margin-left: auto; margin-right: auto;\">
      <tr align=center class=noBorder>
        <td align=center class=noBorder style=\"margin-left: auto; margin-right: auto;\">
          <img src=\"modelica://FCSys/resources/documentation/Subregions/Species/Species/SharePressure.png\">
<br>a:  Pressures of species (A, B, and C) are additive within a phase.
        </td>
        <td align=center class=noBorder style=\"margin-left: auto; margin-right: auto;\">
          <img src=\"modelica://FCSys/resources/documentation/Subregions/Species/Species/ShareVolume.png\">
<br>b:  Volumes of phases (I, II, and III) are additive within a subregion.
        </td>
      </tr>
      <tr align=center class=noBorder style=\"margin-left: auto; margin-right: auto;\">
        <td colspan=2 align=center class=noBorder>Figure 3: Methods of attributing pressure and volume.</td>
      </tr>
    </table>

  <p>**Translational momentum and energy are advected using the <code>semiLinear</code> operator.
  The rate of advection of translational momentum is the
  product of the velocity of the source and the rate of mass
  (<i>m</i> &phi; <i>N&#775;</i>).  The rate of thermal advection is the
  product of the massic enthalpy of the source and the rate of mass
  (<i>m</i> <i>h&#772;</i> <i>N&#775;</i>).  If there multiple sources, then
  their contributions are additive.  If there are multiple sinks, then
  the flow is split on a mass basis.</p>

    <p>Notes regarding the parameters:
    <ul>
    <li>Here (and in the rest of <a href=\"modelica://FCSys\">FCSys</a>), the \"specific\"
    adjective means that the following extensive quantity is divided by particle number.
    (\"Massic\" indicates a quantity divided by mass.)</li>
    <li>In general, if material resistivity, dynamic compressibility, fluidity, or thermal resistivity is zero, then
    it should be set as <code>final</code> so that index reduction may be performed.
    If two <a href=\"modelica://FCSys.Subregions.Species\">Species</a> instances
    are connected through their <code>inert</code> connectors or faces (<code>xNegative</code>,
    <code>xPositive</code>, etc.) and both have zero generalized resistivities for a
    quantity, then index reduction [<a href=\"modelica://FCSys.Subregions.Species\">Mattsson1993B</a>] is necessary.</li>
    <li>Even if an initialization parameter is not selected to be used explicitly,
    it may be used a guess value.</li>
    <li>The <b><i>k</i></b> factor can be used to account for the effects of porosity and tortousity
    on transport.
    It should be changed directly with effective area and inversely with effective length.
    The factor may reflect anisotropic properties; it is a vector with independent components
    for each axis. It affects all of the diffusive transport rates (normal, transverse, and
    thermal) by the same factor.  By default, its components are unity.</li>
    <li>By default, only the x-axis component of translational momentum is included.</li>
    <li>If a state is prescribed, then the
    associated initial condition (IC) will be applied forever.  The
    corresponding conservation equation will not be imposed.
    If <code>consMaterial</code>, <code>consTransX</code>, <code>consTransY</code>, or <code>consTransZ</code> is
    <code>Conservation.IC</code>, then there may be a secondary effect on the energy conservation equation
    and thus temperature.
    In that case, it may be helpful to set <code>consEnergy</code> to <code>Conservation.IC</code> so that
    the energy conservation equation is not imposed.</li>
    <li>If a subregion does not contain any compressible species, then pressure must be prescribed.
    Set <code>consMaterial</code> to <code>Conservation.IC</code> and <code>initMaterial</code>
    to <code>InitScalar.Pressure</code> for one of the species.</li>
    <li>The <code>start</code> values of the initial conditions for pressure and temperature
    (<i>p</i><sub>IC</sub> and <i>T</i><sub>IC</sub>) are the global default pressure and
    temperature (via the <code>outer</code> instance of the <a href=\"modelica://FCSys.Conditions.Environment\">Environment</a> model).
    The <code>start</code> values of the initial conditions for
    other intensive properties (<i>v</i><sub>IC</sub>, <i>h</i><sub>IC</sub>, and
    &mu;<sub>IC</sub>) are related to the initial pressure and temperature
    by the characteristics of the species.  The <code>start</code> value of the
    initial condition for the extensive volume (<i>V</i><sub>IC</sub>) is the volume of the
    subregion, and the <code>start</code> value for particle number (<i>N</i><sub>IC</sub>)
    is related to it via the material characteristics (<code>Data</code>) and the initial pressure and temperature.
    In order to apply other values for any of these initial conditions,
    it may be necessary to do so before translating the model.</li>
    <li>If upstream discretization is not used (<code>upstreamX=false</code>,
    etc.), then the central difference scheme is used.
    <li>If <code>invertEOS</code> is <code>true</code>, then the equation of state is implemented with pressure
    as a function of temperature and specific volume.  Otherwise, specific volume is a function of temperature
    and pressure.</li></p>

    <p>In evaluating the dynamics of a phase, it is typically assumed that all of the species
    exist at the same velocity and temperature.  The translational and thermal time constants
    are usually much shorter than the time span of interest due to the very small coupling
    resistances.  This assumption can be applied in the model by connecting the <code>common</code>
    connectors of the species, which will **

    **    It will cause index reduction during translation.</p>

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
              extent={{-80,80},{80,-80}},
              lineColor={127,127,127},
              pattern=LinePattern.Dash,
              fillColor={225,225,225},
              fillPattern=FillPattern.Solid)}));
    end Species;

    package BaseClasses "Base classes (not generally for direct use)"

      extends Modelica.Icons.BasesPackage;

      type Conservation = enumeration(
          IC "Initial condition imposed forever (no conservation)",
          Static "Static (conservation without storage)",
          Dynamic "Dynamic (conservation with storage)")
        "Options for a conservation equation";
      type InitScalar = enumeration(
          None "No initialization",
          Amount "Prescribed amount",
          AmountSS "Steady-state amount",
          Density "Prescribed density",
          DensitySS "Steady-state density",
          Volume "Prescribed volume",
          VolumeSS "Steady-state volume",
          Pressure "Prescribed pressure",
          PressureSS "Steady-state pressure",
          Temperature "Prescribed temperature",
          TemperatureSS "Steady-state temperature",
          SpecificEnthalpy "Prescribed specific enthalpy",
          SpecificEnthalpySS "Steady-state specific enthalpy",
          PotentialGibbs "Prescribed Gibbs potential",
          PotentialGibbsSS "Steady-state Gibbs potential")
        "Methods of initializing scalar quantities (material and energy)";

      type InitTranslational = enumeration(
          None "No initialization",
          Velocity "Prescribed velocity",
          VelocitySS "Steady-state velocity",
          Current "Prescribed current",
          CurrentSS "Steady-state current")
        "Methods of initializing translational momentum";

    end BaseClasses;

    model Test

      Connectors.ChemicalSpecies chemical(formula="",m=1)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    equation
      chemical.Ndot = chemical.mu;
      actualStream(chemical.sT) = 1;
      chemical.sT = 1;
    end Test;

    model TestTest

      Test test
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    equation
      test.chemical.axis = 1;
    end TestTest;
  end Species;

  model PhaseBoundary
    "Phase boundary (adapter between Amagat and Dalton mixtures)"

    // extends FCSys.BaseClasses.Icons.Names.Top7;

    Connectors.InertAmagat inertAmagat(final n_trans=n_trans)
      "<html>Connector for volume, translational momentum, and thermal energy&mdash;with Amagat's law</html>"
      annotation (Placement(transformation(extent={{60,-80},{80,-60}}),
          iconTransformation(extent={{52,-132},{72,-112}})));
    Connectors.InertDalton inertDalton(final n_trans=n_trans)
      "<html>Connector for volume, translational momentum, and thermal energy&mdash;with Dalton's law</html>"
      annotation (Placement(transformation(extent={{30,-50},{50,-30}}),
          iconTransformation(extent={{25.78,-81.55},{45.78,-61.55}})));

  protected
    outer parameter Integer n_trans
      "Number of components of translational momentum" annotation (
        missingInnerMessage=
          "This model should be used within a subregion model.");

  equation
    // Equal intensive properties
    inertAmagat.phi = inertDalton.phi;
    inertAmagat.T = inertDalton.T;

    // Static balances
    0 = inertAmagat.p + inertDalton.p "Pressure";
    0 = inertAmagat.V + inertDalton.V "Volume";

    // Conservation (without storage)
    zeros(n_trans) = inertAmagat.mPhidot + inertDalton.mPhidot
      "Translational momentum";
    0 = inertAmagat.Qdot + inertDalton.Qdot "Energy";
    annotation (
      Documentation(info="<html><p>This model is essentially an
    adapter between the <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> and
    <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connectors.  Inside a phase,
    Dalton's law is applied.  Outside, Amagat's law is applied.</p>

    <p>See also the documentation in the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
              100,100}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-180,-180},{180,
              180}}), graphics={
          Rectangle(
            extent={{-170,140},{170,180}},
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
            extent={{-170,140},{170,180}},
            textString="%name",
            lineColor={0,0,0})}));
  end PhaseBoundary;

  model Reaction "Electrochemical reaction"
    import FCSys.BaseClasses.Utilities.Chemistry.charge;
    import FCSys.BaseClasses.Utilities.Chemistry.stoich;
    extends FCSys.BaseClasses.Icons.Names.Top2;

    // General parameters
    parameter Axis axis=Axis.x "Axis of the electric field";
    // This is merely propagated to the charged species for their momentum
    // balances.
    parameter Quantities.Capacitance C=U.F "Electrical capacitance";
    parameter Q.Number alpha=0.5
      "<html>Charge transfer coefficient (&alpha;)</html>";
    parameter Integer n_neut(min=1) = 0 "Number of neutral species"
      annotation (Dialog(connectorSizing=true));
    // The default is zero for connectorSizing.
    final parameter Integer n_spec=n_neut + 2 "Number of species";

    // Initialization parameters
    parameter InitElectrical init=InitElectrical.SS "Method of initialization"
      annotation (Evaluate=true, Dialog(tab="Initialization"));
    parameter Quantities.Amount Z_IC(start=0)
      "<html>Initial charge (<i>Z</i><sub>IC</sub>)</html>"
      annotation (Dialog(tab="Initialization"));
    parameter Quantities.Potential w_IC(start=0)
      "<html>Initial potential (<i>w</i><sub>IC</sub>)</html>"
      annotation (Dialog(tab="Initialization"));
    // Note:  These parameters are left enabled even if they aren't used to
    // explicitly initialize the state, since they're used as guess values.

    Quantities.Amount Z(
      min=-Modelica.Constants.inf,
      final start=Z_IC,
      final fixed=false) "Charge";
    Quantities.Potential w(final start=w_IC,final fixed=false)
      "Electrical potential";
    Quantities.Current Ndot "Reaction rate";
    Quantities.MomentumTranslationalSpecific phi[n_trans]
      "Conversion specific mass-velocity product";
    Quantities.PotentialAbsolute sT
      "Conversion specific entropy-temperature product";

    FCSys.Connectors.ChemicalReaction positive(final n_trans=n_trans, final
        axis=axis) "Connector for the charged species on the positive side"
      annotation (Placement(transformation(extent={{90,-10},{110,10}})));
    FCSys.Connectors.ChemicalReaction chemical[n_neut](each final n_trans=
          n_trans) "Connectors for the neutral species"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    FCSys.Connectors.ChemicalReaction negative(final n_trans=n_trans, final
        axis=axis) "Connector for the charged species on the negative side"
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

  protected
    outer parameter Integer n_trans
      "Number of components of translational momentum" annotation (
        missingInnerMessage=
          "This model should be used within a subregion model.");
    final Real n[n_spec]=stoich(cat(
          1,
          {positive.formula},
          chemical.formula,
          {negative.formula})) "Stoichiometric coefficients";
    // This cannot be a constant Integer in Dymola 7.4 since the formulas
    // are propagated through the connectors.

  initial equation
    if init == InitElectrical.SS then
      der(Z) = 0;
    elseif init == InitElectrical.Charge then
      Z = Z_IC;
    elseif init == InitElectrical.Potential then
      w = w_IC;
    end if;

  equation
    // Capacitor
    C*w = Z;

    // Electrochemical potentials
    positive.mu = alpha*(n[2:n_spec - 1]*chemical.mu + w);
    negative.mu = (alpha - 1)*(n[2:n_spec - 1]*chemical.mu + w);

    // Advection
    positive.mPhidot = positive.m*semiLinear(
        positive.Ndot,
        phi,
        positive.phi) "Translational momentum";
    positive.Qdot = semiLinear(
        positive.Ndot,
        sT,
        positive.sT) "Thermal energy";
    for i in 1:n_neut loop
      chemical[i].mPhidot = chemical[i].m*semiLinear(
          chemical[i].Ndot,
          phi,
          chemical[i].phi) "Translational momentum";
      chemical[i].Qdot = semiLinear(
          chemical[i].Ndot,
          sT,
          chemical[i].sT) "Thermal energy";
    end for;
    negative.mPhidot = negative.m*semiLinear(
        negative.Ndot,
        phi,
        negative.phi) "Translational momentum";
    negative.Qdot = semiLinear(
        negative.Ndot,
        sT,
        negative.sT) "Thermal energy";

    // Conservation
    der(Z) = charge(positive.formula)*(positive.Ndot + n[1]*Ndot)
      "Charged species on positive side";
    zeros(n_neut) = chemical.Ndot + n[2:n_spec - 1]*Ndot
      "Neutral species (no storage)";
    der(Z) = -charge(negative.formula)*(negative.Ndot + n[n_spec]*Ndot)
      "Charged species on negative side";
    for ax in 1:n_trans loop
      0 = positive.mPhidot[ax] + sum(chemical.mPhidot[ax]) + negative.mPhidot[
        ax] "Translational momentum (no storage)";
    end for;
    0 = positive.Qdot + sum(chemical.Qdot) + negative.Qdot
      "Energy (no storage)";

    // This model is marked as structurally incomplete because n_neut=0 by
    // default for connectorSizing.
    annotation (
      structurallyIncomplete=true,
      Documentation(info="<html>
    <p>The charge species are stored but without their associated translational momentum
    and thermal energy.  The translational momentum and thermal energy is immediately
    advected from the reactants to the products.</p>
    
    <p>The stoichiometry is determined automatically from the chemical formulas
    of the connected species.  No intermediate species are considered.  Each reaction must be
    completely and uniquely defined by the connected species.
    Otherwise an error message is given.</p>
    </html>"),
      Icon(graphics={
          Line(
            points={{-90,0},{-20,0}},
            color={208,104,0},
            smooth=Smooth.None),
          Line(
            points={{-20,40},{-20,-40}},
            color={208,104,0},
            smooth=Smooth.None,
            thickness=0.5),
          Line(
            points={{20,0},{90,0}},
            color={208,104,0},
            smooth=Smooth.None),
          Line(
            points={{20,40},{20,-40}},
            color={208,104,0},
            smooth=Smooth.None,
            thickness=0.5)}),
      Diagram(graphics));
  end Reaction;

  model Volume "Model to establish a fixed total volume"
    // extends FCSys.BaseClasses.Icons.Names.Top7;

    Connectors.InertAmagat inert(final n_trans=n_trans)
      "Connector for translational momentum and thermal energy, with additivity of volume"
      annotation (Placement(transformation(extent={{60,-80},{80,-60}}),
          iconTransformation(extent={{100,-120},{120,-100}})));

    outer parameter Q.Volume V "Volume" annotation (HideResult=true,
        missingInnerMessage=
          "This model should be used within a subregion model.");
    // Note:  These must be public in Dymola 7.4, so HideResult is set true
    // instead.

  protected
    outer parameter Integer n_trans
      "Number of components of translational momentum" annotation (
        missingInnerMessage=
          "This model should be used within a subregion model.");

  equation
    // Specified volume
    V = inert.V;

    // Conservation (without storage)
    zeros(n_trans) = inert.mPhidot "Translational momentum";
    0 = inert.Qdot "Energy";
    annotation (
      Documentation(info="<html><p>This model uses an <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector that imposes
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
              100,100}}), graphics));
  end Volume;

  package BaseClasses "Base classes (not generally for direct use)"

    extends Modelica.Icons.BasesPackage;

    partial model PartialSubregion
      "Partial model for multi-dimensional and multi-species storage, transport, and exchange"
      import FCSys.BaseClasses.Utilities.cartWrap;
      import FCSys.BaseClasses.Utilities.countTrue;
      import FCSys.BaseClasses.Utilities.enumerate;
      import FCSys.BaseClasses.Utilities.index;
      // extends FCSys.BaseClasses.Icons.Names.Top3;

      // Geometric parameters
      inner parameter Q.Length L[Axis](each min=Modelica.Constants.small) = {U.cm,
        U.cm,U.cm} "<html>Length (<b>L</b>)</html>"
        annotation (Dialog(group="Geometry"));
      final inner parameter Q.Volume V=product(L) "Volume";

      // Assumptions
      // -----------
      // Included components of translational momentum
      parameter Boolean inclTransX=true "X" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));
      parameter Boolean inclTransY=true "Y" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));
      parameter Boolean inclTransZ=true "Z" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));
      //
      // Included faces
      parameter Boolean inclFacesX=true "X" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          tab="Assumptions",
          group="Axes with faces included",
          compact=true));
      parameter Boolean inclFacesY=true "Y" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          tab="Assumptions",
          group="Axes with faces included",
          compact=true));
      parameter Boolean inclFacesZ=true "Z" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          tab="Assumptions",
          group="Axes with faces included",
          compact=true));

      Connectors.FaceBus xPositive if inclFacesX
        "Positive face along the x axis" annotation (Placement(transformation(
              extent={{30,-10},{50,10}}), iconTransformation(extent={{90,-10},{
                110,10}})));
      Connectors.FaceBus xNegative if inclFacesX
        "Negative face along the x axis" annotation (Placement(transformation(
              extent={{-50,-10},{-30,10}}), iconTransformation(extent={{-110,-10},
                {-90,10}})));
      Connectors.FaceBus yPositive if inclFacesY
        "Positive face along the y axis" annotation (Placement(transformation(
              extent={{-10,30},{10,50}}), iconTransformation(extent={{-10,90},{
                10,110}})));
      Connectors.FaceBus yNegative if inclFacesY
        "Negative face along the y axis" annotation (Placement(transformation(
              extent={{-10,-50},{10,-30}}), iconTransformation(extent={{-10,-110},
                {10,-90}})));
      Connectors.FaceBus zPositive if inclFacesZ
        "Positive face along the z axis" annotation (Placement(transformation(
              extent={{-30,-30},{-10,-10}}), iconTransformation(extent={{-60,-60},
                {-40,-40}})));
      Connectors.FaceBus zNegative if inclFacesZ
        "Negative face along the z axis" annotation (Placement(transformation(
              extent={{10,10},{30,30}}), iconTransformation(extent={{40,40},{60,
                60}})));

    protected
      final inner parameter Q.Area A[Axis]={L[cartWrap(axis + 1)]*L[cartWrap(
          axis + 2)] for axis in Axis} "Cross-sectional areas";
      final inner parameter Boolean inclTrans[Axis]={inclTransX,inclTransY,
          inclTransZ}
        "true, if each component of translational momentum is included";
      final inner parameter Integer n_trans=countTrue(inclTrans)
        "Number of components of translational momentum";
      final inner parameter Integer cartAxes[n_trans]=index(inclTrans)
        "Cartesian-axis indices of the components of translational momentum";
      final inner parameter Integer transAxes[Axis]=enumerate(inclTrans)
        "Translational-momentum-component indices of the Cartesian axes";

      Volume volume "Model to establish a fixed total volume"
        annotation (Placement(transformation(extent={{-16,-16},{16,16}})));

      annotation (
        defaultComponentName="subregion",
        Documentation(info="<html>
  <p>This model must be be extended so that models can be added for
  relevant species, phases, and reactions.</p>

  <p>All of the components of translational momentum are included by default.  At least one component must be included.</p>
  </html>"),
        Icon(graphics={
            Line(
              points={{-100,0},{-40,0}},
              color={127,127,127},
              thickness=0.5,
              visible=inclFacesX,
              smooth=Smooth.None),
            Line(
              points={{0,-40},{0,-100}},
              color={127,127,127},
              thickness=0.5,
              visible=inclFacesY,
              smooth=Smooth.None),
            Line(
              points={{40,40},{50,50}},
              color={127,127,127},
              thickness=0.5,
              visible=inclFacesZ,
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
              visible=inclFacesX,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{0,28},{0,-40}},
              color={210,210,210},
              visible=inclFacesY,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{28,0},{100,0}},
              color={127,127,127},
              thickness=0.5,
              visible=inclFacesX,
              smooth=Smooth.None),
            Line(
              points={{0,100},{0,28}},
              color={127,127,127},
              thickness=0.5,
              visible=inclFacesY,
              smooth=Smooth.None),
            Line(
              points={{-12,-12},{40,40}},
              color={210,210,210},
              visible=inclFacesZ,
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
              visible=inclFacesZ,
              smooth=Smooth.None),
            Polygon(
              points={{-40,16},{-16,40},{40,40},{40,-16},{16,-40},{-40,-40},{-40,
                  16}},
              lineColor={127,127,127},
              smooth=Smooth.None),
            Line(
              points={{40,40},{16,16}},
              color={127,127,127},
              smooth=Smooth.None),
            Rectangle(
              extent={{-100,56},{100,96}},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Text(
              extent={{-100,56},{100,96}},
              textString="%name",
              lineColor={0,0,0})}));

    end PartialSubregion;

    type InitElectrical = enumeration(
        None "No initialization",
        SS "Steady state",
        Charge "Charge",
        Potential "Potential")
      "Methods of initializing an electrical capacitor";

  end BaseClasses;
  annotation (Documentation(info="<html>
<p><b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
Copyright 2007&ndash;2012, Georgia Tech Research Corporation.</p>

<p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.ModelicaLicense2\">
FCSys.UsersGuide.ModelicaLicense2</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p>
</html>"));
end Subregions;
