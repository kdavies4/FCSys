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
        "<html>Nafion sulfonate (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>)</html>"
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
        inclTransY=false,
        inclTransZ=false,
        inclFacesY=false,
        inclFacesZ=false,
        inclFacesX=true,
        graphite(
          final 'inclC+'='inclC+',
          final 'incle-'='incle-',
          'C+'(V_IC=subregion.V/4)),
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        liquid(H2O(V_IC=subregion.V/4)),
        ionomer(
          final 'inclC19HF37O5S-'='inclC19HF37O5S-',
          final 'inclH+'='inclH+',
          'C19HF37O5S-'(V_IC=subregion.V/4)))
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      inner FCSys.Conditions.Environment environment(analysis=true)
        annotation (Placement(transformation(extent={{40,40},{60,60}})));
      annotation (
        Placement(transformation(extent={{70,70},{90,90}})),
        experiment(StopTime=10),
        Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.Subregion.mos"));

    end Subregion;

    model SubregionCondensation
      "Test a subregion with evaporation/condensation"

      extends Examples.Subregion(
        inclH2O=true,
        inclH2=false,
        subregion(liquid(inclH2O=inclH2O)));

      annotation (
        experiment(StopTime=2e-07, Tolerance=1e-06),
        Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.SubregionCondensation.mos"),

        Diagram(graphics),
        experimentSetupOutput);

    end SubregionCondensation;

    model SubregionHydration
      "Test a subregion with absorption/desorption between gas and ionomer"

      extends Examples.Subregion(
        'inclC19HF37O5S-'=true,
        inclH2O=true,
        inclH2=false,
        subregion(gas(H2O(p_IC=1.001*U.atm)), ionomer(inclH2O=inclH2O)));
      // In Dymola 7.4, p_IC=1.1*environment.p has no effect on the
      // initial pressure, but p_IC=1.1*U.atm does.
      extends Modelica.Icons.UnderConstruction;
      // **Fix, update plotting script.

      annotation (
        experiment(StopTime=2e-07, Tolerance=1e-06),
        Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.SubregionHydration.mos"),

        Diagram(graphics),
        experimentSetupOutput);

    end SubregionHydration;

    model SubregionHOR
      "<html>Test a subregion with the hydrogen oxidation reaction and the essential species for it (C+, C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S, e<sup>-</sup>, H<sub>2</sub>, and H<sup>+</sup>)</html>"

      extends Examples.Subregion(
        'inclC+'=true,
        'inclC19HF37O5S-'=true,
        'incle-'=true,
        'inclH+'=true,
        inclH2=true,
        subregion(
          gas(H2(
              initMaterial=FCSys.Subregions.Species.BaseClasses.InitScalar.None,

              chemical(Ndot(start=0, fixed=true)),
              consEnergy=FCSys.Subregions.Species.BaseClasses.Conservation.IC)),

          graphite(
            reduceTemp=false,
            'C+'(consEnergy=FCSys.Subregions.Species.BaseClasses.Conservation.IC),

            'e-'(consEnergy=FCSys.Subregions.Species.BaseClasses.Conservation.IC,
                initTransX=FCSys.Subregions.Species.BaseClasses.InitTranslational.None)),

          ionomer(
            reduceTemp=false,
            'C19HF37O5S-'(consEnergy=FCSys.Subregions.Species.BaseClasses.Conservation.IC),

            'H+'(consEnergy=FCSys.Subregions.Species.BaseClasses.Conservation.IC,
                initTransX=FCSys.Subregions.Species.BaseClasses.InitTranslational.None))));

      Conditions.ByConnector.FaceBus.Single.FaceBusIsolated negativeBC(
        graphite(
          'inclC+'=false,
          final 'incle-'='incle-',
          'e-'(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Current
              material(source(final y=0)), redeclare
              Conditions.ByConnector.Face.Single.Translational.Current normal(
                source(y=time*U.A)))),
        ionomer(
          'inclC19HF37O5S-'=false,
          final 'inclH+'='inclH+',
          'H+'(redeclare Conditions.ByConnector.Face.Single.Translational.Force
              normal)),
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Density
              material(source(y=environment.p/environment.T))))) annotation (
          Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={-24,0})));

      Conditions.ByConnector.FaceBus.Single.FaceBusIsolated positiveBC(
        ionomer(
          'inclC19HF37O5S-'=false,
          final 'inclH+'='inclH+',
          'H+'(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Density
              material(source(final y=17842.7*U.C/U.cc)), redeclare
              Conditions.ByConnector.Face.Single.Translational.Velocity normal)),

        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        graphite(
          'inclC+'=false,
          final 'incle-'='incle-',
          'e-'(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Translational.Force
              normal))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={24,0})));

      extends Modelica.Icons.UnderConstruction;

    equation
      connect(negativeBC.face, subregion.xNegative) annotation (Line(
          points={{-20,-1.34539e-15},{-16,-1.34539e-15},{-16,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(subregion.xPositive, positiveBC.face) annotation (Line(
          points={{10,6.10623e-16},{14,6.10623e-16},{14,1.23436e-15},{20,
              1.23436e-15}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        experiment(StopTime=150, Tolerance=1e-08),
        Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.SubregionHOR.mos"),

        experimentSetupOutput);

    end SubregionHOR;

    model SubregionORR
      "<html>Test a subregion with the oxygen reduction reaction and the essential species for it (C+, C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S, e<sup>-</sup>, H<sub>2</sub>O, H<sup>+</sup>, and O<sub>2</sub>)</html>"
      // **Update plotting script.

      extends Examples.SubregionHOR(
        inclH2=false,
        inclH2O=true,
        inclO2=true);
      // **Copy rather than extend.
      annotation (experiment(StopTime=1000, Tolerance=1e-06), Commands(file(
              ensureSimulated=true) =
            "resources/scripts/Dymola/Subregions.Examples.SubregionORR.mos"));

    end SubregionORR;

    model SubregionPipeFlow
      import assert = FCSys.BaseClasses.Utilities.assertEval;
      extends Examples.Subregion(subregion(
          L={100,1,1}*U.mm,
          inclFacesY=true,
          inclFacesZ=true));

      parameter Q.Pressure Deltap=-10*U.Pa
        "<html>Prescribed gas pressure difference (&Delta;<i>p</i>)</html>";

      output Q.Length D=2*subregion.A[Axis.x]/(subregion.L[Axis.y] + subregion.L[
          Axis.z]);
      output Q.Number Re=subregion.gas.H2.phi[Axis.x]*D*subregion.gas.H2.zeta*
          subregion.gas.H2.Data.m*subregion.gas.H2.rho "Reynolds number";
      output Q.Pressure Deltap_Poiseuille=-32*subregion.L[Axis.x]*subregion.gas.H2.phi[
          Axis.x]/(D^2*subregion.gas.H2.zeta)
        "Pressure difference according to Poiseuille's law";

      Conditions.ByConnector.FaceBus.Single.FaceBusIsolated BC1(
        graphite('inclC+'=false, final 'incle-'='incle-'),
        ionomer('inclC19HF37O5S-'=false, final 'inclH+'='inclH+'),
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Translational.Force
              normal, redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Pressure
              material(source(y=environment.p - Deltap/2))),
          H2O(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Translational.Force
              normal, redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Pressure
              material(source(y=environment.p - Deltap/2))),
          N2(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Translational.Force
              normal, redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Pressure
              material(source(y=environment.p - Deltap/2))),
          O2(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Translational.Force
              normal, redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Pressure
              material(source(y=environment.p - Deltap/2))))) annotation (
          Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={-24,0})));

      Conditions.ByConnector.FaceBus.Single.FaceBusIsolated BC2(
        graphite('inclC+'=false, final 'incle-'='incle-'),
        ionomer('inclC19HF37O5S-'=false, final 'inclH+'='inclH+'),
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Translational.Force
              normal, redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Pressure
              material(source(y=environment.p + Deltap/2))),
          H2O(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Translational.Force
              normal, redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Pressure
              material(source(y=environment.p + Deltap/2))),
          N2(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Translational.Force
              normal, redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Pressure
              material(source(y=environment.p + Deltap/2))),
          O2(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Translational.Force
              normal, redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Pressure
              material(source(y=environment.p + Deltap/2))))) annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={24,0})));

      Conditions.ByConnector.FaceBus.Single.FaceBus BC3(
        graphite('inclC+'=false, final 'incle-'='incle-'),
        ionomer('inclC19HF37O5S-'=false, final 'inclH+'='inclH+'),
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Current
              material),
          H2O(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Current
              material),
          N2(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Current
              material),
          O2(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Current
              material))) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=0,
            origin={0,-24})));

      Conditions.ByConnector.FaceBus.Single.FaceBus BC4(
        graphite('inclC+'=false, final 'incle-'='incle-'),
        ionomer('inclC19HF37O5S-'=false, final 'inclH+'='inclH+'),
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Current
              material),
          H2O(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Current
              material),
          N2(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Current
              material),
          O2(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Current
              material))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={0,24})));

      Conditions.ByConnector.FaceBus.Single.FaceBus BC5(
        graphite('inclC+'=false, final 'incle-'='incle-'),
        ionomer('inclC19HF37O5S-'=false, final 'inclH+'='inclH+'),
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Current
              material),
          H2O(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Current
              material),
          N2(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Current
              material),
          O2(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Current
              material))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=315,
            origin={24,24})));
      Conditions.ByConnector.FaceBus.Single.FaceBus BC6(
        graphite('inclC+'=false, final 'incle-'='incle-'),
        ionomer('inclC19HF37O5S-'=false, final 'inclH+'='inclH+'),
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Current
              material),
          H2O(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Current
              material),
          N2(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Current
              material),
          O2(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Current
              material))) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=315,
            origin={-24,-24})));
    equation
      assert(abs(Deltap - (subregion.gas.H2.p_faces[1, Side.p] - subregion.gas.H2.p_faces[
        1, Side.n])) <= Modelica.Constants.small,
        "The actual pressure difference and the prescribed pressure difference are not equal.");

      connect(BC1.face, subregion.xNegative) annotation (Line(
          points={{-20,-1.34539e-15},{-16,-1.34539e-15},{-16,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(BC2.face, subregion.xPositive) annotation (Line(
          points={{20,1.23436e-15},{16,1.23436e-15},{16,6.10623e-16},{10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(subregion.yNegative, BC3.face) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-12.5},{6.10623e-16,-12.5},{
              6.10623e-16,-15},{6.10623e-16,-20},{6.10623e-16,-20}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(BC4.face, subregion.yPositive) annotation (Line(
          points={{6.10623e-16,20},{-3.36456e-22,16},{6.10623e-16,16},{
              6.10623e-16,10}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(BC5.face, subregion.zPositive) annotation (Line(
          points={{21.1716,21.1716},{46,21.1716},{46,-5},{-5,-5}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(BC6.face, subregion.zNegative) annotation (Line(
          points={{-21.1716,-21.1716},{52,-21.1716},{52,5},{5,5}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        Diagram(graphics),
        experiment(StopTime=200, Algorithm="Dassl"),
        experimentSetupOutput,
        Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.SubregionPipeFlow.mos"));
    end SubregionPipeFlow;

    model Subregions
      "<html>Test a one-dimensional array of subregions with an initial pressure difference (C<sup>+</sup> and H<sub>2</sub> included by default)</html>"
      extends Modelica.Icons.Example;

      parameter Integer n_x=0
        "<html>Number of discrete subregions along the x axis, besides the 2 side subregions (<i>n<i><sub>x</sub>)</html>";
      parameter Q.Pressure Deltap_IC=100*U.Pa
        "<html>Initial pressure difference (&Delta;<i>p</i><sub>IC</sub>)</html>";
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

      FCSys.Subregions.Subregion subregion1(
        L={1,1,1}*U.cm,
        inclFacesY=false,
        inclFacesZ=false,
        inclTransY=false,
        inclTransZ=false,
        gas(
          reduceVel=true,
          reduceTemp=true,
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2O(p_IC=environment.p + Deltap_IC/2),
          N2(p_IC=environment.p + Deltap_IC/2),
          O2(p_IC=environment.p + Deltap_IC/2),
          H2(p_IC=environment.p + Deltap_IC/2)),
        graphite(
          reduceTemp=true,
          final 'inclC+'='inclC+',
          final 'incle-'='incle-',
          'C+'(V_IC=subregion1.V/1000)),
        ionomer(
          reduceTemp=true,
          final 'inclC19HF37O5S-'='inclC19HF37O5S-',
          final 'inclH+'='inclH+',
          'C19HF37O5S-'(V_IC=subregion1.V/1000)),
        liquid(H2O(V_IC=subregion1.V/4)))
        annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));

      FCSys.Subregions.Subregion subregions[n_x](
        each L={1,1,1}*U.cm,
        each inclTransY=false,
        each inclTransZ=false,
        each inclFacesY=false,
        each inclFacesZ=false,
        graphite(
          each reduceTemp=true,
          each final 'inclC+'='inclC+',
          each final 'incle-'='incle-',
          'C+'(each V_IC=subregions[1].V/1000)),
        ionomer(
          each reduceTemp=true,
          each final 'inclC19HF37O5S-'='inclC19HF37O5S-',
          each final 'inclH+'='inclH+',
          'C19HF37O5S-'(each V_IC=subregions[1].V/1000)),
        gas(
          each reduceTemp=true,
          each reduceVel=true,
          each final inclH2=inclH2,
          each final inclH2O=inclH2O,
          each final inclN2=inclN2,
          each final inclO2=inclO2,
          H2(p_IC={environment.p + Deltap_IC/2 - i*Deltap_IC/(n_x + 1) for i
                 in 1:n_x}),
          H2O(p_IC=subregions.gas.H2.p_IC),
          N2(p_IC=subregions.gas.H2.p_IC),
          O2(p_IC=subregions.gas.H2.p_IC)),
        liquid(H2O(each V_IC=subregions[1].V/4))) if n_x > 0
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      FCSys.Subregions.Subregion subregion2(
        L={1,1,1}*U.cm,
        inclFacesY=false,
        inclFacesZ=false,
        inclTransY=false,
        inclTransZ=false,
        gas(
          reduceVel=true,
          reduceTemp=true,
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2O(p_IC=environment.p - Deltap_IC/2),
          N2(p_IC=environment.p - Deltap_IC/2),
          O2(p_IC=environment.p - Deltap_IC/2),
          H2(p_IC=environment.p - Deltap_IC/2)),
        graphite(
          reduceTemp=true,
          final 'inclC+'='inclC+',
          final 'incle-'='incle-',
          'C+'(V_IC=subregion2.V/1000)),
        ionomer(
          reduceTemp=true,
          final 'inclC19HF37O5S-'='inclC19HF37O5S-',
          final 'inclH+'='inclH+',
          'C19HF37O5S-'(V_IC=subregion2.V/1000)))
        annotation (Placement(transformation(extent={{20,-10},{40,10}})));

      inner FCSys.Conditions.Environment environment(analysis=true)
        annotation (Placement(transformation(extent={{60,20},{80,40}})));

      Conditions.ByConnector.FaceBus.Single.FaceBusIsolated BC1(
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        graphite('inclC+'=false, final 'incle-'='incle-'),
        ionomer('inclC19HF37O5S-'=false, final 'inclH+'='inclH+')) annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-56,0})));

      Conditions.ByConnector.FaceBus.Single.FaceBusIsolated BC2(
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        graphite('inclC+'=false, final 'incle-'='incle-'),
        ionomer('inclC19HF37O5S-'=false, final 'inclH+'='inclH+')) annotation (
          Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={56,0})));

    equation
      connect(subregion1.xPositive, subregions[1].xNegative) annotation (Line(
          points={{-20,6.10623e-16},{-16,-3.36456e-22},{-16,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      for i in 1:n_x - 1 loop
        connect(subregions[i].xPositive, subregions[i + 1].xNegative)
          "Not shown in the diagram";
      end for;
      if n_x == 0 then
        connect(subregion1.xPositive, subregion2.xNegative)
          "Not shown in the diagram";
      end if;
      connect(subregions[n_x].xPositive, subregion2.xNegative) annotation (Line(
          points={{10,6.10623e-16},{20,-3.36456e-22},{20,6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(BC1.face, subregion1.xNegative) annotation (Line(
          points={{-52,3.65701e-16},{-46,3.65701e-16},{-46,6.10623e-16},{-40,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(subregion2.xPositive, BC2.face) annotation (Line(
          points={{40,6.10623e-16},{46,6.10623e-16},{46,-2.54679e-16},{52,
              -2.54679e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        Placement(transformation(extent={{70,70},{90,90}})),
        experiment(
          StopTime=2,
          Tolerance=1e-06,
          Algorithm="Dassl"),
        Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.Subregions.mos"),
        experimentSetupOutput);

    end Subregions;

    model ElectricalConduction
      "<html>Test a one-dimensional array of subregions with C<sup>+</sup> and e<sup>-</sup></html>"
      // **Make this into an example of Ohm's law
      extends Examples.Subregion(
        'inclC+'=true,
        'incle-'=true,
        inclH2=false,
        subregion(graphite(
            'C+'(consMaterial=Conservation.IC, initMaterial=InitScalar.Pressure),

            reduceTemp=true,
            'e-'(p_IC=4560*U.bar, initTransX=FCSys.Subregions.Species.BaseClasses.InitTranslational.None))));

      output Q.Potential w=(subregion.graphite.'e-'.faces[1, Side.p].mPhidot[
          Orientation.normal]/subregion.graphite.'e-'.faces[1, Side.p].rho -
          subregion.graphite.'e-'.faces[1, Side.n].mPhidot[Orientation.normal]/
          subregion.graphite.'e-'.faces[1, Side.n].rho)/subregion.A[Axis.x]
        "Potential";
      output Q.Current zI=-subregion.graphite.'e-'.Ndot_faces[1, Side.n]
        "Electrical current";
      output Q.ResistanceElectrical R=w/zI "Electrical resistance";
      output Q.ResistanceElectrical R2=subregion.graphite.'e-'.v/subregion.graphite.
          'e-'.mu/U.cm "**Electrical resistance";

      Conditions.ByConnector.FaceBus.Single.FaceBusIsolated BC1(
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        ionomer('inclC19HF37O5S-'=false, final 'inclH+'='inclH+'),
        graphite(
          'inclC+'=false,
          final 'incle-'='incle-',
          'e-'(redeclare
              FCSys.Conditions.ByConnector.Face.Single.Material.Density
              material(source(final y=17842.7*U.C/U.cc)), redeclare
              FCSys.Conditions.ByConnector.Face.Single.Translational.Current
              normal(source(y=U.A))))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-24,0})));
      Conditions.ByConnector.FaceBus.Single.FaceBusIsolated BC2(
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        ionomer('inclC19HF37O5S-'=false, final 'inclH+'='inclH+'),
        graphite(
          'inclC+'=false,
          final 'incle-'='incle-',
          'e-'(material(source(final y=0)),redeclare
              FCSys.Conditions.ByConnector.Face.Single.Translational.Force
              normal))) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={24,0})));
      annotation (
        experiment(StopTime=500, Tolerance=1e-06),
        Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.ElectricalConduction.mos"),

        experimentSetupOutput,
        Diagram(graphics));

    equation
      connect(BC1.face, subregion.xNegative) annotation (Line(
          points={{-20,3.65701e-16},{-16,3.65701e-16},{-16,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(subregion.xPositive, BC2.face) annotation (Line(
          points={{10,6.10623e-16},{16,6.10623e-16},{16,-2.54679e-16},{20,
              -2.54679e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

    end ElectricalConduction;


    model ThermalConduction "Test thermal conduction (through solid)"
      extends Examples.Subregions(
        n_x=8,
        'inclC+'=true,
        inclH2=false,
        subregion1(graphite(T_IC=1.1*environment.T,'C+'(
              consMaterial=Conservation.IC,
              initMaterial=InitScalar.Pressure,
              T(displayUnit="degC")))),
        subregions(graphite('C+'(
              each consMaterial=Conservation.IC,
              each initMaterial=InitScalar.Pressure,
              each T(displayUnit="degC")))),
        subregion2(graphite('C+'(
              consMaterial=Conservation.IC,
              initMaterial=InitScalar.Pressure,
              T(displayUnit="degC")))));

      annotation (
        Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.ThermalConduction.mos"),

        experiment(StopTime=200, Algorithm="Dassl"),
        experimentSetupOutput);

    end ThermalConduction;

    model ThermalConductionConvection
      "Test combined thermal conduction and convection"

      extends Examples.Subregions(
        n_x=8,
        'inclC+'=true,
        inclH2=false,
        inclN2=true,
        subregion1(gas(T_IC=1.1*environment.T,N2(
              p_IC=environment.p,
              T(displayUnit="degC"),
              phi(displayUnit="mm/s"))), graphite(T_IC=1.1*environment.T,'C+'(
                V_IC=0.5*subregion1.V, T(displayUnit="degC")))),
        subregions(gas(N2(
              each p_IC=environment.p,
              each T(displayUnit="degC"),
              phi(each displayUnit="mm/s"))), graphite('C+'(each V_IC=0.5*
                  subregions[1].V, each T(displayUnit="degC")))),
        subregion2(gas(N2(
              p_IC=environment.p,
              T(displayUnit="degC"),
              phi(displayUnit="mm/s"))), graphite('C+'(V_IC=0.5*subregion2.V, T(
                  displayUnit="degC")))));

      annotation (
        Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.ThermalConductionConvection.mos"),

        experiment(
          StopTime=140,
          Tolerance=1e-06,
          Algorithm="Dassl"),
        experimentSetupOutput);

    end ThermalConductionConvection;

    model TestReaction "Test an electrochemical reaction"
      extends Modelica.Icons.Example;
      import FCSys.BaseClasses.Utilities.cartWrap;
      import FCSys.BaseClasses.Utilities.countTrue;
      import FCSys.BaseClasses.Utilities.enumerate;
      import FCSys.BaseClasses.Utilities.index;

      // Assumptions
      // -----------
      // Included components of translational momentum
      parameter Boolean inclTransX=true "X" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));
      parameter Boolean inclTransY=false "Y" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));
      parameter Boolean inclTransZ=false "Z" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));

      Reaction reaction
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      inner Conditions.Environment environment(analysis=true)
        annotation (Placement(transformation(extent={{70,72},{90,92}})));

      ChemicalExchange chemicalExchange[3](n={-2,-2,1}, m={1,1,1}*U.g/U.mol)
        annotation (Placement(transformation(extent={{-30,-30},{-10,-10}})));
      FCSys.Subregions.Depletion layer1(
        transSubstrate=true,
        thermalSubstrate=true,
        side=FCSys.BaseClasses.Side.n,
        redeclare package Data = FCSys.Characteristics.'e-'.Gas)
        annotation (Placement(transformation(extent={{-30,10},{-10,30}})));
      FCSys.Subregions.Depletion layer2(side=FCSys.BaseClasses.Side.p,
          redeclare package Data = FCSys.Characteristics.'H+'.Gas)
        annotation (Placement(transformation(extent={{30,10},{10,30}})));
      replaceable Conditions.ByConnector.ChemicalSpecies.Current current(
        redeclare Modelica.Blocks.Sources.Ramp source(height=100*U.A, duration=
              100),
        final inclTransX=inclTransX,
        final inclTransY=inclTransY,
        final inclTransZ=inclTransZ,
        sT=3000*U.K)
        annotation (Placement(transformation(extent={{-60,-60},{-40,-40}})));

      replaceable Conditions.ByConnector.ChemicalSpecies.Potential potential[2]
        (
        each sT=0,
        each final inclTransX=inclTransX,
        each final inclTransY=inclTransY,
        each final inclTransZ=inclTransZ)
        annotation (Placement(transformation(extent={{-60,-30},{-40,-10}})));
      Conditions.ByConnector.Inert.Inert catalyst(
        thermal(source(y=360*U.K)),
        final inclTransX=inclTransX,
        final inclTransY=inclTransY,
        final inclTransZ=inclTransZ)
        annotation (Placement(transformation(extent={{-10,30},{10,50}})));

      Conditions.ByConnector.Face.Single.Face face1 annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-44,20})));
      Conditions.ByConnector.Face.Single.Face face2 annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={44,20})));

    protected
      final inner parameter Boolean inclTrans[Axis]={inclTransX,inclTransY,
          inclTransZ}
        "true, if each component of translational momentum is included";
      final inner parameter Integer n_trans=countTrue(inclTrans)
        "Number of components of translational momentum";
      final inner parameter Integer cartTrans[:]=index(inclTrans)
        "Cartesian-axis indices of the components of translational momentum";

    equation
      connect(chemicalExchange[1].reaction, reaction.chemical) annotation (Line(
          points={{-18,-20},{6.10623e-16,-20},{6.10623e-16,6.10623e-16}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(chemicalExchange[2].reaction, reaction.chemical) annotation (Line(
          points={{-18,-20},{6.10623e-16,-20},{6.10623e-16,6.10623e-16}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(chemicalExchange[3].reaction, reaction.chemical) annotation (Line(
          points={{-18,-20},{6.10623e-16,-20},{6.10623e-16,6.10623e-16}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(layer2.chemical, reaction.chemical) annotation (Line(
          points={{20,20},{20,6.10623e-16},{6.10623e-16,6.10623e-16}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(layer1.chemical, reaction.chemical) annotation (Line(
          points={{-20,20},{-20,6.10623e-16},{6.10623e-16,6.10623e-16}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(current.chemical, chemicalExchange[3].species) annotation (Line(
          points={{-50,-54},{-50,-60},{-30,-60},{-30,-20},{-22,-20}},
          color={255,195,38},
          smooth=Smooth.None));

      connect(potential.chemical, chemicalExchange[1:2].species) annotation (
          Line(
          points={{-50,-24},{-50,-32},{-30,-32},{-30,-20},{-22,-20}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(catalyst.inert, layer1.inert) annotation (Line(
          points={{6.10623e-16,36},{0,20},{-10,20}},
          color={47,107,251},
          smooth=Smooth.None));
      connect(catalyst.inert, layer2.inert) annotation (Line(
          points={{6.10623e-16,36},{0,20},{10,20}},
          color={47,107,251},
          smooth=Smooth.None));
      connect(face1.face, layer1.face) annotation (Line(
          points={{-40,20},{-30,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(layer2.face, face2.face) annotation (Line(
          points={{30,20},{40,20}},
          color={127,127,127},
          smooth=Smooth.None));
      annotation (
        experiment(StopTime=120),
        Commands(file=
              "resources/scripts/Dymola/Subregions.Examples.TestReaction.mos"),

        Diagram(graphics));

    end TestReaction;

  end Examples;

  model Subregion "Subregion with all phases"
    extends FCSys.Subregions.BaseClasses.EmptySubregion(final hasSpecies=gas.n_spec
           + graphite.n_spec + ionomer.n_spec + liquid.n_spec > 0, volume(
          n_phases=4));

    Phases.Gas gas(
      final n_faces=n_faces,
      final inclHOR=inclHOR,
      final inclORR=inclORR) "Gas" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-10,-10},
              {10,10}})));

    Phases.Graphite graphite(
      final n_faces=n_faces,
      final inclHOR=inclHOR,
      final inclORR=inclORR) "Graphite" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-10,-10},
              {10,10}})));

    Phases.Ionomer ionomer(
      final n_faces=n_faces,
      final inclHOR=inclHOR,
      final inclORR=inclORR) "Ionomer" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-10,-10},
              {10,10}})));

    Phases.Liquid liquid(final n_faces=n_faces) "Liquid" annotation (Dialog(
          group="Phases (click to edit)"), Placement(transformation(extent={{-10,
              -10},{10,10}})));

    parameter Boolean inclHOR=graphite.'incle-' and ionomer.'inclH+' and gas.inclH2
      "Hydrogen oxidation" annotation (Dialog(group="Included reactions",
          compact=true), choices(__Dymola_checkBox=true));
    parameter Boolean inclORR=graphite.'incle-' and ionomer.'inclH+' and gas.inclO2
         and gas.inclH2O "Oxygen reduction" annotation (Dialog(group=
            "Included reactions", compact=true), choices(__Dymola_checkBox=true));

    Reaction reaction if inclHOR or inclORR "Electrochemical reaction"
      annotation (Placement(transformation(extent={{-23,17},{-3,37}})));

  protected
    Connectors.PhysicalBusInternal physical if gas.inclH2O or ionomer.inclH2O
       or liquid.inclH2O "Connector for phase change" annotation (Placement(
          transformation(extent={{-37,3},{-17,23}}), iconTransformation(extent=
              {{-64,36},{-44,56}})));

  equation
    // Phase change
    connect(gas.physical, physical) annotation (Line(
        points={{-5.8,4},{-27,13}},
        color={239,142,1},
        smooth=Smooth.None));
    connect(ionomer.physical, physical) annotation (Line(
        points={{-5.8,4},{-27,13}},
        color={239,142,1},
        smooth=Smooth.None));
    connect(liquid.physical, physical) annotation (Line(
        points={{-5.8,4},{-27,13}},
        color={38,196,52},
        smooth=Smooth.None,
        thickness=0.5));

    // Chemical reactions
    connect(reaction.chemical, gas.chemical) annotation (Line(
        points={{-13,27},{-3.8,7.6}},
        color={255,195,38},
        smooth=Smooth.None));
    connect(reaction.chemical, graphite.chemical) annotation (Line(
        points={{-13,27},{-3.8,7.6}},
        color={255,195,38},
        smooth=Smooth.None));
    connect(reaction.chemical, ionomer.chemical) annotation (Line(
        points={{-13,27},{-3.8,7.6}},
        color={255,195,38},
        smooth=Smooth.None));

    // Phases
    // ------
    // Gas
    connect(gas.inertDalton, volume.inertDalton[1]) annotation (Line(
        points={{8,-8},{11,-11}},
        color={47,107,251},
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
        points={{6.10623e-16,10},{0,10},{0,10},{-4.87687e-22,40},{-4.87687e-22,
            40},{5.55112e-16,40}},
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
    connect(graphite.inertDalton, volume.inertDalton[2]) annotation (Line(
        points={{8,-8},{11,-11}},
        color={47,107,251},
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
        points={{6.10623e-16,10},{0,10},{0,10},{-4.87687e-22,40},{-4.87687e-22,
            40},{5.55112e-16,40}},
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
    connect(ionomer.inertDalton, volume.inertDalton[3]) annotation (Line(
        points={{8,-8},{11,-11}},
        color={47,107,251},
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
    connect(liquid.inertDalton, volume.inertDalton[4]) annotation (Line(
        points={{8,-8},{11,-11}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(liquid.xNegative, xNegative.liquid) annotation (Line(
        points={{-8,6.10623e-16},{-24,0},{-40,5.55112e-16}},
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
        points={{6.10623e-16,10},{0,10},{-4.87687e-22,40},{5.55112e-16,40}},
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
    annotation (Documentation(info="<html>
<p>Please see the documentation of the
   <a href=\"modelica://FCSys.Subregions.BaseClasses.EmptySubregion\">EmptySubregion</a> model.</p></html>"));
  end Subregion;

  model SubregionIonomerOnly "Subregion with only the ionomer phase"
    extends FCSys.Subregions.BaseClasses.EmptySubregion(final hasSpecies=
          ionomer.n_spec > 0, volume(n_phases=1));

    Phases.Ionomer ionomer(
      final n_faces=n_faces,
      inclH2O=true,
      final inclHOR=false,
      final inclORR=false) "Ionomer" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-10,-10},
              {10,10}})));

  equation
    // Ionomer
    connect(ionomer.inertDalton, volume.inertDalton[1]) annotation (Line(
        points={{8,-8},{11,-11}},
        color={47,107,251},
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
      defaultComponentName="subregion",
      Documentation(info="<html>
<p>Please see the documentation of the
   <a href=\"modelica://FCSys.Subregions.BaseClasses.EmptySubregion\">EmptySubregion</a> model.</p></html>"),

      Diagram(graphics));
  end SubregionIonomerOnly;

  model SubregionNoIonomer "Subregion with all phases except ionomer"
    extends FCSys.Subregions.BaseClasses.EmptySubregion(final hasSpecies=gas.n_spec
           + graphite.n_spec + liquid.n_spec > 0, volume(n_phases=3));

    Phases.Gas gas(
      final n_faces=n_faces,
      final inclHOR=false,
      final inclORR=false) "Gas" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-10,-10},
              {10,10}})));
    Phases.Graphite graphite(
      final n_faces=n_faces,
      final inclHOR=false,
      final inclORR=false) "Graphite" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-10,-10},
              {10,10}})));
    Phases.Liquid liquid(final n_faces=n_faces) "Liquid" annotation (Dialog(
          group="Phases (click to edit)"), Placement(transformation(extent={{-10,
              -10},{10,10}})));

  protected
    Connectors.PhysicalBusInternal physical "Connector for phase change"
      annotation (Placement(transformation(extent={{-30,10},{-10,30}}),
          iconTransformation(extent={{-76,16},{-56,36}})));

  equation
    // Phase change
    connect(physical, gas.physical);
    connect(physical, liquid.physical) annotation (Line(
        points={{-20,20},{-5.8,4}},
        color={38,196,52},
        smooth=Smooth.None,
        thickness=0.5));

    // Phases
    // ------
    // Gas
    connect(gas.inertDalton, volume.inertDalton[1]) annotation (Line(
        points={{8,-8},{11,-11}},
        color={47,107,251},
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
    connect(graphite.inertDalton, volume.inertDalton[2]) annotation (Line(
        points={{8,-8},{11,-11}},
        color={47,107,251},
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
    connect(liquid.inertDalton, volume.inertDalton[3]) annotation (Line(
        points={{8,-8},{11,-11}},
        color={47,107,251},
        smooth=Smooth.None));
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
    annotation (defaultComponentName="subregion", Documentation(info="<html>
<p>Please see the documentation of the
   <a href=\"modelica://FCSys.Subregions.BaseClasses.EmptySubregion\">EmptySubregion</a> model.</p></html>"));
  end SubregionNoIonomer;

  package Phases "Phases or mixtures of species"
    extends Modelica.Icons.Package;

    model Gas "Gas phase"
      import FCSys.BaseClasses.Utilities.countTrue;
      extends FCSys.Subregions.Phases.BaseClasses.EmptyPhase(final n_spec=
            countTrue({inclH2,inclH2O,inclN2,inclO2}));

      // Conditionally include species.
      parameter Boolean inclH2=false "<html>Hydrogen (H<sub>2</sub>)</html>"
        annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));

      replaceable Species.H2.Gas.Fixed H2(final n_faces) if inclH2
        constrainedby Species.Species(
        n_faces=n_faces,
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

      parameter Boolean inclH2O=true "<html>Water (H<sub>2</sub>O)</html>"
        annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.H2O.Gas.Fixed H2O(final n_faces) if inclH2O
        constrainedby Species.Species(
        n_faces=n_faces,
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
      replaceable Species.N2.Gas.Fixed N2(final n_faces) if inclN2
        constrainedby Species.Species(
        n_faces=n_faces,
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
      replaceable Species.O2.Gas.Fixed O2(final n_faces) if inclO2
        constrainedby Species.Species(
        n_faces=n_faces,
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

      ChemicalExchange HOR(n=1, final m=H2.Data.m) if inclHOR
        "Interface with the hydrogen oxidation reaction"
        annotation (Placement(transformation(extent={{-10,40},{-30,60}})));
      ChemicalExchange ORR[2](n={-1,2}, final m={O2.Data.m,H2O.Data.m}) if
        inclORR "Interfaces with the oxygen reduction reaction"
        annotation (Placement(transformation(extent={{-10,20},{-30,40}})));
      Connectors.ChemicalReaction chemical if inclHOR or inclORR
        "Connector for chemical reactions" annotation (Placement(transformation(
              extent={{-50,40},{-30,60}}), iconTransformation(extent={{-48,66},
                {-28,86}})));
      Connectors.PhysicalBus physical if inclH2O "Connector for phase change"
        annotation (Placement(transformation(extent={{-37,3},{-17,23}}),
            iconTransformation(extent={{-68,30},{-48,50}})));

      parameter Boolean inclHOR=false "Hydrogen oxidation" annotation (
        HideResult=true,
        Dialog(group="Included reactions", compact=true),
        choices(__Dymola_checkBox=true));
      parameter Boolean inclORR=false "Oxygen reduction" annotation (
        HideResult=true,
        Dialog(group="Included reactions", compact=true),
        choices(__Dymola_checkBox=true));
      // These can't be outer parameters in Dymola 7.4.

    equation
      // Phase change
      connect(physical.H2O, H2O.physical) annotation (Line(
          points={{-27,13},{-6.9,3.8}},
          smooth=Smooth.None,
          color={38,196,52}));

      // Reactions
      // ---------
      // HOR
      connect(H2.chemical, HOR.species) annotation (Line(
          points={{-3.9,7},{-18,50}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(HOR.reaction, chemical) annotation (Line(
          points={{-22,50},{-40,50}},
          color={255,195,38},
          smooth=Smooth.None));
      // ORR
      connect(O2.chemical, ORR[1].species) annotation (Line(
          points={{-3.9,7},{-18,30}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(H2O.chemical, ORR[2].species) annotation (Line(
          points={{-3.9,7},{-18,30}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(ORR[1].reaction, chemical) annotation (Line(
          points={{-22,30},{-30,30},{-30,50},{-40,50}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(ORR[2].reaction, chemical) annotation (Line(
          points={{-22,30},{-30,30},{-30,50},{-40,50}},
          color={255,195,38},
          smooth=Smooth.None));

      // H2
      // --
      // Diffusive exchange
      connect(H2.inert.translational, common.translational) annotation (Line(
          points={{6.9,-3.8},{26.67,-13.33}},
          color={11,43,197},
          smooth=Smooth.None));
      connect(H2.inert.thermal, common.thermal) annotation (Line(
          points={{6.9,-3.8},{26.67,-13.33}},
          color={11,43,197},
          smooth=Smooth.None));
      connect(H2.inertDalton, inertDalton) annotation (Line(
          points={{3.9,-7},{13.33,-26.67}},
          color={11,43,197},
          smooth=Smooth.None));
      // Transport
      connect(H2.faces[facesCart[Axis.x], Side.n], xNegative.H2) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.faces[facesCart[Axis.x], Side.p], xPositive.H2) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},
              {40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.faces[facesCart[Axis.y], Side.n], yNegative.H2) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{6.10623e-16,-40},{5.55112e-16,-40}},

          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.faces[facesCart[Axis.y], Side.p], yPositive.H2) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2.faces[facesCart[Axis.z], Side.n], zNegative.H2) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2.faces[facesCart[Axis.z], Side.p], zPositive.H2) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // H2O
      // ---
      // Diffusive exchange
      connect(H2O.inert.translational, common.translational) annotation (Line(
          points={{6.9,-3.8},{26.67,-13.33}},
          color={11,43,197},
          smooth=Smooth.None));
      connect(H2O.inert.thermal, common.thermal) annotation (Line(
          points={{6.9,-3.8},{26.67,-13.33}},
          color={11,43,197},
          smooth=Smooth.None));
      connect(H2O.inertDalton, inertDalton) annotation (Line(
          points={{3.9,-7},{13.33,-26.67}},
          color={11,43,197},
          smooth=Smooth.None));
      // Transport
      connect(H2O.faces[facesCart[Axis.x], Side.n], xNegative.H2O) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[facesCart[Axis.x], Side.p], xPositive.H2O) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},
              {40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[facesCart[Axis.y], Side.n], yNegative.H2O) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{6.10623e-16,-40},{5.55112e-16,-40}},

          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[facesCart[Axis.y], Side.p], yPositive.H2O) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[facesCart[Axis.z], Side.n], zNegative.H2O) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.faces[facesCart[Axis.z], Side.p], zPositive.H2O) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // N2
      // --
      // Diffusive exchange
      connect(N2.inert.translational, common.translational) annotation (Line(
          points={{6.9,-3.8},{26.67,-13.33}},
          color={11,43,197},
          smooth=Smooth.None));
      connect(N2.inert.thermal, common.thermal) annotation (Line(
          points={{6.9,-3.8},{26.67,-13.33}},
          color={11,43,197},
          smooth=Smooth.None));
      connect(N2.inertDalton, inertDalton) annotation (Line(
          points={{3.9,-7},{13.33,-26.67}},
          color={11,43,197},
          smooth=Smooth.None));

      // Transport
      connect(N2.faces[facesCart[Axis.x], Side.n], xNegative.N2) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.faces[facesCart[Axis.x], Side.p], xPositive.N2) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},
              {40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.faces[facesCart[Axis.y], Side.n], yNegative.N2) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{6.10623e-16,-40},{5.55112e-16,-40}},

          color={127,127,127},
          smooth=Smooth.None));

      connect(N2.faces[facesCart[Axis.y], Side.p], yPositive.N2) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.faces[facesCart[Axis.z], Side.n], zNegative.N2) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(N2.faces[facesCart[Axis.z], Side.p], zPositive.N2) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      // O2
      // --
      // Diffusive exchange
      connect(O2.inert.translational, common.translational) annotation (Line(
          points={{6.9,-3.8},{26.67,-13.33}},
          color={11,43,197},
          smooth=Smooth.None));
      connect(O2.inert.thermal, common.thermal) annotation (Line(
          points={{6.9,-3.8},{26.67,-13.33}},
          color={47,107,251},
          smooth=Smooth.None));
      connect(O2.inertDalton, inertDalton) annotation (Line(
          points={{3.9,-7},{13.33,-26.67}},
          color={11,43,197},
          smooth=Smooth.None));

      // Transport
      connect(O2.faces[facesCart[Axis.x], Side.n], xNegative.O2) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.faces[facesCart[Axis.x], Side.p], xPositive.O2) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(O2.faces[facesCart[Axis.y], Side.n], yNegative.O2) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{5.82867e-16,6.10623e-16},{
              5.82867e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.faces[facesCart[Axis.y], Side.p], yPositive.O2) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.faces[facesCart[Axis.z], Side.n], zNegative.O2) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(O2.faces[facesCart[Axis.z], Side.p], zPositive.O2) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      annotation (Documentation(info="<html>
<p>Please see the documentation of the <a href=\"modelica://FCSys.Subregions.Phases.BaseClasses.EmptyPhase\">EmptyPhase</a> model.</p></html>"));
    end Gas;

    model Graphite "Graphite phase"
      import assert = FCSys.BaseClasses.Utilities.assertEval;
      import FCSys.BaseClasses.Utilities.countTrue;
      extends FCSys.Subregions.Phases.BaseClasses.EmptyPhase(
        final n_spec=countTrue({'inclC+','incle-'}),
        initVelX=not 'inclC+',
        initVelY=not 'inclC+',
        initVelZ=not 'inclC+');

      // Conditionally include species.
      parameter Boolean 'inclC+'='incle-'
        "<html>Carbon plus (C<sup>+</sup>)</html>" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));

      replaceable Species.'C+'.Graphite.Fixed 'C+'(final n_faces) if 'inclC+'
        constrainedby Species.Species(
        n_faces=n_faces,
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

      replaceable Species.'e-'.Graphite.Fixed 'e-'(final n_faces) if 'incle-'
        constrainedby Species.Species(
        n_faces=n_faces,
        initMaterial=if 'inclC+' then InitScalar.None else InitScalar.Volume,
        initTransX=if (initVelX or 'inclC+') and reduceVel then
            InitTranslational.None else InitTranslational.Velocity,
        initTransY=if (initVelY or 'inclC+') and reduceVel then
            InitTranslational.None else InitTranslational.Velocity,
        initTransZ=if (initVelZ or 'inclC+') and reduceVel then
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

      ChemicalExchange HOR(n=-2, final m='e-'.Data.m) if inclHOR
        "Interface with the hydrogen oxidation reaction"
        annotation (Placement(transformation(extent={{-10,40},{-30,60}})));
      ChemicalExchange ORR(n=-4, final m='e-'.Data.m) if inclORR
        "Interface with the oxygen reduction reaction"
        annotation (Placement(transformation(extent={{-10,20},{-30,40}})));
      Connectors.ChemicalReaction chemical if inclHOR or inclORR
        "Connector for chemical reactions" annotation (Placement(transformation(
              extent={{-50,40},{-30,60}}), iconTransformation(extent={{-48,66},
                {-28,86}})));
      Depletion depletionHOR(
        A=A[Axis.x],
        transSubstrate=true,
        thermalSubstrate=true,
        side=FCSys.BaseClasses.Side.p,
        redeclare package Data = FCSys.Characteristics.'e-'.Gas) if inclHOR
        "Depletion region for the HOR"
        annotation (Placement(transformation(extent={{38,60},{18,80}})));

      Depletion depletionORR(
        A=A[Axis.x],
        transSubstrate=true,
        thermalSubstrate=true,
        side=FCSys.BaseClasses.Side.n,
        redeclare package Data = FCSys.Characteristics.'e-'.Gas) if inclORR
        "Depletion region for the ORR"
        annotation (Placement(transformation(extent={{-24,60},{-4,80}})));

      parameter Boolean inclHOR=false "Hydrogen oxidation" annotation (
        HideResult=true,
        Dialog(group="Included reactions", compact=true),
        choices(__Dymola_checkBox=true));
      parameter Boolean inclORR=false "Oxygen reduction" annotation (
        HideResult=true,
        Dialog(group="Included reactions", compact=true),
        choices(__Dymola_checkBox=true));
      // These can't be outer parameters in Dymola 7.4.

    initial equation
      assert(not (inclHOR and inclORR),
        "The HOR and ORR cannot be included simultaneously.");

    equation
      // Reactions
      // ---------
      // HOR
      connect('e-'.chemical, HOR.species) annotation (Line(
          points={{-3.9,7},{-18,50}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(HOR.reaction, chemical) annotation (Line(
          points={{-22,50},{-40,50}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(depletionHOR.face, 'e-'.faces[1, Side.p]) annotation (Line(
          points={{38,70},{50,70},{50,5.55112e-16},{5.55112e-16,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(depletionHOR.inert, 'C+'.inert) annotation (Line(
          points={{18,70},{6.9,70},{6.9,-3.8}},
          color={47,107,251},
          smooth=Smooth.None));
      connect(HOR.reaction, chemical) annotation (Line(
          points={{-22,50},{-40,50}},
          color={255,195,38},
          smooth=Smooth.None));
      // ORR
      connect('e-'.chemical, ORR.species) annotation (Line(
          points={{-3.9,7},{-18,30}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(ORR.reaction, chemical) annotation (Line(
          points={{-22,30},{-30,30},{-30,50},{-40,50}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(depletionORR.face, 'e-'.faces[1, Side.n]) annotation (Line(
          points={{-24,70},{-50,70},{-50,5.55112e-16},{5.55112e-16,5.55112e-16}},

          color={127,127,127},
          smooth=Smooth.None));

      connect(depletionORR.inert, 'C+'.inert) annotation (Line(
          points={{-4,70},{6.9,70},{6.9,-3.8}},
          color={47,107,251},
          smooth=Smooth.None));
      connect(ORR.reaction, chemical) annotation (Line(
          points={{-22,30},{-30,30},{-30,50},{-40,50}},
          color={255,195,38},
          smooth=Smooth.None));

      // C+
      // --
      // Exchange
      connect('C+'.inert.translational, common.translational) annotation (Line(
          points={{6.9,-3.8},{26.67,-13.33}},
          color={11,43,197},
          smooth=Smooth.None));
      connect('C+'.inert.thermal, common.thermal) annotation (Line(
          points={{6.9,-3.8},{26.67,-13.33}},
          color={11,43,197},
          smooth=Smooth.None));
      connect('C+'.inertDalton, inertDalton) annotation (Line(
          points={{3.9,-7},{13.33,-26.67}},
          color={11,43,197},
          smooth=Smooth.None));
      // Transport
      connect('C+'.faces[facesCart[Axis.x], Side.n], xNegative.'C+')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('C+'.faces[facesCart[Axis.x], Side.p], xPositive.'C+')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},
              {40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('C+'.faces[facesCart[Axis.y], Side.n], yNegative.'C+')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{6.10623e-16,-40},{5.55112e-16,-40}},

          color={127,127,127},
          smooth=Smooth.None));

      connect('C+'.faces[facesCart[Axis.y], Side.p], yPositive.'C+')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('C+'.faces[facesCart[Axis.z], Side.n], zNegative.'C+')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('C+'.faces[facesCart[Axis.z], Side.p], zPositive.'C+')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // e-
      // --
      // Exchange
      connect('e-'.inert.translational, common.translational) annotation (Line(
          points={{6.9,-3.8},{26.67,-13.33}},
          color={11,43,197},
          smooth=Smooth.None));
      connect('e-'.inert.thermal, common.thermal) annotation (Line(
          points={{6.9,-3.8},{26.67,-13.33}},
          color={47,107,251},
          smooth=Smooth.None));
      connect('e-'.inertDalton, inertDalton) annotation (Line(
          points={{3.9,-7},{13.33,-26.67}},
          color={11,43,197},
          smooth=Smooth.None));
      // Transport
      connect('e-'.faces[facesCart[Axis.x], Side.n], xNegative.'e-')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.faces[facesCart[Axis.x], Side.p], xPositive.'e-')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},
              {40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.faces[facesCart[Axis.y], Side.n], yNegative.'e-')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{6.10623e-16,-40},{5.55112e-16,-40}},

          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.faces[facesCart[Axis.y], Side.p], yPositive.'e-')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('e-'.faces[facesCart[Axis.z], Side.n], zNegative.'e-')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('e-'.faces[facesCart[Axis.z], Side.p], zPositive.'e-')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(depletionORR.chemical, chemical) annotation (Line(
          points={{-14,70},{-14,60},{-30,60},{-30,50},{-40,50}},
          color={255,195,38},
          smooth=Smooth.None));

      connect(depletionHOR.chemical, chemical) annotation (Line(
          points={{28,70},{28,60},{-30,60},{-30,50},{-40,50}},
          color={255,195,38},
          smooth=Smooth.None));
      annotation (Documentation(info="<html>
    <p>See <a href=\"modelica://FCSys.Subregions.Species.'e-'.Graphite.Fixed\">Species.'e-'.Graphite.Fixed</a> for assumptions.
    For more information, see the
 <a href=\"modelica://FCSys.Subregions.Phases.BaseClasses.EmptyPhase\">EmptyPhase</a> model.</p></html>"),
          Diagram(graphics));
    end Graphite;

    model Ionomer "Ionomer phase"
      import assert = FCSys.BaseClasses.Utilities.assertEval;
      import FCSys.BaseClasses.Utilities.countTrue;
      extends FCSys.Subregions.Phases.BaseClasses.EmptyPhase(
        final n_spec=countTrue({'inclC19HF37O5S-','inclH+',inclH2O}),
        initVelX=not 'inclC19HF37O5S-',
        initVelY=not 'inclC19HF37O5S-',
        initVelZ=not 'inclC19HF37O5S-');

      // Conditionally include species.
      parameter Boolean 'inclC19HF37O5S-'='inclH+'
        "<html>Nafion sulfonate minus (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>)</html>"
        annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));

      replaceable Species.'C19HF37O5S-'.Ionomer.Fixed 'C19HF37O5S-'(final
          n_faces) if 'inclC19HF37O5S-' constrainedby Species.Species(
        n_faces=n_faces,
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

      replaceable Species.'H+'.Ionomer.Fixed 'H+'(final n_faces) if 'inclH+'
        constrainedby Species.Species(
        n_faces=n_faces,
        initMaterial=if 'inclC19HF37O5S-' then InitScalar.None else InitScalar.Volume,

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

      replaceable Species.H2O.Ionomer.Fixed H2O(final n_faces) if inclH2O
        constrainedby Species.Species(
        n_faces=n_faces,
        initTransX=if (initVelX or 'inclC19HF37O5S-') and reduceVel then
            InitTranslational.None else InitTranslational.Velocity,
        initTransY=if (initVelY or 'inclC19HF37O5S-') and reduceVel then
            InitTranslational.None else InitTranslational.Velocity,
        initTransZ=if (initVelZ or 'inclC19HF37O5S-') and reduceVel then
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

      Connectors.PhysicalBus physical if inclH2O "Connector for phase change"
        annotation (Placement(transformation(extent={{-37,3},{-17,23}}),
            iconTransformation(extent={{-68,30},{-48,50}})));

      ChemicalExchange HOR(n=-2, final m='H+'.Data.m) if inclHOR
        "Interface with the hydrogen oxidation reaction"
        annotation (Placement(transformation(extent={{-10,40},{-30,60}})));
      ChemicalExchange ORR(n=-4, final m='H+'.Data.m) if inclORR
        "Interface with the oxygen reduction reaction"
        annotation (Placement(transformation(extent={{-10,20},{-30,40}})));
      Connectors.ChemicalReaction chemical if inclHOR or inclORR
        "Connector for chemical reactions" annotation (Placement(transformation(
              extent={{-50,40},{-30,60}}), iconTransformation(extent={{-48,66},
                {-28,86}})));
      Depletion depletionHOR(
        A=A[Axis.x],
        transSubstrate=false,
        thermalSubstrate=false,
        side=FCSys.BaseClasses.Side.n,
        redeclare package Data = FCSys.Characteristics.'H+'.Gas) if inclHOR
        "Depletion region for the HOR"
        annotation (Placement(transformation(extent={{-24,60},{-4,80}})));

      Depletion depletionORR(
        A=A[Axis.x],
        transSubstrate=false,
        thermalSubstrate=false,
        side=FCSys.BaseClasses.Side.p,
        redeclare package Data = FCSys.Characteristics.'H+'.Gas) if inclORR
        "Depletion region for the ORR"
        annotation (Placement(transformation(extent={{38,60},{18,80}})));

      parameter Boolean inclHOR=false "Hydrogen oxidation" annotation (
        HideResult=true,
        Dialog(group="Included reactions", compact=true),
        choices(__Dymola_checkBox=true));
      parameter Boolean inclORR=false "Oxygen reduction" annotation (
        HideResult=true,
        Dialog(group="Included reactions", compact=true),
        choices(__Dymola_checkBox=true));
      // These can't be outer parameters in Dymola 7.4.

    initial equation
      assert(not (inclHOR and inclORR),
        "The HOR and ORR cannot be included simultaneously.");

    equation
      // Phase change
      connect(physical.H2O, H2O.physical) annotation (Line(
          points={{-27,13},{-6.9,3.8}},
          color={38,196,52},
          smooth=Smooth.None));

      // Reactions
      // ---------
      // HOR
      connect('H+'.chemical, HOR.species) annotation (Line(
          points={{-3.9,7},{-18,50}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(HOR.reaction, chemical) annotation (Line(
          points={{-22,50},{-40,50}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(depletionHOR.face, 'H+'.faces[1, Side.n]) annotation (Line(
          points={{-24,70},{-50,70},{-50,0},{5.55112e-16,0},{5.55112e-16,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(depletionHOR.inert, 'C19HF37O5S-'.inert) annotation (Line(
          points={{-4,70},{6.9,70},{6.9,-3.8}},
          color={47,107,251},
          smooth=Smooth.None));
      connect(HOR.reaction, chemical) annotation (Line(
          points={{-22,50},{-40,50}},
          color={255,195,38},
          smooth=Smooth.None));
      // ORR
      connect('H+'.chemical, ORR.species) annotation (Line(
          points={{-3.9,7},{-18,30}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(ORR.reaction, chemical) annotation (Line(
          points={{-22,30},{-30,30},{-30,50},{-40,50}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(depletionORR.face, 'H+'.faces[1, Side.p]) annotation (Line(
          points={{38,70},{50,70},{50,5.55112e-16},{5.55112e-16,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(depletionORR.inert, 'C19HF37O5S-'.inert) annotation (Line(
          points={{18,70},{6.9,70},{6.9,-3.8}},
          color={47,107,251},
          smooth=Smooth.None));
      connect(ORR.reaction, chemical) annotation (Line(
          points={{-22,30},{-30,30},{-30,50},{-40,50}},
          color={255,195,38},
          smooth=Smooth.None));

      // C19HF37O5S-
      // -----------
      // Exchange
      connect('C19HF37O5S-'.inert.translational, common.translational)
        annotation (Line(
          points={{6.9,-3.8},{26.67,-13.33}},
          color={11,43,197},
          smooth=Smooth.None));
      connect('C19HF37O5S-'.inert.thermal, common.thermal) annotation (Line(
          points={{6.9,-3.8},{26.67,-13.33}},
          color={11,43,197},
          smooth=Smooth.None));
      connect('C19HF37O5S-'.inertDalton, inertDalton) annotation (Line(
          points={{3.9,-7},{4,-8},{13.33,-26.67}},
          color={11,43,197},
          smooth=Smooth.None));
      // Transport
      connect('C19HF37O5S-'.faces[facesCart[Axis.x], Side.n], xNegative.
        'C19HF37O5S-') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},
              {-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('C19HF37O5S-'.faces[facesCart[Axis.x], Side.p], xPositive.
        'C19HF37O5S-') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},
              {40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('C19HF37O5S-'.faces[facesCart[Axis.y], Side.n], yNegative.
        'C19HF37O5S-') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{6.10623e-16,-40},{5.55112e-16,-40}},

          color={127,127,127},
          smooth=Smooth.None));

      connect('C19HF37O5S-'.faces[facesCart[Axis.y], Side.p], yPositive.
        'C19HF37O5S-') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('C19HF37O5S-'.faces[facesCart[Axis.z], Side.n], zNegative.
        'C19HF37O5S-') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('C19HF37O5S-'.faces[facesCart[Axis.z], Side.p], zPositive.
        'C19HF37O5S-') annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // 'H+'
      // ----
      // Exchange
      connect('H+'.inert.translational, common.translational) annotation (Line(
          points={{6.9,-3.8},{26.67,-13.33}},
          color={11,43,197},
          smooth=Smooth.None));
      connect('H+'.inert.thermal, common.thermal) annotation (Line(
          points={{6.9,-3.8},{26.67,-13.33}},
          color={11,43,197},
          smooth=Smooth.None));
      connect('H+'.inertDalton, inertDalton) annotation (Line(
          points={{3.9,-7},{4,-8},{13.33,-26.67}},
          color={11,43,197},
          smooth=Smooth.None));
      // Transport
      connect('H+'.faces[facesCart[Axis.x], Side.n], xNegative.'H+')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},
              {-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.faces[facesCart[Axis.x], Side.p], xPositive.'H+')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},
              {40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.faces[facesCart[Axis.y], Side.n], yNegative.'H+')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{6.10623e-16,-40},{5.55112e-16,-40}},

          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.faces[facesCart[Axis.y], Side.p], yPositive.'H+')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect('H+'.faces[facesCart[Axis.z], Side.n], zNegative.'H+')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect('H+'.faces[facesCart[Axis.z], Side.p], zPositive.'H+')
        annotation (Line(
          points={{6.10623e-16,6.10623e-16},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));

      // H2O
      // ---
      // Exchange
      connect(H2O.inert.translational, common.translational) annotation (Line(
          points={{6.9,-3.8},{26.67,-13.33}},
          color={11,43,197},
          smooth=Smooth.None));
      connect(H2O.inert.thermal, common.thermal) annotation (Line(
          points={{6.9,-3.8},{26.67,-13.33}},
          color={47,107,251},
          smooth=Smooth.None));
      connect(H2O.inertDalton, inertDalton) annotation (Line(
          points={{3.9,-7},{13.33,-26.67}},
          color={11,43,197},
          smooth=Smooth.None));
      // Transport
      connect(H2O.faces[facesCart[Axis.x], Side.n], xNegative.H2O) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{-20,-4.87687e-22},{-20,5.55112e-16},
              {-40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[facesCart[Axis.x], Side.p], xPositive.H2O) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{20,-4.87687e-22},{20,5.55112e-16},
              {40,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[facesCart[Axis.y], Side.n], yNegative.H2O) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{6.10623e-16,-40},{5.55112e-16,-40}},

          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[facesCart[Axis.y], Side.p], yPositive.H2O) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[facesCart[Axis.z], Side.n], zNegative.H2O) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.faces[facesCart[Axis.z], Side.p], zPositive.H2O) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(depletionORR.chemical, chemical) annotation (Line(
          points={{28,70},{28,60},{-30,60},{-30,50},{-40,50}},
          color={255,195,38},
          smooth=Smooth.None));

      connect(depletionHOR.chemical, chemical) annotation (Line(
          points={{-14,70},{-14,60},{-30,60},{-30,50},{-40,50}},
          color={255,195,38},
          smooth=Smooth.None));
      annotation (Documentation(info="<html>
    <p>See <a href=\"modelica://FCSys.Subregions.Species.'H+'.Ionomer.Fixed\">Species.'H+'.Ionomer.Fixed</a> for assumptions.
    For more information, see the
 <a href=\"modelica://FCSys.Subregions.Phases.BaseClasses.EmptyPhase\">EmptyPhase</a> model.</p></html>"));
    end Ionomer;

    model Liquid "Liquid phase"
      extends FCSys.Subregions.Phases.BaseClasses.EmptyPhase(final n_spec=if
            inclH2O then 1 else 0);

      // Conditionally include species.
      parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
        annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Species.H2O.Liquid.Fixed H2O(final n_faces) if inclH2O
        constrainedby Species.Species(
        n_faces=n_faces,
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

      Connectors.PhysicalBus physical if inclH2O "Connector for phase change"
        annotation (Placement(transformation(extent={{-37,3},{-17,23}}),
            iconTransformation(extent={{-68,30},{-48,50}})));

    equation
      // Phase change
      connect(physical.H2O, H2O.physical) annotation (Line(
          points={{-27,13},{-6.9,3.8}},
          color={38,196,52},
          smooth=Smooth.None));

      // H2O
      // ---
      // Exchange
      connect(H2O.inert.translational, common.translational) annotation (Line(
          points={{6.9,-3.8},{26.67,-13.33}},
          color={11,43,197},
          smooth=Smooth.None));
      connect(H2O.inert.thermal, common.thermal) annotation (Line(
          points={{6.9,-3.8},{26.67,-13.33}},
          color={47,107,251},
          smooth=Smooth.None));
      connect(H2O.inertDalton, inertDalton) annotation (Line(
          points={{3.9,-7},{13.33,-26.67}},
          color={11,43,197},
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

      connect(H2O.faces[facesCart[Axis.y], Side.n], yNegative.H2O) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{6.10623e-16,-40},{5.55112e-16,-40}},

          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[facesCart[Axis.y], Side.p], yPositive.H2O) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{-4.87687e-22,20},{5.55112e-16,20},
              {5.55112e-16,40}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(H2O.faces[facesCart[Axis.z], Side.n], zNegative.H2O) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{20,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(H2O.faces[facesCart[Axis.z], Side.p], zPositive.H2O) annotation (
          Line(
          points={{6.10623e-16,6.10623e-16},{-20,-20}},
          color={127,127,127},
          smooth=Smooth.None));
      annotation (Documentation(info="<html><p>Please see the documentation of the
 <a href=\"modelica://FCSys.Subregions.Phases.BaseClasses.EmptyPhase\">EmptyPhase</a> model.</p></html>"));
    end Liquid;

    package BaseClasses "Base classes (generally not for direct use)"
      extends Modelica.Icons.BasesPackage;
      model EmptyPhase "Model for a phase with no species or reactions"
        import FCSys.BaseClasses.Utilities.index;
        // extends FCSys.BaseClasses.Icons.Names.Middle;

        parameter Integer n_spec(start=0) "Number of species"
          annotation (HideResult=true);

        // Geometry
        parameter Integer n_faces(min=1, max=3)
          "<html>Number of pairs of faces (<i>n</i><sub>faces</sub>)</html>"
          annotation (Dialog(group="Geometry"),HideResult=true);
        // This can't be an outer parameter in Dymola 7.4.
        parameter Q.NumberAbsolute k[Axis](
          each min=Modelica.Constants.small,
          each final nominal=1) = {1,1,1} if n_spec > 0
          "<html>Area fill factor (<b><i>k</i></b>)</html>"
          annotation (Dialog(group="Geometry"));

        // Assumptions
        parameter Boolean reduceVel=false "Same velocity for all species"
          annotation (Dialog(tab="Assumptions", enable=n_spec > 1), choices(
              __Dymola_checkBox=true));
        parameter Boolean reduceTemp=false "Same temperature for all species"
          annotation (Dialog(tab="Assumptions", enable=n_spec > 1), choices(
              __Dymola_checkBox=true));

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

        Connectors.InertDalton inertDalton(final n_trans=n_trans) annotation (
            Placement(transformation(extent={{3.33,-36.67},{23.33,-16.67}}),
              iconTransformation(extent={{70,-90},{90,-70}})));
        Connectors.FaceBus xPositive if inclFaces[Axis.x] and n_spec > 0
          "Positive face along the x axis" annotation (Placement(transformation(
                extent={{30,-10},{50,10}}), iconTransformation(extent={{70,-10},
                  {90,10}})));
        Connectors.FaceBus xNegative if inclFaces[Axis.x] and n_spec > 0
          "Negative face along the x axis" annotation (Placement(transformation(
                extent={{-50,-10},{-30,10}}), iconTransformation(extent={{-90,-10},
                  {-70,10}})));
        Connectors.FaceBus yPositive if inclFaces[Axis.y] and n_spec > 0
          "Positive face along the y axis" annotation (Placement(transformation(
                extent={{-10,30},{10,50}}), iconTransformation(extent={{-10,90},
                  {10,110}})));
        Connectors.FaceBus yNegative if inclFaces[Axis.y] and n_spec > 0
          "Negative face along the y axis" annotation (Placement(transformation(
                extent={{-10,-50},{10,-30}}), iconTransformation(extent={{-10,-94},
                  {10,-74}})));
        Connectors.FaceBus zPositive if inclFaces[Axis.z] and n_spec > 0
          "Positive face along the z axis" annotation (Placement(transformation(
                extent={{-30,-30},{-10,-10}}), iconTransformation(extent={{-90,
                  -90},{-70,-70}})));
        Connectors.FaceBus zNegative if inclFaces[Axis.z] and n_spec > 0
          "Negative face along the z axis" annotation (Placement(transformation(
                extent={{10,10},{30,30}}), iconTransformation(extent={{40,40},{
                  60,60}})));

        outer parameter Q.Length L[Axis] if n_spec > 0 "Length" annotation (
            HideResult=true,missingInnerMessage="This model should be used within a subregion model.
");
        outer parameter Q.Area A[Axis] if n_spec > 0 "Cross-sectional area"
          annotation (HideResult=true,missingInnerMessage="This model should be used within a subregion model.
");
        // Note:  These must be public in Dymola 7.4, so HideResult is used.

      protected
        final inner parameter Q.Length Lprime[Axis]=k .* A ./ L if n_spec > 0
          "Effective cross-sectional area per length";
        outer parameter Integer n_trans
          "Number of components of translational momentum" annotation (
            missingInnerMessage="This model should be used within a subregion model.
");
        outer parameter Integer cartTrans[:]
          "Cartesian-axis indices of the components of translational momentum"
          annotation (missingInnerMessage="This model should be used within a subregion model.
");
        outer parameter Integer facesCart[:]
          "Face-pair indices of the Cartesian axes" annotation (
            missingInnerMessage="This model should be used within a subregion model.
");
        outer parameter Boolean inclTrans[Axis]
          "true, if each component of translational momentum is included"
          annotation (missingInnerMessage="This model should be used within a subregion model.
");
        outer parameter Boolean inclFaces[Axis]
          "true, if each pairs of faces is included" annotation (
            missingInnerMessage="This model should be used within a subregion model.
");

        outer Conditions.Environment environment "Environmental conditions";
        // This component is conditional to prevent a mathematical singularity
        // when two or more empty phases (without any species included) are
        // connected.

        Connectors.InertInternal common(
          final n_trans=n_trans,
          final inclTranslational=reduceVel,
          final inclThermal=reduceTemp,
          translational(phi(
              each stateSelect=StateSelect.prefer,
              final start=phi_IC[cartTrans],
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
        // Note:  It would be simpler to use {initVelX, initVelY, initVelZ}[cartTrans]
        // for the fixed attribute of inert.translational, but Dymola 7.4 refuses to
        // accept it.

      equation
        if n_spec == 0 then
          inertDalton.V = 0;
          inertDalton.mPhidot = zeros(n_trans);
          inertDalton.Qdot = 0;
        end if;
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="phase",
          Documentation(info="<html><p>If one of the species has <code>consEnergy = Conservation.IC</code> and
    <code>reduceTemp</code> is <code>true</code>, then
    <code>initTemp</code> should be set to <code>false</code>.
    Likewise, if one of the species has <code>consTransX = Conservation.IC</code>,
    <code>consTransY = Conservation.IC</code>, or <code>consTransZ = Conservation.IC</code>, then
    <code>initVelX</code> and <code>reduceVel</code> is <code>true</code>, <code>initVelY</code>, or <code>initVelZ</code> should
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
                lineColor={0,0,0})}));
      end EmptyPhase;

    end BaseClasses;
    annotation (Documentation(info="<html><p>The graphite, ionomer, and
liquid phases can only be used with a compressible phase (gas).</p></html>"));

  end Phases;

  package Species
    "Models for single-species storage, transport, and exchange of material, translational momentum, and energy"
    extends Modelica.Icons.Package;
    package 'C+' "C"
      extends Modelica.Icons.Package;
      package Graphite "<html>C<sup>+</sup> graphite</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Correlations with adjustment factors"
          extends SpeciesSolid(
            redeclare replaceable package Data = Characteristics.'C+'.Graphite,

            nu=k_nu*Data.nu(T, v),
            theta=k_theta*Data.theta(T, v));

          // Note:  In Dymola 7.4,
          // "redeclare replaceable package Data = FCSys.Characteristics.C.Graphite"
          // must be used instead of
          // "redeclare replaceable FCSys.Characteristics.C.Graphite Data" in
          // order for this model to pass its check.  This applies to the other
          // species models too.

          parameter Q.NumberAbsolute k_nu(final nominal=1) = 1
            "<html>Adjustment factor for thermal independity (<i>k</i><sub>&nu;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'C+'",
            Documentation(info="<html><p>Please see the documentation of the
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
            defaultComponentName="'C+'",
            Documentation(info=
                  "<html><p>Please see the documentation of the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"),

            Icon(graphics));

        end Correlated;

        model Fixed "Fixed properties"

          extends SpeciesSolid(
            redeclare replaceable package Data = Characteristics.'C+'.Graphite
                (
                n_c=0,
                T_lim_c={0,Modelica.Constants.inf},
                b_c=[935*U.J*Data.m/(U.kg*U.K)],
                B_c=[Data.Deltah0_f - (935*U.J*Data.m/U.kg)*298.15, 154.663*U.J
                    /(U.mol*U.K) - (935*U.J*Data.m/(U.kg*U.K))*ln(298.15*U.K)]),

            redeclare final parameter Q.Mobility mu=0,
            redeclare parameter Q.TimeAbsolute nu=Data.nu(),
            redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(11.1*U.W));

          // In Dymola 7.4, the specific heat capacity must be entered explicitly in B_c
          // (i.e., 935*U.J*Data.m/(U.kg*U.K) instead of Data.b_c[1, 1]).

          // Note:  Parameter expressions (e.g., nu=Data.nu(environment.T)) are not
          // used here since they would render the parameters unadjustable in Dymola
          // 7.4.  This also applies to the other species.

          // See the documentation layer for a table of values for the specific heat
          // capacity and thermal resistivity.
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'C+'",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>The specific heat capacity is fixed (independent of temperature).</li>
    <li>The thermal independity and thermal resistivity are fixed (e.g., independent of temperature).</li>
    <li>Mobility is zero (by default).</li>
    </ol></p>

   <p>The default isobaric specific heat capacity (<code>b_c = [935*U.J*Data.m/(U.kg*U.K)]</code>)
   and thermal
   resistivity (<code>theta = U.m*U.K/(11.1*U.W)</code>) are for graphite fiber epoxy (25% vol)
   composite (with heat flow parallel to the fibers) at 300 K
   [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 909].
   The integration offset for specific entropy is set such that
   the specific entropy is 154.663 J/(mol&middot;K) at 25 &deg;C and <i>p</i><sup>o</sup> (1 atm).
   This is the value from Table B in [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>].
   Additional thermal data is listed in <a href=\"#Tab1\">Table 1</a>.</p>

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
          extends SpeciesSolid(
            redeclare replaceable package Data =
                Characteristics.'C19HF37O5S-'.Ionomer,
            nu=k_nu*Data.nu(T, v),
            theta=k_theta*Data.theta(T, v));

          parameter Q.NumberAbsolute k_nu(final nominal=1) = 1
            "<html>Adjustment factor for thermal independity (<i>k</i><sub>&nu;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
            "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'C19HF37O5S-'",
            Documentation(info=
                  "<html><p>Please see the documentation of the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"),

            Diagram(graphics));

        end Calibrated;

        model Correlated "Correlated properties"
          extends SpeciesSolid(redeclare replaceable package Data =
                Characteristics.'C19HF37O5S-'.Ionomer);
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'C19HF37O5S-'",
            Documentation(info=
                  "<html><p>Please see the documentation of the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Correlated;

        model Fixed "Fixed properties"
          extends SpeciesSolid(
            redeclare replaceable package Data =
                Characteristics.'C19HF37O5S-'.Ionomer,
            redeclare parameter Q.TimeAbsolute nu=Data.nu(),
            redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(0.16*U.W));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'C19HF37O5S-'",
            Documentation(info="<html><p>Assumptions:
    <ol>
    <li>The thermal independity and thermal resistivity are fixed (e.g., independent of temperature)</li>
    </ol></p>

    <p>The default thermal resistivity (<code>theta=U.m*U.K/(0.16*U.W)</code>) is of dry
  Nafion 115 [<a href=\"modelica://FCSys.UsersGuide.References\">Kandlikar2009</a>, p. 1277].</p>

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
          extends SpeciesIsochoric(
            redeclare replaceable package Data = Characteristics.'e-'.Gas,
            final tauprime=0,
            mu=k_mu*Data.mu(T, v),
            nu=k_nu*Data.nu(T, v),
            beta=k_beta*Data.beta(T, v),
            zeta=k_zeta*Data.zeta(T, v),
            theta=k_theta*Data.theta(T, v),
            initMaterial=InitScalar.Density);

          parameter Q.NumberAbsolute k_mu(final nominal=1) = 1
            "<html>Adjustment factor for mobility (<i>k</i><sub>&mu;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_nu(final nominal=1) = 1
            "<html>Adjustment factor for thermal independity (<i>k</i><sub>&nu;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_eta(final nominal=1) = 1
            "<html>Adjustment factor for material resistivity (<i>k</i><sub>&eta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
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
            Documentation(info="<html><p>Assumptions:<ol>
          <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change would be
          governed by other configurations.</li>
    </ol></p>

    <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Calibrated;

        model Correlated "Correlated properties"
          extends SpeciesIsochoric(
            redeclare replaceable package Data = Characteristics.'e-'.Gas,
            final tauprime=0,
            initMaterial=InitScalar.Density);

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'e-'",
            Documentation(info="<html><p>Assumptions:<ol>
          <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change would be
          governed by other configurations.</li>
    </ol></p>

    <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare replaceable package Data = Characteristics.'e-'.Gas,
            final tauprime=0,
            redeclare parameter Q.Mobility mu=1e-18*Data.mu(),
            redeclare parameter Q.TimeAbsolute nu=Data.nu(),
            redeclare final parameter Q.Mobility eta=0,
            redeclare final parameter Q.Fluidity beta=0,
            redeclare parameter Q.Fluidity zeta=Data.zeta(),
            redeclare parameter Q.ResistivityThermal theta=Data.theta(),
            initMaterial=InitScalar.Density);
          //    redeclare parameter Q.Mobility eta=Data.eta(),
          //    redeclare parameter Q.Fluidity beta=1e-5*Data.beta(),
          // **set rho_IC, do the same for Correlated and Calibrated.
          // **temp factor on beta
          // add eta to other e- and H+ models.
          annotation (
            group="Material properties",
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'e-'",
            Documentation(info="<html><p>Assumptions:<ol>
          <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change would be
          governed by other configurations.</li>
    </ol></p>

    <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"),

            Diagram(graphics));

        end Fixed;

      end Graphite;

    end 'e-';

    package 'H+' "<html>H<sup>+</sup></html>"
      extends Modelica.Icons.Package;
      package Ionomer "<html>H<sup>+</sup> in ionomer</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Correlations with adjustment factors"
          extends SpeciesIsochoric(
            redeclare replaceable package Data = Characteristics.'H+'.Gas,
            initMaterial=InitScalar.Density,
            final tauprime=0,
            mu=k_mu*Data.mu(T, v),
            nu=k_nu*Data.nu(T, v),
            beta=k_beta*Data.beta(T, v),
            zeta=k_zeta*Data.zeta(T, v),
            theta=k_theta*Data.theta(T, v));

          parameter Q.NumberAbsolute k_mu(final nominal=1) = 1
            "<html>Adjustment factor for mobility (<i>k</i><sub>&mu;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_nu(final nominal=1) = 1
            "<html>Adjustment factor for thermal independity (<i>k</i><sub>&nu;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
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
  (not C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S) for charge neutrality.</li>
          <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change would be
          governed by other configurations.</li>
    </ol></p>

  <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Calibrated;

        model Correlated "Correlated properties"
          extends SpeciesIsochoric(
            redeclare replaceable package Data = Characteristics.'H+'.Gas,
            final tauprime=0,
            initMaterial=InitScalar.Density);
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="'H+'",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>The density of H<sup>+</sup> is equal to that of
  C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup> or approximately 1.912 M.  Note that
  this is greater than that measured by Spry and Fayer (see
  <a href=\"modelica://FCSys.Characteristics.'H+'.Ionomer\">Characteristics.'H+'.Ionomer</a>), but it
  simplifies the model by requiring only C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>
  (not C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S) for charge neutrality.</li>
          <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change would be
          governed by other configurations.</li>
    </ol></p>

    <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare replaceable package Data = Characteristics.'H+'.Gas,
            final tauprime=0,
            redeclare parameter Q.Mobility mu=Data.mu(),
            redeclare parameter Q.TimeAbsolute nu=Data.nu(),
            redeclare final parameter Q.Mobility eta=0,
            redeclare final parameter Q.Fluidity beta=0,
            redeclare parameter Q.Fluidity zeta=1/(5.3e-6*U.Pa*U.s),
            redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(0.1661*U.W),

            initMaterial=InitScalar.Density);
          // **set rho_IC, do the same for Correlated and Calibrated.

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
  (not C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S) for charge neutrality.</li>
          <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change would be
          governed by other configurations.</li>
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
                  b_v=[1], n_v={-1,0}),
            tauprime=k_tauprime*Data.tauprime(T, v),
            mu=k_mu*Data.mu(T, v),
            nu=k_nu*Data.nu(T, v),
            eta=k_eta*Data.eta(T, v),
            beta=k_beta*Data.beta(T, v),
            zeta=k_zeta*Data.zeta(T, v),
            theta=k_theta*Data.theta(T, v));
          parameter Q.NumberAbsolute k_tauprime(final nominal=1) = 1
            "<html>Adjustment factor for the phase change interval (<i>k</i><sub>&tau;&prime;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_mu(final nominal=1) = 1
            "<html>Adjustment factor for mobility (<i>k</i><sub>&mu;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_nu(final nominal=1) = 1
            "<html>Adjustment factor for thermal independity (<i>k</i><sub>&nu;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_eta(final nominal=1) = 1
            "<html>Adjustment factor for material resistivity (<i>k</i><sub>&eta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
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

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Calibrated;

        model Correlated "Correlated properties"
          extends Species(redeclare replaceable package Data =
                Characteristics.H2.Gas (b_v=[1], n_v={-1,0}));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.H2.Gas (
                  b_v=[1], n_v={-1,0}),
            redeclare parameter Q.TimeAbsolute tauprime=Data.tauprime(),
            redeclare parameter Q.Mobility mu=Data.mu(),
            redeclare parameter Q.TimeAbsolute nu=Data.nu(),
            redeclare parameter Q.ResistivityMaterial eta=Data.eta(),
            redeclare parameter Q.Fluidity beta=Data.beta(),
            redeclare parameter Q.Fluidity zeta=1/(8.96e-6*U.Pa*U.s),
            redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(0.183*U.W));
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
<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"),

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
            tauprime=k_tauprime*Data.tauprime(T, v),
            mu=k_mu*Data.mu(T, v),
            nu=k_nu*Data.nu(T, v),
            eta=k_eta*Data.eta(T, v),
            beta=k_beta*Data.beta(T, v),
            zeta=k_zeta*Data.zeta(T, v),
            theta=k_theta*Data.theta(T, v));
          parameter Q.NumberAbsolute k_tauprime(final nominal=1) = 1
            "<html>Adjustment factor for the phase change interval (<i>k</i><sub>&tau;&prime;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_mu(final nominal=1) = 1
            "<html>Adjustment factor for mobility (<i>k</i><sub>&mu;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_nu(final nominal=1) = 1
            "<html>Adjustment factor for thermal independity (<i>k</i><sub>&nu;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_eta(final nominal=1) = 1
            "<html>Adjustment factor for material resistivity (<i>k</i><sub>&eta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
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

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Calibrated;

        model Correlated "Correlated properties"
          extends Species(redeclare replaceable package Data =
                Characteristics.H2O.Gas (b_v=[1], n_v={-1,0}));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.H2O.Gas
                (b_v=[1], n_v={-1,0}),
            redeclare parameter Q.TimeAbsolute tauprime=Data.tauprime(),
            redeclare parameter Q.Mobility mu=Data.mu(),
            redeclare parameter Q.TimeAbsolute nu=Data.nu(),
            redeclare parameter Q.ResistivityMaterial eta=Data.eta(),
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
  </table></ul></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Fixed;

      end Gas;

      package Ionomer "<html>H<sub>2</sub>O in ionomer</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Correlations with adjustment factors"
          extends Species(
            redeclare replaceable package Data = Characteristics.H2O.Ionomer,
            final tauprime=0,
            mu=k_mu*Data.mu(T, v),
            nu=k_nu*Data.nu(T, v),
            eta=k_eta*Data.eta(T, v),
            beta=k_beta*Data.beta(T, v),
            zeta=k_zeta*Data.zeta(T, v),
            theta=k_theta*Data.theta(T, v));

          parameter Q.NumberAbsolute k_mu(final nominal=1) = 1
            "<html>Adjustment factor for mobility (<i>k</i><sub>&mu;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_nu(final nominal=1) = 1
            "<html>Adjustment factor for thermal independity (<i>k</i><sub>&nu;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_eta(final nominal=1) = 1
            "<html>Adjustment factor for material resistivity (<i>k</i><sub>&eta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
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
        <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change is governed by the
        other configurations (e.g., gas).</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Calibrated;

        model Correlated "Correlated properties"
          extends Species(redeclare replaceable package Data =
                Characteristics.H2O.Ionomer (b_v=[1], n_v={-1,0}),final
              tauprime=0);
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
    <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change is governed by the
   other configurations (e.g., gas).</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare replaceable package Data = Characteristics.H2O.Ionomer (
                  b_v=[1], n_v={-1,0}),
            final tauprime=0,
            redeclare parameter Q.Mobility mu=Data.mu(),
            redeclare parameter Q.TimeAbsolute nu=Data.nu(),
            redeclare parameter Q.ResistivityMaterial eta=Data.eta(),
            redeclare parameter Q.Fluidity beta=Data.beta(),
            redeclare parameter Q.Fluidity zeta=Data.zeta(),
            redeclare parameter Q.ResistivityThermal theta=Data.theta());

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
        <li>The generalized resistivities (&beta;, &zeta;, &theta;) are fixed (e.g., independent of temperature).</li>
        <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change is governed by the
        other configurations (e.g., gas).</li>
    </ol></p></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Fixed;

      end Ionomer;

      package Liquid "<html>H<sub>2</sub>O liquid</html>"
        extends Modelica.Icons.Package;
        model Calibrated "Correlations with adjustment factors"
          extends FCSys.Subregions.Species.SpeciesIsochoric(
            redeclare replaceable package Data = Characteristics.H2O.Liquid,
            final tauprime=0,
            mu=k_mu*Data.mu(T, v),
            nu=k_nu*Data.nu(T, v),
            beta=k_beta*Data.beta(T, v),
            zeta=k_zeta*Data.zeta(T, v),
            theta=k_theta*Data.theta(T, v));

          parameter Q.NumberAbsolute k_mu(final nominal=1) = 1
            "<html>Adjustment factor for mobility (<i>k</i><sub>&mu;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_nu(final nominal=1) = 1
            "<html>Adjustment factor for thermal independity (<i>k</i><sub>&nu;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
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
        <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change is governed by the
        other configurations (e.g., gas).</li>
    </ol></p>
         <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Calibrated;

        model Correlated "Correlated properties"
          extends FCSys.Subregions.Species.SpeciesIsochoric(redeclare
              replaceable package Data = Characteristics.H2O.Liquid, final
              tauprime=0);
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
        <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change is governed by the
        other configurations (e.g., gas).</li>
    </ol></p>
         <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Correlated;

        model Fixed "Fixed properties"
          extends FCSys.Subregions.Species.SpeciesIsochoric(
            redeclare replaceable package Data = Characteristics.H2O.Liquid,
            final tauprime=0,
            redeclare parameter Q.Mobility mu=Data.mu(),
            redeclare parameter Q.TimeAbsolute nu=Data.nu(),
            redeclare parameter Q.Fluidity beta=Data.beta(),
            redeclare parameter Q.Fluidity zeta=1/(855e-6*U.Pa*U.s),
            redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(0.613*U.W));

          // See the documentation for tables of values.
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="H2O",
            Documentation(info="<html><p>Assumptions:<ol>
        <li>The generalized resistivities (&beta;, &zeta;, &theta;) are fixed (e.g., independent of temperature).</li>
        <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change is governed by the
        other configurations (e.g., gas).</li>
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

  <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

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
                  b_v=[1], n_v={-1,0}),
            tauprime=k_tauprime*Data.tauprime(T, v),
            mu=k_mu*Data.mu(T, v),
            nu=k_nu*Data.nu(T, v),
            eta=k_eta*Data.eta(T, v),
            beta=k_beta*Data.beta(T, v),
            zeta=k_zeta*Data.zeta(T, v),
            theta=k_theta*Data.theta(T, v));

          parameter Q.NumberAbsolute k_tauprime(final nominal=1) = 1
            "<html>Adjustment factor for the phase change interval (<i>k</i><sub>&tau;&prime;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_mu(final nominal=1) = 1
            "<html>Adjustment factor for mobility (<i>k</i><sub>&mu;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_nu(final nominal=1) = 1
            "<html>Adjustment factor for thermal independity (<i>k</i><sub>&nu;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_eta(final nominal=1) = 1
            "<html>Adjustment factor for material resistivity (<i>k</i><sub>&eta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
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

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Calibrated;

        model Correlated "Correlated properties"
          extends Species(redeclare replaceable package Data =
                Characteristics.N2.Gas (b_v=[1], n_v={-1,0}));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="N2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Correlated;

        model Fixed "Fixed properties"
          import FCSys.BaseClasses.Utilities.Polynomial;

          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.N2.Gas (
                b_v=[1],
                n_v={-1,0},
                n_c=0,
                T_lim_c={0,Modelica.Constants.inf},
                b_c=[1041*U.J*Data.m/(U.kg*U.K)],
                B_c=[Data.Deltah0_f - (1041*U.J*Data.m/U.kg)*298.15, 191.610*U.J
                    /(U.mol*U.K) - (1041*U.J*Data.m/(U.kg*U.K))*ln(298.15*U.K)]),

            redeclare parameter Q.TimeAbsolute tauprime=Data.tauprime(),
            redeclare parameter Q.Mobility mu=Data.mu(),
            redeclare parameter Q.TimeAbsolute nu=Data.nu(),
            redeclare parameter Q.ResistivityMaterial eta=Data.eta(),
            redeclare parameter Q.Fluidity beta=Data.beta(),
            redeclare parameter Q.Fluidity zeta=1/(17.82e-6*U.Pa*U.s),
            redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(25.9e-3*U.W));

          // See the documentation for a table of values.
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="N2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
    <li>Fixed specific heat capacity (independent of temperature)</ol>
        <li>The generalized resistivities (&beta;, &zeta;, &theta;) are fixed (e.g., independent of temperature).</li>

    </ol></p>

<p>The default specific heat capacity (<code>b_c=[1.041e3*U.J*Data.m/(U.kg*U.K)]</code>) and resistivities
(<code>zeta=1/(17.82e-6*U.Pa*U.s)</code> and <code>theta=U.m*U.K/(25.9e-3*U.W))</code>) are based on data of gas at 1 atm and
  300 K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 920].
   The integration offset for specific entropy is set such that
   the specific entropy is 191.610 J/(mol&middot;K) at 25 &deg;C and <i>p</i><sup>o</sup> (1 bar).
   This is the value from Table B in [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>].
   Additional thermal data is listed in <a href=\"#Tab1\">Table 1</a>.  <a href=\"#Tab2\">Table 2</a> lists
  values of the material resistivity or self diffusion coefficient.</p>

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
<tr><td>800</td><td>1.220e3 [sic]</td><td>1/349.1e-7</td><td>1/54.8e-3</td></tr>
<tr><td>900</td><td>1.146e3</td><td>1/375.3e-7</td><td>1/59.7e-3</td></tr>
<tr><td>1000</td><td>1.167e3</td><td>1/399.9e-7</td><td>1/64.7e-3</td></tr>
<tr><td>1100</td><td>1.187e3</td><td>1/423.2e-7</td><td>1/70.0e-3</td></tr>
<tr><td>1200</td><td>1.204e3</td><td>1/445.3e-7</td><td>1/75.8e-3</td></tr>
<tr><td>1300</td><td>1.219e3</td><td>1/466.2e-7</td><td>1/81.0e-3</td></tr>
  </table>

<br>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"top\" id=\"Tab2\">Table 2: Material resistivity of N<sub>2</sub> gas at 1 atm
  [<a href=\"modelica://FCSys.UsersGuide.References\">Present1958</a>, p. 263]</caption>
  <tr>
      <th valign=\"middle\"><code>T<br>/U.K</code></th>
      <th width=1><code>eta*U.s<br>/U.cm^2</code></th>
    </tr>
<tr><td>77.7</td><td>0.0168</td></tr>
<tr><td>194.7</td><td>0.104</td></tr>
<tr><td>273.2</td><td>0.185</td></tr>
<tr><td>353.2</td><td>0.287</td></tr>
  </table></p>

  <p>The fluidity of air at 15.0 &deg;C and 1 atm is given by
       <code>zeta=1/(17.8e-6*U.Pa*U.s)</code>
   (<a href=\"http://en.wikipedia.org/wiki/Viscosity\">http://en.wikipedia.org/wiki/Viscosity</a>).</p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

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
                  b_v=[1], n_v={-1,0}),
            tauprime=k_tauprime*Data.tauprime(T, v),
            mu=k_mu*Data.mu(T, v),
            nu=k_nu*Data.nu(T, v),
            eta=k_eta*Data.eta(T, v),
            beta=k_beta*Data.beta(T, v),
            zeta=k_zeta*Data.zeta(T, v),
            theta=k_theta*Data.theta(T, v));

          parameter Q.NumberAbsolute k_tauprime(final nominal=1) = 1
            "<html>Adjustment factor for the phase change interval (<i>k</i><sub>&tau;&prime;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_mu(final nominal=1) = 1
            "<html>Adjustment factor for mobility (<i>k</i><sub>&mu;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_nu(final nominal=1) = 1
            "<html>Adjustment factor for thermal independity (<i>k</i><sub>&nu;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
          parameter Q.NumberAbsolute k_eta(final nominal=1) = 1
            "<html>Adjustment factor for material resistivity (<i>k</i><sub>&eta;</sub>)</html>"
            annotation (Dialog(group="Material properties"));
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

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Calibrated;

        model Correlated "Correlated properties"
          extends Species(redeclare replaceable package Data =
                Characteristics.O2.Gas (b_v=[1], n_v={-1,0}));
          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="O2",
            Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

        end Correlated;

        model Fixed "Fixed properties"
          extends Species(
            redeclare replaceable package Data = FCSys.Characteristics.O2.Gas (
                  b_v=[1], n_v={-1,0}),
            redeclare parameter Q.TimeAbsolute tauprime=Data.tauprime(),
            redeclare parameter Q.Mobility mu=Data.mu(),
            redeclare parameter Q.TimeAbsolute nu=Data.nu(),
            redeclare parameter Q.ResistivityMaterial eta=Data.eta(),
            redeclare parameter Q.Fluidity beta=Data.beta(),
            redeclare parameter Q.Fluidity zeta=1/(20.72e-6*U.Pa*U.s),
            redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(26.8e-3*U.W));

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
  300 K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 920&ndash;921].
  <a href=\"#Tab1\">Table 1</a> lists the properties at other temperatures. <a href=\"#Tab2\">Table 2</a> lists
  values of the material resistivity or self diffusion coefficient.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"top\" id=\"Tab1\">Table 1: Properties of O<sub>2</sub> gas at 1 atm
  [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 920&ndash;921]</caption>
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

<br>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"top\" id=\"Tab2\">Table 2: Material resistivity of O<sub>2</sub> gas at 1 atm
  [<a href=\"modelica://FCSys.UsersGuide.References\">Present1958</a>, p. 263]</caption>
  <tr>
      <th valign=\"middle\"><code>T<br>/U.K</code></th>
      <th width=1><code>eta*U.s<br>/U.cm^2</code></th>
    </tr>
<tr><td>77.7</td><td>0.0153</td></tr>
<tr><td>194.7</td><td>0.104</td></tr>
<tr><td>273.2</td><td>0.187</td></tr>
<tr><td>353.2</td><td>0.301</td></tr>
  </table></p>

<p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"),

            Icon(graphics));

        end Fixed;

      end Gas;

    end O2;

    model SpeciesSolid
      "<html><a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model for a solid (inert and zero velocity)</html>"
      extends SpeciesIsochoric(
        final upstreamX=false,
        final upstreamY=false,
        final upstreamZ=false,
        final phi_IC=zeros(3),
        final I_IC,
        final consTransX=Conservation.IC,
        final consTransY=Conservation.IC,
        final consTransZ=Conservation.IC,
        final tauprime=0,
        final beta=1,
        final zeta=1);
      // Note:  beta and zeta don't matter as long as they are nonzero.
      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="species",
        Documentation(info="<html><p>Assumptions:<ol>
  <li>Zero dynamic compressibility (&rArr; uniform velocity in the axial direction)</li>
  <li>Zero fluidity (&rArr; no shearing)</li></ol></p>

  <p>Usually, conditions should be applied to specify the velocity (typically <b>0</b>).  In a group of connected solid species
  of a single type (instances of a model derived from this one), there should be exactly one equation to specify the velocity
  along each Cartesian axis.
  For example, the x-axis velocity may be given by setting <code>consTransX</code> to <code>Conservation.IC</code> in one of the instances
  (x-axis velocity will be <code>phi_IC[Axis.x]</code> for all time).  Alternatively, the velocity at one of the faces on the
  outside of the
  group could be specified using a boundary condition from the <a href=\"modelica://FCSys.Conditions\">Conditions</a> package.</p>

  <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

    end SpeciesSolid;

    model SpeciesIsochoric
      "<html><a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model for an isochoric material</html>"
      extends Species(
        initMaterial=InitScalar.Volume,
        invertEOS=false,
        final eta=1);

      // Note:  Pressure, which is the default material IC for the base model,
      // can't be used to initialize an incompressible species.
      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="species",
        Documentation(info="<html>
  <p>Please see the documentation of the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));

    end SpeciesIsochoric;

    model Species
      "Model to exchange, transport, and store the material, momentum, and energy of one species"
      import FCSys.BaseClasses.Utilities.cartWrap;
      import FCSys.BaseClasses.Utilities.inSign;
      import FCSys.BaseClasses.Utilities.Delta;
      import FCSys.BaseClasses.Utilities.Sigma;
      import assert = FCSys.BaseClasses.Utilities.assertEval;
      extends FCSys.BaseClasses.Icons.Names.Top4;

      // Geometry
      parameter Integer n_faces(
        min=1,
        max=3) = 1
        "<html>Number of pairs of faces (<i>n</i><sub>faces</sub>)</html>"
        annotation (Dialog(group="Geometry"),HideResult=true);
      // Note:  This can't be an outer parameter in Dymola 7.4.

      // Material properties
      replaceable package Data = Characteristics.BaseClasses.Characteristic
        constrainedby Characteristics.BaseClasses.Characteristic
        "Characteristic data" annotation (
        Dialog(group="Material properties"),
        choicesAllMatching=true,
        __Dymola_choicesFromPackage=true,
        Placement(transformation(extent={{-60,40},{-40,60}}),
            iconTransformation(extent={{-10,90},{10,110}})));
      Q.TimeAbsolute tauprime(nominal=1e-6*U.s) = Data.tauprime(T, v)
        "<html>Phase change interval (&tau;&prime;)</html>"
        annotation (Dialog(group="Material properties"));
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

      // Assumptions
      // -----------
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
      //
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
      // TODO:  If the static option isn't useful, remove it and go back
      // to Boolean setMaterial, setTransX, etc.
      //
      // Flow conditions
      parameter Q.NumberAbsolute Nu_Phi[Axis]={4,4,4}
        "<html>Translational Nusselt numbers (<b><i>Nu</i><sub>&Phi;</sub></b>)</html>"
        annotation (Dialog(tab="Assumptions",group="Flow conditions"));
      parameter Q.NumberAbsolute Nu_Q=3.66
        "<html>Thermal Nusselt number (<i>Nu</i><sub><i>Q</i></sub>)</html>"
        annotation (Dialog(tab="Assumptions",group="Flow conditions"));

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
      parameter Q.Potential h_IC(start=Data.h(T_IC, p_IC),displayUnit="kJ/mol")
        "<html>Initial specific enthalpy (<i>h</i><sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization", group="Material and energy"));
      parameter Q.Potential g_IC(start=Data.g(T_IC, p_IC),displayUnit="kJ/mol")
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

      // Advanced parameters
      parameter Boolean invertEOS=true "Invert the equation of state"
        annotation (Dialog(
          Evaluate=true,
          tab="Advanced",
          compact=true), choices(__Dymola_checkBox=true));

      // Preferred states
      Q.Amount N(
        min=Modelica.Constants.small,
        nominal=4*U.C,
        final start=N_IC,
        final fixed=false,
        stateSelect=StateSelect.prefer) "Particle number";
      // Note:  The start value for this variable (and others below) isn't fixed
      // because the related initial condition is applied in the initial
      // equation section.
      Q.Velocity phi[n_trans](
        each nominal=10*U.cm/U.s,
        final start=phi_IC[cartTrans],
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
      // StateSelect.never avoids dynamic state selection of this variable and others
      // below in Dymola 7.4.
      Q.Mass M(
        nominal=1e-3*U.g,
        final start=Data.m*N_IC,
        stateSelect=StateSelect.never) "Mass";
      Q.Volume V(
        nominal=U.cc,
        final start=V_IC,
        final fixed=false,
        stateSelect=StateSelect.never) "Volume";
      Q.VolumeSpecific v(
        nominal=U.cc/(4*U.C),
        final start=1/rho_IC,
        final fixed=false,
        stateSelect=StateSelect.never) "Specific volume";
      Q.Potential h(
        nominal=U.V,
        final start=h_IC,
        final fixed=false,
        stateSelect=StateSelect.never) "Specific enthalpy";
      Q.NumberAbsolute s(
        nominal=10,
        final start=(h_IC - g_IC)/T_IC,
        stateSelect=StateSelect.never) "Specific entropy";
      // StateSelect.never avoids dynamic state selection Dymola 7.4.
      Q.Current I[n_trans](
        each nominal=U.A,
        final start=I_IC[cartTrans],
        each final fixed=false) "Current";
      Q.Current Ndot_faces[n_faces, Side](nominal=U.A, final start=outerProduct(
            I_IC[cartFaces], {1,-1}))
        "Total current into the faces (advection and diffusion)";
      Q.PressureAbsolute p_faces[n_faces, Side](each nominal=U.atm, each start=
            p_IC) "Thermodynamic pressures at the faces";

      // Auxiliary variables (for analysis)
      // ----------------------------------
      // Misc. properties and conditions
      output Q.Density rho(stateSelect=StateSelect.never) = 1/v if environment.analysis
        "Density";
      output Q.Potential g(stateSelect=StateSelect.never) = chemical.mu if
        environment.analysis "Electrochemical potential";
      output Q.Amount S(stateSelect=StateSelect.never) = N*s if environment.analysis
        "Entropy";
      output Q.PressureAbsolute q[n_trans](each stateSelect=StateSelect.never)
         = Data.m*phi .* I ./ (2*A[cartTrans]) if environment.analysis
        "Dynamic pressure";
      output Q.CapacityThermalSpecific c_p(stateSelect=StateSelect.never) =
        Data.c_p(T, p) if environment.analysis
        "Isobaric specific heat capacity";
      output Q.CapacityThermalSpecific c_v(stateSelect=StateSelect.never) =
        Data.c_v(T, p) if environment.analysis
        "Isochoric specific heat capacity";
      output Q.PressureReciprocal kappa(stateSelect=StateSelect.never) =
        Data.kappa(T, p) if environment.analysis "Isothermal compressibility";
      //
      // Time constants (only for the axes with translational momentum included; others
      // are infinite)
      output Q.TimeAbsolute tau_NE(
        stateSelect=StateSelect.never,
        start=U.s) = kappa*tauprime*exp((chemical.mu - g0)/T)*T/v if
        environment.analysis "Time constant for phase change";
      output Q.TimeAbsolute tau_PhiE(
        stateSelect=StateSelect.never,
        start=U.s) = Data.m*mu if environment.analysis
        "Time constant for translational exchange";
      output Q.TimeAbsolute tau_QE(
        stateSelect=StateSelect.never,
        start=U.s) = c_p*nu if environment.analysis
        "Time constant for thermal exchange";
      output Q.TimeAbsolute tau_NT[n_faces](
        each stateSelect=StateSelect.never,
        each start=U.s) = fill(V*eta/2, n_faces) ./ Lprime[cartFaces] if
        environment.analysis "Time constants for material transport";
      output Q.TimeAbsolute tau_PhiT_perp[n_faces](
        each stateSelect=StateSelect.never,
        each start=U.s) = fill(M*beta/2, n_faces) ./ Lprime[cartFaces] if
        environment.analysis
        "Time constants for normal translational transport";
      output Q.TimeAbsolute tau_PhiT_para[n_faces](
        each stateSelect=StateSelect.never,
        each start=U.s) = fill(M*zeta/2, n_faces) ./ Lprime[cartFaces] if
        environment.analysis
        "Time constants for transverse translational transport";
      output Q.TimeAbsolute tau_QT[n_faces](
        each stateSelect=StateSelect.never,
        each start=U.s) = fill(N*c_v*theta/2, n_faces) ./ Lprime[cartFaces] if
        environment.analysis "Time constants for thermal transport";
      //
      // Peclet numbers (only for the axes with translational momentum included; others
      // are zero)
      output Q.Number Pe_N[n_trans](each stateSelect=StateSelect.never) = eta*v
        *I ./ Lprime[cartTrans] if environment.analysis
        "Material Peclet numbers";
      output Q.Number Pe_Phi_perp[n_trans](each stateSelect=StateSelect.never)
         = beta*Data.m*I ./ Lprime[cartTrans] if environment.analysis
        "Normal translational Peclet numbers";
      output Q.Number Pe_Phi_para[n_trans](each stateSelect=StateSelect.never)
         = zeta*Data.m*I ./ Lprime[cartTrans] if environment.analysis
        "Transverse translational Peclet numbers";
      output Q.Number Pe_Q[n_trans](each stateSelect=StateSelect.never) = theta
        *Data.c_v(T, p)*I ./ Lprime[cartTrans] if environment.analysis
        "Thermal Peclet numbers";
      //
      // Bulk flow rates
      output Q.Force mphiI[n_trans, n_trans](each stateSelect=StateSelect.never)
         = outerProduct(I, Data.m*phi) if environment.analysis
        "Bulk rate of translational advection (1st index: transport axis, 2nd index: translational component)";
      output Q.Power hI[n_trans](each stateSelect=StateSelect.never) = h*I if
        environment.analysis "Bulk enthalpy flow rate";
      //
      // Translational momentum balance
      output Q.Force Ma[n_trans](each stateSelect=StateSelect.never) = M*(der(
        phi)/U.s + environment.a[cartTrans]) + N*Data.z*environment.E[cartTrans]
        if environment.analysis
        "Acceleration force (including acceleration due to body forces)";
      output Q.Force f_thermo[n_trans](each stateSelect=StateSelect.never) = {(
        if inclFaces[cartTrans[i]] then -Delta(p_faces[facesCart[cartTrans[i]],
        :])*A[cartTrans[i]] else 0) for i in 1:n_trans} "Thermodynamic force";
      output Q.Force f_AE[n_trans](each stateSelect=StateSelect.never) = Data.m
        *((actualStream(chemical.phi) - phi) .* chemical.Ndot + (actualStream(
        physical.phi) - phi) .* physical.Ndot) if environment.analysis
        "Acceleration force due to advective exchange";
      output Q.Force f_DE[n_trans](each stateSelect=StateSelect.never) = inert.translational.mPhidot
         + inertDalton.mPhidot if environment.analysis
        "Friction from other configurations (diffusive exchange)";
      output Q.Force f_AT[n_trans](each stateSelect=StateSelect.never) = {sum((
        faces[j, :].phi[cartWrap(cartTrans[i] - cartFaces[j] + 1)] - {phi[i],
        phi[i]})*Ndot_faces[j, :]*Data.m for j in 1:n_faces) for i in 1:n_trans}
        if environment.analysis "Acceleration force due to advective transport";
      output Q.Force f_DT[n_trans](each stateSelect=StateSelect.never) = {sum(
        Sigma(faces[j, :].mPhidot[cartWrap(cartTrans[i] - cartFaces[j] + 1)])
        for j in 1:n_faces) for i in 1:n_trans} if environment.analysis
        "Friction from other subregions (diffusive transport, including bulk viscosity)";
      //
      // Energy balance
      output Q.Power Ndere(stateSelect=StateSelect.never) = (N*T*der(s) + M*phi
        *der(phi))/U.s if environment.analysis
        "Rate of energy storage (internal and kinetic) and boundary work at constant mass";
      // Note that T*der(s) = der(u) + p*der(v).
      output Q.Power Edot_AE(stateSelect=StateSelect.never) = (chemical.mu +
        actualStream(chemical.sT) - h + (actualStream(chemical.phi)*
        actualStream(chemical.phi) - phi*phi)*Data.m/2)*chemical.Ndot + (
        physical.mu + actualStream(physical.sT) - h + (actualStream(physical.phi)
        *actualStream(physical.phi) - phi*phi)*Data.m/2)*physical.Ndot if
        environment.analysis
        "Relative rate of energy (internal, flow, and kinetic) due to phase change and reaction";
      output Q.Power Edot_DE(stateSelect=StateSelect.never) = inert.translational.phi
        *inert.translational.mPhidot + inert.thermal.Qdot + inertDalton.phi*
        inertDalton.mPhidot + inertDalton.Qdot if environment.analysis
        "Rate of diffusion of energy from other configurations";
      output Q.Power Edot_AT(stateSelect=StateSelect.never) = sum((Data.h(faces[
        j, :].T, p_faces[j, :]) - {h,h})*Ndot_faces[j, :] + sum((faces[j, :].phi[
        cartWrap(cartTrans[i] - cartFaces[j] + 1)] .^ 2 - fill(phi[i]^2, 2))*
        Ndot_faces[j, :]*Data.m/2 for i in 1:n_trans) for j in 1:n_faces) if
        environment.analysis
        "Relative rate of energy (internal, flow, and kinetic) due to advective transport";
      output Q.Power Edot_DT(stateSelect=StateSelect.never) = sum(sum(faces[j,
        :].phi[cartWrap(cartTrans[i] - cartFaces[j] + 1)]*faces[j, :].mPhidot[
        cartWrap(cartTrans[i] - cartFaces[j] + 1)] for i in 1:n_trans) for j
         in 1:n_faces) + sum(faces.Qdot) if environment.analysis
        "Rate of diffusion of energy from other subregions";
      // Note:  The structure of the problem should not change if these
      // auxiliary variables are included (hence StateSelect.never).

      Connectors.ChemicalSpecies chemical(
        final n_trans=n_trans,
        mu(start=g_IC, final fixed=false),
        phi(start=phi_IC[cartTrans],each final fixed=false),
        sT(start=h_IC - g_IC, final fixed=false)) "Connector for reactions"
        annotation (Placement(transformation(extent={{-10,10},{10,30}}),
            iconTransformation(extent={{-29,60},{-49,80}})));
      Connectors.Physical physical(
        final formula=Data.formula,
        final n_trans=n_trans,
        mu(start=g_IC, fixed=false),
        phi(final start=phi_IC[cartTrans],each final fixed=false),
        sT(final start=h_IC - g_IC, final fixed=false))
        "Connector for phase change" annotation (Placement(transformation(
              extent={{-30,-10},{-10,10}}), iconTransformation(extent={{-59,28},
                {-79,48}})));
      Connectors.Inert inert(
        final n_trans=n_trans,
        translational(phi(final start=phi_IC[cartTrans], each final fixed=false)),

        thermal(T(final start=T_IC, final fixed=false)))
        "Connector to directly couple velocity or temperature with other species"
        annotation (Placement(transformation(extent={{10,-10},{30,10}}),
            iconTransformation(extent={{59,-28},{79,-48}})));

      Connectors.InertDalton inertDalton(
        final n_trans=n_trans,
        V(
          min=0,
          final start=V_IC,
          final fixed=false),
        p(final start=p_IC, final fixed=false),
        phi(start=phi_IC[cartTrans]),
        T(start=T_IC))
        "Connector for translational and thermal diffusive exchange, with additivity of pressure"
        annotation (Placement(transformation(extent={{-10,-30},{10,-10}}),
            iconTransformation(extent={{29,-60},{49,-80}})));
      Connectors.Face faces[n_faces, Side](
        rho(each start=rho_IC),
        Ndot(start=outerProduct(I_IC[cartFaces], {1,-1})),
        phi(start={fill({phi_IC[cartWrap(cartFaces[i] + orientation - 1)] for
              orientation in Orientation}, 2) for i in 1:n_faces}),
        mPhidot(each start=0),
        T(each start=T_IC),
        Qdot(each start=0))
        "Connectors to transport material, translational momentum, and thermal energy through the boundaries"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}}),
            iconTransformation(extent={{-10,-10},{10,10}})));

      // Geometric parameters
    protected
      outer parameter Q.Length L[Axis] "Length" annotation (missingInnerMessage
          ="This model should be used within a subregion model.
");
      outer parameter Q.Area A[Axis] "Cross-sectional area" annotation (
          missingInnerMessage="This model should be used within a subregion model.
");
      outer parameter Q.Length Lprime[Axis]
        "Effective cross-sectional area per length" annotation (
          missingInnerMessage="This model should be used within a phase model.
      ");
      outer parameter Boolean inclTrans[Axis]
        "true, if each component of translational momentum is included"
        annotation (missingInnerMessage="This model should be used within a subregion model.
");
      outer parameter Boolean inclFaces[Axis]
        "true, if each pair of faces is included" annotation (
          missingInnerMessage="This model should be used within a subregion model.
");
      outer parameter Boolean inclRot[3]
        "true, if each axis of rotation has all its tangential faces included"
        annotation (missingInnerMessage="This model should be used within a subregion model.
");
      // Note:  The size is also Axis, but it can't be specified here due to
      // an error in Dymola 7.4 (failure in check of Phase models).
      outer parameter Integer n_trans
        "Number of components of translational momentum" annotation (
          missingInnerMessage="This model should be used within a subregion model.
");
      outer parameter Integer cartTrans[:]
        "Cartesian-axis indices of the components of translational momentum"
        annotation (missingInnerMessage="This model should be used within a subregion model.
");
      outer parameter Integer cartFaces[:]
        "Cartesian-axis indices of the pairs of faces" annotation (
          missingInnerMessage="This model should be used within a subregion model.
");
      outer parameter Integer cartRot[:]
        "Cartesian-axis indices of the components of rotational momentum"
        annotation (missingInnerMessage="This model should be used within a subregion model.
");
      // Note:  The size of cartTrans, cartFaces, and cartRot is n_trans,
      // but it can't be specified here due to an error in Dymola 7.4.
      outer parameter Integer transCart[3]
        "Translational-momentum-component indices of the Cartesian axes"
        annotation (missingInnerMessage="This model should be used within a subregion model.
");
      // Note:  The size is also Axis, but it can't be specified here due to
      // an error in Dymola 7.4 (failure in check of Phase models).
      outer parameter Integer facesCart[Axis]
        "Face-pair indices of the Cartesian axes" annotation (
          missingInnerMessage="This model should be used within a subregion model.
");
      final parameter Boolean upstream[Axis]={upstreamX,upstreamY,upstreamZ}
        "true, if each Cartesian axis uses upstream discretization"
        annotation (HideResult=true);
      final parameter Conservation consTrans[Axis]={consTransX,consTransY,
          consTransZ} "Formulation of the translational conservation equations"
        annotation (HideResult=true);
      final parameter InitTranslational initTrans[Axis]={initTransX,initTransY,
          initTransZ} "Initialization methods for translational momentum"
        annotation (HideResult=true);

      // Additional aliases (for common terms)
      Q.Potential g0(
        nominal=U.V,
        final start=g_IC,
        final fixed=false,
        stateSelect=StateSelect.never) "Gibbs potential at reference pressure";
      Q.Force faces_mPhidot[n_faces, Side, 2]
        "Directly calculated shear forces";

      outer Conditions.Environment environment "Environmental conditions";

    initial equation
      // Check the initialization methods.
      assert(initMaterial <> initEnergy or initMaterial == InitScalar.None or
        consMaterial == Conservation.Steady or consEnergy == Conservation.Steady,
        "The initialization methods for material and energy must be different (unless None).");
      if not Data.isCompressible then
        assert(initMaterial <> InitScalar.Pressure and initMaterial <>
          InitScalar.PressureSS or consMaterial == Conservation.IC or
          consMaterial == Conservation.Steady, "The material is incompressible,
yet the initialization method for material of " + Data.formula +
          " involves pressure.");
        assert(initEnergy <> InitScalar.Pressure and initEnergy <> InitScalar.PressureSS
           or consEnergy == Conservation.IC or consEnergy == Conservation.Steady,
          "The material is incompressible,
yet the initialization method for energy of " + Data.formula +
          " involves pressure.");
        if not Data.hasThermalExpansion then
          assert(initMaterial <> InitScalar.Density and initMaterial <>
            InitScalar.DensitySS or consMaterial == Conservation.IC or
            consMaterial == Conservation.Steady, "The material is isochoric,
yet the initialization method for material of " + Data.formula +
            " involves density.");
          assert(initEnergy <> InitScalar.Density and initEnergy <> InitScalar.DensitySS
             or consEnergy == Conservation.IC or consEnergy == Conservation.Steady,
            "The material is isochoric,
yet the initialization method for energy of " + Data.formula +
            " involves density.");
        end if;
      end if;

      // Material
      if consMaterial == Conservation.IC then
        // Ensure that a condition is selected since the state is prescribed.
        assert(initMaterial <> InitScalar.None, "The material state of " + Data.formula
           + " is prescribed,
yet its condition is not defined.  Choose any condition besides None.");
      elseif consMaterial == Conservation.Dynamic then
        // Initialize since there's a time-varying state.
        if initMaterial == InitScalar.Amount then
          N = N_IC;
        elseif initMaterial == InitScalar.AmountSS then
          der(N) = 0;
        elseif initMaterial == InitScalar.Density then
          1/v = rho_IC;
        elseif initMaterial == InitScalar.DensitySS then
          der(1/v) = 0;
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
          chemical.mu = g_IC;
        elseif initMaterial == InitScalar.PotentialGibbsSS then
          der(chemical.mu) = 0;
          // Else there's no initial equation since
          // initMaterial == InitScalar.None or
          // consMaterial == Conservation.Steady.
        end if;
      end if;

      // Velocity
      for i in 1:n_trans loop
        if consTrans[cartTrans[i]] == Conservation.IC then
          // Ensure that a condition is selected since the state is
          // prescribed.
          assert(initTrans[cartTrans[i]] <> InitTranslational.None,
            "The state for the " + {"x","y","z"}[cartTrans[i]] +
            "-axis component of translational momentum of " + Data.formula + " is prescribed,
yet its condition is not defined.  Choose any condition besides None.");
        elseif consTrans[cartTrans[i]] == Conservation.Dynamic then
          // Initialize since there's a time-varying state.
          if initTrans[cartTrans[i]] == InitTranslational.Velocity then
            phi[i] = phi_IC[cartTrans[i]];
          elseif initTrans[cartTrans[i]] == InitTranslational.VelocitySS then
            der(phi[i]) = 0;
          elseif initTrans[cartTrans[i]] == InitTranslational.Current then
            I[i] = I_IC[cartTrans[i]];
          elseif initTrans[cartTrans[i]] == InitTranslational.CurrentSS then
            der(I[i]) = 0;
            // Else there's no initial equation since
            // initTrans[cartTrans[i]] == InitTranslational.None or
            // consTrans[cartTrans[i]] == Conservation.Steady.
          end if;
        end if;
      end for;

      // Energy
      if consEnergy == Conservation.IC then
        // Ensure that a condition is selected since the state is prescribed.
        assert(initEnergy <> InitScalar.None, "The energy state of " + Data.formula
           + " is prescribed,
yet its condition is not defined.  Choose any condition besides None.");
      elseif consEnergy == Conservation.Dynamic then
        // Initialize since there's a time-varying state.
        if initEnergy == InitScalar.Amount then
          N = N_IC;
        elseif initEnergy == InitScalar.AmountSS then
          der(N) = 0;
        elseif initMaterial == InitScalar.Density then
          1/v = rho_IC;
        elseif initMaterial == InitScalar.DensitySS then
          der(1/v) = 0;
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
          chemical.mu = g_IC;
        elseif initEnergy == InitScalar.PotentialGibbsSS then
          der(chemical.mu) = 0;
          // Else there's no initial equation since
          // initEnergy == InitScalar.None or
          // consEnergy == Conservation.Steady.
        end if;
      end if;

    equation
      // Aliases (only to clarify and simplify other equations)
      T = inert.thermal.T;
      p = inertDalton.p;
      V = inertDalton.V;
      v*N = V;
      M = Data.m*N;
      phi = inert.translational.phi;
      I .* L[cartTrans] = N*phi;
      p_faces = {{Data.p_Tv(faces[i, side].T, 1/faces[i, side].rho) for side
         in Side} for i in 1:n_faces};
      Ndot_faces = faces.Ndot + faces.rho .* {{faces[i, Side.n].phi[Orientation.normal],
        -faces[i, Side.p].phi[Orientation.normal]}*A[cartFaces[i]] for i in 1:
        n_faces};

      // Thermodynamic correlations
      if invertEOS then
        p = Data.p_Tv(T, v);
      else
        v = Data.v_Tp(T, p);
      end if;
      h = Data.h(T, p);
      s = Data.s(T, p);
      g0 = Data.g(T, Data.p0);

      // Diffusive exchange
      chemical.mu = h - chemical.sT "Reaction (rate equation elsewhere)";
      if tauprime > Modelica.Constants.small then
        tauprime*physical.Ndot = N*(exp((physical.mu - g0)/T) - exp((chemical.mu
           - g0)/T)) "Phase change";
      else
        physical.mu = chemical.mu;
        // This avoids nonlinear equations when tauprime=0.  Dymola 7.4 can't
        // derive it symbolically from the previous equation.
      end if;
      mu*inertDalton.mPhidot = N*(inertDalton.phi - phi)
        "Translational momentum";
      nu*inertDalton.Qdot = N*(inertDalton.T - T) "Thermal energy";

      // Properties upon outflow due to reaction and phase change
      chemical.phi = phi;
      physical.phi = phi;
      chemical.sT = s*T;
      physical.sT = chemical.sT;

      // Diffusive transport
      for i in 1:n_faces loop
        for side in Side loop
          // Material (central difference)
          eta*faces[i, side].Ndot = Lprime[cartFaces[i]]*(faces[i, side].rho -
            1/v)*2;
          // **Try upstream again, update doc.

          // Translational momentum
          beta*faces[i, side].mPhidot[Orientation.normal] = Lprime[cartFaces[i]]
            *(faces[i, side].phi[Orientation.normal] - (if inclTrans[cartFaces[
            i]] then phi[transCart[cartFaces[i]]] else 0))*(if inclTrans[
            cartFaces[i]] and upstream[cartFaces[i]] then 1 + exp(-inSign(side)
            *I[transCart[cartFaces[i]]]*beta*Data.m/(2*Lprime[cartFaces[i]]))
             else 2) "Normal";
          zeta*faces_mPhidot[i, side, Orientation.following - 1] = Nu_Phi[
            cartFaces[i]]*Lprime[cartFaces[i]]*(faces[i, side].phi[Orientation.following]
             - (if inclTrans[cartWrap(cartFaces[i] + 1)] then phi[transCart[
            cartWrap(cartFaces[i] + 1)]] else 0))*(if inclTrans[cartFaces[i]]
             and upstream[cartFaces[i]] then 1 + exp(-inSign(side)*I[transCart[
            cartFaces[i]]]*zeta*Data.m/(2*Lprime[cartFaces[i]])) else 2)
            "1st transverse";
          zeta*faces_mPhidot[i, side, Orientation.preceding - 1] = Nu_Phi[
            cartFaces[i]]*Lprime[cartFaces[i]]*(faces[i, side].phi[Orientation.preceding]
             - (if inclTrans[cartWrap(cartFaces[i] - 1)] then phi[transCart[
            cartWrap(cartFaces[i] - 1)]] else 0))*(if inclTrans[cartFaces[i]]
             and upstream[cartFaces[i]] then 1 + exp(-inSign(side)*I[transCart[
            cartFaces[i]]]*zeta*Data.m/(2*Lprime[cartFaces[i]])) else 2)
            "2nd transverse";

          // Thermal energy
          theta*faces[i, side].Qdot = Nu_Q*Lprime[cartFaces[i]]*(faces[i, side].T
             - T)*(if inclTrans[cartFaces[i]] and upstream[cartFaces[i]] then 1
             + exp(-inSign(side)*I[transCart[cartFaces[i]]]*theta*Data.c_v(T, p)
            /(2*Lprime[cartFaces[i]])) else 2);
        end for;

        // Direct mapping of transverse forces (calculated above)
        if not inclRot[cartWrap(cartFaces[i] - 1)] then
          faces[i, :].mPhidot[Orientation.following] = faces_mPhidot[i, :,
            Orientation.following - 1];
          // Else the force must be mapped for zero torque (below).
        end if;
        if not inclRot[cartWrap(cartFaces[i] + 1)] then
          faces[i, :].mPhidot[Orientation.preceding] = faces_mPhidot[i, :,
            Orientation.preceding - 1];
          // Else the force must be mapped for zero torque (below).
        end if;
      end for;

      // Zero-torque mapping of transverse forces
      for axis in cartRot loop
        4*cat(1,
              faces[facesCart[cartWrap(axis + 1)], :].mPhidot[Orientation.following],
              faces[facesCart[cartWrap(axis - 1)], :].mPhidot[Orientation.preceding])
          = {{3,1,L[cartWrap(axis - 1)]/L[cartWrap(axis + 1)],-L[cartWrap(axis
           - 1)]/L[cartWrap(axis + 1)]},{1,3,-L[cartWrap(axis - 1)]/L[cartWrap(
          axis + 1)],L[cartWrap(axis - 1)]/L[cartWrap(axis + 1)]},{L[cartWrap(
          axis + 1)]/L[cartWrap(axis - 1)],-L[cartWrap(axis + 1)]/L[cartWrap(
          axis - 1)],3,1},{-L[cartWrap(axis + 1)]/L[cartWrap(axis - 1)],L[
          cartWrap(axis + 1)]/L[cartWrap(axis - 1)],1,3}}*cat(
              1,
              faces_mPhidot[facesCart[cartWrap(axis + 1)], :, Orientation.following
             - 1],
              faces_mPhidot[facesCart[cartWrap(axis - 1)], :, Orientation.preceding
             - 1]);
      end for;

      // Material dynamics
      if consMaterial == Conservation.IC then
        // Apply the IC forever (material not conserved).
        if initMaterial == InitScalar.Amount then
          N = N_IC;
        elseif initMaterial == InitScalar.AmountSS then
          der(N) = 0;
        elseif initMaterial == InitScalar.Density then
          1/v = rho_IC;
        elseif initMaterial == InitScalar.DensitySS then
          der(1/v) = 0;
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
          chemical.mu = g_IC;
        else
          // if initMaterial == InitScalar.PotentialGibbsSS then
          der(chemical.mu) = 0;
          // Note:  initMaterial == InitScalar.None can't occur due to an
          // assertion.
        end if;
      else
        (if consMaterial == Conservation.Dynamic then der(N)/U.s else 0) =
          chemical.Ndot + physical.Ndot + sum(Ndot_faces)
          "Material conservation";
      end if;

      // Translational dynamics
      for j in 1:n_trans loop
        if consTrans[cartTrans[j]] == Conservation.IC then
          // Apply the IC forever (translational momentum isn't conserved along
          // this axis).
          if initTrans[cartTrans[j]] == InitTranslational.Velocity then
            phi[j] = phi_IC[cartTrans[j]];
          elseif initTrans[cartTrans[j]] == InitTranslational.VelocitySS then
            der(phi[j]) = 0;
          elseif initTransX == InitTranslational.Current then
            I[j] = I_IC[cartTrans[j]];
          else
            // if initTrans[cartTrans[j]] == InitTranslational.CurrentSS then
            der(I[j]) = 0;
            // Note:  initTrans[cartTrans[j]] == InitTranslational.None can't
            // occur due to an assertion.
          end if;
        else
          M*((if consTrans[cartTrans[j]] == Conservation.Dynamic then der(phi[j])
            /U.s else 0) + environment.a[cartTrans[j]]) + N*Data.z*environment.E[
            cartTrans[j]] + (if inclFaces[cartTrans[j]] then Delta(p_faces[
            facesCart[cartTrans[j]], :])*A[cartTrans[j]] else 0) = Data.m*((
            actualStream(chemical.phi) - phi) .* chemical.Ndot + (actualStream(
            physical.phi) - phi) .* physical.Ndot)[j] + inert.translational.mPhidot[
            j] + inertDalton.mPhidot[j] + sum((faces[i, :].phi[cartWrap(
            cartTrans[j] - cartFaces[i] + 1)] - {phi[j],phi[j]})*Ndot_faces[i,
            :]*Data.m + Sigma(faces[i, :].mPhidot[cartWrap(cartTrans[j] -
            cartFaces[i] + 1)]) for i in 1:n_faces)
            "Conservation of translational momentum";
          // Note:  Dymola 7.4 (Dassl integrator) runs better with this intensive
          // form of the balance (M*der(phi) = ... rather than der(M*phi) = ...).
          // Note:  In Dymola 7.4 it isn't possible to individually index the
          // components of advective exchange, e.g.,
          // (actualStream(chemical.phi) .* chemical.Ndot)[i]
          // must be used instead of
          // actualStream(chemical.phi[i])*chemical.Ndot[i].
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
          1/v = rho_IC;
        elseif initMaterial == InitScalar.DensitySS then
          der(1/v) = 0;
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
          chemical.mu = g_IC;
        else
          // if initEnergy == InitScalar.PotentialGibbsSS then
          der(chemical.mu) = 0;
          // Note:  initEnergy == InitScalar.None can't occur due to an
          // assertion.
        end if;
      else
        (if consEnergy == Conservation.Dynamic then (N*T*der(s) + M*phi*der(phi))
          /U.s else 0) = (chemical.mu + actualStream(chemical.sT) - h + (
          actualStream(chemical.phi)*actualStream(chemical.phi) - phi*phi)*Data.m
          /2)*chemical.Ndot + (physical.mu + actualStream(physical.sT) - h + (
          actualStream(physical.phi)*actualStream(physical.phi) - phi*phi)*Data.m
          /2)*physical.Ndot + inert.translational.phi*inert.translational.mPhidot
           + inert.thermal.Qdot + inertDalton.phi*inertDalton.mPhidot +
          inertDalton.Qdot + sum((Data.h(faces[i, :].T, p_faces[i, :]) - {h,h})
          *Ndot_faces[i, :] + sum((faces[i, :].phi[cartWrap(cartTrans[j] -
          cartFaces[i] + 1)] .^ 2 - fill(phi[j]^2, 2))*Ndot_faces[i, :]*Data.m/
          2 + faces[i, :].phi[cartWrap(cartTrans[j] - cartFaces[i] + 1)]*faces[
          i, :].mPhidot[cartWrap(cartTrans[j] - cartFaces[i] + 1)] for j in 1:
          n_trans) for i in 1:n_faces) + sum(faces.Qdot)
          "Conservation of energy";
      end if;
      annotation (
        defaultComponentPrefixes="replaceable",
        Documentation(info="<html>
    <p>This model is based on the following fixed assumptions:
    <ol>
       <li>All faces are rectangular.
       <li>The material is orthorhombic.  This implies that a gradient which induces diffusion
       along an axis does not induce diffusion along axes orthogonal to it
       [<a href=\"modelica://FCSys.UsersGuide.References\">Bejan2006</a>,
       pp. 691&ndash;692].</li>
       <li>The coordinate system (x, y, z) is aligned with the principle
       axes of transport.  For example if the species is stratified, then the
       layers must be parallel to one of the planes in the rectilinear
       grid.</li>
       <li>The factors that may cause anisotropic behavior (<b><i>k</i></b>)
          are common to material, translational, and thermal transport.</li>
       <li>There is no radiative heat transfer.</li>
       <li>Rotational momentum is not exchanged, transported, or stored.</li>
       <li>For the purpose of the material, translational momentum, and energy balances, the
       cross sectional areas of the faces are assumed to be the full cross-sectional
       areas of the subregion.  If multiple phases are present, then the areas are
       actually smaller.</li>
    </ol>
    Other assumptions are optional via the parameters.</p>

    <p>Figure 1 shows how instances of
    <a href=\"modelica://FCSys.Subregions.Species\">Species</a> models (derived from this
    model) are
    connected within a <a href=\"modelica://FCSys.Subregions\">Subregion</a>.  The
    generalized resistances (<i>R</i>) affect the phase change rate, forces, and heat flow rates
    associated with differences in activity, velocity, and temperature (respectively) between
    each species and a common node.  These exchange processes are diffusive.

    <p>In general, the resistances are included within the
    <a href=\"modelica://FCSys.Subregions.Species\">Species</a> models.  For reactions, however,
    the rate equation is more complex and is included in the
    <a href=\"modelica://FCSys.Subregions.Reaction\">Reaction</a> model.</p>

    <p>Translational momentum and thermal energy are advected as material is exchanged
    due to phase change or reactions.  This occurs at the velocity (&phi;) and specific entropy-temperature
    product (<i>sT</i>) of the reactants (source configurations), where the reactant/product designation
    depends on the current conditions.</p>

    <p>The advective exchange is modeled using <code>stream</code> connectors
    (<a href=\"modelica://FCSys.Connectors.Physical\">Physical</a> and
    <a href=\"modelica://FCSys.Connectors.ChemicalSpecies\">ChemicalSpecies</a>).
  The rate of advection of translational momentum is the
  product of the velocity of the source (&phi;) and the mass flow rate
  (<i>m</i><i>N&#775;</i>).  The rate of thermal advection is the
  specific entropy-temperature product of the source (<i>sT</i>) times the rate of
  material exchange
  (<i>N&#775;</i>).  If there are multiple sources, then
  their contributions are additive.  If there are multiple sinks, then
  translational momentum is split on a mass basis and the thermal stream is split
  on a particle-number basis.</p>

    <p align=center><img src=\"modelica://FCSys/resources/documentation/Subregions/Species/Species/Exchange.png\">
<br>Figure 1:  Exchange of a quantity (translational momentum or thermal energy) among species
    (A, B, and C) within a subregion.</p>

    <p>Figure 2 shows how
    instances of the same type of <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model
    are connected between neighboring instances of a
    <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a>.
    Material, translational momentum, and thermal energy are transported by both advection and diffusion.
    Upstream discretization is applied if it is enabled via the <code>upstreamX</code>,
    etc. parameters.  Similar to exchange, the transport resistances are included inside the
    <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model.</p>

    <p align=center><img src=\"modelica://FCSys/resources/documentation/Subregions/Species/Species/Transport.png\">
<br>Figure 2:  Transport of a quantity associated with the same species
    between subregions (1 and 2).</p>

    <p>The <a href=\"modelica://FCSys.Subregions.Species\">Species</a> instances
    within a <a href=\"modelica://FCSys.Subregions.Phases\">Phase</a> are combined by Dalton's law of
    partial pressures (see the
    <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector), as shown
    in Figure 3a.  The pressures are additive, and each species is assumed to exist at the
    total extensive volume of the phase.  Within a <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a>,
    the <a href=\"modelica://FCSys.Subregions.Phases\">Phases</a> are combined by Amagat's law of partial volumes
    (see the <a href=\"modelica://FCSys.Subregions.Volume\">Volume</a> model), as shown
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

    <p>Notes regarding the parameters:
    <ul>
    <li>Here (and in the rest of <a href=\"modelica://FCSys\">FCSys</a>), the <i>specific</i>
    adjective means that the following extensive quantity is divided by particle number.
    (<i>Massic</i> indicates a quantity divided by mass.)</li>
    <li>In general, if material resistivity, dynamic compressibility, fluidity, or thermal resistivity is zero, then
    it should be set as <code>final</code> so that index reduction may be performed.
    If two <a href=\"modelica://FCSys.Subregions.Species\">Species</a> instances
    are connected through their <code>inertDalton</code> connectors or faces
    and both have zero generalized resistivities for a
    quantity, then index reduction [<a href=\"modelica://FCSys.Subregions.Species\">Mattsson1993B</a>] is necessary.</li>
    <li>Even if an initialization parameter is not selected for explicit use,
    it may be used a guess value.</li>
    <li>The area fill factor (<b><i>k</i></b>) can be used to account for the effects of porosity and tortousity
    on the rate of transport.  It may reflect anisotropic properties, since it is a vector with independent components
    for each axis.
    By default, its components are unity.  The area fill factor should be adjusted directly with effective
    area and inversely with effective length.
    It affects all of the diffusive transport rates (material, translational, and
    thermal) by the same factor.</li>
    <li>If <code>Conservation.IC</code> is used for a state (via
    <code>consMaterial</code>, <code>consTransX</code>, <code>consTransY</code>,
    <code>consTransZ</code>, or <code>consEnergy</code>),
    then the associated initial condition (IC) will be applied forever instead of the
    corresponding conservation equation.
    If <code>consMaterial</code>, <code>consTransX</code>, <code>consTransY</code>, or <code>consTransZ</code> is
    <code>Conservation.IC</code>, then there may be a secondary effect on the energy conservation equation
    and thus temperature.
    In that case, it may help to set <code>consEnergy</code> to <code>Conservation.IC</code> so that
    the energy conservation equation is not imposed.</li>
    <li>If <code>consTransX</code>, <code>consTransY</code>, or <code>consTransZ</code> is
    <code>Conservation.Steady</code>, then the derivative of the corresponding component of velocity
    is treated as zero and removed from the translational momentum balance.  If <code>consEnergy</code> is
    <code>Conservation.Steady</code>, then <code>T*der(s) + M*phi*der(phi)</code> is treated as
    zero and removed from the energy balance.</li>
    <li>If a component of velocity is not included (via the outer <code>inclTrans[:]</code> parameter
    which maps to <code>{inclTransX, inclTransY, inclTransZ}</code> in the
    <a href=\"modelica://FCSys.Subregions.BaseClasses.EmptySubregion\">Subregion</a> model), then it
    is taken to be zero in each translational transport equation.  However, the corresponding forces
    in the <code>faces</code> connector array are not included in the momentum or energy balances.
    If it is necessary to set a component of velocity to zero but still include it in the energy balance, then
    set the corresponding component of <code>phi_IC</code> to zero and <code>consTransX<code>,
    <code>consTransY<code>, or <code>consTransZ</code> to <code>Conservation.IC</code>.</li>
    <li>If a subregion does not contain any compressible species, then pressure must be prescribed.
    Set <code>consMaterial</code> to <code>Conservation.IC</code> and <code>initMaterial</code>
    to <code>InitScalar.Pressure</code> for one of the species.</li>
    <li>The <code>start</code> values of the initial conditions for pressure and temperature
    (<i>p</i><sub>IC</sub> and <i>T</i><sub>IC</sub>) are the global default pressure and
    temperature (via the <code>outer</code> instance of the <a href=\"modelica://FCSys.Conditions.Environment\">Environment</a> model).
    The <code>start</code> values of the initial conditions for
    other intensive properties (&rho;<sub>IC</sub>, <i>h</i><sub>IC</sub>, and
    <i>g</i><sub>IC</sub>) are related to the initial pressure and temperature
    by the characteristics of the species.  The <code>start</code> value of the
    initial condition for the extensive volume (<i>V</i><sub>IC</sub>) is the volume of the
    subregion.  The <code>start</code> value for particle number (<i>N</i><sub>IC</sub>)
    is related to it via the material characteristics and the initial pressure and temperature.
    In order to apply other values for any of these initial conditions,
    it may be necessary to do so before translating the model.</li>
    <li>Upstream discretization may be applied to translational and thermal transport
    using (<code>upstreamX=true</code>, etc.).  Otherwise, the central difference
    scheme is used.  The central difference scheme
    is always used for material diffusion.</li>
    <li>If <code>invertEOS</code> is <code>true</code>, then the equation of state is implemented with pressure
    as a function of temperature and specific volume.  Otherwise, specific volume is a function of temperature
    and pressure.</li>
    </p>

    <p>In the <code>faces</code> connector array, the transverse translational flow (<i>m</i>&Phi;dot) is only the 
    force due to diffusion.  Translational advection is calculated from the velocity and the material current.
    The thermal flow (<i>Q&#775;</i>) is only the rate of heat transfer due to diffusion.  The advection of 
    thermal energy is determined from the thermodynamic state at the boundary and the material current.</p>

    <p>In evaluating the dynamics of a phase, it is typically assumed that all of the species
    exist at the same velocity and temperature.  The translational and thermal time constants
    are usually much shorter than the time span of interest due to the very small coupling
    resistances.  If this is the case, connect the <code>inert</code>
    connectors of the species.  This will reduce the index of the problem.</p>

    <p>For the variables that relate to transport,
    the first index is the axis and the second index is the side.  The sides
    are ordered from negative to positive, according to the
    <a href=\"modelica://FCSys.BaseClasses.Side\">Side</a> enumeration.
    Velocity and force are additionally indexed by
    the orientation of the momentum with respect to the face.
    The orientations are ordered in Cartesian space starting with the normal axis,
    according to the
    <a href=\"modelica://FCSys.BaseClasses.Orientation\">Orientation</a> enumeration.</p>
    </html>"),
        Diagram(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.1)),
        Icon(graphics={
            Rectangle(
              extent={{-98,80},{98,120}},
              fillPattern=FillPattern.Solid,
              fillColor={255,255,255},
              pattern=LinePattern.None),
            Ellipse(
              extent={{-80,80},{80,-80}},
              lineColor={127,127,127},
              pattern=LinePattern.Dash,
              fillColor={225,225,225},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-98,80},{98,120}},
              textString="%name",
              lineColor={0,0,0})}));
    end Species;

    package BaseClasses "Base classes (generally not for direct use)"

      extends Modelica.Icons.BasesPackage;

      type Conservation = enumeration(
          IC "Initial condition imposed forever (no conservation)",
          Steady "Steady (conservation with steady state)",
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
          Current "Prescribed advective current",
          CurrentSS "Steady-state advective current")
        "Methods of initializing translational momentum";

    end BaseClasses;

  end Species;

  model Reaction "Electrochemical reaction"
    extends FCSys.BaseClasses.Icons.Names.Top2;

    parameter Q.Area A=100*U.cm^2 "Area";
    parameter Q.CurrentAreicAbsolute J_0=1e-3*U.A/U.cm^2
      "<html>Exchange current density (<i>J</i><sub>0</sub>)</html>";
    parameter Q.NumberAbsolute alpha(max=1) = 0.5
      "<html>Charge transfer coefficient (&alpha;)</html>";

    Q.Number Pe(start=0) "Peclet number";

    Connectors.ChemicalReaction chemical(final n_trans=n_trans)
      "Common connector for the reaction" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}}), iconTransformation(extent={{-10,-10},{
              10,10}})));

  protected
    outer parameter Integer n_trans
      "Number of components of translational momentum" annotation (
        missingInnerMessage="This model should be used within a subregion model.
");

  equation
    Pe*chemical.T = chemical.mu - 25.87*U.V;
    // **temp offset
    chemical.Ndot = J_0*A*(exp(alpha*Pe) - exp((alpha - 1)*Pe))
      "Butler-Volmer equation (reaction rate)";

    // Conservation (without storage)
    zeros(n_trans) = chemical.mPhidot "Translational momentum";
    0 = chemical.Qdot_A "Advected thermal energy";
    0 = chemical.Qdot_D + 0*(chemical.mu - 25.87*U.V)*chemical.Ndot
      "Diffused thermal energy (balances heat generation)";
    // **temp offset, zero
    annotation (Documentation(info="<html><p>This model establishes the rate of an electrochemical reaction
  using the Butler-Volmer equation.  It includes a static energy balance with heat generation.
  The heat is rejected to <code>chemical.Qdot_D</code>, independently of the
  thermal stream from the reactants to the products (<code>chemical.Qdot_A</code>).</p>

    <p>The exchange current density (<i>J</i><sub>0</sub>) is the exchange current per unit geometric area (not per
    unit of catalyst surface area).</p>

    <p></p></html>"), Icon(graphics={Ellipse(
            extent={{-40,40},{40,-40}},
            lineColor={127,127,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.Dash)}));
  end Reaction;

  model ChemicalExchange
    "Model to impose stoichiometry and properties on a species in a reaction"
    extends FCSys.BaseClasses.Icons.Names.Top1;

    parameter Q.Number n "Stoichiometric coefficient";
    parameter Q.MassSpecific m "Specific mass";

    // Auxiliary variables (for analysis)
    output Q.Velocity phi_actualStream[n_trans](each stateSelect=StateSelect.never)
       = actualStream(species.phi) if environment.analysis
      "Velocity of the actual stream";
    output Q.PotentialAbsolute sT_actualStream(stateSelect=StateSelect.never)
       = actualStream(species.sT) if environment.analysis
      "Specific entropy-temperature product of the actual stream";

    Connectors.ChemicalReaction reaction(final n_trans=n_trans)
      "Reaction-side connector" annotation (Placement(transformation(extent={{
              10,-10},{30,10}}), iconTransformation(extent={{10,-10},{30,10}})));
    Connectors.ChemicalSpecies species(final n_trans=n_trans)
      "Species-side connector" annotation (Placement(transformation(extent={{-30,
              -10},{-10,10}}), iconTransformation(extent={{-30,-10},{-10,10}})));

  protected
    outer parameter Integer n_trans
      "Number of components of translational momentum" annotation (
        missingInnerMessage="This model should be used within a subregion model.
   ");
    outer Conditions.Environment environment "Environmental conditions";

  equation
    // Intensive properties
    reaction.mu = n*species.mu;
    reaction.phi = species.phi;
    reaction.sT = species.sT;

    // Conservation (without storage)
    0 = species.Ndot + n*reaction.Ndot "Material";
    zeros(n_trans) = reaction.mPhidot + m*actualStream(species.phi)*species.Ndot
      "Translational momentum";
    0 = reaction.Qdot_A + actualStream(species.sT)*species.Ndot
      "Advected thermal energy";
    0 = reaction.Qdot_D + (reaction.phi - actualStream(species.phi))*reaction.mPhidot
      /2 "Diffused thermal energy (balances heat generation)";
    annotation (Documentation(info="<html><p>This model is essentially an adapter between the
  <a href=\"modelica://FCSys.Connectors.ChemicalSpecies\">ChemicalSpecies</a> and the
  <a href=\"modelica://FCSys.Connectors.ChemicalReaction\">ChemicalReaction</a> connectors.
  It should be instantiated once for each species in
  a reaction.</p></html>"), Icon(graphics={Line(
            points={{-10,0},{10,0}},
            color={255,195,38},
            smooth=Smooth.None), Text(
            extent={{-100,-20},{100,-40}},
            lineColor={127,127,127},
            textString="%n")}));
  end ChemicalExchange;

  model Depletion "Electrochemical depletion region"
    import FCSys.BaseClasses.Utilities.inSign;
    //extends FCSys.BaseClasses.Icons.Names.Top2;

    parameter Q.Area A=U.cm^2 "Area" annotation (Dialog(group="Geometry"));
    parameter Side side
      "Side of the subregion which is depleted of the carrier"
      annotation (Dialog(group="Geometry"));
    replaceable package Data = Characteristics.BaseClasses.Characteristic
      constrainedby Characteristics.BaseClasses.CharacteristicEOS
      "Characteristic data" annotation (
      Dialog(group="Material properties"),
      choicesAllMatching=true,
      __Dymola_choicesFromPackage=true,
      Placement(transformation(extent={{-60,40},{-40,60}}), iconTransformation(
            extent={{-10,90},{10,110}})));
    parameter Boolean transSubstrate=false
      "Pass translational momentum through the substrate" annotation (choices(
          __Dymola_checkBox=true), Dialog(tab="Assumptions", compact=true));
    parameter Boolean thermalSubstrate=false
      "Reject heat generated from the reaction into the substrate" annotation (
        choices(__Dymola_checkBox=true), Dialog(tab="Assumptions", compact=true));

    Connectors.ChemicalReaction chemical(final n_trans=n_trans)
      "Chemical connector" annotation (Placement(transformation(extent={{-10,-10},
              {10,10}}), iconTransformation(extent={{-10,-10},{10,10}})));
    Connectors.Inert inert(final n_trans=n_trans) "Connection to the substrate"
      annotation (Placement(transformation(extent={{30,-10},{50,10}}),
          iconTransformation(extent={{90,-10},{110,10}})));
    Connectors.Face face "Boundary of the minority region" annotation (
        Placement(transformation(extent={{-50,-10},{-30,10}}),
          iconTransformation(extent={{-110,-10},{-90,10}})));

  protected
    outer parameter Integer n_trans
      "Number of components of translational momentum" annotation (
        missingInnerMessage="This model should be used within a subregion model.
   ");
    outer parameter Integer cartTrans[:]
      "Cartesian-axis indices of the components of translational momentum"
      annotation (missingInnerMessage="This model should be used within a subregion model.
   ");

  equation
    // Electrical potential and electrostatic force
    Data.z*face.rho*A*chemical.mu = face.mPhidot[1];

    // Translational conditions
    if transSubstrate then
      chemical.phi = inert.translational.phi
        "Products produced at the velocity of the substrate";
    else
      chemical.mPhidot = zeros(n_trans)
        "Translational momentum passed directly from reactants to products";
    end if;

    // Thermal conditions
    if thermalSubstrate then
      chemical.T = inert.thermal.T
        "Substrate sets the temperature of the reaction and receives heat generated by the reaction";
    else
      chemical.Qdot_D = 0 "Substrate is adiabatic w.r.t. reaction";
    end if;
    chemical.Qdot_A = 0 "No thermal energy into the advected stream";

    // Conservation at the far edge of the depletion region (without storage)
    0 = face.Ndot - inSign(side)*face.rho*A*face.phi[1] "Material";
    0 = face.mPhidot[1] - inSign(side)*Data.p_Tv(face.T, 1/face.rho)*A
      "Normal translational momentum";
    zeros(2) = face.mPhidot[2:3] "Transverse translational momentum";
    0 = face.Qdot "Energy";

    // Conservation at the reaction site (without storage)
    for i in 1:n_trans loop
      0 = chemical.mPhidot[i] + inert.translational.mPhidot[i] + (if cartTrans[
        i] == 0 then 0 else face.mPhidot[cartTrans[i]])
        "Translational momentum";
      // The conditional expression is necessary to pass the check since
      // cartTrans is empty by default.
    end for;
    0 = chemical.mu*chemical.Ndot + chemical.phi*chemical.mPhidot + inert.thermal.Qdot
       + chemical.Qdot_D + inert.translational.phi*inert.translational.mPhidot
      "Energy (excluding terms which are zero above)";
    annotation (Documentation(info="<html>
    <p>This model introduces the electrical potential associated with the force on a minority charge carrier in half
    of an electrochemical double layer.  It should be instantiated for each of the two minority regions.  The
    electrical potential is added to the net chemical potential
    of the reaction, resulting in the electrochemical potential.</p>

    <p>If <code>transSubstrate</code> is <code>true</code>, then the translational momentum of the
    reactants is passed to the substrate through the <code>inert</code>
    connector and the products are produced at the velocity of the substrate (typically
    zero).  If it is <code>false</code>, then translational momentum is passed
    directly from the reactants to the products.  If <code>transSubstrate</code> is
    <code>true</code> for both of the depletion regions, then the velocities of the
    substrates will be coupled and this may cause a structural singularity.</p>

    <p>If <code>thermalSubstrate</code> is <code>true</code>, then the reaction occurs at the temperature
    of the substrate and the heat generated by the reaction is rejected to the substrate.  If it is
    <code>false</code>, then the substrate and the reaction are thermally independent.
    Typically, <code>thermalSubstrate</code> should be <code>true</code> for one but not both of the depletion regions.
    It it is <code>true</code> for both of the depletion regions, then the temperatures of the
    substrates will be coupled and this may cause a structural singularity.
    If it is <code>false</code> for both of the depletion regions, then the reaction temperature will be
    undefined.</p>
    </html>"), Icon(graphics={
          Line(
            points={{-90,0},{-20,0}},
            color={127,127,127},
            smooth=Smooth.None),
          Line(
            points={{-20,30},{-20,-30}},
            color={127,127,127},
            smooth=Smooth.None,
            thickness=0.5),
          Line(
            points={{20,0},{90,0}},
            color={47,107,251},
            smooth=Smooth.None),
          Line(
            points={{20,30},{20,-30}},
            color={47,107,251},
            smooth=Smooth.None,
            thickness=0.5),
          Text(
            extent={{-120,40},{120,80}},
            textString="%name",
            lineColor={0,0,0})}));
  end Depletion;

  model Volume "Model to establish a fixed total volume"
    // extends FCSys.BaseClasses.Icons.Names.Top7;

    outer parameter Q.Volume V "Volume" annotation (missingInnerMessage="This model should be used within a subregion model.
", HideResult=true);
    // The check fails if this is protected in Dymola 7.4, so
    // HideResult=true is used instead.
    parameter Integer n_phases(min=1) = 0 "Number of phases"
      annotation (Dialog(connectorSizing=true));
    Connectors.InertDalton inertDalton[n_phases](each n_trans=n_trans)
      "Connector for translational momentum and thermal energy" annotation (
        Placement(transformation(extent={{60,-80},{80,-60}}),
          iconTransformation(extent={{100,-120},{120,-100}})));

  protected
    outer parameter Integer n_trans
      "Number of components of translational momentum" annotation (
        missingInnerMessage="This model should be used within a subregion model.
   ");

  equation
    if n_phases > 0 then
      if n_phases > 1 then
        // Equal properties (KVL)
        inertDalton[1:n_phases - 1].p = inertDalton[2:n_phases].p "Pressures";
        inertDalton[1:n_phases - 1].phi = inertDalton[2:n_phases].phi
          "Velocities";
        inertDalton[1:n_phases - 1].T = inertDalton[2:n_phases].T
          "Temperatures";
      end if;

      // Conservation (KCL)
      V = sum(inertDalton.V) "Volume (fixed)";
      zeros(n_trans) = sum(inertDalton[i].mPhidot for i in 1:n_phases)
        "Translational momentum (no storage)";
      0 = sum(inertDalton.Qdot) "Energy (no storage)";
    end if;
    annotation (
      Documentation(info="<html><p>This model combines the mixtures connected to its <code>inertDalton[:]</code>
    array using additivity of volume.  An entry is added to this array each time a connection is made
    (using <code>connectorSizing</code>).  Each entry establishes a network of
    <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connectors.
    The species connected to each network are combined using additivity of pressure.</p>

    <p>The concept of \"additivity of volume\" is defined by
    <a href=\"http://en.wikipedia.org/wiki/Amagat's_law\">Amagat's law of partial volumes</a>, which
    states that the partial extensive volumes of the components of a mixture sum to the total
    extensive volume of the mixture [<a href=\"modelica://FCSys.UsersGuide.References\">Bejan2006</a>, p. 194].
    The specific or molar volumes of the species are each evaluated at the temperature and the total pressure of the
    mixture.</p>

    <p>This concept loses its physical meaning once the species are mixed [<a href=\"modelica://FCSys.UsersGuide.References\">Woo1995</a>].
    If the species are truly mixed, then it is impossible to distinguish their particles and thus
    determine their partial volumes.
    Therefore, the concept is only used for distinct phases within the same subregion&mdash;not
    for species within a phase.
    For example, if a system contains a solid phase and a gas phase, then it is assumed that the
    partial volumes of the mixtures are additive and the mixtures exist at the same pressure.  Within
    a phase, the species are mixed according to Dalton's law of partial pressures (see the
    <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector).</p>

    <p>Please also see the
    the documentation for the
    <a href=\"modelica://FCSys.Connectors\">Connectors</a> package.</p></html>"),

      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-160,-160},{160,
              160}}), graphics={Polygon(
            points={{-160,60},{-60,160},{160,160},{160,-60},{60,-160},{-160,-160},
                {-160,60}},
            lineColor={127,127,127},
            smooth=Smooth.None,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.Dash), Text(
            extent={{-160,112},{160,152}},
            textString="%name",
            lineColor={0,0,0})}),
      Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
              100,100}}), graphics));
  end Volume;

  package BaseClasses "Base classes (generally not for direct use)"

    extends Modelica.Icons.BasesPackage;

    partial model EmptySubregion
      "Base model for multi-dimensional, multi-species storage, transport, and exchange"
      import FCSys.BaseClasses.Utilities.cartWrap;
      import FCSys.BaseClasses.Utilities.countTrue;
      import FCSys.BaseClasses.Utilities.enumerate;
      import FCSys.BaseClasses.Utilities.index;
      // extends FCSys.BaseClasses.Icons.Names.Top3;

      // Geometric parameters
      inner parameter Q.Length L[Axis](each min=Modelica.Constants.small) = {U.cm,
        U.cm,U.cm} "<html>Length (<b>L</b>)</html>"
        annotation (Dialog(group="Geometry"));
      final inner parameter Q.Area A[Axis]={L[cartWrap(axis + 1)]*L[cartWrap(
          axis + 2)] for axis in Axis} "Cross-sectional areas";
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

      Connectors.FaceBus xNegative if inclFacesX
        "Negative face along the x axis" annotation (Placement(transformation(
              extent={{-50,-10},{-30,10}}), iconTransformation(extent={{-110,-10},
                {-90,10}})));
      Connectors.FaceBus xPositive if inclFacesX
        "Positive face along the x axis" annotation (Placement(transformation(
              extent={{30,-10},{50,10}}), iconTransformation(extent={{90,-10},{
                110,10}})));
      Connectors.FaceBus yNegative if inclFacesY
        "Negative face along the y axis" annotation (Placement(transformation(
              extent={{-10,-50},{10,-30}}), iconTransformation(extent={{-10,-110},
                {10,-90}})));
      Connectors.FaceBus yPositive if inclFacesY
        "Positive face along the y axis" annotation (Placement(transformation(
              extent={{-10,30},{10,50}}), iconTransformation(extent={{-10,90},{
                10,110}})));
      Connectors.FaceBus zNegative if inclFacesZ
        "Negative face along the z axis" annotation (Placement(transformation(
              extent={{10,10},{30,30}}), iconTransformation(extent={{40,40},{60,
                60}})));
      Connectors.FaceBus zPositive if inclFacesZ
        "Positive face along the z axis" annotation (Placement(transformation(
              extent={{-30,-30},{-10,-10}}), iconTransformation(extent={{-60,-60},
                {-40,-40}})));

    protected
      parameter Boolean hasSpecies "true, if any species are included";
      final inner parameter Boolean inclTrans[Axis]={inclTransX,inclTransY,
          inclTransZ}
        "true, if each component of translational momentum is included";
      final inner parameter Boolean inclFaces[Axis]={inclFacesX,inclFacesY,
          inclFacesZ} "true, if each pairs of faces is included";
      final inner parameter Boolean inclRot[Axis]={inclFacesY and inclFacesZ,
          inclFacesZ and inclFacesX,inclFacesX and inclFacesY}
        "true, if each axis of rotation has all its tangential faces included";
      final inner parameter Integer n_trans=countTrue(inclTrans)
        "Number of components of translational momentum";
      final inner parameter Integer n_faces=countTrue(inclFaces)
        "Number of pairs of faces";
      final inner parameter Integer cartTrans[n_trans]=index(inclTrans)
        "Cartesian-axis indices of the components of translational momentum";
      final inner parameter Integer cartFaces[n_faces]=index(inclFaces)
        "Cartesian-axis indices of the pairs of faces";
      final inner parameter Integer cartRot[:]=index(inclRot)
        "Cartesian-axis indices of the components of rotational momentum";
      final inner parameter Integer transCart[Axis]=enumerate(inclTrans)
        "Translational-momentum-component indices of the Cartesian axes";
      final inner parameter Integer facesCart[Axis]=enumerate(inclFaces)
        "Face-pair indices of the Cartesian axes";

      Volume volume if hasSpecies "Model to establish a fixed total volume"
        annotation (Placement(transformation(extent={{-16,-16},{16,16}})));
      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="subregion",
        Documentation(info="<html>
  <p>At least one component of translational momentum must be included.
  All of the components are included by default.</p>

    <p>At least one pair of faces must be included.
  All of the faces are included by default.</p>

  <p>This model should be extended to include the appropriate phases and reactions.</p>
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
            Text(
              extent={{-100,56},{100,96}},
              textString="%name",
              lineColor={0,0,0})}));

    end EmptySubregion;

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
