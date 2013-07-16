within FCSys;
package Regions "3D arrays of discrete, interconnected subregions"
  package Examples "Examples"
    extends Modelica.Icons.ExamplesPackage;
    model FPToFP "Test one flow plate to the other"

      extends Modelica.Icons.Example;
      extends Modelica.Icons.UnderConstruction;

      parameter Q.Length L_y[:]={U.m}
        "<html>Lengths along the channel (<i>L</i><sub>y</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      parameter Q.Length L_z[:]={5*U.mm}
        "<html>Lengths across the channel (<i>L</i><sub>z</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      final parameter Integer n_y=size(L_y, 1)
        "Number of regions along the channel" annotation (HideResult=true);
      final parameter Integer n_z=size(L_z, 1)
        "Number of regions across the channel" annotation (HideResult=true);
      AnFPs.AnFP anFP(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{-70,-10},{-50,10}})));

      AnGDLs.AnGDL anGDL(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));

      AnCLs.AnCL anCL(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));

      PEMs.PEM PEM(
        final L_y=L_y,
        final L_z=L_z,
        Subregion(ionomer('H+'(initTransX=InitTranslational.None))))
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      CaCLs.CaCL caCL(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{10,-10},{30,10}})));

      CaGDLs.CaGDL caGDL(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{30,-10},{50,10}})));

      CaFPs.CaFP caFP(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{50,-10},{70,10}})));

      Conditions.ByConnector.FaceBus.Single.FaceBusEfforts BC1[n_y, n_z](each
          graphite('incle-'=true, 'e-'(redeclare function normalSpec =
                Conditions.ByConnector.Face.Single.TranslationalNormal.currentDensity,
              redeclare Modelica.Blocks.Sources.Ramp normalSource(
              height=-U.A/U.cm^2,
              duration=100.1,
              startTime=0.1)))) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={-84,0})));
      Conditions.ByConnector.FaceBus.Single.FaceBusEfforts BC2[n_y, n_z](each
          graphite('incle-'=true, 'e-'(redeclare function normalSpec =
                Conditions.ByConnector.Face.Single.TranslationalNormal.force)))
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={84,0})));
      Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC3[anFP.n_x, n_z](
          each gas(inclH2=true, inclH2O=true)) annotation (Placement(
            transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-60,-24})));

      Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC4[anFP.n_x, n_z](
          each gas(inclH2=true, inclH2O=true)) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-60,24})));

      Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC5[caFP.n_x, n_z](
          each gas(
          inclH2O=true,
          inclN2=true,
          inclO2=true)) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={60,-24})));
      Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC6[caFP.n_x, n_z](
          each gas(
          inclH2O=true,
          inclN2=true,
          inclO2=true)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={60,24})));
      inner Conditions.Environment environment(analysis=true)
        "Environmental conditions"
        annotation (Placement(transformation(extent={{70,52},{90,72}})));

    equation
      connect(anFP.xPositive, anGDL.xNegative) annotation (Line(
          points={{-50,6.10623e-16},{-50,6.10623e-16}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
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
          points={{10,6.10623e-16},{10,6.10623e-16},{10,6.10623e-16},{10,
              6.10623e-16}},
          color={0,0,240},
          smooth=Smooth.None,
          thickness=0.5));

      connect(caCL.xPositive, caGDL.xNegative) annotation (Line(
          points={{30,6.10623e-16},{30,6.10623e-16}},
          color={0,0,240},
          smooth=Smooth.None,
          thickness=0.5));

      connect(caGDL.xPositive, caFP.xNegative) annotation (Line(
          points={{50,6.10623e-16},{50,6.10623e-16}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));

      connect(BC3.face, anFP.yNegative) annotation (Line(
          points={{-60,-20},{-60,-10}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(BC4.face, anFP.yPositive) annotation (Line(
          points={{-60,20},{-60,10}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(BC5.face, caFP.yNegative) annotation (Line(
          points={{60,-20},{60,-10}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      connect(BC6.face, caFP.yPositive) annotation (Line(
          points={{60,20},{60,10}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        Commands(file="Resources/Scripts/Dymola/Regions.Examples.FPToFP.mos"
            "Regions.Examples.FPToFP.mos"),
        experiment(
          StopTime=20,
          Tolerance=1e-06,
          Algorithm="Dassl"),
        Diagram(graphics));
      connect(BC1.face, anFP.xNegative) annotation (Line(
          points={{-80,-1.34539e-15},{-76,-1.34539e-15},{-76,6.10623e-16},{-70,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(caFP.xPositive, BC2.face) annotation (Line(
          points={{70,6.10623e-16},{76,6.10623e-16},{76,1.23436e-15},{80,
              1.23436e-15}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

    end FPToFP;

    model GDLToGDL "Test one GDL to the other"

      extends Modelica.Icons.Example;
      extends Modelica.Icons.UnderConstruction;

      parameter Q.Length L_y[:]={U.m}
        "<html>Lengths along the channel (<i>L</i><sub>y</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      parameter Q.Length L_z[:]={5*U.mm}
        "<html>Lengths across the channel (<i>L</i><sub>z</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      final parameter Integer n_y=size(L_y, 1)
        "Number of regions along the channel" annotation (HideResult=true);
      final parameter Integer n_z=size(L_z, 1)
        "Number of regions across the channel" annotation (HideResult=true);
      AnGDLs.AnGDL anGDL(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));
      AnCLs.AnCL anCL(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));

      PEMs.PEM PEM(
        final L_y=L_y,
        final L_z=L_z,
        Subregion(ionomer('H+'(initTransX=InitTranslational.None))))
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      CaCLs.CaCL caCL(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{10,-10},{30,10}})));

      CaGDLs.CaGDL caGDL(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{30,-10},{50,10}})));

      Conditions.ByConnector.FaceBus.Single.FaceBusEfforts BC1[n_y, n_z](each
          gas(
          inclH2=true,
          inclH2O=true,
          H2(materialSource(y=(environment.p - environment.p_H2O)/environment.T)),

          H2O(materialSource(y=environment.p_H2O/environment.T))),each graphite(
            'incle-'=true, 'e-'(redeclare function normalSpec =
                Conditions.ByConnector.Face.Single.TranslationalNormal.currentDensity,
              redeclare Modelica.Blocks.Sources.Ramp normalSource(
              height=-U.A/U.cm^2,
              duration=100.1,
              startTime=0.1)))) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={-64,0})));

      Conditions.ByConnector.FaceBus.Single.FaceBusEfforts BC2[n_y, n_z](each
          gas(
          inclH2O=true,
          inclO2=true,
          H2O(materialSource(y=environment.p_H2O/environment.T)),
          O2(materialSource(y=environment.p_O2/environment.T))), each graphite(
            'incle-'=true, 'e-'(redeclare function normalSpec =
                Conditions.ByConnector.Face.Single.TranslationalNormal.force)))
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={64,0})));

      inner Conditions.Environment environment(analysis=true)
        "Environmental conditions"
        annotation (Placement(transformation(extent={{70,50},{90,70}})));

    equation
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
          points={{10,6.10623e-16},{10,6.10623e-16},{10,6.10623e-16}},
          color={0,0,240},
          smooth=Smooth.None,
          thickness=0.5));

      connect(caCL.xPositive, caGDL.xNegative) annotation (Line(
          points={{30,6.10623e-16},{30,6.10623e-16}},
          color={0,0,240},
          smooth=Smooth.None,
          thickness=0.5));

      annotation (
        Commands(file="Resources/Scripts/Dymola/Regions.Examples.GDLToGDL.mos"
            "Regions.Examples.GDLToGDL.mos"),
        experiment(
          StopTime=30,
          Tolerance=1e-06,
          Algorithm="Dassl"),
        Diagram(graphics));
      connect(BC1.face, anGDL.xNegative) annotation (Line(
          points={{-60,-1.34539e-15},{-56,-1.34539e-15},{-56,6.10623e-16},{-50,
              6.10623e-16}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));

      connect(caGDL.xPositive, BC2.face) annotation (Line(
          points={{50,6.10623e-16},{56,6.10623e-16},{56,1.23436e-15},{60,
              1.23436e-15}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));

    end GDLToGDL;

    model CLToCL "Test one catalyst layer to the other"

      extends Modelica.Icons.Example;
      extends Modelica.Icons.UnderConstruction;

      parameter Q.Length L_y[:]={U.m}
        "<html>Lengths along the channel (<i>L</i><sub>y</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      parameter Q.Length L_z[:]={5*U.mm}
        "<html>Lengths across the channel (<i>L</i><sub>z</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      final parameter Integer n_y=size(L_y, 1)
        "Number of regions along the channel" annotation (HideResult=true);
      final parameter Integer n_z=size(L_z, 1)
        "Number of regions across the channel" annotation (HideResult=true);

      AnCLs.AnCL anCL(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));

      PEMs.PEM PEM(
        final L_y=L_y,
        final L_z=L_z,
        Subregion(ionomer('H+'(initTransX=InitTranslational.None))))
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      CaCLs.CaCL caCL(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{10,-10},{30,10}})));

      Conditions.ByConnector.FaceBus.Single.FaceBusEfforts BC1[n_y, n_z](each
          gas(
          inclH2=true,
          inclH2O=true,
          H2(materialSource(y=(environment.p - environment.p_H2O)/environment.T)),

          H2O(materialSource(y=environment.p_H2O/environment.T))),each graphite(
            'incle-'=true, 'e-'(redeclare function normalSpec =
                Conditions.ByConnector.Face.Single.TranslationalNormal.force,
              redeclare Modelica.Blocks.Sources.Ramp normalSource(
              height=-U.A/U.cm^2,
              duration=100.1,
              offset=-0.001*U.A/U.cm^2,
              startTime=0.1)))) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={-44,0})));

      Conditions.ByConnector.FaceBus.Single.FaceBusEfforts BC2[n_y, n_z](each
          gas(
          inclH2O=true,
          inclN2=true,
          inclO2=true,
          H2O(materialSource(y=environment.p_H2O/environment.T)),
          N2(materialSource(y=(environment.p - environment.p_H2O - environment.p_O2)
                  /environment.T)),
          O2(materialSource(y=environment.p_O2/environment.T))), each graphite(
            'incle-'=true, 'e-'(redeclare function normalSpec =
                Conditions.ByConnector.Face.Single.TranslationalNormal.force)))
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={44,0})));
      inner Conditions.Environment environment(analysis=true)
        "Environmental conditions"
        annotation (Placement(transformation(extent={{70,50},{90,70}})));

    equation
      connect(anCL.xPositive, PEM.xNegative) annotation (Line(
          points={{-10,6.10623e-16},{-10,6.10623e-16}},
          color={240,0,0},
          smooth=Smooth.None,
          thickness=0.5));
      connect(PEM.xPositive, caCL.xNegative) annotation (Line(
          points={{10,6.10623e-16},{10,6.10623e-16},{10,6.10623e-16},{10,
              6.10623e-16}},
          color={0,0,240},
          smooth=Smooth.None,
          thickness=0.5));

      annotation (Commands(file=
              "Resources/Scripts/Dymola/Regions.Examples.CLToCL.mos"
            "Regions.Examples.CLToCL.mos"), experiment(
          StopTime=25,
          Tolerance=1e-06,
          Algorithm="Dassl"));
      connect(BC1.face, anCL.xNegative) annotation (Line(
          points={{-40,-1.34539e-15},{-40,0},{-30,0},{-30,6.10623e-16}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(BC2.face, caCL.xPositive) annotation (Line(
          points={{40,1.23436e-15},{40,6.10623e-16},{30,6.10623e-16}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
    end CLToCL;

    model AnFP "Test the anode flow plate"

      extends Modelica.Icons.Example;

      parameter Q.Length L_y[:]={U.m}
        "<html>Lengths along the channel (<i>L</i><sub>y</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      parameter Q.Length L_z[:]={5*U.mm}
        "<html>Lengths across the channel (<i>L</i><sub>z</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      final parameter Integer n_y=size(L_y, 1)
        "Number of regions along the channel" annotation (HideResult=true);
      final parameter Integer n_z=size(L_z, 1)
        "Number of regions across the channel" annotation (HideResult=true);

      AnFPs.AnFP anFP(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{-70,-10},{-50,10}})));

      Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC1[n_y, n_z](each
          graphite('inclC+'=true, 'incle-'=true)) annotation (Placement(
            transformation(
            extent={{10,-10},{-10,10}},
            rotation=90,
            origin={-84,0})));
      Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC2[n_y, n_z](
        each gas(inclH2=true, inclH2O=true),
        each graphite('incle-'=true),
        each liquid(inclH2O=true)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-36,0})));

      Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC3[anFP.n_x, n_z](
          each gas(inclH2=true, inclH2O=true)) annotation (Placement(
            transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-60,-24})));
      Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC4[anFP.n_x, n_z](
          each gas(inclH2=true, inclH2O=true)) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-60,24})));
      inner Conditions.Environment environment(analysis=true)
        "Environmental conditions"
        annotation (Placement(transformation(extent={{70,50},{90,70}})));

    equation
      connect(BC1.face, anFP.xNegative) annotation (Line(
          points={{-80,2.54679e-16},{-80,6.10623e-16},{-70,6.10623e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(BC2.face, anFP.xPositive) annotation (Line(
          points={{-40,1.23436e-15},{-40,6.10623e-16},{-50,6.10623e-16}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));

      connect(BC3.face, anFP.yNegative) annotation (Line(
          points={{-60,-20},{-60,-10}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(BC4.face, anFP.yPositive) annotation (Line(
          points={{-60,20},{-60,10}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (experiment(Tolerance=1e-06), Commands(file(ensureSimulated=
                true) = "Resources/Scripts/Dymola/Regions.Examples.AnFP.mos"));
    end AnFP;

    model AnGDL "Test the anode gas diffusion layer"

      extends Modelica.Icons.Example;

      parameter Q.Length L_y[:]={U.m}
        "<html>Lengths along the channel (<i>L</i><sub>y</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      parameter Q.Length L_z[:]={5*U.mm}
        "<html>Lengths across the channel (<i>L</i><sub>z</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      final parameter Integer n_y=size(L_y, 1)
        "Number of regions along the channel" annotation (HideResult=true);
      final parameter Integer n_z=size(L_z, 1)
        "Number of regions across the channel" annotation (HideResult=true);
      AnGDLs.AnGDL anGDL(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));

      Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC1[n_y, n_z](
        each gas(inclH2=true, inclH2O=true),
        each graphite('incle-'=true),
        each liquid(inclH2O=true)) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=90,
            origin={-64,0})));

      Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC2[n_y, n_z](
        each gas(inclH2=true, inclH2O=true),
        each graphite('incle-'=true),
        each liquid(inclH2O=true)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-16,0})));
      inner Conditions.Environment environment(analysis=true)
        "Environmental conditions"
        annotation (Placement(transformation(extent={{70,50},{90,70}})));

    equation
      connect(BC1.face, anGDL.xNegative) annotation (Line(
          points={{-60,2.54679e-16},{-56,2.54679e-16},{-56,6.10623e-16},{-50,
              6.10623e-16}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));

      connect(anGDL.xPositive, BC2.face) annotation (Line(
          points={{-30,6.10623e-16},{-26,6.10623e-16},{-26,1.23436e-15},{-20,
              1.23436e-15}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (experiment(Tolerance=1e-06, StopTime=10), Commands(file(
              ensureSimulated=true) =
            "Resources/Scripts/Dymola/Regions.Examples.AnGDL.mos"));

    end AnGDL;

    model AnCL "Test the anode catalyst layer"

      extends Modelica.Icons.Example;

      parameter Q.Length L_y[:]={U.m}
        "<html>Lengths along the channel (<i>L</i><sub>y</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      parameter Q.Length L_z[:]={5*U.mm}
        "<html>Lengths across the channel (<i>L</i><sub>z</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      final parameter Integer n_y=size(L_y, 1)
        "Number of regions along the channel" annotation (HideResult=true);
      final parameter Integer n_z=size(L_z, 1)
        "Number of regions across the channel" annotation (HideResult=true);

      /*
  output Q.Potential w[1, n_y, n_z]=-anCL.subregions.reaction.chemical.mu/2
    "Reaction potential";
  output Q.Current zI[1, n_y, n_z]=anCL.subregions.graphite.'e-'.chemical.Ndot
    "Electrical current due to reaction";
  output Q.Number zJ_Apercm2[1, n_y, n_z]=(zI*U.cm^2) ./ (anCL.subregions[:, :,
      :].A[Axis.x]*U.A) "Electrical current density, in A/cm2";
  output Q.Power Qdot=sum(anCL.subregions.graphite.'C+'.Edot_DE)
    "Rate of heat generation";
  output Q.Power P=sum(w .* zI) "Electrical power";
  output Q.NumberAbsolute eta=P/(P + Qdot) "Efficiency";
  */
      // **Fix Qdot, P, eta
      inner Conditions.Environment environment(analysis=true)
        "Environmental conditions"
        annotation (Placement(transformation(extent={{70,50},{90,70}})));

      AnCLs.AnCL anCL(
        final L_y=L_y,
        final L_z=L_z,
        Subregion(ionomer('H+'(initMaterial=InitScalar.None))))
        annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));

      Conditions.ByConnector.FaceBus.Single.FaceBusEfforts BC1[n_y, n_z](each
          gas(
          inclH2=true,
          inclH2O=true,
          H2(materialSource(y=(environment.p - environment.p_H2O)/environment.T)),

          H2O(materialSource(y=environment.p_H2O/environment.T))),each graphite(
            'incle-'=true, 'e-'(redeclare function normalSpec =
                Conditions.ByConnector.Face.Single.TranslationalNormal.currentDensity,
              redeclare Modelica.Blocks.Sources.Ramp normalSource(
              height=-U.A/U.cm^2,
              duration=100.1,
              startTime=0.1)))) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={-44,0})));

      Conditions.ByConnector.FaceBus.Single.FaceBusEfforts BC2[n_y, n_z](each
          ionomer('inclH+'=true, 'H+'(redeclare function normalSpec =
                Conditions.ByConnector.Face.Single.TranslationalNormal.force)))
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={4,0})));

    equation
      connect(BC1.face, anCL.xNegative) annotation (Line(
          points={{-40,-1.34539e-15},{-40,6.10623e-16},{-30,6.10623e-16}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));

      connect(anCL.xPositive, BC2.face) annotation (Line(
          points={{-10,6.10623e-16},{-2,6.10623e-16},{-2,1.23436e-15},{
              6.66134e-16,1.23436e-15}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                {100,100}}), graphics),
        experiment(StopTime=110, Tolerance=1e-06),
        Commands(file="Resources/Scripts/Dymola/Regions.Examples.AnCL.mos"
            "Regions.Examples.AnCL.mos"),
        experimentSetupOutput);
    end AnCL;

    model PEM "Test the proton exchange membrane"

      extends Modelica.Icons.Example;

      parameter Q.Length L_y[:]={U.m}
        "<html>Lengths along the channel (<i>L</i><sub>y</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      parameter Q.Length L_z[:]={5*U.mm}
        "<html>Lengths across the channel (<i>L</i><sub>z</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      final parameter Integer n_y=size(L_y, 1)
        "Number of regions along the channel" annotation (HideResult=true);
      final parameter Integer n_z=size(L_z, 1)
        "Number of regions across the channel" annotation (HideResult=true);

      PEMs.PEM PEM(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC1[n_y, n_z](each
          ionomer('inclH+'=true)) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=90,
            origin={-24,0})));

      Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC2[n_y, n_z](each
          ionomer('inclH+'=true)) annotation (Placement(transformation(
            extent={{10,10},{-10,-10}},
            rotation=90,
            origin={24,0})));
      inner Conditions.Environment environment(analysis=true)
        "Environmental conditions"
        annotation (Placement(transformation(extent={{70,50},{90,70}})));

    equation
      connect(BC1.face, PEM.xNegative) annotation (Line(
          points={{-20,2.54679e-16},{-16,2.54679e-16},{-16,6.10623e-16},{-10,
              6.10623e-16}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));

      connect(PEM.xPositive, BC2.face) annotation (Line(
          points={{10,6.10623e-16},{14,6.10623e-16},{14,-3.65701e-16},{20,
              -3.65701e-16}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (experiment(Tolerance=1e-06, StopTime=10), Commands(file=
              "Resources/Scripts/Dymola/Regions.Examples.PEM.mos"
            "Regions.Examples.PEM.mos"));
    end PEM;

    model CaCL "Test the cathode catalyst layer"

      extends Modelica.Icons.Example;

      parameter Q.Length L_y[:]={U.m}
        "<html>Lengths along the channel (<i>L</i><sub>y</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      parameter Q.Length L_z[:]={5*U.mm}
        "<html>Lengths across the channel (<i>L</i><sub>z</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      final parameter Integer n_y=size(L_y, 1)
        "Number of regions along the channel" annotation (HideResult=true);
      final parameter Integer n_z=size(L_z, 1)
        "Number of regions across the channel" annotation (HideResult=true);

      /*
  output Q.Potential w[1, n_y, n_z]=-caCL.subregions.reaction.chemical.mu/2 
    "Reaction potential";
  output Q.Current zI[1, n_y, n_z]=caCL.subregions.graphite.'e-'.chemical.Ndot 
    "Electrical current due to reaction";
  output Q.Number zJ_Apercm2[1, n_y, n_z]=(zI*U.cm^2) ./ (caCL.subregions[:, :,
      :].A[Axis.x]*U.A) "Electrical current density, in A/cm2";
  output Q.Power Qdot=sum(caCL.subregions.graphite.'C+'.Edot_DE) 
    "Rate of heat generation";
  output Q.Power P=sum(w .* zI) "Electrical power";
  output Q.NumberAbsolute eta=P/(P + Qdot) "Efficiency";
  // **Fix Qdot, P, eta
  
  // **Add and use variable for w in Subregion model.
*/

      CaCLs.CaCL caCL(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{10,-10},{30,10}})));
      Conditions.ByConnector.FaceBus.Single.FaceBusEfforts BC1[n_y, n_z](each
          ionomer('inclH+'=true, 'H+'(redeclare function normalSpec =
                Conditions.ByConnector.Face.Single.TranslationalNormal.force)))
        annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={-24,0})));

      Conditions.ByConnector.FaceBus.Single.FaceBusEfforts BC2[n_y, n_z](each
          gas(
          inclH2O=true,
          inclO2=true,
          H2O(materialSource(y=environment.p_H2O/environment.T)),
          O2(materialSource(y=environment.p_O2/environment.T))), each graphite(
            'incle-'=true, 'e-'(redeclare function normalSpec =
                Conditions.ByConnector.Face.Single.TranslationalNormal.currentDensity,
              redeclare Modelica.Blocks.Sources.Ramp normalSource(
              height=-U.A/U.cm^2,
              duration=100.1,
              startTime=0.1)))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={44,0})));

      inner Conditions.Environment environment(analysis=true)
        "Environmental conditions"
        annotation (Placement(transformation(extent={{70,50},{90,70}})));

    equation
      connect(BC1.face, caCL.xNegative) annotation (Line(
          points={{-20,-1.34539e-15},{10,-1.34539e-15},{10,6.10623e-16}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));

      connect(BC2.face, caCL.xPositive) annotation (Line(
          points={{40,1.23436e-15},{40,6.10623e-16},{30,6.10623e-16}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                {100,100}}), graphics),
        experiment(
          StopTime=110,
          Tolerance=1e-06,
          Algorithm="Dassl"),
        Commands(file="Resources/Scripts/Dymola/Regions.Examples.CaCL.mos"
            "Regions.Examples.CaCL.mos"),
        experimentSetupOutput);
    end CaCL;

    model CaGDL "Test the cathode gas diffusion layer"

      extends Modelica.Icons.Example;

      parameter Q.Length L_y[:]={U.m}
        "<html>Lengths along the channel (<i>L</i><sub>y</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      parameter Q.Length L_z[:]={5*U.mm}
        "<html>Lengths across the channel (<i>L</i><sub>z</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      final parameter Integer n_y=size(L_y, 1)
        "Number of regions along the channel" annotation (HideResult=true);
      final parameter Integer n_z=size(L_z, 1)
        "Number of regions across the channel" annotation (HideResult=true);

      CaGDLs.CaGDL caGDL(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{30,-10},{50,10}})));
      Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC1[n_y, n_z](
        each gas(
          inclH2O=true,
          inclN2=true,
          inclO2=true),
        each graphite('incle-'=true),
        each liquid(inclH2O=true)) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=90,
            origin={16,0})));

      Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC2[n_y, n_z](
        each gas(
          inclH2O=true,
          inclN2=true,
          inclO2=true),
        each graphite('incle-'=true),
        each liquid(inclH2O=true)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={64,0})));
      inner Conditions.Environment environment(analysis=true)
        "Environmental conditions"
        annotation (Placement(transformation(extent={{70,50},{90,70}})));

    equation
      connect(BC1.face, caGDL.xNegative) annotation (Line(
          points={{20,2.54679e-16},{20,6.10623e-16},{30,6.10623e-16}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));

      connect(BC2.face, caGDL.xPositive) annotation (Line(
          points={{60,1.23436e-15},{50,1.23436e-15},{50,6.10623e-16}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (experiment(StopTime=10, Tolerance=1e-06), Commands(file(
              ensureSimulated=true) =
            "Resources/Scripts/Dymola/Regions.Examples.CaGDL.mos"));

    end CaGDL;

    model CaFP "Test the cathode flow plate"

      extends Modelica.Icons.Example;

      parameter Q.Length L_y[:]={U.m}
        "<html>Lengths along the channel (<i>L</i><sub>y</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      parameter Q.Length L_z[:]={5*U.mm}
        "<html>Lengths across the channel (<i>L</i><sub>z</sub>)</html>"
        annotation (Dialog(group="Geometry"));
      final parameter Integer n_y=size(L_y, 1)
        "Number of regions along the channel" annotation (HideResult=true);
      final parameter Integer n_z=size(L_z, 1)
        "Number of regions across the channel" annotation (HideResult=true);

      CaFPs.CaFP caFP(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{50,-10},{70,10}})));

      Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC1[n_y, n_z](
        each gas(
          inclH2O=true,
          inclN2=true,
          inclO2=true),
        each graphite('incle-'=true),
        each liquid(inclH2O=true)) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=90,
            origin={36,0})));

      Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC2[n_y, n_z](each
          graphite('inclC+'=true, 'incle-'=true)) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={84,0})));
      Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC3[caFP.n_x, n_z](
          each gas(
          inclH2O=true,
          inclN2=true,
          inclO2=true)) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={60,-24})));
      Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC4[caFP.n_x, n_z](
          each gas(
          inclH2O=true,
          inclN2=true,
          inclO2=true)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={60,24})));
      inner Conditions.Environment environment(analysis=true)
        "Environmental conditions"
        annotation (Placement(transformation(extent={{70,50},{90,70}})));

    equation
      connect(BC1.face, caFP.xNegative) annotation (Line(
          points={{40,2.54679e-16},{42,2.54679e-16},{42,6.10623e-16},{50,
              6.10623e-16}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));

      connect(caFP.xPositive, BC2.face) annotation (Line(
          points={{70,6.10623e-16},{80,6.10623e-16},{80,1.23436e-15}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));

      connect(BC3.face, caFP.yNegative) annotation (Line(
          points={{60,-20},{60,-10}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      connect(BC4.face, caFP.yPositive) annotation (Line(
          points={{60,20},{60,10}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (experiment(StopTime=10, Tolerance=1e-06), Commands(file(
              ensureSimulated=true) =
            "Resources/Scripts/Dymola/Regions.Examples.CaFP.mos"));
    end CaFP;

  end Examples;
  extends Modelica.Icons.Package;
  import Modelica.Media.IdealGases.Common.SingleGasesData;

  package AnFPs "Anode flow plates"
    extends Modelica.Icons.Package;

    model AnFP "Anode flow plate"
      import FCSys.BaseClasses.Utilities.Polynomial;
      // extends FCSys.BaseClasses.Icons.Names.Top4;
      import FCSys.Characteristics.'C+'.Graphite;

      extends Region(
        L_x={8*U.mm},
        L_y={2*U.m},
        L_z={5*U.mm},
        final inclFacesX=true,
        final inclFacesY=true,
        inclFacesZ=false,
        redeclare replaceable model Subregion =
            FCSys.Subregions.SubregionNoIonomer (
            inclTransX=true,
            inclTransY=true,
            inclTransZ=false,
            gas(
              reduceVel=true,
              reduceTemp=true,
              k_T=fill(epsilon^(3/2), 3),
              inclH2=true,
              inclH2O=true,
              H2(
                initTransX=InitTranslational.None,
                initTransY=InitTranslational.None,
                initTransZ=InitTranslational.None,
                initEnergy=InitScalar.None,
                p_IC=environment.p - environment.p_H2O),
              H2O(consTransX=Conservation.IC, p_IC=environment.p_H2O)),
            graphite(
              reduceTemp=true,
              k_T=fill((1 - epsilon)^(3/2), 3),
              'inclC+'=true,
              'incle-'=true,
              'C+'(V_IC=V - epsilonV, theta=U.m*U.K/(95*U.W)),
              'e-'(
                consTransY=Conservation.IC,
                initTransX=InitTranslational.None,
                initEnergy=InitScalar.None,
                rho_IC=17842.7*U.C/U.cc,
                sigma=U.S/(1.470e-3*U.cm))),
            liquid(
              inclH2O=true,
              k_T=fill(epsilon^(3/2), 3),
              H2O(consTransX=Conservation.IC,V_IC=Modelica.Constants.eps*U.cc))))
        annotation (IconMap(primitivesVisible=false));

      // See the documentation layer of Subregions.Phases.BaseClasses.EmptyPhase
      // regarding the settings of k_T for each phase.

      // TODO:  Find a better solution for liquid.H2O.V_IC than
      // Modelica.Constants.eps and apply it to other layers too.

      parameter Q.NumberAbsolute epsilon(nominal=1) = 0.1
        "<html>Volumetric porosity (&epsilon;)</html>"
        annotation (group="Geometry");

      // Auxiliary variables (for analysis)
      output Q.NumberAbsolute RH[n_x, n_y, n_z](
        each stateSelect=StateSelect.never,
        each displayUnit="%") =
        Modelica.Media.Air.MoistAir.saturationPressureLiquid(subregions.gas.H2O.T
        /U.K) ./ subregions.gas.dalton.p if environment.analysis and
        hasSubregions "Relative humidity";

    protected
      final parameter Q.Volume epsilonV=epsilon*V "Fluid volume";

      outer Conditions.Environment environment "Environmental conditions";

      // Thermal resistivity of some other flowplate materials [Incropera2002,
      // pp. 905 & 907]:
      //                                Stainless steel
      //                                ---------------------------------------------------------------------------------
      //             Aluminium (pure)   AISI 302             AISI 304              AISI 316              AISI 347
      //             ----------------   ------------------   ------------------    -----------------    -----------------
      //                        theta              theta                theta                 theta                theta
      //             c_p*U.kg   *U.W    c_p*U.kg   *U.W      c_p*U.kg   *U.W       c_p*U.kg   *U.W      c_p*U.kg   *U.W
      //             *U.K       /(U.m   *U.K       /(U.m     *U.K       /(U.m      *U.K       /(U.m     *U.K       /(U.m
      //     T/K     /(U.J*m)   *U.K)   /(U.J*m)   *U.K)     /(U.J*m)   *U.K)      /(U.J*m)   *U.K)     /(U.J*m)   *U.K)
      //     ----    --------   -----   --------   ------    --------   -------    --------   ------    --------   ------
      //     100     482        1/302                        272        1/9.2
      //     200     798        1/237                        402        1/12.6
      //     300     903        1/237   480        1/15.1    477        1/14.9     468        1/13.4    480        1/14.2
      //     400     949        1/240   512        1/17.3    515        1/16.6     504        1/15.2    513        1/15.8
      //     600     1033       1/231   559        1/20.0    557        1/19.8     550        1/18.3    559        1/18.9
      //     800     1146       1/218   585        1/22.8    582        1/22.6     576        1/21.3    585        1/21.9
      //     1000                       606        1/25.4    611        1/25.4     602        1/24.2    606        1/24.7

      // Electrical resistivities:
      //     Aluminium
      //     (http://en.wikipedia.org/wiki/Electrical_resistivity, 2008):
      //         2.82e-8 ohm.m
      //     Graphite (http://hypertextbook.com/facts/2004/AfricaBelgrave.shtml):
      //         7.837e-6 to 41e-6 ohm.m
      //     Copper (http://en.wikipedia.org/wiki/Electrical_resistivity, 2008):
      //         1.72e-8 ohm.m
      //     Stainless steel AISI 304
      //     (http://hypertextbook.com/facts/2006/UmranUgur.shtml, 2008):
      //         6.897e-7 ohm.m
      annotation (Documentation(info="<html>
<p>This model represents the anode flow plate of a PEMFC.
The x axis extends from the anode to the cathode.
The y axis extends along the length of the channel. The model is
bidirectional, so that either <code>yNegative</code> or <code>yPositive</code> can be
used as the inlet. The z axis extends across the width of the channel.</p>

<p>By default, the cross-sectional area in the xy plane is 100 cm<sup>2</sup>.</p>

<p>See <a href=\"modelica://FCSys.Species.'C+'.Graphite.Fixed\">Subregions.Species.'C+'.Graphite.Fixed</a>
regarding the default specific heat capacity.  The default thermal resistivity
of the carbon (<code>theta = U.m*U.K/(95*U.W)</code>) and the
electrical conductivity (<code>sigma=U.S/(1.470e-3*U.cm)</code>)
are that of Entegris/Poco Graphite AXF-5Q
[<a href=\"modelica://FCSys.UsersGuide.References\">Entegris2012</a>].
There is additional data in the
text layer of this model.</p>

<p>For more information, see the
    <a href=\"modelica://FCSys.Regions.Region\">Region</a> model.</p>
</html>"), Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.1), graphics={
            Rectangle(
              extent={{-100,60},{100,100}},
              fillColor={255,255,255},
              visible=not inclFacesY,
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Rectangle(
              extent={{-76.648,66.211},{-119.073,52.0689}},
              fillPattern=FillPattern.HorizontalCylinder,
              rotation=45,
              fillColor={135,135,135},
              origin={111.017,77.3801},
              pattern=LinePattern.None,
              lineColor={95,95,95}),
            Rectangle(
              extent={{-20,40},{0,-60}},
              lineColor={95,95,95},
              fillPattern=FillPattern.VerticalCylinder,
              fillColor={135,135,135}),
            Polygon(
              points={{20,0},{42,0},{42,80},{-42,80},{-42,0},{-20,0},{-20,40},{
                  0,60},{20,60},{20,0}},
              smooth=Smooth.None,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Polygon(
              points={{20,0},{42,0},{42,-80},{-42,-80},{-42,0},{-20,0},{-20,-60},
                  {0,-60},{20,-40},{20,0}},
              smooth=Smooth.None,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Rectangle(extent={{-20,40},{0,-60}}, lineColor={0,0,0}),
            Polygon(
              points={{-20,40},{0,60},{20,60},{0,40},{-20,40}},
              lineColor={0,0,0},
              smooth=Smooth.None),
            Polygon(
              points={{20,60},{0,40},{0,-60},{20,-40},{20,60}},
              lineColor={0,0,0},
              fillColor={95,95,95},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{16,48},{4,36},{4,32},{14,42},{14,36},{4,26},{4,12},{16,
                  24},{16,28},{6,18},{6,24},{16,34},{16,48}},
              smooth=Smooth.None,
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None,
              lineColor={0,0,0}),
            Polygon(
              points={{16,28},{4,16},{4,12},{14,22},{14,16},{4,6},{4,-8},{16,4},
                  {16,8},{6,-2},{6,4},{16,14},{16,28}},
              smooth=Smooth.None,
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None,
              lineColor={0,0,0}),
            Polygon(
              points={{16,8},{4,-4},{4,-8},{14,2},{14,-4},{4,-14},{4,-28},{16,-16},
                  {16,-12},{6,-22},{6,-16},{16,-6},{16,8}},
              smooth=Smooth.None,
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None,
              lineColor={0,0,0}),
            Polygon(
              points={{16,-12},{4,-24},{4,-28},{14,-18},{14,-24},{4,-34},{4,-48},
                  {16,-36},{16,-32},{6,-42},{6,-36},{16,-26},{16,-12}},
              smooth=Smooth.None,
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None,
              lineColor={0,0,0}),
            Line(
              points={{10,0},{100,0}},
              color={240,0,0},
              visible=inclFacesX,
              thickness=0.5),
            Line(
              points={{-20,0},{-100,0}},
              color={127,127,127},
              visible=inclFacesX,
              thickness=0.5),
            Line(
              points={{0,-60},{0,-100}},
              color={240,0,0},
              visible=inclFacesY,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{20,20},{50,50}},
              color={253,52,56},
              visible=inclFacesZ,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{-50,-50},{-10,-10}},
              color={253,52,56},
              visible=inclFacesZ,
              smooth=Smooth.None,
              thickness=0.5),
            Ellipse(
              extent={{-4,52},{4,48}},
              lineColor={135,135,135},
              fillColor={240,0,0},
              visible=inclFacesY,
              fillPattern=FillPattern.Sphere),
            Line(
              points={{0,100},{0,50}},
              color={240,0,0},
              visible=inclFacesY,
              smooth=Smooth.None,
              thickness=0.5),
            Text(
              extent={{-100,60},{100,100}},
              textString="%name",
              visible=not inclFacesY,
              lineColor={0,0,0})}));

    end AnFP;

  end AnFPs;

  package AnGDLs "Anode gas diffusion layers"
    extends Modelica.Icons.Package;

    model AnGDL "Anode gas diffusion layer"
      // extends FCSys.BaseClasses.Icons.Names.Top4;

      // Note:  Extensions of AnGDL should be placed directly in the AnGDLs
      // package rather than subpackages (e.g., by manufacturer) so that
      // __Dymola_choicesFromPackage can be used.  Dymola 7.4 launches the
      // parameter dialogs too slowly when __Dymola_choicesAllMatching is
      // used.

      extends Region(
        L_x={0.3*U.mm},
        L_y={2*U.m},
        L_z={5*U.mm},
        final inclFacesX=true,
        inclFacesY=false,
        inclFacesZ=false,
        redeclare replaceable model Subregion =
            FCSys.Subregions.SubregionNoIonomer (
            inclTransX=true,
            inclTransY=false,
            inclTransZ=false,
            gas(
              reduceVel=true,
              reduceTemp=true,
              k_T=fill(epsilon^(3/2), 3),
              inclH2=true,
              inclH2O=true,
              H2(
                initTransX=InitTranslational.None,
                initTransY=InitTranslational.None,
                initTransZ=InitTranslational.None,
                initEnergy=InitScalar.None,
                p_IC=environment.p - environment.p_H2O),
              H2O(consTransX=Conservation.IC, p_IC=environment.p_H2O)),
            graphite(
              reduceTemp=true,
              k_T=fill((1 - epsilon)^(3/2), 3),
              'inclC+'=true,
              'incle-'=true,
              'C+'(V_IC=V - epsilonV, theta=U.m*U.K/(1.18*U.W)),
              'e-'(
                sigma=40*U.S/(12*U.cm),
                initTransX=InitTranslational.None,
                initEnergy=InitScalar.None)),
            liquid(
              inclH2O=true,
              k_T=fill(epsilon^(3/2), 3),
              H2O(V_IC=Modelica.Constants.eps*U.cc,consTransX=Conservation.IC))))
        annotation (IconMap(primitivesVisible=false));
      // **temp excluded C+
      // See the documentation layer of Subregions.Phases.BaseClasses.EmptyPhase
      // regarding the settings of k_T for each phase.

      parameter Q.NumberAbsolute epsilon(nominal=1) = 0.4
        "<html>Volumetric porosity (&epsilon;)</html>"
        annotation (group="Geometry");

      // Auxiliary variables (for analysis)
      output Q.NumberAbsolute RH[n_x, n_y, n_z](
        each stateSelect=StateSelect.never,
        each displayUnit="%") =
        Modelica.Media.Air.MoistAir.saturationPressureLiquid(subregions.gas.H2O.T
        /U.K) ./ subregions.gas.dalton.p if environment.analysis and
        hasSubregions "Relative humidity";

    protected
      final parameter Q.Volume epsilonV=epsilon*V "Fluid volume";

      outer Conditions.Environment environment "Environmental conditions";
      annotation (Documentation(info="<html>
<p>This model represents the anode gas diffusion layer of a PEMFC.
The x axis extends from the anode to the cathode.
The y axis extends along the length of the channel and
the z axis extends across the width of the channel.</p>

<p>By default, the cross-sectional area in the xy plane is 100 cm<sup>2</sup>.</p>

<p>The default porosity (<code>epsilon = 0.4</code>) is
of a compressed GDL according to [<a href=\"modelica://FCSys.UsersGuide.References\">Bernardi1992</a>, p. 2483, Table 3].
  The default thermal conductivity of the carbon (<code>theta = U.m*U.K/(1.18*U.W)</code>)
   represents a compressed SGL Sigracet 10 BA gas diffusion layer
  [<a href=\"modelica://FCSys.UsersGuide.References\">Nitta2008</a>].  The default
  electrical conductivity
  is also for SGL Carbon Group Sigracet&reg; 10 BA [<a href=\"modelica://FCSys.UsersGuide.References\">SGL2007</a>].</p>

<p>For more information, see the
    <a href=\"modelica://FCSys.Regions.Region\">Region</a> model.</p></html>"),
          Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.1), graphics={
            Rectangle(
              extent={{-100,60},{100,100}},
              fillColor={255,255,255},
              visible=not inclFacesY,
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Rectangle(
              extent={{-78.7855,18.6813},{-50.5004,-23.7455}},
              lineColor={64,64,64},
              fillColor={127,127,127},
              rotation=-45,
              fillPattern=FillPattern.VerticalCylinder,
              origin={42.5001,11.0805}),
            Rectangle(
              extent={{-40,40},{0,-60}},
              lineColor={64,64,64},
              fillColor={127,127,127},
              fillPattern=FillPattern.VerticalCylinder),
            Polygon(
              points={{20,0},{42,0},{42,80},{-42,80},{-42,0},{-20,0},{-20,40},{
                  0,60},{20,60},{20,0}},
              smooth=Smooth.None,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Polygon(
              points={{20,0},{42,0},{42,-80},{-42,-80},{-42,0},{-20,0},{-20,-60},
                  {0,-60},{20,-40},{20,0}},
              smooth=Smooth.None,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Polygon(
              points={{0,40},{20,60},{20,-40},{0,-60},{0,40}},
              lineColor={0,0,0},
              smooth=Smooth.None,
              fillPattern=FillPattern.Solid,
              fillColor={64,64,64}),
            Rectangle(extent={{-20,40},{0,-60}}, lineColor={0,0,0}),
            Polygon(
              points={{0,60},{20,60},{0,40},{-20,40},{0,60}},
              lineColor={0,0,0},
              smooth=Smooth.None),
            Line(
              points={{-20,0},{-100,0}},
              color={240,0,0},
              visible=inclFacesX,
              thickness=0.5),
            Line(
              points={{10,0},{100,0}},
              color={240,0,0},
              visible=inclFacesX,
              thickness=0.5),
            Line(
              points={{0,-60},{0,-100}},
              color={253,52,56},
              visible=inclFacesY,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{0,100},{0,50}},
              color={253,52,56},
              visible=inclFacesY,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{-50,-50},{-10,-10}},
              color={253,52,56},
              visible=inclFacesZ,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{20,20},{50,50}},
              color={253,52,56},
              visible=inclFacesZ,
              smooth=Smooth.None,
              thickness=0.5),
            Text(
              extent={{-100,60},{100,100}},
              textString="%name",
              visible=not inclFacesY,
              lineColor={0,0,0})}));

    end AnGDL;

    model Sigracet10BA "<html>SGL Carbon Group Sigracet&reg; 10 BA</html>"
      extends AnGDL(
        L_x={0.400*U.mm},
        epsilon=0.88,
        Subregion(graphite('e-'(sigma=40*U.S/(12*U.cm)))));
      // Additional properties not incorporated [SGL2007]:
      //     Diffusivity:  L = 0.400 mm, p = 0.85 m/s (for air) => D = P*L = 340 mm2/s
      //     Density:  (85 g/m2)/(0.400 mm)/0.88 = 212.5 kg/m3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">SGL2007</a>].</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end Sigracet10BA;

    model Sigracet10BB "<html>SGL Carbon Group Sigracet&reg; 10 BB</html>"
      extends AnGDL(
        L_x={0.420*U.mm},
        epsilon=0.84,
        Subregion(graphite('e-'(sigma=42*U.S/(15*U.cm)))));
      // Additional properties not incorporated [SGL2007]:
      //     Diffusivity:  L = 0.420 mm, p = 0.03 m/s (for air) => D = P*L = 12.6 mm2/s
      //     Density:  (125 g/m2)/(0.420 mm) = 297.62 kg/m3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">SGL2007</a>].</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end Sigracet10BB;

    model Sigracet10BC "<html>SGL Carbon Group Sigracet&reg; 10 BC</html>"
      extends AnGDL(
        L_x={0.420*U.mm},
        epsilon=0.82,
        Subregion(graphite('e-'(sigma=42*U.S/(16*U.cm)))));
      // Additional properties not incorporated [SGL2007]:
      //     Diffusivity:  L = 0.420 mm, p = 0.0145 m/s (for air) => D = P*L = 6.09 mm2/s
      //     Density:  (135 g/m2)/(0.420 mm) = 321.43 kg/m3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">SGL2007</a>].</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end Sigracet10BC;

    model Sigracet24BA "<html>SGL Carbon Group Sigracet&reg; 24 BA</html>"
      extends AnGDL(
        L_x={0.190*U.mm},
        epsilon=0.84,
        Subregion(graphite('e-'(sigma=19*U.S/(10*U.cm)))));
      // Additional properties not incorporated [SGL2004]:
      //     Diffusivity:  L = 0.190 mm, p = 0.30 m/s (for air) => D = P*L = 57 mm2/s
      //     Density:  (54 g/m2)/(0.190 mm) = 284.21 kg/m3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">SGL2004</a>].</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end Sigracet24BA;

    model Sigracet24BC "<html>SGL Carbon Group Sigracet&reg; 24 BC</html>"
      extends AnGDL(
        L_x={0.235*U.mm},
        epsilon=0.76,
        Subregion(graphite('e-'(sigma=23.5*U.S/(11*U.cm)))));
      // Additional properties not incorporated [SGL2004]:
      //     Diffusivity:  L = 0.235 mm, p = 0.0045 m/s (for air) => D = P*L = 1.0575 mm2/s
      //     Density:  (100 g/m2)/(0.235 mm) = 425.53 kg/m3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">SGL2004</a>].</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end Sigracet24BC;

    model Sigracet25BA "<html>SGL Carbon Group Sigracet&reg; 25 BA</html>"
      extends AnGDL(
        L_x={0.190*U.mm},
        epsilon=0.88,
        Subregion(graphite('e-'(sigma=19*U.S/(10*U.cm)))));
      // Additional properties not incorporated [SGL2004]:
      //     Diffusivity:  L = 0.190 mm, p = 0.90 m/s (for air) => D = P*L = 171 mm2/s
      //     Density:  (40 g/m2)/(0.190 mm) = 210.53 kg/m3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">SGL2004</a>].</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end Sigracet25BA;

    model Sigracet25BC "<html>SGL Carbon Group Sigracet&reg; 25 BC</html>"
      extends AnGDL(
        L_x={0.235*U.mm},
        epsilon=0.80,
        Subregion(graphite('e-'(sigma=23.5*U.S/(12*U.cm)))));
      // Additional properties not incorporated [SGL2004]:
      //     Diffusivity:  L = 0.235 mm, p = 0.008 m/s (for air) => D = P*L = 1.88 mm2/s
      //     Density:  (86 g/m2)/(0.235 mm) = 365.96 kg/m3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">SGL2004</a>].</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end Sigracet25BC;

    model TorayTGPH030 "Toray Industries TGP-H-030"
      extends AnGDL(
        L_x={0.11*U.mm},
        epsilon=0.80,
        Subregion(graphite('C+'(theta=U.m*U.K/(1.7*U.W)),'e-'(sigma=U.S/(0.80*U.mm)))));
      // Additional properties not incorporated [Toray2010]:
      //     Diffusivity:  L = 0.110 mm, P/p = 2500 ml.mm/(cm2.hr.mmAq) = 0.70814e-3 m/(s.kPa)
      //         => D = P*L = 7.89e-6 m2/s, assuming p = 101.325 kPa
      //     Bulk density:  0.44 g/cm3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">Toray2010</a>].  The electrical conductivity (&sigma;) is based on the
  through-plane value of resistivity.  The thermal conductivity is not listed but is
  taken to be the same as for TGP-H-060, TGP-H-090, and TGP-H-120.</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end TorayTGPH030;

    model TorayTGPH060 "Toray Industries TGP-H-060"
      extends AnGDL(
        L_x={0.19*U.mm},
        epsilon=0.78,
        Subregion(graphite('C+'(theta=U.m*U.K/(1.7*U.W)),'e-'(sigma=U.S/(0.80*U.mm)))));

      // Additional properties not incorporated [Toray2010]:
      //     Diffusivity: L = 0.190 mm, P/p = 1900 ml.mm/(cm2.hr.mmAq) = 0.53818e-3 m/(s.kPa)
      //         => D = P*L = 10.36e-6 m2/s, assuming p = 101.325 kPa
      //     Bulk density:  0.44 g/cm3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">Toray2010</a>].  The electrical conductivity (&sigma;)  is based on the
  through-plane value of resistivity.  The thermal conductivity is through the plane at
  room temperature.</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end TorayTGPH060;

    model TorayTGPH090 "Toray Industries TGP-H-090"
      extends AnGDL(
        L_x={0.28*U.mm},
        epsilon=0.78,
        Subregion(graphite('C+'(theta=U.m*U.K/(1.7*U.W)),'e-'(sigma=U.S/(0.80*U.mm)))));

      // Additional properties not incorporated [Toray2010]:
      //     Diffusivity: L = 0.280 mm, P/p = 1700 ml.mm/(cm2.hr.mmAq) = 0.48153e-3 m/(s.kPa)
      //         => D = P*L = 13.66e-6 m2/s, assuming p = 101.325 kPa
      //     Bulk density:  0.44 g/cm3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">Toray2010</a>].  The electrical conductivity  (&sigma;) is based on the
  through-plane value of resistivity.  The thermal conductivity is through the plane at
  room temperature.</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end TorayTGPH090;

    model TorayTGPH120 "Toray Industries TGP-H-030"
      extends AnGDL(
        L_x={0.37*U.mm},
        epsilon=0.78,
        Subregion(graphite('C+'(theta=U.m*U.K/(1.7*U.W)),'e-'(sigma=U.S/(0.80*U.mm)))));

      // Additional properties not incorporated [Toray2010]:
      //     Diffusivity: L = 0.370 mm, P/p = 1500 ml.mm/(cm2.hr.mmAq) = 0.42488e-3 m/(s.kPa)
      //         => D = P*L = 15.93e-6 m2/s, assuming p = 101.325 kPa
      //     Bulk density:  0.44 g/cm3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">Toray2010</a>].  The electrical conductivity  (&sigma;) is based on the
  through-plane value of resistivity.  The thermal conductivity is through the plane at
  room temperature.</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end TorayTGPH120;

  end AnGDLs;

  package AnCLs "Anode catalyst layers"
    extends Modelica.Icons.Package;

    model AnCL "Anode catalyst layer"
      // extends FCSys.BaseClasses.Icons.Names.Top4;

      extends Region(
        L_x={28.7*U.micro*U.m},
        L_y={2*U.m},
        L_z={5*U.mm},
        final inclFacesX=true,
        inclFacesY=false,
        inclFacesZ=false,
        redeclare replaceable model Subregion = Subregions.Subregion (
            inclTransX=true,
            inclTransY=false,
            inclTransZ=false,
            gas(
              reduceVel=true,
              reduceTemp=true,
              inclH2=true,
              inclH2O=true,
              k_T=fill(epsilon^(3/2), 3),
              H2(
                initTransX=InitTranslational.None,
                initTransY=InitTranslational.None,
                initTransZ=InitTranslational.None,
                initEnergy=InitScalar.None,
                p_IC=environment.p - environment.p_H2O),
              H2O(consTransX=Conservation.IC, p_IC=environment.p_H2O)),
            graphite(
              reduceTemp=true,
              'inclC+'=true,
              'incle-'=true,
              k_T=fill((0.5*(1 - epsilon))^(3/2), 3),
              'C+'(V_IC=0.5*(V - epsilonV), theta=U.m*U.K/(1.18*U.W)),
              'e-'(
                V_IC=0.5*(V - epsilonV),
                N(stateSelect=StateSelect.prefer),
                initTransX=InitTranslational.None,
                initEnergy=InitScalar.None,
                phi(each stateSelect=StateSelect.default))),
            ionomer(
              reduceTemp=true,
              'inclC19HF37O5S-'=true,
              'inclH+'=true,
              inclH2O=false,
              k_T=fill((0.5*(1 - epsilon))^(3/2), 3),
              'C19HF37O5S-'(V_IC=0.5*(V - epsilonV)),
              'H+'(
                initMaterial=InitScalar.None,
                initTransX=InitTranslational.None,
                initEnergy=InitScalar.None),
              H2O(initEnergy=InitScalar.None, initMaterial=InitScalar.None)),
            liquid(
              inclH2O=false,
              k_T=fill(epsilon^(3/2), 3),
              H2O(V_IC=Modelica.Constants.eps*U.cc, consTransX=Conservation.IC)),

            reaction(J_0=0.5*U.A/U.cm^2))) annotation (IconMap(
            primitivesVisible=false));

      // **temp 0 reaction rate
      // **e-: rho_IC=17842.7*U.C/U.cc,
      // initMaterial=InitScalar.None,
      // chemical(Ndot(start=0, fixed=true))
      // **e-: sigma=40*U.S/(12*U.cm),
      // **H+: mu=0.083*U.S/(0.95*U.M*U.cm),

      // TODO:  Initialize for zero reaction rate
      // (initMaterial=InitScalar.ReactionRate?) for this and CaCL.

      // See the documentation layer of Subregions.Phases.BaseClasses.EmptyPhase
      // regarding the settings of k_T for each phase.

      parameter Q.NumberAbsolute epsilon(nominal=1) = 0.25
        "<html>Volumetric porosity (&epsilon;)</html>"
        annotation (group="Geometry");
      parameter Q.NumberAbsolute lambda_IC=14
        "<html>Initial molar ratio of H<sub>2</sub>O to SO<sub>3</sub><sup>-</sup> (&lambda;<sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization"));

      // Auxiliary variables (for analysis)
      output Q.NumberAbsolute RH[n_x, n_y, n_z](
        each stateSelect=StateSelect.never,
        each displayUnit="%") =
        Modelica.Media.Air.MoistAir.saturationPressureLiquid(subregions.gas.H2O.T
        /U.K) ./ subregions.gas.dalton.p if environment.analysis and
        hasSubregions "Relative humidity";

    protected
      final parameter Q.Volume epsilonV=epsilon*V "Fluid volume";

      outer Conditions.Environment environment "Environmental conditions";

    initial equation
      //  subregions.ionomer.H2O.N = lambda_IC*subregions.ionomer.'C19HF37O5S-'.N;
      annotation (
        defaultComponentPrefixes="replaceable",
        Documentation(info="<html>
<p>This model represents the anode catalyst layer of a PEMFC.
The x axis extends from the anode to the cathode.
The y axis extends along the length of the channel and
the z axis extends across the width of the channel.</p>

<p>By default, the cross-sectional area in the xy plane is 100 cm<sup>2</sup>.</p>

<p>The default thickness (<code>L_x = {28.7*U.micro*U.m}</code>) is from [<a href=\"modelica://FCSys.UsersGuide.References\">Gurau1998</a>].
The default thermal conductivity of the carbon (<code>theta = U.m*U.K/(1.18*U.W)</code>)
   represents a compressed SGL Sigracet 10 BA gas diffusion layer
  [<a href=\"modelica://FCSys.UsersGuide.References\">Nitta2008</a>].  The default
  electronical conductivity (<code>sigma=40*U.S/(12*U.cm)</code>)
  is for SGL Carbon Group Sigracet&reg; 10 BA (see <a href=\"modelica://FCSys.Regions.AnGDLs.Sigracet10BA\">AnGDLs.Sigracet10BA</a>).
  The default
  protonic mobility (<code>mu=0.083*U.S/(0.95*U.M*U.cm)</code>)
  is for DuPont<sup>TM</sup> Nafion&reg; N-112 (see <a href=\"modelica://FCSys.Regions.PEMs.DuPontN112\">PEMs.DuPontN112</a>).</p>

<p>Default assumptions (may be adjusted):
<ol>
<li>Half of the solid is graphite and half is ionomer (by volume).</li>
<li>There is no liquid water.</li>
</ol></p>

<p>For more information, see the
    <a href=\"modelica://FCSys.Regions.Region\">Region</a> model.</p>
</html>"),
        Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.1), graphics={
            Rectangle(
              extent={{-100,60},{100,100}},
              fillColor={255,255,255},
              visible=not inclFacesY,
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Rectangle(
              extent={{-22.6085,-62.7355},{-56.551,-79.7156}},
              lineColor={64,64,64},
              rotation=45,
              fillColor={127,127,127},
              fillPattern=FillPattern.HorizontalCylinder,
              origin={-34.3741,132.347}),
            Rectangle(
              extent={{-105.39,79.1846},{-139.33,70.6991}},
              lineColor={0,0,0},
              fillColor={200,200,200},
              rotation=45,
              fillPattern=FillPattern.HorizontalCylinder,
              origin={148.513,82.5291}),
            Polygon(
              points={{-14,40},{6,60},{14,60},{-6,40},{-14,40}},
              fillPattern=FillPattern.HorizontalCylinder,
              smooth=Smooth.None,
              fillColor={0,0,0},
              pattern=LinePattern.None),
            Rectangle(
              extent={{-6,40},{6,-60}},
              lineColor={0,0,0},
              fillColor={200,200,200},
              fillPattern=FillPattern.VerticalCylinder),
            Rectangle(
              extent={{-38,40},{-14,-60}},
              lineColor={64,64,64},
              fillColor={127,127,127},
              fillPattern=FillPattern.VerticalCylinder),
            Rectangle(
              extent={{-14,40},{-6,-60}},
              fillPattern=FillPattern.Solid,
              fillColor={0,0,0},
              pattern=LinePattern.None),
            Polygon(
              points={{-20,0},{-20,40},{0,60},{20,60},{20,0},{42,0},{42,80},{-42,
                  80},{-42,0},{-20,0}},
              fillPattern=FillPattern.HorizontalCylinder,
              smooth=Smooth.None,
              fillColor={255,255,255},
              pattern=LinePattern.None),
            Polygon(
              points={{-20,0},{-20,-60},{0,-60},{20,-40},{20,0},{42,0},{42,-80},
                  {-42,-80},{-42,0},{-20,0}},
              fillPattern=FillPattern.HorizontalCylinder,
              smooth=Smooth.None,
              fillColor={255,255,255},
              pattern=LinePattern.None),
            Polygon(points={{0,60},{20,60},{0,40},{-20,40},{0,60}}, lineColor={
                  0,0,0}),
            Rectangle(
              extent={{-20,40},{0,-60}},
              pattern=LinePattern.None,
              lineColor={0,0,0}),
            Polygon(
              points={{20,60},{0,40},{0,-60},{20,-40},{20,60}},
              lineColor={0,0,0},
              fillColor={200,200,200},
              fillPattern=FillPattern.Solid),
            Line(
              points={{-20,0},{-100,0}},
              color={240,0,0},
              visible=inclFacesX,
              thickness=0.5),
            Line(
              points={{10,0},{100,0}},
              color={240,0,0},
              visible=inclFacesX,
              thickness=0.5),
            Line(
              points={{0,-60},{0,-100}},
              color={253,52,56},
              visible=inclFacesY,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{0,100},{0,50}},
              color={253,52,56},
              visible=inclFacesY,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{-50,-50},{-10,-10}},
              color={253,52,56},
              visible=inclFacesZ,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{20,20},{50,50}},
              color={253,52,56},
              visible=inclFacesZ,
              smooth=Smooth.None,
              thickness=0.5),
            Text(
              extent={{-100,60},{100,100}},
              textString="%name",
              visible=not inclFacesY,
              lineColor={0,0,0})}));
    end AnCL;

    model AnCGDL "Integrated anode catalyst/gas diffusion layer"
      extends AnCLs.AnCL(L_x={(28.7*U.micro*U.m + 0.3*U.mm)});
      annotation (Documentation(info="<html><p>The default thickness is the total thickness of
  <a href=\"modelica://FCSys.Regions.AnCLs.AnCL\">AnCL</a> and
  <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a>.</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.AnCLs.AnCL\">AnCL</a> model.</p></html>"));

    end AnCGDL;
  end AnCLs;

  package PEMs "Proton exchange membranes"
    extends Modelica.Icons.Package;

    model PEM "Proton exchange membrane"
      // extends FCSys.BaseClasses.Icons.Names.Top4;

      // Note:  Extensions of PEM should be placed directly in the PEMs
      // package rather than subpackages (e.g., by manufacturer) so that
      // __Dymola_choicesFromPackage can be used.  Dymola 7.4 launches the
      // parameter dialogs too slowly when __Dymola_choicesAllMatching is
      // used.

      extends Region(
        L_x={100*U.micro*U.m},
        L_y={2*U.m},
        L_z={5*U.mm},
        final inclFacesX=true,
        inclFacesY=false,
        inclFacesZ=false,
        redeclare replaceable model Subregion =
            FCSys.Subregions.SubregionIonomerOnly (
            inclTransX=true,
            inclTransY=false,
            inclTransZ=false,
            void=true,
            ionomer(
              reduceTemp=true,
              'inclC19HF37O5S-'=true,
              'inclH+'=true,
              inclH2O=true,
              'C19HF37O5S-'(initMaterial=InitScalar.Pressure),
              'H+'(mu=0.083*U.S/(0.95*U.M*U.cm), initEnergy=InitScalar.None),
              H2O(initEnergy=InitScalar.None)))) annotation (IconMap(
            primitivesVisible=false));
      //initMaterial=InitScalar.None,

      parameter Q.NumberAbsolute lambda_IC=14
        "<html>Initial molar ratio of H<sub>2</sub>O to SO<sub>3</sub><sup>-</sup> (&lambda;<sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization"));

      // TODO: Map electrical resistivity and EOD ratio to mobility of ionomer and protons,
      // given the mobility of H2O.

      // Auxiliary variables (for analysis)
      output Q.Number lambda[n_x, n_y, n_z]=subregions.ionomer.H2O.N ./
          subregions.ionomer.'C19HF37O5S-'.N if environment.analysis
        "Molar ratio of H2O to SO3-";

    protected
      outer Conditions.Environment environment "Environmental conditions";

    initial equation
      //**subregions.ionomer.H2O.N = lambda_IC*subregions.ionomer.'C19HF37O5S-'.N;
      annotation (
        defaultComponentName="PEM",
        Documentation(info="<html>
<p>This model represents the proton exchange membrane of a PEMFC.
The x axis extends from the anode to the cathode.
The y axis extends along the length of the channel and
the z axis extends across the width of the channel.</p>

<p>By default, the cross-sectional area in the xy plane is 100 cm<sup>2</sup>.</p>

<p>The default
  protonic mobility (<code>mu=0.083*U.S/(0.95*U.M*U.cm)</code>)
  is for DuPont<sup>TM</sup> Nafion&reg; N-112 (see <a href=\"modelica://FCSys.Regions.PEMs.DuPontN112\">DuPontN112</a>).</p>

<p>Assumptions:<ol>
<li>There are no pores in the PEM.  All H<sub>2</sub>O is absorbed into the ionomer itself.</li>
</ul></p>

<p>For more information, see the
    <a href=\"modelica://FCSys.Regions.Region\">Region</a> model.</p>
</html>"),
        Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.1), graphics={
            Rectangle(
              extent={{-100,60},{100,100}},
              fillColor={255,255,255},
              visible=not inclFacesY,
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Rectangle(
              extent={{-99.092,-21.1179},{-84.9489,-63.5448}},
              lineColor={200,200,200},
              fillColor={255,255,255},
              rotation=-45,
              fillPattern=FillPattern.VerticalCylinder,
              origin={95.001,14.864}),
            Rectangle(
              extent={{-20,40},{0,-60}},
              lineColor={200,200,200},
              fillColor={255,255,255},
              fillPattern=FillPattern.VerticalCylinder),
            Polygon(
              points={{20,0},{42,0},{42,80},{-42,80},{-42,0},{-20,0},{-20,40},{
                  0,60},{20,60},{20,0}},
              smooth=Smooth.None,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Polygon(
              points={{20,0},{42,0},{42,-80},{-42,-80},{-42,0},{-20,0},{-20,-60},
                  {0,-60},{20,-40},{20,0}},
              smooth=Smooth.None,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Polygon(
              points={{0,40},{20,60},{20,-40},{0,-60},{0,40}},
              lineColor={0,0,0},
              smooth=Smooth.None,
              fillPattern=FillPattern.Solid,
              fillColor={200,200,200}),
            Rectangle(extent={{-20,40},{0,-60}}, lineColor={0,0,0}),
            Polygon(
              points={{0,60},{20,60},{0,40},{-20,40},{0,60}},
              lineColor={0,0,0},
              smooth=Smooth.None),
            Line(
              points={{-20,0},{-100,0}},
              color={240,0,0},
              visible=inclFacesX,
              thickness=0.5),
            Line(
              points={{10,0},{100,0}},
              color={0,0,240},
              visible=inclFacesX,
              thickness=0.5),
            Line(
              points={{0,-60},{0,-100}},
              color={127,127,127},
              visible=inclFacesY,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{0,100},{0,50}},
              color={127,127,127},
              visible=inclFacesY,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{-50,-50},{-10,-10}},
              color={127,127,127},
              visible=inclFacesZ,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{20,20},{50,50}},
              color={127,127,127},
              visible=inclFacesZ,
              smooth=Smooth.None,
              thickness=0.5),
            Text(
              extent={{-100,60},{100,100}},
              textString="%name",
              visible=not inclFacesY,
              lineColor={0,0,0})}));
    end PEM;

    model DuPontN112 "<html>DuPont<sup>TM</sup> Nafion&reg; N-112</html>"
      extends PEM(L_x={51*U.micro*U.m},Subregion(ionomer('H+'(mu=0.083*U.S/(
                  0.95*U.M*U.cm)))));
      // Additional properties not incorporated [DuPont2004N]:
      //     Density:  (100 g/m2)/(51 um) = 1.9608 g/cm3
      annotation (defaultComponentName="PEM", Documentation(info="<html>
  <p>The default value of protonic mobility (<code>mu=0.083*U.S/(0.95*U.M*U.cm)</code>)
  is based on the protonic conductivity  of
  DuPont<sup>TM</sup> Nafion&reg; N-112 (0.083 S/cm)
  [<a href=\"modelica://FCSys.UsersGuide.References\">DuPont2004N</a>]
  and the density of protons set by
  <a href=\"modelica://FCSys.Characteristics.'H+'.Ionomer\">FCSys.Characteristics.'H+'.Ionomer</a>
  (0.95 M).</p>

  <p>For more information, see the
          <a href=\"modelica://FCSys.Regions.PEMs.PEM\">PEM</a> model.</p></html>"));

    end DuPontN112;

    model DuPontN115 "<html>DuPont<sup>TM</sup> Nafion&reg; N-115</html>"
      extends PEM(L_x={127*U.micro*U.m},Subregion(ionomer('H+'(mu=0.083*U.S/(
                  0.95*U.M*U.cm)))));
      // Additional properties not incorporated [DuPont2004N] and [DuPont2005]:
      //     Density:  (250 g/m2)/(127 um) = 1.9685 g/cm3
      annotation (defaultComponentName="PEM", Documentation(info="<html>
  <p>The default value of protonic mobility (<code>mu=0.083*U.S/(0.95*U.M*U.cm)</code>)
  is based on the protonic conductivity  of
  DuPont<sup>TM</sup> Nafion&reg; N-115 (0.083 S/cm)
  [<a href=\"modelica://FCSys.UsersGuide.References\">DuPont2004N</a>]
  and the density of protons set by
  <a href=\"modelica://FCSys.Characteristics.'H+'.Ionomer\">FCSys.Characteristics.'H+'.Ionomer</a>
  (0.95 M).</p>

  <p>For more information, see the
          <a href=\"modelica://FCSys.Regions.PEMs.PEM\">PEM</a> model.</p></html>"));

    end DuPontN115;

    model DuPontN117 "<html>DuPont<sup>TM</sup> Nafion&reg; N-117</html>"
      extends PEM(L_x={183*U.micro*U.m},Subregion(ionomer('H+'(mu=0.083*U.S/(
                  0.95*U.M*U.cm)))));
      // Additional properties not incorporated [DuPont2004N] and [DuPont2005]:
      //     Density:  (360 g/m2)/(183 um) = 1.9672 g/cm3
      annotation (defaultComponentName="PEM", Documentation(info="<html>
  <p>The default value of protonic mobility (<code>mu=0.083*U.S/(0.95*U.M*U.cm)</code>)
  is based on the protonic conductivity  of
  DuPont<sup>TM</sup> Nafion&reg; N-117 (0.083 S/cm)
  [<a href=\"modelica://FCSys.UsersGuide.References\">DuPont2004N</a>]
  and the density of protons set by
  <a href=\"modelica://FCSys.Characteristics.'H+'.Ionomer\">FCSys.Characteristics.'H+'.Ionomer</a>
  (0.95 M).</p>

  <p>For more information, see the
          <a href=\"modelica://FCSys.Regions.PEMs.PEM\">PEM</a> model.</p></html>"));

    end DuPontN117;

    model DuPontNE1110 "<html>DuPont<sup>TM</sup> Nafion&reg; NE-1110</html>"
      extends PEM(L_x={254*U.micro*U.m},Subregion(ionomer('H+'(mu=0.083*U.S/(
                  0.95*U.M*U.cm)))));
      // Additional properties not incorporated [DuPont2004N] and [DuPont2005]:
      //     Density:  (500 g/m2)/(254 um) = 1.9685 g/cm3
      annotation (defaultComponentName="PEM", Documentation(info="<html>
  <p>The default value of protonic mobility (<code>mu=0.083*U.S/(0.95*U.M*U.cm)</code>)
  is based on the protonic conductivity  of
  DuPont<sup>TM</sup> Nafion&reg; NE-1110 (0.083 S/cm)
  [<a href=\"modelica://FCSys.UsersGuide.References\">DuPont2004N</a>]
  and the density of protons set by
  <a href=\"modelica://FCSys.Characteristics.'H+'.Ionomer\">FCSys.Characteristics.'H+'.Ionomer</a>
  (0.95 M).</p>

  <p>For more information, see the
          <a href=\"modelica://FCSys.Regions.PEMs.PEM\">PEM</a> model.</p></html>"));

    end DuPontNE1110;

    model DuPontNE1135 "<html>DuPont<sup>TM</sup> Nafion&reg; NE-1135</html>"
      extends PEM(L_x={89*U.micro*U.m},Subregion(ionomer('H+'(mu=0.083*U.S/(
                  0.95*U.M*U.cm)))));
      // Additional properties not incorporated [DuPont2004N] and [DuPont2005]:
      //     Density:  (190 g/m2)/(89 um) = 2.1348 g/cm3
      annotation (defaultComponentName="PEM", Documentation(info="<html>
  <p>The default value of protonic mobility (<code>mu=0.083*U.S/(0.95*U.M*U.cm)</code>)
  is based on the protonic conductivity  of
  DuPont<sup>TM</sup> Nafion&reg; NE-1135 (0.083 S/cm)
  [<a href=\"modelica://FCSys.UsersGuide.References\">DuPont2004N</a>]
  and the density of protons set by
  <a href=\"modelica://FCSys.Characteristics.'H+'.Ionomer\">FCSys.Characteristics.'H+'.Ionomer</a>
  (0.95 M).</p>

  <p>For more information, see the
          <a href=\"modelica://FCSys.Regions.PEMs.PEM\">PEM</a> model.</p></html>"));

    end DuPontNE1135;

    model DuPontNRE211 "<html>DuPont<sup>TM</sup> Nafion&reg; NRE-1110</html>"
      extends PEM(L_x={25.4*U.micro*U.m},Subregion(ionomer('H+'(mu=0.083*U.S/(
                  0.95*U.M*U.cm)))));
      // Additional properties not incorporated [DuPont2004NRE]:
      //     Density:  (50 g/m2)/(25.4 um) = 1.9685 g/cm3
      annotation (defaultComponentName="PEM", Documentation(info="<html>
  <p>The default value of protonic mobility (<code>mu=0.083*U.S/(0.95*U.M*U.cm)</code>)
  is based on the protonic conductivity listed for the
  DuPont<sup>TM</sup> Nafion&reg; N and NE series (0.083 S/cm)
  [<a href=\"modelica://FCSys.UsersGuide.References\">DuPont2004N</a>]
  and the density of protons set by
  <a href=\"modelica://FCSys.Characteristics.'H+'.Ionomer\">FCSys.Characteristics.'H+'.Ionomer</a>
  (0.95 M).  The conductivity is not listed for DuPont<sup>TM</sup> Nafion&reg; NRE-211
  in [<a href=\"modelica://FCSys.UsersGuide.References\">DuPont2004NRE</a>].</p>

  <p>For more information, see the
          <a href=\"modelica://FCSys.Regions.PEMs.PEM\">PEM</a> model.</p></html>"));

    end DuPontNRE211;

    model DuPontNRE212 "<html>DuPont<sup>TM</sup> Nafion&reg; NRE-1110</html>"
      extends PEM(L_x={50.8*U.micro*U.m},Subregion(ionomer('H+'(mu=0.083*U.S/(
                  0.95*U.M*U.cm)))));
      // Additional properties not incorporated [DuPont2004NRE]:
      //     Density:  (100 g/m2)/(50.8 um) = 1.9685 g/cm3
      annotation (defaultComponentName="PEM", Documentation(info="<html>
  <p>The default value of protonic mobility (<code>mu=0.083*U.S/(0.95*U.M*U.cm)</code>)
  is based on the protonic conductivity listed for the
  DuPont<sup>TM</sup> Nafion&reg; N and NE series (0.083 S/cm)
  [<a href=\"modelica://FCSys.UsersGuide.References\">DuPont2004N</a>]
  and the density of protons set by
  <a href=\"modelica://FCSys.Characteristics.'H+'.Ionomer\">FCSys.Characteristics.'H+'.Ionomer</a>
  (0.95 M).  The conductivity is not listed for DuPont<sup>TM</sup> Nafion&reg; NRE-212
  in [<a href=\"modelica://FCSys.UsersGuide.References\">DuPont2004NRE</a>].</p>

  <p>For more information, see the
  <a href=\"modelica://FCSys.Regions.PEMs.PEM\">PEM</a> model.</p></html>"));

    end DuPontNRE212;

  end PEMs;

  package CaCLs "Cathode catalyst layers"
    extends Modelica.Icons.Package;

    model CaCL "Cathode catalyst layer"
      // extends FCSys.BaseClasses.Icons.Names.Top4;

      extends Region(
        L_x={28.7*U.micro*U.m},
        L_y={2*U.m},
        L_z={5*U.mm},
        final inclFacesX=true,
        inclFacesY=false,
        inclFacesZ=false,
        redeclare replaceable model Subregion = Subregions.Subregion (
            inclTransX=true,
            inclTransY=false,
            inclTransZ=false,
            gas(
              reduceVel=true,
              reduceTemp=true,
              k_T=fill(epsilon^(3/2), 3),
              inclH2O=true,
              inclN2=true,
              inclO2=true,
              H2O(p_IC=environment.p_H2O, consTransX=Conservation.IC),
              N2(
                initTransX=InitTranslational.None,
                initTransY=InitTranslational.None,
                initTransZ=InitTranslational.None,
                initEnergy=InitScalar.None,
                p_IC=(1 - environment.n_O2)*(environment.p - environment.p_H2O)),

              O2(
                initTransX=InitTranslational.None,
                initTransY=InitTranslational.None,
                initTransZ=InitTranslational.None,
                initEnergy=InitScalar.None,
                p_IC=environment.n_O2*(environment.p - environment.p_H2O))),
            graphite(
              reduceTemp=true,
              k_T=fill((0.5*(1 - epsilon))^(3/2), 3),
              'inclC+'=true,
              'incle-'=true,
              'C+'(V_IC=0.5*(V - epsilonV),theta=U.m*U.K/(1.18*U.W)),
              'e-'(
                initEnergy=InitScalar.None,
                V_IC=0.5*(V - epsilonV),
                initTransX=InitTranslational.None,
                sigma=U.S/(0.012*U.cm))),
            ionomer(
              reduceTemp=true,
              k_T=fill((0.5*(1 - epsilon))^(3/2), 3),
              'inclC19HF37O5S-'=true,
              inclH2O=true,
              'inclH+'=true,
              'C19HF37O5S-'(V_IC=0.5*(V - epsilonV)),
              'H+'(initEnergy=InitScalar.None, mu=0.083*U.S/(0.95*U.M*U.cm)),
              H2O(initEnergy=InitScalar.None, initMaterial=InitScalar.None)),
            liquid(k_T=fill(epsilon^(3/2), 3)))) annotation (IconMap(
            primitivesVisible=false));

      // See the documentation layer of Subregions.Phases.BaseClasses.EmptyPhase
      // regarding the settings of k_T for each phase.

      parameter Q.NumberAbsolute epsilon(nominal=1) = 0.25
        "<html>Volumetric porosity (&epsilon;)</html>"
        annotation (group="Geometry");
      parameter Q.NumberAbsolute lambda_IC=14
        "<html>Initial molar ratio of H<sub>2</sub>O to SO<sub>3</sub><sup>-</sup> (&lambda;<sub>IC</sub>)</html>"
        annotation (Dialog(tab="Initialization"));

      // Auxiliary variables (for analysis)
      output Q.NumberAbsolute RH[n_x, n_y, n_z](
        each stateSelect=StateSelect.never,
        each displayUnit="%") =
        Modelica.Media.Air.MoistAir.saturationPressureLiquid(subregions.gas.H2O.T
        /U.K) ./ subregions.gas.dalton.p if environment.analysis and
        hasSubregions "Relative humidity";
      output Q.Number n_O2[n_x, n_y, n_z](each stateSelect=StateSelect.never)
         = subregions.gas.O2.N ./ (subregions.gas.N2.N + subregions.gas.O2.N)
        if environment.analysis and hasSubregions "Dry-gas concentration of O2";

    protected
      final parameter Q.Volume epsilonV=epsilon*V "Fluid volume";

      outer Conditions.Environment environment "Environmental conditions";

    initial equation
      subregions.ionomer.H2O.N = lambda_IC*subregions.ionomer.'C19HF37O5S-'.N;
      annotation (Documentation(info="<html>
<p>This model represents the cathode catalyst layer of a PEMFC.
The x axis extends from the anode to the cathode.
The y axis extends along the length of the channel and
the z axis extends across the width of the channel.</p>

<p>By default, the cross-sectional area in the xy plane is 100 cm<sup>2</sup>.</p>

<p>The default thickness (<code>L_x = {28.7*U.micro*U.m}</code>) is from [<a href=\"modelica://FCSys.UsersGuide.References\">Gurau1998</a>].
The default thermal conductivity of the carbon (<code>theta = U.m*U.K/(1.18*U.W)</code>)
   represents a compressed SGL Sigracet 10 BA gas diffusion layer
  [<a href=\"modelica://FCSys.UsersGuide.References\">Nitta2008</a>].  The default
  electronic mobility (<code>sigma=40*U.S/(12*U.cm)</code>)
  is for SGL Carbon Group Sigracet&reg; 10 BA (see <a href=\"modelica://FCSys.Regions.CaGDLs.Sigracet10BA\">CAGDLs.Sigracet10BA</a>).
  The default
  protonic mobility (<code>mu=0.083*U.S/(0.95*U.M*U.cm)</code>)
  is for DuPont<sup>TM</sup> Nafion&reg; N-112 (see <a href=\"modelica://FCSys.Regions.PEMs.DuPontN112\">PEMs.DuPontN112</a>).</p>

<p>Default assumptions (may be adjusted):
<ol>
<li>Half of the solid is graphite and half is ionomer (by volume).</li>
<li>There is no liquid water.</li>
</ol></p>

<p>For more information, see the
    <a href=\"modelica://FCSys.Regions.Region\">Region</a> model.</p>
</html>"), Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.1), graphics={
            Rectangle(
              extent={{-100,60},{100,100}},
              fillColor={255,255,255},
              visible=not inclFacesY,
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Rectangle(
              extent={{-21.6329,-68.4511},{-58.4038,-85.4311}},
              lineColor={64,64,64},
              rotation=45,
              fillColor={127,127,127},
              fillPattern=FillPattern.HorizontalCylinder,
              origin={-15.1055,127.699}),
            Rectangle(
              extent={{-105.385,79.1805},{-139.323,70.6948}},
              lineColor={0,0,0},
              fillColor={200,200,200},
              rotation=45,
              fillPattern=FillPattern.HorizontalCylinder,
              origin={130.507,84.5292}),
            Polygon(
              points={{-14,40},{6,60},{14,60},{-6,40},{-14,40}},
              fillPattern=FillPattern.HorizontalCylinder,
              smooth=Smooth.None,
              fillColor={0,0,0},
              pattern=LinePattern.None),
            Rectangle(
              extent={{-26,40},{-14,-60}},
              lineColor={0,0,0},
              fillColor={200,200,200},
              fillPattern=FillPattern.VerticalCylinder),
            Rectangle(
              extent={{-6,40},{18,-60}},
              lineColor={64,64,64},
              fillColor={127,127,127},
              fillPattern=FillPattern.VerticalCylinder),
            Rectangle(
              extent={{-14,40},{-6,-60}},
              fillPattern=FillPattern.Solid,
              fillColor={0,0,0},
              pattern=LinePattern.None),
            Polygon(
              points={{-20,0},{-20,40},{0,60},{20,60},{20,0},{42,0},{42,80},{-42,
                  80},{-42,0},{-20,0}},
              fillPattern=FillPattern.HorizontalCylinder,
              smooth=Smooth.None,
              fillColor={255,255,255},
              pattern=LinePattern.None),
            Polygon(
              points={{-20,0},{-20,-60},{0,-60},{20,-40},{20,0},{42,0},{42,-80},
                  {-42,-80},{-42,0},{-20,0}},
              fillPattern=FillPattern.HorizontalCylinder,
              smooth=Smooth.None,
              fillColor={255,255,255},
              pattern=LinePattern.None),
            Polygon(points={{0,60},{20,60},{0,40},{-20,40},{0,60}}, lineColor={
                  0,0,0}),
            Rectangle(
              extent={{-20,40},{0,-60}},
              pattern=LinePattern.None,
              lineColor={0,0,0}),
            Polygon(
              points={{20,60},{0,40},{0,-60},{20,-40},{20,60}},
              lineColor={0,0,0},
              fillColor={120,120,120},
              fillPattern=FillPattern.Solid),
            Line(
              points={{-20,0},{-100,0}},
              color={0,0,240},
              visible=inclFacesX,
              thickness=0.5),
            Line(
              points={{10,0},{100,0}},
              color={0,0,240},
              visible=inclFacesX,
              thickness=0.5),
            Line(
              points={{0,-60},{0,-98}},
              color={0,0,240},
              visible=inclFacesY,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{0,100},{0,50}},
              color={0,0,240},
              visible=inclFacesY,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{-50,-50},{-10,-10}},
              color={0,0,240},
              visible=inclFacesZ,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{20,20},{50,50}},
              color={0,0,240},
              visible=inclFacesZ,
              smooth=Smooth.None,
              thickness=0.5),
            Text(
              extent={{-100,60},{100,100}},
              textString="%name",
              visible=not inclFacesY,
              lineColor={0,0,0})}));
    end CaCL;

    model CaCGDL "Integrated cathode catalyst/gas diffusion layer"
      extends CaCLs.CaCL(L_x={(28.7*U.micro*U.m + 0.3*U.mm)});
      annotation (Documentation(info="<html><p>The default thickness is the total thickness of
  <a href=\"modelica://FCSys.Regions.CaCLs.CaCL\">CaCL</a> and
  <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a>.</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.CaCLs.CaCL\">CaCL</a> model.</p></html>"));

    end CaCGDL;
  end CaCLs;

  package CaGDLs "Cathode gas diffusion layers"
    extends Modelica.Icons.Package;

    // Note:  Extensions of CaGDL should be placed directly in the CaGDLs
    // package rather than nested packages (e.g., by manufacturer) so that
    // __Dymola_choicesFromPackage can be used.  In Dymola 7.4 the
    // parameter dialogs launch too slowly when __Dymola_choicesAllMatching
    // is used.

    model CaGDL "Cathode gas diffusion layer"
      // extends FCSys.BaseClasses.Icons.Names.Top4;

      // Note:  Extensions of CaGDL should be placed directly in the CaGDLs
      // package rather than subpackages (e.g., by manufacturer) so that
      // __Dymola_choicesFromPackage can be used.  Dymola 7.4 launches the
      // parameter dialogs too slowly when __Dymola_choicesAllMatching is
      // used.

      extends Region(
        L_x={0.3*U.mm},
        L_y={2*U.m},
        L_z={5*U.mm},
        final inclFacesX=true,
        inclFacesY=false,
        inclFacesZ=false,
        redeclare replaceable model Subregion =
            FCSys.Subregions.SubregionNoIonomer (
            inclTransX=true,
            inclTransY=false,
            inclTransZ=false,
            gas(
              reduceVel=true,
              reduceTemp=true,
              k_T=fill(epsilon^(3/2), 3),
              inclH2O=true,
              inclN2=true,
              inclO2=true,
              H2O(p_IC=environment.p_H2O, consTransX=Conservation.IC),
              N2(
                initTransX=InitTranslational.None,
                initTransY=InitTranslational.None,
                initTransZ=InitTranslational.None,
                initEnergy=InitScalar.None,
                p_IC=(1 - environment.n_O2)*(environment.p - environment.p_H2O)),

              O2(
                initTransX=InitTranslational.None,
                initTransY=InitTranslational.None,
                initTransZ=InitTranslational.None,
                initEnergy=InitScalar.None,
                p_IC=environment.n_O2*(environment.p - environment.p_H2O))),
            graphite(
              reduceTemp=true,
              k_T=fill((1 - epsilon)^(3/2), 3),
              'inclC+'=true,
              'incle-'=true,
              'C+'(V_IC=V - epsilonV,theta=U.m*U.K/(1.18*U.W)),
              'e-'(
                rho_IC=17842.7*U.C/U.cc,
                sigma=40*U.S/(12*U.cm),
                initTransX=InitTranslational.None,
                initEnergy=InitScalar.None)),
            liquid(
              k_T=fill(epsilon^(3/2), 3),
              inclH2O=true,
              H2O(V_IC=Modelica.Constants.eps*U.cc, consTransX=Conservation.IC))))
        annotation (IconMap(primitivesVisible=false));

      // See the documentation layer of Subregions.Phases.BaseClasses.EmptyPhase
      // regarding the settings of k_T for each phase.

      parameter Q.NumberAbsolute epsilon(nominal=1) = 0.4
        "<html>Volumetric porosity (&epsilon;)</html>"
        annotation (group="Geometry");

      // Auxiliary variables (for analysis)
      output Q.NumberAbsolute RH[n_x, n_y, n_z](
        each stateSelect=StateSelect.never,
        each displayUnit="%") =
        Modelica.Media.Air.MoistAir.saturationPressureLiquid(subregions.gas.H2O.T
        /U.K) ./ subregions.gas.dalton.p if environment.analysis and
        hasSubregions "Relative humidity";
      output Q.Number n_O2[n_x, n_y, n_z](each stateSelect=StateSelect.never)
         = subregions.gas.O2.N ./ (subregions.gas.N2.N + subregions.gas.O2.N)
        if environment.analysis and hasSubregions "Dry-gas concentration of O2";

    protected
      final parameter Q.Volume epsilonV=epsilon*V "Fluid volume";

      outer Conditions.Environment environment "Environmental conditions";
      annotation (
        Documentation(info="<html>
<p>This model represents the cathode gas diffusion layer of a PEMFC.
The x axis extends from the anode to the cathode.
The y axis extends along the length of the channel and
the z axis extends across the width of the channel.</p>

<p>By default, the cross-sectional area in the xy plane is 100 cm<sup>2</sup>.</p>

<p>The default porosity (<code>epsilon = 0.4</code>) is
of a compressed GDL according to [<a href=\"modelica://FCSys.UsersGuide.References\">Bernardi1992</a>, p. 2483, Table 3].
  The default thermal conductivity of the carbon (<code>theta = U.m*U.K/(1.18*U.W)</code>)
   represents a compressed SGL Sigracet 10 BA gas diffusion layer
  [<a href=\"modelica://FCSys.UsersGuide.References\">Nitta2008</a>].  The default
  electrical conductivity is also
  for SGL Carbon Group Sigracet&reg; 10 BA [<a href=\"modelica://FCSys.UsersGuide.References\">SGL2007</a>].</p>

<p>For more information, see the
    <a href=\"modelica://FCSys.Regions.Region\">Region</a> model.</p>
</html>"),
        Diagram(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.1)),
        Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.1), graphics={
            Rectangle(
              extent={{-100,60},{100,100}},
              fillColor={255,255,255},
              visible=not inclFacesY,
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Rectangle(
              extent={{-78.7855,18.6813},{-50.5004,-23.7455}},
              lineColor={64,64,64},
              fillColor={127,127,127},
              rotation=-45,
              fillPattern=FillPattern.VerticalCylinder,
              origin={52.5001,1.0805}),
            Rectangle(
              extent={{-20,40},{20,-60}},
              lineColor={64,64,64},
              fillColor={127,127,127},
              fillPattern=FillPattern.VerticalCylinder),
            Polygon(
              points={{20,0},{42,0},{42,80},{-42,80},{-42,0},{-20,0},{-20,40},{
                  0,60},{20,60},{20,0}},
              smooth=Smooth.None,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Polygon(
              points={{20,0},{42,0},{42,-80},{-42,-80},{-42,0},{-20,0},{-20,-60},
                  {0,-60},{20,-40},{20,0}},
              smooth=Smooth.None,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Polygon(
              points={{0,40},{20,60},{20,-40},{0,-60},{0,40}},
              lineColor={0,0,0},
              smooth=Smooth.None,
              fillPattern=FillPattern.Solid,
              fillColor={127,127,127}),
            Rectangle(extent={{-20,40},{0,-60}}, lineColor={0,0,0}),
            Polygon(
              points={{0,60},{20,60},{0,40},{-20,40},{0,60}},
              lineColor={0,0,0},
              smooth=Smooth.None),
            Line(
              points={{-20,0},{-100,0}},
              color={0,0,240},
              visible=inclFacesX,
              thickness=0.5),
            Line(
              points={{10,0},{100,0}},
              color={0,0,240},
              visible=inclFacesX,
              thickness=0.5),
            Line(
              points={{0,-60},{0,-100}},
              color={0,0,240},
              visible=inclFacesY,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{0,100},{0,50}},
              color={0,0,240},
              visible=inclFacesY,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{-50,-50},{-10,-10}},
              color={0,0,240},
              visible=inclFacesZ,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{20,20},{50,50}},
              color={0,0,240},
              visible=inclFacesZ,
              smooth=Smooth.None,
              thickness=0.5),
            Text(
              extent={{-100,60},{100,100}},
              textString="%name",
              visible=not inclFacesY,
              lineColor={0,0,0})}));

    end CaGDL;

    model Sigracet10BA "<html>SGL Carbon Group Sigracet&reg; 10 BA</html>"
      extends CaGDL(
        L_x={0.400*U.mm},
        epsilon=0.88,
        Subregion(graphite('e-'(sigma=U.S/(0.012*U.cm)))));
      // Additional properties not incorporated [SGL2007]:
      //     Diffusivity:  L = 0.400 mm, p = 0.85 m/s (for air) => D = P*L = 340 mm2/s
      //     Density:  (85 g/m2)/(0.400 mm)/0.88 = 212.5 kg/m3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">SGL2007</a>].</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end Sigracet10BA;

    model Sigracet10BB "<html>SGL Carbon Group Sigracet&reg; 10 BB</html>"
      extends CaGDL(
        L_x={0.420*U.mm},
        epsilon=0.84,
        Subregion(graphite('e-'(sigma=42*U.S/(15*U.cm)))));
      // Additional properties not incorporated [SGL2007]:
      //     Diffusivity:  L = 0.420 mm, p = 0.03 m/s (for air) => D = P*L = 12.6 mm2/s
      //     Density:  (125 g/m2)/(0.420 mm) = 297.62 kg/m3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">SGL2007</a>].</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end Sigracet10BB;

    model Sigracet10BC "<html>SGL Carbon Group Sigracet&reg; 10 BC</html>"
      extends CaGDL(
        L_x={0.420*U.mm},
        epsilon=0.82,
        Subregion(graphite('e-'(sigma=42*U.S/(16*U.cm)))));
      // Additional properties not incorporated [SGL2007]:
      //     Diffusivity:  L = 0.420 mm, p = 0.0145 m/s (for air) => D = P*L = 6.09 mm2/s
      //     Density:  (135 g/m2)/(0.420 mm) = 321.43 kg/m3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">SGL2007</a>].</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end Sigracet10BC;

    model Sigracet24BA "<html>SGL Carbon Group Sigracet&reg; 24 BA</html>"
      extends CaGDL(
        L_x={0.190*U.mm},
        epsilon=0.84,
        Subregion(graphite('e-'(sigma=19*U.S/(10*U.cm)))));
      // Additional properties not incorporated [SGL2004]:
      //     Diffusivity:  L = 0.190 mm, p = 0.30 m/s (for air) => D = P*L = 57 mm2/s
      //     Density:  (54 g/m2)/(0.190 mm) = 284.21 kg/m3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">SGL2004</a>].</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end Sigracet24BA;

    model Sigracet24BC "<html>SGL Carbon Group Sigracet&reg; 24 BC</html>"
      extends CaGDL(
        L_x={0.235*U.mm},
        epsilon=0.76,
        Subregion(graphite('e-'(sigma=23.5*U.S/(11*U.cm)))));
      // Additional properties not incorporated [SGL2004]:
      //     Diffusivity:  L = 0.235 mm, p = 0.0045 m/s (for air) => D = P*L = 1.0575 mm2/s
      //     Density:  (100 g/m2)/(0.235 mm) = 425.53 kg/m3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">SGL2004</a>].</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end Sigracet24BC;

    model Sigracet25BA "<html>SGL Carbon Group Sigracet&reg; 25 BA</html>"
      extends CaGDL(
        L_x={0.190*U.mm},
        epsilon=0.88,
        Subregion(graphite('e-'(sigma=19*U.S/(10*U.cm)))));
      // Additional properties not incorporated [SGL2004]:
      //     Diffusivity:  L = 0.190 mm, p = 0.90 m/s (for air) => D = P*L = 171 mm2/s
      //     Density:  (40 g/m2)/(0.190 mm) = 210.53 kg/m3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">SGL2004</a>].</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end Sigracet25BA;

    model Sigracet25BC "<html>SGL Carbon Group Sigracet&reg; 25 BC</html>"
      extends CaGDL(
        L_x={0.235*U.mm},
        epsilon=0.80,
        Subregion(graphite('e-'(sigma=23.5*U.S/(12*U.cm)))));
      // Additional properties not incorporated [SGL2004]:
      //     Diffusivity:  L = 0.235 mm, p = 0.008 m/s (for air) => D = P*L = 1.88 mm2/s
      //     Density:  (86 g/m2)/(0.235 mm) = 365.96 kg/m3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">SGL2004</a>].</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end Sigracet25BC;

    model TorayTGPH030 "Toray Industries TGP-H-030"
      extends CaGDL(
        L_x={0.11*U.mm},
        epsilon=0.80,
        Subregion(graphite('C+'(theta=U.m*U.K/(1.7*U.W)),'e-'(sigma=U.S/(0.80*U.mm)))));
      // Additional properties not incorporated [Toray2010]:
      //     Diffusivity:  L = 0.110 mm, P/p = 2500 ml.mm/(cm2.hr.mmAq) = 0.70814e-3 m/(s.kPa)
      //         => D = P*L = 7.89e-6 m2/s, assuming p = 101.325 kPa
      //     Bulk density:  0.44 g/cm3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">Toray2010</a>].  The electronic mobility (<code>mu</code>) is based on the
  through-plane value of resistivity.  The thermal conductivity is not listed but is
  taken to be the same as for TGP-H-060, TGP-H-090, and TGP-H-120.</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end TorayTGPH030;

    model TorayTGPH060 "Toray Industries TGP-H-060"
      extends CaGDL(
        L_x={0.19*U.mm},
        epsilon=0.78,
        Subregion(graphite('C+'(theta=U.m*U.K/(1.7*U.W)),'e-'(sigma=U.S/(0.80*U.mm)))));

      // Additional properties not incorporated [Toray2010]:
      //     Diffusivity: L = 0.190 mm, P/p = 1900 ml.mm/(cm2.hr.mmAq) = 0.53818e-3 m/(s.kPa)
      //         => D = P*L = 10.36e-6 m2/s, assuming p = 101.325 kPa
      //     Bulk density:  0.44 g/cm3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">Toray2010</a>].  The electronic mobility (<code>mu</code>) is based on the
  through-plane value of resistivity.  The thermal conductivity is through the plane at
  room temperature.</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end TorayTGPH060;

    model TorayTGPH090 "Toray Industries TGP-H-090"
      extends CaGDL(
        L_x={0.28*U.mm},
        epsilon=0.78,
        Subregion(graphite('C+'(theta=U.m*U.K/(1.7*U.W)),'e-'(sigma=U.S/(0.80*U.mm)))));

      // Additional properties not incorporated [Toray2010]:
      //     Diffusivity: L = 0.280 mm, P/p = 1700 ml.mm/(cm2.hr.mmAq) = 0.48153e-3 m/(s.kPa)
      //         => D = P*L = 13.66e-6 m2/s, assuming p = 101.325 kPa
      //     Bulk density:  0.44 g/cm3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">Toray2010</a>].  The electronic mobility (<code>mu</code>) is based on the
  through-plane value of resistivity.  The thermal conductivity is through the plane at
  room temperature.</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end TorayTGPH090;

    model TorayTGPH120 "Toray Industries TGP-H-030"
      extends CaGDL(
        L_x={0.37*U.mm},
        epsilon=0.78,
        Subregion(graphite('C+'(theta=U.m*U.K/(1.7*U.W)),'e-'(sigma=U.S/(0.80*U.mm)))));

      // Additional properties not incorporated [Toray2010]:
      //     Diffusivity: L = 0.370 mm, P/p = 1500 ml.mm/(cm2.hr.mmAq) = 0.42488e-3 m/(s.kPa)
      //         => D = P*L = 15.93e-6 m2/s, assuming p = 101.325 kPa
      //     Bulk density:  0.44 g/cm3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">Toray2010</a>].  The electronic mobility (<code>mu</code>) is based on the
  through-plane value of resistivity.  The thermal conductivity is through the plane at
  room temperature.</p>

  <p>For more information, see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end TorayTGPH120;

  end CaGDLs;

  package CaFPs "Cathode flow plates"
    extends Modelica.Icons.Package;

    model CaFP "Cathode flow plate"
      // extends FCSys.BaseClasses.Icons.Names.Top4;

      extends Region(
        L_x={8*U.mm},
        L_y={2*U.m},
        L_z={5*U.mm},
        final inclFacesX=true,
        final inclFacesY=true,
        inclFacesZ=false,
        redeclare replaceable model Subregion =
            FCSys.Subregions.SubregionNoIonomer (
            inclTransX=true,
            inclTransY=true,
            inclTransZ=false,
            gas(
              reduceVel=true,
              reduceTemp=true,
              k_T=fill(epsilon^(3/2), 3),
              inclH2O=true,
              inclN2=true,
              inclO2=true,
              H2O(p_IC=environment.p_H2O, consTransX=Conservation.IC),
              N2(
                initTransX=InitTranslational.None,
                initTransY=InitTranslational.None,
                initTransZ=InitTranslational.None,
                initEnergy=InitScalar.None,
                p_IC=(1 - environment.n_O2)*(environment.p - environment.p_H2O)),

              O2(
                initTransX=InitTranslational.None,
                initTransY=InitTranslational.None,
                initTransZ=InitTranslational.None,
                initEnergy=InitScalar.None,
                p_IC=environment.n_O2*(environment.p - environment.p_H2O))),
            graphite(
              reduceTemp=true,
              k_T=fill((1 - epsilon)^(3/2), 3),
              'inclC+'=true,
              'incle-'=true,
              'C+'(V_IC=V - epsilonV,theta=U.m*U.K/(95*U.W)),
              'e-'(
                consTransY=Conservation.IC,
                initTransX=InitTranslational.None,
                initEnergy=InitScalar.None,
                rho_IC=17842.7*U.C/U.cc,
                sigma=U.S/(1.470e-3*U.cm))),
            liquid(
              k_T=fill(epsilon^(3/2), 3),
              inclH2O=true,
              H2O(consTransX=Conservation.IC,V_IC=Modelica.Constants.eps*U.cc))))
        annotation (IconMap(primitivesVisible=false));

      // See the documentation layer of Subregions.Phases.BaseClasses.EmptyPhase
      // regarding the settings of k_T for each phase.

      parameter Q.NumberAbsolute epsilon(nominal=1) = 0.1
        "<html>Volumetric porosity (&epsilon;)</html>"
        annotation (group="Geometry");

      // Auxiliary variables (for analysis)
      output Q.NumberAbsolute RH[n_x, n_y, n_z](
        each stateSelect=StateSelect.never,
        each displayUnit="%") =
        Modelica.Media.Air.MoistAir.saturationPressureLiquid(subregions.gas.H2O.T
        /U.K) ./ subregions.gas.dalton.p if environment.analysis and
        hasSubregions "Relative humidity";
      output Q.Number n_O2[n_x, n_y, n_z](each stateSelect=StateSelect.never)
         = subregions.gas.O2.N ./ (subregions.gas.N2.N + subregions.gas.O2.N)
        if environment.analysis and hasSubregions "Dry-gas concentration of O2";

    protected
      final parameter Q.Volume epsilonV=epsilon*V "Fluid volume";

      outer Conditions.Environment environment "Environmental conditions";

      // See AnFPs.AnFP for data on additional materials.
      annotation (Documentation(info="<html>
<p>This model represents the cathode flow plate of a PEMFC.
The x axis extends from the anode to the cathode.
The y axis extends along the length of the channel. The model is
bidirectional, so that either <code>yNegative</code> or <code>yPositive</code> can be
used as the inlet. The z axis extends across the width of the channel.</p>

<p>By default, the cross-sectional area in the xy plane is 100 cm<sup>2</sup>.</p>

<p>See <a href=\"modelica://FCSys.Species.'C+'.Graphite.Fixed\">Subregions.Species.'C+'.Graphite.Fixed</a>
regarding the default specific heat capacity.  The default thermal resistivity
of the carbon (<code>theta = U.m*U.K/(95*U.W)</code>) and the
electrical conductivity (<code>sigma=U.S/(1.470e-3*U.cm)</code>)
are that of Entegris/Poco Graphite AXF-5Q
[<a href=\"modelica://FCSys.UsersGuide.References\">Entegris2012</a>].
  There is additional data in the
text layer of the <a href=\"modelica://FCSys.Regions.AnFPs.AnFP\">AnFP</a> model.</p>

<p>For more information, see the
    <a href=\"modelica://FCSys.Regions.Region\">Region</a> model.</p>
</html>"), Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.1), graphics={
            Rectangle(
              extent={{-100,60},{100,100}},
              fillColor={255,255,255},
              visible=not inclFacesY,
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Rectangle(
              extent={{-76.648,66.211},{-119.073,52.0689}},
              fillPattern=FillPattern.HorizontalCylinder,
              rotation=45,
              fillColor={135,135,135},
              origin={111.017,77.3801},
              pattern=LinePattern.None,
              lineColor={95,95,95}),
            Rectangle(
              extent={{-20,40},{0,-60}},
              lineColor={95,95,95},
              fillPattern=FillPattern.VerticalCylinder,
              fillColor={135,135,135}),
            Polygon(
              points={{20,0},{42,0},{42,80},{-42,80},{-42,0},{-20,0},{-20,40},{
                  0,60},{20,60},{20,0}},
              smooth=Smooth.None,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Polygon(
              points={{20,0},{42,0},{42,-80},{-42,-80},{-42,0},{-20,0},{-20,-60},
                  {0,-60},{20,-40},{20,0}},
              smooth=Smooth.None,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),
            Rectangle(extent={{-20,40},{0,-60}}, lineColor={0,0,0}),
            Polygon(
              points={{-20,40},{0,60},{20,60},{0,40},{-20,40}},
              lineColor={0,0,0},
              smooth=Smooth.None),
            Polygon(
              points={{20,60},{0,40},{0,-60},{20,-40},{20,60}},
              lineColor={0,0,0},
              fillColor={95,95,95},
              fillPattern=FillPattern.Solid),
            Line(
              points={{-20,0},{-100,0}},
              color={0,0,240},
              visible=inclFacesX,
              thickness=0.5),
            Line(
              points={{10,0},{100,0}},
              color={127,127,127},
              visible=inclFacesX,
              thickness=0.5),
            Ellipse(
              extent={{-4,52},{4,48}},
              lineColor={135,135,135},
              fillColor={0,0,240},
              visible=inclFacesY,
              fillPattern=FillPattern.Sphere),
            Line(
              points={{0,-60},{0,-100}},
              color={0,0,240},
              visible=inclFacesY,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{0,100},{0,50}},
              color={0,0,240},
              visible=inclFacesY,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{-50,-50},{-10,-10}},
              color={0,0,240},
              visible=inclFacesZ,
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{20,20},{50,50}},
              color={0,0,240},
              visible=inclFacesZ,
              smooth=Smooth.None,
              thickness=0.5),
            Text(
              extent={{-100,60},{100,100}},
              textString="%name",
              visible=not inclFacesY,
              lineColor={0,0,0})}));

    end CaFP;

  end CaFPs;


  partial model Region "Base model for a 3D array of subregions"
    import FCSys.BaseClasses.Utilities.cartWrap;
    // extends FCSys.BaseClasses.Icons.Names.Top3;
    // extends FCSys.BaseClasses.Icons.Names.Top6;

    // Geometric parameters
    parameter Q.Length L_x[:](each final min=Modelica.Constants.small) = {U.cm}
      "<html>Lengths along the x axis (<i>L</i><sub>x</sub>)</html>"
      annotation (Dialog(group="Geometry"));
    parameter Q.Length L_y[:](each final min=Modelica.Constants.small) = {U.cm}
      "<html>Lengths along the y axis (<i>L</i><sub>y</sub>)</html>"
      annotation (Dialog(group="Geometry"));
    parameter Q.Length L_z[:](each final min=Modelica.Constants.small) = {U.cm}
      "<html>Lengths along the z axis (<i>L</i><sub>z</sub>)</html>"
      annotation (Dialog(group="Geometry"));
    final parameter Integer n_x=size(L_x, 1)
      "Number of sets of subregions along the x axis"
      annotation (HideResult=true);
    final parameter Integer n_y=size(L_y, 1)
      "Number of sets of subregions along the y axis"
      annotation (HideResult=true);
    final parameter Integer n_z=size(L_z, 1)
      "Number of sets of subregions along the z axis"
      annotation (HideResult=true);

    // Assumptions about included faces
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

    // Auxiliary parameters (for analysis only)
    final parameter Q.Length L[Axis]={sum(L_x),sum(L_y),sum(L_z)} if
      hasSubregions "Length";
    final parameter Q.Area A[Axis]={L[cartWrap(axis + 1)]*L[cartWrap(axis + 2)]
        for axis in Axis} if hasSubregions "Cross-sectional areas";
    final parameter Q.Volume V=product(L) if hasSubregions "Volume";

    replaceable model Subregion = Subregions.Subregion constrainedby
      FCSys.Subregions.BaseClasses.EmptySubregion(
      final inclFacesX=inclFacesX,
      final inclFacesY=inclFacesY,
      final inclFacesZ=inclFacesZ) "Base subregion model" annotation (
        __Dymola_choicesAllMatching=true);
    // The inclFacesX, inclFacesX, inclFacesX parameters still appear in the
    // parameter dialog even though they are final.
    // TODO:  Find a way to hide them.
    Subregion subregions[n_x, n_y, n_z](final L={{L_x[i_x],L_y[i_y],L_z[i_z]}
          for i_z in 1:n_z, i_y in 1:n_y, i_x in 1:n_x}) if hasSubregions
      "Instances of the subregion model"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    Connectors.FaceBus xNegative[n_y, n_z] if inclFacesX
      "Negative face along the x axis" annotation (Placement(transformation(
            extent={{-50,-10},{-30,10}}), iconTransformation(extent={{-110,-10},
              {-90,10}})));
    Connectors.FaceBus xPositive[n_y, n_z] if inclFacesX
      "Positive face along the x axis" annotation (Placement(transformation(
            extent={{30,-10},{50,10}}), iconTransformation(extent={{90,-10},{
              110,10}})));
    Connectors.FaceBus yNegative[n_x, n_z] if inclFacesY
      "Negative face along the y axis" annotation (Placement(transformation(
            extent={{-10,-50},{10,-30}}), iconTransformation(extent={{-10,-110},
              {10,-90}})));
    Connectors.FaceBus yPositive[n_x, n_z] if inclFacesY
      "Positive face along the y axis" annotation (Placement(transformation(
            extent={{-10,30},{10,50}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.FaceBus zNegative[n_x, n_y] if inclFacesZ
      "Negative face along the z axis" annotation (Placement(transformation(
            extent={{10,10},{30,30}}), iconTransformation(extent={{40,40},{60,
              60}})));
    Connectors.FaceBus zPositive[n_x, n_y] if inclFacesZ
      "Positive face along the z axis" annotation (Placement(transformation(
            extent={{-30,-30},{-10,-10}}), iconTransformation(extent={{-60,-60},
              {-40,-40}})));

  protected
    final parameter Boolean hasSubregions=n_x > 0 and n_y > 0 and n_z > 0
      "true, if there are any subregions";

  equation
    // X axis
    connect(xNegative, subregions[1, :, :].xNegative) annotation (Line(
        points={{-40,5.55112e-16},{-38,5.55112e-16},{-40,0},{-30,0},{-30,
            6.10623e-16},{-10,6.10623e-16}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    for i in 1:n_x - 1 loop
      connect(subregions[i, :, :].xPositive, subregions[i + 1, :, :].xNegative)
        "Connection b/w neighboring subregions (not shown the diagram)";
    end for;
    connect(subregions[n_x, :, :].xPositive, xPositive) annotation (Line(
        points={{10,6.10623e-16},{20,6.10623e-16},{20,0},{30,0},{30,5.55112e-16},
            {40,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    if n_x == 0 then
      connect(xNegative, xPositive)
        "Direct pass-through (not shown the diagram)";
    end if;

    // Y axis
    connect(yNegative, subregions[:, 1, :].yNegative) annotation (Line(
        points={{5.55112e-16,-40},{0,-10},{6.10623e-16,-10}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    for i in 1:n_y - 1 loop
      connect(subregions[:, i, :].yPositive, subregions[:, i + 1, :].yNegative)
        "Connection b/w neighboring subregions (not shown the diagram)";
    end for;
    connect(subregions[:, n_y, :].yPositive, yPositive) annotation (Line(
        points={{6.10623e-16,10},{6.10623e-16,24},{0,24},{0,40},{5.55112e-16,40}},

        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));

    if n_y == 0 then
      connect(yNegative, yPositive)
        "Direct pass-through (not shown the diagram)";
    end if;

    // Z axis
    connect(zNegative, subregions[:, :, 1].zNegative) annotation (Line(
        points={{20,20},{5,5}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    for i in 1:n_z - 1 loop
      connect(subregions[:, :, i].zPositive, subregions[:, :, i + 1].zNegative)
        "Connection b/w neighboring subregions (not shown the diagram)";
    end for;
    connect(zPositive, subregions[:, :, n_z].zPositive) annotation (Line(
        points={{-20,-20},{-5,-5}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    if n_z == 0 then
      connect(zNegative, zPositive)
        "Direct pass-through (not shown the diagram)";
    end if;
    // TODO:  Once primitivesVisible is supported by Modelica tools (not
    // supported by Dymola 7.4 or 2013 FD01), complete the icon of this model.
    // Until then, the icon should be blank so that the fuel cell layer models
    // can remain as they are.
    annotation (
      Documentation(info="<html>
  <p>If <code>L_x</code> is an empty vector (e.g., <code>zeros(0)</code>,
  <code>ones(0)</code>, or <code>fill(1, 0)</code>), then there are no
  subregions along the x axis and the boundaries along the x axis are
  directly connected.  The same applies to the other axes.</li></p>
</html>"),
      Placement(transformation(extent={{-10,-10},{10,10}})),
      Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
              100,100}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Text(
            extent={{-100,120},{100,160}},
            textString="%name",
            visible=inclFacesY,
            lineColor={0,0,0}), Text(
            extent={{-100,56},{100,96}},
            textString="%name",
            visible=not inclFacesY,
            lineColor={0,0,0})}));
  end Region;
  annotation (Documentation(info="<html>
<p><b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
Copyright 2007&ndash;2013, Georgia Tech Research Corporation.</p>

<p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.ModelicaLicense2\">
FCSys.UsersGuide.ModelicaLicense2</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p>
</html>"));

end Regions;
