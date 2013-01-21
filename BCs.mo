within FCSys;
package BCs "Models for boundary conditions"
  extends Modelica.Icons.SourcesPackage;

  package Examples "Examples and tests"
    extends Modelica.Icons.ExamplesPackage;

    model FaceBC "<html>Test the BCs for the face of a subregion</html>"
      extends Modelica.Icons.Example;

      FaceBus.Subregion subregionFaceBC(gas(inclH2O=true, H2O(
            redeclare Face.Material.Current normal(spec(k=0)),
            inviscidX=false,
            inviscidZ=false)))
        annotation (Placement(transformation(extent={{-10,14},{10,34}})));
      Subregions.Subregion subregion(
        L={1,1,1}*U.cm,
        inclReact=false,
        inclFacesX=false,
        inclFacesY=true,
        inclFacesZ=false,
        inclLinX=false,
        inclLinY=true,
        graphite('inclC+'=true, 'C+'(V_IC=0.5*U.cc)),
        gas(inclH2O=true, H2O(
            xNegative(isobaric=true, inviscidY=true),
            xPositive(isobaric=true, inviscidY=true),
            zNegative(inviscidY=true),
            zPositive(inviscidY=true),
            yPositive(
              isobaric=false,
              inviscidX=false,
              inviscidZ=false))))
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      inner BCs.Environment environment
        annotation (Placement(transformation(extent={{30,30},{50,50}})));
    equation
      connect(subregion.yPositive, subregionFaceBC.face) annotation (Line(
          points={{6.10623e-16,10},{-3.36456e-22,16},{6.10623e-16,16},{
              6.10623e-16,20}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        experiment(NumberOfIntervals=5000),
        experimentSetupOutput,
        Commands(file="resources/scripts/Dymola/BCs.Examples.FaceBC.mos"));
    end FaceBC;

    model FaceBCPhases
      "<html>Test the BCs for the face of a subregion with phases</html>"
      extends Modelica.Icons.Example;

      // Geometric parameters
      inner parameter Q.Length L[Axis](each min=Modelica.Constants.small,start=
            ones(3)*U.cm) "<html>Length (<b>L</b>)</html>"
        annotation (Dialog(group="Geometry"));
      final inner parameter Q.Area A[Axis]={L[cartWrap(axis + 1)]*L[cartWrap(
          axis + 2)] for axis in Axis} "Cross-sectional area";

      FaceBus.Phases.Gas phaseFaceBC(
        inclH2O=true,
        H2O(isobaric=false, redeclare Face.Material.Current normal(spec(k=0))),

        axis=FCSys.BaseClasses.Axis.y)
        annotation (Placement(transformation(extent={{-10,14},{10,34}})));

      Subregions.Volume volume
        annotation (Placement(transformation(extent={{-16,-16},{16,16}})));
      FCSys.Subregions.Phases.Gas gas(
        inclReact=false,
        inclLin={false,true,false},
        inclH2=false,
        inclH2O=true,
        H2O(
          xNegative(isobaric=true, inviscidY=true),
          xPositive(isobaric=true, inviscidY=true),
          yPositive(isobaric=false),
          zNegative(inviscidY=true),
          zPositive(inviscidY=true)))
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      inner BCs.Environment environment
        annotation (Placement(transformation(extent={{30,30},{50,50}})));
    equation
      connect(gas.yPositive, phaseFaceBC.face) annotation (Line(
          points={{6.10623e-16,10},{-3.36456e-22,16},{6.10623e-16,16},{
              6.10623e-16,20}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(gas.inert, volume.inert) annotation (Line(
          points={{8,-8},{11,-11}},
          color={72,90,180},
          smooth=Smooth.None));
      annotation (
        experiment,
        experimentSetupOutput,
        Diagram(graphics));
    end FaceBCPhases;

    model Router "<html>Test the <code>Router<code> model</html>"
      extends Modelica.Icons.Example;

      // TODO:  Make this into a meaningful example.
      Router router
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    end Router;

    model AnodeAdapter "<html>Test the <code>'Adapte-'</code> model</html>"

      extends Modelica.Icons.Example;
      extends Modelica.Icons.UnderConstruction;
      /* **fails simulation:
  "Model error - power: (1/subregion.graphite.'e-'.rho_face[1, 2]) ** (-1.66666666666667) = (-0.000669271) ** (-1.66667)

  Non-linear solver will attempt to handle this problem."
  */
      inner Modelica.Fluid.System system(T_ambient=293.15 + 5)
        annotation (Placement(transformation(extent={{40,70},{60,90}})));
      inner BCs.Environment environment(T=350*U.K)
        annotation (Placement(transformation(extent={{70,72},{90,92}})));
      Subregions.SubregionNoIonomer subregion(
        L={1,1,1}*U.cm,
        inclReact=false,
        inclFacesY=false,
        inclFacesZ=false,
        gas(inclH2=true, inclH2O=true),
        graphite('inclC+'=true, 'incle-'=true),
        liquid(inclH2O=true))
        annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));
      Adapters.Anode anodeAdapter(redeclare package LiquidMedium =
            Modelica.Media.CompressibleLiquids.LinearColdWater)
        annotation (Placement(transformation(extent={{-20,-10},{0,10}})));
      Modelica.Electrical.Analog.Basic.Ground ground
        annotation (Placement(transformation(extent={{0,-40},{20,-20}})));

      Modelica.Fluid.Vessels.ClosedVolume gasVolume(
        use_portsData=false,
        nPorts=1,
        V=1e-6,
        use_HeatTransfer=true,
        redeclare
          Modelica.Fluid.Vessels.BaseClasses.HeatTransfer.IdealHeatTransfer
          HeatTransfer,
        redeclare package Medium = Adapters.Media.AnodeGas,
        medium(p(fixed=true),X(each fixed=true)))
        annotation (Placement(transformation(extent={{30,20},{50,40}})));
      Modelica.Fluid.Vessels.ClosedVolume liquidVolume(
        nPorts=1,
        use_HeatTransfer=true,
        redeclare
          Modelica.Fluid.Vessels.BaseClasses.HeatTransfer.IdealHeatTransfer
          HeatTransfer,
        V=0.5e-6,
        use_portsData=false,
        redeclare package Medium =
            Modelica.Media.CompressibleLiquids.LinearColdWater,
        medium(p(fixed=true),T(fixed=true)))
        annotation (Placement(transformation(extent={{70,20},{90,40}})));

    equation
      connect(ground.p, anodeAdapter.pin) annotation (Line(
          points={{10,-20},{10,2},{-2,2}},
          color={0,0,255},
          smooth=Smooth.None));

      connect(subregion.xPositive, anodeAdapter.face) annotation (Line(
          points={{-30,6.10623e-16},{-24,-3.36456e-22},{-24,6.10623e-16},{-18,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(gasVolume.heatPort, anodeAdapter.heatPort) annotation (Line(
          points={{30,30},{20,30},{20,-2},{-2,-2}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(gasVolume.ports[1], anodeAdapter.gasPort) annotation (Line(
          points={{40,20},{40,6},{-2,6}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(liquidVolume.heatPort, anodeAdapter.heatPort) annotation (Line(
          points={{70,30},{60,30},{60,-2},{-2,-2}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(anodeAdapter.liquidPort, liquidVolume.ports[1]) annotation (Line(
          points={{-2,-6},{80,-6},{80,20}},
          color={0,127,255},
          smooth=Smooth.None));
      annotation (
        experiment(StopTime=2e-10),
        experimentSetupOutput,
        Commands(file="resources/scripts/Dymola/BCs.Examples.Adapteminus.mos"),

        Diagram(graphics));
    end AnodeAdapter;
  end Examples;

  package Adapters
    "<html>Adapters to <a href=\"modelica://Modelica\">Package Modelica</a></html>"
    extends Modelica.Icons.Package;

    model Anode
      "<html>Adapter between <a href=\"modelica://Modelica\">Modelica</a> and the face connector of a <a href=\"modelica://FCSys.Assemblies.Cells.Cell\">Cell</a>, <a href=\"modelica://FCSys.Regions.Region\">Region</a>, or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a></html>"
      extends FCSys.BaseClasses.Icons.Names.Top4;

      replaceable package GasMedium = Media.AnodeGas constrainedby
        Modelica.Media.Interfaces.PartialMedium "Medium model for the gas"
        annotation (choicesAllMatching=true, Dialog(group="Material properties"));
      replaceable package LiquidMedium =
          Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby
        Modelica.Media.Interfaces.PartialMedium "Medium model for the liquid"
        annotation (choicesAllMatching=true, Dialog(group="Material properties"));

      FCSys.Connectors.FaceBus face
        "Multi-species connector for linear momentum and heat" annotation (
          Placement(transformation(extent={{-90,-10},{-70,10}}),
            iconTransformation(extent={{-90,-10},{-70,10}})));
      Modelica.Fluid.Interfaces.FluidPort_b gasPort(redeclare final package
          Medium = GasMedium) "Modelica fluid port for the gas" annotation (
          Placement(transformation(extent={{70,50},{90,70}}),
            iconTransformation(extent={{70,50},{90,70}})));
      Modelica.Electrical.Analog.Interfaces.NegativePin pin
        "Modelica electrical pin" annotation (Placement(transformation(extent={
                {70,10},{90,30}}), iconTransformation(extent={{70,10},{90,30}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
        "Modelica heat port" annotation (Placement(transformation(extent={{70,-30},
                {90,-10}}),iconTransformation(extent={{70,-30},{90,-10}})));
      Modelica.Fluid.Interfaces.FluidPort_b liquidPort(redeclare final package
          Medium = LiquidMedium) "Modelica fluid port for the liquid"
        annotation (Placement(transformation(extent={{70,-70},{90,-50}}),
            iconTransformation(extent={{70,-70},{90,-50}})));
      Phases.AnodeGas gas(redeclare final package Medium = GasMedium)
        "Gas subadapter"
        annotation (Placement(transformation(extent={{-10,30},{10,50}})));
      Phases.Graphite graphite "Graphite subadapter"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      Phases.Liquid liquid(redeclare final package Medium = LiquidMedium)
        "Liquid subadapter"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
    equation
      connect(gas.face, face.gas) annotation (Line(
          points={{-8,40},{-40,40},{-40,5.55112e-16},{-80,5.55112e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(gasPort, gas.fluidPort) annotation (Line(
          points={{80,60},{50,60},{50,36},{8,36}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(gas.heatPort, heatPort) annotation (Line(
          points={{8,40},{30,40},{30,-20},{80,-20}},
          color={191,0,0},
          smooth=Smooth.None));

      connect(graphite.face, face.graphite) annotation (Line(
          points={{-8,6.10623e-16},{-40,6.10623e-16},{-40,5.55112e-16},{-80,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None,
          thickness=0.5));

      connect(graphite.pin, pin) annotation (Line(
          points={{8,4},{50,4},{50,20},{80,20}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(graphite.heatPort, heatPort) annotation (Line(
          points={{8,6.10623e-16},{30,6.10623e-16},{30,-20},{80,-20}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(liquid.face, face.liquid) annotation (Line(
          points={{-8,-40},{-40,-40},{-40,5.55112e-16},{-80,5.55112e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(liquidPort, liquid.fluidPort) annotation (Line(
          points={{80,-60},{50,-60},{50,-44},{8,-44}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(liquid.heatPort, heatPort) annotation (Line(
          points={{8,-40},{30,-40},{30,-20},{80,-20}},
          color={191,0,0},
          smooth=Smooth.None));

      annotation (Icon(graphics={Line(
                  points={{0,60},{0,-60}},
                  color={0,0,0},
                  smooth=Smooth.None,
                  pattern=LinePattern.Dash,
                  thickness=0.5),Line(
                  points={{0,0},{-80,0}},
                  color={127,127,127},
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{0,20},{80,20}},
                  color={0,0,255},
                  smooth=Smooth.None),Line(
                  points={{0,-20},{80,-20}},
                  color={191,0,0},
                  smooth=Smooth.None),Line(
                  points={{0,60},{80,60}},
                  color={0,127,255},
                  smooth=Smooth.None),Line(
                  points={{0,-60},{80,-60}},
                  color={0,127,255},
                  smooth=Smooth.None)}));
    end Anode;

    model Cathode
      "<html>Adapter between <a href=\"modelica://Modelica\">Modelica</a> and the face connector of a <a href=\"modelica://FCSys.Assemblies.Cells.Cell\">Cell</a>, <a href=\"modelica://FCSys.Regions.Region\">Region</a>, or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a></html>"
      extends FCSys.BaseClasses.Icons.Names.Top4;

      replaceable package GasMedium = Adapters.Media.CathodeGas constrainedby
        Modelica.Media.Interfaces.PartialMedium "Medium model for the gas"
        annotation (choicesAllMatching=true, Dialog(group="Material properties"));
      replaceable package LiquidMedium =
          Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby
        Modelica.Media.Interfaces.PartialMedium "Medium model for the liquid"
        annotation (choicesAllMatching=true, Dialog(group="Material properties"));

      FCSys.Connectors.FaceBus face
        "Multi-species connector for linear momentum and heat" annotation (
          Placement(transformation(extent={{-90,-10},{-70,10}}),
            iconTransformation(extent={{-90,-10},{-70,10}})));
      Modelica.Fluid.Interfaces.FluidPort_b gasPort(redeclare final package
          Medium = GasMedium) "Modelica fluid port for the gas" annotation (
          Placement(transformation(extent={{70,50},{90,70}}),
            iconTransformation(extent={{70,50},{90,70}})));
      Modelica.Electrical.Analog.Interfaces.NegativePin pin
        "Modelica electrical pin" annotation (Placement(transformation(extent={
                {70,10},{90,30}}), iconTransformation(extent={{70,10},{90,30}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
        "Modelica heat port" annotation (Placement(transformation(extent={{70,-30},
                {90,-10}}),iconTransformation(extent={{70,-30},{90,-10}})));
      Modelica.Fluid.Interfaces.FluidPort_b liquidPort(redeclare final package
          Medium = LiquidMedium) "Modelica fluid port for the liquid"
        annotation (Placement(transformation(extent={{70,-70},{90,-50}}),
            iconTransformation(extent={{70,-70},{90,-50}})));

      Phases.CathodeGas gas(redeclare final package Medium = GasMedium)
        "Gas subadapter"
        annotation (Placement(transformation(extent={{-10,30},{10,50}})));
      Phases.Graphite graphite "Graphite subadapter"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Phases.Liquid liquid(redeclare final package Medium = LiquidMedium)
        "Liquid subadapter"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
    equation
      connect(gas.face, face.gas) annotation (Line(
          points={{-8,40},{-40,40},{-40,5.55112e-16},{-80,5.55112e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(gasPort, gas.fluidPort) annotation (Line(
          points={{80,60},{50,60},{50,36},{8,36}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(gas.heatPort, heatPort) annotation (Line(
          points={{8,40},{30,40},{30,-20},{80,-20}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(graphite.face, face.graphite) annotation (Line(
          points={{-8,6.10623e-16},{-40,6.10623e-16},{-40,5.55112e-16},{-80,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None,
          thickness=0.5));

      connect(graphite.pin, pin) annotation (Line(
          points={{8,4},{50,4},{50,20},{80,20}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(graphite.heatPort, heatPort) annotation (Line(
          points={{8,6.10623e-16},{30,6.10623e-16},{30,-20},{80,-20}},
          color={191,0,0},
          smooth=Smooth.None));
      connect(liquid.face, face.liquid) annotation (Line(
          points={{-8,-40},{-40,-40},{-40,5.55112e-16},{-80,5.55112e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(liquidPort, liquid.fluidPort) annotation (Line(
          points={{80,-60},{50,-60},{50,-44},{8,-44}},
          color={0,127,255},
          smooth=Smooth.None));
      connect(liquid.heatPort, heatPort) annotation (Line(
          points={{8,-40},{30,-40},{30,-20},{80,-20}},
          color={191,0,0},
          smooth=Smooth.None));

      annotation (Icon(graphics={Line(
                  points={{0,60},{0,-60}},
                  color={0,0,0},
                  smooth=Smooth.None,
                  pattern=LinePattern.Dash,
                  thickness=0.5),Line(
                  points={{0,0},{-80,0}},
                  color={127,127,127},
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{0,20},{80,20}},
                  color={0,0,255},
                  smooth=Smooth.None),Line(
                  points={{0,-20},{80,-20}},
                  color={191,0,0},
                  smooth=Smooth.None),Line(
                  points={{0,60},{80,60}},
                  color={0,127,255},
                  smooth=Smooth.None),Line(
                  points={{0,-60},{80,-60}},
                  color={0,127,255},
                  smooth=Smooth.None)}));
    end Cathode;

    package Phases "Adapters for material phases"
      extends Modelica.Icons.Package;

      model AnodeGas
        "<html>Adapter for PEMFC anode gas between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
        extends BaseClasses.PartialPhase;

        replaceable package Medium = Media.AnodeGas constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model (Modelica)"
          annotation (choicesAllMatching=true, Dialog(group=
                "Material properties"));

        Species.FluidNonionic H2(redeclare package Medium =
              Modelica.Media.IdealGases.SingleGases.H2 (referenceChoice=
                  Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                excludeEnthalpyOfFormation=false), redeclare package Data =
              FCSys.Characteristics.H2.Gas)
          annotation (Placement(transformation(extent={{-10,10},{10,30}})));
        Species.FluidNonionic H2O(redeclare package Data =
              FCSys.Characteristics.H2O.Gas (referenceChoice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                excludeEnthalpyOfFormation=false), redeclare final package
            Medium = Modelica.Media.IdealGases.SingleGases.H2O)
          annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));
        Junctions.Junction2 junction
          annotation (Placement(transformation(extent={{60,-50},{40,-30}})));

        Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final package
            Medium = Medium) "Modelica fluid port" annotation (Placement(
              transformation(extent={{70,-50},{90,-30}}), iconTransformation(
                extent={{70,-50},{90,-30}})));

      equation
        // H2
        connect(H2.face.normal, face.H2.normal) annotation (Line(
            points={{-8,20},{-40,20},{-40,5.55112e-16},{-80,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(H2.face.thermal, face.H2.thermal) annotation (Line(
            points={{-8,20},{-40,20},{-40,5.55112e-16},{-80,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(H2.heatPort, heatPort) annotation (Line(
            points={{8,20},{60,20},{60,5.55112e-16},{80,5.55112e-16}},
            color={191,0,0},
            smooth=Smooth.None));
        connect(H2.fluidPort, junction.purePort1) annotation (Line(
            points={{8,16},{34,16},{34,-36.2},{42,-36.2}},
            color={0,127,255},
            smooth=Smooth.None));

        // H2O
        connect(H2O.face.normal, face.H2O.normal) annotation (Line(
            points={{-8,-20},{-40,-20},{-40,5.55112e-16},{-80,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(H2O.face.thermal, face.H2O.thermal) annotation (Line(
            points={{-8,-20},{-40,-20},{-40,5.55112e-16},{-80,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(H2O.heatPort, heatPort) annotation (Line(
            points={{8,-20},{60,-20},{60,5.55112e-16},{80,5.55112e-16}},
            color={191,0,0},
            smooth=Smooth.None));
        connect(H2O.fluidPort, junction.purePort2) annotation (Line(
            points={{8,-24},{26,-24},{26,-44},{42,-44}},
            color={0,127,255},
            smooth=Smooth.None));

        // Mixture
        connect(junction.mixturePort, fluidPort) annotation (Line(
            points={{58,-40},{80,-40}},
            color={0,127,255},
            smooth=Smooth.None));
        annotation (Icon(graphics={Line(
                      points={{0,20},{0,-60}},
                      color={0,0,0},
                      smooth=Smooth.None,
                      pattern=LinePattern.Dash,
                      thickness=0.5),Line(
                      points={{0,0},{-80,0}},
                      color={127,127,127},
                      smooth=Smooth.None,
                      thickness=0.5),Line(
                      points={{0,0},{80,0}},
                      color={191,0,0},
                      smooth=Smooth.None),Line(
                      points={{0,-40},{80,-40}},
                      color={0,127,255},
                      smooth=Smooth.None)}));
      end AnodeGas;

      model CathodeGas
        "<html>Adapter for PEMFC cathode gas between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
        extends BaseClasses.PartialPhase;

        replaceable package Medium = Media.CathodeGas constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model (Modelica)"
          annotation (choicesAllMatching=true, Dialog(group=
                "Material properties"));

        Junctions.Junction3 junction(
          redeclare package Medium1 = Modelica.Media.IdealGases.SingleGases.H2O
              (referenceChoice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                excludeEnthalpyOfFormation=false),
          redeclare package Medium2 = Modelica.Media.IdealGases.SingleGases.N2
              (referenceChoice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                excludeEnthalpyOfFormation=false),
          redeclare package Medium3 = Modelica.Media.IdealGases.SingleGases.O2
              (referenceChoice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                excludeEnthalpyOfFormation=false),
          redeclare package MixtureMedium = Medium)
          annotation (Placement(transformation(extent={{60,-50},{40,-30}})));

        Species.FluidNonionic H2O(redeclare package Data =
              FCSys.Characteristics.H2O.Gas, redeclare final package Medium =
              Modelica.Media.IdealGases.SingleGases.H2O (referenceChoice=
                  Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                excludeEnthalpyOfFormation=false))
          annotation (Placement(transformation(extent={{-10,10},{10,30}})));
        Species.FluidNonionic N2(redeclare package Data =
              FCSys.Characteristics.N2.Gas, redeclare final package Medium =
              Modelica.Media.IdealGases.SingleGases.N2 (referenceChoice=
                  Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                excludeEnthalpyOfFormation=false))
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
        Species.FluidNonionic O2(redeclare package Data =
              FCSys.Characteristics.O2.Gas, redeclare final package Medium =
              Modelica.Media.IdealGases.SingleGases.O2 (referenceChoice=
                  Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
                excludeEnthalpyOfFormation=false))
          annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));

        Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final package
            Medium = Medium) "Modelica fluid port" annotation (Placement(
              transformation(extent={{70,-50},{90,-30}}), iconTransformation(
                extent={{70,-50},{90,-30}})));
      equation

        // H2O
        connect(H2O.face.normal, face.H2O.normal) annotation (Line(
            points={{-8,20},{-40,20},{-40,5.55112e-16},{-80,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(H2O.face.thermal, face.H2O.thermal) annotation (Line(
            points={{-8,20},{-40,20},{-40,5.55112e-16},{-80,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(H2O.fluidPort, junction.purePort1) annotation (Line(
            points={{8,16},{34,16},{34,-36.2},{42,-36.2}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(H2O.heatPort, heatPort) annotation (Line(
            points={{8,20},{60,20},{60,5.55112e-16},{80,5.55112e-16}},
            color={191,0,0},
            smooth=Smooth.None));

        // N2
        connect(N2.face.normal, face.N2.normal) annotation (Line(
            points={{-8,6.10623e-16},{-40,6.10623e-16},{-40,5.55112e-16},{-80,
                5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));

        connect(N2.face.thermal, face.N2.thermal) annotation (Line(
            points={{-8,6.10623e-16},{-40,6.10623e-16},{-40,5.55112e-16},{-80,
                5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));

        connect(N2.fluidPort, junction.purePort2) annotation (Line(
            points={{8,-4},{30,-4},{30,-40},{42,-40}},
            color={0,127,255},
            smooth=Smooth.None));

        connect(N2.heatPort, heatPort) annotation (Line(
            points={{8,6.10623e-16},{60,6.10623e-16},{60,5.55112e-16},{80,
                5.55112e-16}},
            color={191,0,0},
            smooth=Smooth.None));

        // O2
        connect(O2.face.normal, face.O2.normal) annotation (Line(
            points={{-8,-20},{-40,-20},{-40,5.55112e-16},{-80,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(O2.face.thermal, face.O2.thermal) annotation (Line(
            points={{-8,-20},{-40,-20},{-40,5.55112e-16},{-80,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(O2.fluidPort, junction.purePort3) annotation (Line(
            points={{8,-24},{26,-24},{26,-44},{42,-44}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(O2.heatPort, heatPort) annotation (Line(
            points={{8,-20},{60,-20},{60,5.55112e-16},{80,5.55112e-16}},
            color={191,0,0},
            smooth=Smooth.None));

        // Mixture
        connect(junction.mixturePort, fluidPort) annotation (Line(
            points={{58,-40},{80,-40}},
            color={0,127,255},
            smooth=Smooth.None));
        annotation (Icon(graphics={Line(
                      points={{0,20},{0,-60}},
                      color={0,0,0},
                      smooth=Smooth.None,
                      pattern=LinePattern.Dash,
                      thickness=0.5),Line(
                      points={{0,0},{-80,0}},
                      color={127,127,127},
                      smooth=Smooth.None,
                      thickness=0.5),Line(
                      points={{0,0},{80,0}},
                      color={191,0,0},
                      smooth=Smooth.None),Line(
                      points={{0,-40},{80,-40}},
                      color={0,127,255},
                      smooth=Smooth.None)}));
      end CathodeGas;

      model Graphite
        "<html>Adapter for graphite between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
        extends BaseClasses.PartialPhase;

        Species.'e-' 'e-'(redeclare package Data =
              FCSys.Characteristics.'e-'.Graphite)
          annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));

        Species.Solid 'C+'(redeclare package Data =
              FCSys.Characteristics.'C+'.Graphite)
          annotation (Placement(transformation(extent={{-10,10},{10,30}})));
        Modelica.Electrical.Analog.Interfaces.NegativePin pin
          "Modelica electrical pin" annotation (Placement(transformation(extent
                ={{70,30},{90,50}}), iconTransformation(extent={{70,30},{90,50}})));
      equation
        // C
        connect('C+'.face.thermal, face.C.thermal) annotation (Line(
            points={{-8,20},{-40,20},{-40,5.55112e-16},{-80,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect('C+'.heatPort, heatPort) annotation (Line(
            points={{8,20},{40,20},{40,5.55112e-16},{80,5.55112e-16}},
            color={191,0,0},
            smooth=Smooth.None));

        // e-
        connect('e-'.face.normal, face.'e-'.normal) annotation (Line(
            points={{-8,-20},{-40,-20},{-40,5.55112e-16},{-80,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));

        connect('e-'.face.thermal, face.'e-'.thermal) annotation (Line(
            points={{-8,-20},{-40,-20},{-40,5.55112e-16},{-80,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect('e-'.heatPort, heatPort) annotation (Line(
            points={{8,-20},{40,-20},{40,5.55112e-16},{80,5.55112e-16}},
            color={191,0,0},
            smooth=Smooth.None));

        connect('e-'.pin, pin) annotation (Line(
            points={{8,-16},{60,-16},{60,40},{80,40}},
            color={0,0,255},
            smooth=Smooth.None));

        annotation (Icon(graphics={Line(
                      points={{0,60},{0,-20}},
                      color={0,0,0},
                      smooth=Smooth.None,
                      pattern=LinePattern.Dash,
                      thickness=0.5),Line(
                      points={{0,0},{-80,0}},
                      color={127,127,127},
                      smooth=Smooth.None,
                      thickness=0.5),Line(
                      points={{0,40},{80,40}},
                      color={0,0,255},
                      smooth=Smooth.None),Line(
                      points={{0,0},{80,0}},
                      color={191,0,0},
                      smooth=Smooth.None)}));
      end Graphite;

      model Liquid
        "<html>Adapter for liquid between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
        extends BaseClasses.PartialPhase;

        replaceable package Medium =
            Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby
          Modelica.Media.Interfaces.PartialPureSubstance
          "Medium model (Modelica)" annotation (choicesAllMatching=true, Dialog(
              group="Material properties"));

        Species.FluidNonionic H2O(redeclare package Data =
              FCSys.Characteristics.H2O.Liquid, redeclare final package Medium
            = Medium)
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

        Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final package
            Medium = Medium) "Modelica fluid port" annotation (Placement(
              transformation(extent={{70,-50},{90,-30}}), iconTransformation(
                extent={{70,-50},{90,-30}})));

      equation
        // H2O
        connect(H2O.face.normal, face.H2.normal) annotation (Line(
            points={{-8,6.10623e-16},{-54,6.10623e-16},{-54,5.55112e-16},{-80,
                5.55112e-16}},
            color={0,0,0},
            smooth=Smooth.None));

        connect(H2O.face.thermal, face.H2.thermal) annotation (Line(
            points={{-8,6.10623e-16},{-54,6.10623e-16},{-54,5.55112e-16},{-80,
                5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));

        connect(H2O.heatPort, heatPort) annotation (Line(
            points={{8,6.10623e-16},{40,6.10623e-16},{40,5.55112e-16},{80,
                5.55112e-16}},
            color={191,0,0},
            smooth=Smooth.None));

        connect(H2O.fluidPort, fluidPort) annotation (Line(
            points={{8,-4},{40,-4},{40,-40},{80,-40}},
            color={0,127,255},
            smooth=Smooth.None));

        annotation (Icon(graphics={Line(
                      points={{0,20},{0,-60}},
                      color={0,0,0},
                      smooth=Smooth.None,
                      pattern=LinePattern.Dash,
                      thickness=0.5),Line(
                      points={{0,0},{-80,0}},
                      color={127,127,127},
                      smooth=Smooth.None,
                      thickness=0.5),Line(
                      points={{0,-40},{80,-40}},
                      color={0,127,255},
                      smooth=Smooth.None)}));
      end Liquid;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;

        partial model PartialPhase
          "<html>Partial adapter for a phase between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
          extends FCSys.BaseClasses.Icons.Names.Top3;

          Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
            "Modelica heat port" annotation (Placement(transformation(extent={{
                    70,-10},{90,10}}), iconTransformation(extent={{70,-10},{90,
                    10}})));
          FCSys.Connectors.FaceBus face "FCSys face connector" annotation (
              Placement(transformation(extent={{-90,-10},{-70,10}}),
                iconTransformation(extent={{-90,-10},{-70,10}})));

          annotation (Icon(graphics={Line(
                          points={{0,0},{-80,0}},
                          color={127,127,127},
                          smooth=Smooth.None,
                          thickness=0.5),Line(
                          points={{0,0},{80,0}},
                          color={191,0,0},
                          smooth=Smooth.None)}));
        end PartialPhase;
      end BaseClasses;
    end Phases;

    package Species "Adapters for single species"
      extends Modelica.Icons.Package;

      model 'e-'
        "<html>Adapter to connect e<sup>-</sup> between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a> (electrical and heat only)</html>"

        extends BaseClasses.PartialSpecies(redeclare
            FCSys.Characteristics.'e-'.Graphite Data, face(final iosbaric=false));

        Modelica.Electrical.Analog.Interfaces.NegativePin pin
          "Modelica electrical pin" annotation (Placement(transformation(extent
                ={{70,30},{90,50}}), iconTransformation(extent={{70,30},{90,50}})));

      equation
        // Efforts
        Data.g(face.thermal.T, Data.p_Tv(face.thermal.T, 1/face.normal.rho)) =
          Data.z*pin.v*U.V "Electrical potential";

        // Conservation (no storage)
        0 = face.normal.Ndot + pin.i*U.A/Data.z "Material";

        annotation (
          Documentation(info="<html><p>For additional information, see the
    <a href=\"modelica://FCSys.BCs.Adapters.Species.BaseClasses.PartialSpecies\">
    PartialSpecies</a> model.</p>
    </html>"),
          Icon(graphics={Line(
                      points={{0,40},{80,40}},
                      color={0,0,255},
                      smooth=Smooth.None),Line(
                      points={{0,60},{0,-20}},
                      color={0,0,0},
                      smooth=Smooth.None,
                      pattern=LinePattern.Dash)}),
          Diagram(graphics));
      end 'e-';

      model FluidNonionic
        "<html>Adapter to connect a single nonionic fluid species between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"

        extends BaseClasses.PartialSpecies(face(final isobaric=false));

        replaceable package Medium = Modelica.Media.IdealGases.SingleGases.H2O
          constrainedby Modelica.Media.Interfaces.PartialPureSubstance
          "Medium model (Modelica)" annotation (choicesAllMatching=true, Dialog(
              group="Material properties"));

        Medium.BaseProperties medium "Base properties of the fluid";

        Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final package
            Medium = Medium) "Modelica fluid port" annotation (Placement(
              transformation(extent={{70,-50},{90,-30}}), iconTransformation(
                extent={{70,-50},{90,-30}})));

      equation
        // Thermodynamic state and properties
        medium.p = fluidPort.p;
        medium.T = heatPort.T;
        medium.Xi = ones(Medium.nXi)/Medium.nXi;

        // Efforts
        face.normal.rho = (medium.d/medium.MM)*U.mol/U.m^3;
        medium.h = fluidPort.h_outflow;

        // Conservation (no storage)
        0 = face.normal.Ndot + (fluidPort.m_flow/medium.MM)*U.mol/U.s
          "Material";

        // See the partial model for additional equations.
        annotation (
          Documentation(info="<html><p>The electrical connector (<code>pin</code>) is only included
    if the species is ionic.
    </p>
    <p>For additional information, see the
    <a href=\"modelica://FCSys.BCs.Adapters.Species.BaseClasses.PartialSpecies\">
    PartialSpecies</a> model.</p>
    </html>"),
          Icon(graphics={Line(
                      points={{0,-40},{80,-40}},
                      color={0,127,255},
                      smooth=Smooth.None),Line(
                      points={{0,20},{0,-60}},
                      color={0,0,0},
                      smooth=Smooth.None,
                      pattern=LinePattern.Dash)}),
          Diagram(graphics));
      end FluidNonionic;

      model Solid
        "<html>Adapter to connect a single solid species between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a> (heat only)</html>"

        extends BaseClasses.PartialSpecies(face(final isobaric=true));

        annotation (Documentation(info="<html><p>For additional information, see the
    <a href=\"modelica://FCSys.BCs.Adapters.Species.BaseClasses.PartialSpecies\">
    PartialSpecies</a> model.</p>
    </html>"), Icon(graphics={Line(
                      points={{0,20},{0,-20}},
                      color={0,0,0},
                      smooth=Smooth.None,
                      pattern=LinePattern.Dash)}));
      end Solid;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;

        partial model PartialSpecies
          "<html>Partial single-species adapter between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
          extends FCSys.BaseClasses.Icons.Names.Top3;

          replaceable package Data =
              FCSys.Characteristics.BaseClasses.Characteristic
            "Characteristic data (FCSys)" annotation (
            Dialog(group="Material properties"),
            __Dymola_choicesAllMatching=true,
            Placement(transformation(extent={{-60,40},{-40,60}}),
                iconTransformation(extent={{-10,90},{10,110}})));

          FCSys.Connectors.Face face(
            axis=Axis.x,
            final inviscidX=true,
            final inviscidY=true,
            final inviscidZ=true)
            "Connector for material and heat of a single species (linear momentum excluded)"
            annotation (Placement(transformation(extent={{-90,-10},{-70,10}}),
                iconTransformation(extent={{-90,-10},{-70,10}})));

          // Note:  The axis doesn't matter since transverse linear momentum
          // isn't included.

          Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
            "Modelica heat port" annotation (Placement(transformation(extent={{
                    70,-10},{90,10}}), iconTransformation(extent={{70,-10},{90,
                    10}})));

        equation
          // Efforts
          face.thermal.T = heatPort.T*U.K "Temperature";

          // Conservation (no storage)
          0 = face.thermal.Qdot + heatPort.Q_flow*U.W "Energy";
          // Note:  The enthalpy, kinetic energy, and electric work terms each
          // cancel since there's no material storage and the thermodynamic state
          // and electrical potential is continuous across the junction.

          annotation (
            defaultComponentName="species",
            Documentation(info="<html><p>Note that shear force is not included.</p>
  </html>"),
            Icon(graphics={Line(
                          points={{0,0},{-80,0}},
                          color={127,127,127},
                          smooth=Smooth.None),Line(
                          points={{0,0},{80,0}},
                          color={191,0,0},
                          smooth=Smooth.None)}));
        end PartialSpecies;
      end BaseClasses;
    end Species;

    package Junctions
      "<html><a href=\"modelica://Modelica\">Modelica</a> junctions between pure substances and their mixtures</html>"
      extends Modelica.Icons.Package;

      model Junction2 "Junction between two pure substances and their mixture"
        extends BaseClasses.PartialJunction;

        replaceable package Medium1 = Modelica.Media.IdealGases.SingleGases.H2
            (referenceChoice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
              excludeEnthalpyOfFormation=false) constrainedby
          Modelica.Media.Interfaces.PartialPureSubstance
          "<html>Medium model for the 1<sup>st</sup> pure substance</html>"
          annotation (choicesAllMatching=true, Dialog(group=
                "Material properties"));
        replaceable package Medium2 = Modelica.Media.IdealGases.SingleGases.H2O
            (referenceChoice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
              excludeEnthalpyOfFormation=false) constrainedby
          Modelica.Media.Interfaces.PartialPureSubstance
          "<html>Medium model for the 2<sup>nd</sup> pure substance</html>"
          annotation (choicesAllMatching=true, Dialog(group=
                "Material properties"));

        Modelica.Fluid.Interfaces.FluidPort_b purePort1(redeclare final package
            Medium = Medium1) "Fluid port for the 1st pure substance"
          annotation (Placement(transformation(extent={{70,28},{90,48}}),
              iconTransformation(extent={{70,28},{90,48}})));
        Modelica.Fluid.Interfaces.FluidPort_b purePort2(redeclare final package
            Medium = Medium2) "Fluid port for the 2nd pure substance"
          annotation (Placement(transformation(extent={{70,-50},{90,-30}}),
              iconTransformation(extent={{70,-50},{90,-30}})));

      initial equation

        // Check the number and names of substances
        assert(MixtureMedium.nS == 2,
          "The mixture medium must have exactly two substances.");
        assert(MixtureMedium.substanceNames[1] == Medium1.substanceNames[1],
          "The first substance of the mixture medium (MixtureMedium) is \"" +
          MixtureMedium.substanceNames[1] + "\",
but the first pure substance is \"" + Medium1.substanceNames[1] + "\".");
        assert(MixtureMedium.substanceNames[2] == Medium2.substanceNames[1],
          "The second substance of the mixture medium (MixtureMedium) is \"" +
          MixtureMedium.substanceNames[2] + "\",
but the second pure substance is \"" + Medium2.substanceNames[1] + "\".");

        // Check the extra properties.
        assert(MixtureMedium.nC == Medium1.nC and MixtureMedium.nC == Medium2.nC,
          "The media must all have the same number of extra properties.");
        for i in 1:MixtureMedium.nC loop
          assert(MixtureMedium.extraPropertiesNames[i] == Medium1.extraPropertiesNames[
            i], "Extra property #" + String(i) +
            " of the mixture medium (MixtureMedium) is \"" + MixtureMedium.extraPropertiesNames[
            i] + "\",
but that of the first pure substance (Medium1) is \"" + Medium1.extraPropertiesNames[
            i] + "\".");
          assert(MixtureMedium.extraPropertiesNames[i] == Medium2.extraPropertiesNames[
            i], "Extra property #" + String(i) +
            " of the mixture medium (MixtureMedium) is \"" + MixtureMedium.extraPropertiesNames[
            i] + "\",
but that of the second pure substance (Medium2) is \"" + Medium2.extraPropertiesNames[
            i] + "\".");
        end for;

      equation
        // Dalton's law (additivity of pressure)
        mixturePort.p = purePort1.p + purePort2.p;

        // Streams
        // -------
        // Enthalpy
        purePort1.h_outflow = inStream(mixturePort.h_outflow);
        purePort2.h_outflow = inStream(mixturePort.h_outflow);
        mixturePort.h_outflow = X*{inStream(purePort1.h_outflow),inStream(
          purePort2.h_outflow)};
        //
        // Extra properties
        purePort1.C_outflow = inStream(mixturePort.C_outflow);
        purePort2.C_outflow = inStream(mixturePort.C_outflow);
        mixturePort.C_outflow = X*{inStream(purePort1.C_outflow),inStream(
          purePort2.C_outflow)};

        // Mass conservation (no storage)
        0 = X[1]*mixturePort.m_flow + purePort1.m_flow "Substance 1";
        0 = X[2]*mixturePort.m_flow + purePort2.m_flow "Substance 2";

        annotation (
          defaultComponentName="junction",
          Documentation(info="<html><p>
  Assumptions:
  <ol>
  <li>The mixing is ideal.  If the pure substances are being combined,
  then the massic enthalpy of the mixture is the mass-weighted sum of the pure substances.
  If the mixture is being split, then each of the pure substances receives fluid
  at the same massic enthalpy.
  </li>
  <li>The mixture observes Dalton's law.  The pressure of the mixture is the sum
  of the pressures of the pure substances.
  </li>
  </ol>
  </p></html>"),
          Icon(graphics={Line(
                      points={{0,-40},{80,-40}},
                      color={0,127,255},
                      smooth=Smooth.None),Line(
                      points={{0,40},{80,40}},
                      color={0,127,255},
                      smooth=Smooth.None)}));
      end Junction2;

      model Junction3
        "Junction between three pure substances and their mixture"
        extends BaseClasses.PartialJunction(redeclare replaceable package
            MixtureMedium = Media.CathodeGas);

        replaceable package Medium1 = Modelica.Media.IdealGases.SingleGases.H2O
            (referenceChoice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
              excludeEnthalpyOfFormation=false) constrainedby
          Modelica.Media.Interfaces.PartialPureSubstance
          "<html>Medium model for the 1<sup>st</sup> pure substance</html>"
          annotation (choicesAllMatching=true, Dialog(group=
                "Material properties"));
        replaceable package Medium2 = Modelica.Media.IdealGases.SingleGases.N2
            (referenceChoice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
              excludeEnthalpyOfFormation=false) constrainedby
          Modelica.Media.Interfaces.PartialPureSubstance
          "<html>Medium model for the 2<sup>nd</sup> pure substance</html>"
          annotation (choicesAllMatching=true, Dialog(group=
                "Material properties"));
        replaceable package Medium3 = Modelica.Media.IdealGases.SingleGases.O2
            (referenceChoice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,
              excludeEnthalpyOfFormation=false) constrainedby
          Modelica.Media.Interfaces.PartialPureSubstance
          "<html>Medium model for the 3<sup>rd</sup> pure substance</html>"
          annotation (choicesAllMatching=true, Dialog(group=
                "Material properties"));

        Modelica.Fluid.Interfaces.FluidPort_b purePort1(redeclare final package
            Medium = Medium1) "Fluid port for the 1st pure substance"
          annotation (Placement(transformation(extent={{70,28},{90,48}}),
              iconTransformation(extent={{70,28},{90,48}})));
        Modelica.Fluid.Interfaces.FluidPort_b purePort2(redeclare final package
            Medium = Medium2) "Fluid port for the 2nd pure substance"
          annotation (Placement(transformation(extent={{70,-10},{90,10}}),
              iconTransformation(extent={{70,-10},{90,10}})));
        Modelica.Fluid.Interfaces.FluidPort_b purePort3(redeclare final package
            Medium = Medium3) "Fluid port for the 3rd pure substance"
          annotation (Placement(transformation(extent={{70,-50},{90,-30}}),
              iconTransformation(extent={{70,-50},{90,-30}})));

      initial equation

        // Check the number and names of substances
        assert(MixtureMedium.nS == 3,
          "The mixture medium must have exactly three substances.");
        assert(MixtureMedium.substanceNames[1] == Medium1.substanceNames[1],
          "The first substance of the mixture medium (MixtureMedium) is \"" +
          MixtureMedium.substanceNames[1] + "\",
but the first pure substance is \"" + Medium1.substanceNames[1] + "\".");
        assert(MixtureMedium.substanceNames[2] == Medium2.substanceNames[1],
          "The second substance of the mixture medium (MixtureMedium) is \"" +
          MixtureMedium.substanceNames[2] + "\",
but the second pure substance is \"" + Medium2.substanceNames[1] + "\".");
        assert(MixtureMedium.substanceNames[3] == Medium3.substanceNames[1],
          "The second substance of the mixture medium (MixtureMedium) is \"" +
          MixtureMedium.substanceNames[3] + "\",
but the third pure substance is \"" + Medium2.substanceNames[1] + "\".");

        // Check the extra properties.
        assert(MixtureMedium.nC == Medium1.nC and MixtureMedium.nC == Medium2.nC
           and MixtureMedium.nC == Medium3.nC,
          "The media must all have the same number of extra properties.");
        for i in 1:MixtureMedium.nC loop
          assert(MixtureMedium.extraPropertiesNames[i] == Medium1.extraPropertiesNames[
            i], "Extra property #" + String(i) +
            " of the mixture medium (MixtureMedium) is \"" + MixtureMedium.extraPropertiesNames[
            i] + "\",
but that of the first pure substance (Medium1) is \"" + Medium1.extraPropertiesNames[
            i] + "\".");
          assert(MixtureMedium.extraPropertiesNames[i] == Medium2.extraPropertiesNames[
            i], "Extra property #" + String(i) +
            " of the mixture medium (MixtureMedium) is \"" + MixtureMedium.extraPropertiesNames[
            i] + "\",
but that of the second pure substance (Medium2) is \"" + Medium2.extraPropertiesNames[
            i] + "\".");
          assert(MixtureMedium.extraPropertiesNames[i] == Medium3.extraPropertiesNames[
            i], "Extra property #" + String(i) +
            " of the mixture medium (MixtureMedium) is \"" + MixtureMedium.extraPropertiesNames[
            i] + "\",
but that of the third pure substance (Medium3) is \"" + Medium3.extraPropertiesNames[
            i] + "\".");
        end for;

      equation
        // Dalton's law (additivity of pressure)
        mixturePort.p = purePort1.p + purePort2.p + purePort3.p;

        // Streams
        // -------
        // Enthalpy
        purePort1.h_outflow = inStream(mixturePort.h_outflow);
        purePort2.h_outflow = inStream(mixturePort.h_outflow);
        purePort3.h_outflow = inStream(mixturePort.h_outflow);
        mixturePort.h_outflow = X*{inStream(purePort1.h_outflow),inStream(
          purePort2.h_outflow),inStream(purePort3.h_outflow)};
        //
        // Extra properties
        purePort1.C_outflow = inStream(mixturePort.C_outflow);
        purePort2.C_outflow = inStream(mixturePort.C_outflow);
        purePort3.C_outflow = inStream(mixturePort.C_outflow);
        mixturePort.C_outflow = X*{inStream(purePort1.C_outflow),inStream(
          purePort2.C_outflow),inStream(purePort3.C_outflow)};

        // Mass conservation (no storage)
        0 = X[1]*mixturePort.m_flow + purePort1.m_flow "Substance 1";
        0 = X[2]*mixturePort.m_flow + purePort2.m_flow "Substance 2";
        0 = X[3]*mixturePort.m_flow + purePort3.m_flow "Substance 3";

        annotation (
          defaultComponentName="junction",
          Documentation(info="<html><p>
  Assumptions:
  <ol>
  <li>The mixing is ideal.  If the pure substances are being combined,
  then the massic enthalpy of the mixture is the mass-weighted sum of the pure substances.
  If the mixture is being split, then each of the pure substances receives fluid
  at the same massic enthalpy.
  </li>
  <li>The mixture observes Dalton's law.  The pressure of the mixture is the sum
  of the pressures of the pure substances.
  </li>
  </ol>
  </p></html>"),
          Icon(graphics={Line(
                      points={{0,-40},{80,-40}},
                      color={0,127,255},
                      smooth=Smooth.None),Line(
                      points={{0,40},{80,40}},
                      color={0,127,255},
                      smooth=Smooth.None),Line(
                      points={{6,0},{80,0}},
                      color={0,127,255},
                      smooth=Smooth.None)}));
      end Junction3;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialJunction
          "Partial model for a junction between pure substances and their mixture"
          extends FCSys.BaseClasses.Icons.Names.Top3;

          replaceable package MixtureMedium = Media.AnodeGas constrainedby
            Modelica.Media.Interfaces.PartialMedium
            "Medium model for the mixture" annotation (choicesAllMatching=true,
              Dialog(group="Material properties"));

          Modelica.Fluid.Interfaces.FluidPort_a mixturePort(redeclare final
              package Medium = MixtureMedium) "Fluid port for the mixture"
            annotation (Placement(transformation(extent={{-90,-10},{-70,10}}),
                iconTransformation(extent={{-90,-10},{-70,10}})));

          SI.MassFraction X[MixtureMedium.nX]
            "Mass fractions within the mixture";

        equation
          // Mass fractions
          X = if MixtureMedium.fixedX then MixtureMedium.reference_X else if
            MixtureMedium.reducedX then cat(
                    1,
                    inStream(mixturePort.Xi_outflow),
                    1 - sum(X[1:MixtureMedium.nXi])) else inStream(mixturePort.Xi_outflow);
          X = if MixtureMedium.reducedX then cat(
                    1,
                    mixturePort.Xi_outflow,
                    1 - sum(X[1:MixtureMedium.nXi])) else mixturePort.Xi_outflow;

          annotation (defaultComponentName="junction", Icon(graphics={Line(
                          points={{-80,0},{0,0}},
                          color={0,127,255},
                          smooth=Smooth.None),Line(
                          points={{0,-40},{0,40}},
                          color={0,127,255},
                          smooth=Smooth.None),Ellipse(
                          extent={{-6,6},{6,-6}},
                          lineColor={0,127,255},
                          fillColor={255,255,255},
                          fillPattern=FillPattern.Solid)}));
        end PartialJunction;
      end BaseClasses;
    end Junctions;

    package Media
      "<html><a href=\"modelica://Modelica.Media\">Modelica media</a> models to interface with the <a href=\"modelica://FCSys.Assemblies.Cells.Cell\">fuel cell</a></html>"
      extends Modelica.Icons.MaterialPropertiesPackage;

      package AnodeGas "Gas mixture for PEMFC anode (H2 and H2O)"
        extends Modelica.Media.IdealGases.Common.MixtureGasNasa(
          mediumName="AnodeGas",
          data={Modelica.Media.IdealGases.Common.SingleGasesData.H2,Modelica.Media.IdealGases.Common.SingleGasesData.H2O},

          fluidConstants={Modelica.Media.IdealGases.Common.FluidData.H2,
              Modelica.Media.IdealGases.Common.FluidData.H2O},
          substanceNames={"H2","H2O"},
          reference_X=fill(1/nX, nX),
          referenceChoice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,

          excludeEnthalpyOfFormation=false);

        annotation (Documentation(info="<html>

</html>"));
      end AnodeGas;

      package CathodeGas "Gas mixture for PEMFC cathode (H2O, N2, and O2)"
        extends Modelica.Media.IdealGases.Common.MixtureGasNasa(
          mediumName="CathodeGas",
          data={Modelica.Media.IdealGases.Common.SingleGasesData.H2O,Modelica.Media.IdealGases.Common.SingleGasesData.N2,
              Modelica.Media.IdealGases.Common.SingleGasesData.O2},
          fluidConstants={Modelica.Media.IdealGases.Common.FluidData.H2O,
              Modelica.Media.IdealGases.Common.FluidData.N2,Modelica.Media.IdealGases.Common.FluidData.O2},

          substanceNames={"H2O","N2","O2"},
          reference_X=fill(1/nX, nX),
          referenceChoice=Modelica.Media.Interfaces.PartialMedium.Choices.ReferenceEnthalpy.ZeroAt25C,

          excludeEnthalpyOfFormation=false);

        annotation (Documentation(info="<html>

</html>"));
      end CathodeGas;
    end Media;
  end Adapters;

  package TestStands "Test stands"
    extends Modelica.Icons.Package;

    model TestProfile "Cell test profile"
      extends Modelica.Icons.Example;
      extends BaseClasses.PartialTestStandNoIO;

      annotation (structurallyIncomplete=true, Diagram(coordinateSystem(
              preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
            graphics));
    end TestProfile;

    model Replay
      "Regenerate signals recorded from HNEI's Greenlight FC test stand"
      extends FCSys.BaseClasses.Icons.Blocks.ContinuousShort;

      parameter Integer n(final min=1) = 1 "index of the data set";
      parameter Boolean terminateMaxTime=true
        "Terminate at maximum time of source data"
        annotation (Dialog(compact=true), choices(__Dymola_checkBox=true));

      Modelica.Blocks.Sources.CombiTimeTable combiTimeTable(
        tableOnFile=true,
        extrapolation=Modelica.Blocks.Types.Extrapolation.HoldLastPoint,
        tableName="data" + String(n),
        fileName="FCSys/test/LOOCV/data.mat",
        columns=2:19,
        smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments)
        "Block to load and replay data" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,80})));
      FCSys.Connectors.RealOutputBus y "Output signals as a bus" annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-100}),iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-40})));

    protected
      Modelica.Blocks.Math.Add sumAnMFC annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-14,30})));
      Modelica.Blocks.Math.Add sumCaMFC annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={14,30})));

      Modelica.Blocks.Math.Gain unit1(k=U.V) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-152,0})));
      Modelica.Blocks.Math.Gain unit2(k=U.'%') annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-134,-30})));
      Modelica.Blocks.Math.Gain unit3(k=U.'%') annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-114,0})));
      FCSys.BCs.BaseClasses.RealFunction unit4(y=unit4.u*U.kPa + U.atm)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-94,-30})));
      FCSys.BCs.BaseClasses.RealFunction unit5(y=unit5.u*U.kPa + U.atm)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-74,0})));
      FCSys.BCs.BaseClasses.RealFunction unit6(y=unit6.u*U.kPa + U.atm)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-54,-30})));
      FCSys.BCs.BaseClasses.RealFunction unit7(y=unit7.u*U.kPa + U.atm)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-34,0})));
      Modelica.Blocks.Math.Gain unit8(k=U.L/U.min) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-14,-30})));
      Modelica.Blocks.Math.Gain unit9(k=U.L/U.min) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={14,0})));
      FCSys.BCs.BaseClasses.RealFunction unit10(y=(unit10.u + 273.15)*U.K)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={32,-30})));
      FCSys.BCs.BaseClasses.RealFunction unit11(y=(unit11.u + 273.15)*U.K)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={52,0})));
      FCSys.BCs.BaseClasses.RealFunction unit12(y=(unit12.u + 273.15)*U.K)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={72,-30})));
      FCSys.BCs.BaseClasses.RealFunction unit13(y=(unit13.u + 273.15)*U.K)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={92,0})));
      FCSys.BCs.BaseClasses.RealFunction unit14(y=(unit14.u + 273.15)*U.K)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={112,-30})));
      FCSys.BCs.BaseClasses.RealFunction unit15(y=(unit15.u + 273.15)*U.K)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={132,0})));
      Modelica.Blocks.Math.Gain unit16(k=U.A) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={152,-30})));

      FCSys.Connectors.RealOutputInternal Deltamu(final unit="l2.m/(N.T2)")
        "CVM Cell 1 Voltage" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-152,-60})));
      FCSys.Connectors.RealOutputInternal RHAnFPNegX(
        final unit="1",
        displayUnit="%",
        final min=0) "Anode inlet RH" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-134,-60})));
      FCSys.Connectors.RealOutputInternal RHCaFPNegX(
        final unit="1",
        displayUnit="%",
        final min=0) "Cathode inlet RH" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-114,-60})));
      FCSys.Connectors.RealOutputInternal p_anFPNegY(final unit="m/(l.T2)")
        "Pressure anode inlet" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-94,-60})));
      FCSys.Connectors.RealOutputInternal p_anFPPosY(final unit="m/(l.T2)",
          final min=0) "Pressure anode outlet" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-74,-60})));
      FCSys.Connectors.RealOutputInternal p_caFPNegY(final unit="m/(l.T2)")
        "Pressure cathode inlet" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-54,-60})));
      FCSys.Connectors.RealOutputInternal p_caFPPosY(final unit="m/(l.T2)",
          final min=0) "Pressure anode outlet" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-34,-60})));
      FCSys.Connectors.RealOutputInternal Vdot_anFPNegY_H2(final unit="l3/T")
        "Flow anode H2 MFC" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-14,-60})));
      FCSys.Connectors.RealOutputInternal Vdot_caFPNegY_air(final unit="l3/T")
        "Flow cathode H2 MFC" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={14,-60})));
      FCSys.Connectors.RealOutputInternal T_anFPNegY(
        final unit="l2.m/(N.T2)",
        displayUnit="K",
        final min=0) "Temperature anode inlet" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={32,-60})));
      FCSys.Connectors.RealOutputInternal T_anFPPosY(
        final unit="l2.m/(N.T2)",
        displayUnit="K",
        final min=0) "Temperature anode outlet" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={52,-60})));
      FCSys.Connectors.RealOutputInternal T_caFPNegY(
        final unit="l2.m/(N.T2)",
        displayUnit="K",
        final min=0) "Temperature cathode inlet" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={72,-60})));
      FCSys.Connectors.RealOutputInternal T_caFPPosY(
        final unit="l2.m/(N.T2)",
        displayUnit="K",
        final min=0) "Temperature cathode outlet" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={92,-60})));
      FCSys.Connectors.RealOutputInternal T_anFPX(
        final unit="l2.m/(N.T2)",
        displayUnit="K",
        final min=0) "Temperature end plate anode" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={112,-60})));
      FCSys.Connectors.RealOutputInternal T_caFPX(
        final unit="l2.m/(N.T2)",
        displayUnit="K",
        final min=0) "Temperature end plate cathode" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={132,-60})));
      FCSys.Connectors.RealOutputInternal J(final unit="N/T") "Measured load"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={152,-60})));

    equation
      //  Terminate as desired
      if terminateMaxTime then
        when time > combiTimeTable.t_max then
          terminate("The end of the data has been reached.");
        end when;
      end if;

      // Connections from source to unit conversion
      connect(unit1.u, combiTimeTable.y[1]) annotation (Line(
          points={{-152,12},{-152,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(unit2.u, combiTimeTable.y[2]) annotation (Line(
          points={{-134,-18},{-134,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(unit3.u, combiTimeTable.y[3]) annotation (Line(
          points={{-114,12},{-114,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(unit4.u, combiTimeTable.y[4]) annotation (Line(
          points={{-94,-20},{-94,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(unit5.u, combiTimeTable.y[5]) annotation (Line(
          points={{-74,10},{-74,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(unit6.u, combiTimeTable.y[6]) annotation (Line(
          points={{-54,-20},{-54,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(unit7.u, combiTimeTable.y[7]) annotation (Line(
          points={{-34,10},{-34,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(sumAnMFC.u1, combiTimeTable.y[8]) annotation (Line(
          points={{-8,42},{-8,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(sumAnMFC.u2, combiTimeTable.y[9]) annotation (Line(
          points={{-20,42},{-20,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(sumCaMFC.u1, combiTimeTable.y[10]) annotation (Line(
          points={{20,42},{20,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(sumCaMFC.u2, combiTimeTable.y[11]) annotation (Line(
          points={{8,42},{8,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(unit10.u, combiTimeTable.y[12]) annotation (Line(
          points={{32,-20},{32,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(unit11.u, combiTimeTable.y[13]) annotation (Line(
          points={{52,10},{52,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(unit12.u, combiTimeTable.y[14]) annotation (Line(
          points={{72,-20},{72,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(unit13.u, combiTimeTable.y[15]) annotation (Line(
          points={{92,10},{92,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(unit14.u, combiTimeTable.y[16]) annotation (Line(
          points={{112,-20},{112,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(unit15.u, combiTimeTable.y[17]) annotation (Line(
          points={{132,10},{132,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(unit16.u, combiTimeTable.y[18]) annotation (Line(
          points={{152,-18},{152,50},{-1.44329e-15,50},{-1.44329e-15,69}},
          color={0,0,127},
          smooth=Smooth.None));

      // Connections from unit conversion to internal outputs
      connect(Deltamu, unit1.y) annotation (Line(
          points={{-152,-60},{-152,-11}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(RHAnFPNegX, unit2.y) annotation (Line(
          points={{-134,-60},{-134,-41}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(RHCaFPNegX, unit3.y) annotation (Line(
          points={{-114,-60},{-114,-11}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p_anFPNegY, unit4.y) annotation (Line(
          points={{-94,-60},{-94,-40}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p_anFPPosY, unit5.y) annotation (Line(
          points={{-74,-60},{-74,-10}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p_caFPNegY, unit6.y) annotation (Line(
          points={{-54,-60},{-54,-40}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p_caFPPosY, unit7.y) annotation (Line(
          points={{-34,-60},{-34,-10}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(Vdot_anFPNegY_H2, unit8.y) annotation (Line(
          points={{-14,-60},{-14,-50.5},{-14,-41},{-14,-41}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(Vdot_caFPNegY_air, unit9.y) annotation (Line(
          points={{14,-60},{14,-35.5},{14,-11},{14,-11}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_anFPNegY, unit10.y) annotation (Line(
          points={{32,-60},{32,-50},{32,-40},{32,-40}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_anFPPosY, unit11.y) annotation (Line(
          points={{52,-60},{52,-10}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_caFPNegY, unit12.y) annotation (Line(
          points={{72,-60},{72,-40}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_caFPPosY, unit13.y) annotation (Line(
          points={{92,-60},{92,-10}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_anFPX, unit14.y) annotation (Line(
          points={{112,-60},{112,-40}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_caFPX, unit15.y) annotation (Line(
          points={{132,-60},{132,-10}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(J, unit16.y) annotation (Line(
          points={{152,-60},{152,-41}},
          color={0,0,127},
          smooth=Smooth.None));

      // Summations
      connect(sumAnMFC.y, unit8.u) annotation (Line(
          points={{-14,19},{-14,9.75},{-14,9.75},{-14,0.5},{-14,-18},{-14,-18}},

          color={0,0,127},
          smooth=Smooth.None));

      connect(sumCaMFC.y, unit9.u) annotation (Line(
          points={{14,19},{14,17.25},{14,17.25},{14,15.5},{14,12},{14,12}},
          color={0,0,127},
          smooth=Smooth.None));

      // Connections from internal outputs to public output
      connect(Deltamu, y.Deltamu) annotation (Line(
          points={{-152,-60},{-152,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(RHAnFPNegX, y.RHAnFPNegX) annotation (Line(
          points={{-134,-60},{-134,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(RHCaFPNegX, y.RHCaFPNegX) annotation (Line(
          points={{-114,-60},{-114,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(p_anFPNegY, y.p_anFPNegY) annotation (Line(
          points={{-94,-60},{-94,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(p_anFPPosY, y.p_anFPPosY) annotation (Line(
          points={{-74,-60},{-74,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(p_caFPNegY, y.p_caFPNegY) annotation (Line(
          points={{-54,-60},{-54,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(p_caFPPosY, y.p_caFPPosY) annotation (Line(
          points={{-34,-60},{-34,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(Vdot_anFPNegY_H2, y.Vdot_anFPNegY_H2) annotation (Line(
          points={{-14,-60},{-14,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(Vdot_caFPNegY_air, y.Vdot_caFPNegY_air) annotation (Line(
          points={{14,-60},{14,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_anFPNegY, y.T_anFPNegY) annotation (Line(
          points={{32,-60},{32,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_anFPPosY, y.T_anFPPosY) annotation (Line(
          points={{52,-60},{52,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_caFPNegY, y.T_caFPNegY) annotation (Line(
          points={{72,-60},{72,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_caFPPosY, y.T_caFPPosY) annotation (Line(
          points={{92,-60},{92,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(T_anFPX, y.T_anFPX) annotation (Line(
          points={{112,-60},{112,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(T_caFPX, y.T_caFPX) annotation (Line(
          points={{132,-60},{132,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(J, y.J) annotation (Line(
          points={{152,-60},{152,-80},{5.55112e-16,-80},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      annotation (
        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-180,-100},
                {180,100}}),graphics),
        experiment(StopTime=15481, Algorithm="Euler"),
        experimentSetupOutput);
    end Replay;

    package BaseClasses "Base classes (not for direct use)"
      extends Modelica.Icons.BasesPackage;
      partial model PartialTestStand "Partial cell test stand"
        extends FCSys.BaseClasses.Icons.Names.Top9;
        final parameter Integer n_x_an=1
          "<html>Number of subregions along the through-cell axis in anode FP (<i>n</i><sub>x an</sub>)</html>";
        final parameter Integer n_x_ca=1
          "<html>Number of subregions along the through-cell axis in anode FP (<i>n</i><sub>x ca</sub>)</html>";
        final parameter Integer n_y=1
          "<html>Number of subregions along the channel (<i>n</i><sub>y</sub>)</html>";
        final parameter Integer n_z=1
          "<html>Number of subregions across the channel (<i>n</i><sub>z</sub>)</html>";

        parameter Boolean inclIO=false "true, if input and output are included"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true));

        FCSys.Connectors.FaceBus anEnd[n_y, n_z] "Anode end plate" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-160,0}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-160,0})));
        FCSys.Connectors.FaceBus caEnd[n_y, n_z] "Cathode end plate"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={160,0}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={162,0})));
        FCSys.Connectors.FaceBus anSource[n_x_an, n_z] "Anode source"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-40,-160}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-40,-160})));
        FCSys.Connectors.FaceBus anSink[n_x_an, n_z] "Anode sink" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-40,160}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-40,160})));
        FCSys.Connectors.FaceBus caSource[n_x_ca, n_z] "Cathode source"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={40,-160}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={40,160})));
        FCSys.Connectors.FaceBus caSink[n_x_ca, n_z] "Cathode sink" annotation
          (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={40,160}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={40,-160})));

        FaceBus.SubregionClosed anEndBC[n_y, n_z](each final axis=FCSys.BaseClasses.Axis.x,
            each graphite('inclC+'=true, 'incle-'=true)) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-136,0})));
        FaceBus.SubregionClosed caEndBC[n_y, n_z](each final axis=FCSys.BaseClasses.Axis.x,
            each graphite('inclC+'=true, 'incle-'=true)) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={136,0})));
        FaceBus.SubregionClosed anSourceBC[n_x_an, n_z](each final axis=FCSys.BaseClasses.Axis.y,
            each gas(inclH2=true, inclH2O=true)) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-40,-136})));
        FaceBus.SubregionClosed anSinkBC[n_x_an, n_z](each final axis=FCSys.BaseClasses.Axis.y,
            each gas(inclH2=true, inclH2O=true)) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-40,136})));
        FaceBus.SubregionClosed caSourceBC[n_x_ca, n_z](each final axis=FCSys.BaseClasses.Axis.y,
            each gas(
            inclH2O=true,
            inclN2=true,
            inclO2=true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={40,-136})));
        FaceBus.SubregionClosed caSinkBC[n_x_ca, n_z](each final axis=FCSys.BaseClasses.Axis.y,
            each gas(
            inclH2O=true,
            inclN2=true,
            inclO2=true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={40,136})));
        FCSys.Connectors.RealInputBus u[n_y, n_z] if inclIO annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=315,
              origin={-160,160})));
        FCSys.Connectors.RealOutputBus y[n_y, n_z] if inclIO annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=315,
              origin={160,-160})));
        replaceable FaceBusDifferential.Subregion current[n_y, n_z](each final
            axis=FCSys.BaseClasses.Axis.x, graphite(
            'inclC+'=true,
            'C+'(isobaric=true),
            'incle-'=true,
            'e-'(isobaric=false))) if inclIO constrainedby
          FaceBusDifferential.Subregion(graphite(
            'inclC+'=true,
            'C+'(isobaric=true),
            'incle-'=true,
            'e-'(isobaric=false)))
          annotation (Placement(transformation(extent={{-140,20},{-120,40}})));

        replaceable FCSys.Sensors.FaceBusDifferential.Subregion voltage[n_y,
          n_z](each final axis=FCSys.BaseClasses.Axis.x) if inclIO
          annotation (Placement(transformation(extent={{120,-40},{140,-20}})));
      equation

        connect(anSourceBC.face, anSource) annotation (Line(
            points={{-40,-140},{-40,-160}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));

        connect(anSinkBC.face, anSink) annotation (Line(
            points={{-40,140},{-40,160}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));

        connect(caSourceBC.face, caSource) annotation (Line(
            points={{40,-140},{40,-160}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));

        connect(caSinkBC.face, caSink) annotation (Line(
            points={{40,140},{40,160}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));

        connect(caEndBC.face, caEnd) annotation (Line(
            points={{140,3.65701e-16},{140,5.55112e-16},{160,5.55112e-16}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(current.positive, caEnd) annotation (Line(
            points={{-120,30},{150,30},{150,5.55112e-16},{160,5.55112e-16}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));

        connect(current.negative, anEnd) annotation (Line(
            points={{-140,30},{-150,30},{-150,5.55112e-16},{-160,5.55112e-16}},

            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));

        connect(u, current.u) annotation (Line(
            points={{-160,160},{-130,130},{-130,34}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(voltage.negative, anEndBC.face) annotation (Line(
            points={{120,-30},{-150,-30},{-150,1.23436e-15},{-140,1.23436e-15}},

            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));

        connect(caEnd, voltage.positive) annotation (Line(
            points={{160,5.55112e-16},{150,5.55112e-16},{150,-30},{140,-30}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(voltage.y, y) annotation (Line(
            points={{130,-40},{130,-140},{160,-160}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        connect(anEndBC.face, anEnd) annotation (Line(
            points={{-140,1.23436e-15},{-140,5.55112e-16},{-160,5.55112e-16}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (
          defaultComponentName="testStand",
          Diagram(coordinateSystem(preserveAspectRatio=true,extent={{-160,-160},
                  {160,160}}), graphics),
          Icon(coordinateSystem(preserveAspectRatio=true,extent={{-160,-160},{
                  160,160}}), graphics={Rectangle(
                      extent={{-160,160},{160,-160}},
                      lineColor={191,191,191},
                      fillColor={255,255,255},
                      fillPattern=FillPattern.Backward),Rectangle(extent={{-160,
                160},{160,-160}}, lineColor={0,0,0})}));
      end PartialTestStand;

      partial model PartialTestStandNoIO
        "Partial cell test stand without inputs/outputs"
        extends FCSys.BaseClasses.Icons.Names.Top9;

        final parameter Integer n_x_an=1
          "<html>Number of subregions along the through-cell axis in anode FP (<i>n</i><sub>x an</sub>)</html>";
        final parameter Integer n_x_ca=1
          "<html>Number of subregions along the through-cell axis in anode FP (<i>n</i><sub>x ca</sub>)</html>";
        final parameter Integer n_y=1
          "<html>Number of subregions along the channel (<i>n</i><sub>y</sub>)</html>";
        final parameter Integer n_z=1
          "<html>Number of subregions across the channel (<i>n</i><sub>z</sub>)</html>";

        FaceBus.SubregionClosed anEnd[n_y, n_z](each final axis=FCSys.BaseClasses.Axis.x,
            each graphite(
            'inclC+'=true,
            'incle-'=true,
            'e-'(redeclare Face.Material.Current normal(redeclare
                  Modelica.Blocks.Sources.Ramp spec(height=U.A,duration=50)))))
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-30,0})));
        FaceBus.SubregionClosed caEnd[n_y, n_z](each final axis=FCSys.BaseClasses.Axis.x,
            each graphite(
            'inclC+'=true,
            'incle-'=true,
            'e-'(redeclare Face.Material.Current normal(redeclare
                  Modelica.Blocks.Sources.Ramp spec(height=U.A,duration=50)))))
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={30,0})));
        FaceBus.SubregionClosed anSource[n_x_an, n_z](each final axis=FCSys.BaseClasses.Axis.y,
            each gas(inclH2=true, inclH2O=true)) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-20,-30})));
        FaceBus.SubregionClosed anSink[n_x_an, n_z](each final axis=FCSys.BaseClasses.Axis.y,
            each gas(inclH2=true, inclH2O=true)) annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-20,30})));
        FaceBus.SubregionClosed caSource[n_x_ca, n_z](each final axis=FCSys.BaseClasses.Axis.y,
            each gas(
            inclH2O=true,
            inclN2=true,
            inclO2=true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={20,-30})));
        FaceBus.SubregionClosed caSink[n_x_ca, n_z](each final axis=FCSys.BaseClasses.Axis.y,
            each gas(
            inclH2O=true,
            inclN2=true,
            inclO2=true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={20,30})));

        inner Environment environment
          annotation (Placement(transformation(extent={{50,20},{70,40}})));
      equation

        annotation (
          defaultComponentName="testStand",
          Diagram(coordinateSystem(preserveAspectRatio=true,extent={{-100,-100},
                  {100,100}}), graphics),
          Icon(coordinateSystem(preserveAspectRatio=true,extent={{-160,-160},{
                  160,160}}), graphics));
      end PartialTestStandNoIO;
    end BaseClasses;
  end TestStands;

  package ChemicalBus
    "<html>BCs for the <a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a> connector</html>"
    extends Modelica.Icons.Package;

    model Gas "BC for gas"

      extends BaseClasses.NullPhase;

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
      Chemical.Species H2(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ,
        redeclare FCSys.Characteristics.H2.Gas Data) if inclH2 "Model"
        annotation (Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclH2), Placement(transformation(extent={{-10,-10},{10,10}})));

      parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      Chemical.Species H2O(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ,
        redeclare FCSys.Characteristics.H2O.Gas Data) if inclH2O "Model"
        annotation (Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,10}})));

      parameter Boolean inclN2=false "<html>Nitrogen (N<sub>2</sub>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      Chemical.Species N2(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ,
        redeclare FCSys.Characteristics.N2.Gas Data) if inclN2 "Model"
        annotation (Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclN2), Placement(transformation(extent={{-10,-10},{10,10}})));

      parameter Boolean inclO2=false "<html>Oxygen (O<sub>2</sub>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      Chemical.Species O2(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ,
        redeclare FCSys.Characteristics.O2.Gas Data) if inclO2 "Model"
        annotation (Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclO2), Placement(transformation(extent={{-10,-10},{10,10}})));

    equation
      // Note:  It would be helpful if Modelica allowed elements of expandable
      // connectors to be named by the contents of a string variable and the
      // name of an instance of a model was accessible through a string (like
      // %name is expanded to be the name of the instance of the model).  Then,
      // the connection equations that follow could be generic.

      // H2
      connect(H2.chemical, chemical.H2) annotation (Line(
          points={{-5.08852e-16,-4},{5.55112e-16,-40}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(u.H2, H2.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,14},{-5.08852e-16,14},{-5.08852e-16,
              4}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      // H2O
      connect(H2O.chemical, chemical.H2O) annotation (Line(
          points={{-5.08852e-16,-4},{5.55112e-16,-40}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(u.H2O, H2O.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,14},{-5.08852e-16,14},{-5.08852e-16,
              4}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      // N2
      connect(N2.chemical, chemical.N2) annotation (Line(
          points={{-5.08852e-16,-4},{5.55112e-16,-40}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(u.N2, N2.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,14},{-5.08852e-16,14},{-5.08852e-16,
              4}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      // O2
      connect(O2.chemical, chemical.O2) annotation (Line(
          points={{-5.08852e-16,-4},{5.55112e-16,-40}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(u.O2, O2.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,14},{-5.08852e-16,14},{-5.08852e-16,
              4}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));
    end Gas;

    model Graphite "BC for graphite"

      extends BaseClasses.NullPhase;

      // Conditionally include species.
      parameter Boolean Error "<html>Carbon plus (C<sup>+</sup>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      Chemical.Species 'C+'(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ,
        redeclare FCSys.Characteristics.'C+'.Graphite Data) if '' + ERROR
        annotation (Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable='' + ERROR annotation (
            Evaluate=true,
            HideResult=true,
            choices(__Dymola_checkBox=true),
            Dialog(
              group="Species",
              __Dymola_descriptionLabel=true,
              __Dymola_joinNext=true))));
      Chemical.Species 'e-'(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ,
        redeclare FCSys.Characteristics.'e-'.Graphite Data) if 'incle-' "Model"
        annotation (Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable='incle-'), Placement(transformation(extent={{-10,-10},{10,10}})));

    equation
      // C+
      connect('C+'.chemical, chemical.'C+') annotation (Line(
          points={{-5.08852e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(u.'C+', 'C+'.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,14},{-5.08852e-16,14},{-5.08852e-16,
              4}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      // e-
      connect('e-'.chemical, chemical.'e-') annotation (Line(
          points={{-5.08852e-16,-4},{5.55112e-16,-40}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(u.'e-', 'e-'.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,14},{-5.08852e-16,14},{-5.08852e-16,
              4}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));
    end Graphite;

    model Ionomer "BC for ionomer"

      extends BaseClasses.NullPhase;

      // Conditionally include species.
      parameter Boolean 'inclC19HF37O5S-'=false
        "<html>Nafion sulfonate (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      Chemical.Species 'C19HF37O5S-'(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ,
        redeclare FCSys.Characteristics.'C19HF37O5S-'.Ionomer Data) if
        'inclC19HF37O5S-' "Model" annotation (Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable='inclC19HF37O5S-'), Placement(transformation(extent={{-10,-10},
                {10,10}})));
      parameter Boolean 'inclH+'=false "<html>Protons (H<sup>+</sup>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      Chemical.Species 'H+'(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ,
        redeclare FCSys.Characteristics.'H+'.Ionomer Data) if 'inclH+' "Model"
        annotation (Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable='inclH+'), Placement(transformation(extent={{-10,-10},{10,10}})));
      parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      Chemical.Species H2O(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ,
        redeclare FCSys.Characteristics.H2O.Gas Data) if inclH2O "Model"
        annotation (Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,10}})));

    equation
      // C19HF37O5S-
      connect('C19HF37O5S-'.chemical, chemical.'C19HF37O5S-') annotation (Line(
          points={{-5.08852e-16,-4},{5.55112e-16,-40}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(u.'C19HF37O5S-', 'C19HF37O5S-'.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,14},{-5.08852e-16,14},{-5.08852e-16,
              4}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      // H+
      connect('H+'.chemical, chemical.'H+') annotation (Line(
          points={{-5.08852e-16,-4},{5.55112e-16,-40}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(u.'H+', 'H+'.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,14},{-5.08852e-16,14},{-5.08852e-16,
              4}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      // H2O
      connect(H2O.chemical, chemical.H2O) annotation (Line(
          points={{-5.08852e-16,-4},{5.55112e-16,-40}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(u.H2O, H2O.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,14},{-5.08852e-16,14},{-5.08852e-16,
              4}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));
    end Ionomer;

    model Liquid "BC for liquid"

      extends BaseClasses.NullPhase;

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
      Chemical.Species H2O(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ,
        redeclare FCSys.Characteristics.H2O.Liquid Data) if inclH2O "Model"
        annotation (Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,10}})));

    equation
      // H2O
      connect(H2O.chemical, chemical.H2O) annotation (Line(
          points={{-5.08852e-16,-4},{5.55112e-16,-40}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(u.H2O, H2O.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,14},{-5.08852e-16,14},{-5.08852e-16,
              4}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));
    end Liquid;

    package BaseClasses "Base classes (not for direct use)"
      extends Modelica.Icons.BasesPackage;

      model NullPhase "Empty BC for a phase (no species)"
        extends FCSys.BaseClasses.Icons.BCs.Single;

        parameter Boolean inclLinX=true "X" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(group="Axes with linear momentum included", compact=true));
        parameter Boolean inclLinY=false "Y" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(group="Axes with linear momentum included", compact=true));
        parameter Boolean inclLinZ=false "Z" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(group="Axes with linear momentum included", compact=true));

        FCSys.Connectors.ChemicalBus chemical
          "Multi-species connector for material"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

        FCSys.Connectors.RealInputBus u "Bus of inputs to specify conditions"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));
      end NullPhase;
    end BaseClasses;
  end ChemicalBus;

  package Chemical
    "<html>BCs for the <a href=\"modelica://FCSys.Connectors.ChemicalOutput\">ChemicalOutput</a> connector</html>"
    extends Modelica.Icons.Package;

    model Species
      "<html>BCs for the <a href=\"modelica://FCSys.Connectors.ChemicalOutput\">ChemicalOutput</a> connector</html>"

      import FCSys.BCs.Chemical.BaseClasses.BCTypeMaterial;
      import FCSys.BCs.Chemical.BaseClasses.BCTypeMechanical;
      import FCSys.BCs.Chemical.BaseClasses.BCTypeFluid;
      extends FCSys.BaseClasses.Icons.BCs.Single;

      replaceable package Data =
          FCSys.Characteristics.BaseClasses.Characteristic constrainedby
        FCSys.Characteristics.BaseClasses.Characteristic "Characteristic data"
        annotation (
        __Dymola_choicesAllMatching=true,
        Dialog(group="Material properties"),
        Placement(transformation(extent={{-60,40},{-40,60}}),
            iconTransformation(extent={{-10,90},{10,110}})));

      // Material
      parameter BCTypeMaterial material=BCTypeMaterial.PotentialPerTemperature
        "Type of BC"
        annotation (Dialog(group="Material", __Dymola_descriptionLabel=true));
      parameter Boolean internalMaterial=true "Use internal specification"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Material",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Modelica.Blocks.Sources.Constant materialSpec(k(start=Data.g(
              298.15*U.K)/(298.15*U.K))) if internalMaterial constrainedby
        Modelica.Blocks.Interfaces.SO "Internal specification" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Material",
          __Dymola_descriptionLabel=true,
          enable=internalMaterial),
        Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-90,10})));

      // X component of linear momentum
      parameter Boolean inclLinX=true "Include" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="X component of linear momentum",
          __Dymola_descriptionLabel=true,
          compact=true));
      parameter BCTypeMechanical linXBC=BCTypeMechanical.Velocity "Type of BC"
        annotation (Dialog(
          group="X component of linear momentum",
          enable=inclLinX,
          __Dymola_descriptionLabel=true));
      parameter Boolean internalLinX=true "Use internal specification"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="X component of linear momentum",
          enable=inclLinX,
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Modelica.Blocks.Sources.Constant linXSpec(k(start=0)) if
        inclLinX and internalLinX constrainedby Modelica.Blocks.Interfaces.SO
        "Internal specification" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="X component of linear momentum",
          __Dymola_descriptionLabel=true,
          enable=inclLinX and internalLinX),
        Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-50,10})));

      // Y component of linear momentum
      parameter Boolean inclLinY=false "Include" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Y component of linear momentum",
          __Dymola_descriptionLabel=true,
          compact=true));
      parameter BCTypeMechanical linYBC=BCTypeMechanical.Velocity "Type of BC"
        annotation (Dialog(
          group="Y component of linear momentum",
          enable=inclLinY,
          __Dymola_descriptionLabel=true));
      parameter Boolean internalLinY=true "Use internal specification"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Y component of linear momentum",
          enable=inclLinY,
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Modelica.Blocks.Sources.Constant linYSpec(k(start=0)) if
        inclLinY and internalLinY constrainedby Modelica.Blocks.Interfaces.SO
        "Internal specification" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Y component of linear momentum",
          enable=inclLinY and internalLinY,
          __Dymola_descriptionLabel=true,
          enable=internalLinY),
        Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={30,10})));

      // Z component of linear momentum
      parameter Boolean inclLinZ=false "Include" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Z component of linear momentum",
          __Dymola_descriptionLabel=true,
          compact=true));
      parameter BCTypeMechanical linZBC=BCTypeMechanical.Velocity "Type of BC"
        annotation (Dialog(
          group="Z component of linear momentum",
          enable=inclLinZ,
          __Dymola_descriptionLabel=true));
      parameter Boolean internalLinZ=true "Use internal specification"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Z component of linear momentum",
          enable=inclLinZ,
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Modelica.Blocks.Sources.Constant linZSpec(k(start=0)) if
        inclLinZ and internalLinZ constrainedby Modelica.Blocks.Interfaces.SO
        "Internal specification" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Z component of linear momentum",
          enable=inclLinZ and internalLinZ,
          __Dymola_descriptionLabel=true,
          enable=internalLinZ),
        Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-10,10})));

      // Enthalpy
      parameter BCTypeMechanical fluidBC=BCTypeFluid.EnthalpyMassic
        "Type of BC"
        annotation (Dialog(group="Enthalpy", __Dymola_descriptionLabel=true));
      parameter Boolean internalFluid=true "Use internal specification"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Enthalpy",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Modelica.Blocks.Sources.Constant fluidSpec(k(start=Data.h()/
              Data.m)) if internalFluid constrainedby
        Modelica.Blocks.Interfaces.SO "Internal specification" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Enthalpy",
          __Dymola_descriptionLabel=true,
          enable=internalFluid),
        Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={70,10})));

      FCSys.Connectors.ChemicalOutput chemical(
        final n_lin=n_lin,
        final formula=Data.formula,
        final m=Data.m)
        "Single-species connector for material, with advection of linear momentum and enthalpy"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
      FCSys.Connectors.RealInputBus u "Bus of inputs to specify conditions"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,60}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40})));

    protected
      final parameter Integer n_lin=countTrue({inclLinX,inclLinY,inclLinZ})
        "Number of components of linear momentum" annotation (Evaluate=true);

      FCSys.Connectors.RealInputInternal u_normal(final unit=if material ==
            BCTypeMaterial.PotentialPerTemperature then "1" else "N/T") if not
        internalMaterial "Material signal" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-70,30}),iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40})));
      FCSys.Connectors.RealInputInternal u_lin_x(final unit="l/T") if not
        internalLinX and inclLinX
        "Signal for the x component of linear momentum" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-30,30}),iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40})));
      FCSys.Connectors.RealInputInternal u_lin_y(final unit="l/T") if not
        internalLinY and inclLinY
        "Signal for the y component of linear momentum" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={10,30}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40})));
      FCSys.Connectors.RealInputInternal u_lin_z(final unit="l/T") if not
        internalLinZ and inclLinZ
        "Signal for the z component of linear momentum" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={50,30}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40})));
      FCSys.Connectors.RealInputInternal u_fluid(final unit="l2/T2") if not
        internalFluid "Fluid signal" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={90,30}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40})));
      FCSys.Connectors.RealInputInternal u_material_int
        "Internal material signal" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-70,-20}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40})));
      FCSys.Connectors.RealInputInternal u_lin_int[n_lin]
        "Internal mechanical signal" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={10,-20}),iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40})));
      FCSys.Connectors.RealInputInternal u_fluid_int "Internal fluid signal"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={90,-20}),iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40})));

    equation
      // Material
      if material == BCTypeMaterial.PotentialPerTemperature then
        chemical.muPerT = u_material_int;
      else
        chemical.Ndot = u_material_int;
      end if;

      // X component of linear momentum
      if inclLinX then
        //  if linXBC == BCTypeMechanical.Velocity then
        chemical.phi[1] = u_lin_int[1];
        //  end if;
      end if;

      // Y component of linear momentum
      if inclLinY then
        //  if linYBC == BCTypeMechanical.Velocity then
        chemical.phi[2] = u_lin_int[2];
        //  end if;
      end if;

      // Z component of linear momentum
      if inclLinZ then
        //  if linZBC == BCTypeMechanical.Velocity then
        chemical.phi[3] = u_lin_int[3];
        //  end if;
      end if;

      // Enthalpy
      //  if fluidBC== BCTypeFluid.EnthalpyMassic then
      chemical.hbar = u_fluid_int;
      //  end if;

      connect(u.normal, u_material) annotation (Line(
          points={{5.55112e-16,60},{5.55112e-16,40},{-70,40},{-70,30}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(u.linX, u_lin_x) annotation (Line(
          points={{5.55112e-16,60},{5.55112e-16,40},{-30,40},{-30,30}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(u.linY, u_lin_y) annotation (Line(
          points={{5.55112e-16,60},{5.55112e-16,40},{10,40},{10,30}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(u.linZ, u_lin_z) annotation (Line(
          points={{5.55112e-16,60},{5.55112e-16,40},{50,40},{50,30}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(u.fluid, u_fluid) annotation (Line(
          points={{5.55112e-16,60},{5.55112e-16,40},{90,40},{90,30}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(u_material_int, u_material) annotation (Line(
          points={{-70,-20},{-70,30}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u_lin_x, u_lin_int[1]) annotation (Line(
          points={{-30,30},{-30,-8},{10,-8},{10,-20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u_lin_y, u_lin_int[2]) annotation (Line(
          points={{10,30},{10,-20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u_lin_z, u_lin_int[3]) annotation (Line(
          points={{50,30},{50,-8},{10,-8},{10,-20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u_fluid, u_fluid_int) annotation (Line(
          points={{90,30},{90,-20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(materialSpec.y, u_material_int) annotation (Line(
          points={{-90,-1},{-90,-8},{-70,-8},{-70,-20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(linXSpec.y, u_lin_int[1]) annotation (Line(
          points={{-50,-1},{-50,-8},{10,-8},{10,-20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(linYSpec.y, u_lin_int[2]) annotation (Line(
          points={{30,-1},{30,-8},{10,-8},{10,-20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(linZSpec.y, u_lin_int[3]) annotation (Line(
          points={{-10,-1},{-10,-8},{10,-8},{10,-20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(fluidSpec.y, u_fluid_int) annotation (Line(
          points={{70,-1},{70,-8},{90,-8},{90,-20}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-120,
                -100},{120,100}}), graphics), Icon(coordinateSystem(
              preserveAspectRatio=true, extent={{-120,-100},{120,100}}),
            graphics));
    end Species;

    package BaseClasses "Base classes (not for direct use)"
      extends Modelica.Icons.BasesPackage;
      type BCTypeMaterial = enumeration(
          PotentialPerTemperature
            "Prescribed quotient of electrochemical potential and temperature",

          Current "Prescribed current") "Types of material BCs";

      type BCTypeMechanical = enumeration(
          Velocity "Prescribed velocity") "Types of mechanical BCs";
      type BCTypeFluid = enumeration(
          EnthalpyMassic "Prescribed massic enthalpy") "Types of fluid BCs";
    end BaseClasses;
  end Chemical;

  package InertAmagat
    "<html>BCs for the <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector</html>"
    extends Modelica.Icons.Package;

    model Phase
      "<html>BC for the <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Phases\">Phase</a> model</html>"
      extends FCSys.BaseClasses.Icons.BCs.Single;

      // Volume
      replaceable Volume.Volume volumeBC(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ) constrainedby Volume.BaseClasses.PartialBC
        "Condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Volume",__Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-90,-14},{-70,6}})));

      // X component of linear momentum
      parameter Boolean inclLinX=true "Include" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="X component of linear momentum",
          __Dymola_descriptionLabel=true,
          compact=true));
      replaceable Mechanical.Force linXBC(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ) if inclLinX constrainedby
        Mechanical.BaseClasses.PartialBC "Condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="X component of linear momentum",
          enable=inclLinX,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-50,-14},{-30,6}})));

      // Y component of linear momentum
      parameter Boolean inclLinY=true "Include" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Y component of linear momentum",
          __Dymola_descriptionLabel=true,
          compact=true));
      replaceable Mechanical.Force linYBC(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ) if inclLinY constrainedby
        Mechanical.BaseClasses.PartialBC "Condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Y component of linear momentum",
          enable=inclLinY,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-10,-14},{10,6}})));

      // Z component of linear momentum
      parameter Boolean inclLinZ=true "Include" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Z component of linear momentum",
          __Dymola_descriptionLabel=true,
          compact=true));
      replaceable Mechanical.Force linZBC(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ) if inclLinZ constrainedby
        Mechanical.BaseClasses.PartialBC "Condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Z component of linear momentum",
          enable=inclLinZ,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{30,-14},{50,6}})));

      // Heat
      replaceable Thermal.HeatFlowRate thermal(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ) constrainedby Thermal.BaseClasses.PartialBC
        "Condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Heat",__Dymola_descriptionLabel=true),
        Placement(transformation(extent={{70,-14},{90,6}})));

      FCSys.Connectors.RealInputBus u "Input bus for external signal sources"
        annotation (HideResult=not (internalVolume or internalLinX or
            internalLinY or internalLinZ or internalThermal), Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={2,40})));

      FCSys.Connectors.InertAmagat inert(final n_lin=countTrue({inclLinX,
            inclLinY,inclLinZ}))
        "Single-species connector for linear momentum and heat, with additivity of pressure"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}}),
            iconTransformation(extent={{-10,-50},{10,-30}})));

    equation
      // Volume
      connect(volumeBC.inert, inert) annotation (Line(
          points={{-80,-8},{-80,-20},{5.55112e-16,-20},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(u.volume, volumeBC.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,20},{-80,20},{-80,6.66134e-16}},

          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // X component of linear momentum
      connect(linXBC.inert, inert) annotation (Line(
          points={{-40,-8},{-40,-20},{5.55112e-16,-20},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(u.linX, linXBC.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,20},{-40,20},{-40,6.66134e-16}},

          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // Y component of linear momentum
      connect(linYBC.inert, inert) annotation (Line(
          points={{6.10623e-16,-8},{6.10623e-16,-20},{5.55112e-16,-20},{
              5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(u.linY, linYBC.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,6.66134e-16},{6.10623e-16,
              6.66134e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // Z component of linear momentum
      connect(linZBC.inert, inert) annotation (Line(
          points={{40,-8},{40,-20},{5.55112e-16,-20},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(u.linZ, linZBC.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,20},{40,20},{40,6.66134e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // Heat
      connect(thermal.inert, inert) annotation (Line(
          points={{80,-8},{80,-20},{5.55112e-16,-20},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(u.thermal, thermal.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,20},{80,20},{80,6.66134e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
    end Phase;

    package Volume "BCs for additivity of volume"
      extends Modelica.Icons.Package;

      model Pressure "Prescribed pressure"
        extends BaseClasses.PartialBC(
          final bCType=BaseClasses.BCType.Pressure,
          u(final unit="m/(l.T2)"),
          spec(k(start=U.atm)));

      equation
        inert.p = u_final;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="volumeBC");
      end Pressure;

      model Volume "Prescribed volume"
        extends BaseClasses.PartialBC(
          final bCType=BaseClasses.BCType.Volume,
          u(final unit="l3"),
          spec(k(start=U.cc)));
      equation
        inert.V = u_final;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="volumeBC");
      end Volume;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialBC "Partial model for a BC for volume"
          extends FCSys.BCs.InertAmagat.BaseClasses.PartialBC;

          constant BCType bCType "Type of BC";
          // Note:  This is included so that the type of BC is recorded with the
          // results.

        equation
          inert.mPhidot = zeros(n_lin) "No force";
          inert.Qdot = 0 "Adiabatic";
          annotation (defaultComponentName="volumeBC");
        end PartialBC;

        type BCType = enumeration(
            Volume "Volume",
            Pressure "Pressure") "Types of BCs";
      end BaseClasses;
    end Volume;

    package Mechanical "Mechanical BCs"
      extends Modelica.Icons.Package;
      model Velocity "Prescribed velocity"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Velocity,
            u(final unit="l/T"));
      equation
        inert.phi[linAxes[axis]] = u_final;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="mechanicalBC");
      end Velocity;

      model Force "Prescribed force"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Force, u(
              final unit="l.m/T2"));
      equation
        inert.mPhidot[linAxes[axis]] = u_final;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="mechanicalBC");
      end Force;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialBC "Partial model for a mechanical BC"
          extends FCSys.BCs.InertAmagat.BaseClasses.PartialBC;

          parameter Axis axis=Axis.x "Axis" annotation (HideResult=true);

          constant BCType bCType "Type of BC";
          // Note:  This is included so that the type of BC is recorded with the
          // results.

        protected
          final parameter Integer cartAxes[n_lin]=index({inclLinX,inclLinY,
              inclLinZ})
            "Cartesian-axis indices of the axes of linear momentum";
          final parameter Integer linAxes[Axis]=enumerate({inclLinX,inclLinY,
              inclLinZ})
            "Linear momentum component indices of the Cartesian axes";

        equation
          inert.V = 0 "No volume";
          for i in 1:n_lin loop
            if cartAxes[i] <> axis then
              inert.mPhidot[i] = 0 "No force along the other axes";
            end if;
          end for;
          inert.Qdot = 0 "Adiabatic";

          annotation (defaultComponentName="mechanicalBC");
        end PartialBC;

        type BCType = enumeration(
            Velocity "Velocity",
            Force "Force") "Types of BCs";
      end BaseClasses;
    end Mechanical;

    package Thermal "Thermal BCs"
      extends Modelica.Icons.Package;

      model Temperature "Prescribed temperature"
        extends BaseClasses.PartialBC(
          final bCType=BaseClasses.BCType.Temperature,
          u(final unit="l2.m/(N.T2)", displayUnit="K"),
          spec(k(start=298.15*U.K)));
      equation
        inert.T = u_final;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="thermal");
      end Temperature;

      model HeatFlowRate "Prescribed heat flow rate"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.HeatFlowRate,
            u(final unit="l2.m/T3"));
      equation
        inert.Qdot = u_final;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="thermal");
      end HeatFlowRate;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialBC "Partial model for a thermal BC"
          extends FCSys.BCs.InertAmagat.BaseClasses.PartialBC;

          constant BCType bCType "Type of BC";
          // Note:  This is included so that the type of BC is recorded with the
          // results.

        equation
          inert.V = 0 "No volume";
          inert.mPhidot = zeros(n_lin) "No force";
          annotation (defaultComponentName="thermal");
        end PartialBC;

        type BCType = enumeration(
            Temperature "Temperature",
            HeatFlowRate "Heat flow rate") "Types of BCs";
      end BaseClasses;
    end Thermal;

    package BaseClasses "Base classes (not for direct use)"
      extends Modelica.Icons.BasesPackage;
      partial model PartialBC "Partial model for a BC"
        extends FCSys.BaseClasses.Icons.BCs.Single;

        parameter Boolean inclLinX=true "X" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(group="Axes with linear momentum included", compact=true));
        parameter Boolean inclLinY=false "Y" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(group="Axes with linear momentum included", compact=true));
        parameter Boolean inclLinZ=false "Z" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(group="Axes with linear momentum included", compact=true));

        parameter Boolean internal=true "Use internal specification"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(__Dymola_descriptionLabel=true));

        replaceable Modelica.Blocks.Sources.Constant spec if internal
          constrainedby Modelica.Blocks.Interfaces.SO "Internal specification"
          annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(__Dymola_descriptionLabel=true, enable=internal),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-30,30})));

        FCSys.Connectors.RealInput u if not internal "Value of BC" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));
        FCSys.Connectors.InertAmagat inert(final n_lin=n_lin)
          "Connector for linear momentum and heat, with additivity of volume"
          annotation (HideResult=true, Placement(transformation(extent={{-10,-50},
                  {10,-30}})));

      protected
        final parameter Integer n_lin=countTrue({inclLinX,inclLinY,inclLinZ})
          "Number of components of linear momentum" annotation (Evaluate=true);

        FCSys.Connectors.RealInputInternal u_final "Final value of BC"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-10})));

      equation
        connect(u, u_final) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,27.5},{5.55112e-16,27.5},{
                5.55112e-16,15},{5.55112e-16,-10},{5.55112e-16,-10}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(spec.y, u_final) annotation (Line(
            points={{-30,19},{-30,10},{5.55112e-16,10},{5.55112e-16,-10}},
            color={0,0,127},
            smooth=Smooth.None));
      end PartialBC;
    end BaseClasses;
  end InertAmagat;

  package InertDalton
    "<html>BCs for the <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector</html>"
    extends Modelica.Icons.Package;

    model Species
      "<html>BC for the <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model</html>"

      extends FCSys.BaseClasses.Icons.BCs.Single;

      // Pressure
      replaceable Pressure.Volume pressureBC(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ) constrainedby Pressure.BaseClasses.PartialBC
        "Condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Pressure", __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-90,-14},{-70,6}})));

      // X component of linear momentum
      parameter Boolean inclLinX=true "Include" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="X component of linear momentum",
          __Dymola_descriptionLabel=true,
          compact=true));
      replaceable Mechanical.Force linXBC(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ) if inclLinX constrainedby
        Mechanical.BaseClasses.PartialBC "Condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="X component of linear momentum",
          enable=inclLinX,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-50,-14},{-30,6}})));

      // Y component of linear momentum
      parameter Boolean inclLinY=true "Include" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Y component of linear momentum",
          __Dymola_descriptionLabel=true,
          compact=true));
      replaceable Mechanical.Force linYBC(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ) if inclLinY constrainedby
        Mechanical.BaseClasses.PartialBC "Condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Y component of linear momentum",
          enable=inclLinY,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-10,-14},{10,6}})));

      // Z component of linear momentum
      parameter Boolean inclLinZ=true "Include" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Z component of linear momentum",
          __Dymola_descriptionLabel=true,
          compact=true));
      replaceable Mechanical.Force linZBC(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ) if inclLinZ constrainedby
        Mechanical.BaseClasses.PartialBC "Condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Z component of linear momentum",
          enable=inclLinZ,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{30,-14},{50,6}})));

      // Heat
      replaceable Thermal.HeatFlowRate thermal(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ) constrainedby Thermal.BaseClasses.PartialBC
        "Condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Heat", __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{70,-14},{90,6}})));

      FCSys.Connectors.RealInputBus u "Input bus for external signal sources"
        annotation (HideResult=not (internalVol or internalLinX or internalLinY
             or internalLinZ or internalThermal), Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40})));

      FCSys.Connectors.InertDalton inert(final n_lin=countTrue({inclLinX,
            inclLinY,inclLinZ}))
        "Single-species connector for linear momentum and heat, with additivity of pressure"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}}),
            iconTransformation(extent={{-10,-50},{10,-30}})));

    equation
      // Pressure
      connect(pressureBC.inert, inert) annotation (Line(
          points={{-80,-8},{-80,-20},{5.55112e-16,-20},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(u.pressure, pressureBC.u) annotation (Line(
          points={{5.55112e-16,40},{0,40},{0,20},{-80,20},{-80,6.66134e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // X component of linear momentum
      connect(linXBC.inert, inert) annotation (Line(
          points={{-40,-8},{-40,-20},{5.55112e-16,-20},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(u.linX, linXBC.u) annotation (Line(
          points={{5.55112e-16,40},{0,40},{0,20},{-40,20},{-40,6.66134e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // Y component of linear momentum
      connect(linYBC.inert, inert) annotation (Line(
          points={{6.10623e-16,-8},{6.10623e-16,-20},{5.55112e-16,-20},{
              5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(u.linY, linYBC.u) annotation (Line(
          points={{5.55112e-16,40},{0,40},{0,6.66134e-16},{6.10623e-16,
              6.66134e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // Z component of linear momentum
      connect(linZBC.inert, inert) annotation (Line(
          points={{40,-8},{40,-20},{5.55112e-16,-20},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(u.linZ, linZBC.u) annotation (Line(
          points={{5.55112e-16,40},{0,40},{0,20},{40,20},{40,6.66134e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // Heat
      connect(thermal.inert, inert) annotation (Line(
          points={{80,-8},{80,-20},{5.55112e-16,-20},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(u.thermal, thermal.u) annotation (Line(
          points={{5.55112e-16,40},{0,40},{0,20},{80,20},{80,6.66134e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
    end Species;

    package Pressure "BCs for additivity of pressure"
      extends Modelica.Icons.Package;

      model Volume "Prescribed volume"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Volume, u(
              final unit="l3"));
      equation
        inert.V = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="pressureBC");
      end Volume;

      model Pressure "Prescribed pressure"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Pressure,
            u(final unit="m/(l.T2)"));
      equation
        inert.p = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="pressureBC");
      end Pressure;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialBC "Partial model for a BC for pressure"
          extends FCSys.BCs.InertDalton.BaseClasses.PartialBC;

          constant BCType bCType "Type of BC";
          // Note:  This is included so that the type of BC is recorded with the
          // results.
        equation
          inert.mPhidot = zeros(n_lin) "No force";
          inert.Qdot = 0 "No heat flow";
          annotation (defaultComponentName="pressureBC");
        end PartialBC;

        type BCType = enumeration(
            Volume "Volume",
            Pressure "Pressure") "Types of BCs";
      end BaseClasses;
    end Pressure;

    package Mechanical "Mechanical BCs"
      extends Modelica.Icons.Package;
      model Velocity "Prescribed velocity"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Velocity,
            u(final unit="l/T"));
      equation
        inert.phi[linAxes[axis]] = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="mechanicalBC");
      end Velocity;

      model Force "Prescribed force"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Force, u(
              final unit="l.m/T2"));
      equation
        inert.mPhidot[linAxes[axis]] = u_final;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="mechanicalBC");
      end Force;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialBC "Partial model for a mechanical BC"
          extends FCSys.BCs.InertDalton.BaseClasses.PartialBC;

          parameter Axis axis=Axis.x "Axis" annotation (HideResult=true);

          constant BCType bCType "Type of BC";
          // Note:  This is included so that the type of BC is recorded with the
          // results.

        protected
          final parameter Integer cartAxes[n_lin]=index({inclLinX,inclLinY,
              inclLinZ})
            "Cartesian-axis indices of the axes of linear momentum";
          final parameter Integer linAxes[Axis]=enumerate({inclLinX,inclLinY,
              inclLinZ})
            "Linear momentum component indices of the Cartesian axes";

        equation
          inert.p = 0 "No pressure";
          for i in 1:n_lin loop
            if cartAxes[i] <> axis then
              inert.mPhidot[i] = 0 "No force along the other axes";
            end if;
          end for;
          inert.Qdot = 0 "Adiabatic";

          annotation (defaultComponentName="mechanicalBC");
        end PartialBC;

        type BCType = enumeration(
            Velocity "Velocity",
            Force "Force") "Types of BCs";
      end BaseClasses;
    end Mechanical;

    package Thermal "BCs for heat"
      extends Modelica.Icons.Package;

      model Temperature "Prescribed temperature"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Temperature,
            u(final unit="l2.m/(N.T2)", displayUnit="K"));

      equation
        inert.T = u_final;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="thermal");
      end Temperature;

      model HeatFlowRate "Prescribed heat flow rate"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.HeatFlowRate,
            u(final unit="l2.m/T3"));

      equation
        inert.Qdot = u_final;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="thermal");
      end HeatFlowRate;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialBC "Partial model for a thermal BC"
          extends FCSys.BCs.InertDalton.BaseClasses.PartialBC;

          constant BCType bCType "Type of BC";
          // Note:  This is included so that the type of BC is recorded with the
          // results.

        equation
          inert.p = 0 "No pressure";
          inert.mPhidot = zeros(n_lin) "No force";
          annotation (defaultComponentName="thermal");
        end PartialBC;

        type BCType = enumeration(
            Temperature "Temperature",
            HeatFlowRate "Heat flow rate") "Types of BCs";
      end BaseClasses;
    end Thermal;

    package BaseClasses "Base classes (not for direct use)"
      extends Modelica.Icons.BasesPackage;
      partial model PartialBC "Partial model for a BC"
        extends FCSys.BaseClasses.Icons.BCs.Single;

        parameter Boolean inclLinX=true "X" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(group="Axes with linear momentum included", compact=true));
        parameter Boolean inclLinY=false "Y" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(group="Axes with linear momentum included", compact=true));
        parameter Boolean inclLinZ=false "Z" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(group="Axes with linear momentum included", compact=true));

        parameter Boolean internal=true "Use internal specification"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(__Dymola_descriptionLabel=true));

        replaceable Modelica.Blocks.Sources.Constant spec if internal
          constrainedby Modelica.Blocks.Interfaces.SO "Internal specification"
          annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(__Dymola_descriptionLabel=true, enable=internal),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-30,30})));

        FCSys.Connectors.RealInput u if not internal "Value of BC" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));
        FCSys.Connectors.InertDalton inert(final n_lin=n_lin)
          "Connector for linear momentum and heat, with additivity of pressure"
          annotation (HideResult=true, Placement(transformation(extent={{-10,-50},
                  {10,-30}})));

      protected
        final parameter Integer n_lin=countTrue({inclLinX,inclLinY,inclLinZ})
          "Number of components of linear momentum" annotation (Evaluate=true);

        FCSys.Connectors.RealInputInternal u_final "Final value of BC"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-10})));

      equation
        connect(u, u_final) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,27.5},{5.55112e-16,27.5},{
                5.55112e-16,15},{5.55112e-16,-10},{5.55112e-16,-10}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(spec.y, u_final) annotation (Line(
            points={{-30,19},{-30,10},{5.55112e-16,10},{5.55112e-16,-10}},
            color={0,0,127},
            smooth=Smooth.None));
      end PartialBC;
    end BaseClasses;
  end InertDalton;

  package FaceBus
    "<html>BCs for a <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector</html>"
    extends Modelica.Icons.Package;

    model Subregion
      "<html>BC for a face of a <a href=\"modelica://FCSys.Regions.Region\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with efforts by default</html>"
      extends FCSys.BaseClasses.Icons.BCs.Single;

      parameter Axis axis=Axis.x "Axis normal to the face";

      Phases.Gas gas(final axis=axis) "Gas" annotation (Dialog(group="Phases",
            __Dymola_descriptionLabel=true), Placement(transformation(extent={{
                -10,-10},{10,10}})));

      Phases.Graphite graphite(final axis=axis) "Graphite" annotation (Dialog(
            group="Phases", __Dymola_descriptionLabel=true), Placement(
            transformation(extent={{-10,-10},{10,10}})));

      Phases.Ionomer ionomer(final axis=axis) "Ionomer" annotation (Dialog(
            group="Phases", __Dymola_descriptionLabel=true), Placement(
            transformation(extent={{-10,-10},{10,10}})));

      Phases.Liquid liquid(final axis=axis) "Liquid" annotation (Dialog(group=
              "Phases", __Dymola_descriptionLabel=true), Placement(
            transformation(extent={{-10,-10},{10,10}})));

      FCSys.Connectors.FaceBus face
        "Connector for linear momentum and heat of multiple species"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}}),
            iconTransformation(extent={{-10,-50},{10,-30}})));
      FCSys.Connectors.RealInputBus u "Bus of inputs to specify conditions"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40})));

    equation
      connect(gas.face, face.gas) annotation (Line(
          points={{6.10623e-16,-4},{7.35523e-16,-4},{8.60423e-16,-40},{
              5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(graphite.face, face.graphite) annotation (Line(
          points={{6.10623e-16,-4},{7.35523e-16,-4},{8.60423e-16,-40},{
              5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(ionomer.face, face.ionomer) annotation (Line(
          points={{6.10623e-16,-4},{8.60423e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));
      connect(liquid.face, face.liquid) annotation (Line(
          points={{6.10623e-16,-4},{8.60423e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));
      connect(u.gas, gas.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
              6.10623e-16,4}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(u.graphite, graphite.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
              6.10623e-16,4}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(u.ionomer, ionomer.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
              6.10623e-16,4}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(u.liquid, liquid.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
              6.10623e-16,4}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="subregionFaceBC",
        Diagram(graphics));
    end Subregion;

    model SubregionFlow
      "<html>BC for a face of a <a href=\"modelica://FCSys.Subregions.Subregion\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with flows by default</html>"
      extends FaceBus.Subregion(
        gas(
          H2(
            redeclare replaceable Face.Material.Current normal(spec(k=0)),
            redeclare replaceable Face.Mechanical.Force transverseX,
            redeclare replaceable Face.Mechanical.Force transverseY,
            redeclare replaceable Face.Mechanical.Force transverseZ,
            redeclare replaceable Face.Thermal.HeatFlowRate thermal(spec(k=0))),

          H2O(
            redeclare replaceable Face.Material.Current normal(spec(k=0)),
            redeclare replaceable Face.Mechanical.Force transverseX,
            redeclare replaceable Face.Mechanical.Force transverseY,
            redeclare replaceable Face.Mechanical.Force transverseZ,
            redeclare replaceable Face.Thermal.HeatFlowRate thermal(spec(k=0))),

          N2(
            redeclare replaceable Face.Material.Current normal(spec(k=0)),
            redeclare replaceable Face.Mechanical.Force transverseX,
            redeclare replaceable Face.Mechanical.Force transverseY,
            redeclare replaceable Face.Mechanical.Force transverseZ,
            redeclare replaceable Face.Thermal.HeatFlowRate thermal(spec(k=0))),

          O2(
            redeclare replaceable Face.Material.Current normal(spec(k=0)),
            redeclare replaceable Face.Mechanical.Force transverseX,
            redeclare replaceable Face.Mechanical.Force transverseY,
            redeclare replaceable Face.Mechanical.Force transverseZ,
            redeclare replaceable Face.Thermal.HeatFlowRate thermal(spec(k=0)))),

        graphite('C+'(
            redeclare replaceable Face.Material.Current normal(spec(k=0)),
            redeclare replaceable Face.Mechanical.Force transverseX,
            redeclare replaceable Face.Mechanical.Force transverseY,
            redeclare replaceable Face.Mechanical.Force transverseZ,
            redeclare replaceable Face.Thermal.HeatFlowRate thermal(spec(k=0))),
            'e-'(
            redeclare replaceable Face.Material.Current normal(spec(k=0)),
            redeclare replaceable Face.Mechanical.Force transverseX,
            redeclare replaceable Face.Mechanical.Force transverseY,
            redeclare replaceable Face.Mechanical.Force transverseZ,
            redeclare replaceable Face.Thermal.HeatFlowRate thermal(spec(k=0)))),

        ionomer);

      redeclare replaceable Face.Mechanical.Force transverseX;
      redeclare replaceable Face.Mechanical.Force Error;

      annotation (defaultComponentPrefixes="replaceable",defaultComponentName=
            "subregionFaceBC");
    end SubregionFlow;

    model SubregionClosed
      "<html>BC for a face of a <a href=\"modelica://FCSys.Subregions.Subregion\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with efforts except closed by default</html>"
      extends FaceBus.Subregion(
        gas(
          H2(isobaric=true),
          H2O(isobaric=true),
          N2(isobaric=true),
          O2(isobaric=true)),
        graphite('C+'(isobaric=true), 'e-'(thermoOpt=ThermoOpt.ClosedDiabatic)),

        ionomer(Error),
        H2O(isobaric=true),
        'H+'(isobaric=true));

      annotation (defaultComponentPrefixes="replaceable",defaultComponentName=
            "subregionFaceBC");
    end SubregionClosed;

    model SubregionClosedAdiabatic
      "<html>BC for a face of a <a href=\"modelica://FCSys.Subregions.Subregion\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with flows except zero shear velocities by default</html>"
      extends FaceBus.Subregion(
        gas(
          H2(isobaric=true),
          H2O(isobaric=true),
          N2(isobaric=true),
          O2(isobaric=true)),
        graphite('C+'(isobaric=true), 'e-'(thermoOpt=ThermoOpt.ClosedAdiabatic)),

        ionomer(Error),
        H2O(isobaric=true),
        'H+'(isobaric=true));

      annotation (defaultComponentPrefixes="replaceable",defaultComponentName=
            "subregionFaceBC");
    end SubregionClosedAdiabatic;

    package Phases
      "<html>BCs for the <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Phase\">Phase</a> model (multi-species)</html>"
      extends Modelica.Icons.Package;

      model Gas "BC for gas"

        extends BaseClasses.NullPhase;

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
        Face.Species H2(final axis=axis) if inclH2 "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2), Placement(transformation(extent={{-10,-10},{10,10}})));

        parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species H2O(final axis=axis) if inclH2O "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,10}})));

        parameter Boolean inclN2=false "<html>Nitrogen (N<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species N2(final axis=axis) if inclN2 "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclN2), Placement(transformation(extent={{-10,-10},{10,10}})));

        parameter Boolean inclO2=false "<html>Oxygen (O<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species O2(final axis=axis) if inclO2 "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclO2), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // H2
        connect(H2.face.normal, face.H2.normal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.transverseX, face.H2.transverseX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.transverseY, face.H2.transverseY) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.transverseZ, face.H2.transverseZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.thermal, face.H2.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.H2, H2.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2O
        connect(H2O.face.normal, face.H2O.normal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.transverseX, face.H2O.transverseX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.transverseY, face.H2O.transverseY) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.transverseZ, face.H2O.transverseZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.thermal, face.H2O.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.H2O, H2O.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // N2
        connect(N2.face.normal, face.N2.normal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.transverseX, face.N2.transverseX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.transverseY, face.N2.transverseY) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.transverseZ, face.N2.transverseZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.thermal, face.N2.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.N2, N2.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // O2
        connect(O2.face.normal, face.O2.normal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.transverseX, face.O2.transverseX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.transverseY, face.O2.transverseY) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.transverseZ, face.O2.transverseZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.thermal, face.O2.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.O2, O2.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (Diagram(graphics));
      end Gas;

      model Graphite "BC for graphite"

        extends BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean HideResult=true;
        parameter Boolean choices(__Dymola_checkBox=true);
        parameter Boolean Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true);
        Face.Species 'C+'(final axis=axis, isobaric=true) if '' + ERROR
          annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='' + ERROR annotation (
              Evaluate=true,
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true))));
        Face.Species 'e-'(final axis=axis) if 'incle-' "Model" annotation (
            Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='incle-'), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // C+
        connect('C+'.face.normal, face.'C+'.normal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('C+'.face.transverseX, face.'C+'.transverseX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('C+'.face.transverseY, face.'C+'.transverseY) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('C+'.face.transverseZ, face.'C+'.transverseZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('C+'.face.thermal, face.'C+'.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'C+', 'C+'.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // e-
        connect('e-'.face.normal, face.'e-'.normal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.transverseX, face.'e-'.transverseX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.transverseY, face.'e-'.transverseY) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.transverseZ, face.'e-'.transverseZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.thermal, face.H2.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'e-', 'e-'.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

      end Graphite;

      model Ionomer "BC for ionomer"

        extends BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean 'inclC19HF37O5S-'=false
          "<html>Nafion sulfonate minus (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species 'C19HF37O5S-'(final axis=axis,isobaric=true) if
          'inclC19HF37O5S-' "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclC19HF37O5S-'), Placement(transformation(extent={{-10,-10},
                  {10,10}})));

        parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species H2O(final axis=axis) if inclH2O "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,10}})));

        parameter Boolean 'inclH+'=false "<html>Protons (H<sup>+</sup>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species 'H+'(final axis=axis) if 'inclH+' "Model" annotation (
            Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclH+'), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // C19HF37O5S-
        connect('C19HF37O5S-'.face.normal, face.'C19HF37O5S-'.normal)
          annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('C19HF37O5S-'.face.transverseX, face.'C19HF37O5S-'.transverseX)
          annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('C19HF37O5S-'.face.transverseY, face.'C19HF37O5S-'.transverseY)
          annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('C19HF37O5S-'.face.transverseZ, face.'C19HF37O5S-'.transverseZ)
          annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('C19HF37O5S-'.face.thermal, face.'C19HF37O5S-'.thermal)
          annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'C19HF37O5S-', 'C19HF37O5S-'.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H+
        connect('H+'.face.normal, face.'H+'.normal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.transverseX, face.'H+'.transverseX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.transverseY, face.'H+'.transverseY) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.transverseZ, face.'H+'.transverseZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.thermal, face.'H+'.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'H+', 'H+'.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

      end Ionomer;

      model Liquid "BC for liquid"

        extends BaseClasses.NullPhase;

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
        Face.Species H2O(final axis=axis) if inclH2O "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // H2O
        connect(H2O.face.normal, face.H2O.normal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.transverseX, face.H2O.transverseX) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.transverseY, face.H2O.transverseY) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.transverseZ, face.H2O.transverseZ) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.thermal, face.H2O.thermal) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.H2O, H2O.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
      end Liquid;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        model NullPhase "Empty BC for a phase (no species)"
          extends FCSys.BaseClasses.Icons.BCs.Single;

          parameter Axis axis "Axis material to the face";

          FCSys.Connectors.FaceBus face
            "Multi-species connector for linear momentum and heat"
            annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
          FCSys.Connectors.RealInputBus u "Bus of inputs to specify conditions"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,40})));
        end NullPhase;
      end BaseClasses;
    end Phases;

    annotation (Documentation(info="<html><p>Since the connectors in
<a href=\"modelica://FCSys\">FCSys</a> are hierarchical
(see the <a href=\"modelica://FCSys.Connectors\">Connectors</a> package),
the models for the boundary conditions must be as well.  A
<a href=\"modelica://FCSys.Connectors.Face\">Face</a>
connector (<a href=\"modelica://FCSys.Connectors.Face\">Face</a>,
<a href=\"modelica://FCSys.Connectors.Face\">Face</a>, or
<a href=\"modelica://FCSys.Connectors.Face\">Face</a>)
is used in <a href=\"modelica://FCSys.Subregions.Species\">Species</a> models,
and there is a corresponding <a href=\"modelica://FCSys.BCs.Face.Species.Species\">Species
boundary condition</a> model in this package. The
<a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a>
connector is used in <a href=\"modelica://FCSys.Subregions.Phases\">Phase</a> models,
and there are corresponding <a href=\"modelica://FCSys.BCs.Face.Phases\">Phase
boundary condition</a> models.  The
<a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector is nested once more
in models such as the <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a>,
and there is a corresponding <a href=\"modelica://FCSys.BCs.Face.Subregion\">Subregion
boundary condition</a> model.
</p></html>"));
  end FaceBus;

  package Face
    "<html>BCs for a <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"
    model Species
      "<html>BC for a face of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"
      extends BaseClasses.PartialSpecies;

      parameter Axis axis "Axis normal to the face";

      // X-axis linear momentum
      parameter Boolean inviscidX=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="X-axis linear momentum (if not normal)",
          enable=axis <> 1,
          __Dymola_descriptionLabel=true));
      // Note:  Dymola 7.4 doesn't recognize enumerations in the dialog enable
      // option, e.g.,
      //     enable=axis <> Axis.x.
      // Therefore, the values of the enumerations are specified numerically.
      replaceable Mechanical.Velocity transverseX(spec(k=0)) if axis <> Axis.x
         and not inviscidX constrainedby Mechanical.BaseClasses.PartialBC
        "Type of condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="X-axis linear momentum (if not normal)",
          enable=axis <> 1 and not inviscidX,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-40,-14},{-20,6}})));

      // Y-axis linear momentum
      parameter Boolean inviscidY=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="Y-axis linear momentum (if not normal)",
          enable=axis <> 2,
          __Dymola_descriptionLabel=true));
      replaceable Mechanical.Velocity transverseY(spec(k=0)) if axis <> Axis.y
         and not inviscidY constrainedby Mechanical.BaseClasses.PartialBC
        "Type of condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Y-axis linear momentum (if not normal)",
          enable=axis <> 2 and not inviscidY,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-10,-14},{10,6}})));

      // Z-axis linear momentum
      parameter Boolean inviscidZ=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="Z-axis linear momentum (if not normal)",
          enable=axis <> 3,
          __Dymola_descriptionLabel=true));
      replaceable Mechanical.Velocity transverseZ(spec(k=0)) if axis <> Axis.z
         and not inviscidZ constrainedby Mechanical.BaseClasses.PartialBC
        "Type of condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Z-axis linear momentum (if not normal)",
          enable=axis <> 3 and not inviscidZ,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{20,-14},{40,6}})));

      FCSys.Connectors.Face face(
        final isobaric=isobaric,
        final inviscidX=inviscidX,
        final inviscidY=inviscidY,
        final inviscidZ=inviscidZ)
        "Single-species connector for  linear momentum and heat" annotation (
          Placement(transformation(extent={{-10,-50},{10,-30}}),
            iconTransformation(extent={{-10,-50},{10,-30}})));

    equation
      // Material
      connect(material.normal, face.normal) annotation (Line(
          points={{-60,-8},{-60,-20},{5.55112e-16,-20},{5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));

      // X-axis linear momentum
      connect(transverseX.mechanical, face.transverseX) annotation (Line(
          points={{-30,-8},{-30,-20},{5.55112e-16,-20},{5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(u.transverseX, transverseX.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,20},{-30,20},{-30,6.66134e-16}},

          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // Y-axis linear momentum
      connect(transverseY.mechanical, face.transverseY) annotation (Line(
          points={{6.10623e-16,-8},{6.10623e-16,-20},{0,-20},{0,-40},{
              5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));

      connect(u.transverseY, transverseY.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,20},{6.10623e-16,20},{
              6.10623e-16,6.66134e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // Z-axis linear momentum
      connect(transverseZ.mechanical, face.transverseZ) annotation (Line(
          points={{30,-8},{30,-20},{5.55112e-16,-20},{5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(u.transverseZ, transverseZ.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,20},{30,20},{30,6.66134e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // Heat
      connect(thermal.thermal, face.thermal) annotation (Line(
          points={{60,-8},{60,-20},{5.55112e-16,-20},{5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));

      annotation (Diagram(graphics));
    end Species;
    extends Modelica.Icons.Package;

    package Material "Material BCs"
      extends Modelica.Icons.Package;

      model Density "Prescribed density"

        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Density,
            u(final unit="N/l3"));

      equation
        material.rho = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="material");
      end Density;

      model Current "Prescribed current"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Current,
            u(final unit="N/T"));

      equation
        material.Ndot = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="material");
      end Current;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialBC "Partial model for a material BC"
          extends Face.BaseClasses.PartialBC;

          constant BCType bCType "Type of BC";
          // Note:  This is included so that the type of BC is recorded with the
          // results.

          FCSys.Connectors.BaseClasses.Normal material
            "Material subconnector for the face"
            annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
          annotation (defaultComponentName="material");
        end PartialBC;

        type BCType = enumeration(
            Density "Density",
            Current "Current") "Types of BCs";
      end BaseClasses;
    end Material;

    package Mechanical "Mechanical BCs"
      extends Modelica.Icons.Package;

      model Velocity "Prescribed velocity"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Velocity,
            u(final unit="l/T"));

      equation
        mechanical.phi = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="mechanicalBC");
      end Velocity;

      model Force "Prescribed shear force"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Force, u(
              final unit="l.m/T2"));

      equation
        mechanical.mPhidot = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="mechanicalBC");
      end Force;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialBC "Partial model for a mechanical BC"
          extends Face.BaseClasses.PartialBC;

          constant Mechanical.BaseClasses.BCType bCType "Type of BC";
          // Note:  This is included so that the type of BC is recorded with the
          // results.

          FCSys.Connectors.BaseClasses.Transverse mechanical
            "Mechanical subconnector for the face"
            annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
          annotation (defaultComponentName="mechanicalBC");
        end PartialBC;

        type BCType = enumeration(
            Velocity "Velocity",
            Force "Shear force") "Types of BCs";
      end BaseClasses;
    end Mechanical;

    package Thermal "Thermal BCs"
      extends Modelica.Icons.Package;

      model Temperature "Prescribed temperature"
        extends Thermal.BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Temperature,
            u(final unit="l2.m/(N.T2)", displayUnit="K"));

      equation
        thermal.T = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="thermal");
      end Temperature;

      model HeatFlowRate "Prescribed heat flow rate"
        extends Thermal.BaseClasses.PartialBC(final bCType=BaseClasses.BCType.HeatFlowRate,
            u(final unit="l2.m/T3"));

      equation
        thermal.Qdot = u_final;
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="thermal",
          Diagram(graphics));
      end HeatFlowRate;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialBC "Partial model for a thermal BC"
          extends Face.BaseClasses.PartialBC;

          constant BCType bCType "Type of BC";
          // Note:  This is included so that the type of BC is recorded with the
          // results.

          FCSys.Connectors.Thermal thermal "Thermal subconnector for the face"
            annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

          annotation (defaultComponentName="thermal");
        end PartialBC;

        type BCType = enumeration(
            Temperature "Temperature",
            HeatFlowRate "Heat flow rate") "Types of BCs";
      end BaseClasses;
    end Thermal;

    package BaseClasses "Base classes (not for direct use)"
      extends Modelica.Icons.BasesPackage;

      model PartialSpecies
        "<html>Partial BC for a face of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"
        extends FCSys.BaseClasses.Icons.BCs.Single;

        // Material
        parameter Boolean isobaric=false "Isobaric condition" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(group="Normal linear momentum",compact=true));
        replaceable Material.Density normal(spec(k(start=U.atm/(298.15*U.K))))
          if not isobaric constrainedby Material.BaseClasses.PartialBC
          "Type of condition" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="Normal linear momentum",
            enable=not isobaric,
            __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{-70,-14},{-50,6}})));
        // Note:  In Dymola 7.4, the value of k must be specified here instead
        // of at the lower level (e.g., Material.Density) so that the spec
        // subcomponent can be replaced by blocks that don't contain the
        // parameter k.

        // Heat
        replaceable Thermal.Temperature thermal(spec(k(start=298.15*U.K)))
          constrainedby Thermal.BaseClasses.PartialBC "Type of condition"
          annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Heat", __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{50,-14},{70,6}})));

        FCSys.Connectors.RealInputBus u "Input bus for external signal sources"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));

      equation
        // Material
        connect(u.normal, material.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,20},{-60,20},{-60,6.66134e-16}},

            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        // Heat
        connect(u.thermal, thermal.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,20},{60,20},{60,6.66134e-16}},

            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

      end PartialSpecies;

      partial model PartialBC "Partial model for a BC"
        extends FCSys.BaseClasses.Icons.BCs.Single;

        parameter Boolean internal=true "Use internal specification"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(__Dymola_descriptionLabel=true));

        replaceable Modelica.Blocks.Sources.Constant spec if internal
          constrainedby Modelica.Blocks.Interfaces.SO "Internal specification"
          annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(__Dymola_descriptionLabel=true, enable=internal),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-30,30})));

        FCSys.Connectors.RealInput u if not internal "Value of BC" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));

      protected
        FCSys.Connectors.RealInputInternal u_final "Final value of BC"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-10})));

      equation
        connect(u, u_final) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,27.5},{5.55112e-16,27.5},{
                5.55112e-16,15},{5.55112e-16,-10},{5.55112e-16,-10}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(spec.y, u_final) annotation (Line(
            points={{-30,19},{-30,10},{5.55112e-16,10},{5.55112e-16,-10}},
            color={0,0,127},
            smooth=Smooth.None));
      end PartialBC;
    end BaseClasses;
  end Face;

  package FaceBusDifferential
    "<html>BCs for a pair of <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connectors</html>"
    extends Modelica.Icons.Package;

    model Subregion
      "<html>BC for faces of a <a href=\"modelica://FCSys.Regions.Region\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with efforts by default</html>"
      extends FCSys.BaseClasses.Icons.BCs.Double;

      parameter Axis axis=Axis.x "Axis normal to the face";

      replaceable Phases.Gas gas(final axis=axis) "Gas" annotation (Dialog(
            group="Phases", __Dymola_descriptionLabel=true), Placement(
            transformation(extent={{-10,-10},{10,10}})));

      replaceable Phases.Graphite graphite(final axis=axis) "Graphite"
        annotation (Dialog(group="Phases", __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{-10,-10},{10,10}})));

      replaceable Phases.Ionomer ionomer(final axis=axis) "Ionomer" annotation
        (Dialog(group="Phases", __Dymola_descriptionLabel=true), Placement(
            transformation(extent={{-10,-10},{10,10}})));

      FCSys.Connectors.FaceBus negative annotation (Placement(transformation(
              extent={{-110,-10},{-90,10}}), iconTransformation(extent={{-110,-10},
                {-90,10}})));
      FCSys.Connectors.RealInputBus u annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40})));
      FCSys.Connectors.FaceBus positive annotation (Placement(transformation(
              extent={{90,-10},{110,10}}), iconTransformation(extent={{90,-10},
                {110,10}})));
    equation
      connect(gas.negative, negative.gas) annotation (Line(
          points={{-10,6.10623e-16},{-100,0},{-100,5.55112e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(graphite.negative, negative.graphite) annotation (Line(
          points={{-10,6.10623e-16},{-100,6.10623e-16},{-100,5.55112e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(ionomer.negative, negative.ionomer) annotation (Line(
          points={{-10,6.10623e-16},{-100,5.55112e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));
      connect(gas.positive, positive.gas) annotation (Line(
          points={{10,6.10623e-16},{20,0},{100,5.55112e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(graphite.positive, positive.graphite) annotation (Line(
          points={{10,6.10623e-16},{20,6.10623e-16},{20,5.55112e-16},{100,
              5.55112e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(ionomer.positive, positive.ionomer) annotation (Line(
          points={{10,6.10623e-16},{100,5.55112e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));
      connect(u.gas, gas.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,4},{6.10623e-16,4}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(u.graphite, graphite.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,4},{6.10623e-16,4}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(u.ionomer, ionomer.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,4},{6.10623e-16,4}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="subregionFaceBC",
        Diagram(graphics));
    end Subregion;

    model SubregionFlow
      "<html>BC for faces of a <a href=\"modelica://FCSys.Regions.Region\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with flows by default</html>"
      extends Subregion(
        gas(
          H2(
            redeclare replaceable FaceDifferential.Material.Current material,
            redeclare replaceable FaceDifferential.Mechanical.Force transverseX,

            redeclare replaceable FaceDifferential.Mechanical.Force transverseY,

            redeclare replaceable FaceDifferential.Mechanical.Force transverseZ,

            redeclare replaceable FaceDifferential.Thermal.HeatFlowRate thermal),

          H2O(
            redeclare replaceable FaceDifferential.Material.Current material,
            redeclare replaceable FaceDifferential.Mechanical.Force transverseX,

            redeclare replaceable FaceDifferential.Mechanical.Force transverseY,

            redeclare replaceable FaceDifferential.Mechanical.Force transverseZ,

            redeclare replaceable FaceDifferential.Thermal.HeatFlowRate thermal),

          N2(
            redeclare replaceable FaceDifferential.Material.Current material,
            redeclare replaceable FaceDifferential.Mechanical.Force transverseX,

            redeclare replaceable FaceDifferential.Mechanical.Force transverseY,

            redeclare replaceable FaceDifferential.Mechanical.Force transverseZ,

            redeclare replaceable FaceDifferential.Thermal.HeatFlowRate thermal),

          O2(
            redeclare replaceable FaceDifferential.Material.Current material,
            redeclare replaceable FaceDifferential.Mechanical.Force transverseX,

            redeclare replaceable FaceDifferential.Mechanical.Force transverseY,

            redeclare replaceable FaceDifferential.Mechanical.Force transverseZ,

            redeclare replaceable FaceDifferential.Thermal.HeatFlowRate thermal)),

        graphite('C+'(
            redeclare replaceable FaceDifferential.Material.Current material,
            redeclare replaceable FaceDifferential.Mechanical.Force transverseX,

            redeclare replaceable FaceDifferential.Mechanical.Force transverseY,

            redeclare replaceable FaceDifferential.Mechanical.Force transverseZ,

            redeclare replaceable FaceDifferential.Thermal.HeatFlowRate thermal),
            'e-'(
            redeclare replaceable FaceDifferential.Material.Current material,
            redeclare replaceable FaceDifferential.Mechanical.Force transverseX,

            redeclare replaceable FaceDifferential.Mechanical.Force transverseY,

            redeclare replaceable FaceDifferential.Mechanical.Force transverseZ,

            redeclare replaceable FaceDifferential.Thermal.HeatFlowRate thermal)),

        ionomer(
          redeclare replaceable FaceDifferential.Mechanical.Force transverseX,
          redeclare replaceable FaceDifferential.Mechanical.Force transverseY,
          redeclare replaceable FaceDifferential.Mechanical.Force transverseZ,
          redeclare replaceable FaceDifferential.Thermal.HeatFlowRate thermal),

        H2O(
          redeclare replaceable FaceDifferential.Material.Current material,
          redeclare replaceable FaceDifferential.Mechanical.Force transverseX,
          redeclare replaceable FaceDifferential.Mechanical.Force transverseY,
          redeclare replaceable FaceDifferential.Mechanical.Force transverseZ,
          redeclare replaceable FaceDifferential.Thermal.HeatFlowRate thermal),

        'H+'(
          redeclare replaceable FaceDifferential.Material.Current material,
          redeclare replaceable FaceDifferential.Mechanical.Force transverseX,
          redeclare replaceable FaceDifferential.Mechanical.Force transverseY,
          redeclare replaceable FaceDifferential.Mechanical.Force transverseZ,
          redeclare replaceable FaceDifferential.Thermal.HeatFlowRate thermal));

      annotation (defaultComponentPrefixes="replaceable",defaultComponentName=
            "subregionFaceBC");
    end SubregionFlow;

    model SubregionClosed
      "<html>BC for faces of a <a href=\"modelica://FCSys.Regions.Region\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with efforts except closed by default</html>"
      extends Subregion(
        gas(
          H2(isobaric=true),
          H2O(isobaric=true),
          N2(isobaric=true),
          O2(isobaric=true)),
        graphite('C+'(isobaric=true), 'e-'(thermoOpt=ThermoOpt.ClosedDiabatic)),

        ionomer(Error),
        H2O(isobaric=true),
        'H+'(isobaric=true));

      annotation (defaultComponentPrefixes="replaceable",defaultComponentName=
            "subregionFaceBC");
    end SubregionClosed;

    model SubregionClosedAdiabatic
      "<html>BC for faces of a <a href=\"modelica://FCSys.Regions.Region\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with flows except zero shear velocities by default</html>"
      extends Subregion(
        gas(
          H2(isobaric=true),
          H2O(isobaric=true),
          N2(isobaric=true),
          O2(isobaric=true)),
        graphite('C+'(isobaric=true), 'e-'(thermoOpt=ThermoOpt.ClosedAdiabatic)),

        ionomer(Error),
        H2O(isobaric=true),
        'H+'(isobaric=true));

      annotation (defaultComponentPrefixes="replaceable",defaultComponentName=
            "subregionFaceBC");
    end SubregionClosedAdiabatic;

    package Phases
      "<html>BCs for the <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Phase\">Phase</a> model (multi-species)</html>"
      extends Modelica.Icons.Package;

      model Gas "BC for gas"

        extends BaseClasses.NullPhase;

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
        FCSys.BCs.FaceDifferential.Species H2 if inclH2 "Model" annotation (
            Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2), Placement(transformation(extent={{-10,-10},{10,10}})));

        parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FaceDifferential.Species H2O if inclH2O "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,10}})));

        parameter Boolean inclN2=false "<html>Nitrogen (N<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FaceDifferential.Species N2 if inclN2 "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclN2), Placement(transformation(extent={{-10,-10},{10,10}})));

        parameter Boolean inclO2=false "<html>Oxygen (O<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FaceDifferential.Species O2 if inclO2 "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclO2), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // H2
        connect(H2.negative, negative.H2) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive, positive.H2) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.H2, H2.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2O
        connect(H2O.negative, negative.H2O) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive, positive.H2O) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.H2O, H2O.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // N2
        connect(N2.negative, negative.N2) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive, positive.N2) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.N2, N2.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // O2
        connect(O2.negative, negative.O2) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive, positive.O2) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.O2, O2.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
      end Gas;

      model Graphite "BC for graphite"

        extends BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean 'inclC+'=false
          "<html>Carbon plus (C<sup>+</sup>)</html>" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FaceDifferential.Species 'C+'(isobaric=true) if '' + ERROR annotation (
            Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='' + ERROR annotation (
              Evaluate=true,
              HideResult=true,
              choices(__Dymola_checkBox=true),
              Dialog(
                group="Species",
                __Dymola_descriptionLabel=true,
                __Dymola_joinNext=true))));
        FaceDifferential.Species 'e-' if 'incle-' "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='incle-'), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // C+
        connect('C+'.negative, negative.'C+') annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('C+'.positive, positive.'C+') annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'C+', 'C+'.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // e-
        connect('e-'.negative, negative.'e-') annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive, positive.'e-') annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'e-', 'e-'.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

      end Graphite;

      model Ionomer "BC for ionomer"

        extends BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean 'inclC19HF37O5S-'=false
          "<html>Nafion sulfonate minus (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FaceDifferential.Species 'C19HF37O5S-'(isobaric=true) if
          'inclC19HF37O5S-' "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclC19HF37O5S-'), Placement(transformation(extent={{-10,-10},
                  {10,10}})));

        parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FaceDifferential.Species H2O if inclH2O "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,10}})));

        parameter Boolean 'inclH+'=false "<html>Protons (H<sup>+</sup>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FaceDifferential.Species 'H+' if 'inclH+' "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclH+'), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // C19HF37O5S-
        connect('C19HF37O5S-'.negative, negative.'C19HF37O5S-') annotation (
            Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('C19HF37O5S-'.positive, positive.'C19HF37O5S-') annotation (
            Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'C19HF37O5S-', 'C19HF37O5S-'.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2O
        connect(H2O.negative, negative.H2O) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive, positive.H2O) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.H2O, H2O.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H+
        connect('H+'.negative, negative.'H+') annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive, positive.'H+') annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'H+', 'H+'.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

      end Ionomer;

      model Liquid "BC for liquid"

        extends BaseClasses.NullPhase;

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
        FaceDifferential.Species H2O if inclH2O "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // H2O
        connect(H2O.negative, negative.H2O) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive, positive.H2O) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.H2O, H2O.u) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,14},{6.10623e-16,14},{
                6.10623e-16,4}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
      end Liquid;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        model NullPhase "Empty BC for a phase (no species)"
          extends FCSys.BaseClasses.Icons.BCs.Double;

          parameter Axis axis=Axis.x "Axis normal to the face";

          FCSys.Connectors.FaceBus negative
            "Multi-species connector for linear momentum and heat" annotation (
              Placement(transformation(extent={{-110,-10},{-90,10}})));
          FCSys.Connectors.FaceBus positive
            "Multi-species connector for linear momentum and heat"
            annotation (Placement(transformation(extent={{90,-10},{110,10}})));

          FCSys.Connectors.RealInputBus u "Bus of inputs to specify conditions"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,40})));
        end NullPhase;
      end BaseClasses;
    end Phases;

    annotation (Documentation(info="<html><p>The hierarchy of these
 boundary condition models is similar to that of the models in the
 <a href=\"modelica://FCSys.BCs.FaceBus\">BCs.FaceBus</a> package.
 For more information, please see the documentation in that package.</p></html>"));
  end FaceBusDifferential;

  package FaceDifferential
    "<html>BCs for a pair of <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connectors</html>"

    extends Modelica.Icons.Package;

    model Species
      "<html>BC for a pair of faces of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"
      extends BaseClasses.PartialSpecies;

      parameter Axis axis "Axis normal to the face";

      // X-axis linear momentum
      parameter Boolean inviscidX=false "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="X-axis linear momentum",
          enable=axis <> 1,
          __Dymola_descriptionLabel=true));
      // Note:  Dymola 7.4 doesn't recognize enumerations in the dialog enable
      // option, e.g.,
      //     enable=axis <> Axis.x.
      // Therefore, the values of the enumerations are specified numerically.
      replaceable Mechanical.Velocity transverseX(spec(k=0)) if axis <> Axis.x
         and not inviscidX constrainedby Mechanical.BaseClasses.PartialBC
        "Type of condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="X-axis linear momentum",
          enable=not inviscidX and axis <> 1,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-40,-10},{-20,10}})));

      // Y-axis linear momentum
      parameter Boolean inviscidY=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="Y-axis linear momentum",
          enable=axis <> 2,
          __Dymola_descriptionLabel=true));
      replaceable Mechanical.Velocity transverseY(spec(k=0)) if axis <> Axis.y
         and not inviscidY constrainedby Mechanical.BaseClasses.PartialBC
        "Type of condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Y-axis linear momentum",
          enable=not inviscidY and axis <> 2,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-10,-20},{10,0}})));

      // Z-axis linear momentum
      parameter Boolean inviscidZ=false "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="Z-axis linear momentum",
          enable=axis <> 3,
          __Dymola_descriptionLabel=true));

      replaceable Mechanical.Velocity transverseZ(spec(k=0)) if axis <> Axis.z
         and not inviscidZ constrainedby Mechanical.BaseClasses.PartialBC
        "Type of condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Z-axis linear momentum",
          enable=not inviscidZ and axis <> 2,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{20,-30},{40,-10}})));

      FCSys.Connectors.Face negative(
        final axis=axis,
        final isobaric=isobaric,
        final inviscidX=inviscidX,
        final inviscidY=inviscidY,
        final inviscidZ=inviscidZ) "Negative face" annotation (Placement(
            transformation(extent={{-110,-40},{-90,-20}}), iconTransformation(
              extent={{-110,-10},{-90,10}})));
      FCSys.Connectors.Face positive(
        final axis=axis,
        final inviscidX=inviscidX,
        final inviscidY=inviscidY,
        final inviscidZ=inviscidZ) "Positive face" annotation (Placement(
            transformation(extent={{90,-40},{110,-20}}), iconTransformation(
              extent={{90,-10},{110,10}})));

    equation
      // Material
      connect(material.negative, negative.normal) annotation (Line(
          points={{-70,10},{-80,10},{-80,-30},{-100,-30}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(material.positive, positive.normal) annotation (Line(
          points={{-50,10},{90,10},{90,-30},{100,-30}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));

      // X-axis linear momentum
      connect(transverseX.negative, negative.transverseX) annotation (Line(
          points={{-40,6.10623e-16},{-80,6.10623e-16},{-80,-30},{-100,-30}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(transverseX.positive, positive.transverseX) annotation (Line(
          points={{-20,6.10623e-16},{90,6.10623e-16},{90,-30},{100,-30}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(u.transverseX, transverseX.u) annotation (Line(
          points={{5.55112e-16,40},{0,40},{0,20},{-30,20},{-30,4}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // Y-axis linear momentum
      connect(transverseY.negative, negative.transverseY) annotation (Line(
          points={{-10,-10},{-80,-10},{-80,-30},{-100,-30}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(transverseY.positive, positive.transverseY) annotation (Line(
          points={{10,-10},{90,-10},{90,-30},{100,-30}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(u.transverseY, transverseY.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,10},{6.10623e-16,10},{
              6.10623e-16,-6}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // Z-axis linear momentum
      connect(transverseZ.negative, negative.transverseZ) annotation (Line(
          points={{20,-20},{-80,-20},{-80,-30},{-100,-30}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(transverseZ.positive, positive.transverseZ) annotation (Line(
          points={{40,-20},{90,-20},{90,-30},{100,-30}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(u.transverseZ, transverseZ.u) annotation (Line(
          points={{5.55112e-16,40},{0,40},{0,20},{30,20},{30,-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // Heat
      connect(thermal.negative, negative.thermal) annotation (Line(
          points={{50,-30},{-26,-30},{-26,-30},{-100,-30}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(thermal.positive, positive.thermal) annotation (Line(
          points={{70,-30},{86,-30},{86,-30},{100,-30}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));

      annotation (Icon(graphics));
    end Species;

    package Material "Material BCs"
      extends Modelica.Icons.Package;

      model Density "Prescribed density difference, with material conservation"

        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Density,
            u(final unit="N/L3"));

      equation
        negative.rho - positive.rho = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="material");
      end Density;

      model Current "Prescribed current, with material conservation"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Current,
            u(final unit="N/T"));

      equation
        negative.Ndot = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="material");
      end Current;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialBC "Partial model for a material BC"
          extends FaceDifferential.BaseClasses.PartialBC;

          constant BCType bCType "Type of BC";
          // Note:  This is included so that the type of BC is recorded with the
          // results.

          FCSys.Connectors.BaseClasses.Normal negative
            "Material connector for the negative face" annotation (Placement(
                transformation(extent={{-110,-10},{-90,10}})));
          FCSys.Connectors.BaseClasses.Normal positive
            "Material connector for the positive face"
            annotation (Placement(transformation(extent={{90,-10},{110,10}})));
        equation
          0 = negative.Ndot + positive.Ndot
            "Material conservation (no storage)";
          annotation (defaultComponentName="material");
        end PartialBC;

        type BCType = enumeration(
            Density "Density difference",
            Current "Current") "Types of BCs";
      end BaseClasses;
    end Material;

    package Mechanical "Mechanical BCs"
      extends Modelica.Icons.Package;

      model Velocity
        "Prescribed shear velocity, with conservation of linear momentum"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Velocity,
            u(final unit="l/T"));
      equation
        negative.phi - positive.phi = u_final;

        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="mechanicalBC");
      end Velocity;

      model Force
        "Prescribed shear force, with conservation of linear momentum"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Force, u(
              final unit="l.m/T2"));

      equation
        negative.mPhidot = u_final;

        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="mechanicalBC");
      end Force;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialBC "Partial model for a BC for linear momentum"
          extends FaceDifferential.BaseClasses.PartialBC;

          constant BCType bCType "Type of BC";
          // Note:  This is included so that the type of BC is recorded with the
          // results.

          FCSys.Connectors.BaseClasses.Transverse negative
            "Linear momentum connector for the negative face" annotation (
              Placement(transformation(extent={{-110,-10},{-90,10}})));
          FCSys.Connectors.BaseClasses.Transverse positive
            "Linear momentum connector for the positive face"
            annotation (Placement(transformation(extent={{90,-10},{110,10}})));

        equation
          0 = negative.mPhidot + positive.mPhidot
            "Conservation of linear momentum (no storage)";
          annotation (defaultComponentName="mechanicalBC");
        end PartialBC;

        type BCType = enumeration(
            Velocity "Shear velocity",
            Force "Shear force") "Types of BCs";
      end BaseClasses;
    end Mechanical;

    package Thermal "Thermal BCs"
      extends Modelica.Icons.Package;

      model Temperature
        "Prescribed temperature difference, with energy conservation"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.Temperature,
            redeclare FCSys.Connectors.RealInput u(final unit="l2.m/(N.T2)",
              displayUnit="K"));

      equation
        negative.T - positive.T = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="thermal");
      end Temperature;

      model HeatRate "Prescribed heat flow rate, with energy conservation"
        extends BaseClasses.PartialBC(final bCType=BaseClasses.BCType.HeatRate,
            redeclare FCSys.Connectors.RealInput u(final unit="l2.m/T3"));

      equation
        negative.Qdot = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="thermal");
      end HeatRate;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialBC "Partial model for a BC for heat"
          extends FaceDifferential.BaseClasses.PartialBC;

          constant BCType bCType "Type of BC";
          // Note:  This is included so that the type of BC is recorded with the
          // results.

          FCSys.Connectors.Thermal negative
            "Heat connector for the negative face" annotation (Placement(
                transformation(extent={{-110,-10},{-90,10}})));
          FCSys.Connectors.Thermal positive
            "Heat connector for the positive face"
            annotation (Placement(transformation(extent={{90,-10},{110,10}})));

        equation
          0 = negative.Qdot + positive.Qdot "Energy conservation (no storage)";
          annotation (defaultComponentName="thermal");
        end PartialBC;

        type BCType = enumeration(
            Temperature "Temperature difference",
            HeatRate "Heat flow rate") "Types of BCs";
      end BaseClasses;
    end Thermal;

    package BaseClasses "Base classes (not for direct use)"
      extends Modelica.Icons.BasesPackage;

      partial model PartialSpecies
        "<html>Partial BC for a pair of faces of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"
        extends FCSys.BaseClasses.Icons.BCs.Double;

        // Material
        parameter Boolean isobaric=false "Isobaric condition" annotation (
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(group="Normal linear momentum",compact=true));
        replaceable Material.Density normal(spec(k(start=4*U.C/U.cc))) if not
          isobaric constrainedby Material.BaseClasses.PartialBC
          "Type of condition" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(
            group="Normal linear momentum",
            enable=not isobaric,
            __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{-70,0},{-50,20}})));
        // Note:  In Dymola 7.4, the value of k must be specified here instead
        // of at the lower level (e.g., Material.Density) so that the spec
        // subcomponent can be replaced by blocks that don't contain the
        // parameter k.

        // Heat
        replaceable Thermal.Temperature thermal(spec(k(start=298.15*U.K)))
          constrainedby Thermal.BaseClasses.PartialBC "Type of condition"
          annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Heat", __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{50,-40},{70,-20}})));

        FCSys.Connectors.RealInputBus u "Input bus for external signal sources"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));

      equation
        // Material
        connect(u.normal, material.u) annotation (Line(
            points={{5.55112e-16,40},{0,40},{0,20},{-60,20},{-60,14}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        // Heat
        connect(u.thermal, thermal.u) annotation (Line(
            points={{5.55112e-16,40},{0,40},{0,20},{60,20},{60,-26}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        annotation (Diagram(graphics));
      end PartialSpecies;

      partial model PartialBC "Partial model for a BC"
        extends FCSys.BaseClasses.Icons.BCs.Double;

        parameter Boolean internal=true "Use internal specification"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(__Dymola_descriptionLabel=true));

        replaceable Modelica.Blocks.Sources.Constant spec(k=0) if internal
          constrainedby Modelica.Blocks.Interfaces.SO "Internal specification"
          annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(__Dymola_descriptionLabel=true, enable=internal),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-30,30})));

        FCSys.Connectors.RealInput u if not internal "Value of BC" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));

      protected
        FCSys.Connectors.RealInputInternal u_final "Final value of BC"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-10})));

      equation
        connect(u, u_final) annotation (Line(
            points={{5.55112e-16,40},{5.55112e-16,27.5},{5.55112e-16,27.5},{
                5.55112e-16,15},{5.55112e-16,-10},{5.55112e-16,-10}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(spec.y, u_final) annotation (Line(
            points={{-30,19},{-30,10},{5.55112e-16,10},{5.55112e-16,-10}},
            color={0,0,127},
            smooth=Smooth.None));
      end PartialBC;
    end BaseClasses;
  end FaceDifferential;

  model Router "Connect two pairs of faces to pass through or cross over"
    extends FCSys.BaseClasses.Icons.Names.Top3;
    parameter Boolean crossOver=false "Cross over (otherwise, pass through)"
      annotation (Evaluate=true,choices(__Dymola_checkBox=true));
    FCSys.Connectors.FaceBus negative1 "Negative face 1" annotation (Placement(
          transformation(extent={{-90,-50},{-70,-30}}, rotation=0),
          iconTransformation(extent={{-90,-50},{-70,-30}})));
    FCSys.Connectors.FaceBus positive1 "Positive face 1" annotation (Placement(
          transformation(extent={{70,-50},{90,-30}}, rotation=0),
          iconTransformation(extent={{70,-50},{90,-30}})));
    FCSys.Connectors.FaceBus negative2 "Negative face 2" annotation (Placement(
          transformation(extent={{-90,30},{-70,50}}, rotation=0),
          iconTransformation(extent={{-90,30},{-70,50}})));
    FCSys.Connectors.FaceBus positive2 "Positive face 2" annotation (Placement(
          transformation(extent={{70,30},{90,50}}, rotation=0),
          iconTransformation(extent={{70,30},{90,50}})));
  equation
    if crossOver then
      connect(negative1, positive2) annotation (Line(
          points={{-80,-40},{-40,-40},{0,0},{40,40},{80,40}},
          color={127,127,127},
          smooth=Smooth.Bezier,
          thickness=0.5,
          pattern=LinePattern.Dash));
      connect(negative2, positive1) annotation (Line(
          points={{-80,40},{-80,40},{-40,40},{0,0},{40,-40},{80,-40}},
          color={127,127,127},
          smooth=Smooth.Bezier,
          thickness=0.5,
          pattern=LinePattern.Dash));
    else
      // Pass-through
      connect(negative1, positive1) annotation (Line(
          points={{-80,-40},{-80,-40},{80,-40}},
          color={127,127,127},
          thickness=0.5));
      connect(negative2, positive2) annotation (Line(
          points={{-80,40},{-80,40},{80,40}},
          color={127,127,127},
          thickness=0.5));
    end if;
    annotation (
      Documentation(info="<html>
<p>This model acts as a connection switch.
It has a single parameter, <code>crossOver</code>.  If <code>crossOver</code> is
set to <code>false</code>, then
the router is in the pass-through mode.
In pass-through mode, shown by Figure 1a,
<code>negative1</code> is connected to <code>positive1</code> and <code>negative2</code>
is connected to <code>positive2</code>.
If <code>crossOver</code> is set to true, then the router is in cross-over mode.
In cross-over mode, shown by Figure 1b, <code>negative1</code> is connected to <code>positive2</code>
and <code>negative2</code> is
connected to <code>positive1</code>. The only equations are
those generated by the model's <code>connect</code> statements.</p>

    <table border=\"0\" cellspacing=\"0\" cellpadding=\"2\" align=center>
      <tr align=center>
        <td align=center width=120>
          <img src=\"modelica://FCSys/resources/documentation/BCs/Router/PassThrough.png\">
<br><b>a:</b>  Pass-through
        </td>
        <td align=center>
          <img src=\"modelica://FCSys/resources/documentation/BCs/Router/CrossOver.png\">
<br><b>b:</b>  Cross-over
        </td>
      </tr>
      <tr align=center>
        <td colspan=2 align=center><b>Figure 1:</b> Modes of connection.</td>
      </tr>
    </table>
</html>"),
      Icon(graphics={Line(
              points={{-80,40},{-40,40},{0,0},{40,-40},{80,-40}},
              color={127,127,127},
              thickness=0.5,
              visible=crossOver,
              smooth=Smooth.Bezier),Line(
              points={{-80,40},{80,40}},
              color={127,127,127},
              visible=not crossOver,
              smooth=Smooth.None,
              thickness=0.5),Line(
              points={{-80,-40},{80,-40}},
              color={127,127,127},
              visible=not crossOver,
              smooth=Smooth.None,
              thickness=0.5),Line(
              points={{-80,-40},{-40,-40},{0,0},{40,40},{80,40}},
              color={127,127,127},
              thickness=0.5,
              visible=crossOver,
              smooth=Smooth.Bezier)}),
      Diagram(graphics));
  end Router;

  record Environment "Environmental properties for a model"
    extends FCSys.BaseClasses.Icons.Names.Top3;

    // Store the values of the base constants and units.
    final constant U.Bases.Base base=U.base "Base constants and units";

    parameter Boolean analysis=true "Include optional variables for analysis"
      annotation (choices(__Dymola_checkBox=true));

    parameter Q.PressureAbsolute p(nominal=U.atm) = 1*U.atm "Pressure";
    parameter Q.TemperatureAbsolute T(nominal=298.15*U.K) = 298.15*U.K
      "Temperature";
    parameter Q.NumberAbsolute RH(displayUnit="%") = 1 "Relative humidity";
    parameter Q.NumberAbsolute x_O2_dry(
      final max=1,
      displayUnit="%") = 0.208
      "<html>Dry gas O<sub>2</sub> fraction (<i>y</i><sub>O2 dry</sub>)</html>";
    // Value from http://en.wikipedia.org/wiki/Oxygen
    parameter Q.Acceleration a[FCSys.BaseClasses.Axis]=zeros(3)
      "Acceleration of the reference frame";

    final parameter Q.NumberAbsolute x_H2O(
      final max=1,
      displayUnit="%") = 0.2
      "<html>Gas H<sub>2</sub>O fraction (<i>y</i><sub>H2O</sub>)</html>";
    // TODO:  Cast this in terms of relative humidity.

    annotation (
      defaultComponentPrefixes="inner",
      missingInnerMessage="Your model is using an outer \"environment\" record, but an inner \"environment\" record is not defined.
For simulation, specify global default settings by dragging FCSys.BCs.Environment into your model.
The default global default settings will be used for the current simulation.",
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Rectangle(
              extent={{-120,60},{120,100}},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None),Text(
              extent={{-120,60},{120,100}},
              textString="%name",
              lineColor={0,0,0}),Rectangle(
              extent={{-80,60},{80,-100}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.Dash),Rectangle(
              extent={{-70,50},{70,-98}},
              lineColor={255,255,255},
              fillPattern=FillPattern.HorizontalCylinder,
              fillColor={170,170,255}),Rectangle(
              extent={{-72,-60},{72,-98}},
              fillPattern=FillPattern.Solid,
              fillColor={255,255,255},
              pattern=LinePattern.None,
              lineColor={0,0,0}),Line(points={{-70,-60},{70,-60}}, color={0,0,0}),
            Line(points={{-40,-20},{-10,-50},{40,0}}, color={0,0,0}),Ellipse(
              extent={{30,10},{50,-10}},
              pattern=LinePattern.None,
              lineColor={255,255,255},
              fillColor={240,0,0},
              fillPattern=FillPattern.Sphere),Line(points={{-66,-90},{-36,-60}},
            color={0,0,0}),Line(points={{2,-90},{32,-60}}, color={0,0,0}),Line(
            points={{36,-90},{66,-60}}, color={0,0,0}),Line(points={{-32,-90},{
            -2,-60}}, color={0,0,0}),Rectangle(
              extent={{70,50},{76,-60}},
              fillPattern=FillPattern.Solid,
              fillColor={255,255,255},
              pattern=LinePattern.None,
              lineColor={0,0,0}),Rectangle(
              extent={{-76,50},{-70,-60}},
              fillPattern=FillPattern.Solid,
              fillColor={255,255,255},
              pattern=LinePattern.None,
              lineColor={0,0,0})}),
      Diagram(graphics));
  end Environment;

  package BaseClasses "Base classes (not for direct use)"
    extends Modelica.Icons.BasesPackage;

    block RealFunction
      "<html>Set an output signal according to a <code>Real</code> function of an input</html>"

      extends FCSys.BaseClasses.Icons.Names.Top2;
      //extends FCSys.BaseClasses.Icons.Blocks.ContinuousShort;

      FCSys.Connectors.RealInput u "Value of Real input" annotation (Placement(
            transformation(extent={{-110,-10},{-90,10}}), iconTransformation(
              extent={{-110,-10},{-90,10}})));
      FCSys.Connectors.RealOutput y=0.0*u
        "<html>Value of <code>Real</code> output</html>" annotation (Dialog(
            group="Time varying output signal"), Placement(transformation(
              extent={{90,-10},{110,10}}, rotation=0), iconTransformation(
              extent={{90,-10},{110,10}})));

      annotation (
        Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            grid={2,2}),graphics={Rectangle(
                  extent={{-100,40},{100,-40}},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  lineColor={0,0,0}),Text(
                  extent={{-100,-10},{100,10}},
                  lineColor={127,127,127},
                  textString="%y")}),
        Diagram(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            grid={2,2}),graphics),
        Documentation(info="<html>
<p>
The (time-varying) <code>Real</code> output signal of this block can be defined in its
parameter menu via variable <i>y</i>. The purpose is to support the
easy definition of <code>Real</code> expressions in a block diagram. For example,
in the y-menu the definition <code>if time &gt; 1 then 1 else 0</code> can be given in order
to produce an output signal that transitions from zero to one at time one.  Note that
<code>time</code> is a built-in variable that is always
accessible and represents the \"model time.\"
Variable <i>y</i> is both a variable and a connector.
Variable <i>u</i> is too, and it may be used in the expression for <i>y</i>.
</p>
<p>This block has been adapted from <a href=\"modelica://Modelica.Blocks.Sources.RealExpression\">Modelica.Blocks.Sources.RealExpression</a>.</p>
</html>"));
    end RealFunction;
  end BaseClasses;

end BCs;
