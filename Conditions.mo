within FCSys;
package Conditions "Models to impose and measure operating conditions"
  extends Modelica.Icons.Package;

  package Examples "Examples"
    extends Modelica.Icons.ExamplesPackage;

    model FaceCondition
      "<html>Test the Conditions for the face of a subregion</html>"
      extends Modelica.Icons.Example;

      FaceBus.Subregion subregionFaceCondition(gas(inclH2O=true, H2O(redeclare
              Face.Normal.CurrentAreic normal(spec(k=0)))))
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
        gas(inclH2O=true))
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      inner Conditions.Environment environment
        annotation (Placement(transformation(extent={{30,30},{50,50}})));

    equation
      connect(subregion.yPositive, subregionFaceCondition.face) annotation (
          Line(
          points={{6.10623e-16,10},{-3.36456e-22,16},{6.10623e-16,16},{
              6.10623e-16,20}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (experiment(NumberOfIntervals=5000), Commands(file=
              "resources/scripts/Dymola/Conditions.Examples.FaceCondition.mos"));
    end FaceCondition;

    model FaceConditionPhases
      "<html>Test the Conditions for the face of a subregion with phases</html>"
      extends Modelica.Icons.Example;

      // Geometric parameters
      inner parameter Q.Length L[Axis](each min=Modelica.Constants.small,start=
            ones(3)*U.cm) "<html>Length (<b>L</b>)</html>"
        annotation (Dialog(group="Geometry"));
      final inner parameter Q.Area A[Axis]={L[cartWrap(axis + 1)]*L[cartWrap(
          axis + 2)] for axis in Axis} "Cross-sectional area";

      FaceBus.Phases.Gas phaseFaceCondition(inclH2O=true, H2O(redeclare
            Face.Normal.CurrentAreic normal(spec(k=0))))
        annotation (Placement(transformation(extent={{-10,14},{10,34}})));

      Subregions.Volume volume
        annotation (Placement(transformation(extent={{-16,-16},{16,16}})));
      FCSys.Subregions.Phases.Gas gas(
        inclReact=false,
        inclLin={false,true,false},
        inclH2=false,
        inclH2O=true)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      inner Conditions.Environment environment
        annotation (Placement(transformation(extent={{30,30},{50,50}})));

    equation
      connect(gas.yPositive, phaseFaceCondition.face) annotation (Line(
          points={{6.10623e-16,10},{-3.36456e-22,16},{6.10623e-16,16},{
              6.10623e-16,20}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(gas.inert, volume.inert) annotation (Line(
          points={{8,-8},{11,-11}},
          color={72,90,180},
          smooth=Smooth.None));
      annotation (experiment);
    end FaceConditionPhases;

    model Router "<html>Test the <code>Router<code> model</html>"
      extends Modelica.Icons.Example;

      // TODO:  Make this into a meaningful example.
      Conditions.Router router
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

    end Router;

    model AnodeAdapter "<html>Test the <code>'Adapte-'</code> model</html>"

      extends Modelica.Icons.Example;
      extends Modelica.Icons.UnderConstruction;
      // **fails check
      inner Modelica.Fluid.System system(T_ambient=293.15 + 5)
        annotation (Placement(transformation(extent={{40,70},{60,90}})));
      inner Conditions.Environment environment(T=350*U.K)
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
      annotation (experiment(StopTime=2e-10), Commands(file=
              "resources/scripts/Dymola/Conditions.Examples.Adapteminus.mos"));
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
                {70,10},{90,30}}),iconTransformation(extent={{70,10},{90,30}})));
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
                {70,10},{90,30}}),iconTransformation(extent={{70,10},{90,30}})));
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
        connect(H2.face, face.H2) annotation (Line(
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
        connect(H2O.face, face.H2O) annotation (Line(
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
        connect(H2O.face, face.H2O) annotation (Line(
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
        connect(N2.face, face.N2) annotation (Line(
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
        connect(O2.face, face.O2) annotation (Line(
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
        connect('C+'.face, face.'C+') annotation (Line(
            points={{-8,20},{-40,20},{-40,5.55112e-16},{-80,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect('C+'.heatPort, heatPort) annotation (Line(
            points={{8,20},{40,20},{40,5.55112e-16},{80,5.55112e-16}},
            color={191,0,0},
            smooth=Smooth.None));

        // e-
        connect('e-'.face, face.'e-') annotation (Line(
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
        connect(H2O.face, face.H2) annotation (Line(
            points={{-8,6.10623e-16},{-54,6.10623e-16},{-54,5.55112e-16},{-80,
                5.55112e-16}},
            color={0,0,0},
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

      package BaseClasses "Base classes (not generally for direct use)"
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
                  thickness=0.5), Line(
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
            FCSys.Characteristics.'e-'.Graphite Data);

        parameter Q.Area A=U.cm^2 "Area of the interface";
        parameter Side side=Side.n
          "Side of the interface w.r.t. this component";

        Modelica.Electrical.Analog.Interfaces.NegativePin pin
          "Modelica electrical pin" annotation (Placement(transformation(extent
                ={{70,30},{90,50}}), iconTransformation(extent={{70,30},{90,50}})));

      equation
        // Efforts
        Data.g(face.T, inSign(side)*face.mPhidot_0/A) = Data.z*pin.v*U.V
          "Electrical potential";

        // Conservation (no storage)
        0 = A*face.J + pin.i*U.A/Data.z "Material";
        annotation (Documentation(info="<html><p>For additional information, see the
    <a href=\"modelica://FCSys.Conditions.Adapters.Species.BaseClasses.PartialSpecies\">
    PartialSpecies</a> model.</p>
    </html>"), Icon(graphics={Line(
                      points={{0,40},{80,40}},
                      color={0,0,255},
                      smooth=Smooth.None),Line(
                      points={{0,60},{0,-20}},
                      color={0,0,0},
                      smooth=Smooth.None,
                      pattern=LinePattern.Dash)}));
      end 'e-';

      model FluidNonionic
        "<html>Adapter to connect a single nonionic fluid species between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"

        extends BaseClasses.PartialSpecies;

        parameter Q.Area A=U.cm^2 "Area of the interface";
        parameter Side side=Side.n
          "Side of the interface w.r.t. this component";

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
        face.mPhidot_0 = inSign(side)*A*fluidPort.p*U.Pa;
        medium.h = fluidPort.h_outflow;

        // Conservation (no storage)
        0 = face.J*A + (fluidPort.m_flow/medium.MM)*U.mol/U.s "Material";

        // See the partial model for additional equations.
        annotation (Documentation(info="<html><p>The electrical connector (<code>pin</code>) is only included
    if the species is ionic.
    </p>
    <p>For additional information, see the
    <a href=\"modelica://FCSys.Conditions.Adapters.Species.BaseClasses.PartialSpecies\">
    PartialSpecies</a> model.</p>
    </html>"), Icon(graphics={Line(
                points={{0,-40},{80,-40}},
                color={0,127,255},
                smooth=Smooth.None), Line(
                points={{0,20},{0,-60}},
                color={0,0,0},
                smooth=Smooth.None,
                pattern=LinePattern.Dash)}));
      end FluidNonionic;

      model Solid
        "<html>Adapter to connect a single solid species between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a> (heat only)</html>"

        extends BaseClasses.PartialSpecies;

      equation
        face.J = 0 "Closed";
        annotation (Documentation(info="<html><p>For additional information, see the
    <a href=\"modelica://FCSys.Conditions.Adapters.Species.BaseClasses.PartialSpecies\">
    PartialSpecies</a> model.</p>
    </html>"), Icon(graphics={Line(
                      points={{0,20},{0,-20}},
                      color={0,0,0},
                      smooth=Smooth.None,
                      pattern=LinePattern.Dash)}));
      end Solid;

      package BaseClasses "Base classes (not generally for direct use)"
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

          FCSys.Connectors.Face face
            "Connector for linear momentum and heat of a single species"
            annotation (Placement(transformation(extent={{-90,-10},{-70,10}}),
                iconTransformation(extent={{-90,-10},{-70,10}})));

          Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPort
            "Modelica heat port" annotation (Placement(transformation(extent={{
                    70,-10},{90,10}}), iconTransformation(extent={{70,-10},{90,
                    10}})));

        equation
          // Efforts
          face.T = heatPort.T*U.K "Temperature";

          // Conservation (no storage)
          face.mPhidot = {0,0} "Transverse linear momentum";
          0 = face.Qdot + heatPort.Q_flow*U.W "Energy";
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
                  smooth=Smooth.None), Line(
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

      package BaseClasses "Base classes (not generally for direct use)"
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
      FCSys.Conditions.BaseClasses.RealFunction unit4(y=unit4.u*U.kPa + U.atm)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-94,-30})));
      FCSys.Conditions.BaseClasses.RealFunction unit5(y=unit5.u*U.kPa + U.atm)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-74,0})));
      FCSys.Conditions.BaseClasses.RealFunction unit6(y=unit6.u*U.kPa + U.atm)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-54,-30})));
      FCSys.Conditions.BaseClasses.RealFunction unit7(y=unit7.u*U.kPa + U.atm)
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
      FCSys.Conditions.BaseClasses.RealFunction unit10(y=(unit10.u + 273.15)*U.K)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={32,-30})));
      FCSys.Conditions.BaseClasses.RealFunction unit11(y=(unit11.u + 273.15)*U.K)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={52,0})));
      FCSys.Conditions.BaseClasses.RealFunction unit12(y=(unit12.u + 273.15)*U.K)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={72,-30})));
      FCSys.Conditions.BaseClasses.RealFunction unit13(y=(unit13.u + 273.15)*U.K)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={92,0})));
      FCSys.Conditions.BaseClasses.RealFunction unit14(y=(unit14.u + 273.15)*U.K)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={112,-30})));
      FCSys.Conditions.BaseClasses.RealFunction unit15(y=(unit15.u + 273.15)*U.K)
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
      annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-180,
                -100},{180,100}}), graphics), experiment(StopTime=15481,
            Algorithm="Euler"));
    end Replay;

    package BaseClasses "Base classes (not generally for direct use)"
      extends Modelica.Icons.BasesPackage;
      partial model PartialTestStand "Partial cell test stand"
        extends FCSys.BaseClasses.Icons.Names.Top9;
        extends Modelica.Icons.UnderConstruction;
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

        FaceBus.SubregionFlow anEndCondition[n_y, n_z](each graphite('inclC+'=
                true, 'incle-'=true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-136,0})));
        FaceBus.SubregionFlow caEndCondition[n_y, n_z](each graphite('inclC+'=
                true, 'incle-'=true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={136,0})));
        FaceBus.SubregionFlow anSourceCondition[n_x_an, n_z](each gas(inclH2=
                true, inclH2O=true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-40,-136})));
        FaceBus.SubregionFlow anSinkCondition[n_x_an, n_z](each gas(inclH2=true,
              inclH2O=true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-40,136})));
        FaceBus.SubregionFlow caSourceCondition[n_x_ca, n_z](each gas(
            inclH2O=true,
            inclN2=true,
            inclO2=true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={40,-136})));
        FaceBus.SubregionFlow caSinkCondition[n_x_ca, n_z](each gas(
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
        replaceable FCSys.Conditions.FaceBusPair.Subregion current[n_y, n_z](
            graphite('inclC+'=true, 'incle-'=true)) if inclIO constrainedby
          FCSys.Conditions.FaceBusPair.Subregion(graphite('inclC+'=true,
              'incle-'=true))
          annotation (Placement(transformation(extent={{-140,20},{-120,40}})));

      equation
        connect(anSourceCondition.face, anSource) annotation (Line(
            points={{-40,-140},{-40,-160}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));

        connect(anSinkCondition.face, anSink) annotation (Line(
            points={{-40,140},{-40,160}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));

        connect(caSourceCondition.face, caSource) annotation (Line(
            points={{40,-140},{40,-160}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));

        connect(caSinkCondition.face, caSink) annotation (Line(
            points={{40,140},{40,160}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));

        connect(caEndCondition.face, caEnd) annotation (Line(
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
        connect(voltage.negative, anEndCondition.face) annotation (Line(
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

        connect(anEndCondition.face, anEnd) annotation (Line(
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

        FaceBus.SubregionFlow anEnd[n_y, n_z](each graphite(
            'inclC+'=true,
            'incle-'=true,
            'e-'(redeclare Face.Normal.CurrentAreic normal(redeclare
                  Modelica.Blocks.Sources.Ramp spec(height=U.A/U.cm^2,duration=
                      50))))) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-30,0})));
        FaceBus.SubregionFlow caEnd[n_y, n_z](each graphite(
            'inclC+'=true,
            'incle-'=true,
            'e-'(redeclare Face.Normal.CurrentAreic normal(redeclare
                  Modelica.Blocks.Sources.Ramp spec(height=U.A/U.cm^2,duration=
                      50))))) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={30,0})));
        FaceBus.SubregionFlow anSource[n_x_an, n_z](each gas(inclH2=true,
              inclH2O=true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={-20,-30})));
        FaceBus.SubregionFlow anSink[n_x_an, n_z](each gas(inclH2=true, inclH2O
              =true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-20,30})));
        FaceBus.SubregionFlow caSource[n_x_ca, n_z](each gas(
            inclH2O=true,
            inclN2=true,
            inclO2=true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=180,
              origin={20,-30})));
        FaceBus.SubregionFlow caSink[n_x_ca, n_z](each gas(
            inclH2O=true,
            inclN2=true,
            inclO2=true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={20,30})));

        inner Environment environment
          annotation (Placement(transformation(extent={{50,20},{70,40}})));
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
    "<html>Conditions for the <a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a> connector</html>"
    extends Modelica.Icons.Package;

    model Gas "Condition for gas"

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

    model Graphite "Condition for graphite"

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
      Chemical.Species 'C+'(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ,
        redeclare FCSys.Characteristics.'C+'.Graphite Data) if 'inclC+' "Model"
        annotation (Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable='inclC+'), Placement(transformation(extent={{-10,-10},{10,10}})));

      parameter Boolean 'incle-'=false "<html>Electrons (e<sup>-</sup>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
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

    model Ionomer "Condition for ionomer"

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

    model Liquid "Condition for liquid"

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

    package BaseClasses "Base classes (not generally for direct use)"
      extends Modelica.Icons.BasesPackage;

      model NullPhase "Empty condition for a phase (no species)"
        extends FCSys.Conditions.BaseClasses.Icons.Single;

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
    "<html>Conditions for the <a href=\"modelica://FCSys.Connectors.ChemicalOutput\">ChemicalOutput</a> connector</html>"
    extends Modelica.Icons.Package;

    model Species
      "<html>Conditions for the <a href=\"modelica://FCSys.Connectors.ChemicalOutput\">ChemicalOutput</a> connector</html>"

      import FCSys.Conditions.Chemical.BaseClasses.ConditionTypeMaterial;
      import FCSys.Conditions.Chemical.BaseClasses.ConditionTypeMechanical;
      import FCSys.Conditions.Chemical.BaseClasses.ConditionTypeFluid;
      extends FCSys.Conditions.BaseClasses.Icons.Single;

      replaceable package Data =
          FCSys.Characteristics.BaseClasses.Characteristic constrainedby
        FCSys.Characteristics.BaseClasses.Characteristic "Characteristic data"
        annotation (
        __Dymola_choicesAllMatching=true,
        Dialog(group="Material properties"),
        Placement(transformation(extent={{-60,40},{-40,60}}),
            iconTransformation(extent={{-10,90},{10,110}})));

      // Material
      parameter ConditionTypeMaterial material=ConditionTypeMaterial.PotentialPerTemperature
        "Type of condition"
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
      parameter ConditionTypeMechanical xCondition=ConditionTypeMechanical.Velocity
        "Type of condition" annotation (Dialog(
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
      replaceable Modelica.Blocks.Sources.Constant xSpec(k(start=0)) if
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
      parameter ConditionTypeMechanical yCondition=ConditionTypeMechanical.Velocity
        "Type of condition" annotation (Dialog(
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
      replaceable Modelica.Blocks.Sources.Constant ySpec(k(start=0)) if
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
      parameter ConditionTypeMechanical zCondition=ConditionTypeMechanical.Velocity
        "Type of condition" annotation (Dialog(
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
      replaceable Modelica.Blocks.Sources.Constant zSpec(k(start=0)) if
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
      parameter ConditionTypeMechanical fluidCondition=ConditionTypeFluid.EnthalpyMassic
        "Type of condition"
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

      FCSys.Connectors.RealInputInternal u_material(final unit=if material ==
            ConditionTypeMaterial.PotentialPerTemperature then "1" else "N/T")
        if not internalMaterial "Material signal" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-70,30}),iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40})));
      FCSys.Connectors.RealInputInternal u_x(final unit="l/T") if not
        internalLinX and inclLinX
        "Signal for the x component of linear momentum" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-30,30}),iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40})));
      FCSys.Connectors.RealInputInternal u_y(final unit="l/T") if not
        internalLinY and inclLinY
        "Signal for the y component of linear momentum" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={10,30}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40})));
      FCSys.Connectors.RealInputInternal u_z(final unit="l/T") if not
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
      FCSys.Connectors.RealInputInternal u_int[n_lin]
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
      if material == ConditionTypeMaterial.PotentialPerTemperature then
        chemical.muPerT = u_material_int;
      else
        chemical.Ndot = u_material_int;
      end if;

      // X component of linear momentum
      if inclLinX then
        //  if xCondition == ConditionTypeMechanical.Velocity then
        chemical.phi[1] = u_int[1];
        //  end if;
      end if;

      // Y component of linear momentum
      if inclLinY then
        //  if yCondition == ConditionTypeMechanical.Velocity then
        chemical.phi[2] = u_int[2];
        //  end if;
      end if;

      // Z component of linear momentum
      if inclLinZ then
        //  if zCondition == ConditionTypeMechanical.Velocity then
        chemical.phi[3] = u_int[3];
        //  end if;
      end if;

      // Enthalpy
      //  if fluidCondition == ConditionTypeFluid.EnthalpyMassic then
      chemical.hbar = u_fluid_int;
      //  end if;

      connect(u.material, u_material) annotation (Line(
          points={{5.55112e-16,60},{5.55112e-16,40},{-70,40},{-70,30}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(u.x, u_x) annotation (Line(
          points={{5.55112e-16,60},{5.55112e-16,40},{-30,40},{-30,30}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(u.y, u_y) annotation (Line(
          points={{5.55112e-16,60},{5.55112e-16,40},{10,40},{10,30}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(u.z, u_z) annotation (Line(
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
      connect(u_x, u_int[1]) annotation (Line(
          points={{-30,30},{-30,-8},{10,-8},{10,-20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u_y, u_int[2]) annotation (Line(
          points={{10,30},{10,-20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u_z, u_int[3]) annotation (Line(
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
      connect(xSpec.y, u_int[1]) annotation (Line(
          points={{-50,-1},{-50,-8},{10,-8},{10,-20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(ySpec.y, u_int[2]) annotation (Line(
          points={{30,-1},{30,-8},{10,-8},{10,-20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(zSpec.y, u_int[3]) annotation (Line(
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

    package Material "Material Conditions"
      extends Modelica.Icons.Package;

      model Current "Impose current"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.CurrentAreic,
            u(final unit="N/(l2.T)"));

      equation
        face.J = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="normal");
      end Current;

      model Force "Impose normal force"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Force,
            u(final unit="l.m/T2"));

      equation
        face.mPhidot_0 = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="normal");
      end Force;

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition "Partial model for a normal condition"
          extends Face.BaseClasses.PartialCondition;

          constant ConditionType conditionType "Type of condition";
          // Note:  This is included so that the type of condition is recorded with the
          // results.

        equation
          face.mPhidot = {0,0} "No transverse forces";
          face.Qdot = 0 "Adiabatic";
          annotation (defaultComponentName="normal");
        end PartialCondition;

        type ConditionType = enumeration(
            CurrentAreic "Areic current",
            Force "Force") "Types of Conditions";

      end BaseClasses;

    end Material;

    package BaseClasses "Base classes (not generally for direct use)"
      extends Modelica.Icons.BasesPackage;
      type ConditionTypeMaterial = enumeration(
          PotentialPerTemperature
            "Impose quotient of electrochemical potential and temperature",
          Current "Impose current") "Types of material Conditions";

      type ConditionTypeMechanical = enumeration(
          Velocity "Impose velocity",
          Force "Impose force") "Types of mechanical Conditions";
      type ConditionTypeFluid = enumeration(
          EnthalpyMassic "Impose massic enthalpy",
          EnthalpyRate "Impose enthalpy flow rate") "Types of fluid Conditions";

    end BaseClasses;

  end Chemical;

  package Inert
    "<html>Conditions for the <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> or <a href=\"modelica://FCSys.Connectors.InertInternal\">InertInternal</a> connector</html>"
    extends Modelica.Icons.Package;

  end Inert;

  package InertAmagat
    "<html>Conditions for the <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector</html>"
    extends Modelica.Icons.Package;

    model Phase
      "<html>Condition for the <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Phases\">Phase</a> model</html>"
      extends FCSys.Conditions.BaseClasses.Icons.Single;

      // Volume
      replaceable Volume.Volume volumeCondition(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ) constrainedby
        Volume.BaseClasses.PartialCondition "Condition" annotation (
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
      replaceable Mechanical.Force linXCondition(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ) if inclLinX constrainedby
        Mechanical.BaseClasses.PartialCondition "Condition" annotation (
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
      replaceable Mechanical.Force linYCondition(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ) if inclLinY constrainedby
        Mechanical.BaseClasses.PartialCondition "Condition" annotation (
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
      replaceable Mechanical.Force linZCondition(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ) if inclLinZ constrainedby
        Mechanical.BaseClasses.PartialCondition "Condition" annotation (
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
        final inclLinZ=inclLinZ) constrainedby
        Thermal.BaseClasses.PartialCondition "Condition" annotation (
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
      connect(volumeCondition.inert, inert) annotation (Line(
          points={{-80,-8},{-80,-20},{5.55112e-16,-20},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(u.volume, volumeCondition.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,20},{-80,20},{-80,6.66134e-16}},

          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // X component of linear momentum
      connect(linXCondition.inert, inert) annotation (Line(
          points={{-40,-8},{-40,-20},{5.55112e-16,-20},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(u.linX, linXCondition.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,20},{-40,20},{-40,6.66134e-16}},

          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // Y component of linear momentum
      connect(linYCondition.inert, inert) annotation (Line(
          points={{6.10623e-16,-8},{6.10623e-16,-20},{5.55112e-16,-20},{
              5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(u.linY, linYCondition.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,6.66134e-16},{6.10623e-16,
              6.66134e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // Z component of linear momentum
      connect(linZCondition.inert, inert) annotation (Line(
          points={{40,-8},{40,-20},{5.55112e-16,-20},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(u.linZ, linZCondition.u) annotation (Line(
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

    package Volume "Conditions for additivity of volume"
      extends Modelica.Icons.Package;

      model Pressure "Impose pressure"
        extends BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.Pressure,
          u(final unit="m/(l.T2)"),
          spec(k(start=U.atm)));

      equation
        inert.p = u_final;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="volumeCondition");
      end Pressure;

      model Volume "Impose volume"
        extends BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.Volume,
          u(final unit="l3"),
          spec(k(start=U.cc)));

      equation
        inert.V = u_final;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="volumeCondition");
      end Volume;

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition
          "Partial model of a condition for volume"
          extends FCSys.Conditions.InertAmagat.BaseClasses.PartialCondition;

          constant ConditionType conditionType "Type of condition";
          // Note:  This is included so that the type of condition is recorded with the
          // results.

        equation
          inert.mPhidot = zeros(n_lin) "No force";
          inert.Qdot = 0 "Adiabatic";
          annotation (defaultComponentName="volumeCondition");
        end PartialCondition;

        type ConditionType = enumeration(
            Volume "Volume",
            Pressure "Pressure") "Types of Conditions";

      end BaseClasses;

    end Volume;

    package Mechanical "Mechanical Conditions"
      extends Modelica.Icons.Package;
      model Velocity "Impose velocity"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Velocity,
            u(final unit="l/T"));

      equation
        inert.phi[linAxes[axis]] = u_final;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="mechanicalCondition");
      end Velocity;

      model Force "Impose force"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Force,
            u(final unit="l.m/T2"));

      equation
        inert.mPhidot[linAxes[axis]] = u_final;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="mechanicalCondition");
      end Force;

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition
          "Partial model for a mechanical condition"
          extends FCSys.Conditions.InertAmagat.BaseClasses.PartialCondition;

          parameter Axis axis=Axis.x "Axis" annotation (HideResult=true);

          constant ConditionType conditionType "Type of condition";
          // Note:  This is included so that the type of condition is recorded with the
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
          annotation (defaultComponentName="mechanicalCondition");
        end PartialCondition;

        type ConditionType = enumeration(
            Velocity "Velocity",
            Force "Force") "Types of Conditions";

      end BaseClasses;

    end Mechanical;

    package Thermal "Thermal conditions"
      extends Modelica.Icons.Package;

      model Temperature "Impose temperature"
        extends BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.Temperature,
          u(final unit="l2.m/(N.T2)", displayUnit="K"),
          spec(k(start=298.15*U.K)));

      equation
        inert.T = u_final;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="thermal");
      end Temperature;

      model HeatFlowRate "Impose heat flow rate"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.HeatFlowRate,
            u(final unit="l2.m/T3"));

      equation
        inert.Qdot = u_final;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="thermal");
      end HeatFlowRate;

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition "Partial model for a thermal condition"
          extends FCSys.Conditions.InertAmagat.BaseClasses.PartialCondition;

          constant ConditionType conditionType "Type of condition";
          // Note:  This is included so that the type of condition is recorded with the
          // results.

        equation
          inert.V = 0 "No volume";
          inert.mPhidot = zeros(n_lin) "No force";
          annotation (defaultComponentName="thermal");
        end PartialCondition;

        type ConditionType = enumeration(
            Temperature "Temperature",
            HeatFlowRate "Heat flow rate") "Types of Conditions";

      end BaseClasses;

    end Thermal;

    package BaseClasses "Base classes (not generally for direct use)"
      extends Modelica.Icons.BasesPackage;
      partial model PartialCondition "Partial model of a condition"
        extends FCSys.Conditions.BaseClasses.Icons.Single;

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

        FCSys.Connectors.RealInput u if not internal "Value of condition"
          annotation (Placement(transformation(
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

        FCSys.Connectors.RealInputInternal u_final "Final value of condition"
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

      end PartialCondition;

    end BaseClasses;

  end InertAmagat;

  package InertDalton
    "<html>Conditions for the <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector</html>"
    extends Modelica.Icons.Package;

    model Species
      "<html>Condition for the <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model</html>"

      extends FCSys.Conditions.BaseClasses.Icons.Single;

      // Pressure
      replaceable Pressure.Volume pressureCondition(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ) constrainedby
        Pressure.BaseClasses.PartialCondition "Condition" annotation (
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
      replaceable Mechanical.Force linXCondition(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ) if inclLinX constrainedby
        Mechanical.BaseClasses.PartialCondition "Condition" annotation (
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
      replaceable Mechanical.Force linYCondition(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ) if inclLinY constrainedby
        Mechanical.BaseClasses.PartialCondition "Condition" annotation (
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
      replaceable Mechanical.Force linZCondition(
        final inclLinX=inclLinX,
        final inclLinY=inclLinY,
        final inclLinZ=inclLinZ) if inclLinZ constrainedby
        Mechanical.BaseClasses.PartialCondition "Condition" annotation (
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
        final inclLinZ=inclLinZ) constrainedby
        Thermal.BaseClasses.PartialCondition "Condition" annotation (
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
      connect(pressureCondition.inert, inert) annotation (Line(
          points={{-80,-8},{-80,-20},{5.55112e-16,-20},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(u.pressure, pressureCondition.u) annotation (Line(
          points={{5.55112e-16,40},{0,40},{0,20},{-80,20},{-80,6.66134e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // X component of linear momentum
      connect(linXCondition.inert, inert) annotation (Line(
          points={{-40,-8},{-40,-20},{5.55112e-16,-20},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(u.linX, linXCondition.u) annotation (Line(
          points={{5.55112e-16,40},{0,40},{0,20},{-40,20},{-40,6.66134e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // Y component of linear momentum
      connect(linYCondition.inert, inert) annotation (Line(
          points={{6.10623e-16,-8},{6.10623e-16,-20},{5.55112e-16,-20},{
              5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(u.linY, linYCondition.u) annotation (Line(
          points={{5.55112e-16,40},{0,40},{0,6.66134e-16},{6.10623e-16,
              6.66134e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      // Z component of linear momentum
      connect(linZCondition.inert, inert) annotation (Line(
          points={{40,-8},{40,-20},{5.55112e-16,-20},{5.55112e-16,-40}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(u.linZ, linZCondition.u) annotation (Line(
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

    package Pressure "Conditions for additivity of pressure"
      extends Modelica.Icons.Package;

      model Volume "Impose volume"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Volume,
            u(final unit="l3"));

      equation
        inert.V = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="pressureCondition");
      end Volume;

      model Pressure "Impose pressure"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Pressure,
            u(final unit="m/(l.T2)"));

      equation
        inert.p = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="pressureCondition");
      end Pressure;

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition
          "Partial model of a condition for pressure"
          extends FCSys.Conditions.InertDalton.BaseClasses.PartialCondition;

          constant ConditionType conditionType "Type of condition";
          // Note:  This is included so that the type of condition is recorded with the
          // results.

        equation
          inert.mPhidot = zeros(n_lin) "No force";
          inert.Qdot = 0 "No heat flow";
          annotation (defaultComponentName="pressureCondition");
        end PartialCondition;

        type ConditionType = enumeration(
            Volume "Volume",
            Pressure "Pressure") "Types of Conditions";

      end BaseClasses;

    end Pressure;

    package Mechanical "Mechanical Conditions"
      extends Modelica.Icons.Package;
      model Velocity "Impose velocity"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Velocity,
            u(final unit="l/T"));

      equation
        inert.phi[linAxes[axis]] = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="mechanicalCondition");
      end Velocity;

      model Force "Impose force"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Force,
            u(final unit="l.m/T2"));

      equation
        inert.mPhidot[linAxes[axis]] = u_final;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="mechanicalCondition");
      end Force;

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition
          "Partial model for a mechanical condition"
          extends FCSys.Conditions.InertDalton.BaseClasses.PartialCondition;

          parameter Axis axis=Axis.x "Axis" annotation (HideResult=true);

          constant ConditionType conditionType "Type of condition";
          // Note:  This is included so that the type of condition is recorded with the
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
          annotation (defaultComponentName="mechanicalCondition");
        end PartialCondition;

        type ConditionType = enumeration(
            Velocity "Velocity",
            Force "Force") "Types of Conditions";

      end BaseClasses;

    end Mechanical;

    package Thermal "Conditions for heat"
      extends Modelica.Icons.Package;

      model Temperature "Impose temperature"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Temperature,
            u(final unit="l2.m/(N.T2)", displayUnit="K"));

      equation
        inert.T = u_final;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="thermal");
      end Temperature;

      model HeatFlowRate "Impose heat flow rate"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.HeatFlowRate,
            u(final unit="l2.m/T3"));

      equation
        inert.Qdot = u_final;
        annotation (defaultComponentPrefixes="replaceable",defaultComponentName
            ="thermal");
      end HeatFlowRate;

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition "Partial model for a thermal condition"
          extends FCSys.Conditions.InertDalton.BaseClasses.PartialCondition;

          constant ConditionType conditionType "Type of condition";
          // Note:  This is included so that the type of condition is recorded with the
          // results.

        equation
          inert.p = 0 "No pressure";
          inert.mPhidot = zeros(n_lin) "No force";
          annotation (defaultComponentName="thermal");
        end PartialCondition;

        type ConditionType = enumeration(
            Temperature "Temperature",
            HeatFlowRate "Heat flow rate") "Types of Conditions";

      end BaseClasses;

    end Thermal;

    package BaseClasses "Base classes (not generally for direct use)"
      extends Modelica.Icons.BasesPackage;
      partial model PartialCondition "Partial model of a condition"
        extends FCSys.Conditions.BaseClasses.Icons.Single;

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

        FCSys.Connectors.RealInput u if not internal "Value of condition"
          annotation (Placement(transformation(
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

        FCSys.Connectors.RealInputInternal u_final "Final value of condition"
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

      end PartialCondition;

    end BaseClasses;

  end InertDalton;

  package FaceBus
    "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, e.g., of a <a href=\"modelica://FCSys.Regions.Region\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model</html>"

    extends Modelica.Icons.Package;

    model Subregion
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, with efforts by default</html>"

      extends FCSys.Conditions.BaseClasses.Icons.Single;

      Phases.Gas gas "Gas" annotation (Dialog(group="Phases",
            __Dymola_descriptionLabel=true), Placement(transformation(extent={{
                -10,-10},{10,10}})));

      Phases.Graphite graphite "Graphite" annotation (Dialog(group="Phases",
            __Dymola_descriptionLabel=true), Placement(transformation(extent={{
                -10,-10},{10,10}})));

      Phases.Ionomer ionomer "Ionomer" annotation (Dialog(group="Phases",
            __Dymola_descriptionLabel=true), Placement(transformation(extent={{
                -10,-10},{10,10}})));

      Phases.Liquid liquid "Liquid" annotation (Dialog(group="Phases",
            __Dymola_descriptionLabel=true), Placement(transformation(extent={{
                -10,-10},{10,10}})));

      FCSys.Connectors.FaceBus face
        "Connector for linear momentum and heat of multiple species"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}}),
            iconTransformation(extent={{-10,-50},{10,-30}})));
      FCSys.Connectors.RealInputBus u "Bus of inputs to specify conditions"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-100,0})));

      FCSys.Connectors.RealOutputBus y "Output bus of measurements" annotation
        (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={100,0})));

    equation
      // Gas
      connect(gas.face, face.gas) annotation (Line(
          points={{6.10623e-16,-4},{8.60423e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(u.gas, gas.u) annotation (Line(
          points={{-100,5.55112e-16},{-10,5.55112e-16},{-10,6.10623e-16}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(gas.y, y.gas) annotation (Line(
          points={{10,6.10623e-16},{100,6.10623e-16},{100,5.55112e-16}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      // Graphite
      connect(graphite.face, face.graphite) annotation (Line(
          points={{6.10623e-16,-4},{8.60423e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(u.graphite, graphite.u) annotation (Line(
          points={{-100,5.55112e-16},{-10,5.55112e-16},{-10,6.10623e-16}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(graphite.y, y.graphite) annotation (Line(
          points={{10,6.10623e-16},{100,6.10623e-16},{100,5.55112e-16}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      // Ionomer
      connect(ionomer.face, face.ionomer) annotation (Line(
          points={{6.10623e-16,-4},{-4.87687e-22,-4},{-4.87687e-22,-40},{
              5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));
      connect(u.ionomer, ionomer.u) annotation (Line(
          points={{-100,5.55112e-16},{-10,5.55112e-16},{-10,6.10623e-16}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(ionomer.y, y.ionomer) annotation (Line(
          points={{10,6.10623e-16},{100,6.10623e-16},{100,5.55112e-16}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      // Liquid
      connect(liquid.face, face.liquid) annotation (Line(
          points={{6.10623e-16,-4},{-4.87687e-22,-4},{-4.87687e-22,-40},{
              5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));
      connect(u.liquid, liquid.u) annotation (Line(
          points={{-100,5.55112e-16},{-10,5.55112e-16},{-10,6.10623e-16}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(liquid.y, y.liquid) annotation (Line(
          points={{10,6.10623e-16},{100,6.10623e-16},{100,5.55112e-16}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (defaultComponentPrefixes="replaceable", defaultComponentName=
            "subregion");
    end Subregion;

    model SubregionFlow
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, with flows by default</html>"

      extends FaceBus.Subregion(
        gas(
          H2(
            redeclare replaceable Face.Normal.Force normal,
            redeclare replaceable Face.Transverse.Force transverse1,
            redeclare replaceable Face.Transverse.Force transverse2,
            redeclare replaceable Face.Thermal.HeatRate thermal(source(k=0))),
          H2O(
            redeclare replaceable Face.Normal.Force normal,
            redeclare replaceable Face.Transverse.Force transverse1,
            redeclare replaceable Face.Transverse.Force transverse2,
            redeclare replaceable Face.Thermal.HeatRate thermal(source(k=0))),
          N2(
            redeclare replaceable Face.Normal.Force normal,
            redeclare replaceable Face.Transverse.Force transverse1,
            redeclare replaceable Face.Transverse.Force transverse2,
            redeclare replaceable Face.Thermal.HeatRate thermal(source(k=0))),
          O2(
            redeclare replaceable Face.Normal.Force normal,
            redeclare replaceable Face.Transverse.Force transverse1,
            redeclare replaceable Face.Transverse.Force transverse2,
            redeclare replaceable Face.Thermal.HeatRate thermal(source(k=0)))),

        graphite('C+'(
            redeclare replaceable Face.Normal.Force normal,
            redeclare replaceable Face.Transverse.Force transverse1,
            redeclare replaceable Face.Transverse.Force transverse2,
            redeclare replaceable Face.Thermal.HeatRate thermal(source(k=0))),
            'e-'(
            redeclare replaceable Face.Normal.Force normal,
            redeclare replaceable Face.Transverse.Force transverse1,
            redeclare replaceable Face.Transverse.Force transverse2,
            redeclare replaceable Face.Thermal.HeatRate thermal(source(k=0)))),

        ionomer(
          'C19HF37O5S-'(
            redeclare replaceable Face.Normal.Force normal,
            redeclare replaceable Face.Transverse.Force transverse1,
            redeclare replaceable Face.Transverse.Force transverse2,
            redeclare replaceable Face.Thermal.HeatRate thermal(source(k=0))),
          'H+'(
            redeclare replaceable Face.Normal.Force normal,
            redeclare replaceable Face.Transverse.Force transverse1,
            redeclare replaceable Face.Transverse.Force transverse2,
            redeclare replaceable Face.Thermal.HeatRate thermal(source(k=0))),
          H2O(
            redeclare replaceable Face.Normal.Force normal,
            redeclare replaceable Face.Transverse.Force transverse1,
            redeclare replaceable Face.Transverse.Force transverse2,
            redeclare replaceable Face.Thermal.HeatRate thermal(source(k=0)))),

        liquid(H2O(
            redeclare replaceable Face.Normal.Force normal,
            redeclare replaceable Face.Transverse.Force transverse1,
            redeclare replaceable Face.Transverse.Force transverse2,
            redeclare replaceable Face.Thermal.HeatRate thermal(source(k=0)))));

      annotation (defaultComponentPrefixes="replaceable",defaultComponentName=
            "subregion");

    end SubregionFlow;

    package Phases
      "<html>Conditions for the <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Phase\">Phase</a> model (multi-species)</html>"
      extends Modelica.Icons.Package;

      model Gas "Condition for gas"

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
        Face.Species H2 if inclH2 "Conditions" annotation (Dialog(
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
        Face.Species H2O if inclH2O "Conditions" annotation (Dialog(
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
        Face.Species N2 if inclN2 "Conditions" annotation (Dialog(
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
        Face.Species O2 if inclO2 "Conditions" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclO2), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // H2
        connect(H2.face, face.H2) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.H2, H2.u) annotation (Line(
            points={{-100,5.55112e-16},{-100,0},{-10,0},{-10,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(H2.y, y.H2) annotation (Line(
            points={{10,6.10623e-16},{10,0},{100,0},{100,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2O
        connect(H2O.face, face.H2O) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.H2O, H2O.u) annotation (Line(
            points={{-100,5.55112e-16},{-100,0},{-10,0},{-10,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(H2O.y, y.H2O) annotation (Line(
            points={{10,6.10623e-16},{10,0},{100,0},{100,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // N2
        connect(N2.face, face.N2) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.N2, N2.u) annotation (Line(
            points={{-100,5.55112e-16},{-100,0},{-10,0},{-10,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(N2.y, y.N2) annotation (Line(
            points={{10,6.10623e-16},{10,0},{100,0},{100,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // O2
        connect(O2.face, face.O2) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.O2, O2.u) annotation (Line(
            points={{-100,5.55112e-16},{-100,0},{-10,0},{-10,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(O2.y, y.O2) annotation (Line(
            points={{10,6.10623e-16},{10,0},{100,0},{100,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (Diagram(graphics));
      end Gas;

      model Graphite "Condition for graphite"

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
        Face.Species 'C+' if 'inclC+' "Conditions" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclC+'), Placement(transformation(extent={{-10,-10},{10,10}})));

        parameter Boolean 'incle-'=false
          "<html>Electrons (e<sup>-</sup>)</html>" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species 'e-' if 'incle-' "Conditions" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='incle-'), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // C+
        connect('C+'.face, face.'C+') annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'C+', 'C+'.u) annotation (Line(
            points={{-100,5.55112e-16},{-100,0},{-10,0},{-10,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect('C+'.y, y.'C+') annotation (Line(
            points={{10,6.10623e-16},{10,0},{100,0},{100,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // e-
        connect('e-'.face, face.'e-') annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'e-', 'e-'.u) annotation (Line(
            points={{-100,5.55112e-16},{-100,0},{-10,0},{-10,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect('e-'.y, y.'e-') annotation (Line(
            points={{10,6.10623e-16},{10,0},{100,0},{100,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

      end Graphite;

      model Ionomer "Condition for ionomer"

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
        Face.Species 'C19HF37O5S-' if 'inclC19HF37O5S-' "Conditions"
          annotation (Dialog(
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
        Face.Species 'H+' if 'inclH+' "Conditions" annotation (Dialog(
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
        Face.Species H2O if inclH2O "Conditions" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // C19HF37O5S-
        connect('C19HF37O5S-'.face, face.'C19HF37O5S-') annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'C19HF37O5S-', 'C19HF37O5S-'.u) annotation (Line(
            points={{-100,5.55112e-16},{-100,0},{-10,0},{-10,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect('C19HF37O5S-'.y, y.'C19HF37O5S-') annotation (Line(
            points={{10,6.10623e-16},{10,0},{100,0},{100,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H+
        connect('H+'.face, face.'H+') annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'H+', 'H+'.u) annotation (Line(
            points={{-100,5.55112e-16},{-100,0},{-10,0},{-10,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect('H+'.y, y.'H+') annotation (Line(
            points={{10,6.10623e-16},{10,0},{100,0},{100,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2O
        connect(H2O.face, face.H2O) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.H2O, H2O.u) annotation (Line(
            points={{-100,5.55112e-16},{-100,0},{-10,0},{-10,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(H2O.y, y.H2O) annotation (Line(
            points={{10,6.10623e-16},{10,0},{100,0},{100,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

      end Ionomer;

      model Liquid "Condition for liquid"

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
        Face.Species H2O if inclH2O "Conditions" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // H2O
        connect(H2O.face, face.H2O) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.H2O, H2O.u) annotation (Line(
            points={{-100,5.55112e-16},{-100,0},{-10,0},{-10,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(H2O.y, y.H2O) annotation (Line(
            points={{10,6.10623e-16},{10,0},{100,0},{100,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

      end Liquid;

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;
        model NullPhase "Empty condition for a phase (no species)"
          extends FCSys.Conditions.BaseClasses.Icons.Single;

          FCSys.Connectors.FaceBus face
            "Multi-species connector for linear momentum and heat"
            annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
          FCSys.Connectors.RealInputBus u
            "Input bus for values of imposed conditions" annotation (Placement(
                transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-100,0}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-100,0})));

          FCSys.Connectors.RealOutputBus y "Output bus of measurements"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={100,0}),iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={100,0})));

        end NullPhase;

      end BaseClasses;

    end Phases;

  end FaceBus;

  package Face
    "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"
    model Species
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connector</html>"

      extends FCSys.Conditions.BaseClasses.Icons.Single;

      replaceable Normal.CurrentAreic normal(source(k=0)) constrainedby
        Normal.BaseClasses.PartialCondition "Normal condition" annotation (
          __Dymola_choicesFromPackage=true, Placement(transformation(extent={{-58,
                10},{-38,30}})));
      replaceable Transverse.Velocity transverse1(source(k=0),final orientation
          =Orientation.preceding) constrainedby
        Transverse.BaseClasses.PartialCondition
        "<html>1<sup>st</sup> transverse condition</html>" annotation (
          __Dymola_choicesFromPackage=true, Placement(transformation(extent={{-26,
                -2},{-6,18}})));
      replaceable Transverse.Velocity transverse2(source(k=0),final orientation
          =Orientation.following) constrainedby
        Transverse.BaseClasses.PartialCondition
        "<html>2<sup>nd</sup> transverse condition</html>" annotation (
          __Dymola_choicesFromPackage=true, Placement(transformation(extent={{6,
                -18},{26,2}})));
      replaceable Thermal.Temperature thermal(source(k=298.15*U.K))
        constrainedby Thermal.BaseClasses.PartialCondition "Thermal condition"
        annotation (__Dymola_choicesFromPackage=true, Placement(transformation(
              extent={{38,-30},{58,-10}})));
      // Note:  In Dymola 7.4, the value of k must be specified here instead
      // of at the lower level (e.g., Thermal.Temperature) so that the source
      // subcomponent can be replaced by blocks that don't contain the
      // parameter k.

      FCSys.Connectors.Face face
        "Single-species connector for linear momentum and heat" annotation (
          Placement(transformation(extent={{-10,-50},{10,-30}}),
            iconTransformation(extent={{-10,-50},{10,-30}})));

      FCSys.Connectors.RealInputBus u
        "Input bus for values of imposed conditions" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-100,0}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-100,0})));

      FCSys.Connectors.RealOutputBus y "Output bus of measurements" annotation
        (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={100,0}),iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={100,0})));

    equation
      // Normal
      connect(normal.face, face) annotation (Line(
          points={{-48,16},{-48,-30},{0,-30},{0,-40},{5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(u.normal, normal.u) annotation (Line(
          points={{-100,5.55112e-16},{-80,5.55112e-16},{-80,20},{-58,20}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(normal.y, y.normal) annotation (Line(
          points={{-38,20},{80,20},{80,0},{100,0},{100,5.55112e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{6,3},{6,3}}));

      // 1st transverse
      connect(transverse1.face, face) annotation (Line(
          points={{-16,4},{-16,-30},{0,-30},{0,-40},{5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(u.transverse1, transverse1.u) annotation (Line(
          points={{-100,5.55112e-16},{-80,5.55112e-16},{-80,8},{-26,8}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(transverse1.y, y.transverse1) annotation (Line(
          points={{-6,8},{80,8},{80,5.55112e-16},{100,5.55112e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{6,3},{6,3}}));

      // 2nd transverse
      connect(transverse2.face, face) annotation (Line(
          points={{16,-12},{16,-30},{0,-30},{0,-40},{5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(u.transverse2, transverse2.u) annotation (Line(
          points={{-100,5.55112e-16},{-80,5.55112e-16},{-80,-8},{6,-8}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(transverse2.y, y.transverse2) annotation (Line(
          points={{26,-8},{80,-8},{80,0},{100,0},{100,5.55112e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{6,3},{6,3}}));

      // Thermal
      connect(thermal.face, face) annotation (Line(
          points={{48,-24},{48,-30},{0,-30},{0,-40},{5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(u.thermal, thermal.u) annotation (Line(
          points={{-100,5.55112e-16},{-80,5.55112e-16},{-80,-20},{38,-20}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(thermal.y, y.thermal) annotation (Line(
          points={{58,-20},{80,-20},{80,0},{100,0},{100,5.55112e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{6,3},{6,3}}));
      annotation (Icon(graphics));
    end Species;
    extends Modelica.Icons.Package;

    package Normal "Normal mechanical Conditions"
      extends Modelica.Icons.Package;

      model CurrentAreic "Impose areic current (measure normal force)"
        extends BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.CurrentAreic,
          u(final unit="N/(l2.T)"),
          final y(unit="l.m/T2") = face.mPhidot_0);

      equation
        face.J = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="normal");
      end CurrentAreic;

      model Force "Impose normal force (measure areic current)"
        extends BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.Force,
          u(final unit="l.m/T2"),
          final y(unit="N/(l2.T)") = face.J);

      equation
        face.mPhidot_0 = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="normal");
      end Force;

      model Custom "Apply condition to a custom expression"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Custom,
            y=face.mPhidot_0);

        Real x=face.J "Expression to which the condition is applied"
          annotation (Dialog(group="Specification"));

      equation
        x = u_final;
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="normal",
          Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>face.J</code> and/or <code>face.mPhidot_0</code>.</p></html>"));
      end Custom;

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition "Partial model for a normal condition"
          extends Face.BaseClasses.PartialCondition;

          constant ConditionType conditionType "Type of condition";
          // Note:  This is included so that the type of condition is recorded with the
          // results.

        equation
          face.mPhidot = {0,0} "No transverse forces";
          face.Qdot = 0 "Adiabatic";
          annotation (defaultComponentName="normal");
        end PartialCondition;

        type ConditionType = enumeration(
            CurrentAreic "Impose areic current (measure force)",
            Force "Impose force (measure areic current)",
            Custom "Custom expression") "Types of conditions";

      end BaseClasses;

    end Normal;

    package Transverse "Transverse mechanical Conditions"
      extends Modelica.Icons.Package;

      model Velocity "Impose velocity (measure shear force)"
        extends BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.Velocity,
          u(final unit="l/T"),
          final y(unit="l.m/T2") = face.mPhidot[orientation]);

      equation
        face.phi[orientation] = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="transverse");
      end Velocity;

      model Force "Impose shear force (measure velocity)"
        extends BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.Force,
          u(final unit="l.m/T2"),
          final y(unit="l/T") = face.phi[orientation]);

      equation
        face.mPhidot[orientation] = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="transverse");
      end Force;

      model Custom "Apply condition to a custom expression"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Custom,
            y=face.mPhidot[orientation]);

        Real x=face.phi[orientation]
          "Expression to which the condition is applied"
          annotation (Dialog(group="Specification"));

      equation
        x = u_final;
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="normal",
          Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>face.phi[orientation]</code> and/or <code>face.mPhidot[orientation]</code>.</p></html>"));
      end Custom;

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition
          "Partial model for a transverse mechanical condition"
          extends Face.BaseClasses.PartialCondition;

          parameter Orientation orientation=Orientation.preceding
            "Orientation of linear momentum";

          constant ConditionType conditionType "Type of condition";
          // Note:  This is included so that the type of condition is recorded with the
          // results.

        equation
          face.mPhidot_0 = 0 "No normal force";
          face.mPhidot[mod1(orientation + 1, 2)] = 0
            "No force along the other transverse axis";
          face.Qdot = 0 "Adiabatic";
          annotation (defaultComponentName="transverse");
        end PartialCondition;

        type ConditionType = enumeration(
            Velocity "Impose velocity (measure shear force)",
            Force "Impose shear force (measure velocity)",
            Custom "Custom expression") "Types of conditions";

      end BaseClasses;

    end Transverse;

    package Thermal "Thermal Conditions"
      extends Modelica.Icons.Package;

      model Temperature "Impose temperature (measure heat flow rate)"
        extends Thermal.BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.Temperature,
          u(final unit="l2.m/(N.T2)", displayUnit="K"),
          final y(unit="l2.m/T3") = face.Qdot);

      equation
        face.T = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="thermal");
      end Temperature;

      model HeatRate "Impose heat flow rate (measure temperature)"
        extends Thermal.BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.HeatRate,
          u(final unit="l2.m/T3"),
          final y(
            final unit="l2.m/(N.T2)",
            displayUnit="K") = face.T);

      equation
        face.Qdot = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="thermal");
      end HeatRate;

      model Custom "Apply condition to a custom expression"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Custom,
            y=face.Qdot);

        Real x=face.T "Expression to which the condition is applied"
          annotation (Dialog(group="Specification"));

      equation
        x = u_final;
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="normal",
          Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>face.T</code> and/or <code>face.Qdot</code>.</p></html>"));
      end Custom;

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition "Partial model for a thermal condition"
          extends Face.BaseClasses.PartialCondition;

          constant ConditionType conditionType "Type of condition";
          // Note:  This is included so that the type of condition is recorded with the
          // results.

        equation
          face.mPhidot_0 = 0 "No normal force";
          face.mPhidot = {0,0} "No transverse forces";
          annotation (defaultComponentName="thermal");
        end PartialCondition;

        type ConditionType = enumeration(
            Temperature "Impose temperature (measure heat flow rate)",
            HeatRate "Impose heat flow rate (measure temperature)",
            Custom "Custom expression") "Types of conditions";

      end BaseClasses;

    end Thermal;

    package BaseClasses "Base classes (not generally for direct use)"
      extends Modelica.Icons.BasesPackage;

      partial model PartialCondition
        "Partial model to specify and measure conditions on a connector"
        extends FCSys.Conditions.BaseClasses.Icons.Single;

        parameter Boolean internal=true "Use internal condition" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(group="Specification"));

        replaceable Modelica.Blocks.Sources.Constant source if internal
          constrainedby Modelica.Blocks.Interfaces.SO
          "Source of internal condition" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Specification",enable=internal),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-70,30})));

        FCSys.Connectors.RealInput u if not internal
          "Value of imposed condition" annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-100,0})));

        FCSys.Connectors.RealOutput y "Expression of measurement" annotation (
            Dialog(group="Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={100,0}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={100,0})));

        FCSys.Connectors.Face face
          "Connector to transport linear momentum and heat of a single species"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

      protected
        FCSys.Connectors.RealOutputInternal u_final
          "Final value of imposed condition" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=0,
              origin={-20,0})));

      equation
        connect(u, u_final) annotation (Line(
            points={{-100,5.55112e-16},{-62,-4.87687e-22},{-62,5.55112e-16},{-20,
                5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(source.y, u_final) annotation (Line(
            points={{-59,30},{-40,30},{-40,5.55112e-16},{-20,5.55112e-16}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation (Icon(graphics));
      end PartialCondition;

    end BaseClasses;

  end Face;

  package FaceBusPair
    "<html>Conditions for a pair of <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connectors</html>"
    extends Modelica.Icons.Package;

    model Subregion
      "<html>Condition for faces of a <a href=\"modelica://FCSys.Regions.Region\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with efforts by default</html>"
      extends FCSys.Conditions.BaseClasses.Icons.Double;

      replaceable Phases.Gas gas "Gas" annotation (Dialog(group="Phases",
            __Dymola_descriptionLabel=true), Placement(transformation(extent={{
                -10,-10},{10,10}})));

      replaceable Phases.Graphite graphite "Graphite" annotation (Dialog(group=
              "Phases", __Dymola_descriptionLabel=true), Placement(
            transformation(extent={{-10,-10},{10,10}})));

      replaceable Phases.Ionomer ionomer "Ionomer" annotation (Dialog(group=
              "Phases", __Dymola_descriptionLabel=true), Placement(
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
      annotation (defaultComponentPrefixes="replaceable", defaultComponentName=
            "subregion");
    end Subregion;

    model SubregionFlow
      "<html>Condition for faces of a <a href=\"modelica://FCSys.Regions.Region\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model, with flows by default</html>"
      extends Subregion(
        gas(
          H2(
            redeclare replaceable FaceDifferential.Normal.Force normal,
            redeclare replaceable FaceDifferential.Transverse.Force transverse1,

            redeclare replaceable FaceDifferential.Transverse.Force transverse2,

            redeclare replaceable FaceDifferential.Thermal.HeatRate thermal),
          H2O(
            redeclare replaceable FaceDifferential.Normal.Force normal,
            redeclare replaceable FaceDifferential.Transverse.Force transverse1,

            redeclare replaceable FaceDifferential.Transverse.Force transverse2,

            redeclare replaceable FaceDifferential.Thermal.HeatRate thermal),
          N2(
            redeclare replaceable FaceDifferential.Normal.Force normal,
            redeclare replaceable FaceDifferential.Transverse.Force transverse1,

            redeclare replaceable FaceDifferential.Transverse.Force transverse2,

            redeclare replaceable FaceDifferential.Thermal.HeatRate thermal),
          O2(
            redeclare replaceable FaceDifferential.Normal.Force normal,
            redeclare replaceable FaceDifferential.Transverse.Force transverse1,

            redeclare replaceable FaceDifferential.Transverse.Force transverse2,

            redeclare replaceable FaceDifferential.Thermal.HeatRate thermal)),
        graphite('C+'(
            redeclare replaceable FaceDifferential.Normal.Force normal,
            redeclare replaceable FaceDifferential.Transverse.Force transverse1,

            redeclare replaceable FaceDifferential.Transverse.Force transverse2,

            redeclare replaceable FaceDifferential.Thermal.HeatRate thermal),
            'e-'(
            redeclare replaceable FaceDifferential.Normal.Force normal,
            redeclare replaceable FaceDifferential.Transverse.Force transverse1,

            redeclare replaceable FaceDifferential.Transverse.Force transverse2,

            redeclare replaceable FaceDifferential.Thermal.HeatRate thermal)),
        ionomer(
          'C19HF37O5S-'(
            redeclare replaceable FaceDifferential.Normal.Force normal,
            redeclare replaceable FaceDifferential.Transverse.Force transverse1,

            redeclare replaceable FaceDifferential.Transverse.Force transverse2,

            redeclare replaceable FaceDifferential.Thermal.HeatRate thermal),
          H2O(
            redeclare replaceable FaceDifferential.Normal.Force normal,
            redeclare replaceable FaceDifferential.Transverse.Force transverse1,

            redeclare replaceable FaceDifferential.Transverse.Force transverse2,

            redeclare replaceable FaceDifferential.Thermal.HeatRate thermal),
          'H+'(
            redeclare replaceable FaceDifferential.Normal.Force normal,
            redeclare replaceable FaceDifferential.Transverse.Force transverse1,

            redeclare replaceable FaceDifferential.Transverse.Force transverse2,

            redeclare replaceable FaceDifferential.Thermal.HeatRate thermal)));

      annotation (defaultComponentPrefixes="replaceable", defaultComponentName=
            "subregion");

    end SubregionFlow;

    package Phases
      "<html>Conditions for the <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Phase\">Phase</a> model (multi-species)</html>"
      extends Modelica.Icons.Package;

      model Gas "Condition for gas"

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
        FCSys.Conditions.FacePair.Species H2 if inclH2 "Model" annotation (
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
        FCSys.Conditions.FacePair.Species H2O if inclH2O "Model" annotation (
            Dialog(
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
        FCSys.Conditions.FacePair.Species N2 if inclN2 "Model" annotation (
            Dialog(
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
        FCSys.Conditions.FacePair.Species O2 if inclO2 "Model" annotation (
            Dialog(
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

      model Graphite "Condition for graphite"

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
        FCSys.Conditions.FacePair.Species 'C+' if 'inclC+' "Model" annotation (
            Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclC+'), Placement(transformation(extent={{-10,-10},{10,10}})));

        parameter Boolean 'incle-'=false
          "<html>Electrons (e<sup>-</sup>)</html>" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.Conditions.FacePair.Species 'e-' if 'incle-' "Model" annotation (
            Dialog(
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

      model Ionomer "Condition for ionomer"

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
        FCSys.Conditions.FacePair.Species 'C19HF37O5S-' if 'inclC19HF37O5S-'
          "Model" annotation (Dialog(
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
        FCSys.Conditions.FacePair.Species 'H+' if 'inclH+' "Model" annotation (
            Dialog(
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
        FCSys.Conditions.FacePair.Species H2O if inclH2O "Model" annotation (
            Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,10}})));

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

      model Liquid "Condition for liquid"

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
        FCSys.Conditions.FacePair.Species H2O if inclH2O "Model" annotation (
            Dialog(
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

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;
        model NullPhase "Empty condition for a phase (no species)"
          extends FCSys.Conditions.BaseClasses.Icons.Double;

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
 <a href=\"modelica://FCSys.Conditions.FaceBus\">Conditions.FaceBus</a> package.
 For more information, please see the documentation in that package.</p></html>"));

  end FaceBusPair;

  package FaceBusPair2
    "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, e.g., of a <a href=\"modelica://FCSys.Regions.Region\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model</html>"

    extends Modelica.Icons.Package;

    model Subregion
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, with efforts by default</html>"

      extends FCSys.Conditions.BaseClasses.Icons.Single;

      Phases.Gas gas "Gas" annotation (Dialog(group="Phases",
            __Dymola_descriptionLabel=true), Placement(transformation(extent={{
                -10,-10},{10,10}})));

      Phases.Graphite graphite "Graphite" annotation (Dialog(group="Phases",
            __Dymola_descriptionLabel=true), Placement(transformation(extent={{
                -10,-10},{10,10}})));

      Phases.Ionomer ionomer "Ionomer" annotation (Dialog(group="Phases",
            __Dymola_descriptionLabel=true), Placement(transformation(extent={{
                -10,-10},{10,10}})));

      Phases.Liquid liquid "Liquid" annotation (Dialog(group="Phases",
            __Dymola_descriptionLabel=true), Placement(transformation(extent={{
                -10,-10},{10,10}})));

      FCSys.Connectors.FaceBus face
        "Connector for linear momentum and heat of multiple species"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}}),
            iconTransformation(extent={{-10,-50},{10,-30}})));
      FCSys.Connectors.RealInputBus u "Bus of inputs to specify conditions"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-100,0})));

      FCSys.Connectors.RealOutputBus y "Output bus of measurements" annotation
        (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={100,0})));

    equation
      // Gas
      connect(gas.face, face.gas) annotation (Line(
          points={{6.10623e-16,-4},{8.60423e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(u.gas, gas.u) annotation (Line(
          points={{-100,5.55112e-16},{-10,5.55112e-16},{-10,6.10623e-16}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(gas.y, y.gas) annotation (Line(
          points={{10,6.10623e-16},{100,6.10623e-16},{100,5.55112e-16}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      // Graphite
      connect(graphite.face, face.graphite) annotation (Line(
          points={{6.10623e-16,-4},{8.60423e-16,-40},{5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(u.graphite, graphite.u) annotation (Line(
          points={{-100,5.55112e-16},{-10,5.55112e-16},{-10,6.10623e-16}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(graphite.y, y.graphite) annotation (Line(
          points={{10,6.10623e-16},{100,6.10623e-16},{100,5.55112e-16}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      // Ionomer
      connect(ionomer.face, face.ionomer) annotation (Line(
          points={{6.10623e-16,-4},{-4.87687e-22,-4},{-4.87687e-22,-40},{
              5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));
      connect(u.ionomer, ionomer.u) annotation (Line(
          points={{-100,5.55112e-16},{-10,5.55112e-16},{-10,6.10623e-16}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(ionomer.y, y.ionomer) annotation (Line(
          points={{10,6.10623e-16},{100,6.10623e-16},{100,5.55112e-16}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      // Liquid
      connect(liquid.face, face.liquid) annotation (Line(
          points={{6.10623e-16,-4},{-4.87687e-22,-4},{-4.87687e-22,-40},{
              5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));
      connect(u.liquid, liquid.u) annotation (Line(
          points={{-100,5.55112e-16},{-10,5.55112e-16},{-10,6.10623e-16}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(liquid.y, y.liquid) annotation (Line(
          points={{10,6.10623e-16},{100,6.10623e-16},{100,5.55112e-16}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (defaultComponentPrefixes="replaceable", defaultComponentName=
            "subregion");
    end Subregion;

    model SubregionFlow
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, with flows by default</html>"

      extends Subregion(
        gas(
          H2(
            redeclare replaceable Face.Normal.Force normal,
            redeclare replaceable Face.Transverse.Force transverse1,
            redeclare replaceable Face.Transverse.Force transverse2,
            redeclare replaceable Face.Thermal.HeatRate thermal(source(k=0))),
          H2O(
            redeclare replaceable Face.Normal.Force normal,
            redeclare replaceable Face.Transverse.Force transverse1,
            redeclare replaceable Face.Transverse.Force transverse2,
            redeclare replaceable Face.Thermal.HeatRate thermal(source(k=0))),
          N2(
            redeclare replaceable Face.Normal.Force normal,
            redeclare replaceable Face.Transverse.Force transverse1,
            redeclare replaceable Face.Transverse.Force transverse2,
            redeclare replaceable Face.Thermal.HeatRate thermal(source(k=0))),
          O2(
            redeclare replaceable Face.Normal.Force normal,
            redeclare replaceable Face.Transverse.Force transverse1,
            redeclare replaceable Face.Transverse.Force transverse2,
            redeclare replaceable Face.Thermal.HeatRate thermal(source(k=0)))),

        graphite('C+'(
            redeclare replaceable Face.Normal.Force normal,
            redeclare replaceable Face.Transverse.Force transverse1,
            redeclare replaceable Face.Transverse.Force transverse2,
            redeclare replaceable Face.Thermal.HeatRate thermal(source(k=0))),
            'e-'(
            redeclare replaceable Face.Normal.Force normal,
            redeclare replaceable Face.Transverse.Force transverse1,
            redeclare replaceable Face.Transverse.Force transverse2,
            redeclare replaceable Face.Thermal.HeatRate thermal(source(k=0)))),

        ionomer(
          'C19HF37O5S-'(
            redeclare replaceable Face.Normal.Force normal,
            redeclare replaceable Face.Transverse.Force transverse1,
            redeclare replaceable Face.Transverse.Force transverse2,
            redeclare replaceable Face.Thermal.HeatRate thermal(source(k=0))),
          'H+'(
            redeclare replaceable Face.Normal.Force normal,
            redeclare replaceable Face.Transverse.Force transverse1,
            redeclare replaceable Face.Transverse.Force transverse2,
            redeclare replaceable Face.Thermal.HeatRate thermal(source(k=0))),
          H2O(
            redeclare replaceable Face.Normal.Force normal,
            redeclare replaceable Face.Transverse.Force transverse1,
            redeclare replaceable Face.Transverse.Force transverse2,
            redeclare replaceable Face.Thermal.HeatRate thermal(source(k=0)))),

        liquid(H2O(
            redeclare replaceable Face.Normal.Force normal,
            redeclare replaceable Face.Transverse.Force transverse1,
            redeclare replaceable Face.Transverse.Force transverse2,
            redeclare replaceable Face.Thermal.HeatRate thermal(source(k=0)))));

      annotation (defaultComponentPrefixes="replaceable",defaultComponentName=
            "subregion");

    end SubregionFlow;

    package Phases
      "<html>Conditions for the <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Phase\">Phase</a> model (multi-species)</html>"
      extends Modelica.Icons.Package;

      model Gas "Condition for gas"

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
        Face.Species H2 if inclH2 "Conditions" annotation (Dialog(
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
        Face.Species H2O if inclH2O "Conditions" annotation (Dialog(
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
        Face.Species N2 if inclN2 "Conditions" annotation (Dialog(
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
        Face.Species O2 if inclO2 "Conditions" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclO2), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // H2
        connect(H2.face, face.H2) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.H2, H2.u) annotation (Line(
            points={{-100,5.55112e-16},{-100,0},{-10,0},{-10,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(H2.y, y.H2) annotation (Line(
            points={{10,6.10623e-16},{10,0},{100,0},{100,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2O
        connect(H2O.face, face.H2O) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.H2O, H2O.u) annotation (Line(
            points={{-100,5.55112e-16},{-100,0},{-10,0},{-10,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(H2O.y, y.H2O) annotation (Line(
            points={{10,6.10623e-16},{10,0},{100,0},{100,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // N2
        connect(N2.face, face.N2) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.N2, N2.u) annotation (Line(
            points={{-100,5.55112e-16},{-100,0},{-10,0},{-10,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(N2.y, y.N2) annotation (Line(
            points={{10,6.10623e-16},{10,0},{100,0},{100,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // O2
        connect(O2.face, face.O2) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.O2, O2.u) annotation (Line(
            points={{-100,5.55112e-16},{-100,0},{-10,0},{-10,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(O2.y, y.O2) annotation (Line(
            points={{10,6.10623e-16},{10,0},{100,0},{100,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (Diagram(graphics));
      end Gas;

      model Graphite "Condition for graphite"

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
        Face.Species 'C+' if 'inclC+' "Conditions" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclC+'), Placement(transformation(extent={{-10,-10},{10,10}})));

        parameter Boolean 'incle-'=false
          "<html>Electrons (e<sup>-</sup>)</html>" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species 'e-' if 'incle-' "Conditions" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='incle-'), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // C+
        connect('C+'.face, face.'C+') annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'C+', 'C+'.u) annotation (Line(
            points={{-100,5.55112e-16},{-100,0},{-10,0},{-10,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect('C+'.y, y.'C+') annotation (Line(
            points={{10,6.10623e-16},{10,0},{100,0},{100,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // e-
        connect('e-'.face, face.'e-') annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'e-', 'e-'.u) annotation (Line(
            points={{-100,5.55112e-16},{-100,0},{-10,0},{-10,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect('e-'.y, y.'e-') annotation (Line(
            points={{10,6.10623e-16},{10,0},{100,0},{100,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

      end Graphite;

      model Ionomer "Condition for ionomer"

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
        Face.Species 'C19HF37O5S-' if 'inclC19HF37O5S-' "Conditions"
          annotation (Dialog(
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
        Face.Species 'H+' if 'inclH+' "Conditions" annotation (Dialog(
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
        Face.Species H2O if inclH2O "Conditions" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // C19HF37O5S-
        connect('C19HF37O5S-'.face, face.'C19HF37O5S-') annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'C19HF37O5S-', 'C19HF37O5S-'.u) annotation (Line(
            points={{-100,5.55112e-16},{-100,0},{-10,0},{-10,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect('C19HF37O5S-'.y, y.'C19HF37O5S-') annotation (Line(
            points={{10,6.10623e-16},{10,0},{100,0},{100,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H+
        connect('H+'.face, face.'H+') annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.'H+', 'H+'.u) annotation (Line(
            points={{-100,5.55112e-16},{-100,0},{-10,0},{-10,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect('H+'.y, y.'H+') annotation (Line(
            points={{10,6.10623e-16},{10,0},{100,0},{100,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2O
        connect(H2O.face, face.H2O) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.H2O, H2O.u) annotation (Line(
            points={{-100,5.55112e-16},{-100,0},{-10,0},{-10,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(H2O.y, y.H2O) annotation (Line(
            points={{10,6.10623e-16},{10,0},{100,0},{100,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

      end Ionomer;

      model Liquid "Condition for liquid"

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
        Face.Species H2O if inclH2O "Conditions" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // H2O
        connect(H2O.face, face.H2O) annotation (Line(
            points={{6.10623e-16,-4},{1.16573e-15,-40},{5.55112e-16,-40}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(u.H2O, H2O.u) annotation (Line(
            points={{-100,5.55112e-16},{-100,0},{-10,0},{-10,6.10623e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(H2O.y, y.H2O) annotation (Line(
            points={{10,6.10623e-16},{10,0},{100,0},{100,5.55112e-16}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

      end Liquid;

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;
        model NullPhase "Empty condition for a phase (no species)"
          extends FCSys.Conditions.BaseClasses.Icons.Single;

          FCSys.Connectors.FaceBus face
            "Multi-species connector for linear momentum and heat"
            annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
          FCSys.Connectors.RealInputBus u
            "Input bus for values of imposed conditions" annotation (Placement(
                transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-100,0}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={-100,0})));

          FCSys.Connectors.RealOutputBus y "Output bus of measurements"
            annotation (Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={100,0}),iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=0,
                origin={100,0})));

        end NullPhase;

      end BaseClasses;

    end Phases;

  end FaceBusPair2;

  package FacePair
    "<html>Conditions for a pair of <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connectors</html>"

    extends Modelica.Icons.Package;

    model Species
      "<html>Condition for a pair of faces of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"
      extends FCSys.Conditions.BaseClasses.Icons.Double;

      // Normal
      replaceable Normal.CurrentAreic normal constrainedby
        Normal.BaseClasses.PartialCondition "Normal condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Linear momentum", __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-70,0},{-50,20}})));

      // 1st transverse
      replaceable Transverse.Velocity transverse1(spec(k=0),final orientation=
            Orientation.preceding) constrainedby
        Transverse.BaseClasses.PartialCondition
        "<html>1<sup>st</sup> transverse condition</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Linear momentum", __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-26,-14},{-6,6}})));

      // 2nd transverse
      replaceable Transverse.Velocity transverse2(spec(k=0),final orientation=
            Orientation.following) constrainedby
        Transverse.BaseClasses.PartialCondition
        "<html>2<sup>nd</sup> transverse condition</html>" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Linear momentum", __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{6,-28},{26,-8}})));

      // Thermal
      replaceable Thermal.Temperature thermal(spec(k(start=298.15*U.K)))
        constrainedby Thermal.BaseClasses.PartialCondition "Type of condition"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(group="Heat", __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{50,-40},{70,-20}})));
      // Note:  In Dymola 7.4, the value of k must be specified here instead
      // of at the lower level (e.g., Thermal.Temperature) so that the spec
      // subcomponent can be replaced by blocks that don't contain the
      // parameter k.
      FCSys.Connectors.RealInputBus u "Input bus for external signal sources"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,40})));

      FCSys.Connectors.Face negative "Negative face" annotation (Placement(
            transformation(extent={{-110,-40},{-90,-20}}), iconTransformation(
              extent={{-110,-10},{-90,10}})));
      FCSys.Connectors.Face positive "Positive face" annotation (Placement(
            transformation(extent={{90,-40},{110,-20}}), iconTransformation(
              extent={{90,-10},{110,10}})));

    equation
      // Normal
      connect(u.normal, normal.u) annotation (Line(
          points={{5.55112e-16,40},{0,40},{0,20},{-60,20},{-60,14}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(normal.negative, negative) annotation (Line(
          points={{-70,10},{-80,10},{-80,-30},{-100,-30}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(normal.positive, positive) annotation (Line(
          points={{-50,10},{90,10},{90,-30},{100,-30}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));

      // 1st transverse
      connect(u.transverse1, transverse1.u) annotation (Line(
          points={{5.55112e-16,40},{0,40},{0,20},{-16,20},{-16,6.66134e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(transverse1.negative, negative) annotation (Line(
          points={{-26,-4},{-80,-4},{-80,-30},{-100,-30}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(transverse1.positive, positive) annotation (Line(
          points={{-6,-4},{90,-4},{90,-30},{100,-30}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));

      // 2nd transverse
      connect(u.transverse2, transverse2.u) annotation (Line(
          points={{5.55112e-16,40},{5.55112e-16,20},{16,20},{16,-14}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));

      connect(transverse2.negative, negative) annotation (Line(
          points={{6,-18},{-80,-18},{-80,-30},{-100,-30}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(transverse2.positive, positive) annotation (Line(
          points={{26,-18},{90,-18},{90,-30},{100,-30}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));

      // Thermal
      connect(u.thermal, thermal.u) annotation (Line(
          points={{5.55112e-16,40},{0,40},{0,20},{60,20},{60,-26}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(thermal.negative, negative) annotation (Line(
          points={{50,-30},{-26,-30},{-26,-30},{-100,-30}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(thermal.positive, positive) annotation (Line(
          points={{70,-30},{86,-30},{86,-30},{100,-30}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      annotation (Icon(graphics));
    end Species;

    package Normal "Normal mechanical Conditions"
      extends Modelica.Icons.Package;

      model CurrentAreic "Impose areic current, with material conservation"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.CurrentAreic,
            u(final unit="N/(l2.T)"));

      equation
        negative.J - positive.J = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="normal");
      end CurrentAreic;

      model Force "Impose normal force, with conservation of linear momentum"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Force,
            u(final unit="l.m/T2"));

      equation
        negative.mPhidot_0 = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="transverse");
      end Force;

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition
          "Partial model for a normal mechanical condition"

          extends FacePair.BaseClasses.PartialCondition;

          constant ConditionType conditionType "Type of condition";
          // Note:  This is included so that the type of condition is recorded with the
          // results.

        equation
          negative.mPhidot = positive.mPhidot;
          negative.Qdot = positive.Qdot;
          annotation (defaultComponentName="normal");
        end PartialCondition;

        type ConditionType = enumeration(
            CurrentAreic "Areic current",
            Force "Force") "Types of Conditions";

      end BaseClasses;

    end Normal;

    package Transverse "Transverse mechanical Conditions"
      extends Modelica.Icons.Package;

      model Velocity
        "Impose shear velocity, with conservation of linear momentum"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Velocity,
            u(final unit="l/T"));

      equation
        negative.phi[orientation] - positive.phi[orientation] = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="transverse");
      end Velocity;

      model Force "Impose shear force, with conservation of linear momentum"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Force,
            u(final unit="l.m/T2"));

      equation
        negative.mPhidot[orientation] = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="transverse");
      end Force;

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition
          "Partial model of a condition for linear momentum"

          extends FacePair.BaseClasses.PartialCondition;

          constant ConditionType conditionType "Type of condition";
          // Note:  This is included so that the type of condition is recorded with the
          // results.

          parameter Orientation orientation=Orientation.preceding
            "Orientation of linear momentum";

        equation
          negative.mPhidot_0 = positive.mPhidot_0;
          negative.mPhidot[mod1(orientation + 1, 2)] = positive.mPhidot[mod1(
            orientation + 1, 2)];
          negative.Qdot = positive.Qdot;
          annotation (defaultComponentName="transverse");
        end PartialCondition;

        type ConditionType = enumeration(
            Velocity "Shear velocity",
            Force "Shear force") "Types of Conditions";

      end BaseClasses;

    end Transverse;

    package Thermal "Thermal Conditions"
      extends Modelica.Icons.Package;

      model Temperature
        "Impose temperature difference, with energy conservation"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Temperature,
            redeclare FCSys.Connectors.RealInput u(final unit="l2.m/(N.T2)",
              displayUnit="K"));

      equation
        negative.T - positive.T = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="thermal");
      end Temperature;

      model HeatRate "Impose heat flow rate, with energy conservation"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.HeatRate,
            redeclare FCSys.Connectors.RealInput u(final unit="l2.m/T3"));

      equation
        negative.Qdot = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="thermal");
      end HeatRate;

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition "Partial model of a condition for heat"

          extends FacePair.BaseClasses.PartialCondition;

          constant ConditionType conditionType "Type of condition";
          // Note:  This is included so that the type of condition is recorded with the
          // results.

        equation
          negative.mPhidot_0 = positive.mPhidot_0;
          negative.mPhidot = positive.mPhidot;
          annotation (defaultComponentName="thermal");
        end PartialCondition;

        type ConditionType = enumeration(
            Temperature "Temperature difference",
            HeatRate "Heat flow rate") "Types of Conditions";

      end BaseClasses;

    end Thermal;

    package BaseClasses "Base classes (not generally for direct use)"
      extends Modelica.Icons.BasesPackage;

      partial model PartialCondition "Partial model of a condition"
        extends FCSys.Conditions.BaseClasses.Icons.Double;

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

        FCSys.Connectors.RealInput u if not internal "Value of condition"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));
        FCSys.Connectors.Face negative "Connector for the negative face"
          annotation (Placement(transformation(extent={{-110,-10},{-90,10}}),
              iconTransformation(extent={{-110,-10},{-90,10}})));
        FCSys.Connectors.Face positive "Connector for the positive face"
          annotation (Placement(transformation(extent={{90,-10},{110,10}}),
              iconTransformation(extent={{90,-10},{110,10}})));

      protected
        FCSys.Connectors.RealInputInternal u_final "Final value of condition"
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

        // Conservation (no storage)
        0 = negative.mPhidot_0 + positive.mPhidot_0 "Normal linear momentum";
        {0,0} = negative.mPhidot + positive.mPhidot
          "Transverse linear momentum";
        0 = negative.Qdot + positive.Qdot "Energy conservation (no storage)";
        annotation (Diagram(graphics));
      end PartialCondition;

    end BaseClasses;

  end FacePair;

  package FacePair2
    "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"
    model Species
      "<html>Conditions for a <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connector</html>"

      extends FCSys.Conditions.BaseClasses.Icons.Single;

      replaceable Normal.CurrentAreic normal(source(k=0)) constrainedby
        Normal.BaseClasses.PartialCondition "Normal condition" annotation (
          __Dymola_choicesFromPackage=true, Placement(transformation(extent={{-58,
                10},{-38,30}})));
      replaceable Transverse.Velocity transverse1(source(k=0),final orientation
          =Orientation.preceding) constrainedby
        Transverse.BaseClasses.PartialCondition
        "<html>1<sup>st</sup> transverse condition</html>" annotation (
          __Dymola_choicesFromPackage=true, Placement(transformation(extent={{-26,
                -2},{-6,18}})));
      replaceable Transverse.Velocity transverse2(source(k=0),final orientation
          =Orientation.following) constrainedby
        Transverse.BaseClasses.PartialCondition
        "<html>2<sup>nd</sup> transverse condition</html>" annotation (
          __Dymola_choicesFromPackage=true, Placement(transformation(extent={{6,
                -18},{26,2}})));
      replaceable Thermal.Temperature thermal(source(k=298.15*U.K))
        constrainedby Thermal.BaseClasses.PartialCondition "Thermal condition"
        annotation (__Dymola_choicesFromPackage=true, Placement(transformation(
              extent={{38,-30},{58,-10}})));
      // Note:  In Dymola 7.4, the value of k must be specified here instead
      // of at the lower level (e.g., Thermal.Temperature) so that the source
      // subcomponent can be replaced by blocks that don't contain the
      // parameter k.

      FCSys.Connectors.Face face
        "Single-species connector for linear momentum and heat" annotation (
          Placement(transformation(extent={{-10,-50},{10,-30}}),
            iconTransformation(extent={{-10,-50},{10,-30}})));

      FCSys.Connectors.RealInputBus u
        "Input bus for values of imposed conditions" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-100,0}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-100,0})));

      FCSys.Connectors.RealOutputBus y "Output bus of measurements" annotation
        (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={100,0}),iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={100,0})));

    equation
      // Normal
      connect(normal.face, face) annotation (Line(
          points={{-48,16},{-48,-30},{0,-30},{0,-40},{5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(u.normal, normal.u) annotation (Line(
          points={{-100,5.55112e-16},{-80,5.55112e-16},{-80,20},{-58,20}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(normal.y, y.normal) annotation (Line(
          points={{-38,20},{80,20},{80,0},{100,0},{100,5.55112e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{6,3},{6,3}}));

      // 1st transverse
      connect(transverse1.face, face) annotation (Line(
          points={{-16,4},{-16,-30},{0,-30},{0,-40},{5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(u.transverse1, transverse1.u) annotation (Line(
          points={{-100,5.55112e-16},{-80,5.55112e-16},{-80,8},{-26,8}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(transverse1.y, y.transverse1) annotation (Line(
          points={{-6,8},{80,8},{80,5.55112e-16},{100,5.55112e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{6,3},{6,3}}));

      // 2nd transverse
      connect(transverse2.face, face) annotation (Line(
          points={{16,-12},{16,-30},{0,-30},{0,-40},{5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(u.transverse2, transverse2.u) annotation (Line(
          points={{-100,5.55112e-16},{-80,5.55112e-16},{-80,-8},{6,-8}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(transverse2.y, y.transverse2) annotation (Line(
          points={{26,-8},{80,-8},{80,0},{100,0},{100,5.55112e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{6,3},{6,3}}));

      // Thermal
      connect(thermal.face, face) annotation (Line(
          points={{48,-24},{48,-30},{0,-30},{0,-40},{5.55112e-16,-40}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(u.thermal, thermal.u) annotation (Line(
          points={{-100,5.55112e-16},{-80,5.55112e-16},{-80,-20},{38,-20}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%first",
          index=-1,
          extent={{-6,3},{-6,3}}));
      connect(thermal.y, y.thermal) annotation (Line(
          points={{58,-20},{80,-20},{80,0},{100,0},{100,5.55112e-16}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{6,3},{6,3}}));
      annotation (Icon(graphics));
    end Species;
    extends Modelica.Icons.Package;

    package Normal "Normal mechanical Conditions"
      extends Modelica.Icons.Package;

      model CurrentAreic "Impose areic current (measure normal force)"
        extends BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.CurrentAreic,
          u(final unit="N/(l2.T)"),
          final y(unit="l.m/T2") = face.mPhidot_0);

      equation
        face.J = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="normal");
      end CurrentAreic;

      model Force "Impose normal force (measure areic current)"
        extends BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.Force,
          u(final unit="l.m/T2"),
          final y(unit="N/(l2.T)") = face.J);

      equation
        face.mPhidot_0 = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="normal");
      end Force;

      model Custom "Apply condition to a custom expression"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Custom,
            y=face.mPhidot_0);

        Real x=face.J "Expression to which the condition is applied"
          annotation (Dialog(group="Imposition"));

      equation
        x = u_final;
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="normal",
          Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>face.J</code> and/or <code>face.mPhidot_0</code>.</p></html>"));
      end Custom;

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition "Partial model for a normal condition"
          extends FCSys.Conditions.FacePair2.BaseClasses.PartialCondition;

          constant ConditionType conditionType "Type of condition";
          // Note:  This is included so that the type of condition is recorded with the
          // results.

        equation
          face.mPhidot = {0,0} "No transverse forces";
          face.Qdot = 0 "Adiabatic";
          annotation (defaultComponentName="normal");
        end PartialCondition;

        type ConditionType = enumeration(
            CurrentAreic "Impose areic current (measure force)",
            Force "Impose force (measure areic current)",
            Custom "Custom expression") "Types of conditions";

      end BaseClasses;

    end Normal;

    package Transverse "Transverse mechanical Conditions"
      extends Modelica.Icons.Package;

      model Velocity "Impose velocity (measure shear force)"
        extends BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.Velocity,
          u(final unit="l/T"),
          final y(unit="l.m/T2") = face.mPhidot[orientation]);

      equation
        face.phi[orientation] = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="transverse");
      end Velocity;

      model Force "Impose shear force (measure velocity)"
        extends BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.Force,
          u(final unit="l.m/T2"),
          final y(unit="l/T") = face.phi[orientation]);

      equation
        face.mPhidot[orientation] = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="transverse");
      end Force;

      model Custom "Apply condition to a custom expression"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Custom,
            y=face.mPhidot[orientation]);

        Real x=face.phi[orientation]
          "Expression to which the condition is applied"
          annotation (Dialog(group="Imposition"));

      equation
        x = u_final;
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="normal",
          Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>face.phi[orientation]</code> and/or <code>face.mPhidot[orientation]</code>.</p></html>"));
      end Custom;

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition
          "Partial model for a transverse mechanical condition"
          extends FCSys.Conditions.FacePair2.BaseClasses.PartialCondition;

          parameter Orientation orientation=Orientation.preceding
            "Orientation of linear momentum";

          constant ConditionType conditionType "Type of condition";
          // Note:  This is included so that the type of condition is recorded with the
          // results.

        equation
          face.mPhidot_0 = 0 "No normal force";
          face.mPhidot[mod1(orientation + 1, 2)] = 0
            "No force along the other transverse axis";
          face.Qdot = 0 "Adiabatic";
          annotation (defaultComponentName="transverse");
        end PartialCondition;

        type ConditionType = enumeration(
            Velocity "Impose velocity (measure shear force)",
            Force "Impose shear force (measure velocity)",
            Custom "Custom expression") "Types of conditions";

      end BaseClasses;

    end Transverse;

    package Thermal "Thermal Conditions"
      extends Modelica.Icons.Package;

      model Temperature "Impose temperature (measure heat flow rate)"
        extends Thermal.BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.Temperature,
          u(final unit="l2.m/(N.T2)", displayUnit="K"),
          final y(unit="l2.m/T3") = face.Qdot);

      equation
        face.T = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="thermal");
      end Temperature;

      model HeatRate "Impose heat flow rate (measure temperature)"
        extends Thermal.BaseClasses.PartialCondition(
          final conditionType=BaseClasses.ConditionType.HeatRate,
          u(final unit="l2.m/T3"),
          final y(
            final unit="l2.m/(N.T2)",
            displayUnit="K") = face.T);

      equation
        face.Qdot = u_final;
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="thermal");
      end HeatRate;

      model Custom "Apply condition to a custom expression"
        extends BaseClasses.PartialCondition(final conditionType=BaseClasses.ConditionType.Custom,
            y=face.Qdot);

        Real x=face.T "Expression to which the condition is applied"
          annotation (Dialog(group="Imposition"));

      equation
        x = u_final;
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="normal",
          Documentation(info="<html><p>The expression to which the condition is applied (<code>x</code>)
    must involve <code>face.T</code> and/or <code>face.Qdot</code>.</p></html>"));
      end Custom;

      package BaseClasses "Base classes (not generally for direct use)"
        extends Modelica.Icons.BasesPackage;
        partial model PartialCondition "Partial model for a thermal condition"
          extends FCSys.Conditions.FacePair2.BaseClasses.PartialCondition;

          constant ConditionType conditionType "Type of condition";
          // Note:  This is included so that the type of condition is recorded with the
          // results.

        equation
          face.mPhidot_0 = 0 "No normal force";
          face.mPhidot = {0,0} "No transverse forces";
          annotation (defaultComponentName="thermal");
        end PartialCondition;

        type ConditionType = enumeration(
            Temperature "Impose temperature (measure heat flow rate)",
            HeatRate "Impose heat flow rate (measure temperature)",
            Custom "Custom expression") "Types of conditions";

      end BaseClasses;

    end Thermal;

    package BaseClasses "Base classes (not generally for direct use)"
      extends Modelica.Icons.BasesPackage;

      partial model PartialCondition "Partial model of a condition"
        extends FCSys.Conditions.BaseClasses.Icons.Double;

        parameter Boolean internal=true "Use internal condition" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(group="Imposition"));

        replaceable Modelica.Blocks.Sources.Constant source if internal
          constrainedby Modelica.Blocks.Interfaces.SO
          "Source of internal condition" annotation (
          __Dymola_choicesFromPackage=true,
          Dialog(group="Imposition",enable=internal),
          Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={30,20})));

        FCSys.Connectors.RealInput u if not internal
          "Value of imposed condition" annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,40})));

        FCSys.Connectors.RealOutput y "Expression of measurement" annotation (
            Dialog(group="Measurement"), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-40}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-40})));

        FCSys.Connectors.Face negative "Negative-side connector"
          annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
        FCSys.Connectors.Face positive "Positive-side connector"
          annotation (Placement(transformation(extent={{90,-10},{110,10}})));

      protected
        FCSys.Connectors.RealOutputInternal u_final
          "Final value of imposed condition" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-12})));

      equation
        connect(u, u_final) annotation (Line(
            points={{5.55112e-16,40},{-4.87687e-22,14},{5.55112e-16,14},{
                5.55112e-16,-12}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(source.y, u_final) annotation (Line(
            points={{30,9},{30,0},{0,0},{0,-12},{5.55112e-16,-12}},
            color={0,0,127},
            smooth=Smooth.None));
        annotation (Icon(graphics));
      end PartialCondition;

    end BaseClasses;

  end FacePair2;

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
    annotation (Documentation(info="<html>
<p>This model acts as a connection switch.
It has a single parameter, <code>crossOver</code>.</p>

<p>If <code>crossOver</code> is
set to <code>false</code>, then
the router will be in the pass-through mode.  In that case,
<code>negative1</code> is connected to <code>positive1</code> and <code>negative2</code>
is connected to <code>positive2</code>, as shown by Figure 1a.</p>

<p>If <code>crossOver</code> is set to <code>true</code>, then the router will be in cross-over mode.  In that case, <code>negative1</code> is connected to <code>positive2</code>
and <code>negative2</code> is
connected to <code>positive1</code>, as shown by Figure 1b.</p>

    <table border=\"0\" cellspacing=\"0\" cellpadding=\"2\" align=center>
      <tr align=center>
        <td align=center width=120>
          <img src=\"modelica://FCSys/resources/documentation/Conditions/Router/PassThrough.png\">
<br><b>a:</b>  Pass-through
        </td>
        <td align=center>
          <img src=\"modelica://FCSys/resources/documentation/Conditions/Router/CrossOver.png\">
<br><b>b:</b>  Cross-over
        </td>
      </tr>
      <tr align=center>
        <td colspan=2 align=center><b>Figure 1:</b> Modes of connection.</td>
      </tr>
    </table>
</html>"), Icon(graphics={
          Line(
            points={{-80,40},{-40,40},{0,0},{40,-40},{80,-40}},
            color={127,127,127},
            thickness=0.5,
            visible=crossOver,
            smooth=Smooth.Bezier),
          Line(
            points={{-80,40},{80,40}},
            color={127,127,127},
            visible=not crossOver,
            smooth=Smooth.None,
            thickness=0.5),
          Line(
            points={{-80,-40},{80,-40}},
            color={127,127,127},
            visible=not crossOver,
            smooth=Smooth.None,
            thickness=0.5),
          Line(
            points={{-80,-40},{-40,-40},{0,0},{40,40},{80,40}},
            color={127,127,127},
            thickness=0.5,
            visible=crossOver,
            smooth=Smooth.Bezier)}));
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
For simulation, specify global default settings by dragging FCSys.Conditions.Environment into your model.
The default global default settings will be used for the current simulation.",
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={
          Rectangle(
            extent={{-120,60},{120,100}},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Text(
            extent={{-120,60},{120,100}},
            textString="%name",
            lineColor={0,0,0}),
          Rectangle(
            extent={{-80,60},{80,-100}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.Dash),
          Rectangle(
            extent={{-70,50},{70,-98}},
            lineColor={255,255,255},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={170,170,255}),
          Rectangle(
            extent={{-72,-60},{72,-98}},
            fillPattern=FillPattern.Solid,
            fillColor={255,255,255},
            pattern=LinePattern.None,
            lineColor={0,0,0}),
          Line(points={{-70,-60},{70,-60}}, color={0,0,0}),
          Line(points={{-40,-20},{-10,-50},{40,0}}, color={0,0,0}),
          Ellipse(
            extent={{30,10},{50,-10}},
            pattern=LinePattern.None,
            lineColor={255,255,255},
            fillColor={240,0,0},
            fillPattern=FillPattern.Sphere),
          Line(points={{-66,-90},{-36,-60}}, color={0,0,0}),
          Line(points={{2,-90},{32,-60}}, color={0,0,0}),
          Line(points={{36,-90},{66,-60}}, color={0,0,0}),
          Line(points={{-32,-90},{-2,-60}}, color={0,0,0}),
          Rectangle(
            extent={{70,50},{76,-60}},
            fillPattern=FillPattern.Solid,
            fillColor={255,255,255},
            pattern=LinePattern.None,
            lineColor={0,0,0}),
          Rectangle(
            extent={{-76,50},{-70,-60}},
            fillPattern=FillPattern.Solid,
            fillColor={255,255,255},
            pattern=LinePattern.None,
            lineColor={0,0,0})}));

  end Environment;

  package BaseClasses "Base classes (not generally for direct use)"
    extends Modelica.Icons.BasesPackage;

    package Icons "Icons for conditions"
      extends Modelica.Icons.Package;
      partial class Double "Icon for a two-connector boundary condition"
        //extends Names.Middle;
        annotation (Icon(graphics={
              Rectangle(
                extent={{-100,40},{100,-40}},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid,
                pattern=LinePattern.None),
              Line(
                points={{-100,40},{100,40}},
                pattern=LinePattern.None,
                smooth=Smooth.None),
              Line(
                points={{-100,-40},{-100,40}},
                color={0,0,0},
                smooth=Smooth.None,
                pattern=LinePattern.Dash),
              Text(
                extent={{-150,-20},{150,20}},
                textString="%name",
                lineColor={0,0,0}),
              Line(
                points={{-100,-40},{100,-40}},
                pattern=LinePattern.None,
                smooth=Smooth.None),
              Line(
                points={{100,-40},{100,40}},
                color={0,0,0},
                smooth=Smooth.None,
                pattern=LinePattern.Dash)}));

      end Double;

      partial class Single "Icon for a single-connector boundary condition"
        //extends Names.Middle;
        annotation (Icon(graphics={
              Rectangle(
                extent={{-100,40},{100,-40}},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid,
                pattern=LinePattern.None),
              Line(
                points={{-100,-40},{-100,40},{100,40},{100,-40}},
                pattern=LinePattern.None,
                smooth=Smooth.None),
              Line(
                points={{-100,-40},{100,-40}},
                color={0,0,0},
                smooth=Smooth.None,
                pattern=LinePattern.Dash),
              Text(
                extent={{-100,-20},{100,20}},
                textString="%name",
                lineColor={0,0,0})}));

      end Single;

    end Icons;

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
              lineColor={0,0,0}), Text(
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
  annotation (Documentation(info="<html><p>**Since the connectors in
<a href=\"modelica://FCSys\">FCSys</a> are hierarchical
(see the <a href=\"modelica://FCSys.Connectors\">Connectors</a> package),
the models for the boundary conditions must be as well.  A
<a href=\"modelica://FCSys.Connectors.Face\">Face</a>
connector (<a href=\"modelica://FCSys.Connectors.Face\">Face</a>,
<a href=\"modelica://FCSys.Connectors.Face\">Face</a>, or
<a href=\"modelica://FCSys.Connectors.Face\">Face</a>)
is used in <a href=\"modelica://FCSys.Subregions.Species\">Species</a> models,
and there is a corresponding <a href=\"modelica://FCSys.Conditions.Face.Species.Species\">Species
boundary condition</a> model in this package. The
<a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a>
connector is used in <a href=\"modelica://FCSys.Subregions.Phases\">Phase</a> models,
and there are corresponding <a href=\"modelica://FCSys.Conditions.Face.Phases\">Phase
boundary condition</a> models.  The
<a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector is nested once more
in models such as the <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a>,
and there is a corresponding <a href=\"modelica://FCSys.Conditions.Face.Subregion\">Subregion
boundary condition</a> model.
</p>

<p>**The <a href=\"modelica://FCSys.Conditions.Chemical\">Chemical</a>,
  <a href=\"modelica://FCSys.Conditions.ChemicalBus\">ChemicalBus</a>, <a href=\"modelica://FCSys.Conditions.InertAmagat\">Inert</a>,
  <a href=\"modelica://FCSys.Conditions.InertAmagat\">InertAmagat</a>,
  <a href=\"modelica://FCSys.Conditions.InertDalton\">InertDalton</a>, <a href=\"modelica://FCSys.Conditions.Face\">Face</a>, and
  <a href=\"modelica://FCSys.Conditions.FaceBus\">FaceBus</a> packages contain models to impose boundary conditions on the
  connectors with the same names (<a href=\"modelica://FCSys.Connectors.ChemicalInput\">ChemicalInput</a> or
  <a href=\"modelica://FCSys.Connectors.ChemicalOutput\">ChemicalOutput</a>, <a href=\"modelica://FCSys.Connectors.Inert\">Inert</a> or
  <a href=\"modelica://FCSys.Connectors.InertInternal\">InertInternal</a>, <a href=\"modelica://FCSys.Conditions.InertAmagat\">InertAmagat</a>,
  <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a>, <a href=\"modelica://FCSys.Conditions.Face\">Face</a>, and
  <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a>).  The <a href=\"modelica://FCSys.Conditions.FaceDifferential\">FacePair</a> and
  <a href=\"modelica://FCSys.Conditions.FaceBusPair\">FaceBusPair</a> packages contain models for pairs of <a href=\"modelica://FCSys.Conditions.Face\">Face</a> and
  <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connectors.

</p></html>"));

end Conditions;
